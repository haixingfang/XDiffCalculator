# -*- coding: utf-8 -*-
# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Deterministic crystallography, geometry, orientation and workflow checks."""
from __future__ import division

import math
import json
import os
import random
import tempfile
import threading
import unittest
import gemmi

from app import DEFAULTS, calculate, settings_for_disk, settings_from_disk, example_cifs, intensity_color, ring_segments
from diffraction import cif_blocks, crystal_metadata, load_crystal, reflections, powder_lines, simulate_spots, polarization_factor, lorentz_factor
from geometry import HC_KEV_ANGSTROM, detector_geometry, maximum_two_theta, omega_solutions, detector_hit, matvec
from orientations import random_orientations, validate_rotation, load_orientations, save_orientations

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class CalculatorTests(unittest.TestCase):
    """Check ideal analytical limits rather than comparing only with the implementation."""

    @classmethod
    def setUpClass(cls):
        """Prepare one BCC and one FCC crystal and an asymmetric-pitch detector."""
        cls.fe = load_crystal(os.path.join(BASE, 'examples', 'Fe_bcc.cif'))
        cls.al = load_crystal(os.path.join(BASE, 'examples', 'Al_fcc.cif'))
        cls.fe_records = reflections(cls.fe, 40, 25)
        cls.detector = detector_geometry(100, 2001, 1501, 150, 200, 1000, 750)

    def test_cif_block_and_missing_block(self):
        """Structure blocks must be explicit and errors must not fall back silently."""
        path = os.path.join(BASE, 'examples', 'Fe_bcc.cif')
        self.assertEqual(cif_blocks(path), ['Fe_bcc'])
        with self.assertRaises(ValueError):
            load_crystal(path, 'absent')

    def test_crystal_metadata(self):
        """Resolve known cubic space groups and non-cubic cell-only metadata."""
        for name, sgno in [('Fe_bcc.cif', 229), ('Al_fcc.cif', 225)]:
            metadata = crystal_metadata(os.path.join(BASE, 'examples', name))
            self.assertEqual(metadata['sgno'], sgno)
            self.assertEqual(metadata['cell'][3:], [90, 90, 90])
        text = ('data_first\n_cell_length_a 4\n_cell_length_b 5\n_cell_length_c 6\n'
                '_cell_angle_alpha 80\n_cell_angle_beta 95\n_cell_angle_gamma 110\n_space_group_name_H-M_alt "P 1"\n')
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, 'cells.cif')
            with open(path, 'w') as stream:
                stream.write(text + text.replace('data_first', 'data_second').replace('_cell_length_a 4', '_cell_length_a 7'))
            with self.assertRaises(ValueError):
                crystal_metadata(path)
            self.assertEqual(crystal_metadata(path, 'first'), dict(cell=[4, 5, 6, 80, 95, 110], sgno=1))
            self.assertEqual(crystal_metadata(path, 'second')['cell'][0], 7)
            with self.assertRaises(ValueError):
                crystal_metadata(path, 'missing')

    def test_bcc_extinctions_spacing_and_multiplicity(self):
        """BCC 110 has multiplicity 12 and d=a/sqrt(2); odd h+k+l is absent."""
        self.assertTrue(all(sum(record['hkl']) % 2 == 0 for record in self.fe_records))
        first = powder_lines(self.fe_records)[0]
        self.assertEqual(first['hkl'], (1, 1, 0))
        self.assertEqual(first['multiplicity'], 12)
        self.assertAlmostEqual(first['d'], self.fe.cell.a/math.sqrt(2), places=12)

    def test_bcc_absolute_structure_factor(self):
        """Special-position occupancy must give F110=2*fFe, not 96*fFe."""
        record = next(item for item in self.fe_records if item['hkl'] == (1, 1, 0))
        atomic = gemmi.Element('Fe').it92.calculate_sf(1/(4*record['d']**2))
        self.assertAlmostEqual(record['f2']/(4*atomic**2), 1, places=6)

    def test_fcc_structure_factor(self):
        """FCC 111 has multiplicity 8 and F=4*fAl."""
        records = reflections(self.al, 40, 20)
        self.assertTrue(all(len(set(abs(value) % 2 for value in item['hkl'])) == 1 for item in records))
        first = powder_lines(records)[0]
        self.assertEqual(first['hkl'], (1, 1, 1))
        self.assertEqual(first['multiplicity'], 8)
        atomic = gemmi.Element('Al').it92.calculate_sf(1/(4*first['d']**2))
        self.assertAlmostEqual(first['sum_f2']/(8*16*atomic**2), 1, places=6)
        self.assertEqual(first['index'], 0)
        self.assertEqual(first['hkl'], (1, 1, 1))

    def test_occupancy_and_displacement(self):
        """Half occupancy reduces F squared by four; positive U reduces intensity."""
        crystal = load_crystal(os.path.join(BASE, 'examples', 'Fe_bcc.cif'))
        calc = gemmi.StructureFactorCalculatorX(crystal.cell)
        initial = abs(calc.calculate_sf_from_small_structure(crystal, (1, 1, 0)))**2
        crystal.sites[0].occ *= 0.5
        half = abs(calc.calculate_sf_from_small_structure(crystal, (1, 1, 0)))**2
        self.assertAlmostEqual(half/initial, 0.25, places=12)
        crystal.sites[0].u_iso = 0.02
        self.assertLess(abs(calc.calculate_sf_from_small_structure(crystal, (1, 1, 0)))**2, half)

    def test_triclinic_reciprocal_basis(self):
        """The inverse-transpose direct basis must reproduce Gemmi's d spacings."""
        cell = gemmi.UnitCell(3, 4, 5, 73, 84, 108)
        frac = cell.frac.mat.tolist()
        basis = [[frac[j][i] for j in range(3)] for i in range(3)]
        for hkl in [(1, 2, -3), (0, 1, 0), (-2, 1, 4)]:
            g = matvec(basis, hkl)
            self.assertAlmostEqual(math.sqrt(sum(x*x for x in g)), 1/cell.calculate_d(hkl), places=12)

    def test_random_rotations_and_reproducibility(self):
        """Uniform rotations are proper, reproducible, and approximately isotropic."""
        matrices = random_orientations(4000, 14)
        self.assertEqual(matrices[:4], random_orientations(4, 14))
        for matrix in matrices:
            validate_rotation(matrix, 1e-12)
        for axis in range(3):
            mean = sum(matrix[axis][0] for matrix in matrices)/len(matrices)
            moment = sum(matrix[axis][0]**2 for matrix in matrices)/len(matrices)
            self.assertLess(abs(mean), 0.035)
            self.assertLess(abs(moment - 1/3), 0.035)

    def test_invalid_rotation(self):
        """Reject reflections, scale errors and NaN rather than repairing matrices."""
        for matrix in [[[1, 0, 0], [0, 1, 0], [0, 0, -1]], [[2, 0, 0], [0, 1, 0], [0, 0, 1]], [[float('nan'), 0, 0], [0, 1, 0], [0, 0, 1]]]:
            with self.assertRaises(ValueError):
                validate_rotation(matrix)

    def test_orientation_csv_roundtrip(self):
        """Exported row-major matrices can be loaded without convention changes."""
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, 'orientations.csv')
            matrices = random_orientations(5, 11)
            save_orientations(path, matrices)
            loaded = load_orientations(path)
            for a, b in zip(matrices, loaded):
                self.assertLess(max(abs(a[i][j] - b[i][j]) for i in range(3) for j in range(3)), 1e-13)

    def test_analytic_bragg_and_rotation_sign(self):
        """g=(1,0,0), lambda=1 diffracts at 120/240 degrees with opposite Y."""
        roots = omega_solutions((1, 0, 0), 1, 0, 360)
        self.assertEqual(len(roots), 2)
        self.assertAlmostEqual(roots[0], 120)
        self.assertAlmostEqual(roots[1], 240)
        detector = detector_geometry(1, 1001, 1001, 10, 10, 500, 500)
        hit = detector_hit((1, 0, 0), roots[0], 1, detector)
        self.assertGreater(hit[0], detector['cx'])
        self.assertAlmostEqual(hit[1], detector['cy'])

    def test_bragg_energy_and_ring_radius_randomized(self):
        """Random solutions satisfy elastic scattering and r=L*tan(2theta)."""
        rng = random.Random(18)
        wavelength = HC_KEV_ANGSTROM/40
        checked = 0
        for unused in range(300):
            g = tuple(rng.uniform(-1, 1) for axis in range(3))
            norm2 = sum(x*x for x in g)
            for omega in omega_solutions(g, wavelength, -80, 280):
                radians = math.radians(omega)
                gx = math.cos(radians)*g[0] - math.sin(radians)*g[1]
                self.assertAlmostEqual(2*wavelength*gx + wavelength*wavelength*norm2, 0, places=12)
                hit = detector_hit(g, omega, wavelength, self.detector)
                if hit:
                    radius = math.hypot((hit[0] - self.detector['cx'])*self.detector['py'], (hit[1] - self.detector['cy'])*self.detector['pz'])
                    tth = 2*math.asin(wavelength*math.sqrt(norm2)/2)
                    self.assertAlmostEqual(radius, self.detector['distance']*math.tan(tth), places=10)
                    if g[2] > 0:
                        self.assertLess(hit[1], self.detector['cy'])
                    checked += 1
        self.assertGreater(checked, 200)

    def test_boundaries_and_impossible_reflections(self):
        """Frame windows are half-open, periodic roots are unique, axial g is absent."""
        self.assertEqual(len(omega_solutions((1, 0, 0), 1, 120, 240)), 1)
        self.assertEqual(omega_solutions((0, 0, 1), 1, 0, 360), [])
        self.assertEqual(omega_solutions((10, 0, 0), 1, 0, 360), [])
        self.assertEqual(len(omega_solutions((2, 0, 0), 1, 0, 360)), 1)

    def test_frames_and_no_multiplicity_in_spots(self):
        """One explicit reflection creates two events, not a powder-multiplied set."""
        record = dict(hkl=(1, 0, 0), g=(1, 0, 0), d=1, two_theta=60, f2=7)
        detector = detector_geometry(1, 1001, 1001, 10, 10, 500, 500)
        identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        events = simulate_spots([record], [identity], HC_KEV_ANGSTROM, detector, 0, 360, 120)
        self.assertEqual([event['frame'] for event in events], [1, 2])
        self.assertEqual([event['f2'] for event in events], [7, 7])

    def test_detector_validation_and_coverage(self):
        """Validate detector geometry including finite coordinates and square corners."""
        detector = detector_geometry(100, 1000, 1000, 100, 100, 499.5, 499.5)
        self.assertAlmostEqual(maximum_two_theta(detector), math.degrees(math.atan(math.sqrt(2)/2)), places=12)
        with self.assertRaises(ValueError):
            detector_geometry(0, 10, 10, 100, 100, 5, 5)
        with self.assertRaises(ValueError):
            detector_geometry(100, 1.5, 10, 100, 100, 5, 5)

    def test_polarization_limits_and_average(self):
        """Thomson limits: horizontal scattering is suppressed by horizontal polarization."""
        self.assertAlmostEqual(polarization_factor(60, 0, 1), 0.25)
        self.assertAlmostEqual(polarization_factor(60, math.pi/2, 1), 1)
        self.assertAlmostEqual(polarization_factor(60, 0, -1), 1)
        self.assertAlmostEqual(polarization_factor(60, 0, 0), 0.625)
        self.assertEqual(polarization_factor(60, 0, None), 1)
        mean = sum(polarization_factor(60, math.radians(i), 0.95) for i in range(360))/360
        self.assertAlmostEqual(mean, polarization_factor(60, polarization=0.95))
        for bad in (2, -2, float('nan')):
            with self.assertRaises(ValueError):
                polarization_factor(60, 0, bad)

    def test_lorentz_geometry(self):
        """Check powder and equatorial rotation limits and an off-equator ray."""
        self.assertAlmostEqual(lorentz_factor(60), math.cos(math.pi/6)/0.75)
        self.assertAlmostEqual(lorentz_factor(60, 0), 2/math.sqrt(3))
        self.assertAlmostEqual(lorentz_factor(60, math.pi/3), 4/math.sqrt(3))
        with self.assertRaises(ValueError):
            lorentz_factor(60, math.pi/2)

    def test_corrected_intensity_workflows(self):
        """Corrections change weights, preserve geometry, and multiply rather than divide."""
        settings = dict(DEFAULTS, mode='Random grains', number='3', max_tth='20')
        corrected = calculate(settings)
        disabled = calculate(dict(settings, polarization=''))
        unpolarized = calculate(dict(settings, polarization='0'))
        self.assertEqual(len(corrected['events']), len(disabled['events']))
        for spot, no_p, unpol in zip(corrected['events'], disabled['events'], unpolarized['events']):
            self.assertEqual((spot['omega'], spot['col'], spot['row']), (no_p['omega'], no_p['col'], no_p['row']))
            self.assertAlmostEqual(spot['intensity'], spot['f2']*spot['polarization']*spot['lorentz'])
            self.assertAlmostEqual(no_p['intensity'], no_p['f2']*no_p['lorentz'])
            self.assertAlmostEqual(unpol['polarization'], (1 + math.cos(math.radians(spot['two_theta']))**2)/2)
            self.assertAlmostEqual(spot['lorentz'], 1/abs(math.sin(math.radians(spot['two_theta']))*math.cos(spot['azimuth'])))
        self.assertEqual(corrected['intensity_max'], max(spot['intensity'] for spot in corrected['events']))
        for line in corrected['lines']:
            self.assertAlmostEqual(line['intensity'], line['sum_f2']*line['polarization']*line['lorentz'])

    def test_ring_colours(self):
        """Horizontal and vertical ring sectors differ; guides do not imply intensity."""
        result = calculate(DEFAULTS)
        line = result['lines'][0]
        segments = list(ring_segments(line, result['detector'], result['polarization'], result['intensity_max']))
        self.assertEqual(len(segments), 360)
        self.assertNotEqual(segments[0][-1], segments[90][-1])
        guides = list(ring_segments(line, result['detector'], result['polarization'], result['intensity_max'], True))
        self.assertEqual(len(set(segment[-1] for segment in guides)), 1)
        self.assertEqual(intensity_color(0), '#1e3791')
        self.assertEqual(intensity_color(1), '#fff046')

    def test_cancellation(self):
        """Cancellation must stop a calculation, not return a partial result."""
        cancel = threading.Event()
        cancel.set()
        with self.assertRaises(RuntimeError):
            reflections(self.fe, 40, 20, cancel)

    def test_full_workflows(self):
        """Powder and broad-scan random/imported grain workflows produce results."""
        for mode in ['Powder rings', 'Random grains', 'Imported orientations']:
            settings = dict(DEFAULTS, mode=mode, number='3', max_tth='20', end='180', orientation_file=os.path.join(BASE, 'examples', 'orientations.csv'))
            result = calculate(settings)
            self.assertTrue(result['lines'])
            if mode != 'Powder rings':
                self.assertTrue(result['events'])
        settings = dict(DEFAULTS, energy='nan')
        with self.assertRaises(ValueError):
            calculate(settings)

    def test_portable_example_settings(self):
        """Packaged examples must not retain temporary single-file extraction paths."""
        saved = settings_for_disk(DEFAULTS)
        self.assertEqual(saved['cif'], '@example/Fe_bcc.cif')
        restored = settings_from_disk(saved)
        self.assertEqual(restored['cif'], DEFAULTS['cif'])
        with self.assertRaises(ValueError):
            settings_from_disk(dict(cif='@example/../bad.cif'))

    def test_detector_presets(self):
        """Verify full EIGER image dimensions and user-specified effective camera pitches."""
        with open(os.path.join(BASE, 'detectors.json'), encoding='utf-8') as stream:
            presets = json.load(stream)
        expected = {'EIGER2 X CdTe 4M': (2068, 2162, 75), 'EIGER2 X CdTe 16M': (4148, 4362, 75),
                    'marana_20x_0.55um': (2048, 2048, 0.55), 'marana_10x_1.1um': (2048, 2048, 1.1),
                    'marana_7.5x_1.47um': (2048, 2048, 1.47), 'marana_5x_2.2um': (2048, 2048, 2.2),
                    'Hamamatsu_2x_3.25um': (2048, 2048, 3.25)}
        for name, (columns, rows, pitch) in expected.items():
            preset = presets[name]
            self.assertEqual((preset['columns'], preset['rows'], preset['pixel_y'], preset['pixel_z']), (columns, rows, pitch, pitch))
            detector = detector_geometry(100, columns, rows, pitch, pitch, (columns - 1)/2, (rows - 1)/2)
            self.assertAlmostEqual(detector['py'], pitch/1000, places=12)
            self.assertAlmostEqual(detector['columns']*detector['py'], columns*pitch/1000, places=12)

    def test_example_discovery(self):
        """Discover every CIF, not just the initial Fe/Al examples; preserve filename case."""
        names = example_cifs()
        expected = sorted([name for name in os.listdir(os.path.join(BASE, 'examples'))
                           if name.lower().endswith('.cif') and os.path.isfile(os.path.join(BASE, 'examples', name))], key=str.lower)
        self.assertEqual(names, expected)
        self.assertIn('Fe_bcc.cif', names)
        self.assertNotIn('orientations.csv', names)
        self.assertEqual(example_cifs(os.path.join(BASE, 'does_not_exist')), [])

    def test_new_example_settings(self):
        """Every example path survives serialization without opening or altering its CIF."""
        for name in example_cifs():
            settings = dict(DEFAULTS, cif=os.path.join(BASE, 'examples', name))
            encoded = settings_for_disk(settings)
            self.assertEqual(encoded['cif'], '@example/' + name)
            self.assertEqual(os.path.normpath(settings_from_disk(encoded)['cif']), os.path.normpath(settings['cif']))


if __name__ == '__main__':
    unittest.main()
