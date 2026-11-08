# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Analytical Laue, harmonic, spectrum, export and monochromatic consistency checks."""
import csv
import json
import math
import os
import tempfile
import threading
import unittest
from app import DEFAULTS, calculate, export_laue
from geometry import detector_geometry, HC_KEV_ANGSTROM, omega_solutions, detector_hit
from laue import simulate_laue_spots, group_laue_spots, weight_laue_spots
from spectrum import load_spectrum, evaluate_spectrum, source_spectrum, TUBE_VOLTAGES

IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


class LaueTests(unittest.TestCase):
    def setUp(self):
        self.detector = detector_geometry(100, 4000, 4000, 100, 100, 1999.5, 1999.5)
        self.records = [dict(hkl=(-1, 3, 0), g=(-1, 3, 0), d=1/math.sqrt(10), f2=4),
                        dict(hkl=(-2, 6, 0), g=(-2, 6, 0), d=1/math.sqrt(40), f2=9)]

    def events(self, **kwargs):
        options = dict(omega=0, detector=self.detector, energy_min=1, energy_max=150)
        options.update(kwargs)
        return simulate_laue_spots(self.records, [IDENTITY], **options)

    def test_analytical_energy_position_and_harmonics(self):
        events = self.events()
        self.assertEqual(len(events), 2)
        self.assertAlmostEqual(events[0]['wavelength_A'], 0.2)
        self.assertAlmostEqual(events[0]['energy_keV'], HC_KEV_ANGSTROM/0.2)
        self.assertAlmostEqual(events[0]['col'], 2749.5)
        self.assertAlmostEqual(events[0]['row'], 1999.5)
        self.assertAlmostEqual(events[1]['energy_keV'], 2*events[0]['energy_keV'])
        self.assertAlmostEqual(events[1]['col'], events[0]['col'])
        weighted = weight_laue_spots(events, dict(kind='flat'), None)
        groups = group_laue_spots(weighted)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]['intensity'], 13)
        self.assertEqual(len(groups[0]['members']), 2)

    def test_mono_geometry_consistency_tilted_nonzero_omega(self):
        detector = detector_geometry(100, 4000, 4000, 100, 100, 1999.5, 1999.5, 5, -8, 9)
        events = self.events(omega=13, detector=detector)
        self.assertTrue(events)
        for record, event in zip(self.records, events):
            hit = detector_hit(record['g'], 13, event['wavelength_A'], detector)
            self.assertAlmostEqual(hit[0], event['col'], places=9)
            self.assertAlmostEqual(hit[1], event['row'], places=9)
            roots = omega_solutions(record['g'], event['wavelength_A'], 0, 30)
            self.assertTrue(any(abs(angle - 13) < 1e-8 for angle in roots))

    def test_band_sign_angle_and_cancellation(self):
        self.assertEqual(len(self.events(energy_max=80)), 1)
        self.assertEqual(self.events(energy_min=1, energy_max=2), [])
        self.assertEqual(self.events(max_two_theta=10), [])
        self.assertEqual(self.events(omega=180), [])
        cancel = threading.Event()
        cancel.set()
        with self.assertRaises(RuntimeError):
            self.events(cancel=cancel)
        for options in (dict(energy_min=80, energy_max=40), dict(omega=float('nan'))):
            with self.assertRaises(ValueError):
                self.events(**options)

    def test_spectrum_interpolation_and_validation(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, 's.csv')
            with open(path, 'w') as stream:
                stream.write('# Example\nenergy_keV,photon_density\n10,2\n20,6\n')
            spectrum = load_spectrum(path)
            self.assertEqual(evaluate_spectrum(spectrum, [9, 10, 15, 20, 21]), [0, 2, 4, 6, 0])
            for rows in ('10,1\n10,2', '10,-1\n20,1', '10,nan\n20,1', '10,0\n20,0'):
                with open(path, 'w') as stream:
                    stream.write('energy_keV,photon_density\n' + rows)
                with self.assertRaises(ValueError):
                    load_spectrum(path)

    def test_reference_spectra(self):
        for target, voltages in TUBE_VOLTAGES.items():
            for voltage in voltages:
                spectrum = source_spectrum(dict(source=target, tube_voltage=voltage))
                self.assertEqual(spectrum['energies'][-1], float(voltage))
                self.assertEqual(evaluate_spectrum(spectrum, [float(voltage) + 1]), [0])
                self.assertTrue(max(spectrum['densities']) > 0)
        with self.assertRaises(ValueError):
            source_spectrum(dict(source='Cr', tube_voltage='120'))

    def test_pipeline_export_and_cache(self):
        settings = dict(DEFAULTS, beam_mode='Polychromatic (Laue)', mode='Single orientation', source='Flat')
        result = calculate(settings)
        self.assertTrue(result['events'])
        self.assertTrue(any(len(group['members']) > 1 for group in result['groups']))
        shifted = calculate(dict(settings, omega='5'))
        self.assertIs(result['reflections'], shifted['reflections'])
        self.assertTrue(all(event['lorentz'] is None for event in result['events']))
        self.assertNotEqual(result['events'][0]['col'], shifted['events'][0]['col'])
        uniform = calculate(dict(settings, laue_weight='Uniform'))
        self.assertTrue(all(group['intensity'] == 1 for group in uniform['groups']))
        with tempfile.TemporaryDirectory() as folder:
            base = os.path.join(folder, 'export')
            export_laue(result, base)
            with open(base + '_laue.csv') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), len(result['events']))
            self.assertIn('energy_keV', rows[0])
            with open(base + '_laue_settings.json') as stream:
                metadata = json.load(stream)
            self.assertEqual(metadata['matrices'], result['matrices'])
            self.assertEqual(len(metadata['euler_bunge_deg'][0]), 3)
            self.assertIn('spectrum', metadata)

    def test_invalid_powder_and_orientation(self):
        with self.assertRaises(ValueError):
            calculate(dict(DEFAULTS, beam_mode='Polychromatic (Laue)'))
        with self.assertRaises(ValueError):
            calculate(dict(DEFAULTS, beam_mode='Polychromatic (Laue)', mode='Single orientation', single_euler='1,181,0'))


if __name__ == '__main__':
    unittest.main()
