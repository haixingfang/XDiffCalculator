# -*- coding: utf-8 -*-
# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""CIF reflection calculation and ideal point-grain rotation diffraction."""
from __future__ import division

import math
import gemmi

from geometry import HC_KEV_ANGSTROM, positive, matvec, maximum_two_theta, omega_solutions, detector_hit
from orientations import validate_rotation


def cif_blocks(path):
    """Return CIF data blocks containing unit-cell lengths."""
    return [block.name for block in gemmi.cif.read_file(path) if block.find_value('_cell_length_a')]


def crystal_metadata(path, block_name=None):
    """Read lattice metadata without requiring valid scattering sites.
    Args:
        path: Local CIF file.
        block_name: Unit-cell block; required when several blocks are available.
    Returns:
        metadata: Cell lengths in angstrom, angles in degrees and sgno (or None).
    """
    document = gemmi.cif.read_file(path)
    names = [block.name for block in document if block.find_value('_cell_length_a')]
    if block_name is None and len(names) == 1:
        block_name = names[0]
    if block_name not in names:
        raise ValueError('Select a valid unit-cell block.')
    structure = gemmi.make_small_structure_from_block(document.find_block(block_name))
    cell = structure.cell
    if not all(math.isfinite(value) for value in cell.parameters) or min(cell.a, cell.b, cell.c) <= 0 or not cell.volume > 0:
        raise ValueError('Invalid unit cell.')
    return dict(cell=list(cell.parameters), sgno=structure.spacegroup.number if structure.spacegroup else None)


def load_crystal(path, block_name=None):
    """Read a small-structure CIF with symmetry, sites, occupancies and ADPs.
    Args:
        path: Local CIF file; atom-site fractional coordinates are required.
        block_name: Selected data block, required if several structures exist.
    Returns:
        structure: Gemmi structure with crystallographic special-position weights.
    """
    document = gemmi.cif.read_file(path)
    names = [block.name for block in document if block.find_value('_cell_length_a')]
    if not names:
        raise ValueError('CIF has no unit-cell block.')
    if block_name is None:
        if len(names) != 1:
            raise ValueError('Select a CIF block: {}.'.format(', '.join(names)))
        block_name = names[0]
    if block_name not in names:
        raise ValueError('Unknown structure block: {}.'.format(block_name))
    structure = gemmi.make_small_structure_from_block(document.find_block(block_name))
    if structure.spacegroup is None:
        raise ValueError('CIF space group could not be resolved; do not assume P1.')
    if not structure.sites:
        raise ValueError('CIF must contain fractional atomic coordinates to calculate structure factors.')
    cell = structure.cell
    if not all(math.isfinite(value) for value in cell.parameters) or min(cell.a, cell.b, cell.c) <= 0 or not cell.volume > 0:
        raise ValueError('Invalid unit cell.')
    for site in structure.sites:
        values = (site.fract.x, site.fract.y, site.fract.z, site.occ, site.u_iso)
        if not all(math.isfinite(value) for value in values) or not 0 <= site.occ <= 1:
            raise ValueError('Invalid coordinates, occupancy or displacement for {}.'.format(site.label))
        if site.element.atomic_number == 0 or site.element.it92 is None:
            raise ValueError('Unsupported atomic scattering factor for {}.'.format(site.label))
    structure.setup_cell_images()
    structure.change_occupancies_to_crystallographic()
    return structure


def reflections(structure, energy, max_two_theta, cancel=None):
    """Calculate ALL signed HKL directions, without powder multiplicity factors.
    Args:
        structure: Validated structure from load_crystal.
        energy: X-ray energy in keV.
        max_two_theta: Scattering-angle limit in degrees, strictly below 90.
        cancel: Optional threading.Event for cancellation.
    Returns:
        records: HKL, reciprocal Cartesian g (1/A, no 2*pi), d, two_theta, f2.
            Uses neutral-atom f0 and CIF ADPs; anomalous scattering is omitted.
    """
    wavelength = HC_KEV_ANGSTROM/positive(energy, 'Energy')
    angle = positive(max_two_theta, 'Maximum 2theta')
    if angle >= 90:
        raise ValueError('Maximum 2theta must be below 90 degrees.')
    dmin = wavelength/(2*math.sin(math.radians(angle/2)))
    result = enumerate_reflections(structure, dmin, cancel)
    for record in result:
        record['two_theta'] = math.degrees(2*math.asin(min(1, wavelength/(2*record['d']))))
    return sorted(result, key=lambda item: (item['two_theta'], item['hkl']))


def enumerate_reflections(structure, dmin, cancel=None):
    """Enumerate signed HKLs without choosing a wavelength.
    Args:
        structure: Validated Gemmi structure from load_crystal.
        dmin: Minimum spacing in angstrom.
        cancel: Optional threading.Event.
    Returns:
        records: HKL, reciprocal vector (no 2*pi), spacing and neutral-atom F squared.
    """
    dmin = positive(dmin, 'Minimum spacing')
    cell = structure.cell
    limits = [int(math.ceil(length/dmin)) for length in (cell.a, cell.b, cell.c)]
    if (2*limits[0] + 1)*(2*limits[1] + 1)*(2*limits[2] + 1) > 2000000:
        raise ValueError('HKL search exceeds two million indices; reduce energy or angular coverage.')
    # Direct basis: a along +x, b in xy, c with positive z. B = A^(-T).
    fractional = cell.frac.mat.tolist()
    basis = [[fractional[j][i] for j in range(3)] for i in range(3)]
    calculator = gemmi.StructureFactorCalculatorX(cell)
    operations = structure.spacegroup.operations()
    result = []
    for h in range(-limits[0], limits[0] + 1):
        if cancel is not None and cancel.is_set():
            raise RuntimeError('Calculation cancelled.')
        for k in range(-limits[1], limits[1] + 1):
            for l in range(-limits[2], limits[2] + 1):
                hkl = (h, k, l)
                if hkl == (0, 0, 0) or operations.is_systematically_absent(hkl):
                    continue
                d = cell.calculate_d(hkl)
                if d < dmin*(1 - 1e-12):
                    continue
                f2 = abs(calculator.calculate_sf_from_small_structure(structure, hkl))**2
                if f2 <= 1e-16:
                    continue
                result.append(dict(hkl=hkl, g=matvec(basis, hkl), d=d, f2=f2))
    return sorted(result, key=lambda item: (-item['d'], item['hkl']))


def polarization_factor(two_theta, azimuth=None, polarization=0.95):
    """Calculate Thomson polarization for azimuth measured from +Y toward +Z.
    Args:
        two_theta: Scattering angle in degrees.
        azimuth: Detector azimuth in radians; None returns the full-ring average.
        polarization: Signed linear polarization degree [-1,1], or None for P=1.
            +1 is horizontal (+Y), -1 vertical (+Z), 0 unpolarized.
    Returns:
        factor: Multiplicative scattering intensity factor, not its inverse.
    """
    if polarization is None:
        return 1.0
    if not math.isfinite(polarization) or not -1 <= polarization <= 1:
        raise ValueError('Polarization must be between -1 and 1, or blank to disable.')
    sin2 = math.sin(math.radians(two_theta))**2
    modulation = 0 if azimuth is None else polarization*math.cos(2*azimuth)
    return 1 - 0.5*sin2*(1 + modulation)


def lorentz_factor(two_theta, azimuth=None):
    """Return powder-profile L or omega-integrated spot L for rotation about +Z.
    Args:
        two_theta: Scattering angle in degrees, strictly between 0 and 90.
        azimuth: Radians from +Y toward +Z; None selects powder-profile L.
    Returns:
        factor: Powder cos(theta)/sin(2theta)^2; spots 1/|sin(2theta)cos(azimuth)|.
            Exact rotation tangencies raise an error; finite rocking width is needed.
    """
    if not math.isfinite(two_theta) or not 0 < two_theta < 90:
        raise ValueError('Lorentz factor requires 0 < 2theta < 90 degrees.')
    angle = math.radians(two_theta)
    if azimuth is None:
        return math.cos(angle/2)/math.sin(angle)**2
    denominator = abs(math.sin(angle)*math.cos(azimuth))
    if denominator < 1e-10:
        raise ValueError('Singular rotation Lorentz factor: reflection tangent to the Ewald sphere. A finite rocking-profile model is required.')
    return 1/denominator


def powder_lines(records, polarization=0.95):
    """Sum coincident signed reflections into powder lines; retain all HKL members.
    Returns:
        lines: Relative powder-profile weights sum(|F|^2)*P_mean*L, not counts.
            multiplicity counts all coincident directions, including accidental overlaps.
    """
    lines = []
    for record in sorted(records, key=lambda item: item['two_theta']):
        if not lines or abs(record['two_theta'] - lines[-1]['two_theta']) > 1e-7:
            lines.append(dict(two_theta=record['two_theta'], d=record['d'], sum_f2=0.0, hkls=[]))
        lines[-1]['sum_f2'] += record['f2']
        lines[-1]['hkls'].append(record['hkl'])
    for index, line in enumerate(lines):
        line['index'] = index
        line['polarization'] = polarization_factor(line['two_theta'], polarization=polarization)
        line['lorentz'] = lorentz_factor(line['two_theta'])
        line['intensity'] = line['sum_f2']*line['polarization']*line['lorentz']
    maximum = max([line['intensity'] for line in lines] or [1])
    for line in lines:
        line['multiplicity'] = len(line['hkls'])
        line['hkl'] = max(line['hkls'])
        line['relative'] = 100*line['intensity']/maximum
    return lines


def simulate_spots(records, matrices, energy, detector, start, end, step, cancel=None, progress=None, polarization=0.95):
    """Compute exact Bragg events for origin-centred grains in finite omega frames.
    Args:
        records: All signed reflections from reflections(), not powder representatives.
        matrices: Crystal Cartesian to sample rotations at omega=0.
        energy: Energy in keV; must match the reflection calculation.
        detector: Tilted or perpendicular geometry from detector_geometry().
        start, end, step: Degrees; increasing [start,end), at most one full turn.
        cancel: Optional cancellation event.
        progress: Optional callback receiving completed and total grain counts.
        polarization: Signed linear polarization degree, or None to disable P.
    Returns:
        events: Grain/HKL/omega/frame/column/row/f2 dictionaries, sorted by frame.
            intensity=|F|^2*P*L; f2 remains uncorrected. No finite rocking profile.
    """
    start, end, step = float(start), float(end), positive(step, 'Omega step')
    if not math.isfinite(start) or not math.isfinite(end) or not 0 < end - start <= 360:
        raise ValueError('Use increasing omega limits spanning at most 360 degrees.')
    if math.ceil((end - start)/step) > 100000:
        raise ValueError('More than 100000 frames requested; increase the omega step.')
    if len(records)*len(matrices) > 20000000:
        raise ValueError('More than 20 million grain/reflection pairs; reduce grain count or angular coverage.')
    wavelength = HC_KEV_ANGSTROM/positive(energy, 'Energy')
    events = []
    for grain, matrix in enumerate(matrices):
        matrix = validate_rotation(matrix)
        for i, record in enumerate(records):
            if i % 256 == 0 and cancel is not None and cancel.is_set():
                raise RuntimeError('Calculation cancelled.')
            g = matvec(matrix, record['g'])
            for omega in omega_solutions(g, wavelength, start, end):
                hit = detector_hit(g, omega, wavelength, detector)
                if hit is not None:
                    angle = math.radians(omega)
                    azimuth = math.atan2(g[2], math.sin(angle)*g[0] + math.cos(angle)*g[1])
                    p = polarization_factor(record['two_theta'], azimuth, polarization)
                    lorentz = lorentz_factor(record['two_theta'], azimuth)
                    frame = int(math.floor((omega - start)/step + 1e-10))
                    frame = min(frame, int(math.ceil((end - start)/step)) - 1)
                    events.append(dict(grain=grain, hkl=record['hkl'], omega=omega, frame=frame,
                                       col=hit[0], row=hit[1], f2=record['f2'], two_theta=record['two_theta'],
                                       azimuth=azimuth, polarization=p, lorentz=lorentz, intensity=record['f2']*p*lorentz))
                    if len(events) > 300000:
                        raise ValueError('More than 300000 spots; reduce grain count or scan range.')
        if progress is not None and (grain % 20 == 0 or grain + 1 == len(matrices)):
            progress(grain + 1, len(matrices))
    return sorted(events, key=lambda event: (event['frame'], event['grain'], event['omega']))
