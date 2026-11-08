# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Lightweight photon-spectrum import and evaluation; no spectrum-generator dependency."""
import bisect
import csv
import hashlib
import math
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
TUBE_VOLTAGES = {'W': ('40', '80', '120'), 'Cr': ('40', '50'), 'Cu': ('40', '50')}


def load_spectrum(file_path):
    """Read photon spectral density, preserving comments as provenance.
    Args:
        file_path: CSV with energy_keV, photon_density header and increasing samples.
            Density may use relative units but must be per keV, not per energy bin.
    Returns:
        spectrum: Energies, densities, comments and SHA256 of the original CSV.
    """
    with open(file_path, 'rb') as stream:
        digest = hashlib.sha256(stream.read()).hexdigest()
    with open(file_path, encoding='utf-8-sig') as stream:
        lines = stream.readlines()
    comments = [line.strip()[1:].strip() for line in lines if line.lstrip().startswith('#')]
    reader = csv.reader(line for line in lines if line.strip() and not line.lstrip().startswith('#'))
    if next(reader, None) != ['energy_keV', 'photon_density']:
        raise ValueError('Spectrum CSV header must be energy_keV,photon_density.')
    energies, densities = [], []
    for row in reader:
        if len(row) != 2:
            raise ValueError('Spectrum rows must have exactly two columns.')
        energy, density = map(float, row)
        if not math.isfinite(energy) or energy <= 0 or not math.isfinite(density) or density < 0:
            raise ValueError('Spectrum requires positive finite energies and nonnegative finite densities.')
        if energies and energy <= energies[-1]:
            raise ValueError('Spectrum energies must be strictly increasing; duplicate bins are not allowed.')
        energies.append(energy)
        densities.append(density)
    if len(energies) < 2 or not any(densities):
        raise ValueError('Spectrum requires at least two samples and some positive flux.')
    return dict(kind='sampled', energies=energies, densities=densities, comments=comments, sha256=digest, path=os.path.abspath(file_path))


def source_spectrum(settings):
    """Select a flat spectrum, bundled reference tube, or measured CSV."""
    source = settings.get('source', 'Flat')
    if source == 'Flat':
        return dict(kind='flat', comments=['Constant relative photon density per keV within the selected band.'])
    if source == 'CSV':
        return load_spectrum(settings.get('spectrum_file', ''))
    if source not in ('W', 'Cr', 'Cu'):
        raise ValueError('Unknown spectrum source.')
    voltage = settings.get('tube_voltage', '80')
    if str(voltage) not in TUBE_VOLTAGES[source]:
        raise ValueError('Available {} reference voltages: {} kV. Import CSV for other conditions.'.format(source, ', '.join(TUBE_VOLTAGES[source])))
    return load_spectrum(os.path.join(ROOT, 'spectra', '{}_{}kV.csv'.format(source, voltage)))


def evaluate_spectrum(spectrum, energies):
    """Evaluate photon density per keV; sampled spectra are zero outside their support."""
    result = []
    for energy in energies:
        if not math.isfinite(energy) or energy <= 0:
            raise ValueError('Photon energies must be positive and finite.')
        if spectrum['kind'] == 'flat':
            result.append(1.0)
            continue
        grid, values = spectrum['energies'], spectrum['densities']
        if energy < grid[0] or energy > grid[-1]:
            result.append(0.0)
            continue
        index = min(len(grid) - 2, max(0, bisect.bisect_right(grid, energy) - 1))
        fraction = (energy - grid[index])/(grid[index + 1] - grid[index])
        result.append(values[index]*(1 - fraction) + values[index + 1]*fraction)
    return result
