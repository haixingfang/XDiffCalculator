# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Exact stationary Laue geometry and explicitly approximate spectral weights."""
import math
import os
from functools import lru_cache
from diffraction import load_crystal, enumerate_reflections, polarization_factor
from geometry import HC_KEV_ANGSTROM, matvec, ray_detector_hit
from orientations import validate_rotation
from spectrum import evaluate_spectrum


@lru_cache(maxsize=2)
def _cached_reflections(path, block, modified_ns, size, dmin):
    """Cache completed HKL enumeration; file revision and spacing are part of the key."""
    return enumerate_reflections(load_crystal(path, block), dmin)


def laue_reflections(file_path, block, energy_max, max_two_theta):
    """Get cached energy-independent HKLs bounded by highest energy and angular coverage."""
    if not math.isfinite(energy_max) or energy_max <= 0 or not 0 < max_two_theta < 90:
        raise ValueError('Require positive maximum energy and 0 < maximum 2theta < 90 degrees.')
    dmin = HC_KEV_ANGSTROM/energy_max/(2*math.sin(math.radians(max_two_theta/2)))
    path = os.path.abspath(file_path)
    stat = os.stat(path)
    return _cached_reflections(path, block, stat.st_mtime_ns, stat.st_size, dmin)


def simulate_laue_spots(records, matrices, omega, detector, energy_min, energy_max, max_two_theta=89.9, cancel=None, progress=None):
    """Solve exact wavelengths at one fixed omega and intersect the existing detector.
    Args:
        records: Signed HKLs with crystal Cartesian g in inverse angstrom, without 2*pi.
        matrices: Crystal-to-sample rotations at omega zero. All grains are at the origin.
        omega: Fixed angle in degrees, positive CCW viewed from +Z.
        detector: Existing detector_geometry dictionary, including tilts.
        energy_min, energy_max: Inclusive positive energy limits in keV.
        max_two_theta: Angular acceptance in degrees, less than 90.
        cancel, progress: Optional cancellation event and grain progress callback.
    Returns:
        events: One record per grain/signed HKL, retaining harmonic contributions.
    """
    if not all(math.isfinite(x) for x in (omega, energy_min, energy_max, max_two_theta)) or not 0 < energy_min < energy_max or not 0 < max_two_theta < 90:
        raise ValueError('Require finite omega, 0 < Emin < Emax and 0 < maximum 2theta < 90.')
    if len(records)*len(matrices) > 20000000:
        raise ValueError('More than 20 million grain/HKL pairs; reduce energy, coverage or grain count.')
    cosine, sine = math.cos(math.radians(omega)), math.sin(math.radians(omega))
    events = []
    for grain, matrix in enumerate(matrices):
        matrix = validate_rotation(matrix)
        for index, record in enumerate(records):
            if index % 256 == 0 and cancel is not None and cancel.is_set():
                raise RuntimeError('Calculation cancelled.')
            x, y, z = matvec(matrix, record['g'])
            gx, gy, gz = cosine*x - sine*y, sine*x + cosine*y, z
            norm2 = gx*gx + gy*gy + gz*gz
            if gx >= 0 or norm2 == 0:
                continue
            wavelength = -2*gx/norm2
            energy = HC_KEV_ANGSTROM/wavelength
            if not energy_min <= energy <= energy_max:
                continue
            ray = (1 + wavelength*gx, wavelength*gy, wavelength*gz)
            two_theta = math.degrees(math.acos(max(-1, min(1, ray[0]))))
            if two_theta > max_two_theta:
                continue
            hit = ray_detector_hit(ray, detector)
            if hit is None:
                continue
            events.append(dict(grain=grain, hkl=record['hkl'], d=record['d'], omega=omega, frame=0, col=hit[0], row=hit[1],
                               energy_keV=energy, wavelength_A=wavelength, two_theta=two_theta, azimuth=math.atan2(gz, gy), f2=record['f2']))
            if len(events) > 300000:
                raise ValueError('More than 300000 Laue reflections; reduce grain count or energy range.')
        if progress is not None:
            progress(grain + 1, len(matrices))
    return events


def weight_laue_spots(events, spectrum, polarization=0.0):
    """Return copies weighted by S_E*F_squared*P, not quantitative Laue intensities.
    Args:
        events: Exact Laue event records.
        spectrum: Photon spectral density per keV from spectrum.py.
        polarization: Signed linear polarization; None disables its correction.
    Returns:
        weighted: Events with spectral_weight and intensity aliases for plotting.
            No rotation Lorentz factor, Laue Jacobian, absorption, anomalous
            scattering or detector response is included. Line bins are interpolated.
    """
    weighted = []
    for event, density in zip(events, evaluate_spectrum(spectrum, [e['energy_keV'] for e in events])):
        p = polarization_factor(event['two_theta'], event['azimuth'], polarization)
        weight = density*event['f2']*p
        weighted.append(dict(event, spectral_density=density, polarization=p, spectral_weight=weight, intensity=weight, lorentz=None))
    return weighted


def group_laue_spots(events):
    """Combine exact harmonic rays per grain; keep constituent records for hover/export.
    Returns:
        groups: Summed weights and weight-averaged energy for each primitive HKL ray.
            Different grains remain separate even if detector coordinates coincide.
    """
    groups = {}
    for event in events:
        h, k, l = event['hkl']
        divisor = math.gcd(math.gcd(abs(h), abs(k)), abs(l))
        key = (event['grain'], h//divisor, k//divisor, l//divisor)
        groups.setdefault(key, []).append(event)
    result = []
    for members in groups.values():
        weight = sum(e['intensity'] for e in members)
        energy = sum(e['energy_keV']*e['intensity'] for e in members)/weight if weight else sum(e['energy_keV'] for e in members)/len(members)
        result.append(dict(members[0], members=members, intensity=weight, energy_keV=energy))
    return result
