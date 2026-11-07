# -*- coding: utf-8 -*-
# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Right-handed laboratory geometry; X beam, Z up, Y outward."""
from __future__ import division

import math

HC_KEV_ANGSTROM = 12.398419843320026


def positive(value, name):
    """Return a finite positive float or raise an informative error."""
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError('{} must be finite and positive.'.format(name))
    return value


def matvec(matrix, vector):
    """Multiply a three-by-three matrix by a Cartesian column vector."""
    return tuple(sum(matrix[i][j]*vector[j] for j in range(3)) for i in range(3))


def detector_geometry(distance, columns, rows, pixel_y_um, pixel_z_um, center_col, center_row, tilt_x=0, tilt_y=0, tilt_z=0):
    """Validate a detector tilted about the fixed direct-beam point (distance,0,0).
    Args:
        distance: Distance along the incident +X beam in mm, not normal distance.
        columns, rows: Positive integer image dimensions.
        pixel_y_um, pixel_z_um: Pixel pitches in micrometres.
        center_col, center_row: Zero-based beam-centre pixel coordinates.
        tilt_x, tilt_y, tilt_z: Degrees, right-handed fixed lab rotations X then Y then Z.
            R=Rz*Ry*Rx. Zero tilt has columns +Y and rows -Z.
    Returns:
        detector: Geometry dictionary; image rows increase downward (-Z).
    """
    if int(columns) != float(columns) or int(rows) != float(rows) or min(int(columns), int(rows)) <= 0:
        raise ValueError('Detector dimensions must be positive integers.')
    cx, cy = float(center_col), float(center_row)
    if not math.isfinite(cx) or not math.isfinite(cy):
        raise ValueError('Beam centre must be finite.')
    angles = tuple(float(value) for value in (tilt_x, tilt_y, tilt_z))
    if not all(math.isfinite(value) for value in angles):
        raise ValueError('Detector tilt angles must be finite.')
    tx, ty, tz = [math.radians(value) for value in angles]
    sx, sy, sz = math.sin(tx), math.sin(ty), math.sin(tz)
    cox, coy, coz = math.cos(tx), math.cos(ty), math.cos(tz)
    rotation = ((coz*coy, coz*sy*sx - sz*cox, coz*sy*cox + sz*sx),
                (sz*coy, sz*sy*sx + coz*cox, sz*sy*cox - coz*sx), (-sy, coy*sx, coy*cox))
    detector = dict(distance=positive(distance, 'Distance'), columns=int(columns), rows=int(rows),
                py=positive(pixel_y_um, 'Horizontal pitch')/1000, pz=positive(pixel_z_um, 'Vertical pitch')/1000,
                cx=cx, cy=cy, tilts=angles, normal=matvec(rotation, (1, 0, 0)),
                col_axis=matvec(rotation, (0, 1, 0)), row_axis=matvec(rotation, (0, 0, -1)))
    if detector['normal'][0] <= 1e-8:
        raise ValueError('Detector must face the sample; reduce the tilt angles.')
    if any(pixel_lab(detector, col, row)[0] <= 0 for col in (-0.5, detector['columns'] - 0.5) for row in (-0.5, detector['rows'] - 0.5)):
        raise ValueError('Tilt places part of the detector behind the sample (X <= 0); increase distance or reduce tilt.')
    return detector


def pixel_lab(detector, col, row):
    """Return a detector pixel's lab coordinates in mm, allowing fractional pixels."""
    u, v = (col - detector['cx'])*detector['py'], (row - detector['cy'])*detector['pz']
    return tuple((detector['distance'] if axis == 0 else 0) + u*detector['col_axis'][axis] + v*detector['row_axis'][axis] for axis in range(3))


def ray_detector_hit(ray, detector, clip=True):
    """Intersect a lab ray from the origin with the detector; optionally clip to its edges."""
    denominator = sum(a*b for a, b in zip(ray, detector['normal']))
    if ray[0] <= 0 or denominator <= 1e-12:
        return None
    length = detector['distance']*detector['normal'][0]/denominator
    delta = tuple(length*ray[axis] - (detector['distance'] if axis == 0 else 0) for axis in range(3))
    col = detector['cx'] + sum(a*b for a, b in zip(delta, detector['col_axis']))/detector['py']
    row = detector['cy'] + sum(a*b for a, b in zip(delta, detector['row_axis']))/detector['pz']
    if not clip or (-0.5 <= col < detector['columns'] - 0.5 and -0.5 <= row < detector['rows'] - 0.5):
        return col, row
    return None


def clip_detector_segment(first, last, detector):
    """Clip a straight ring segment to detector edges before drawing at any zoom."""
    x, y = first
    dx, dy = last[0] - x, last[1] - y
    lower, upper = 0.0, 1.0
    for p, q in ((-dx, x + 0.5), (dx, detector['columns'] - 0.5 - x), (-dy, y + 0.5), (dy, detector['rows'] - 0.5 - y)):
        if p == 0:
            if q < 0:
                return None
            continue
        ratio = q/p
        if p < 0:
            lower = max(lower, ratio)
        else:
            upper = min(upper, ratio)
        if lower > upper:
            return None
    return x + lower*dx, y + lower*dy, x + upper*dx, y + upper*dy


def maximum_two_theta(detector):
    """Return the detector corner's maximum scattering angle in degrees."""
    corners = [pixel_lab(detector, col, row) for col in (-0.5, detector['columns'] - 0.5) for row in (-0.5, detector['rows'] - 0.5)]
    return max(math.degrees(math.atan2(math.hypot(y, z), x)) for x, y, z in corners)


def omega_solutions(g, wavelength, start, end):
    """Solve the elastic Bragg condition for positive Rz rotation.
    Args:
        g: Reciprocal vector at omega=0, in inverse angstroms WITHOUT 2*pi.
        wavelength: Incident wavelength in angstroms.
        start, end: Increasing exposure limits in degrees, half-open [start,end).
    Returns:
        angles: Sorted Bragg angles; positive omega is CCW viewed from +Z.
    """
    gx, gy, gz = g
    radius = math.hypot(gx, gy)
    if radius < 1e-15:
        return []
    ratio = -wavelength*(gx*gx + gy*gy + gz*gz)/(2*radius)
    if abs(ratio) > 1 + 1e-12:
        return []
    phase = math.atan2(-gy, gx)
    delta = math.acos(max(-1, min(1, ratio)))
    result = []
    for root in (math.degrees(phase - delta), math.degrees(phase + delta)):
        first = int(math.floor((start - root)/360)) - 1
        last = int(math.ceil((end - root)/360)) + 1
        for turn in range(first, last + 1):
            angle = root + 360*turn
            if abs(angle - start) < 1e-9:
                angle = start
            if abs(angle - end) < 1e-9:
                angle = end
            if start <= angle < end and not any(abs(angle - old) < 1e-8 for old in result):
                result.append(angle)
    return sorted(result)


def detector_hit(g, omega, wavelength, detector):
    """Return (column,row) for an on-Bragg ray, or None outside the detector."""
    angle = math.radians(omega)
    gx = math.cos(angle)*g[0] - math.sin(angle)*g[1]
    gy = math.sin(angle)*g[0] + math.cos(angle)*g[1]
    dx, dy, dz = 1 + wavelength*gx, wavelength*gy, wavelength*g[2]
    return ray_detector_hit((dx, dy, dz), detector)
