# -*- coding: utf-8 -*-
# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Uniform rotations and Bunge ZXZ Euler-angle input/output in degrees."""
from __future__ import division

import csv
import math
import random


def validate_rotation(matrix, tolerance=1e-5):
    """Validate a proper rotation without silently orthogonalizing it."""
    if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
        raise ValueError('An orientation must be a 3 by 3 matrix.')
    u = [[float(value) for value in row] for row in matrix]
    if not all(math.isfinite(value) for row in u for value in row):
        raise ValueError('Orientation entries must be finite.')
    error = max(abs(sum(u[k][i]*u[k][j] for k in range(3)) - int(i == j)) for i in range(3) for j in range(3))
    det = (u[0][0]*(u[1][1]*u[2][2] - u[1][2]*u[2][1]) - u[0][1]*(u[1][0]*u[2][2] - u[1][2]*u[2][0])
           + u[0][2]*(u[1][0]*u[2][1] - u[1][1]*u[2][0]))
    if error > tolerance or abs(det - 1) > tolerance:
        raise ValueError('Orientation must satisfy U.T U = I and det(U) = +1; error={:.3g}, det={:.6g}.'.format(error, det))
    return u


def random_orientations(number, seed=12345):
    """Generate reproducible Haar-uniform crystal-to-sample rotation matrices.
    Args:
        number: Integer grain count, from 1 to 100000.
        seed: Seed for a local random generator; global state is unchanged.
    Returns:
        matrices: List of 3 by 3 matrices, generated from uniform unit quaternions.
    """
    if int(number) != float(number) or not 1 <= int(number) <= 100000:
        raise ValueError('Number of grains must be an integer in [1, 100000].')
    rng = random.Random(seed)
    matrices = []
    for unused in range(int(number)):
        q = [rng.gauss(0, 1) for k in range(4)]
        norm = math.sqrt(sum(value*value for value in q))
        if norm == 0:
            raise RuntimeError('Zero random quaternion; choose another seed.')
        w, x, y, z = [value/norm for value in q]
        matrices.append([[1 - 2*(y*y + z*z), 2*(x*y - z*w), 2*(x*z + y*w)],
                         [2*(x*y + z*w), 1 - 2*(x*x + z*z), 2*(y*z - x*w)],
                         [2*(x*z - y*w), 2*(y*z + x*w), 1 - 2*(x*x + y*y)]])
    return matrices


def euler_to_matrix(angles):
    """Convert Bunge Euler angles to the crystal-to-sample matrix U.
    Args:
        angles: (phi1, Phi, phi2) in degrees, with Phi in [0,180].
    Returns:
        matrix: U = Rz(phi1) Rx(Phi) Rz(phi2), the transpose of the
            conventional passive Bunge sample-to-crystal matrix.
    """
    if len(angles) != 3:
        raise ValueError('Expected three Bunge Euler angles: phi1, Phi, phi2 in degrees.')
    phi1, Phi, phi2 = [float(value) for value in angles]
    if not all(math.isfinite(value) for value in (phi1, Phi, phi2)) or not 0 <= Phi <= 180:
        raise ValueError('Euler angles must be finite and Phi must be in [0,180] degrees.')
    a, b, c = [math.radians(value) for value in (phi1 % 360, Phi, phi2 % 360)]
    ca, sa, cb, sb, cc, sc = math.cos(a), math.sin(a), math.cos(b), math.sin(b), math.cos(c), math.sin(c)
    return [[ca*cc - sa*cb*sc, -ca*sc - sa*cb*cc, sa*sb],
            [sa*cc + ca*cb*sc, -sa*sc + ca*cb*cc, -ca*sb],
            [sb*sc, sb*cc, cb]]


def matrix_to_euler(matrix):
    """Convert crystal-to-sample U to canonical Bunge ZXZ angles in degrees.
    Args:
        matrix: Proper 3 by 3 rotation matrix.
    Returns:
        angles: phi1/phi2 in [0,360), Phi in [0,180]. At Phi=0 or 180,
            set phi2=0 and retain the equivalent combined rotation in phi1.
            No crystal-symmetry reduction is applied.
    """
    u = validate_rotation(matrix)
    sine = math.hypot(u[2][0], u[2][1])
    Phi = math.atan2(sine, u[2][2])
    if sine < 1e-12:
        phi1, phi2 = math.atan2(u[1][0], u[0][0]), 0.0
    else:
        phi1, phi2 = math.atan2(u[0][2], -u[1][2]), math.atan2(u[2][0], u[2][1])
    return [math.degrees(phi1) % 360, math.degrees(Phi), math.degrees(phi2) % 360]


def load_orientations(path):
    """Read Bunge ZXZ Euler angles (degrees), one phi1/Phi/phi2 triple per row.
    Args:
        path: CSV or whitespace TXT, with optional phi1_deg,Phi_deg,phi2_deg
            or phi1,Phi,phi2 header. Blank lines and # comments are allowed.
    Returns:
        matrices: Crystal-to-sample rotations for the existing diffraction core.
    """
    matrices = []
    headers = [['phi1_deg', 'Phi_deg', 'phi2_deg'], ['phi1', 'Phi', 'phi2']]
    with open(path, 'r', encoding='utf-8-sig') as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            fields = line.replace(',', ' ').split()
            if fields in headers and not matrices:
                continue
            if len(fields) == 9:
                raise ValueError('Orientation line {}: matrix files are no longer accepted. Supply phi1_deg,Phi_deg,phi2_deg (Bunge ZXZ).'.format(line_number))
            try:
                values = [float(value) for value in fields]
                matrices.append(euler_to_matrix(values))
            except ValueError as error:
                raise ValueError('Orientation line {}: {}'.format(line_number, error))
            if len(matrices) > 100000:
                raise ValueError('Orientation file exceeds 100000 grains.')
    if not matrices:
        raise ValueError('No Euler-angle orientations found.')
    return matrices


def save_orientations(path, matrices):
    """Export calculated orientations as Bunge phi1/Phi/phi2 angles in degrees."""
    with open(path, 'w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['phi1_deg', 'Phi_deg', 'phi2_deg'])
        for matrix in matrices:
            writer.writerow(['{:.15g}'.format(value) for value in matrix_to_euler(matrix)])
