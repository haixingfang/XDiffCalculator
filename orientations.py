# -*- coding: utf-8 -*-
# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Uniform rotations and explicit crystal-to-sample matrix import."""
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


def load_orientations(path):
    """Read one row-major U matrix per CSV/whitespace line; optional U11..U33 header."""
    matrices = []
    header = ['U11', 'U12', 'U13', 'U21', 'U22', 'U23', 'U31', 'U32', 'U33']
    with open(path, 'r', encoding='utf-8-sig') as stream:
        for line_number, line in enumerate(stream, 1):
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            fields = line.replace(',', ' ').split()
            if fields == header and not matrices:
                continue
            try:
                values = [float(value) for value in fields]
                if len(values) != 9:
                    raise ValueError('Expected nine matrix entries.')
                matrices.append(validate_rotation([values[0:3], values[3:6], values[6:9]]))
            except ValueError as error:
                raise ValueError('Orientation line {}: {}'.format(line_number, error))
            if len(matrices) > 100000:
                raise ValueError('Orientation file exceeds 100000 grains.')
    if not matrices:
        raise ValueError('No orientation matrices found.')
    return matrices


def save_orientations(path, matrices):
    """Export row-major crystal-to-sample matrices with an explicit CSV header."""
    with open(path, 'w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['U11', 'U12', 'U13', 'U21', 'U22', 'U23', 'U31', 'U32', 'U33'])
        for matrix in matrices:
            writer.writerow(['{:.15g}'.format(value) for row in validate_rotation(matrix) for value in row])
