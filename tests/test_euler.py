# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Bunge convention, singular-angle round trips and user input checks."""
import math
import os
import tempfile
import unittest
from orientations import euler_to_matrix, matrix_to_euler, load_orientations, save_orientations, random_orientations
from app import DEFAULTS, single_orientation, settings_from_disk


class EulerTests(unittest.TestCase):
    def assert_matrix_close(self, first, second, tolerance=1e-10):
        self.assertLess(max(abs(a - b) for row1, row2 in zip(first, second) for a, b in zip(row1, row2)), tolerance)

    def test_known_bunge_rotations(self):
        self.assert_matrix_close(euler_to_matrix([90, 0, 0]), [[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        self.assert_matrix_close(euler_to_matrix([0, 90, 0]), [[1, 0, 0], [0, 0, -1], [0, 1, 0]])
        # Independent composition of fixed-axis elementary matrices.
        def multiply(a, b):
            return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        def rz(deg):
            c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
            return [[c, -s, 0], [s, c, 0], [0, 0, 1]]
        c, s = math.cos(math.radians(42)), math.sin(math.radians(42))
        expected = multiply(multiply(rz(17), [[1, 0, 0], [0, c, -s], [0, s, c]]), rz(73))
        self.assert_matrix_close(euler_to_matrix([17, 42, 73]), expected)

    def test_random_and_singular_round_trips(self):
        matrices = random_orientations(1000, 28)
        for angles in ([30, 0, 80], [30, 180, 80], [10, 0.0000001, 20], [40, 179.9999999, 50]):
            matrices.append(euler_to_matrix(angles))
        for matrix in matrices:
            angles = matrix_to_euler(matrix)
            self.assertTrue(0 <= angles[0] < 360 and 0 <= angles[1] <= 180 and 0 <= angles[2] < 360)
            self.assert_matrix_close(euler_to_matrix(angles), matrix)
        for Phi in (0, 180):
            self.assertEqual(matrix_to_euler(euler_to_matrix([30, Phi, 80]))[2], 0)

    def test_default_matches_previous_ipf_z_example(self):
        matrix = single_orientation(DEFAULTS['single_euler'])
        expected = [[0.99681528, 0, -0.07974522], [-0.00396964, 0.99876026, -0.04962044], [0.07964636, 0.04977897, 0.99557949]]
        self.assert_matrix_close(matrix, expected, 2e-8)
        self.assertAlmostEqual(math.degrees(math.acos(matrix[2][2])), 5.389322, places=6)
        self.assertGreater(matrix[2][0], matrix[2][1])
        self.assertGreater(matrix[2][1], 0)
        self.assertEqual((DEFAULTS['energy_min'], DEFAULTS['energy_max']), ('5', '120'))

    def test_legacy_single_matrix_settings(self):
        settings = settings_from_disk(dict(single_u='0,-1,0,1,0,0,0,0,1'))
        self.assertNotIn('single_u', settings)
        self.assert_matrix_close(single_orientation(settings['single_euler']), [[0, -1, 0], [1, 0, 0], [0, 0, 1]])

    def test_euler_file_and_errors(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, 'euler.csv')
            matrices = random_orientations(10, 12)
            save_orientations(path, matrices)
            with open(path) as stream:
                self.assertEqual(stream.readline().strip(), 'phi1_deg,Phi_deg,phi2_deg')
            for actual, expected in zip(load_orientations(path), matrices):
                self.assert_matrix_close(actual, expected)
            for text in ('# degrees\nphi1 Phi phi2\n0 0 0\n90 0 0', '0,0,0\n90,0,0'):
                with open(path, 'w') as stream:
                    stream.write(text)
                self.assertEqual(len(load_orientations(path)), 2)
            for text in ('1,0,0,0,1,0,0,0,1', '0,nan,0', '0,-1,0', '0,181,0', '0,0', ''):
                with open(path, 'w') as stream:
                    stream.write(text)
                with self.assertRaises(ValueError):
                    load_orientations(path)


if __name__ == '__main__':
    unittest.main()
