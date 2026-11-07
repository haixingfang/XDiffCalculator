# -*- coding: utf-8 -*-
# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Independent analytical checks for tilted-detector geometry and intensity conventions."""
import math
import random
import unittest

from app import DEFAULTS, calculate, ring_segments
from geometry import detector_geometry, pixel_lab, ray_detector_hit, maximum_two_theta, clip_detector_segment


class DetectorTiltTests(unittest.TestCase):
    """Check physical lab geometry independently of the GUI controls."""

    def test_zero_tilt_and_direct_beam(self):
        """Zero tilt reproduces the original projection and tilts preserve the beam point."""
        detector = detector_geometry(150, 2001, 2101, 75, 100, 1000, 1050)
        ray = (1, 0.1, -0.2)
        col, row = ray_detector_hit(ray, detector)
        self.assertAlmostEqual(col, 1000 + 150*0.1/0.075)
        self.assertAlmostEqual(row, 1050 + 150*0.2/0.1)
        for angles in ((0, 0, 0), (10, -15, 20), (-10, 25, -5)):
            tilted = detector_geometry(150, 2001, 2101, 75, 100, 1000, 1050, *angles)
            self.assertEqual(ray_detector_hit((1, 0, 0), tilted), (1000, 1050))
            self.assertEqual(pixel_lab(tilted, 1000, 1050), (150, 0, 0))

    def test_single_axis_signs(self):
        """Analytical single-axis rotations establish angle signs and image directions."""
        sine, cosine = math.sin(math.radians(20)), math.cos(math.radians(20))
        rolled = detector_geometry(150, 4001, 4001, 100, 100, 2000, 2000, 90, 0, 0)
        col, row = ray_detector_hit((cosine, sine, 0), rolled)
        self.assertAlmostEqual(col, 2000)
        self.assertAlmostEqual(row, 2000 + 150*math.tan(math.radians(20))/0.1)
        tilted_y = detector_geometry(150, 4001, 4001, 100, 100, 2000, 2000, 0, 15, 0)
        col, row = ray_detector_hit((cosine, 0, sine), tilted_y)
        self.assertAlmostEqual(row, 2000 - 150*sine/math.cos(math.radians(35))/0.1)
        tilted_z = detector_geometry(150, 4001, 4001, 100, 100, 2000, 2000, 0, 0, 15)
        col, row = ray_detector_hit((cosine, sine, 0), tilted_z)
        self.assertAlmostEqual(col, 2000 + 150*sine/math.cos(math.radians(5))/0.1)

    def test_rotation_order_and_inverse(self):
        """Independent sequential rotations and sampled pixels round-trip through ray intersection."""
        detector = detector_geometry(150, 2068, 2162, 75, 75, 1033.5, 1080.5, 12, -18, 23)
        vector = (0, 1, 0)
        x, y, z = vector
        a, b, c = [math.radians(value) for value in (12, -18, 23)]
        y, z = math.cos(a)*y - math.sin(a)*z, math.sin(a)*y + math.cos(a)*z
        x, z = math.cos(b)*x + math.sin(b)*z, -math.sin(b)*x + math.cos(b)*z
        x, y = math.cos(c)*x - math.sin(c)*y, math.sin(c)*x + math.cos(c)*y
        for actual, expected in zip(detector['col_axis'], (x, y, z)):
            self.assertAlmostEqual(actual, expected)
        rng = random.Random(5)
        maximum = maximum_two_theta(detector)
        for unused in range(100):
            col, row = rng.uniform(0, 2067), rng.uniform(0, 2161)
            point = pixel_lab(detector, col, row)
            recovered = ray_detector_hit(point, detector)
            self.assertAlmostEqual(recovered[0], col, places=9)
            self.assertAlmostEqual(recovered[1], row, places=9)
            angle = math.degrees(math.atan2(math.hypot(point[1], point[2]), point[0]))
            self.assertLessEqual(angle, maximum + 1e-12)

    def test_invalid_and_clipping(self):
        """Reject backward/invalid detector planes and safely clip long ring segments."""
        for angles in ((0, 90, 0), (0, 180, 0), (float('nan'), 0, 0), (0, 89, 0)):
            with self.assertRaises(ValueError):
                detector_geometry(10, 2001, 2001, 100, 100, 1000, 1000, *angles)
        detector = detector_geometry(150, 100, 100, 75, 75, 50, 50)
        self.assertIsNone(ray_detector_hit((-1, 0, 0), detector))
        self.assertIsNone(ray_detector_hit((0, 1, 0), detector))
        self.assertEqual(clip_detector_segment((-100, 50), (200, 50), detector), (-0.5, 50, 99.5, 50))
        self.assertIsNone(clip_detector_segment((-10, -10), (-10, 50), detector))

    def test_tilted_rings_and_lab_intensity(self):
        """Tilted ring vertices satisfy the diffraction cone; tilts do not alter physical spot weights."""
        result = calculate(dict(DEFAULTS, tilt_x='12', tilt_y='-18', tilt_z='23'))
        line = result['lines'][0]
        segments = list(ring_segments(line, result['detector'], result['polarization'], result['intensity_max']))
        self.assertTrue(segments)
        for segment in segments:
            for col, row in ((segment[0], segment[1]), (segment[2], segment[3])):
                point = pixel_lab(result['detector'], col, row)
                self.assertAlmostEqual(math.degrees(math.atan2(math.hypot(point[1], point[2]), point[0])), line['two_theta'], places=8)
        settings = dict(DEFAULTS, mode='Random grains', number='8', end='180')
        original = calculate(settings)
        tilted = calculate(dict(settings, tilt_x='12', tilt_y='-18', tilt_z='23'))
        key = lambda event: (event['grain'], event['hkl'], round(event['omega'], 8))
        lookup = dict((key(event), event) for event in original['events'])
        shared = [event for event in tilted['events'] if key(event) in lookup]
        self.assertGreater(len(shared), 100)
        self.assertTrue(any(abs(event['col'] - lookup[key(event)]['col']) > 1 for event in shared))
        for event in shared:
            reference = lookup[key(event)]
            self.assertAlmostEqual(event['intensity'], reference['intensity'], places=8)
            self.assertAlmostEqual(event['azimuth'], reference['azimuth'], places=12)


if __name__ == '__main__':
    unittest.main()
