"""FAF→IAF 교차고도 보간, 방향, 기체 버퍼와 3D 거리 회귀 검사."""
import importlib.util
from pathlib import Path
import unittest

import pandas as pd
import numpy as np
from shapely.geometry import LineString

spec = importlib.util.spec_from_file_location(
    'faf_height_gap', Path(__file__).with_name('FAF_외곽이동시_1500ft와의_간격.py'))
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class HeightGapTests(unittest.TestCase):
    def setUp(self):
        self.ring = LineString([(5, -5), (15, -5), (15, 5), (5, 5), (5, -5)])

    def route(self, points):
        row = pd.Series({'geometry': LineString(points)})
        return analysis.analyze_route(row, 'test', 1, self.ring, 457.2, .01)

    def test_crossing_height_and_buffer_clearance(self):
        route = self.route([(0, 0, 400), (20, 0, 600)])
        a, b = route['contacts']
        self.assertAlmostEqual(a['route_z_m'], 450)
        self.assertAlmostEqual(a['signed_gap_m'], -7.2)
        self.assertAlmostEqual(a['vertical_clearance_m'], 3.44)
        self.assertAlmostEqual(b['route_z_m'], 550)
        self.assertEqual(a['crossing_role'], 'interior_crossing')
        self.assertFalse(a['vertical_buffer_overlap'])

    def test_reversed_storage_has_same_faf_progress_and_heights(self):
        forward = self.route([(0, 0, 400), (20, 0, 600)])
        reverse = self.route([(20, 0, 600), (0, 0, 400)])
        self.assertEqual(forward['contacts'], reverse['contacts'])

    def test_vertical_overlap_and_exact_3d_distance(self):
        route = self.route([(0, 0, 455), (20, 0, 455)])
        self.assertTrue(route['contacts'][0]['vertical_buffer_overlap'])
        ring_xyz = np.column_stack((np.asarray(self.ring.coords), np.full(5, 457.2)))
        analysis.analyze_buffer(route, ring_xyz, 3.76)
        self.assertAlmostEqual(route['buffer_result']['min_distance_3d_m'], 2.2)
        self.assertTrue(route['buffer_result']['buffer_overlap'])

    def test_endpoint_is_distinct_from_mid_route_crossing(self):
        route = self.route([(5, 0, 400), (20, 0, 600)])
        self.assertEqual(route['contacts'][0]['crossing_role'], 'endpoint_contact')
        self.assertEqual(route['contacts'][1]['crossing_role'], 'interior_crossing')


if __name__ == '__main__':
    unittest.main()
