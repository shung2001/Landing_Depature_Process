"""Verify moving fleets, failed overpacking, OD laps, and Kepler coordinates."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd

from network_capacity import Network
from circulation_capacity import RingTrack, preload, feed, pair_minimum, kepler_csv
from export_kepler_flows import trips


class CirculationChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = json.loads(Path(__file__).with_name('capacity_config.json').read_text(encoding='utf-8'))
        cls.net = Network(cls.cfg)
        cls.tracks = [RingTrack(cls.net, ri) for ri in (0, 1)]

    def test_original_closed_ring_and_direction(self):
        for ri, track in enumerate(self.tracks):
            np.testing.assert_allclose(track.at(np.array([0.])), track.at(np.array([track.length])))
            xy = track.xyz[:, :2] - track.xyz[0, :2]
            area = np.sum(xy[:-1, 0] * xy[1:, 1] - xy[1:, 0] * xy[:-1, 1])
            self.assertGreater(area if ri == 0 else -area, 0)

    def test_fourteen_aircraft_complete_multiple_real_laps(self):
        for track in self.tracks:
            result = preload(track, 100, 14, 'uniform', 42)
            self.assertTrue(result['valid'])
            self.assertGreater(result['laps_per_aircraft'], 2)
            self.assertGreaterEqual(result['min_swept_m'], 100)

    def test_overpacked_ring_rejected(self):
        for track in self.tracks:
            result = preload(track, 1000, track.geometric_max(1000) + 1, 'uniform', 42)
            self.assertFalse(result['valid'])

    def test_continuous_demand_reaches_fourteen_without_stopping(self):
        cfg = dict(self.cfg, duration_s=1200)
        for track in self.tracks:
            result, history, records, states = feed(track, self.net, cfg, 100, 14, 360, 42, 'weighted_entries', True)
            self.assertEqual(result['admitted'], 14)
            self.assertGreater(result['seconds_at_least_10_moving'], 600)
            self.assertGreater(result['end_queue'], 0)
            frame = pd.DataFrame(records).drop_duplicates(['time_s', 'aircraft_id'])
            first = frame[frame.aircraft_id == frame.aircraft_id.iloc[0]].sort_values('time_s')
            speed = np.diff(first.distance_m) / np.diff(first.time_s)
            np.testing.assert_allclose(speed, track.speed)

    def test_swept_pair_minimum_detects_midstep_crossing(self):
        a = np.array([[-100., 0, 0], [0, -100., 0]])
        b = np.array([[100., 0, 0], [0, 100., 0]])
        self.assertAlmostEqual(pair_minimum(a, b), 0)

    def test_extra_laps_stay_on_directed_network(self):
        expanded = Network(dict(self.cfg, arrival_extra_ring_laps=2))
        longer = 0
        for key, original in self.net.paths.items():
            route = expanded.paths[key]
            for a, b in zip(route['path'], route['path'][1:]):
                self.assertEqual(expanded.edges[a].v, expanded.edges[b].u)
            if key[0] == 'arrival' and any(self.net.edges[e].ring >= 0 for e in original['path']):
                self.assertGreater(len(route['xyz']), len(original['xyz']))
                longer += 1
        self.assertGreater(longer, 0)

    def test_kepler_coordinates_timestamp_and_separate_phases(self):
        point = self.tracks[0].xyz[0]
        rows = [dict(time_s=t, aircraft_id=1, x=point[0], y=point[1], z=point[2], ring=0,
                     phase='arrival' if t < 20 else 'departure') for t in (0, 10, 20, 30)]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            frame = kepler_csv(rows, path / 'flow.csv', 'test', 42)
            self.assertTrue(frame.longitude.between(124, 132).all())
            self.assertTrue(frame.latitude.between(33, 39).all())
            np.testing.assert_array_equal(np.diff(frame.timestamp_ms), [10000, 10000, 10000])
            self.assertEqual(trips(frame, path / 'trips.geojson'), 2)
            data = json.loads((path / 'trips.geojson').read_text(encoding='utf-8'))
            self.assertEqual(len(data['features'][0]['geometry']['coordinates'][0]), 4)


if __name__ == '__main__':
    unittest.main(verbosity=2)
