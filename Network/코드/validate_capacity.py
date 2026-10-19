"""Focused checks of routing, continuous collision tests, and aircraft lifecycle."""
import sys
sys.dont_write_bytecode = True
import json
import copy
from pathlib import Path
import unittest
import numpy as np
from network_capacity import (Network, demand, simulate, swept_distances, validate_config,
                              capacity_scenarios, capacity_expansion_history)


class CapacityChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads(Path(__file__).with_name('capacity_config.json').read_text(encoding='utf-8'))
        cls.net = Network(cls.config)

    def test_continuous_crossing_not_only_endpoints(self):
        # Both endpoints are far apart; the paths intersect at mid-step.
        d = swept_distances(np.array([-100., 0, 0]), np.array([100., 0, 0]),
                            np.array([[0., -100, 0]]), np.array([[0., 100, 0]]))
        self.assertAlmostEqual(d[0], 0.)

    def test_crossing_at_different_altitudes(self):
        d = swept_distances(np.array([-100., 0, 0]), np.array([100., 0, 0]),
                            np.array([[0., -100, 152.4]]), np.array([[0., 100, 152.4]]))
        self.assertAlmostEqual(d[0], 152.4)

    def test_all_od_reachable(self):
        self.assertEqual(len(self.net.paths), 64)
        for (kind, name, fid), route in self.net.paths.items():
            start = self.net.terminals[name] if kind == 'arrival' else self.net.ports[fid]
            end = self.net.ports[fid] if kind == 'arrival' else self.net.terminals[name]
            np.testing.assert_array_equal(route['xyz'][0], self.net.nodes[start])
            np.testing.assert_array_equal(route['xyz'][-1], self.net.nodes[end])
            for a, b in zip(route['path'], route['path'][1:]):
                self.assertEqual(self.net.edges[a].v, self.net.edges[b].u)

    def test_ring_directions_and_exact_endpoints(self):
        for ring in (0, 1):
            edges = [e for e in self.net.edges if e.ring == ring]
            area = sum(np.sum(e.xyz[:-1, 0] * e.xyz[1:, 1] - e.xyz[1:, 0] * e.xyz[:-1, 1]) for e in edges)
            self.assertGreater(area if ring == 0 else -area, 0)
        for edge in self.net.edges:
            np.testing.assert_array_equal(edge.xyz[0], self.net.nodes[edge.u])
            np.testing.assert_array_equal(edge.xyz[-1], self.net.nodes[edge.v])

    def test_demand_split_and_reproducibility(self):
        config = dict(self.config, duration_s=3600)
        requests = demand(self.net, config, 100, 42)
        self.assertEqual(requests, demand(self.net, config, 100, 42))
        counts = {name: sum(r['origin'] == name for r in requests) for name in self.net.terminals}
        self.assertEqual(counts['강원'], 10)
        self.assertTrue(all(count == 30 for name, count in counts.items() if name != '강원'))

    def test_landing_before_takeoff_and_pad_capacity(self):
        config = dict(self.config, duration_s=2400, warmup_s=300)
        requests = demand(self.net, config, 30, 42)
        metrics, history, states, tracks, events = simulate(self.net, config, requests, (14, 14),
                                                           100, 'validation', 42, detail=True)
        self.assertEqual(metrics['buffer_violations'], 0)
        self.assertGreaterEqual(metrics['min_swept_separation_m'], 100 - 1e-6)
        self.assertTrue((history.outer <= 14).all())
        self.assertTrue((history.inner <= 14).all())
        self.assertTrue((history.holding <= 8).all())
        occupants = {}
        for event in events:
            fid, aid = event['fid'], event['aircraft_id']
            if event['event'] == 'landed':
                self.assertNotIn(fid, occupants)
                occupants[fid] = aid
            elif event['event'] == 'takeoff':
                self.assertEqual(occupants.pop(fid), aid)
        for state in states:
            if np.isfinite(state['departure_s']):
                self.assertTrue(np.isfinite(state['landed_s']))
                self.assertGreaterEqual(state['departure_s'] - state['landed_s'], config['turnaround_s'])

    def test_replay_preserves_integer_demand(self):
        config = dict(self.config, duration_s=1200)
        original = demand(self.net, config, 18, 42)
        # 6 * 3600 / 1200 is exact, while this nextafter value deliberately
        # represents the downward floating-point error of a stored realised rate.
        realised = np.nextafter(len(original) * 3600 / config['duration_s'], 0)
        replay = demand(self.net, config, realised, 42, count=len(original))
        self.assertEqual(original, replay)

    def test_invalid_measurement_and_separation_inputs(self):
        for change in (dict(step_s=0), dict(warmup_s=self.config['duration_s']),
                       dict(separations_m=[50]), dict(separations_m=[1001]),
                       dict(capacity_checkpoints=[28, 14])):
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate_config(dict(self.config, **change))
        cfg = copy.deepcopy(self.config)
        cfg['speed_defaults_kts']['접근'] = float('nan')
        with self.assertRaises(ValueError):
            validate_config(cfg)

    def test_incremental_expansion_and_geometric_spacing(self):
        history = capacity_expansion_history(self.net, self.config)
        for sep in self.config['separations_m']:
            rows = [r for r in history if r['separation_m'] == sep]
            for a, b in zip(rows, rows[1:]):
                self.assertEqual((b['outer'] - a['outer']) + (b['inner'] - a['inner']), 1)
        for row in capacity_scenarios(self.net, self.config):
            self.assertGreaterEqual(row['outer_uniform_chord_m'], row['separation_m'])
            self.assertGreaterEqual(row['inner_uniform_chord_m'], row['separation_m'])

    def test_pad_reserved_before_final_zone_prevents_circular_wait(self):
        metrics, _, states, _, _ = simulate(self.net, self.config,
            demand(self.net, self.config, 15, 20261019), (14, 14), 100, 'pad_regression', 20261019)
        # This sequence previously stopped after ten landings: an inbound
        # aircraft owned the departure's zone while waiting for its occupied pad.
        self.assertGreater(sum(np.isfinite(s['landed_s']) for s in states), 10)
        self.assertEqual(metrics['buffer_violations'], 0)
        self.assertGreaterEqual(metrics['min_ring_swept_separation_m'], 100 - 1e-6)


if __name__ == '__main__':
    unittest.main(verbosity=2)
