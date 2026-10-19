"""Focused checks of routing, continuous collision tests, and aircraft lifecycle."""
import sys
sys.dont_write_bytecode = True
import json
from pathlib import Path
import unittest
import numpy as np
from network_capacity import Network, demand, simulate, swept_distances


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


if __name__ == '__main__':
    unittest.main(verbosity=2)
