"""Directed 3D network and reproducible, conservative capacity experiment.

python Network/코드/network_capacity.py
python Network/코드/network_capacity.py --quick
Original SHP files are read only. All generated data go to Network/결과물.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
from dataclasses import dataclass
import hashlib
import heapq
import json
import math
import os
from pathlib import Path
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import get_coordinates
from shapely.geometry import LineString, Point
from shapely.ops import substring
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / '자료/최종정리'
OUTPUT = ROOT / 'Network/결과물'
KT = 1852 / 3600
RINGS = ('외곽_2800', '내곽_2000')
LINES = ('접근', 'Expert_Off', 'Interchange_In', 'Interchange_out',
         'IAF_FAF_In_Update', 'IAF_FAF_Out_Update', 'Short_Final')
POINTS = ('접근_start', '각_건물의_최종_버티포트', 'FAF',
          '외곽_IAF_Interchange', '내곽_Interchange', '내곽_IAF', 'Holding_point')


def csv(rows, path):
    pd.DataFrame(rows).to_csv(path, index=False, encoding='utf-8-sig')


def json_text(value):
    """Use JSON null for unobserved statistics, including numpy scalars."""
    def clean(item):
        if isinstance(item, dict):
            return {k: clean(v) for k, v in item.items()}
        if isinstance(item, (list, tuple)):
            return [clean(v) for v in item]
        if isinstance(item, np.generic):
            item = item.item()
        if isinstance(item, float) and not math.isfinite(item):
            return None
        return item
    return json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False)


def length(xyz):
    return np.r_[0., np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))]


def validate_config(cfg):
    """Reject inputs that silently invalidate spacing or measurement windows."""
    for key in ('duration_s', 'step_s', 'aircraft_buffer_m', 'destination_distance_scale_m'):
        if not math.isfinite(cfg[key]) or cfg[key] <= 0:
            raise ValueError(f'{key} must be finite and positive')
    if not 0 <= cfg['warmup_s'] < cfg['duration_s']:
        raise ValueError('warmup_s must be within the simulation horizon')
    for key in ('hover_s', 'turnaround_s'):
        if not math.isfinite(cfg[key]) or cfg[key] < 0:
            raise ValueError(f'{key} must be finite and nonnegative')
    if cfg['replicates'] < 1 or int(cfg['replicates']) != cfg['replicates']:
        raise ValueError('replicates must be a positive integer')
    for key in ('arrival_rates_per_hour', 'capacity_checkpoints', 'separations_m'):
        if not cfg[key] or any(not math.isfinite(v) or v <= 0 for v in cfg[key]):
            raise ValueError(f'{key} must contain positive finite values')
    if any(not 2 * cfg['aircraft_buffer_m'] <= s <= 1000 for s in cfg['separations_m']):
        raise ValueError('separations_m must be between twice aircraft_buffer_m and 1000 m')
    if cfg['capacity_checkpoints'] != sorted(set(cfg['capacity_checkpoints'])) or any(
            int(v) != v or v < 14 for v in cfg['capacity_checkpoints']):
        raise ValueError('capacity_checkpoints must be increasing integers >= 14')
    for layer in (*LINES, *RINGS):
        speed = cfg['speed_defaults_kts'].get(layer, math.nan)
        if not math.isfinite(speed) or speed <= 0:
            raise ValueError(f'{layer}: finite positive default speed required')


@dataclass
class Edge:
    id: int
    u: int
    v: int
    xyz: np.ndarray
    layer: str
    row: int
    speed: float
    speed_source: str
    ring: int = -1

    @property
    def seconds(self):
        return float(length(self.xyz)[-1] / (self.speed * KT))


class Network:
    def __init__(self, config, source=SOURCE):
        validate_config(config)
        self.config = config
        self.source = source
        self.frames, self.audit, self.nodes, self.edges = {}, [], [], []
        self.adj = defaultdict(list)
        self.terminals, self.ports, self.holdings = {}, {}, {}
        for name in (*POINTS, *LINES, *RINGS):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                frame = gpd.read_file(source / f'{name}.shp')
            if frame.crs is None:
                raise ValueError(f'{name}: CRS missing')
            self.frames[name] = frame.to_crs(5186)
            self.audit.append({'layer': name, 'features': len(frame), 'original_crs': str(frame.crs),
                               'fields': '|'.join(frame.columns), 'geometry': '|'.join(frame.geom_type.unique())})
        for name in POINTS:
            for row, record in self.frames[name].iterrows():
                xyz = get_coordinates(record.geometry, include_z=True)
                if len(xyz) != 1:
                    raise ValueError(f'{name}/{row}: one point per feature required')
                point = xyz[0].copy()
                if not np.isfinite(point[2]):
                    point[2] = float(record['Altitude'])
                node = self.node(point)
                if name == '접근_start':
                    self.terminals[str(record['Name'])] = node
                elif name == '각_건물의_최종_버티포트':
                    self.ports[int(record['fid'])] = node
                elif name == 'Holding_point':
                    self.holdings[int(record['id'])] = dict(node=node, capacity=int(record['Capacity']),
                        radius=float(record['Radius']), straight=float(record['straight']),
                        speed=float(record['Speed_kts']), type=record['Type'])
        for name in LINES:
            for row, record in self.frames[name].iterrows():
                xyz = get_coordinates(record.geometry, include_z=True).copy()
                if record.geometry.geom_type != 'LineString' or len(xyz) < 2:
                    raise ValueError(f'{name}/{row}: LineString required')
                if not np.isfinite(xyz[:, 2]).all():
                    if name != 'Expert_Off':
                        raise ValueError(f'{name}/{row}: missing Z')
                    # Expert_Off has only XY; both endpoint levels are 609.6 m.
                    z0 = self.nearest_xy(xyz[0])[2]
                    z1 = self.nearest_xy(xyz[-1])[2]
                    d = length(xyz[:, :2])
                    xyz[:, 2] = np.interp(d, [0, d[-1]], [z0, z1])
                if name == '접근':
                    if self.closest_node(xyz[0]) not in self.terminals.values():
                        xyz = xyz[::-1].copy()
                elif name == 'Expert_Off':
                    if self.closest_node(xyz[0]) in self.terminals.values():
                        xyz = xyz[::-1].copy()
                elif name in ('Interchange_In', 'IAF_FAF_In_Update', 'Short_Final'):
                    if xyz[0, 2] < xyz[-1, 2]:
                        xyz = xyz[::-1].copy()
                else:
                    if xyz[0, 2] > xyz[-1, 2]:
                        xyz = xyz[::-1].copy()
                self.add_edge(xyz, name, row, record)
                if name == 'Short_Final':
                    self.add_edge(xyz[::-1].copy(), name, row, record)
        # Rings are split only at real endpoint/point nodes, never at XY crossings.
        self.ring_lengths = []
        self.ring_geometry = []
        for ri, name in enumerate(RINGS):
            record = self.frames[name].iloc[0]
            xyz = get_coordinates(record.geometry, include_z=True)
            area = np.sum(xyz[:-1, 0] * xyz[1:, 1] - xyz[1:, 0] * xyz[:-1, 1])
            if (ri == 0 and area < 0) or (ri == 1 and area > 0):
                xyz = xyz[::-1].copy()
            curve = LineString(xyz[:, :2])
            self.ring_lengths.append(curve.length)
            center = (xyz[:, :2].min(axis=0) + xyz[:, :2].max(axis=0)) / 2
            self.ring_geometry.append((center, curve.length / (2 * math.pi), xyz[0, 2]))
            stations = []
            for node, point in enumerate(self.nodes):
                if abs(point[2] - xyz[0, 2]) < .1 and curve.distance(Point(point[:2])) < .1:
                    stations.append((curve.project(Point(point[:2])), node))
            stations.sort()
            if len(stations) < 2:
                raise ValueError(f'{name}: too few ring nodes')
            for i, (a, u) in enumerate(stations):
                b, v = stations[(i + 1) % len(stations)]
                if b > a:
                    xy = np.asarray(substring(curve, a, b).coords)
                else:
                    xy = np.vstack([np.asarray(substring(curve, a, curve.length).coords),
                                    np.asarray(substring(curve, 0, b).coords)[1:]])
                arc = np.column_stack([xy, np.full(len(xy), xyz[0, 2])])
                arc[0], arc[-1] = self.nodes[u], self.nodes[v]
                self.add_edge(arc, name, i, record, ring=ri)
        # One racetrack queue is assigned to each external entry, by spatial proximity.
        self.holding_for = {}
        for name, node in self.terminals.items():
            hid = min(self.holdings, key=lambda h: np.linalg.norm(self.nodes[node] - self.nodes[self.holdings[h]['node']]))
            self.holding_for[name] = hid
        if len(set(self.holding_for.values())) != len(self.terminals):
            raise ValueError('Holding-to-terminal assignment is not one-to-one')
        # Conservative, exclusive crossing-corridor reservations. Spatial tests
        # alone permit head-on deadlock; reservations happen before entering.
        self.edge_resources = defaultdict(set)
        self.resource_edges = {}
        self.edge_trees = {}
        self.conflicts = []
        samples = {}
        for edge in self.edges:
            d = length(edge.xyz)
            s = np.linspace(0, d[-1], max(2, math.ceil(d[-1] / 20) + 1))
            samples[edge.id] = np.column_stack([np.interp(s, d, edge.xyz[:, axis]) for axis in range(3)])
            self.edge_trees[edge.id] = cKDTree(samples[edge.id])
        for eid, sample in samples.items():
            tree = cKDTree(sample)
            for other, other_sample in samples.items():
                if other <= eid:
                    continue
                if self.edges[eid].ring >= 0 and self.edges[other].ring >= 0:
                    continue
                closest = float(tree.query(other_sample)[0].min())
                if closest <= 2 * config['aircraft_buffer_m'] + 20:
                    resource = len(self.conflicts)
                    self.conflicts.append(dict(resource_id=resource, edge_a=eid, edge_b=other,
                                               sampled_min_distance_m=closest))
                    self.edge_resources[eid].add(resource)
                    self.edge_resources[other].add(resource)
                    self.resource_edges[resource] = (eid, other)
        self.paths = {}
        self.route_rows = []
        for name, terminal in self.terminals.items():
            for fid, port in self.ports.items():
                for kind, start, end in [('arrival', terminal, port), ('departure', port, terminal)]:
                    path = self.shortest(start, end)
                    if kind == 'arrival':
                        extra = config.get('arrival_extra_ring_laps', 0)
                        if int(extra) != extra or extra < 0:
                            raise ValueError('arrival_extra_ring_laps must be a nonnegative integer')
                        if extra:
                            expanded, visited = [], set()
                            for eid in path:
                                edge = self.edges[eid]
                                if edge.ring >= 0 and edge.ring not in visited:
                                    visited.add(edge.ring)
                                    cycle, node = [], edge.u
                                    while True:
                                        ring_edge = next(e for e in self.adj[node]
                                                         if self.edges[e].ring == edge.ring)
                                        cycle.append(ring_edge)
                                        node = self.edges[ring_edge].v
                                        if node == edge.u:
                                            break
                                        if len(cycle) > len(self.edges):
                                            raise ValueError('Ring is not a closed directed cycle')
                                    expanded.extend(cycle * int(extra))
                                expanded.append(eid)
                            path = expanded
                    key = (kind, name, fid)
                    self.paths[key] = self.sample(path, hover=kind == 'arrival')
                    self.route_rows.append(dict(kind=kind, terminal=name, fid=fid,
                        distance_m=sum(length(self.edges[e].xyz)[-1] for e in path),
                        freeflow_s=sum(self.edges[e].seconds for e in path) + (config['hover_s'] if kind == 'arrival' else 0),
                        edge_ids='|'.join(map(str, path)), layers='|'.join(self.edges[e].layer for e in path)))
        # Overlapping conflict neighbourhoods must be acquired as one zone;
        # otherwise two aircraft can each own half of the same intersection.
        parent = list(range(len(self.conflicts)))
        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for route in self.paths.values():
            for claim in route['claims']:
                if claim:
                    first = next(iter(claim))
                    for resource in claim:
                        parent[find(resource)] = find(first)
        for row in self.conflicts:
            row['zone_id'] = find(row['resource_id'])
        for route in self.paths.values():
            claims = [set(find(resource) for resource in claim) for claim in route['claims']]
            route['claims'] = [frozenset(c) for c in claims]
            route['release'] = {i: claims[i - 1] - claims[i] for i in range(1, len(claims))}
            route['merge_claims'] = {end: frozenset().union(*(
                route['claims'][i] for i in np.flatnonzero(route['merge_end'] == end)))
                for end in route['merge_trees']}
            # Acquire the pad before any zone needed by the final leg. Otherwise
            # an inbound aircraft can own the zone needed by the pad's departure,
            # while itself waiting for that departure to vacate the pad.
            final_zones = frozenset().union(*route['claims'][route['last_ring'] + 1:])
            gate = route['last_ring'] + 1
            for i, claim in enumerate(route['claims']):
                if route['merge'][i] >= 0 and route['ring'][i] < 0:
                    claim = claim | route['merge_claims'][int(route['merge_end'][i])]
                if claim & final_zones:
                    gate = min(gate, i)
                    break
            route['pad_gate'] = max(0, gate)

    def closest_node(self, point):
        distances = np.linalg.norm(np.asarray(self.nodes) - point, axis=1)
        i = int(distances.argmin())
        return i if distances[i] < .1 else -1

    def nearest_xy(self, point):
        distances = np.linalg.norm(np.asarray(self.nodes)[:, :2] - point[:2], axis=1)
        i = int(distances.argmin())
        if distances[i] > .1:
            raise ValueError('2D endpoint lacks a known 3D point')
        return self.nodes[i]

    def node(self, point):
        if not np.isfinite(point).all():
            raise ValueError('Nonfinite coordinate')
        if self.nodes:
            i = self.closest_node(point)
            if i >= 0:
                return i
        self.nodes.append(np.asarray(point).copy())
        return len(self.nodes) - 1

    def add_edge(self, xyz, layer, row, record, ring=-1):
        u, v = self.node(xyz[0]), self.node(xyz[-1])
        xyz[0], xyz[-1] = self.nodes[u], self.nodes[v]
        speed, source = self.config['speed_defaults_kts'][layer], 'assumed_kts'
        for col in record.index:
            if col.lower() in ('speed', 'speed_kts', 'speed_kt', 'speed_knot') and pd.notna(record[col]):
                speed, source = float(record[col]), col + '_assumed_kts'
                break
        if not math.isfinite(speed) or speed <= 0:
            raise ValueError('Speed must be positive')
        edge = Edge(len(self.edges), u, v, xyz, layer, int(row), speed, source, ring)
        self.edges.append(edge)
        self.adj[u].append(edge.id)

    def shortest(self, start, end):
        best, previous, queue = {start: 0.}, {}, [(0., start)]
        while queue:
            cost, node = heapq.heappop(queue)
            if cost != best[node]:
                continue
            if node == end:
                path = []
                while node != start:
                    edge = previous[node]
                    path.append(edge)
                    node = self.edges[edge].u
                return path[::-1]
            for eid in self.adj[node]:
                edge = self.edges[eid]
                new = cost + edge.seconds
                if new < best.get(edge.v, math.inf):
                    best[edge.v], previous[edge.v] = new, eid
                    heapq.heappush(queue, (new, edge.v))
        raise ValueError(f'No directed route: {start} → {end}')

    def sample(self, path, hover=False):
        positions, rings, layers, edges = [], [], [], []
        for eid in path:
            edge = self.edges[eid]
            dist = length(edge.xyz)
            # Each edge occupies an integer number of seconds: speed never exceeds field value.
            n = max(1, math.ceil(edge.seconds / self.config['step_s']))
            sample = np.linspace(0, dist[-1], n + 1)
            xyz = np.column_stack([np.interp(sample, dist, edge.xyz[:, axis]) for axis in range(3)])
            positions.extend(xyz[:-1])
            rings.extend([edge.ring] * n)
            layers.extend([edge.layer] * n)
            edges.extend([eid] * n)
        positions.append(self.nodes[self.edges[path[-1]].v])
        rings.append(-1)
        layers.append('hover' if hover else 'exit')
        edges.append(path[-1])
        if hover:
            for _ in range(math.ceil(self.config['hover_s'] / self.config['step_s'])):
                positions.append(positions[-1])
                rings.append(-1)
                layers.append('hover')
                edges.append(path[-1])
        # Only the spatial conflict neighbourhood is exclusive; reserving a whole
        # corridor could deadlock an inbound merge against an unrelated SID exit.
        claims = [set() for _ in positions]
        edges_array, xyz_array = np.array(edges), np.array(positions)
        for eid in set(edges):
            indices = np.flatnonzero(edges_array == eid)
            for resource in self.edge_resources[eid]:
                a, b = self.resource_edges[resource]
                other = b if eid == a else a
                d = self.edge_trees[other].query(xyz_array[indices])[0]
                affected = indices[d <= 2 * self.config['aircraft_buffer_m'] + 50]
                for k in affected:
                    claims[k].add(resource)
                    if k:
                        claims[k - 1].add(resource)
        claims = [frozenset(c) for c in claims]
        release = {i: claims[i - 1] - claims[i] for i in range(1, len(claims))}
        merge = np.full(len(positions), -1, dtype=int)
        merge_point = np.zeros((len(positions), 3))
        merge_end = np.full(len(positions), -1, dtype=int)
        merge_trees = {}
        for end in range(1, len(positions)):
            if rings[end] < 0 or rings[end - 1] >= 0:
                continue
            ri = rings[end]
            center, radius, z = self.ring_geometry[ri]
            start = end - 1
            while start > 0 and rings[start - 1] < 0:
                start -= 1
            points = xyz_array[start:end]
            distance_to_ring = np.hypot(np.linalg.norm(points[:, :2] - center, axis=1) - radius, points[:, 2] - z)
            close = np.flatnonzero(distance_to_ring <= 150)
            if len(close):
                guard = max(start, start + int(close[0]) - 1)
                merge[guard:end + 1] = ri
                merge_point[guard:end + 1] = positions[end]
                merge_end[guard:end + 1] = end
                merge_trees[end] = cKDTree(xyz_array[guard:end + 1])
        return {'xyz': np.array(positions), 'ring': np.array(rings),
                'layer': np.array(layers), 'edge': np.array(edges), 'path': path,
                'claims': claims, 'release': release,
                'merge': merge, 'merge_point': merge_point,
                'merge_end': merge_end, 'merge_trees': merge_trees,
                'last_ring': max((i for i, r in enumerate(rings) if r >= 0), default=-1)}


def demand(network, config, rate, seed, count=None):
    """Exact 30/30/30/10 proportions in every shuffled block of ten arrivals."""
    rng = np.random.default_rng(seed)
    names = list(network.terminals)
    weights = np.array([.1 if name == '강원' else .3 for name in names])
    assert np.isclose(weights.sum(), 1.)
    fids = sorted(network.ports)
    p1 = network.nodes[network.ports[1]]
    distances = np.array([np.linalg.norm(network.nodes[network.ports[f]][:2] - p1[:2]) for f in fids])
    priorities = np.exp(-distances / config['destination_distance_scale_m'])
    priorities /= priorities.sum()
    # Replay uses the original integer count; realised rate can round down again.
    count = int(rate * config['duration_s'] / 3600) if count is None else int(count)
    # Stratified random times; no accidental Poisson variation between capacity scenarios.
    times = (np.arange(count) + rng.random(count)) * config['duration_s'] / max(count, 1)
    blocks = []
    for _ in range(math.ceil(count / 10)):
        block = np.repeat(names, (weights * 10).astype(int))
        rng.shuffle(block)
        blocks.extend(block)
    destinations = rng.choice(fids, size=count, p=priorities)
    exits = rng.choice(names, size=count, p=weights)
    return [dict(aircraft_id=i + 1, requested_s=float(t), origin=str(blocks[i]),
                 destination_fid=int(destinations[i]), departure_destination=str(exits[i]))
            for i, t in enumerate(times)]


@dataclass
class Aircraft:
    request: dict
    route: dict
    phase: str
    start: float
    index: int = 0
    wait: int = 0
    arrival_finish: float = math.nan
    last_block: str = ''


def swept_distances(p, q, others_p, others_q):
    """Exact closest approach for straight motions within a simulation time step."""
    relative = p - others_p
    velocity = (q - p) - (others_q - others_p)
    norm = np.einsum('ij,ij->i', velocity, velocity)
    at = np.clip(-np.einsum('ij,ij->i', relative, velocity) / np.maximum(norm, 1e-20), 0, 1)
    return np.linalg.norm(relative + at[:, None] * velocity, axis=1)


def simulate(net, cfg, requests, caps, separation, scenario, seed, detail=False):
    dt = cfg['step_s']
    duration, warmup = cfg['duration_s'], cfg['warmup_s']
    threshold = 2 * cfg['aircraft_buffer_m']
    active, ground, pad_lock, reservations, merge_slots = [], {}, {}, {}, {}
    pending = {name: deque() for name in net.terminals}
    holding = {name: deque() for name in net.terminals}
    request_cursor = 0
    events, tracks, history = [], [], []
    landing_times, exit_times, takeoff_times, delays = [], [], [], []
    aircraft_state = {r['aircraft_id']: dict(r, admitted_s=math.nan, landed_s=math.nan,
                      departure_s=math.nan, exited_s=math.nan, air_wait_s=0,
                      outer_distance_m=0., inner_distance_m=0.) for r in requests}
    blocks = defaultdict(int)
    peak_active = peak_total = 0
    peak_ring = [0, 0]
    min_distance = math.inf
    min_ring_distance = math.inf
    arrivals_in_measurement = 0
    peak_ring_moving = [0, 0]
    ring_moving_steps, ring_occupied_steps = [0, 0], [0, 0]
    seconds_ten_moving = [0., 0.]

    def log(t, a, event):
        if detail:
            events.append(dict(time_s=t, aircraft_id=a.request['aircraft_id'], event=event,
                               fid=a.request['destination_fid'], phase=a.phase))

    for tick in range(int(duration / dt) + 1):
        t = tick * dt
        while request_cursor < len(requests) and requests[request_cursor]['requested_s'] <= t:
            r = requests[request_cursor]
            pending[r['origin']].append(r)
            request_cursor += 1
        for name in net.terminals:
            hid = net.holding_for[name]
            capacity = net.holdings[hid]['capacity']
            while pending[name] and len(holding[name]) < capacity:
                holding[name].append(pending[name].popleft())
        # Landing is completed only after the 20 s hover samples have elapsed.
        remaining = []
        for a in active:
            if a.index == len(a.route['xyz']) - 1:
                fid = a.request['destination_fid']
                state = aircraft_state[a.request['aircraft_id']]
                state['air_wait_s'] += a.wait * dt
                for resource, owner in list(reservations.items()):
                    if owner == a.request['aircraft_id']:
                        del reservations[resource]
                merge_slots.pop(a.request['aircraft_id'], None)
                if a.phase == 'arrival':
                    ground[fid] = (a, t + cfg['turnaround_s'])
                    a.arrival_finish = t
                    state['landed_s'] = t
                    landing_times.append(t)
                    delays.append((a.request['requested_s'], t - a.request['requested_s'] - (len(a.route['xyz']) - 1) * dt))
                    log(t, a, 'landed')
                else:
                    state['exited_s'] = t
                    exit_times.append(t)
                    log(t, a, 'exited')
            else:
                remaining.append(a)
        active = remaining
        # Pad remains exclusive until the departing vehicle has cleared 100 m.
        for fid, owner in list(pad_lock.items()):
            if fid in ground:
                continue
            departures = [a for a in active if a.request['aircraft_id'] == owner and a.phase == 'departure']
            if departures and np.linalg.norm(departures[0].route['xyz'][departures[0].index] - net.nodes[net.ports[fid]]) >= threshold:
                del pad_lock[fid]
        # Admission and takeoff are checked against airborne traffic. Ground vehicles
        # occupy pads, but are not included in the airborne spherical-buffer test.
        candidates = []
        for name, queue in holding.items():
            if queue:
                r = queue[0]
                route = net.paths[('arrival', name, r['destination_fid'])]
                candidates.append((0, r['requested_s'], 'arrival', name, r, route))
        for fid, (a, ready) in ground.items():
            if ready <= t:
                r = a.request
                route = net.paths[('departure', r['departure_destination'], fid)]
                candidates.append((1, ready, 'departure', fid, r, route))
        for _, _, phase, source, r, route in sorted(candidates, key=lambda x: (x[0], x[1], x[4]['aircraft_id'])):
            claim = route['claims'][0]
            needs_pad = phase == 'arrival' and route['pad_gate'] == 0
            if needs_pad and r['destination_fid'] in pad_lock:
                continue
            if any(resource in reservations and reservations[resource] != r['aircraft_id'] for resource in claim):
                continue
            if active and min(np.linalg.norm(route['xyz'][0] - a.route['xyz'][a.index]) for a in active) < threshold:
                continue
            if phase == 'departure':
                # Give airborne arrivals near the destination landing segment priority.
                if any(a.phase == 'arrival' and a.request['destination_fid'] == source
                       and a.route['layer'][a.index] in ('Short_Final', 'hover') for a in active):
                    continue
                del ground[source]
                aircraft_state[r['aircraft_id']]['departure_s'] = t
                takeoff_times.append(t)
            else:
                holding[source].popleft()
                aircraft_state[r['aircraft_id']]['admitted_s'] = t
                if t >= warmup:
                    arrivals_in_measurement += 1
            a = Aircraft(r, route, phase, t)
            if needs_pad:
                pad_lock[r['destination_fid']] = r['aircraft_id']
            for resource in claim:
                reservations[resource] = r['aircraft_id']
            active.append(a)
            log(t, a, 'admitted' if phase == 'arrival' else 'takeoff')
        n = len(active)
        moving_count = [0, 0]
        ring_count = [sum(a.route['ring'][a.index] == ri for a in active) for ri in (0, 1)]
        if n:
            old = np.array([a.route['xyz'][a.index] for a in active])
            proposed = old.copy()
            current_ring = np.array([a.route['ring'][a.index] for a in active])
            new_ring = current_ring.copy()
            moves = np.zeros(n, dtype=bool)

            def priority(i):
                a = active[i]
                nxt = min(a.index + 1, len(a.route['xyz']) - 1)
                merge = a.route['merge'][nxt] >= 0 and a.route['ring'][a.index] < 0
                final = a.phase == 'arrival' and a.route['layer'][a.index] in ('Short_Final', 'hover')
                return (0 if final else 1, 0 if merge else 1, a.start, a.request['aircraft_id'])

            for i in sorted(range(n), key=priority):
                a = active[i]
                nxt = min(a.index + 1, len(a.route['xyz']) - 1)
                q = a.route['xyz'][nxt]
                r0, r1 = current_ring[i], int(a.route['ring'][nxt])
                fid, aid = a.request['destination_fid'], a.request['aircraft_id']
                merge_ring = int(a.route['merge'][nxt])
                entry_point = a.route['merge_point'][nxt]
                merge_tree = a.route['merge_trees'].get(int(a.route['merge_end'][nxt]))
                if merge_ring >= 0 and r0 < 0 and aid not in merge_slots:
                    occupied_slots = sum(slot[0] == merge_ring for slot in merge_slots.values())
                    if ring_count[merge_ring] + occupied_slots >= caps[merge_ring]:
                        blocks['ring_capacity'] += 1
                        a.last_block = 'merge_capacity'
                        continue
                    on_ring = new_ring == merge_ring
                    on_ring[i] = False
                    if (on_ring.any() and (np.any(np.linalg.norm(proposed[on_ring] - entry_point, axis=1) < separation) or
                                          np.any(merge_tree.query(proposed[on_ring])[0] < threshold + 20))) or any(
                            slot[0] == merge_ring and (np.linalg.norm(slot[1] - entry_point) < separation or
                                np.min(slot[2].query(merge_tree.data)[0]) < threshold + 20)
                            for slot in merge_slots.values()):
                        blocks['merge_gap'] += 1
                        a.last_block = 'merge_gap'
                        continue
                if r1 >= 0 and any(owner != aid and slot[0] == r1 and
                                   (np.linalg.norm(q - slot[1]) < separation or slot[2].query(q)[0] < threshold + 20)
                                   for owner, slot in merge_slots.items()):
                    blocks['merge_priority'] += 1
                    a.last_block = 'merge_priority:' + str(list(merge_slots))
                    continue
                final_entry = a.phase == 'arrival' and a.route['layer'][nxt] in ('Short_Final', 'hover')
                claim = a.route['claims'][nxt]
                if merge_ring >= 0 and r0 < 0:
                    claim = claim | a.route['merge_claims'][int(a.route['merge_end'][nxt])]
                if any(resource in reservations and reservations[resource] != aid for resource in claim):
                    blocks['corridor_reservation'] += 1
                    a.last_block = 'reservation:' + str({r: reservations[r] for r in claim if r in reservations and reservations[r] != aid})
                    continue
                # Reserve the landing pad at the start of the final corridor,
                # so an arrival never occupies a corridor while awaiting takeoff.
                terminal_entry = a.phase == 'arrival' and nxt >= a.route['pad_gate']
                if terminal_entry and fid in pad_lock and pad_lock[fid] != aid:
                    blocks['pad'] += 1
                    a.last_block = 'pad:' + str(pad_lock[fid])
                    continue
                if final_entry and fid in pad_lock and pad_lock[fid] != aid:
                    blocks['pad'] += 1
                    a.last_block = 'pad:' + str(pad_lock[fid])
                    continue
                if r1 >= 0 and r0 != r1 and ring_count[r1] >= caps[r1]:
                    blocks['ring_capacity'] += 1
                    a.last_block = 'ring_capacity'
                    continue
                others = np.arange(n) != i
                if others.any():
                    distances = swept_distances(old[i], q, old[others], proposed[others])
                    required = np.full(n - 1, threshold)
                    # Ring spacing is a 3D/chord separation for co-ring vehicles.
                    same = ((current_ring[others] == r0) & (r0 >= 0)) | ((new_ring[others] == r1) & (r1 >= 0))
                    required[same] = max(threshold, separation)
                    if np.any(distances < required - 1e-7):
                        blocks['separation'] += 1
                        a.last_block = 'separation:' + str([active[j].request['aircraft_id'] for j in np.flatnonzero(others)[distances < required - 1e-7]])
                        continue
                moves[i], proposed[i], new_ring[i] = True, q, r1
                a.last_block = ''
                if merge_ring >= 0 and r0 < 0:
                    merge_slots[aid] = (merge_ring, entry_point, merge_tree)
                if r1 >= 0:
                    merge_slots.pop(aid, None)
                for resource in claim:
                    reservations[resource] = aid
                if r0 != r1:
                    if r0 >= 0:
                        ring_count[r0] -= 1
                    if r1 >= 0:
                        ring_count[r1] += 1
                if final_entry or terminal_entry:
                    pad_lock[fid] = aid
            # Validate actual simultaneous swept motion, independently of priority.
            for i in range(n - 1):
                d = swept_distances(old[i], proposed[i], old[i + 1:], proposed[i + 1:])
                min_distance = min(min_distance, float(d.min()))
                if np.any(d < threshold - 1e-6):
                    raise AssertionError(f'Buffer overlap at t={t}, aircraft {active[i].request["aircraft_id"]}')
                same = (((current_ring[i + 1:] == current_ring[i]) & (current_ring[i] >= 0)) |
                        ((new_ring[i + 1:] == new_ring[i]) & (new_ring[i] >= 0)))
                if same.any():
                    min_ring_distance = min(min_ring_distance, float(d[same].min()))
                    if np.any(d[same] < separation - 1e-6):
                        raise AssertionError(f'Ring spacing violation at t={t}')
            for i, a in enumerate(active):
                ri = int(current_ring[i])
                if ri >= 0:
                    ring_occupied_steps[ri] += 1
                    if moves[i] and np.linalg.norm(proposed[i] - old[i]) > 1e-8:
                        moving_count[ri] += 1
                        ring_moving_steps[ri] += 1
                        field = 'outer_distance_m' if ri == 0 else 'inner_distance_m'
                        aircraft_state[a.request['aircraft_id']][field] += float(np.linalg.norm(proposed[i] - old[i]))
                if moves[i]:
                    a.index += 1
                    keep = a.route['claims'][a.index]
                    if a.route['merge'][a.index] >= 0 and a.route['ring'][a.index] < 0:
                        keep = keep | a.route['merge_claims'][int(a.route['merge_end'][a.index])]
                    for resource, owner in list(reservations.items()):
                        if owner == a.request['aircraft_id'] and resource not in keep:
                            del reservations[resource]
                else:
                    a.wait += 1
            if detail and tick % 10 == 0:
                for i, a in enumerate(active):
                    tracks.append(dict(time_s=t, aircraft_id=a.request['aircraft_id'], phase=a.phase,
                        x=proposed[i, 0], y=proposed[i, 1], z=proposed[i, 2],
                        edge_id=int(a.route['edge'][a.index]), ring=int(new_ring[i])))
        peak_active = max(peak_active, n)
        peak_total = max(peak_total, n + len(ground) + sum(map(len, holding.values())))
        peak_ring = [max(peak_ring[i], ring_count[i]) for i in (0, 1)]
        for ri in (0, 1):
            peak_ring_moving[ri] = max(peak_ring_moving[ri], moving_count[ri])
            if moving_count[ri] >= 10:
                seconds_ten_moving[ri] += dt
        if tick % 10 == 0:
            history.append(dict(time_s=t, airborne=n, outer=ring_count[0], inner=ring_count[1],
                outer_moving=moving_count[0], inner_moving=moving_count[1],
                grounded=len(ground), holding=sum(map(len, holding.values())),
                external_queue=sum(map(len, pending.values())),
                total_backlog=request_cursor - len(landing_times),
                roundtrip_backlog=request_cursor - len(exit_times),
                landed=len(landing_times), exited=len(exit_times)))
    history = pd.DataFrame(history)
    measured = history[history.time_s >= warmup]
    hours = (duration - warmup) / 3600
    arrivals = sum(warmup <= t <= duration for t in landing_times)
    departures = sum(warmup <= t <= duration for t in takeoff_times)
    exits = sum(warmup <= t <= duration for t in exit_times)
    cutoff = max(0, duration - (1200 if duration >= 2400 else 900))
    mature = [s for s in aircraft_state.values() if warmup <= s['requested_s'] <= cutoff]
    completed = sum(np.isfinite(s['exited_s']) for s in mature)
    completion = completed / len(mature) if mature else math.nan
    ds = []
    for state in mature:
        if np.isfinite(state['exited_s']):
            arrival = net.paths[('arrival', state['origin'], state['destination_fid'])]
            departure = net.paths[('departure', state['departure_destination'], state['destination_fid'])]
            baseline = (len(arrival['xyz']) + len(departure['xyz']) - 2) * dt + cfg['turnaround_s']
            ds.append(state['exited_s'] - state['requested_s'] - baseline)
    for a in active:
        aircraft_state[a.request['aircraft_id']]['air_wait_s'] += a.wait * dt
        aircraft_state[a.request['aircraft_id']].update(final_phase=a.phase, final_index=a.index,
            final_edge=int(a.route['edge'][a.index]), final_block_reason=a.last_block,
            held_zones=str([r for r, owner in reservations.items() if owner == a.request['aircraft_id']]))
    growth = float(np.polyfit(measured.time_s, measured.roundtrip_backlog, 1)[0] * 3600) if len(measured) > 1 else 0
    metrics = dict(scenario=scenario, seed=seed, separation_m=separation,
        outer_capacity=caps[0], inner_capacity=caps[1], requested_arrivals=len(requests),
        arrival_demand_per_h=len(requests) * 3600 / duration,
        landings_per_h=arrivals / hours, departures_per_h=departures / hours,
        external_exits_per_h=exits / hours,
        movements_per_h=(arrivals + departures) / hours,
        peak_airborne=peak_active, peak_system_including_holding_ground=peak_total,
        peak_outer=peak_ring[0], peak_inner=peak_ring[1],
        peak_moving_outer=peak_ring_moving[0], peak_moving_inner=peak_ring_moving[1],
        outer_seconds_at_least_10_moving=seconds_ten_moving[0],
        inner_seconds_at_least_10_moving=seconds_ten_moving[1],
        outer_motion_fraction=ring_moving_steps[0] / max(1, ring_occupied_steps[0]),
        inner_motion_fraction=ring_moving_steps[1] / max(1, ring_occupied_steps[1]),
        end_external_queue=int(history.iloc[-1].external_queue), end_backlog=int(history.iloc[-1].total_backlog),
        completion_fraction=completion, queue_growth_per_h=growth,
        mean_delay_s=float(np.mean(ds)) if ds else math.nan,
        p95_delay_s=float(np.percentile(ds, 95)) if ds else math.nan,
        min_swept_separation_m=min_distance if np.isfinite(min_distance) else math.nan,
        min_ring_swept_separation_m=min_ring_distance if np.isfinite(min_ring_distance) else math.nan,
        mature_aircraft=len(mature), completed_mature_aircraft=completed,
        buffer_violations=0, blocked_pad_steps=blocks['pad'],
        blocked_reservation_steps=blocks['corridor_reservation'],
        blocked_merge_steps=blocks['merge_gap'] + blocks['merge_priority'],
        blocked_ring_steps=blocks['ring_capacity'], blocked_separation_steps=blocks['separation'])
    metrics['stable'] = bool(completion >= cfg['stable_completion_fraction'] and
        metrics['p95_delay_s'] <= cfg['stable_p95_delay_s'] and growth <= cfg['stable_queue_growth_per_hour'])
    return metrics, history, list(aircraft_state.values()), tracks, events


def capacity_expansion_history(net, cfg):
    rng = np.random.default_rng(cfg['seed'])
    rows = []
    for separation in cfg['separations_m']:
        radii = np.array(net.ring_lengths) / (2 * math.pi)
        maximum = np.floor(math.pi / np.arcsin(separation / (2 * radii))).astype(int)
        cap = np.minimum([14, 14], maximum)
        step = 0
        rows.append(dict(separation_m=separation, step=step, changed='initial',
                         outer=int(cap[0]), inner=int(cap[1])))
        targets = sorted(set([14, *cfg['capacity_checkpoints'], int(max(maximum))]))
        for target in targets:
            # Random increments are reproducible and never exceed geometric packing.
            while min(cap) < target and np.any(cap < maximum):
                choices = np.flatnonzero(cap < np.minimum(maximum, target))
                if not len(choices):
                    break
                ri = int(rng.choice(choices))
                cap[ri] += 1
                step += 1
                rows.append(dict(separation_m=separation, step=step,
                    changed='outer' if ri == 0 else 'inner', outer=int(cap[0]), inner=int(cap[1])))
    return rows


def capacity_scenarios(net, cfg):
    rows = []
    history = capacity_expansion_history(net, cfg)
    for separation in cfg['separations_m']:
        radii = np.array(net.ring_lengths) / (2 * math.pi)
        maximum = np.floor(math.pi / np.arcsin(separation / (2 * radii))).astype(int)
        saved = set()
        for target in sorted(set([14, *cfg['capacity_checkpoints'], int(max(maximum))])):
            state = next(r for r in history if r['separation_m'] == separation and
                         r['outer'] >= min(target, maximum[0]) and r['inner'] >= min(target, maximum[1]))
            cap = np.array([state['outer'], state['inner']])
            pair = tuple(map(int, cap))
            if pair in saved:
                continue
            saved.add(pair)
            rows.append(dict(scenario=f'S{separation}_O{pair[0]}_I{pair[1]}', separation_m=separation,
                outer_capacity=pair[0], inner_capacity=pair[1],
                outer_geometric_max=int(maximum[0]), inner_geometric_max=int(maximum[1]),
                outer_mean_spacing_m=net.ring_lengths[0] / pair[0], inner_mean_spacing_m=net.ring_lengths[1] / pair[1],
                outer_uniform_chord_m=2 * radii[0] * math.sin(math.pi / pair[0]),
                inner_uniform_chord_m=2 * radii[1] * math.sin(math.pi / pair[1])))
    return rows


def configure_plot():
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    if 'Malgun Gothic' in fonts:
        plt.rcParams['font.family'] = 'Malgun Gothic'
    plt.rcParams.update({'axes.unicode_minus': False, 'figure.facecolor': '#f7f9fc',
                         'savefig.facecolor': '#f7f9fc', 'font.size': 10})


def init_worker(cfg):
    global WORKER_NET, WORKER_CFG
    WORKER_CFG = cfg
    WORKER_NET = Network(cfg)


def run_batch(job):
    scenarios, rate, seed = job
    net, cfg = WORKER_NET, WORKER_CFG
    requests = demand(net, cfg, rate, seed)
    cached, batch = None, []
    for scenario in scenarios:
        if cached is not None:
            metrics, states = dict(cached[0]), cached[1]
            metrics.update(scenario=scenario['scenario'], outer_capacity=scenario['outer_capacity'],
                           inner_capacity=scenario['inner_capacity'], computation='identical_nonbinding_capacity')
        else:
            metrics, _, states, _, _ = simulate(net, cfg, requests,
                (scenario['outer_capacity'], scenario['inner_capacity']), scenario['separation_m'],
                scenario['scenario'], seed)
            metrics['computation'] = 'simulated'
            if metrics['blocked_ring_steps'] == 0:
                cached = (dict(metrics), states)
        batch.append((metrics, states, rate))
        metrics['target_arrivals_per_h'] = rate
    return batch


def plots(net, results, scenarios, history, out):
    configure_plot()
    colors = {'접근': '#8995a3', 'Expert_Off': '#8995a3', '외곽_2800': '#df8b28', '내곽_2000': '#277eb6',
              'Interchange_In': '#159e91', 'Interchange_out': '#69c2ba',
              'IAF_FAF_In_Update': '#519ab8', 'IAF_FAF_Out_Update': '#d18aaf', 'Short_Final': '#8952ae'}
    origin = np.mean([net.nodes[n] for n in net.ports.values()], axis=0)
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(121)
    ax3 = fig.add_subplot(122, projection='3d')
    seen = set()
    for edge in net.edges:
        xyz = edge.xyz - origin
        xyz[:, :2] /= 1000
        xyz[:, 2] = edge.xyz[:, 2]
        label = edge.layer if edge.layer not in seen else None
        seen.add(edge.layer)
        ax.plot(xyz[:, 0], xyz[:, 1], color=colors[edge.layer], lw=1, label=label)
        ax3.plot(*xyz.T, color=colors[edge.layer], lw=.8)
        if edge.ring >= 0 and edge.row % 3 == 0 and len(xyz) > 2:
            mid = len(xyz) // 2
            ax.annotate('', xy=xyz[min(len(xyz)-1, mid+max(1, len(xyz)//6)), :2],
                        xytext=xyz[mid, :2], arrowprops={'arrowstyle': '-|>', 'color': colors[edge.layer], 'lw': 1.2})
    for name, node in net.terminals.items():
        xy = (net.nodes[node][:2] - origin[:2]) / 1000
        ax.scatter(*xy, s=25, color='black')
        ax.annotate(name, xy, xytext=(5, 5), textcoords='offset points')
    for fid, node in net.ports.items():
        xy = (net.nodes[node][:2] - origin[:2]) / 1000
        ax.scatter(*xy, marker='s', color='#8952ae', s=18)
        ax.annotate(str(fid), xy, fontsize=8)
    ax.set_aspect('equal'); ax.set_xlabel('동서 상대거리 (km)'); ax.set_ylabel('남북 상대거리 (km)')
    ax.set_title('방향성 네트워크: 외곽 반시계 / 내곽 시계')
    ax.legend(fontsize=7, loc='lower left')
    ax3.set_zlabel('원본 고도 (m)'); ax3.set_xlabel('동서 (km)'); ax3.set_ylabel('남북 (km)')
    ax3.set_title('3D 연결 · XY 교차를 접속점으로 처리하지 않음')
    fig.tight_layout(); fig.savefig(out / 'network_3d.png', dpi=180, bbox_inches='tight'); plt.close(fig)
    separations = sorted(results.separation_m.unique())
    columns = min(2, len(separations))
    fig, axes = plt.subplots(math.ceil(len(separations) / columns), columns,
                             figsize=(7 * columns, 5 * math.ceil(len(separations) / columns)), squeeze=False)
    for ax in axes.flat[len(separations):]:
        ax.set_visible(False)
    for sep, ax in zip(separations, axes.flat):
        for scenario, data in results[results.separation_m == sep].groupby('scenario'):
            grouped = data.groupby('arrival_demand_per_h').movements_per_h.agg(['mean', 'min', 'max'])
            ax.plot(grouped.index, grouped['mean'], marker='o', ms=4, label=scenario.split('_', 1)[1])
            ax.fill_between(grouped.index, grouped['min'], grouped['max'], alpha=.12)
        ax.set_title(f'링 분리거리 {sep:g} m')
        ax.set_xlabel('요구 도착량 (대/h)'); ax.set_ylabel('착륙+이륙 (회/h)'); ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle('수요·링 용량별 처리량 (밴드: 반복 실행 최솟값~최댓값)', fontsize=16)
    fig.tight_layout(); fig.savefig(out / 'capacity_comparison.png', dpi=180); plt.close(fig)
    fig, axes = plt.subplots(3, 1, figsize=(13, 10), sharex=True)
    for col in ('airborne', 'outer', 'inner', 'holding'):
        axes[0].plot(history.time_s / 60, history[col], label=col)
    for col in ('total_backlog', 'roundtrip_backlog', 'external_queue'):
        axes[1].plot(history.time_s / 60, history[col], label=col)
    for col in ('landed', 'exited'):
        axes[2].plot(history.time_s / 60, history[col], label=col)
    for ax in axes:
        ax.legend(); ax.grid(alpha=.2)
    axes[0].set_ylabel('대수'); axes[1].set_ylabel('미완료 도착 대수'); axes[2].set_ylabel('누적 완료')
    axes[2].set_xlabel('시뮬레이션 시간 (분)')
    fig.suptitle('선정 시나리오의 점유·대기열·처리량'); fig.tight_layout()
    fig.savefig(out / 'selected_timeline.png', dpi=180); plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('capacity_config.json'))
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--build-only', action='store_true')
    parser.add_argument('--workers', type=int, default=min(6, os.cpu_count() or 1))
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text(encoding='utf-8'))
    if args.quick:
        cfg.update(duration_s=1200, warmup_s=120, replicates=1,
                   separations_m=[100], capacity_checkpoints=[14], arrival_rates_per_hour=[30, 90])
    validate_config(cfg)
    args.output.mkdir(parents=True, exist_ok=True)
    net = Network(cfg)
    csv(net.audit, args.output / 'input_audit.csv')
    csv(net.conflicts, args.output / 'conflict_corridors.csv')
    csv(net.route_rows, args.output / 'OD_routes.csv')
    csv([dict(edge_id=e.id, start_node=e.u, end_node=e.v, layer=e.layer, source_row=e.row,
              length_m=length(e.xyz)[-1], travel_time_s=e.seconds, speed_kts=e.speed,
              speed_source=e.speed_source, ring=e.ring) for e in net.edges], args.output / 'network_edges.csv')
    csv([dict(node_id=i, x=p[0], y=p[1], z=p[2]) for i, p in enumerate(net.nodes)], args.output / 'network_nodes.csv')
    print(f'Network: {len(net.nodes)} nodes, {len(net.edges)} directed edges, {len(net.paths)} OD routes', flush=True)
    if args.build_only:
        return
    scenarios = capacity_scenarios(net, cfg)
    if args.quick:
        scenarios = scenarios[:1]
    csv(scenarios, args.output / 'capacity_scenarios.csv')
    csv(capacity_expansion_history(net, cfg), args.output / 'capacity_expansion_history.csv')
    (args.output / 'effective_config.json').write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding='utf-8')
    results, od_rows = [], []
    total = len(scenarios) * len(cfg['arrival_rates_per_hour']) * cfg['replicates']
    started = time.monotonic()
    jobs = [([s for s in scenarios if s['separation_m'] == sep], rate, cfg['seed'] + rep)
            for sep in cfg['separations_m'] for rate in cfg['arrival_rates_per_hour']
            for rep in range(cfg['replicates'])]
    def collect(batch):
        for metrics, states, rate in batch:
            results.append(metrics)
            for state in states:
                base = dict(scenario=metrics['scenario'], seed=metrics['seed'], demand_per_h=rate, **state)
                od_rows.append(dict(base, movement='arrival', origin_name=state['origin'],
                                   destination_name=str(state['destination_fid'])))
                od_rows.append(dict(base, movement='departure', origin_name=str(state['destination_fid']),
                                   destination_name=state['departure_destination']))
            print(f'[{len(results)}/{total}] {metrics["scenario"]} rate={rate} seed={metrics["seed"]} '
                  f'movements={metrics["movements_per_h"]:.1f}/h backlog={metrics["end_backlog"]} '
                  f'stable={metrics["stable"]} elapsed={time.monotonic()-started:.1f}s', flush=True)
        csv(results, args.output / 'capacity_results.csv')
    if args.workers == 1 or args.quick:
        init_worker(cfg)
        for job in jobs:
            collect(run_batch(job))
    else:
        with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=(cfg,)) as pool:
            futures = [pool.submit(run_batch, job) for job in jobs]
            for future in as_completed(futures):
                collect(future.result())
    od = pd.DataFrame(od_rows)
    inbound = od.movement == 'arrival'
    od['movement_ready_s'] = np.where(inbound, od.requested_s, od.landed_s + cfg['turnaround_s'])
    od['movement_started_s'] = np.where(inbound, od.admitted_s, od.departure_s)
    od['movement_finished_s'] = np.where(inbound, od.landed_s, od.exited_s)
    od['movement_status'] = np.where(od.movement_finished_s.notna(), 'completed',
        np.where(od.movement_started_s.notna(), 'in_progress',
                 np.where(od.movement_ready_s.notna(), 'waiting', 'awaiting_landing')))
    od.to_csv(args.output / 'OD.csv', index=False, encoding='utf-8-sig')
    results = pd.DataFrame(results)
    results = results.sort_values(['separation_m', 'outer_capacity', 'inner_capacity', 'target_arrivals_per_h', 'seed'])
    results.to_csv(args.output / 'capacity_results.csv', index=False, encoding='utf-8-sig')
    grouped = results.groupby(['scenario', 'separation_m', 'outer_capacity', 'inner_capacity', 'arrival_demand_per_h'])
    aggregate = grouped.agg(stable_all_replicates=('stable', 'all'),
        stable_replicates=('stable', 'sum'), movements_mean_per_h=('movements_per_h', 'mean'),
        movements_min_per_h=('movements_per_h', 'min'), movements_max_per_h=('movements_per_h', 'max'),
        peak_airborne_max=('peak_airborne', 'max'), p95_delay_worst_s=('p95_delay_s', 'max'),
        completion_min=('completion_fraction', 'min'), queue_growth_worst_per_h=('queue_growth_per_h', 'max')).reset_index()
    aggregate.to_csv(args.output / 'capacity_summary.csv', index=False, encoding='utf-8-sig')
    viable = aggregate[aggregate.stable_all_replicates]
    chosen = (viable if len(viable) else aggregate).sort_values(
        ['movements_mean_per_h', 'p95_delay_worst_s'], ascending=[False, True]).iloc[0]
    selected = results[(results.scenario == chosen.scenario) &
                       (results.arrival_demand_per_h == chosen.arrival_demand_per_h)].sort_values('seed').iloc[0]
    requests = demand(net, cfg, selected.target_arrivals_per_h, int(selected.seed),
                      count=int(selected.requested_arrivals))
    _, history, states, tracks, events = simulate(net, cfg, requests,
        (int(selected.outer_capacity), int(selected.inner_capacity)), selected.separation_m,
        selected.scenario, int(selected.seed), detail=True)
    history.to_csv(args.output / 'selected_timeline.csv', index=False, encoding='utf-8-sig')
    csv(tracks, args.output / 'selected_trajectories.csv')
    csv(events, args.output / 'selected_events.csv')
    plots(net, results, scenarios, history, args.output)
    summary = dict(selection_status='stable_tested' if len(viable) else 'diagnostic_only_no_stable_condition',
                   selected=selected.to_dict(), selected_aggregate=chosen.to_dict(),
                   max_tested_stable_arrivals_per_h=float(viable.arrival_demand_per_h.max()) if len(viable) else None,
                   max_stable_observed_airborne=int(viable.peak_airborne_max.max()) if len(viable) else None,
                   runs=len(results), reachable_ODs=len(net.paths),
                   max_observed_airborne=int(results.peak_airborne.max()),
                   max_observed_movements_per_h=float(results.movements_per_h.max()),
                   stable_runs=int(results.stable.sum()))
    (args.output / 'summary.json').write_text(json_text(summary), encoding='utf-8')
    print(json_text(summary), flush=True)


if __name__ == '__main__':
    main()
