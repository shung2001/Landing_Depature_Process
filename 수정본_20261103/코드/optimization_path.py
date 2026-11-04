"""Joint clockwise STAR/SID, opposite FAFs, spacing, and ring-boundary scenarios.

Run from any directory: python -B 수정본_20261103/코드/optimization_path.py
Inputs are never modified. Reviewed output: ../결과물/optimization_path_reviewed.
--star-only runs the earlier STAR-only analysis in ../결과물/optimization_path.

Model: EPSG:5186 XY, supplied ring altitude, FAF `alt`, constant climb/descent
gradient along analytic horizontal arclength. STAR 70 -> 60 kt, SID 60 -> 70 kt,
constant tangential acceleration. TAS in still air; no wind/IAS conversion.
RSR is constructed geometrically using the external tangent of two right-turn
circles. STAR heads toward the OLS apex; SID departs away from the OLS apex
and joins the clockwise ring tangent. SID is NOT a reversed STAR geometry.
Turn radius >= v_start**2/(g*tan(bank_limit)); this is a coordinated-turn
research model, NOT a type-specific FAA/NASA approved operating radius.

Every direction x IAF x FAF x radius combination is evaluated. A spatially
diverse shortlist is jointly optimized by multistart block coordinate descent.
Default cases: inside_only, outside_allowed, outside_detour. Opposite FAFs
are assigned to STAR/SID for each pair. XY route intersections have a large
penalty; minimum and pairwise XY spacing deficits are penalized. The 100 m
soft spacing target and weights are configurable. No vertical separation or
departure-time simulation is applied. Landing_Area no-crossing is tried first,
then softened only if full opposite-FAF pair coverage is impossible in the
enumerated individual candidates. No global optimality claim is made.
Landing_Area is a zero-thickness Polygon Z surface, not an occupied volume.

References:
https://www.faa.gov/sites/faa.gov/files/pilots/pilot_handbook.pdf
https://ntrs.nasa.gov/citations/20240005129
https://ompl.kavrakilab.org/classompl_1_1base_1_1DubinsStateSpace.html
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import shutil
import time

import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely import get_coordinates
from shapely.geometry import LineString, Point, Polygon

BASE = Path(__file__).resolve().parents[1]
KNOT = 1852 / 3600
G = 9.80665
TAU = 2 * math.pi


@dataclass
class Config:
    iaf_count: int = 1000
    ring_z: float = 350.0
    start_kt: float = 70.0
    end_kt: float = 60.0
    max_bank: float = 30.0
    max_slope_deg: float = 2.36
    max_ft_nm: float = 250.0
    radii: tuple = ()  # empty -> computed minimum, 350, 500, 750 m
    step_m: float = 12.0
    chord_error_m: float = 0.02
    endpoint_tol_m: float = 0.001
    separation: bool = False
    horizontal_sep_m: float = 100.0
    vertical_sep_m: float = 70.0
    node_neighborhood_m: float = 250.0
    node_weight: float = 20.0
    crossing_weight: float = 10000.0
    overlap_weight: float = 10.0
    clearance_weight: float = 50000.0
    clearance_scale_m: float = 100.0  # smooth penalty scale; not hard separation
    length_weight: float = 1.0  # per km
    bank_weight: float = 1.0
    shortlist_per_pair: int = 25
    multistarts: int = 16
    seed: int = 20261104


def load(path):
    data = gpd.read_file(path)
    if data.empty or data.crs is None:
        raise ValueError(f'Empty input or missing CRS: {path}')
    data = data.to_crs(5186)
    if data.geometry.is_empty.any() or data.geometry.isna().any():
        raise ValueError(f'Empty geometry: {path}')
    return data


def single_point(geometry):
    if geometry.geom_type == 'MultiPoint' and len(geometry.geoms) == 1:
        return geometry.geoms[0]
    if geometry.geom_type == 'Point':
        return geometry
    raise ValueError('Each FAF/IAF feature must contain exactly one point.')


def fingerprint(paths):
    result = {}
    for path in sorted(set(paths)):
        for suffix in ('.shp', '.shx', '.dbf', '.prj', '.cpg'):
            item = path.with_suffix(suffix)
            if item.exists():
                result[str(item.relative_to(BASE))] = hashlib.sha256(item.read_bytes()).hexdigest()
    return result


def prepare(qgis, config):
    """Validate pairs, fit supplied ring, and associate FAF with OLS by XY."""
    iaf_path = qgis / 'IAF_후보군/내곽_꼭짓점.shp'
    sources = [iaf_path, qgis / 'FAF.shp']
    source = load(iaf_path)
    xy = np.unique(np.array([(p.x, p.y) for p in map(single_point, source.geometry)]), axis=0)
    # Fit locally centered coordinates to avoid squared Korean-grid roundoff.
    origin = xy.mean(axis=0)
    local = xy - origin
    fit, *_ = np.linalg.lstsq(np.column_stack((2 * local, np.ones(len(local)))),
                             np.sum(local * local, axis=1), rcond=None)
    center = origin + fit[:2]
    ring_radius = float(np.sqrt(fit[2] + np.dot(fit[:2], fit[:2])))
    residual = float(np.max(np.abs(np.linalg.norm(xy - center, axis=1) - ring_radius)))
    if residual > 0.05:
        raise ValueError(f'Input IAF points are not circular; max fit residual {residual} m.')
    theta0 = math.atan2(single_point(source.geometry.iloc[0]).y - center[1],
                        single_point(source.geometry.iloc[0]).x - center[0])
    theta = theta0 - np.arange(config.iaf_count) * TAU / config.iaf_count
    nodes = center + ring_radius * np.column_stack((np.cos(theta), np.sin(theta)))
    headings = theta - math.pi / 2  # CW tangent, independent of SHP vertex order
    areas = []
    for path in sorted((qgis / 'Landing_Area').glob('*.shp')):
        sources.append(path)
        for _, row in load(path).iterrows():
            geom = row.geometry
            if geom.geom_type != 'Polygon' or not geom.is_valid or len(geom.interiors):
                raise ValueError(f'Expected valid single triangle without holes: {path}')
            xyz = np.asarray(geom.exterior.coords, dtype=float)[:-1]
            if xyz.shape != (3, 3) or not np.isfinite(xyz).all():
                raise ValueError(f'Expected finite triangular Polygon Z: {path}')
            apex = np.array([float(row.center_x), float(row.center_y)])
            idx = int(np.argmin(np.linalg.norm(xyz[:, :2] - apex, axis=1)))
            if np.linalg.norm(xyz[idx, :2] - apex) > .01:
                raise ValueError(f'OLS apex attribute does not match geometry: {path}')
            mid = np.delete(xyz, idx, axis=0).mean(axis=0)
            # Centered plane avoids poorly conditioned intercept on large XY.
            plane_origin = xyz[0].copy()
            gradient = np.linalg.solve(xyz[1:, :2] - xyz[0, :2], xyz[1:, 2] - xyz[0, 2])
            areas.append(dict(pair=path.stem, fid=int(row.fid), geometry=Polygon(xyz[:, :2]),
                              xyz=xyz, apex=xyz[idx], mid=mid,
                              plane_origin=plane_origin, gradient=gradient))
    if not areas:
        raise ValueError('Landing_Area shapefiles are missing.')
    global_faf = load(qgis / 'FAF.shp')
    global_ids = {}
    for _, row in global_faf.iterrows():
        if int(row.id) in global_ids:
            raise ValueError('Duplicate global FAF id.')
        global_ids[int(row.id)] = row
    fafs = []
    pair_paths = sorted((qgis / 'FAF_List').glob('*.shp'))
    if not pair_paths:
        raise ValueError('No FAF_List pairs.')
    for pair, path in enumerate(pair_paths, 1):
        sources.append(path)
        data = load(path)
        if len(data) != 2:
            raise ValueError(f'Expected exactly two alternative FAF features: {path}')
        for _, row in data.iterrows():
            point = single_point(row.geometry)
            endpoint = np.array([point.x, point.y, float(row.alt)])
            if not np.isfinite(endpoint).all():
                raise ValueError(f'Invalid FAF alt/XY: {path}')
            if endpoint[2] > config.ring_z:
                raise ValueError(f'FAF {row.id} above ring; descending STAR impossible.')
            matches = [a for a in areas if np.linalg.norm(a['mid'][:2] - endpoint[:2]) < .01]
            if len(matches) != 1:
                raise ValueError(f'FAF {row.id}: expected one corresponding Landing_Area.')
            area = matches[0]
            if int(row.id) not in global_ids:
                raise ValueError(f'FAF {row.id} missing from master FAF.shp')
            master = global_ids[int(row.id)]
            mp = single_point(master.geometry)
            if math.dist((mp.x, mp.y), endpoint[:2]) > .01 or abs(float(master.alt) - endpoint[2]) > 1e-6:
                raise ValueError(f'FAF_List and FAF.shp differ for {row.id}')
            delta = area['apex'][:2] - endpoint[:2]
            fafs.append(dict(pair=pair, pair_name=path.stem, faf_id=int(row.id), xyz=endpoint,
                             heading=math.atan2(delta[1], delta[0]), area_fid=area['fid'],
                             source_kt=float(row.get('speed_kts', np.nan)),
                             geom_z=float(point.z) if point.has_z else None))
    if len({f['faf_id'] for f in fafs}) != len(fafs):
        raise ValueError('FAF IDs must be unique across pairs.')
    audit = dict(iaf_source_count=len(source), iaf_unique_count=len(xy),
                 generated_iaf_count=len(nodes), ring_center=center.tolist(),
                 ring_radius_m=ring_radius, ring_fit_residual_m=residual,
                 pair_count=len(pair_paths), faf_count=len(fafs), landing_count=len(areas))
    return nodes, headings, fafs, areas, audit, sources


def rsr(start, h0, end, h1, radius):
    """External tangent construction, returns arc/straight/arc lengths in m."""
    c0 = start + radius * np.array([math.sin(h0), -math.cos(h0)])
    c1 = end + radius * np.array([math.sin(h1), -math.cos(h1)])
    delta = c1 - c0
    distance = float(np.linalg.norm(delta))
    if distance < 1e-8:
        return None  # degenerate coincident circles, no unique CSC tangent
    tangent = math.atan2(delta[1], delta[0])
    a = (h0 - tangent) % TAU
    b = (tangent - h1) % TAU
    if a > TAU - 1e-10:
        a = 0.
    if b > TAU - 1e-10:
        b = 0.
    return np.array([radius * a, distance, radius * b])


def sample_rsr(start, h0, end_xyz, radius, lengths, config):
    """Include CSC joints, preserve exact arclength and analytic end headings."""
    step = min(config.step_m, math.sqrt(8 * radius * config.chord_error_m))
    position = np.array(start, dtype=float)
    heading = h0
    coords = [position.copy()]
    distances = [0.]
    cursor = 0.
    for segment, length in enumerate(lengths):
        if length <= 1e-9:
            continue
        ds = np.linspace(0, length, max(1, math.ceil(length / step)) + 1)[1:]
        if segment == 1:
            points = position + ds[:, None] * np.array([math.cos(heading), math.sin(heading)])
        else:
            angle = heading - ds / radius
            points = position + radius * np.column_stack((math.sin(heading) - np.sin(angle),
                                                          np.cos(angle) - math.cos(heading)))
        coords.extend(points)
        distances.extend(cursor + ds)
        position = points[-1]
        if segment != 1:
            heading -= length / radius
        cursor += length
    xy = np.asarray(coords)
    if np.linalg.norm(xy[-1] - end_xyz[:2]) > 1e-6:
        raise RuntimeError('RSR endpoint integration failed.')
    xy[-1] = end_xyz[:2]
    s = np.asarray(distances)
    z = config.ring_z + (end_xyz[2] - config.ring_z) * s / s[-1]
    z[-1] = end_xyz[2]
    return np.column_stack((xy, z)), s


def minimum_abs_linear(values):
    return 0. if min(values) <= 0 <= max(values) else float(min(abs(v) for v in values))


def landing_metrics(xyz, areas, details=False):
    """Intersect individual chords: exact Z interpolation on the exported line.

    Point-only boundary contacts have zero length and are allowed. Nonzero
    overlaps are evaluated even for the route's own selected Landing_Area.
    Plane-vs-route height difference is linear on each clipped segment, so
    endpoint values detect its exact minimum, including zero crossings.
    """
    import shapely
    line = LineString(xyz[:, :2])
    polys = [a['geometry'] for a in areas]
    hit_indices = np.flatnonzero(shapely.intersects(line, polys))
    overlaps, minimums, detail = [], [], []
    segments = None
    for ai in hit_indices:
        area = areas[ai]
        if line.intersection(area['geometry']).length <= 1e-5:
            continue
        if segments is None:
            segments = shapely.linestrings(np.stack((xyz[:-1, :2], xyz[1:, :2]), axis=1))
        indices = np.flatnonzero(shapely.intersects(segments, area['geometry']))
        intersections = shapely.intersection(segments[indices], area['geometry'])
        overlap, clearance, signed_min = 0., math.inf, math.inf
        for si, clipped in zip(indices, intersections):
            if clipped.length <= 1e-7:
                continue
            points = get_coordinates(clipped)
            start, end = xyz[si], xyz[si + 1]
            direction = end[:2] - start[:2]
            t = np.clip((points - start[:2]) @ direction / np.dot(direction, direction), 0, 1)
            route_z = start[2] + t * (end[2] - start[2])
            p0 = area['plane_origin']
            area_z = p0[2] + (points - p0[:2]) @ area['gradient']
            dz = route_z - area_z
            clearance = min(clearance, minimum_abs_linear(dz))
            signed_min = min(signed_min, float(dz.min()))
            overlap += clipped.length
        if overlap <= 1e-5:
            continue
        overlaps.append(overlap)
        minimums.append(clearance)
        if details:
            detail.append(dict(area_pair=area['pair'], area_fid=area['fid'],
                               overlap_m=overlap, min_abs_z=clearance, min_signed_z=signed_min))
    return len(overlaps), sum(overlaps), minimums, detail


def generate(nodes, headings, fafs, areas, config, output):
    v0 = config.start_kt * KNOT
    minimum_radius = v0 ** 2 / (G * math.tan(math.radians(config.max_bank)))
    radii = sorted(set(config.radii or (minimum_radius * 1.000001, 350., 500., 750.)))
    if any(r < minimum_radius - 1e-6 for r in radii):
        raise ValueError(f'Turn radius must be >= {minimum_radius:.6f} m.')
    gradient_limit = min(math.tan(math.radians(config.max_slope_deg)), config.max_ft_nm * .3048 / 1852)
    rows, rejected = [], dict(slope=0, degenerate=0)
    total = len(nodes) * len(fafs) * len(radii)
    started, done = time.monotonic(), 0
    for fi, faf in enumerate(fafs):
        for ni, (node, heading) in enumerate(zip(nodes, headings)):
            for radius in radii:
                done += 1
                lengths = rsr(node, heading, faf['xyz'][:2], faf['heading'], radius)
                if lengths is None:
                    rejected['degenerate'] += 1
                    continue
                length = float(lengths.sum())
                # Conservative exported-chord gradient bound (sagitta tolerance).
                step = min(config.step_m, math.sqrt(8 * radius * config.chord_error_m))
                chord_ratio = np.sinc(step / (2 * radius) / math.pi)
                gradient = (config.ring_z - faf['xyz'][2]) / length
                if gradient / chord_ratio > gradient_limit + 1e-12:
                    rejected['slope'] += 1
                    continue
                xyz, _ = sample_rsr(node, heading, faf['xyz'], radius, lengths, config)
                crossings, overlap, gaps, _ = landing_metrics(xyz, areas)
                bank = math.degrees(math.atan(v0 ** 2 / (G * radius)))
                clearance_cost = sum(math.exp(-gap / config.clearance_scale_m) for gap in gaps)
                base_cost = config.length_weight * length / 1000 + config.bank_weight * bank / config.max_bank
                score = (base_cost + config.crossing_weight * crossings + config.overlap_weight * overlap
                         + config.clearance_weight * clearance_cost)
                rows.append(dict(cand_id=len(rows) + 1, pair=faf['pair'], faf_index=fi,
                                 faf_id=faf['faf_id'], iaf_id=ni + 1, radius_m=radius,
                                 length_m=length, slope_deg=math.degrees(math.atan(gradient)),
                                 ft_nm=gradient * 1852 / .3048, bank_deg=bank,
                                 cross_n=crossings, overlap_m=overlap,
                                 min_gap_m=min(gaps) if gaps else np.nan,
                                 arc1_m=lengths[0], straight_m=lengths[1], arc2_m=lengths[2],
                                 base_cost=base_cost, score=score))
        print(f'Candidates FAF {faf["faf_id"]}: {done}/{total}, '
              f'accepted={len(rows)}, elapsed={time.monotonic()-started:.1f}s', flush=True)
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError(f'No gradient/bank-feasible candidates: {rejected}')
    frame.to_csv(output / 'all_candidate_metrics.csv', index=False, encoding='utf-8-sig')
    return frame, dict(total_combinations=total, hard_feasible=len(rows), rejected=rejected,
                       minimum_radius_m=minimum_radius, radii_m=radii,
                       gradient_limit=gradient_limit, generation_seconds=time.monotonic()-started)


def route_xyz(row, nodes, headings, fafs, config):
    ni = int(row.iaf_id) - 1
    faf = fafs[int(row.faf_index)]
    lengths = np.array([row.arc1_m, row.straight_m, row.arc2_m])
    return sample_rsr(nodes[ni], headings[ni], faf['xyz'], float(row.radius_m), lengths, config)


def clipped_parameter_box(a, b, limit):
    """Box (t,u) clipped by |za(t)-zb(u)| <= limit."""
    vertices = [np.array(p, dtype=float) for p in ((0, 0), (1, 0), (1, 1), (0, 1))]
    base = a[0, 2] - b[0, 2]
    linear = np.array([a[1, 2] - a[0, 2], -(b[1, 2] - b[0, 2])])
    for sign in (1, -1):
        new = []
        for p, q in zip(vertices, vertices[1:] + vertices[:1]):
            fp, fq = sign * (base + linear @ p) - limit, sign * (base + linear @ q) - limit
            if fp <= 0:
                new.append(p)
            if (fp <= 0) != (fq <= 0):
                new.append(p + (q - p) * fp / (fp - fq))
        vertices = new
        if not vertices:
            break
    return vertices


def separation_conflict(xyz_a, xyz_b, horizontal, vertical):
    """Continuous piecewise-linear corridor check; optional, disabled by default.

    Within a pair of segments minimize horizontal squared distance on the
    parameter box clipped by the vertical band. The minimum occurs at the
    unconstrained least-squares point or an edge of this convex polygon.
    """
    import shapely
    line_a, line_b = LineString(xyz_a[:, :2]), LineString(xyz_b[:, :2])
    if line_a.distance(line_b) >= horizontal:
        return False
    seg_a = shapely.linestrings(np.stack((xyz_a[:-1, :2], xyz_a[1:, :2]), axis=1))
    seg_b = shapely.linestrings(np.stack((xyz_b[:-1, :2], xyz_b[1:, :2]), axis=1))
    tree = shapely.STRtree(seg_b)
    indices = tree.query(seg_a, predicate='dwithin', distance=horizontal)
    for i, j in indices.T:
        a, b = xyz_a[i:i+2], xyz_b[j:j+2]
        poly = clipped_parameter_box(a, b, vertical - 1e-8)
        if not poly:
            continue
        origin = a[0, :2] - b[0, :2]
        matrix = np.column_stack((a[1, :2] - a[0, :2], -(b[1, :2] - b[0, :2])))
        unconstrained = np.linalg.lstsq(matrix, -origin, rcond=None)[0]
        if (np.all(unconstrained >= 0) and np.all(unconstrained <= 1)
                and abs(a[0, 2] - b[0, 2] + unconstrained @ np.array([
                    a[1, 2] - a[0, 2], -(b[1, 2] - b[0, 2])])) < vertical
                and np.linalg.norm(origin + matrix @ unconstrained) < horizontal):
            return True
        for p, q in zip(poly, poly[1:] + poly[:1]):
            start, vector = origin + matrix @ p, matrix @ (q - p)
            norm = np.dot(vector, vector)
            t = np.clip(-np.dot(start, vector) / norm, 0, 1) if norm > 0 else 0
            if np.linalg.norm(start + t * vector) < horizontal - 1e-8:
                return True
    return False


def select(frame, nodes, headings, fafs, config, strict):
    """Multistart discrete coordinate descent over ALL eligible candidates.

    Objective: sum unary scores + Gaussian penalty on neighboring IAF nodes.
    A pair contributes one route, never both alternate FAFs. Strict mode uses
    only zero-intersection candidates; impossible groups are explicitly absent.
    Optional separation can further reduce coverage; coverage is prioritized.
    """
    data = frame.loc[frame.cross_n == 0] if strict else frame
    groups = {int(pair): rows.sort_values(['score', 'cand_id']) for pair, rows in data.groupby('pair')}
    if not groups:
        return [], dict(coverage=0, score=None, local_optimum=True, multistarts=0)
    cache, conflicts = {}, {}

    def geometry(index):
        if index not in cache:
            cache[index] = route_xyz(frame.loc[index], nodes, headings, fafs, config)[0]
        return cache[index]

    def conflict(i, j):
        key = tuple(sorted((i, j)))
        if key not in conflicts:
            # Include total chord-error allowance so a curve between vertices
            # cannot evade a requested horizontal corridor separation.
            conflicts[key] = separation_conflict(geometry(i), geometry(j),
                config.horizontal_sep_m + 2 * config.chord_error_m, config.vertical_sep_m)
        return conflicts[key]

    def evaluate(chosen):
        cost = sum(float(frame.loc[i, 'score']) for i in chosen.values())
        ids = list(chosen.values())
        for k, i in enumerate(ids):
            for j in ids[k+1:]:
                d = np.linalg.norm(nodes[int(frame.loc[i, 'iaf_id'])-1] - nodes[int(frame.loc[j, 'iaf_id'])-1])
                cost += config.node_weight * math.exp(-(d/config.node_neighborhood_m)**2)
        return cost

    rng = np.random.default_rng(config.seed)
    best, best_key, converged = {}, (0, math.inf), False
    for run in range(config.multistarts):
        chosen = {}
        order = list(groups)
        if run:
            rng.shuffle(order)
        for pair in order:
            rows = groups[pair]
            if run and not config.separation:
                # Randomized initialization is restricted to the best 100
                # unary candidates; subsequent updates examine the full group.
                chosen[pair] = int(rng.choice(rows.index[:100]))
            else:
                for index in rows.index:
                    if not config.separation or all(not conflict(index, j) for j in chosen.values()):
                        chosen[pair] = index
                        break
        for iteration in range(50):
            changed = False
            for pair in order:
                rows = groups[pair]
                costs = rows.score.to_numpy(copy=True)
                starts = nodes[rows.iaf_id.to_numpy(dtype=int)-1]
                for other, index in chosen.items():
                    if other == pair:
                        continue
                    d = np.linalg.norm(starts - nodes[int(frame.loc[index, 'iaf_id'])-1], axis=1)
                    costs += config.node_weight * np.exp(-(d/config.node_neighborhood_m)**2)
                rank = np.argsort(costs, kind='stable')
                choice = None
                for k in rank:
                    index = int(rows.index[k])
                    if not config.separation or all(not conflict(index, j) for other, j in chosen.items() if other != pair):
                        choice = index
                        break
                if choice is not None and chosen.get(pair) != choice:
                    chosen[pair] = choice
                    changed = True
            if not changed:
                break
        key = (-len(chosen), evaluate(chosen))
        if key < best_key:
            best, best_key = chosen.copy(), key
            converged = not changed
    lower_bound = sum(float(groups[pair].score.min()) for pair in best)
    return [best[pair] for pair in sorted(best)], dict(
        coverage=len(best), score=best_key[1], local_optimum=converged,
        multistarts=config.multistarts, global_optimality_proven=False,
        unary_lower_bound=lower_bound,
        relative_bound_gap=(best_key[1]-lower_bound)/max(abs(best_key[1]), 1e-12))


def write_lines(frame, indices, path, nodes, headings, fafs, config):
    """SHP field names <=10 chars, XYZ verified after reopening."""
    records, lines = [], []
    for index in indices:
        row = frame.loc[index]
        xyz, _ = route_xyz(row, nodes, headings, fafs, config)
        lengths = np.linalg.norm(np.diff(xyz[:, :2], axis=0), axis=1)
        slopes = np.abs(np.diff(xyz[:, 2])) / lengths
        limit = min(math.tan(math.radians(config.max_slope_deg)), config.max_ft_nm * .3048 / 1852)
        if slopes.max() > limit + 1e-8:
            raise RuntimeError('Exported segment descent-gradient validation failed.')
        if float(row.bank_deg) > config.max_bank + 1e-7:
            raise RuntimeError('Exported bank validation failed.')
        records.append({key: row[key] for key in ['cand_id', 'pair', 'faf_id', 'iaf_id', 'radius_m',
            'length_m', 'slope_deg', 'ft_nm', 'bank_deg', 'cross_n', 'overlap_m', 'min_gap_m', 'score']})
        records[-1].update(z_start=config.ring_z, z_end=float(xyz[-1, 2]),
                           v_start=config.start_kt, v_end=config.end_kt, path_type='RSR')
        records[-1]['outside_m'] = float(LineString(xyz[:, :2]).difference(Polygon(nodes)).length)
        lines.append(LineString(xyz))
    if not lines:
        return
    result = gpd.GeoDataFrame(records, geometry=lines, crs=5186)
    result.to_file(path, driver='ESRI Shapefile', engine='pyogrio', encoding='UTF-8')
    saved = gpd.read_file(path)
    if len(saved) != len(lines) or saved.crs.to_epsg() != 5186 or not saved.has_z.all():
        raise RuntimeError(f'Z/CRS/count roundtrip failed: {path}')
    for actual, expected in zip(saved.geometry, lines):
        if not np.allclose(np.asarray(actual.coords), np.asarray(expected.coords), rtol=0, atol=1e-7):
            raise RuntimeError(f'XYZ roundtrip failed: {path}')


def visualize(frame, indices, nodes, headings, fafs, areas, output, config, suffix):
    colors = plt.get_cmap('tab10')
    fig, ax = plt.subplots(figsize=(12, 11))
    for area in areas:
        xy = np.asarray(area['geometry'].exterior.coords)
        ax.fill(xy[:, 0], xy[:, 1], color='tomato', alpha=.15)
        ax.plot(xy[:, 0], xy[:, 1], color='tomato', linewidth=.6)
    ring = np.vstack((nodes, nodes[0]))
    ax.plot(ring[:, 0], ring[:, 1], '--', color='gray',
            label=f'Inner ring: CW, {config.ring_z:g} m, {config.start_kt:g} kt')
    ax.scatter(nodes[:, 0], nodes[:, 1], s=.5, color='gray', alpha=.25)
    for faf in fafs:
        p = faf['xyz']
        ax.scatter(*p[:2], marker='x', color='black', s=20)
    fig3 = plt.figure(figsize=(13, 10))
    ax3 = fig3.add_subplot(111, projection='3d')
    for area in areas:
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        ax3.add_collection3d(Poly3DCollection([area['xyz']], facecolor='tomato', alpha=.12))
    ax3.plot(ring[:, 0], ring[:, 1], np.full(len(ring), config.ring_z), '--', color='gray')
    profiles, axes = plt.subplots(max(len(indices), 1), 1, figsize=(11, max(4, 2.5*len(indices))), squeeze=False)
    vertices = []
    for k, index in enumerate(indices):
        row = frame.loc[index]
        xyz, s = route_xyz(row, nodes, headings, fafs, config)
        color = colors((int(row.pair)-1) % 10)
        ax.plot(xyz[:, 0], xyz[:, 1], color=color, linewidth=1.5,
                label=f'Pair {int(row.pair)} / FAF {int(row.faf_id)} / {row.length_m:.0f} m')
        ax.scatter(*xyz[0, :2], color=color, marker='o', s=25)
        ax.scatter(*xyz[-1, :2], color=color, marker='*', s=100)
        # Direction arrows distinguish travel direction from geometric crossings.
        at = min(len(xyz)-2, max(0, len(xyz)//3))
        ax.annotate('', xy=xyz[at+1, :2], xytext=xyz[at, :2], arrowprops=dict(arrowstyle='->', color=color))
        ax3.plot(*xyz.T, color=color, label=f'Pair {int(row.pair)} / FAF {int(row.faf_id)}')
        axis = axes[k, 0]
        axis.plot(s, xyz[:, 2], color=color)
        axis.set(xlabel='Horizontal arclength (m)', ylabel='Z (m)',
                 title=f'Pair {int(row.pair)}: {row.slope_deg:.3f} deg, R={row.radius_m:.1f} m, '
                       f'bank <= {row.bank_deg:.2f} deg, crossings={int(row.cross_n)}')
        axis.grid(alpha=.25)
        v0, v1 = config.start_kt*KNOT, config.end_kt*KNOT
        speed = np.sqrt(v0*v0 + (v1*v1-v0*v0)*s/s[-1])
        # TAS is along the descending 3D path; s is horizontal arclength.
        distance3d = s * math.hypot(1., (config.ring_z-xyz[-1,2])/s[-1])
        elapsed = 2*distance3d/(v0+speed)
        for j, (point, distance, velocity, seconds) in enumerate(zip(xyz, s, speed, elapsed)):
            vertices.append(dict(pair=int(row.pair), faf_id=int(row.faf_id), vertex=j,
                                 x=point[0], y=point[1], z=point[2], s_m=distance,
                                 s3d_m=distance3d[j], speed_kt=velocity/KNOT, time_s=seconds,
                                 acceleration=(v1*v1-v0*v0)/(2*distance3d[-1])))
    ax.set(aspect='equal', xlabel='EPSG:5186 X (m)', ylabel='EPSG:5186 Y (m)',
           title=f'Inner IAF -> FAF / {suffix} / separation disabled' if not config.separation else suffix)
    ax.ticklabel_format(style='plain', useOffset=False)
    ax.grid(alpha=.2)
    ax.legend(loc='upper left', fontsize=8)
    fig.tight_layout()
    fig.savefig(output / f'map_{suffix}.png', dpi=180)
    ax3.set(xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)', title=f'RSR routes: {suffix}')
    ax3.legend(fontsize=8)
    fig3.savefig(output / f'view3d_{suffix}.png', dpi=160)
    profiles.tight_layout()
    profiles.savefig(output / f'profiles_{suffix}.png', dpi=160)
    pd.DataFrame(vertices).to_csv(output / f'vertices_{suffix}.csv', index=False, encoding='utf-8-sig')
    plt.close('all')


def self_check():
    """Geometry and continuous-clearance regression checks, no source writes."""
    cfg = Config()
    rng = np.random.default_rng(10)
    for _ in range(100):
        start, end = rng.normal(size=(2, 2))*2000
        h0, h1 = rng.uniform(-math.pi, math.pi, 2)
        r = rng.uniform(230, 750)
        lengths = rsr(start, h0, end, h1, r)
        xyz, s = sample_rsr(start, h0, np.r_[end, 250], r, lengths, cfg)
        assert np.linalg.norm(xyz[-1, :2]-end) < 1e-7
        assert abs(((h0-(lengths[0]+lengths[2])/r-h1+math.pi) % TAU)-math.pi) < 1e-9
        assert np.all(np.diff(s) > 0) and xyz[0, 2] == 350 and xyz[-1, 2] == 250
    a = np.array([[0,0,350],[10,0,350]], dtype=float)
    b = np.array([[5,-5,300],[5,5,300]], dtype=float)
    assert separation_conflict(a,b,1,70)
    b[:,2] = 270
    assert not separation_conflict(a,b,1,70)
    b = np.array([[0,2,350],[10,2,350]], dtype=float)
    assert not separation_conflict(a,b,1,70)
    assert separation_conflict(a,b,3,70)
    area = dict(pair='test', fid=1, geometry=Polygon([(0,-1),(10,-1),(0,9)]),
                plane_origin=np.array([0.,0.,0.]), gradient=np.array([0.,0.]))
    metrics = landing_metrics(np.array([[-1.,0.,10.],[11.,0.,10.]]),[area])
    assert metrics[0] == 1 and abs(metrics[1]-9) < 1e-8 and metrics[2] == [10.]
    assert landing_metrics(np.array([[-1.,0.,-6.],[11.,0.,6.]]),[area])[2] == [0.]
    assert landing_metrics(np.array([[-1.,0.,10.],[0.,0.,10.]]),[area])[0] == 0
    print('Self-check passed: 100 RSR endpoints/headings + separation + Landing_Area clipping/clearance.')


def legacy_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qgis', type=Path, default=BASE/'QGIS')
    parser.add_argument('--output', type=Path, default=BASE/'결과물/optimization_path')
    parser.add_argument('--iaf-count', type=int, default=1000)
    parser.add_argument('--ring-z', type=float, default=350.)
    parser.add_argument('--radii', type=float, nargs='+')
    parser.add_argument('--max-bank', type=float, default=30.)
    parser.add_argument('--max-slope-deg', type=float, default=2.36)
    parser.add_argument('--start-kt', type=float, default=70.)
    parser.add_argument('--end-kt', type=float, default=60.)
    parser.add_argument('--enable-separation', action='store_true')
    parser.add_argument('--horizontal-sep', type=float, default=100.)
    parser.add_argument('--vertical-sep', type=float, default=70.)
    parser.add_argument('--self-check', action='store_true')
    parser.add_argument('--reuse-metrics', action='store_true', help='Reuse metrics only when config and input hashes match')
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return
    config = Config(iaf_count=args.iaf_count, ring_z=args.ring_z,
                    radii=tuple(args.radii or ()), max_bank=args.max_bank,
                    max_slope_deg=args.max_slope_deg, start_kt=args.start_kt, end_kt=args.end_kt,
                    separation=args.enable_separation, horizontal_sep_m=args.horizontal_sep,
                    vertical_sep_m=args.vertical_sep)
    values = [config.ring_z, config.start_kt, config.end_kt, config.max_bank,
              config.max_slope_deg, config.horizontal_sep_m, config.vertical_sep_m]
    if not all(math.isfinite(v) for v in values) or config.iaf_count < 2 or not (
        0 < config.end_kt <= config.start_kt and 0 < config.max_bank < 90
        and 0 < config.max_slope_deg < 90 and config.horizontal_sep_m > 0 and config.vertical_sep_m > 0):
        raise ValueError('Invalid speed/bank/slope/count/separation settings.')
    if any(not math.isfinite(r) or r <= 0 for r in config.radii):
        raise ValueError('Radii must be finite and positive.')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    nodes, headings, fafs, areas, audit, sources = prepare(args.qgis, config)
    hashes = fingerprint(sources)
    signature = json.loads(json.dumps(dict(config=asdict(config), inputs=hashes)))
    cache_path = output/'candidate_run.json'
    if args.reuse_metrics:
        old = json.loads(cache_path.read_text(encoding='utf-8'))
        if old['signature'] != signature:
            raise ValueError('Cached input/config differs. Rerun without --reuse-metrics.')
        frame = pd.read_csv(output/'all_candidate_metrics.csv')
        generation = old['generation']
    else:
        frame, generation = generate(nodes, headings, fafs, areas, config, output)
        cache_path.write_text(json.dumps(dict(signature=signature, generation=generation), indent=2,
                                        ensure_ascii=False), encoding='utf-8')
    summary = dict(config=asdict(config), input_audit=audit, generation=generation,
                   separation_applied=config.separation, input_sha256=hashes,
                   speed_model=f'constant tangential deceleration along 3D arclength, {config.start_kt}->{config.end_kt} kt TAS, still air',
                   z_model='ring_z -> FAF alt, same vertical datum assumed; no datum conversion',
                   objective='crossing_weight*n + overlap_weight*m + clearance_weight*sum(exp(-abs_gap/scale)) + length_weight*km + bank_weight*bank/max_bank + node_weight*sum(exp(-(node_distance/neighborhood)**2))',
                   optimality='multistart coordinate descent over enumerated IAF/FAF/radius grid; not global proof',
                   chord_error_m=config.chord_error_m)
    gpd.GeoDataFrame(dict(iaf_id=np.arange(1,len(nodes)+1), heading=np.degrees(headings)%360),
        geometry=[Point(x,y,config.ring_z) for x,y in nodes], crs=5186).to_file(
            output/'iaf_candidates_z.shp', engine='pyogrio', encoding='UTF-8')
    shortlists, counts = [], []
    for pair, rows in frame.groupby('pair'):
        strict_rows = rows[rows.cross_n == 0]
        counts.append(dict(pair=int(pair), feasible=len(rows), strict=len(strict_rows)))
        shortlists.extend(rows.nsmallest(config.shortlist_per_pair, 'score').index.tolist())
        shortlists.extend(strict_rows.nsmallest(5, 'score').index.tolist())
    selected_scenarios = {}
    for strict, suffix in ((True,'strict'), (False,'penalized')):
        indices, selection = select(frame, nodes, headings, fafs, config, strict)
        selected_scenarios[suffix] = indices
        summary[suffix] = selection
        summary[suffix]['missing_pairs'] = sorted(set(range(1,audit['pair_count']+1)) - set(frame.loc[indices,'pair'].astype(int)))
        summary[suffix]['length_m'] = float(frame.loc[indices,'length_m'].sum())
        summary[suffix]['crossings'] = int(frame.loc[indices,'cross_n'].sum())
        frame.loc[indices].to_csv(output/f'selected_{suffix}.csv', index=False, encoding='utf-8-sig')
        write_lines(frame, indices, output/f'routes_{suffix}_z.shp', nodes, headings, fafs, config)
        visualize(frame, indices, nodes, headings, fafs, areas, output, config, suffix)
        details = []
        for index in indices:
            xyz, _ = route_xyz(frame.loc[index], nodes, headings, fafs, config)
            n, length, gaps, items = landing_metrics(xyz, areas, details=True)
            if n != int(frame.loc[index,'cross_n']) or not math.isclose(length,frame.loc[index,'overlap_m'],abs_tol=1e-5):
                raise RuntimeError('Final Landing_Area audit differs from candidate evaluation.')
            details.extend(dict(pair=int(frame.loc[index,'pair']),faf_id=int(frame.loc[index,'faf_id']),**d) for d in items)
        pd.DataFrame(details, columns=['pair','faf_id','area_pair','area_fid','overlap_m','min_abs_z','min_signed_z']).to_csv(
            output/f'landing_crossings_{suffix}.csv', index=False, encoding='utf-8-sig')
        print(f'{suffix}: {selection}', flush=True)
    primary = 'strict' if summary['strict']['coverage'] == audit['pair_count'] else 'penalized'
    summary['primary_scenario'] = primary
    summary['ring_boundary_is_hard_constraint'] = False
    summary['candidate_counts'] = counts
    summary['shortest_gradient_feasible_total_m'] = float(frame.groupby('pair').length_m.min().sum())
    for suffix in selected_scenarios:
        summary[suffix]['length_ratio_to_shortest_independent'] = (
            summary[suffix]['length_m']/summary['shortest_gradient_feasible_total_m'])
    best = selected_scenarios[primary]
    summary['primary_outside_ring_m'] = {}
    boundary = Polygon(nodes)
    for index in best:
        xyz, _ = route_xyz(frame.loc[index], nodes, headings, fafs, config)
        summary['primary_outside_ring_m'][str(int(frame.loc[index,'pair']))] = float(
            LineString(xyz[:,:2]).difference(boundary).length)
    write_lines(frame, best, output/'selected_routes_z.shp', nodes, headings, fafs, config)
    final_fafs = [fafs[int(frame.loc[i,'faf_index'])] for i in best]
    gpd.GeoDataFrame([dict(pair=f['pair'],faf_id=f['faf_id'],alt=f['xyz'][2],speed_kt=config.end_kt,
                           source_kt=f['source_kt']) for f in final_fafs],
                     geometry=[Point(*f['xyz']) for f in final_fafs],crs=5186).to_file(
                         output/'selected_faf_z.shp',engine='pyogrio',encoding='UTF-8')
    indices = sorted(set(shortlists + [i for ids in selected_scenarios.values() for i in ids]))
    write_lines(frame, indices, output/'shortlist_routes_z.shp', nodes, headings, fafs, config)
    if fingerprint(sources) != hashes:
        raise RuntimeError('Input source hashes changed during execution.')
    summary['original_inputs_unchanged'] = True
    (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    report(output, summary, frame.loc[best], fafs)
    print(f'Complete: {output}', flush=True)


def report(output, summary, selected, fafs):
    records = selected[['pair','faf_id','iaf_id','radius_m','length_m','slope_deg','bank_deg','cross_n','overlap_m']]
    table = records.to_html(index=False,float_format=lambda v:f'{v:.3f}')
    text = f'''<!doctype html><html lang="ko"><meta charset="utf-8"><title>Inner IAF–FAF 항로 분석</title>
<style>body{{font-family:Malgun Gothic,sans-serif;max-width:1200px;margin:36px auto;line-height:1.7}}img{{max-width:100%}}table{{border-collapse:collapse}}td,th{{padding:6px;border:1px solid #ccc}}.notice{{background:#fff1d7;padding:16px}}</style>
<h1>Inner IAF → FAF 항로 분석</h1>
<p>1000 IAF × 16 FAF × {len(summary['generation']['radii_m'])} 반경: 총 {summary['generation']['total_combinations']:,}개 조합.
내곽 {summary['config']['ring_z']:g}m, {summary['config']['start_kt']:g}→{summary['config']['end_kt']:g}kt, 시계방향 RSR. FAF 높이는 alt 필드 사용.</p>
<p class="notice">수평·수직 분리 조건 활성화: {summary['separation_applied']}. 비활성화된 경우 서로 교차하는 항로가 포함될 수 있습니다.
선회반경은 무풍·TAS·협조선회 모델이며 기체별 승인 성능이나 운항 안전성 인증을 뜻하지 않습니다.
건물·지형·기상·롤 전이·출발 시각은 분석 범위에 없습니다. 고도 기준은 입력 자료와 링 높이가 같다고 가정합니다.</p>
<p>링은 IAF 출발 위치와 접선 방향을 정하는 데 사용하며, 링 내부에 항로를 가두는 경계 제약은 적용하지 않았습니다.
따라서 일부 항로는 링 외부로 우회합니다. SHP의 outside_m 필드에 링 밖 구간 길이를 기록했습니다.</p>
<p>하강 제한: {summary['config']['max_slope_deg']:g}°와 {summary['config']['max_ft_nm']:g}ft/NM 중 엄격한 값.
뱅크 ≤{summary['config']['max_bank']:g}°. 최소반경 {summary['generation']['minimum_radius_m']:.2f}m.
각 쌍에서 하나의 FAF만 선택합니다. 원본 SHP는 변경하지 않습니다.</p>
<p>IAF→FAF 속도는 고정된 시작·끝 속도 사이에서 3D 경로를 따라 등감속하는 모델로 설정했습니다.
모든 후보에 같은 속도 조건을 적용하므로 별도 속도 선택 최적화는 하지 않습니다.
FAF 이후의 −1m/s² 감속과 35kt 역천이는 이번 분석 범위에 포함하지 않습니다.</p>
<h2>선정 결과</h2><p>기본 결과: {summary['primary_scenario']}. 엄격 회피 {summary['strict']['coverage']}/8쌍,
관통 페널티 적용 {summary['penalized']['coverage']}/8쌍.</p>{table}
<p>Landing_Area 회피와 노드 분산을 무시하고 하강·뱅크 조건만 충족하는 쌍별 최단 항로 합:
{summary['shortest_gradient_feasible_total_m']/1000:.3f}km.
엄격 회피 선정 합 {summary['strict']['length_m']/1000:.3f}km,
관통 페널티 선정 합 {summary['penalized']['length_m']/1000:.3f}km.
연결된 쌍 수가 같을 때 이 값을 비교해 회피로 늘어난 길이를 판단할 수 있습니다.</p>
<p>관통은 Polygon XY 내부와 겹치는 선 길이로 측정하며 끝점 단독 접촉은 허용합니다.
원호는 최대 {summary['chord_error_m']:g}m 현 오차의 폴리라인으로 근사하며, 관통·높이 차는 저장한 폴리라인 기준으로 검사합니다.
모든 16개 Landing_Area를 회피 대상으로 사용합니다(미선택 반대편 FAF의 영역도 포함).
페널티는 관통 영역 수, 중첩 길이, 표면과의 최소 절대 높이 차, 항로 길이, 뱅크각,
선정 IAF 사이의 근접성으로 구성됩니다. 높이 차가 클수록 관통 페널티가 감소합니다.</p>
<p>모든 개별 조합을 평가한 뒤 다중 시작 좌표하강법으로 조합을 선정했습니다.
설정한 반경과 1000개 점에 대한 탐색이며 연속 공간 또는 전체 조합의 전역 최적성은 증명하지 않았습니다.
strict와 penalized 결과를 모두 제공하며, 비효율 판단에 의한 자동 완화는 하지 않습니다.</p>
<h2>엄격 회피</h2><img src="map_strict.png"><h2>관통 페널티 대안</h2><img src="map_penalized.png">
<h2>기본 결과 3D</h2><img src="view3d_{summary['primary_scenario']}.png">
<h2>종단 고도</h2><img src="profiles_{summary['primary_scenario']}.png">
<h2>파일</h2><ul><li>selected_routes_z.shp: 기본 선정 항로 (LineString Z)</li>
<li>routes_strict_z.shp / routes_penalized_z.shp: 시나리오별 선정 항로</li>
<li>shortlist_routes_z.shp: 쌍별 상위 25개 및 엄격 후보·선정 경로</li>
<li>iaf_candidates_z.shp / selected_faf_z.shp: 1000개 IAF 및 선정 FAF</li>
<li>all_candidate_metrics.csv: 모든 하강·뱅크 조건 통과 후보의 지표</li>
<li>vertices_*.csv: 꼭짓점 XYZ·속도·경과시간·가속도</li>
<li>landing_crossings_*.csv: 영역별 중첩 길이 및 최소 수직 높이 차</li>
<li>summary.json: 설정·후보 수·점수·원본 해시·최적성 한계</li></ul>
<h2>설정 변경</h2><pre>python -B 수정본_20261103/코드/optimization_path.py
# 별도 결과 폴더에 분리 조건을 켜는 예시(기본 분석에는 적용 안 함)
python -B 수정본_20261103/코드/optimization_path.py --enable-separation --horizontal-sep 100 --vertical-sep 70 --output 수정본_20261103/결과물/separation_test</pre>
<p>분리 옵션은 두 항로의 임의 위치 쌍에서 수평·수직이 동시에 지정값 미만이면 배제하는 보수적 공간 회랑 검사입니다.
Outer→IF와 Inner Ring의 이격 검사는 이번 범위에 없으며 해당 항로 추가 시 따로 구현해야 합니다.</p>
<p>근거: <a href="https://www.faa.gov/sites/faa.gov/files/pilots/pilot_handbook.pdf">FAA PHAK</a>,
<a href="https://ntrs.nasa.gov/citations/20240005129">NASA UAM 기체 설계 연구</a>,
<a href="https://ompl.kavrakilab.org/classompl_1_1base_1_1DubinsStateSpace.html">OMPL Dubins 모델</a>.</p></html>'''
    (output/'report.html').write_text(text,encoding='utf-8')
    (output/'README.md').write_text(
        '# Inner IAF–FAF 분석\n\nreport.html을 열어 지도와 종단면을 확인하세요.\n\n'
        'QGIS 기본 항로: selected_routes_z.shp (EPSG:5186, Z 있음).\n'
        f'수평·수직 분리 조건 활성화: {summary["separation_applied"]}.\n'
        '엄격 회피와 관통 페널티 대안은 별도 SHP로 제공됩니다.\n'
        '전역 최적성·기체별 운항 승인·실제 충돌 안전성을 보장하는 결과가 아닙니다.\n',encoding='utf-8')


@dataclass
class TrafficConfig:
    """Plan-view spacing is a soft objective, unrelated to Outer/IF Z spacing."""
    spacing_target_m: float = 100.
    route_cross_weight: float = 1_000_000.
    pair_gap_weight: float = 200.
    minimum_gap_weight: float = 5000.
    candidate_limit: int = 36
    angular_sectors: int = 24
    multistarts: int = 80
    max_sweeps: int = 35
    seed: int = 20261104
    outside_min_m: float = 1.
    max_arc_deg: float = 180.
    max_total_turn_deg: float = 180.


def traffic_parameters(row, nodes, headings, fafs, config):
    ni, fi = int(row.iaf_id)-1, int(row.faf_index)
    faf = fafs[fi]
    if row.operation == 'STAR':
        return nodes[ni], headings[ni], faf['xyz'], config
    # New clockwise SID, NOT a reversed clockwise STAR (which would be LSL).
    return (faf['xyz'][:2], faf['heading']+math.pi,
            np.r_[nodes[ni],config.ring_z], replace(config, ring_z=float(faf['xyz'][2])))


def traffic_xyz(row, nodes, headings, fafs, config):
    start, heading, end, sampling = traffic_parameters(row,nodes,headings,fafs,config)
    return sample_rsr(start,heading,end,float(row.radius_m),
                      np.array([row.arc1_m,row.straight_m,row.arc2_m]),sampling)


def maximum_ring_radius(start, heading, radius, lengths, ring_center):
    """Exact radial maximum of RSR arcs and straight, not a vertex-only test."""
    position = np.array(start,dtype=float)
    maximum = float(np.linalg.norm(position-ring_center))
    for segment, length in enumerate(lengths):
        if segment == 1:
            position = position + length*np.array([math.cos(heading),math.sin(heading)])
        else:
            center = position + radius*np.array([math.sin(heading),-math.cos(heading)])
            radial = heading+math.pi/2
            target = math.atan2(*(center-ring_center)[::-1])
            if (radial-target) % TAU <= length/radius+1e-12:
                maximum = max(maximum,float(np.linalg.norm(center-ring_center)+radius))
            angle = heading-length/radius
            position = position+radius*np.array([math.sin(heading)-math.sin(angle),
                                                math.cos(angle)-math.cos(heading)])
            heading = angle
        maximum = max(maximum,float(np.linalg.norm(position-ring_center)))
    return maximum


def traffic_generate(nodes, headings, fafs, areas, config, audit, hashes, output):
    """Reuse validated STAR inventory; calculate all SID RSR combinations anew."""
    signature = json.loads(json.dumps(dict(version=2,config=asdict(config),inputs=hashes)))
    cache = output/'traffic_candidate_run.json'
    csv = output/'all_star_sid_candidates.csv'
    for cache_folder in dict.fromkeys((output,BASE/'결과물/optimization_path_star_sid')):
        candidate_meta = cache_folder/'traffic_candidate_run.json'
        candidate_csv = cache_folder/'all_star_sid_candidates.csv'
        if candidate_meta.exists() and candidate_csv.exists():
            saved = json.loads(candidate_meta.read_text(encoding='utf-8'))
            if saved['signature'] == signature:
                print('Reusing validated STAR/SID candidate inventory.',flush=True)
                if cache_folder!=output:
                    shutil.copy2(candidate_csv,csv)
                    shutil.copy2(candidate_meta,cache)
                return pd.read_csv(csv),saved['generation']
    old = BASE/'결과물/optimization_path'
    old_meta = old/'candidate_run.json'
    star = None
    if old_meta.exists():
        saved = json.loads(old_meta.read_text(encoding='utf-8'))
        expected = json.loads(json.dumps(dict(config=asdict(config),inputs=hashes)))
        if saved['signature'] == expected:
            star = pd.read_csv(old/'all_candidate_metrics.csv')
            print(f'Reusing {len(star)} validated STAR candidates.',flush=True)
    if star is None:
        star, _ = generate(nodes,headings,fafs,areas,config,output)
    star['operation'] = 'STAR'
    center = np.array(audit['ring_center'])
    radial = []
    for row in star.itertuples():
        radial.append(maximum_ring_radius(nodes[int(row.iaf_id)-1],headings[int(row.iaf_id)-1],
            row.radius_m,[row.arc1_m,row.straight_m,row.arc2_m],center))
    star['max_rad_m'] = radial
    minimum_radius = max(config.start_kt,config.end_kt)**2*KNOT**2/(G*math.tan(math.radians(config.max_bank)))
    radii = sorted(set(config.radii or (minimum_radius*1.000001,350.,500.,750.)))
    if any(r < minimum_radius-1e-6 for r in radii):
        raise ValueError('Turn radius below configured bank/speed limit.')
    limit = min(math.tan(math.radians(config.max_slope_deg)),config.max_ft_nm*.3048/1852)
    records, rejected = [],0
    started = time.monotonic()
    for fi,faf in enumerate(fafs):
        for ni,(node,heading) in enumerate(zip(nodes,headings)):
            for radius in radii:
                lengths = rsr(faf['xyz'][:2],faf['heading']+math.pi,node,heading,radius)
                if lengths is None:
                    continue
                length = float(lengths.sum())
                gradient = (config.ring_z-faf['xyz'][2])/length
                step = min(config.step_m,math.sqrt(8*radius*config.chord_error_m))
                if gradient/np.sinc(step/(2*radius)/math.pi) > limit+1e-12:
                    rejected += 1
                    continue
                row = dict(operation='SID',pair=faf['pair'],faf_index=fi,faf_id=faf['faf_id'],
                           iaf_id=ni+1,radius_m=radius,length_m=length,
                           slope_deg=math.degrees(math.atan(gradient)),ft_nm=gradient*1852/.3048,
                           arc1_m=lengths[0],straight_m=lengths[1],arc2_m=lengths[2])
                xyz,_ = traffic_xyz(pd.Series(row),nodes,headings,fafs,config)
                crossings,overlap,gaps,_ = landing_metrics(xyz,areas)
                bank = math.degrees(math.atan((config.start_kt*KNOT)**2/(G*radius)))
                base = config.length_weight*length/1000+config.bank_weight*bank/config.max_bank
                row.update(bank_deg=bank,cross_n=crossings,overlap_m=overlap,
                           min_gap_m=min(gaps) if gaps else np.nan,base_cost=base,
                           score=base+config.crossing_weight*crossings+config.overlap_weight*overlap
                           +config.clearance_weight*sum(math.exp(-g/config.clearance_scale_m) for g in gaps),
                           max_rad_m=maximum_ring_radius(faf['xyz'][:2],faf['heading']+math.pi,
                                                       radius,lengths,center))
                records.append(row)
        print(f'SID FAF {faf["faf_id"]}: accepted={len(records)}, elapsed={time.monotonic()-started:.1f}s',flush=True)
    frame = pd.concat([star,pd.DataFrame(records)],ignore_index=True)
    frame['cand_id'] = np.arange(1,len(frame)+1)
    frame['inside'] = frame.max_rad_m <= audit['ring_radius_m']+1e-7
    frame['detour'] = frame.max_rad_m > audit['ring_radius_m']+1.
    frame.to_csv(csv,index=False,encoding='utf-8-sig')
    generation = dict(star_feasible=len(star),sid_feasible=len(records),sid_slope_rejected=rejected,
                      total_combinations=2*len(nodes)*len(fafs)*len(radii),radii_m=radii,
                      sid_generation_seconds=time.monotonic()-started,
                      minimum_radius_m=minimum_radius)
    cache.write_text(json.dumps(dict(signature=signature,generation=generation),ensure_ascii=False,indent=2),encoding='utf-8')
    return frame,generation


def review_candidate_shapes(frame,nodes,headings,fafs,config,traffic,output):
    """Explicit shape screening; not an FAA/NASA operational approval rule.

    Cap total clockwise change to avoid circle-like detours used merely to
    meet a gradient. Never change altitude/slope/speed/direction to fill a pair.
    Self-intersection is evaluated on the high-resolution exported XY line.
    """
    frame = frame.copy()
    frame['first_deg'] = np.degrees(frame.arc1_m/frame.radius_m)
    frame['last_deg'] = np.degrees(frame.arc2_m/frame.radius_m)
    frame['turn_deg'] = frame.first_deg+frame.last_deg
    frame['shape_ok'] = ((frame.first_deg<=traffic.max_arc_deg+1e-8)
                         &(frame.last_deg<=traffic.max_arc_deg+1e-8)
                         &(frame.turn_deg<=traffic.max_total_turn_deg+1e-8))
    frame['shape_reason'] = np.where(frame.shape_ok,'pass','excessive_turn')
    for index,row in frame[frame.shape_ok].iterrows():
        xyz,_ = traffic_xyz(row,nodes,headings,fafs,config)
        if not LineString(xyz[:,:2]).is_simple:
            frame.loc[index,['shape_ok','shape_reason']] = [False,'self_intersection']
    frame.to_csv(output/'review_candidate_metrics.csv',index=False,encoding='utf-8-sig')
    stats = dict(total_gradient_bank_feasible=len(frame),
                 shape_feasible=int(frame.shape_ok.sum()),
                 excessive_turn_rejected=int((frame.shape_reason=='excessive_turn').sum()),
                 self_intersection_rejected=int((frame.shape_reason=='self_intersection').sum()),
                 max_arc_deg=traffic.max_arc_deg,max_total_turn_deg=traffic.max_total_turn_deg,
                 angle_limit_is_model_assumption_not_regulatory_standard=True)
    print(f'Geometry review: {stats}',flush=True)
    return frame,stats


def traffic_pool(frame, scenario, config, traffic, audit):
    """Spatially diverse shortlist per operation/FAF; all candidates remain in CSV."""
    eligible = frame[frame.shape_ok] if 'shape_ok' in frame else frame
    if scenario == 'inside_only':
        eligible = eligible[eligible.inside]
    elif scenario == 'outside_detour':
        eligible = eligible[eligible.max_rad_m>audit['ring_radius_m']+traffic.outside_min_m]
    pairs = sorted(frame.pair.unique().astype(int))

    def feasible_orientations(data):
        orientations = {}
        for pair in pairs:
            ids = sorted(frame.loc[frame.pair==pair,'faf_id'].unique().astype(int))
            options = []
            for a,b in (ids,ids[::-1]):
                star = data[(data.pair==pair)&(data.operation=='STAR')&(data.faf_id==a)]
                sid = data[(data.pair==pair)&(data.operation=='SID')&(data.faf_id==b)]
                if len(star) and len(sid):
                    options.append((a,b))
            orientations[pair] = options
        return orientations

    strict = eligible[eligible.cross_n==0]
    options = feasible_orientations(strict)
    use_strict = all(options.values())
    strict_pairs = [p for p,v in options.items() if v]
    relaxed_pairs = [p for p,v in options.items() if not v]
    if use_strict:
        eligible = strict
    else:
        # Preserve the original hard Landing_Area rule wherever an opposite
        # STAR/SID pair can meet it. Do NOT globally relax all buildings just
        # because a different building has unavoidable overlap.
        eligible = pd.concat([strict[strict.pair.isin(strict_pairs)],
                              eligible[eligible.pair.isin(relaxed_pairs)]]).sort_index()
        options = feasible_orientations(eligible)
    indices = []
    for _,rows in eligible.groupby(['operation','faf_id']):
        indices.extend(rows.nsmallest(traffic.candidate_limit,'score').index)
        # Avoid a shortlist made entirely of neighboring ring nodes.
        sectors = ((rows.iaf_id-1)*traffic.angular_sectors//config.iaf_count).astype(int)
        for _,sector in rows.groupby(sectors):
            indices.extend(sector.nsmallest(2,'score').index)
    pool = eligible.loc[sorted(set(indices))].copy().reset_index(drop=True)
    return pool,options,dict(landing_mode='strict' if use_strict else 'pairwise_strict_then_penalized',
                              eligible_candidates=len(eligible),shortlisted_candidates=len(pool),
                              strict_pairs=[int(p) for p in strict_pairs],
                              penalized_pairs=[int(p) for p in relaxed_pairs],
                              strict_missing_pairs=[int(p) for p,v in feasible_orientations(strict).items() if not v],
                              missing_pairs=[int(p) for p,v in options.items() if not v])


def traffic_pair_matrices(pool,nodes,headings,fafs,config,traffic,cache_path=None,cache_key=None):
    """Exact exported-polyline distances/intersections, vectorized by GEOS.

    Simplify only redundant collinear points with 1e-8 m tolerance. Distance
    matrix is symmetric and includes endpoints. No Z separation is applied.
    """
    import shapely
    from concurrent.futures import ThreadPoolExecutor, as_completed
    if cache_path is not None and cache_path.exists():
        with np.load(cache_path,allow_pickle=False) as saved:
            if saved['key'].item()==cache_key:
                distances,crossing = saved['distance'].copy(),saved['crossing'].copy()
                gap = np.maximum(0.,1.-distances/traffic.spacing_target_m)**2
                np.fill_diagonal(gap,0.)
                print('Reusing validated route-spacing matrix.',flush=True)
                return None,distances,crossing,gap
    lines = np.empty(len(pool),dtype=object)
    for i,row in pool.iterrows():
        xyz,_ = traffic_xyz(row,nodes,headings,fafs,config)
        lines[i] = LineString(xyz[:,:2]).simplify(1e-8,preserve_topology=False)
    n = len(pool)
    distances = np.zeros((n,n),dtype=float)
    crossing = np.zeros((n,n),dtype=bool)
    def block(first):
        result = []
        for i in range(first,min(first+32,n)):
            # GEOS calls release Python's GIL. Workers read immutable geometry
            # objects and return arrays; only the main thread mutates matrices.
            hits = shapely.intersects(lines[i],lines[i+1:])
            values = np.zeros(len(hits))
            if (~hits).any():
                values[~hits] = shapely.distance(lines[i],lines[i+1:][~hits])
            result.append((i,values,hits))
        return result
    completed,next_update = 0,250
    with ThreadPoolExecutor(max_workers=4) as executor:
        pending = [executor.submit(block,i) for i in range(0,n,32)]
        for future in as_completed(pending):
            values = future.result()
            for i,distance,hits in values:
                distances[i,i+1:],crossing[i,i+1:] = distance,hits
            completed += len(values)
            if completed>=next_update:
                print(f'Route spacing matrix: {completed}/{n}',flush=True)
                next_update += 250
    distances += distances.T
    crossing |= crossing.T
    # Smooth deficit, with zero penalty at/above the requested soft target.
    gap = np.maximum(0.,1.-distances/traffic.spacing_target_m)**2
    np.fill_diagonal(gap,0.)
    if cache_path is not None:
        np.savez_compressed(cache_path,key=cache_key,distance=distances,crossing=crossing)
    return lines,distances,crossing,gap


def traffic_select(pool,options,distances,crossings,gaps,nodes,traffic,config,
                   initial_indices=None,forbid_crossings=False,minimum_spacing_m=0.,forbidden_pairs=None):
    """Joint opposite-FAF assignment and route choice; multistart block descent.

    Each block is a building pair: STAR uses one FAF, SID the other. Both
    opposite-FAF assignments and ALL STAR x SID candidate combinations within
    the block are evaluated exactly against other selected routes. This is
    still a shortlist block-descent heuristic, not a global optimizer.
    """
    groups = {(op,int(fid)):np.asarray(rows.index,dtype=int)
              for (op,fid),rows in pool.groupby(['operation','faf_id'])}
    unary = pool.score.to_numpy(dtype=float)
    starts = nodes[pool.iaf_id.to_numpy(dtype=int)-1]
    node_distance = np.linalg.norm(starts[:,None,:]-starts[None,:,:],axis=2)
    additive = (traffic.route_cross_weight*crossings.astype(float)+traffic.pair_gap_weight*gaps
                +config.node_weight*np.exp(-(node_distance/config.node_neighborhood_m)**2))
    pairs = [p for p in sorted(options) if options[p]]
    rng = np.random.default_rng(traffic.seed)

    def objective(ids):
        if not ids:
            return math.inf
        ids = np.asarray(ids,dtype=int)
        sub = np.ix_(ids,ids)
        minimum = float(np.min(distances[sub][np.triu_indices(len(ids),1)])) if len(ids)>1 else traffic.spacing_target_m
        return (float(unary[ids].sum()+np.triu(additive[sub],1).sum())
                +traffic.minimum_gap_weight*max(0.,1.-minimum/traffic.spacing_target_m)**2)

    def best_single(candidates,others):
        if len(others):
            sub = np.ix_(candidates,others)
            costs = unary[candidates]+additive[sub].sum(axis=1)
            mins = distances[sub].min(axis=1)
            if len(others)>1:
                oo = distances[np.ix_(others,others)]
                mins = np.minimum(mins,oo[np.triu_indices(len(others),1)].min())
            costs += traffic.minimum_gap_weight*np.maximum(0.,1.-mins/traffic.spacing_target_m)**2
        else:
            costs = unary[candidates]
        return int(candidates[np.argmin(costs)])

    best, best_cost = {},math.inf
    for run in range(traffic.multistarts):
        chosen = {}
        order = list(pairs)
        if run:
            rng.shuffle(order)
        if initial_indices is not None:
            for pair in order:
                star = [i for i in initial_indices if int(pool.loc[i,'pair'])==pair and pool.loc[i,'operation']=='STAR']
                sid = [i for i in initial_indices if int(pool.loc[i,'pair'])==pair and pool.loc[i,'operation']=='SID']
                if len(star)!=1 or len(sid)!=1:
                    raise ValueError('Invalid joint-network initial solution.')
                chosen[pair] = (star[0],sid[0])
        for pair in ([] if initial_indices is not None else order):
            orientation = options[pair][int(rng.integers(len(options[pair])))] if run else options[pair][0]
            star,sid = groups['STAR',orientation[0]],groups['SID',orientation[1]]
            others = [i for values in chosen.values() for i in values]
            si = best_single(star,others)
            di = best_single(sid,others+[si])
            chosen[pair] = (si,di)
        current = objective([i for values in chosen.values() for i in values])
        for sweep in range(traffic.max_sweeps):
            improved = False
            for pair in order:
                others = [i for p,values in chosen.items() if p!=pair for i in values]
                block,block_cost = chosen[pair],current
                for a,b in options[pair]:
                    star,sid = groups['STAR',a],groups['SID',b]
                    cstar,csid = unary[star].copy(),unary[sid].copy()
                    dmin = distances[np.ix_(star,sid)].copy()
                    if others:
                        so,do = np.ix_(star,others),np.ix_(sid,others)
                        cstar += additive[so].sum(axis=1)
                        csid += additive[do].sum(axis=1)
                        dmin = np.minimum(dmin,distances[so].min(axis=1)[:,None])
                        dmin = np.minimum(dmin,distances[do].min(axis=1)[None,:])
                        if len(others)>1:
                            oo = distances[np.ix_(others,others)]
                            dmin = np.minimum(dmin,oo[np.triu_indices(len(others),1)].min())
                    costs = (cstar[:,None]+csid[None,:]+additive[np.ix_(star,sid)]
                             +traffic.minimum_gap_weight*np.maximum(0.,1.-dmin/traffic.spacing_target_m)**2)
                    if forbid_crossings:
                        hard_conflicts = crossings if forbidden_pairs is None else forbidden_pairs
                        invalid = hard_conflicts[np.ix_(star,sid)].copy()
                        if minimum_spacing_m:
                            invalid |= dmin<minimum_spacing_m-1e-7
                        if others:
                            invalid |= hard_conflicts[np.ix_(star,others)].any(axis=1)[:,None]
                            invalid |= hard_conflicts[np.ix_(sid,others)].any(axis=1)[None,:]
                        costs[invalid] = np.inf
                        if not np.isfinite(costs).any():
                            continue
                    i,j = np.unravel_index(np.argmin(costs),costs.shape)
                    si,di = int(star[i]),int(sid[j])
                    cost = objective(others+[si,di])
                    if cost < block_cost-1e-8:
                        block,block_cost = (si,di),cost
                if block != chosen[pair]:
                    chosen[pair],current = block,block_cost
                    improved = True
            if not improved:
                break
        if current < best_cost:
            best,best_cost = chosen.copy(),current
    indices = [i for p in sorted(best) for i in best[p]]
    if not indices:
        return [],dict(score=None,selected_pairs=0,selected_routes=0)
    sub = np.ix_(indices,indices)
    upper = np.triu_indices(len(indices),1)
    return indices,dict(score=best_cost,selected_pairs=len(best),selected_routes=len(indices),
                        route_crossing_pairs=int(crossings[sub][upper].sum()),
                        minimum_distance_m=float(distances[sub][upper].min()),
                        below_target_pairs=int((distances[sub][upper]<traffic.spacing_target_m).sum()),
                        global_optimality_proven=False,multistarts=traffic.multistarts)


def traffic_zero_crossing_search(pool,options,crossings,seed_indices,seconds=20.):
    """Time-bounded MILP feasibility with lazy route-conflict constraints.

    A feasible incumbent is always checked against the COMPLETE conflict
    matrix. An incomplete/time-limited search never certifies infeasibility.
    Infeasibility of a relaxation certifies only this spatial shortlist.
    """
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    active = [p for p in sorted(options) if options[p]]
    n = len(pool)
    base_rows = []
    for pair in active:
        for operation in ('STAR','SID'):
            base_rows.append(pool.index[(pool.pair==pair)&(pool.operation==operation)].to_numpy(dtype=int))
        for faf_id in pool.loc[pool.pair==pair,'faf_id'].unique():
            base_rows.append(pool.index[(pool.pair==pair)&(pool.faf_id==faf_id)].to_numpy(dtype=int))
    br = np.concatenate([np.full(len(ids),i,dtype=int) for i,ids in enumerate(base_rows)])
    bc = np.concatenate(base_rows)
    edges = set()
    def add_conflicts(indices):
        for i in indices:
            for j in np.flatnonzero(crossings[i]):
                if i!=j:
                    edges.add((min(int(i),int(j)),max(int(i),int(j))))
    add_conflicts(seed_indices)
    upper = pool.pair.isin(active).to_numpy(dtype=float)
    objective = pool.score.to_numpy(dtype=float)
    objective = (objective-objective.min())/max(1.,float(np.ptp(objective)))
    started = time.monotonic()
    attempts = 0
    while time.monotonic()-started<seconds:
        edge_array = np.array(sorted(edges),dtype=int).reshape(-1,2)
        er = np.repeat(np.arange(len(edge_array))+len(base_rows),2)
        rows = np.r_[br,er]
        columns = np.r_[bc,edge_array.ravel()]
        matrix = coo_matrix((np.ones(len(rows)),(rows,columns)),shape=(len(base_rows)+len(edge_array),n)).tocsc()
        lower = np.r_[np.ones(len(base_rows)),np.full(len(edge_array),-np.inf)]
        high = np.ones(len(base_rows)+len(edge_array))
        remaining = seconds-(time.monotonic()-started)
        if remaining<=0:
            break
        solved = milp(objective,integrality=np.ones(n),bounds=Bounds(np.zeros(n),upper),
                      constraints=LinearConstraint(matrix,lower,high),
                      options={'time_limit':min(5.,remaining),'mip_rel_gap':.05})
        attempts += 1
        if solved.status==2:
            return None,dict(status='infeasible_in_shortlist',attempts=attempts,
                             seconds=time.monotonic()-started,not_a_full_candidate_infeasibility_proof=True)
        if solved.x is None:
            continue
        chosen = np.flatnonzero(solved.x>.5).tolist()
        if len(chosen)!=2*len(active) or any(sum(solved.x[ids]>.5)!=1 for ids in base_rows):
            continue
        sub = crossings[np.ix_(chosen,chosen)]
        if not np.triu(sub,1).any():
            return chosen,dict(status='zero_crossing_found',attempts=attempts,
                               seconds=time.monotonic()-started)
        old_count = len(edges)
        add_conflicts(chosen)
        if len(edges)==old_count:
            raise RuntimeError('MILP incumbent violates already-added route conflicts.')
    return None,dict(status='time_limit_no_zero_crossing_incumbent',attempts=attempts,
                     seconds=time.monotonic()-started,infeasibility_proven=False)


def traffic_noncrossing_subset(pool,options,crossings,seconds=20.,distances=None,
                              minimum_spacing_m=0.,required_pairs=None):
    """Maximize complete opposite-FAF pairs subject to zero XY intersections.

    Every selected pair has both an Approach and a Departure. Building binary
    variables allow omission to expose the coverage cost of no crossings;
    no altitude/gradient/direction/shape condition is relaxed. Optimality,
    when proved, applies only to this explicitly saved candidate shortlist.
    """
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    active = [p for p in sorted(options) if options[p]]
    n,k = len(pool),len(active)
    if not k:
        return [],dict(status='no_eligible_pairs',coverage_optimal_in_shortlist=True,selected_pairs=0)
    row_indices,column_indices,values = [],[],[]
    row_count = 0
    for pi,pair in enumerate(active):
        groups = [pool.index[(pool.pair==pair)&(pool.operation==op)].to_numpy(dtype=int) for op in ('STAR','SID')]
        groups += [pool.index[(pool.pair==pair)&(pool.faf_id==fid)].to_numpy(dtype=int)
                   for fid in pool.loc[pool.pair==pair,'faf_id'].unique()]
        for ids in groups:
            row_indices.extend([row_count]*(len(ids)+1))
            column_indices.extend([*ids,n+pi])
            values.extend([*np.ones(len(ids)),-1.])
            row_count += 1
    if required_pairs is not None:
        row_indices.extend([row_count]*k)
        column_indices.extend(range(n,n+k))
        values.extend(np.ones(k))
        row_count += 1
    equality_rows = row_count
    conflicts = crossings if not minimum_spacing_m else (crossings|(distances<minimum_spacing_m-1e-7))
    ii,jj = np.where(np.triu(conflicts,1))
    pair = pool.pair.to_numpy()
    op = pool.operation.to_numpy()
    fid = pool.faf_id.to_numpy()
    possible = (np.isin(pair[ii],active)&np.isin(pair[jj],active)
                &~((pair[ii]==pair[jj])&((op[ii]==op[jj])|(fid[ii]==fid[jj]))))
    ii,jj = ii[possible],jj[possible]
    er = np.repeat(np.arange(len(ii))+row_count,2)
    ec = np.column_stack((ii,jj)).ravel()
    matrix = coo_matrix((np.r_[values,np.ones(len(ec))],
                          (np.r_[row_indices,er],np.r_[column_indices,ec])),
                         shape=(row_count+len(ii),n+k)).tocsc()
    lower = np.r_[np.zeros(equality_rows),np.full(len(ii),-np.inf)]
    high = np.r_[np.zeros(equality_rows),np.ones(len(ii))]
    if required_pairs is not None:
        lower[equality_rows-1] = high[equality_rows-1] = required_pairs
    unary = pool.score.to_numpy(dtype=float)
    objective = np.r_[unary/max(float(np.max(np.abs(unary))),1.)*.0001,np.full(k,-1000.)]
    upper = np.r_[pool.pair.isin(active).to_numpy(dtype=float),np.ones(k)]
    started = time.monotonic()
    solved = milp(objective,integrality=np.ones(n+k),bounds=Bounds(np.zeros(n+k),upper),
                  constraints=LinearConstraint(matrix,lower,high),
                  options={'time_limit':seconds,'mip_rel_gap':.0001})
    selected = [] if solved.x is None else np.flatnonzero(solved.x[:n]>.5).tolist()
    if selected:
        if np.triu(crossings[np.ix_(selected,selected)],1).any():
            raise RuntimeError('Noncrossing subset solver returned intersecting routes.')
        for p,rows in pool.loc[selected].groupby('pair'):
            if len(rows)!=2 or set(rows.operation)!={'STAR','SID'} or rows.faf_id.nunique()!=2:
                raise RuntimeError('Noncrossing subset violates complete opposite-FAF pairing.')
    return selected,dict(status=str(solved.message),selected_pairs=len(selected)//2,
                         status_code=int(solved.status),minimum_spacing_m=minimum_spacing_m,
                         coverage_optimal_in_shortlist=bool(solved.status==0),
                         not_a_full_candidate_optimality_proof=True,seconds=time.monotonic()-started,
                         crossing_constraints=len(ii))


def traffic_maxmin_spacing(pool,options,crossings,distances,fafs,seed_indices,
                           seconds=20.,initial_goal_m=100.,tolerance_m=2.):
    """Bounded max-min spacing search while preserving zero crossings/coverage.

    An exact endpoint-distance bound is derived over all eligible building
    subsets of the retained cardinality. MILP probes impose a distance floor.
    A timeout without an incumbent NEVER reduces the upper bound.
    """
    from itertools import combinations
    k = len(seed_indices)//2
    if not k:
        return seed_indices,dict(status='no_selected_pairs')
    active = [p for p in sorted(options) if options[p]]
    by_pair = {p:np.array([f['xyz'][:2] for f in fafs if f['pair']==p]) for p in active}
    endpoint_bound = 0.
    for group in combinations(active,k):
        xy = np.concatenate([by_pair[p] for p in group])
        d = np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=2)
        np.fill_diagonal(d,np.inf)
        endpoint_bound = max(endpoint_bound,float(d.min()))
    def minimum(ids):
        d = distances[np.ix_(ids,ids)]
        return float(d[np.triu_indices(len(ids),1)].min())
    best = list(seed_indices)
    low,high = minimum(best),endpoint_bound
    started = time.monotonic()
    probes = []
    unknown = False
    while high-low>tolerance_m and time.monotonic()-started<seconds-1.:
        if not probes and low<initial_goal_m<high:
            goal = initial_goal_m
        else:
            goal = (low+high)/2
        remaining = seconds-(time.monotonic()-started)
        ids,status = traffic_noncrossing_subset(pool,options,crossings,min(5.,remaining),
            distances=distances,minimum_spacing_m=goal,required_pairs=k)
        probes.append(dict(goal_m=goal,**status))
        if len(ids)==2*k:
            actual = minimum(ids)
            if actual<goal-1e-6:
                raise RuntimeError('Spacing feasibility probe returned an invalid incumbent.')
            best,low = ids,max(low,actual)
        elif status.get('status_code')==2:
            high = goal
        else:
            unknown = True
            break
    return best,dict(status='bounded_maxmin_search',retained_pairs=k,
                     achieved_minimum_m=low,upper_bound_m=high,endpoint_upper_bound_m=endpoint_bound,
                     tolerance_m=tolerance_m,maxmin_optimal_to_tolerance_in_shortlist=(not unknown and high-low<=tolerance_m),
                     seconds=time.monotonic()-started,probes=probes)


def traffic_strict_layout(pool,options,strict_pairs,distances,crossings,gaps,nodes,headings,fafs,
                          areas,audit,config,traffic,directory,scenario,seconds=20.):
    """Preferred hard-constraint frontier: no Landing_Area overlap or XY crossing."""
    allowed = {p:v for p,v in options.items() if p in strict_pairs}
    ids,solver = traffic_noncrossing_subset(pool,allowed,crossings,seconds)
    if ids:
        ids,spacing = traffic_maxmin_spacing(pool,allowed,crossings,distances,fafs,ids,
            seconds=seconds,initial_goal_m=traffic.spacing_target_m)
        selected_pairs = set(pool.loc[ids,'pair'].astype(int))
        selected_options = {p:v for p,v in allowed.items() if p in selected_pairs}
        ids,selection = traffic_select(pool,selected_options,distances,crossings,gaps,nodes,traffic,config,
            initial_indices=ids,forbid_crossings=True,
            minimum_spacing_m=max(0.,spacing.get('achieved_minimum_m',0.)-1e-6))
    else:
        selection = dict(selected_pairs=0,selected_routes=0,route_crossing_pairs=0)
        spacing = dict(status='no_selected_pairs')
    result,pairs = traffic_export(pool,ids,directory,nodes,headings,fafs,areas,audit,config,traffic,scenario)
    if not result.empty and (result.cross_n.astype(int).sum()!=0 or any(r['intersects'] for r in pairs)):
        raise RuntimeError('Strict layout violated Landing_Area or route-intersection hard constraint.')
    return dict(**selection,solver=solver,spacing_search=spacing,
                missing_pairs=sorted(set(range(1,audit['pair_count']+1))-set(result.pair.astype(int))),
                length_m=float(result.length_m.sum()),landing_crossings=int(result.cross_n.sum()),
                exact_crossing_pairs=sum(r['intersects'] for r in pairs),
                exact_min_distance_m=min((r['min_xy_m'] for r in pairs),default=None))


def write_connection_points(routes,directory):
    """Write actual connected endpoints, preserving each route's associations."""
    records,points = [],[]
    for _,row in routes.iterrows():
        xyz = np.asarray(row.geometry.coords)
        is_approach = row.operation=='STAR'
        for role,index in (('IAF',0 if is_approach else -1),('FAF',-1 if is_approach else 0)):
            point = xyz[index]
            speed = row.v_start if index==0 else row.v_end
            records.append(dict(point_id=len(records)+1,route_id=int(row.cand_id),pair=int(row.pair),
                                operation=row.operation,role=role,iaf_id=int(row.iaf_id),
                                faf_id=int(row.faf_id),alt=float(point[2]),speed_kt=float(speed)))
            points.append(Point(*point))
    result = gpd.GeoDataFrame(pd.DataFrame(records,columns=['point_id','route_id','pair','operation',
        'role','iaf_id','faf_id','alt','speed_kt']),geometry=points,crs=routes.crs)
    for name,data in [('connected_points_z',result),('IAF_points_z',result[result.role=='IAF']),
                      ('FAF_points_z',result[result.role=='FAF'])]:
        path = directory/f'{name}.shp'
        data.to_file(path,engine='pyogrio',encoding='UTF-8',geometry_type='Point Z')
        saved = gpd.read_file(path)
        if len(saved)!=len(data) or not saved.has_z.all():
            raise RuntimeError('Connected endpoint PointZ count/Z verification failed.')
        for actual,expected in zip(saved.geometry,data.geometry):
            if not np.allclose(actual.coords[0],expected.coords[0],rtol=0,atol=1e-7):
                raise RuntimeError('Connected endpoint XYZ verification failed.')
    return result


def traffic_export(pool,indices,directory,nodes,headings,fafs,areas,audit,config,traffic,scenario):
    """Export separate STAR/SID + combined PolyLineZ, vertices, exact diagnostics."""
    directory.mkdir(parents=True,exist_ok=True)
    records,geometries,vertex_rows,landing_rows = [],[],[],[]
    center,ring_radius = np.array(audit['ring_center']),audit['ring_radius_m']
    boundary = Polygon(nodes)
    for index in indices:
        row = pool.loc[index]
        xyz,s = traffic_xyz(row,nodes,headings,fafs,config)
        is_star = row.operation=='STAR'
        faf = fafs[int(row.faf_index)]
        ni = int(row.iaf_id)-1
        expected_start = np.r_[nodes[ni],config.ring_z] if is_star else faf['xyz']
        expected_end = faf['xyz'] if is_star else np.r_[nodes[ni],config.ring_z]
        if not np.allclose(xyz[0],expected_start,rtol=0,atol=1e-7) or not np.allclose(xyz[-1],expected_end,rtol=0,atol=1e-7):
            raise RuntimeError('STAR/SID endpoint XYZ validation failed.')
        dz = np.diff(xyz[:,2])
        if (is_star and dz.max()>1e-8) or (not is_star and dz.min()<-1e-8):
            raise RuntimeError('Incorrect STAR descent / SID climb direction.')
        slope = np.abs(dz)/np.linalg.norm(np.diff(xyz[:,:2],axis=0),axis=1)
        limit = min(math.tan(math.radians(config.max_slope_deg)),config.max_ft_nm*.3048/1852)
        if slope.max()>limit+1e-8 or row.bank_deg>config.max_bank+1e-7:
            raise RuntimeError('STAR/SID slope or bank limit violation.')
        start,heading,end,_ = traffic_parameters(row,nodes,headings,fafs,config)
        end_heading = heading-(row.arc1_m+row.arc2_m)/row.radius_m
        wanted_heading = faf['heading'] if is_star else headings[ni]
        if abs((end_heading-wanted_heading+math.pi)%TAU-math.pi)>1e-8:
            raise RuntimeError('FAF alignment or clockwise ring-merge heading mismatch.')
        turn_deg = math.degrees((row.arc1_m+row.arc2_m)/row.radius_m)
        if (turn_deg>traffic.max_total_turn_deg+1e-7
                or math.degrees(max(row.arc1_m,row.arc2_m)/row.radius_m)>traffic.max_arc_deg+1e-7
                or not LineString(xyz[:,:2]).is_simple):
            raise RuntimeError('Loop-free route shape validation failed.')
        analytic_radial = maximum_ring_radius(start,heading,row.radius_m,
            [row.arc1_m,row.straight_m,row.arc2_m],center)
        if scenario=='inside_only' and analytic_radial>ring_radius+1e-7:
            raise RuntimeError('Analytic arc exits ring in inside-only scenario.')
        line = LineString(xyz)
        count,overlap,gaps,items = landing_metrics(xyz,areas,details=True)
        if count!=int(row.cross_n):
            raise RuntimeError('Landing_Area crossing roundtrip mismatch.')
        rec = {key:row[key] for key in ['cand_id','pair','faf_id','iaf_id','radius_m','length_m',
                                        'slope_deg','bank_deg','cross_n','overlap_m']}
        rec.update(operation=row.operation,z_start=xyz[0,2],z_end=xyz[-1,2],
                   max_rad_m=analytic_radial,
                   outside_m=(0. if analytic_radial<=ring_radius+1e-7 else
                              LineString(xyz[:,:2]).difference(boundary).length),
                   path_type='RSR',turn_deg=turn_deg,v_start=config.start_kt if is_star else config.end_kt,
                   v_end=config.end_kt if is_star else config.start_kt)
        records.append(rec)
        geometries.append(line)
        v0,v1 = rec['v_start']*KNOT,rec['v_end']*KNOT
        s3d = s*math.hypot(1.,(xyz[-1,2]-xyz[0,2])/s[-1])
        speed = np.sqrt(v0*v0+(v1*v1-v0*v0)*s/s[-1])
        for j,point in enumerate(xyz):
            vertex_rows.append(dict(operation=row.operation,pair=int(row.pair),faf_id=int(row.faf_id),
                                    vertex=j,x=point[0],y=point[1],z=point[2],s_m=s[j],s3d_m=s3d[j],
                                    speed_kt=speed[j]/KNOT,time_s=2*s3d[j]/(v0+speed[j])))
        landing_rows.extend(dict(operation=row.operation,pair=int(row.pair),**item) for item in items)
    columns = ['cand_id','pair','faf_id','iaf_id','radius_m','length_m','slope_deg','bank_deg',
               'cross_n','overlap_m','operation','z_start','z_end','max_rad_m','outside_m','path_type','turn_deg','v_start','v_end']
    result = gpd.GeoDataFrame(pd.DataFrame(records,columns=columns),geometry=geometries,crs=5186)
    for filename,data in [('combined_routes_z.shp',result),
                          ('Approach_z.shp',result[result.operation=='STAR']),
                          ('Departure_z.shp',result[result.operation=='SID'])]:
        path = directory/filename
        data.to_file(path,engine='pyogrio',encoding='UTF-8',geometry_type='LineString Z')
        saved = gpd.read_file(path)
        if len(saved)!=len(data) or saved.crs.to_epsg()!=5186 or not saved.has_z.all():
            raise RuntimeError('SHP Z/CRS/count validation failed.')
        for actual,expected in zip(saved.geometry,data.geometry):
            if not np.allclose(np.asarray(actual.coords),np.asarray(expected.coords),rtol=0,atol=1e-7):
                raise RuntimeError('SHP XYZ validation failed.')
    result.drop(columns='geometry').to_csv(directory/'selected_routes.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(vertex_rows,columns=['operation','pair','faf_id','vertex','x','y','z','s_m','s3d_m','speed_kt','time_s']).to_csv(
        directory/'vertices.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(landing_rows,columns=['operation','pair','area_pair','area_fid','overlap_m','min_abs_z','min_signed_z']).to_csv(
        directory/'landing_crossings.csv',index=False,encoding='utf-8-sig')
    # Exact pairwise audit on full exported polylines, including all endpoints.
    pairs,crossing_points = [],[]
    for i in range(len(result)):
        for j in range(i+1,len(result)):
            a,b = result.iloc[i],result.iloc[j]
            la,lb = LineString(np.asarray(a.geometry.coords)[:,:2]),LineString(np.asarray(b.geometry.coords)[:,:2])
            distance = la.distance(lb)
            hit = la.intersection(lb)
            pairs.append(dict(route_a=i+1,operation_a=a.operation,pair_a=int(a.pair),
                              route_b=j+1,operation_b=b.operation,pair_b=int(b.pair),
                              min_xy_m=distance,intersects=not hit.is_empty,
                              intersection_type=hit.geom_type,below_target=distance<traffic.spacing_target_m))
            if not hit.is_empty:
                if hit.geom_type=='Point':
                    points = [hit]
                elif hit.geom_type=='MultiPoint':
                    points = list(hit.geoms)
                else:
                    points = [Point(p) for p in get_coordinates(hit)]
                crossing_points.extend(points)
    pd.DataFrame(pairs,columns=['route_a','operation_a','pair_a','route_b','operation_b','pair_b',
                               'min_xy_m','intersects','intersection_type','below_target']).to_csv(
        directory/'route_pair_spacing.csv',index=False,encoding='utf-8-sig')
    gpd.GeoDataFrame(pd.DataFrame([dict(pair=int(r.pair),operation=r.operation,faf_id=int(r.faf_id),alt=float(fafs[int(pool.loc[idx,'faf_index'])]['xyz'][2]))
        for idx,r in zip(indices,result.itertuples())],columns=['pair','operation','faf_id','alt']),
        geometry=[Point(*fafs[int(pool.loc[idx,'faf_index'])]['xyz']) for idx in indices],crs=5186).to_file(
            directory/'assigned_faf_z.shp',engine='pyogrio',encoding='UTF-8',geometry_type='Point Z')
    write_connection_points(result,directory)
    traffic_plot(result,areas,nodes,crossing_points,directory,scenario,config)
    if config.separation:
        raise ValueError('Use legacy --star-only for old corridor separation; traffic model currently optimizes XY spacing.')
    for pair,rows in result.groupby('pair'):
        if set(rows.operation)!= {'STAR','SID'} or rows.faf_id.nunique()!=2:
            raise RuntimeError('Opposite FAF STAR/SID pairing validation failed.')
    return result,pairs


def traffic_plot(result,areas,nodes,crossing_points,output,scenario,config):
    colors = plt.get_cmap('tab10')
    ring = np.vstack([nodes,nodes[0]])
    for operation in ('combined','STAR','SID'):
        fig,ax = plt.subplots(figsize=(13,12))
        for area in areas:
            xy = np.asarray(area['geometry'].exterior.coords)
            ax.fill(*xy.T,color='tomato',alpha=.10)
            ax.plot(*xy.T,color='tomato',lw=.5,alpha=.6)
        ax.plot(*ring.T,'--',color='black',lw=1,label=f'Inner ring CW, {config.ring_z:g} m')
        for _,row in result.iterrows():
            if operation!='combined' and row.operation!=operation:
                continue
            xyz = np.asarray(row.geometry.coords)
            color = colors((int(row.pair)-1)%10)
            ax.plot(*xyz[:,:2].T,color=color,ls='-' if row.operation=='STAR' else '--',lw=1.5,
                    label=f'{row.operation} P{int(row.pair)} F{int(row.faf_id)} ({row.length_m:.0f} m)')
            # Place arrows on arc/straight interior and distinguish ring/FAF.
            k = min(len(xyz)-2,len(xyz)//2)
            ax.annotate('',xy=xyz[k+1,:2],xytext=xyz[k,:2],arrowprops=dict(arrowstyle='->',color=color))
            faf = xyz[-1] if row.operation=='STAR' else xyz[0]
            ax.scatter(*faf[:2],marker='*' if row.operation=='STAR' else '^',s=65,color=color)
        if operation=='combined' and crossing_points:
            ax.scatter([p.x for p in crossing_points],[p.y for p in crossing_points],marker='x',s=60,
                       color='red',zorder=10,label='Route intersections')
        ax.set(aspect='equal',xlabel='EPSG:5186 X (m)',ylabel='EPSG:5186 Y (m)',
               title=f'{scenario}: {operation} / solid STAR, dashed SID')
        ax.ticklabel_format(style='plain',useOffset=False)
        ax.legend(fontsize=7,loc='upper left',ncol=2)
        ax.grid(alpha=.2)
        fig.tight_layout()
        fig.savefig(output/f'map_{operation.lower()}.png',dpi=160)
        plt.close(fig)
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(13,10))
    ax = fig.add_subplot(111,projection='3d')
    all_xyz = []
    for area in areas:
        ax.add_collection3d(Poly3DCollection([area['xyz']],facecolor='tomato',alpha=.10))
        all_xyz.extend(area['xyz'])
    for _,row in result.iterrows():
        xyz = np.asarray(row.geometry.coords)
        all_xyz.extend(xyz)
        ax.plot(*xyz.T,color=colors((int(row.pair)-1)%10),ls='-' if row.operation=='STAR' else '--')
    xyz = np.asarray(all_xyz)
    ax.auto_scale_xyz(xyz[:,0],xyz[:,1],xyz[:,2])
    ax.set(xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)',title=f'{scenario}: STAR descent / SID climb')
    fig.savefig(output/'view3d.png',dpi=160)
    plt.close(fig)
    fig,axes = plt.subplots(4,2,figsize=(14,13),squeeze=False)
    for pair,rows in result.groupby('pair'):
        ax = axes.flat[int(pair)-1]
        for _,row in rows.iterrows():
            xyz = np.asarray(row.geometry.coords)
            distance = np.r_[0,np.cumsum(np.linalg.norm(np.diff(xyz[:,:2],axis=0),axis=1))]
            ax.plot(distance,xyz[:,2],ls='-' if row.operation=='STAR' else '--',label=f'{row.operation} F{int(row.faf_id)}')
        ax.set(title=f'Pair {int(pair)}',xlabel='Horizontal distance (m)',ylabel='Z (m)')
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(output/'profiles.png',dpi=150)
    plt.close(fig)


def traffic_report(output,summary,tables):
    labels = dict(inside_only='링 내부 제한',outside_allowed='링 밖 우회 허용',outside_detour='링 밖 우회 필수')
    review_notice = ('<p class="note">형태 수정 검토 중: FAF 부근의 큰 원형 선회가 포함된 기존 기하 후보입니다. '
                     '현재 SHP 분리는 완료했으며 항로 모양은 아직 변경하지 않았습니다. '
                     '이 결과를 최종 Approach/Departure 설계안으로 해석하지 마세요.</p>'
                     if summary.get('geometry_revision_pending',False) else '')
    body = []
    for scenario,info in summary['scenarios'].items():
        strict = info.get('strict_noncrossing',{})
        sg = strict.get('exact_min_distance_m')
        sg_label = f'{sg:.2f}m' if sg is not None else '선정 항로 없음'
        strict_html = (f'''<h3>불관용 기준 비교안: Landing_Area 관통 0건 + 항로 교차 0건</h3>
<p>FAF {strict['selected_pairs']}쌍 / Approach·Departure 합계 {strict['selected_routes']}개.
최소 항로 간격 {sg_label}. 미연결 쌍: {strict['missing_pairs']}.
선별 후보 안에서 연결 쌍 수 최적성 입증: {strict['solver']['coverage_optimal_in_shortlist']}.</p>
<p><a href="{scenario}/strict_noncrossing/Approach_z.shp">불관용 Approach SHP</a> ·
<a href="{scenario}/strict_noncrossing/Departure_z.shp">불관용 Departure SHP</a></p>
<img src="{scenario}/strict_noncrossing/map_combined.png">'''
                       if strict else '')
        gap_label = f"{info['exact_min_distance_m']:.2f}m" if info['exact_min_distance_m'] is not None else '선정 항로 없음'
        subset = info.get('noncrossing_subset',{})
        subset_gap = subset.get('exact_min_distance_m')
        subset_gap_label = f'{subset_gap:.2f}m' if subset_gap is not None else '선정 항로 없음'
        subset_html = (f'''<h3>Landing_Area 페널티 허용: 교차 0건 부분 구성</h3><p>연결한 FAF 쌍 {subset['selected_pairs']}개,
Approach/Departure 합계 {subset['selected_routes']}개. 교차 {subset['exact_crossing_pairs']}쌍,
최소 간격 {subset_gap_label}. Landing_Area 관통 {subset['landing_crossings']}건. 미연결 쌍: {subset['missing_pairs']}.
선별 후보 안에서 연결 쌍 수 최적성 입증: {subset['solver']['coverage_optimal_in_shortlist']}.</p>
<p><a href="{scenario}/noncrossing_subset/Approach_z.shp">교차 0건 Approach SHP</a> ·
<a href="{scenario}/noncrossing_subset/Departure_z.shp">교차 0건 Departure SHP</a></p>
<details><summary>교차 0건 지도</summary><img src="{scenario}/noncrossing_subset/map_combined.png"></details>'''
                       if subset else '')
        table = tables[scenario][['operation','pair','faf_id','iaf_id','length_m','slope_deg','bank_deg','outside_m']].to_html(
            index=False,float_format=lambda v:f'{v:.2f}')
        body.append(f'''<h2>{labels[scenario]}</h2>{strict_html}<h3>연결 범위 우선 비교안 — 교차가 남으면 최종 무교차 항로가 아님</h3>
<p>STAR {info['star_count']}개 / SID {info['sid_count']}개.
항로 교차 {info['exact_crossing_pairs']}쌍, 최소 XY 간격 {gap_label},
길이 합 {info['length_m']/1000:.2f}km, 링 밖 우회 {info['outside_routes']}개.
Landing_Area 관통 {info['landing_crossings']}건 ({info['landing_mode']}).</p>
<p>선정하지 못한 쌍: {info['missing_pairs']}. 후보 수 {info['eligible_candidates']:,}, 조합 탐색용 후보 {info['shortlisted_candidates']:,}.</p>
<p>교차 0건 추가 탐색: {info.get('zero_crossing_search',{}).get('status','not_run')}.
시간 제한으로 해를 찾지 못한 경우 교차 없는 해가 존재하지 않는다는 뜻은 아닙니다.</p>
<p><a href="{scenario}/Approach_z.shp">Approach 전용 SHP</a> ·
<a href="{scenario}/Departure_z.shp">Departure 전용 SHP</a> ·
<a href="{scenario}/route_pair_spacing.csv">항로 간격·교차 CSV</a></p>
<details><summary>통합 비교 지도 (Approach와 Departure SHP는 별도)</summary><img src="{scenario}/map_combined.png"></details>
<h3>Approach</h3><img src="{scenario}/map_star.png"><h3>Departure</h3><img src="{scenario}/map_sid.png">
<details><summary>3D / 고도 / 선정 표</summary><img src="{scenario}/view3d.png"><img src="{scenario}/profiles.png">{table}</details>
{subset_html}''')
    html = f'''<!doctype html><html lang="ko"><meta charset="utf-8"><title>STAR/SID 항로 비교</title>
<style>body{{font-family:Malgun Gothic,sans-serif;max-width:1300px;margin:32px auto;line-height:1.7}}img{{max-width:100%}}table{{border-collapse:collapse}}td,th{{padding:5px;border:1px solid #ccc}}.note{{background:#fff2d8;padding:16px}}</style>
<h1>Inner Ring STAR / 반대편 FAF SID 공동 선정</h1>
{review_notice}
<p>형태 재검토: 원형 우회로 구배를 맞추던 기존 해를 제거했습니다.
총 선회량 ≤{summary['traffic_config'].get('max_total_turn_deg',180):g}°, 개별 원호 ≤{summary['traffic_config'].get('max_arc_deg',180):g}°,
자기교차 금지를 추가했습니다. 이 각도 제한은 이번 모델의 형태 필터이며 FAA/NASA의 공통 승인 기준이 아닙니다.
초기 조건 대조와 수정 근거는 <a href="design_review.md">design_review.md</a>에 기록했습니다.</p>
<p>쌍마다 STAR 하나와 반대편 FAF의 SID 하나를 함께 선정합니다. STAR는 링→FAF, SID는 FAF→링으로
각각 새 시계방향 RSR을 계산합니다. SID를 STAR의 좌표 역순으로 만들지 않습니다.</p>
<p>링 {summary['config']['ring_z']:g}m, STAR {summary['config']['start_kt']:g}→{summary['config']['end_kt']:g}kt,
SID {summary['config']['end_kt']:g}→{summary['config']['start_kt']:g}kt. 상승·하강 모두
{summary['config']['max_slope_deg']:g}°와 {summary['config']['max_ft_nm']:g}ft/NM 중 엄격한 제한 적용.
뱅크 ≤{summary['config']['max_bank']:g}°. 속도는 무풍 TAS로 가정하고 FAF alt를 사용합니다.</p>
<p class="note">항로 간 거리는 전체 중심선의 XY 최소 거리이며 끝점도 포함합니다.
고정된 16개 FAF의 최소 간격이 {summary['fixed_faf_spacing_upper_bound_m']:.2f}m라서 모든 항로 간격을
{summary['traffic_config']['spacing_target_m']:g}m 이상으로 만드는 것은 불가능할 수 있습니다.
간격 목표는 연성 조건이고, 항로 교차를 강하게 벌점 처리합니다. 수직 분리로 XY 교차를 면제하지 않습니다.</p>
<p>inside_only는 원호와 직선의 해석적 최대 반경을 검사하여 링 밖으로 나가지 못하게 합니다.
outside_allowed는 동일 후보에 링 밖 우회를 추가로 허용합니다.
outside_detour는 각 항로가 링 반경보다 1m 이상 벗어나는 후보만 사용하는 추가 비교안입니다.
이전 STAR 단독 결과는 그대로 보존합니다.</p>
{''.join(body)}
<h2>선정 방법과 해석</h2><p>먼저 각 FAF 쌍에서 Landing_Area 관통 없는 STAR/SID 구성이 가능한지 검사합니다.
가능한 쌍은 관통 금지를 유지하고, 불가능한 쌍에만 관통 개수·길이·표면 높이 차 페널티를 적용합니다.
반경·구배·시계방향·링 내부 경계 같은 기하 조건은 완화하지 않습니다.
해당 조건으로 구성할 수 없는 쌍은 보고서와 JSON에 누락으로 명시합니다.</p>
<p>개별 후보 전수 평가 후, 각 운항방향·FAF의 점수 상위 후보와 링 구간별 분산 후보를 모아
STAR/SID 쌍 단위 다중 시작 좌표하강법으로 공동 선정합니다. 교차 쌍 수, 최소 간격 부족,
전체 항로 쌍의 간격 부족, 노드 집중, 항로 길이·뱅크각이 목적함수에 들어갑니다.
전체 후보 조합의 전역 최적성이나 ‘해가 존재하지 않음’을 증명하지 않습니다.
교차가 남으면 제한 시간 내 MILP로 교차 없는 조합을 추가로 찾습니다. 이 탐색이 성공하면 교차 0건을 유지하며 간격을 재개선합니다.
시나리오별 독립 휴리스틱 결과이므로 우회 허용의 수치가 내부 제한보다 항상 우수하다는 보장은 없습니다.</p>
<p>원호는 최대 {summary['config']['chord_error_m']:g}m 현 오차로 저장하고, 최종 XY 교차·간격은 원본 해상도의
저장 폴리라인으로 다시 검사합니다. 간격 분석은 실제 비행 시각을 고려한 충돌 분석이 아닙니다.
기체별 성능 승인·바람·지형·건물 장애물·롤 전이는 포함하지 않습니다.</p>
<h2>재실행</h2><pre>python -B 수정본_20261103/코드/optimization_path.py
# 간격 목표와 교차 페널티 변경 (기하 후보 캐시는 재사용)
python -B 수정본_20261103/코드/optimization_path.py --spacing-target 150 --route-cross-weight 2000000
# 이전 STAR 단독 분석 모드
python -B 수정본_20261103/코드/optimization_path.py --star-only --reuse-metrics</pre>
<p>모든 후보 지표: all_star_sid_candidates.csv. 설정·해시·선정 요약: summary.json.
각 폴더의 vertices.csv에는 Z·속도·시간, route_pair_spacing.csv에는 STAR–STAR, SID–SID,
STAR–SID 모든 조합의 최소 간격과 교차 여부가 있습니다.</p></html>'''
    (output/'report.html').write_text(html,encoding='utf-8')


def traffic_self_check():
    self_check()
    config = Config()
    rng = np.random.default_rng(62)
    for _ in range(100):
        start,end = rng.normal(size=(2,2))*1500
        h0,h1 = rng.uniform(-math.pi,math.pi,2)
        radius = float(rng.uniform(230,750))
        lengths = rsr(start,h0,end,h1,radius)
        xyz,_ = sample_rsr(start,h0,np.r_[end,350.],radius,lengths,replace(config,ring_z=230.))
        analytic = maximum_ring_radius(start,h0,radius,lengths,np.array([0.,0.]))
        sampled = np.linalg.norm(xyz[:,:2],axis=1).max()
        assert sampled <= analytic+1e-7 and analytic-sampled <= config.chord_error_m+.001
        assert np.all(np.diff(xyz[:,2])>=0)
    print('Traffic self-check passed: ascending SID and analytic ring-boundary extrema.',flush=True)


def write_design_review(output,summary):
    """Trace initial requirements, corrections, assumptions and remaining limits."""
    text = '''# 초기 조건 대조 및 설계 재검토

기존 결과는 수학적으로 허용된 RSR 후보를 골랐지만, 큰 원형 선회를 제한하지 않아
구배를 맞추기 위해 거의 한 바퀴 도는 항로가 선택되었다. 또 일부 FAF 쌍에서
Landing_Area 회피가 불가능하다는 이유로 다른 쌍까지 관통 금지를 완화했다.
두 부분을 수정했다. 이 문서는 실제 비행절차 승인 문서가 아니다.

| 요구사항 | 이번 처리 |
|---|---|
| 이번 범위는 내곽 링과 FAF | 외곽 링·Interchange·Outer→IF는 계산하지 않음 |
| IF 제외 | 새로운 IF나 중간 항법 지점을 추가하지 않음 |
| 내곽·외곽 높이 350m | 이번 내곽 시작/합류 높이 350m 고정. 불리해도 임의로 낮추지 않음 |
| FAF 고도는 alt | 기하 Z나 평균 259.3m로 대체하지 않음. 입력 수직 기준과 동일하다고 가정 |
| 링 후보 1000개 | 공급된 내곽 꼭짓점으로 확인한 반경 2000m 링에 1000개 균등 배치 |
| 내곽 속도 70kt, FAF 60kt | Approach 70→60kt, Departure 60→70kt. 원본의 80kt 속성은 사용하지 않음 |
| CSC, 내곽 시계방향 RSR | 양방향 각각 새 RSR 생성. Approach 역순을 Departure로 사용하지 않음 |
| FAF 진입·이탈 방향 | Approach는 FAF→해당 OLS 꼭짓점 방향으로 끝나고, Departure는 그 반대 방향으로 출발 |
| 링 접선 | Approach 시작과 Departure 끝이 시계방향 링 접선에 일치하는지 해석적으로 검증 |
| 최대 하강 2.36°·250ft/NM | 더 엄격한 250ft/NM 적용. 동일 원리 요청에 따라 Departure 상승에도 같은 수치 사용 |
| 뱅크 최대 30°, 최소 선회반경 | 70kt 최대 속도에서 약 229m 이상. 229/350/500/750m 후보 비교 |
| 낮은 뱅크 선호 | 길이와 함께 연성 비용에 반영 |
| FAA/NASA 기준 | 사용한 반경은 협조선회 관계식의 모델값. 기종별 승인 반경이라고 주장하지 않음 |
| Landing_Area 우선 회피 | 각 FAF 쌍에 대해 회피 가능한 경우 hard 유지. 불가능한 쌍만 관통 개수·길이·수직 높이 차 penalty |
| 한 쌍의 FAF 중 하나를 STAR에 사용 | 다른 FAF를 반드시 SID에 사용. 한 FAF를 양쪽에 중복 배정하지 않음 |
| 항로 간 간격 최대한 확보·교차 최소화 | 모든 Approach/Departure 조합의 XY 교차와 최소 간격을 함께 평가. 100m는 변경 가능한 목표이며 hard 보장값 아님 |
| 이전 수직 70m 조건 제외 | Outer→IF와 Inner Ring에 관한 조건이므로 이번 항로 쌍에 재적용하지 않음 |
| 노드 집중 억제 | 링 노드와 주변 노드의 근접성에 penalty |
| 링 내부/외부 비교 | 내부는 원호의 해석적 최대 반경을 검사. 외부 허용은 경계 제약만 해제. 외부 필수는 별도 보조 비교 |
| Approach/Departure SHP 분리 | 각 시나리오에 Approach_z.shp와 Departure_z.shp 별도 저장, 각각 Z 검증 |
| 조건을 만족하지 못한 경우 | 누락 FAF 쌍을 표시. 고도·구배·회전 방향을 임의로 완화하거나 원형 우회로 채우지 않음 |

## 이번 피드백에 따른 추가 형태 필터

총 시계방향 선회량 180° 이하, 각 원호 180° 이하, 자기교차 금지를 기본으로 적용한다.
이는 한 바퀴에 가까운 우회 해를 제거하기 위한 **모델 선택**이며 초기 사용자 수치나
FAA/NASA의 공통 규정이라고 주장하지 않는다. --max-turn-deg와 --max-arc-deg로 조정할 수 있다.
350m와 FAF alt 차이를 이 필터 안에서 해소할 길이가 없으면 그 조합은 제외한다.

## 범위와 남은 한계

FAF 이후 Landing Area→버티포트의 -1m/s² 감속과 35kt 역천이는 처음에 정한
IAF→FAF 및 반대편 FAF→링 범위 밖이므로 이번 경로에 끼워 넣지 않았다.
별도 IF나 마지막 직선 구간을 덧붙이면 요청한 CSC 형태 자체가 바뀌므로 임의 추가하지 않았다.
현재 FAF에서 방향 접선은 맞추지만, 실제 롤 전이·2.36°에서 8°로의 종단 전이·기체 성능은 모델링하지 않는다.
무풍 TAS 가정이며 KIAS/TAS 변환, 실제 운항 시각, 지형·건물 장애물 분석도 포함하지 않는다.
따라서 이 결과는 조건별 기하 후보 분석이며 실제 운항 가능한 절차의 승인이나 안전성 보증이 아니다.

불관용 기준 비교안(strict_noncrossing)은 Landing_Area 관통과 항로 교차를 모두 금지한 채
연결 가능한 FAF 쌍 수를 먼저 최대화하고, 그 수를 유지하며 최소 간격을 개선한다.
Landing_Area 페널티 허용 부분 구성(noncrossing_subset)과 전체 연결 비교안은 별도 파일로 남긴다.
서로 다른 조건의 결과를 하나의 최종 항로처럼 합치지 않는다.

개별 IAF/FAF/반경 후보를 평가한 뒤 공간적으로 분산된 후보 집합에서 조합을 최적화한다.
전역 최적성은 증명하지 않는다. 교차 0건 MILP 탐색이 시간 제한에 도달하면 이를 불가능성 증명으로 쓰지 않는다.
교차가 남으면 최종 CSV와 지도에 그대로 표시한다.
'''
    text += '\n## 실제 실행 요약\n\n```json\n'+json.dumps(dict(
        geometry_review=summary['geometry_review'],scenarios=summary['scenarios']),ensure_ascii=False,indent=2)+'\n```\n'
    (output/'design_review.md').write_text(text,encoding='utf-8')


def main():
    import sys
    if '--inside-new-options' in sys.argv:
        sys.argv.remove('--inside-new-options')
        from optimization_new_options import main as run_new_options
        run_new_options()
        return
    if '--star-only' in sys.argv:
        sys.argv.remove('--star-only')
        legacy_main()
        return
    parser = argparse.ArgumentParser(description='Joint clockwise STAR/SID; inside and outside ring cases.')
    parser.add_argument('--qgis',type=Path,default=BASE/'QGIS')
    parser.add_argument('--output',type=Path,default=BASE/'결과물/optimization_path_reviewed')
    parser.add_argument('--iaf-count',type=int,default=1000)
    parser.add_argument('--ring-z',type=float,default=350.)
    parser.add_argument('--radii',type=float,nargs='+')
    parser.add_argument('--max-bank',type=float,default=30.)
    parser.add_argument('--max-slope-deg',type=float,default=2.36)
    parser.add_argument('--spacing-target',type=float,default=100.)
    parser.add_argument('--route-cross-weight',type=float,default=1_000_000.)
    parser.add_argument('--gap-weight',type=float,default=200.)
    parser.add_argument('--minimum-gap-weight',type=float,default=5000.)
    parser.add_argument('--max-turn-deg',type=float,default=180.,help='Total RSR turn cap; explicit shape filter, not a regulatory value')
    parser.add_argument('--max-arc-deg',type=float,default=180.,help='Individual arc turn cap')
    parser.add_argument('--candidate-limit',type=int,default=36)
    parser.add_argument('--multistarts',type=int,default=80)
    parser.add_argument('--crossing-search-seconds',type=float,default=20.,
                        help='Extra MILP search for a fully nonintersecting shortlist network; 0 disables')
    parser.add_argument('--subset-search-seconds',type=float,default=20.,
                        help='Time limit per zero-crossing subset coverage/spacing search stage')
    parser.add_argument('--self-check',action='store_true')
    parser.add_argument('--generate-only',action='store_true')
    args = parser.parse_args()
    if args.self_check:
        traffic_self_check()
        return
    config = Config(iaf_count=args.iaf_count,ring_z=args.ring_z,radii=tuple(args.radii or ()),
                    max_bank=args.max_bank,max_slope_deg=args.max_slope_deg)
    traffic = TrafficConfig(spacing_target_m=args.spacing_target,route_cross_weight=args.route_cross_weight,
                            pair_gap_weight=args.gap_weight,minimum_gap_weight=args.minimum_gap_weight,
                            candidate_limit=args.candidate_limit,multistarts=args.multistarts,
                            max_total_turn_deg=args.max_turn_deg,max_arc_deg=args.max_arc_deg)
    values = [args.ring_z,args.max_bank,args.max_slope_deg,args.spacing_target,args.route_cross_weight,
              args.gap_weight,args.minimum_gap_weight,args.max_turn_deg,args.max_arc_deg,*config.radii]
    if (not all(math.isfinite(v) and v>0 for v in values) or args.max_bank>=90 or args.max_slope_deg>=90
            or args.iaf_count<2 or args.candidate_limit<1 or args.multistarts<1
            or not math.isfinite(args.crossing_search_seconds) or args.crossing_search_seconds<0
            or not math.isfinite(args.subset_search_seconds) or args.subset_search_seconds<=0):
        raise ValueError('Invalid configuration values.')
    output = args.output.resolve()
    output.mkdir(parents=True,exist_ok=True)
    nodes,headings,fafs,areas,audit,sources = prepare(args.qgis,config)
    hashes = fingerprint(sources)
    frame,generation = traffic_generate(nodes,headings,fafs,areas,config,audit,hashes,output)
    if args.generate_only:
        return
    frame,shape_review = review_candidate_shapes(frame,nodes,headings,fafs,config,traffic,output)
    fafxy = np.array([f['xyz'][:2] for f in fafs])
    fixed = np.linalg.norm(fafxy[:,None,:]-fafxy[None,:,:],axis=2)
    np.fill_diagonal(fixed,np.inf)
    summary = dict(config=asdict(config),traffic_config=asdict(traffic),input_audit=audit,
                   generation=generation,input_sha256=hashes,fixed_faf_spacing_upper_bound_m=float(fixed.min()),
                   joint_opposite_faf_assignment=True,vertical_separation_applied=False,
                   crossing_search_seconds=args.crossing_search_seconds,geometry_review=shape_review,
                   subset_search_seconds=args.subset_search_seconds,
                   geometry_revision_pending=False,scenarios={})
    tables = {}
    for scenario in ('inside_only','outside_allowed','outside_detour'):
        print(f'Optimizing {scenario}...',flush=True)
        pool,options,info = traffic_pool(frame,scenario,config,traffic,audit)
        if pool.empty or not any(options.values()):
            indices = []
            distances = np.zeros((len(pool),len(pool)))
            crossings = np.zeros((len(pool),len(pool)),dtype=bool)
            gaps = np.zeros_like(distances)
            selection = dict(score=None,selected_pairs=0,selected_routes=0,route_crossing_pairs=0,
                             minimum_distance_m=None,below_target_pairs=0,global_optimality_proven=False,
                             status='no_opposite_faf_pair_under_hard_constraints')
        else:
            matrix_key = hashlib.sha256(json.dumps(dict(version=1,config=asdict(config),inputs=hashes,
                candidate_ids=pool.cand_id.astype(int).tolist()),sort_keys=True).encode('utf-8')).hexdigest()
            _,distances,crossings,gaps = traffic_pair_matrices(pool,nodes,headings,fafs,config,traffic,
                output/f'spacing_matrix_{scenario}.npz',matrix_key)
            indices,selection = traffic_select(pool,options,distances,crossings,gaps,nodes,traffic,config)
        zero_info = dict(status='heuristic_already_has_zero_crossings' if not selection['route_crossing_pairs'] else 'disabled')
        if selection['route_crossing_pairs'] and args.crossing_search_seconds:
            print(f'{scenario}: searching for a fully nonintersecting network...',flush=True)
            zero,zero_info = traffic_zero_crossing_search(pool,options,crossings,indices,args.crossing_search_seconds)
            if zero is not None:
                indices,selection = traffic_select(pool,options,distances,crossings,gaps,nodes,traffic,config,
                    initial_indices=zero,forbid_crossings=True)
        info['zero_crossing_search'] = zero_info
        if any(options.values()) and not pool.empty:
            print(f'{scenario}: maximizing coverage with zero route intersections...',flush=True)
            subset_ids,subset_solver = traffic_noncrossing_subset(pool,options,crossings,args.subset_search_seconds)
        else:
            subset_ids,subset_solver = [],dict(status='no_eligible_pairs',coverage_optimal_in_shortlist=True)
        if subset_ids:
            subset_ids,subset_spacing_search = traffic_maxmin_spacing(pool,options,crossings,distances,fafs,
                subset_ids,seconds=args.subset_search_seconds,initial_goal_m=traffic.spacing_target_m)
            spacing_floor = max(0.,subset_spacing_search.get('achieved_minimum_m',0.)-1e-6)
            active_subset = set(pool.loc[subset_ids,'pair'].astype(int))
            subset_options = {p:values for p,values in options.items() if p in active_subset}
            subset_ids,subset_selection = traffic_select(pool,subset_options,distances,crossings,gaps,nodes,traffic,config,
                initial_indices=subset_ids,forbid_crossings=True,minimum_spacing_m=spacing_floor)
            if subset_selection['selected_pairs']==selection['selected_pairs'] and selection['route_crossing_pairs']:
                indices,selection = subset_ids,subset_selection
        else:
            subset_selection = dict(selected_pairs=0,selected_routes=0,route_crossing_pairs=0)
            subset_spacing_search = dict(status='no_selected_pairs')
        folder = output/scenario
        result,pair_rows = traffic_export(pool,indices,folder,nodes,headings,fafs,areas,audit,config,traffic,scenario)
        subset_result,subset_spacing = traffic_export(pool,subset_ids,folder/'noncrossing_subset',
            nodes,headings,fafs,areas,audit,config,traffic,scenario)
        if subset_spacing and min(r['min_xy_m'] for r in subset_spacing)<spacing_floor-1e-6:
            raise RuntimeError('Exported subset failed its optimized spacing floor.')
        info['noncrossing_subset'] = dict(**subset_selection,solver=subset_solver,spacing_search=subset_spacing_search,
            missing_pairs=sorted(set(range(1,audit['pair_count']+1))-set(subset_result.pair.astype(int))),
            length_m=float(subset_result.length_m.sum()),landing_crossings=int(subset_result.cross_n.sum()),
            exact_crossing_pairs=sum(r['intersects'] for r in subset_spacing),
            exact_min_distance_m=min((r['min_xy_m'] for r in subset_spacing),default=None))
        info['strict_noncrossing'] = traffic_strict_layout(pool,options,info['strict_pairs'],distances,
            crossings,gaps,nodes,headings,fafs,areas,audit,config,traffic,folder/'strict_noncrossing',scenario,
            seconds=args.subset_search_seconds)
        pool['selected'] = pool.index.isin(indices)
        pool.to_csv(folder/'candidate_shortlist.csv',index=False,encoding='utf-8-sig')
        info.update(selection,star_count=int((result.operation=='STAR').sum()),
                    sid_count=int((result.operation=='SID').sum()),length_m=float(result.length_m.sum()),
                    landing_crossings=int(result.cross_n.sum()),outside_routes=int((result.outside_m>0.05).sum()),
                    exact_crossing_pairs=sum(r['intersects'] for r in pair_rows),
                    exact_min_distance_m=min((r['min_xy_m'] for r in pair_rows),default=None),
                    median_pair_distance_m=float(np.median([r['min_xy_m'] for r in pair_rows])) if pair_rows else None,
                    below_spacing_target_pairs=sum(r['below_target'] for r in pair_rows),
                    roundtrip_and_analytic_boundary_validated=True)
        summary['scenarios'][scenario] = info
        tables[scenario] = result
        print(f'{scenario}: {info}',flush=True)
    if fingerprint(sources)!=hashes:
        raise RuntimeError('Input source hashes changed during execution.')
    summary['original_inputs_unchanged'] = True
    (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    write_design_review(output,summary)
    traffic_report(output,summary,tables)
    print(f'Complete: {output}',flush=True)


if __name__ == '__main__':
    main()
