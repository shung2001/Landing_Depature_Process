"""워크플로에 따른 연구용 최종 OLS 선정 (공식 승인/운항 안전 판정 아님).

실행: python -B 코드/polygon/최종_OLS_선정.py
검증: python -B 코드/polygon/최종_OLS_선정.py --self-test
기존 결과의 PNG만 갱신: python -B 코드/polygon/최종_OLS_선정.py --png-only

기존 Polygon을 건물·중심·정규화 3D 형상으로 중복 제거한다.
135도 이상인 PASS/PASS, PASS/REVIEW, REVIEW/REVIEW 쌍을 우선 탐색한다.
쌍이 없으면 조건부/부적합을 포함해 최선의 단일 후보를 반드시 배정한다.
기존 host_index를 후보 건물 소유관계로 사용하고 원본 버티포트 ID를 보존한다.
같은 판정 조합을 유지하는 선택지 안에서 건물 간 XY 중첩 면적 합을 공동 최소화한다.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import Transformer
from shapely.geometry import Polygon, Point, MultiPolygon

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / '자료/안전성/송파구'
TOL = 1e-6  # 수치적 접촉 공차(m), 운항 Separation 기준이 아님
AREA_TOL = 1e-6  # 수치적 면적 공차(m²)


def sample_ground(geometries, crs, dem):
    points = gpd.GeoSeries(geometries, crs=crs).representative_point().to_crs(dem.crs)
    values = []
    for value in dem.sample([(p.x, p.y) for p in points], indexes=1, masked=True):
        cell = value[0]
        values.append(np.nan if np.ma.is_masked(cell) else float(cell))
    return np.asarray(values)


def surface_plane(surface):
    xyz = np.asarray(surface.exterior.coords)[:-1]
    if xyz.shape != (3, 3) or not np.isfinite(xyz).all():
        raise ValueError('유한한 세 꼭짓점의 Polygon Z가 필요합니다.')
    origin = xyz[0, :2]
    coefficients = np.linalg.solve(np.column_stack([xyz[:, :2]-origin, np.ones(3)]), xyz[:, 2])
    return origin, coefficients


def evaluate(surface, footprints, roofs, tree):
    """교차 영역 모든 꼭짓점에서 평면고를 계산한다 (구멍/복합/선/점 포함)."""
    xy = shapely.force_2d(surface)
    origin, coefficients = surface_plane(surface)
    indices = np.asarray(tree.query(xy, predicate='intersects'), dtype=int)
    known = indices[np.isfinite(roofs[indices])]
    unknown = indices[~np.isfinite(roofs[indices])]
    clearances = []
    for index in known:
        overlap = xy.intersection(footprints[index])
        coords = shapely.get_coordinates(overlap)
        elevations = (coords-origin) @ coefficients[:2] + coefficients[2]
        clearances.append(float(elevations.min()-roofs[index]))
    values = np.asarray(clearances)
    penetrations = known[values < -TOL]
    contacts = known[np.abs(values) <= TOL]
    minimum = float(values.min()) if len(values) else None
    limiting = int(known[values.argmin()]) if len(values) else None
    return dict(result='FAIL' if len(penetrations) else 'REVIEW' if len(unknown)+len(contacts) else 'PASS',
                intersecting_buildings=len(indices), known_penetrations=len(penetrations),
                touching_buildings=len(contacts), unknown_buildings=len(unknown),
                min_clearance_m=minimum, limiting_building_row=limiting,
                unknown_building_rows=json.dumps(sorted(unknown.tolist())),
                penetration_building_rows=json.dumps(sorted(penetrations.tolist())),
                touching_building_rows=json.dumps(sorted(contacts.tolist())))


def consecutive_sectors(angles):
    """1도 샘플을 연결한다. 359→0 연속 및 360개 전체 PASS 처리."""
    accepted = {int(a) for a in angles}
    if not accepted:
        return []
    if accepted == set(range(360)):
        return [list(range(360))]
    starts = sorted(a for a in accepted if (a-1) % 360 not in accepted)
    groups = []
    for start in starts:
        group, angle = [], start
        while angle in accepted:
            group.append(angle)
            angle = (angle+1) % 360
        groups.append(group)
    return groups


def check_inputs(summary, candidate_dir):
    """서로 다른 ID는 같은 좌표여도 병합하지 않고 각 ID의 모든 후보 속성을 확인."""
    if summary.id.duplicated().any() or summary.empty:
        raise ValueError('요약표의 버티포트 ID는 고유해야 합니다.')
    numeric = summary[['longitude', 'latitude', 'agl_m', 'ground_z_m', 'apex_z_m', 'ols_slope_deg']].to_numpy(float)
    if not np.isfinite(numeric).all() or not (summary.agl_m > 0).all():
        raise ValueError('유효한 중심점과 시작 높이가 필요합니다.')
    if not np.allclose(summary.ols_slope_deg, 8) or not np.allclose(summary.apex_z_m, summary.ground_z_m+summary.agl_m):
        raise ValueError('현재 실행은 기존 8도/DEM+AGL 설계를 요구합니다.')
    for row in summary.itertuples():
        path = candidate_dir / str(int(row.id)) / f'Polygon_{int(row.id)}.shp'
        if not path.exists():
            if row.polygon_count == 0:
                continue
            raise FileNotFoundError(path)
        frame = gpd.read_file(path)
        expected = dict(heli_id=row.id, apex_lon=row.longitude, apex_lat=row.latitude,
                        apex_agl=row.agl_m, ground_z=row.ground_z_m, apex_z=row.apex_z_m,
                        height_m=1000, width_m=573, slope_deg=8)
        for field, value in expected.items():
            if not np.allclose(pd.to_numeric(frame[field], errors='raise'), value, rtol=0, atol=1e-6):
                raise ValueError(f'{path}: {field} 불일치')
        if not frame.z_ref.eq('DEM+AGL').all():
            raise ValueError(f'{path}: 높이 기준 불일치')


def add_tables(path, tables):
    """비공간 표도 GeoPackage attributes로 등록하여 QGIS에서 읽을 수 있게 한다."""
    with sqlite3.connect(path) as connection:
        for name, frame in tables.items():
            frame.to_sql(name, connection, if_exists='replace', index=False)
            connection.execute('DELETE FROM gpkg_contents WHERE table_name=?', (name,))
            connection.execute("INSERT INTO gpkg_contents(table_name,data_type,identifier,description,last_change) VALUES (?, 'attributes', ?, '', ?)",
                               (name, name, datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S.%fZ')))


def write_layer(records, path, layer, geometry_type):
    import pyogrio
    if records:
        frame = gpd.GeoDataFrame(records, crs=5186)
    else:
        frame = gpd.GeoDataFrame({'vp_id': pd.Series(dtype='int64'), 'geometry': gpd.GeoSeries([], crs=5186)}, crs=5186)
    pyogrio.write_dataframe(frame, path, layer=layer, driver='GPKG', geometry_type=geometry_type)


def self_test():
    # 중심점은 높은 쪽에 있어도 footprint 경계의 낮은 표면고로 침범을 잡아야 함.
    s = Polygon([(0, 0, 10), (10, 0, 20), (0, 10, 10)])
    b = np.array([Polygon([(1, 1), (8, 1), (8, 2), (1, 2)])], dtype=object)
    result = evaluate(s, b, np.array([12.]), shapely.STRtree(b))
    assert result['result'] == 'FAIL' and np.isclose(result['min_clearance_m'], -1)
    assert evaluate(s, b, np.array([11.]), shapely.STRtree(b))['touching_buildings'] == 1
    assert evaluate(s, b, np.array([np.nan]), shapely.STRtree(b))['result'] == 'REVIEW'
    assert evaluate(s, b, np.array([9.]), shapely.STRtree(b))['result'] == 'PASS'
    far = np.array([Point(50, 50).buffer(1)], dtype=object)
    r = evaluate(s, far, np.array([np.nan]), shapely.STRtree(far))
    assert r['result'] == 'PASS' and r['min_clearance_m'] is None
    assert consecutive_sectors([358, 359, 0, 1, 20]) == [[20], [358, 359, 0, 1]]
    assert consecutive_sectors([]) == [] and len(consecutive_sectors(range(360))[0]) == 360
    # 면적이 아닌 꼭짓점만의 접촉도 평가한다.
    touch = np.array([Polygon([(10, 0), (11, 0), (11, 1)])], dtype=object)
    assert evaluate(s, touch, np.array([20.]), shapely.STRtree(touch))['result'] == 'REVIEW'
    print('SELF-TEST OK: 부분 교차, 접촉, 결측, 비중첩, 순환 sector, 점 접촉')



def angle_distance(a, b):
    return abs((a-b+180) % 360-180)


def geometry_key(record, include_owner=True):
    """링 시작점/방향을 정규화하되 좌표 반올림은 하지 않는다. Z 차이도 보존."""
    geom = record['geometry']
    definition = (record['center_x'], record['center_y'], record['z0_dem_agl_m'],
                  record['slope_deg'], record['length_m'], record['start_width_m'],
                  record['end_width_m'], record['surface_model_id'])
    return ((record['building_id'],) if include_owner else ()) + definition + (
        shapely.to_wkb(shapely.normalize(geom), output_dimension=3),)


def deduplicate(records):
    unique, mapping, keys = [], [], {}
    for record in records:
        key = geometry_key(record)
        duplicate = key in keys
        if not duplicate:
            keys[key] = record
            unique.append(record.copy())
        representative = keys[key]
        mapping.append(dict(building_id=record['building_id'], source_vp_id=record['source_vp_id'],
            original_polygon_id=record['polygon_id'], representative_polygon_id=representative['polygon_id'],
            is_duplicate=duplicate, source_file=record['source_file'], source_row=record['source_row']))
    return unique, mapping


def clearance_value(record):
    value = record['min_clearance_m']
    return float('inf') if value is None or pd.isna(value) else float(value)


def mark_near_duplicates(records, tolerance=0.01):
    """같은 설계의 1cm 이내 유사 형상은 플래그만 붙이고 제거하지 않는다."""
    matches = {r['polygon_id']: [] for r in records}
    groups = {}
    for record in records:
        key = tuple(record[k] for k in ('building_id', 'slope_deg', 'length_m',
            'start_width_m', 'end_width_m', 'surface_model_id'))
        groups.setdefault(key, []).append(record)
    for group in groups.values():
        geometries = np.array([r['geometry'] for r in group], dtype=object)
        for i, a in enumerate(group[:-1]):
            distances = shapely.hausdorff_distance(a['geometry'], geometries[i+1:])
            for j in np.flatnonzero(distances <= tolerance) + i + 1:
                b = group[int(j)]
                if abs(a['z0_dem_agl_m'] - b['z0_dem_agl_m']) <= tolerance and np.hypot(
                    a['center_x']-b['center_x'], a['center_y']-b['center_y']) <= tolerance:
                    matches[a['polygon_id']].append(b['polygon_id'])
                    matches[b['polygon_id']].append(a['polygon_id'])
    for record in records:
        record['near_duplicate_flag'] = bool(matches[record['polygon_id']])
        record['near_duplicate_polygon_ids'] = json.dumps(matches[record['polygon_id']])


def single_key(record):
    """작은 키 우선. PASS > REVIEW > FAIL, 마지막 동률은 방위각과 ID로 결정."""
    common = (record['azimuth_deg'], record['polygon_id'])
    if record['status'] == 'PASS':
        return (0, -clearance_value(record), *common)
    if record['status'] == 'REVIEW':
        return (1, record['unknown_count'], record['touching_buildings'],
                -clearance_value(record), *common)
    return (2, record['worst_penetration_m'], record['penetration_count'],
            record['unknown_count'], *common)


def choose_candidates(records):
    """단일 최적을 먼저 고정하지 않고 135도 이상 모든 조합 탐색. FAIL 포함 쌍 제외."""
    eligible = [r for r in records if r['status'] != 'FAIL']
    best, best_key = None, None
    for i, a in enumerate(eligible):
        for b in eligible[i+1:]:
            diff = angle_distance(a['azimuth_deg'], b['azimuth_deg'])
            if diff < 135:
                continue
            passes = int(a['status']=='PASS') + int(b['status']=='PASS')
            clearances = [clearance_value(a), clearance_value(b)]
            key = (2-passes, a['unknown_count']+b['unknown_count'],
                   a['touching_buildings']+b['touching_buildings'],
                   -min(clearances), -max(clearances), -diff,
                   min(a['azimuth_deg'], b['azimuth_deg']), max(a['azimuth_deg'], b['azimuth_deg']))
            if best_key is None or key < best_key:
                best, best_key = [a, b], key
    if best:
        return sorted(best, key=single_key), 'SELECTED'
    return ([min(records, key=single_key)] if records else []), 'NOT_FOUND'


def overlap_matrix(records):
    """서로 다른 소유자 사이의 모든 후보 쌍을 정확한 XY 교차 면적으로 비교."""
    geometries = shapely.force_2d(np.array([r['geometry'] for r in records], dtype=object))
    matrix = np.zeros((len(records), len(records)), dtype=float)
    if not len(records):
        return matrix
    owners = np.array([r['building_id'] for r in records])
    left, right = shapely.STRtree(geometries).query(geometries, predicate='intersects')
    keep = (left < right) & (owners[left] != owners[right])
    left, right = left[keep], right[keep]
    for start in range(0, len(left), 20000):
        a, b = left[start:start+20000], right[start:start+20000]
        areas = shapely.area(shapely.intersection(geometries[a], geometries[b]))
        areas[areas <= AREA_TOL] = 0
        matrix[a, b] = areas
        matrix[b, a] = areas
        print(f'중첩 면적 계산 {min(start+20000,len(left))}/{len(left)}', flush=True)
    return matrix


def admissible_bundles(records, owners):
    """판정 조합을 건물별 최상위 등급으로 고정한다. 기존 품질 순서로 정렬."""
    groups, baseline = {}, {}
    by_polygon = {r['polygon_id']: i for i, r in enumerate(records)}
    for owner in owners:
        building_id = owner['building_id']
        indices = [i for i,r in enumerate(records) if r['building_id']==building_id]
        options = [records[i] for i in indices]
        if not options:
            continue
        original, _ = choose_candidates(options)
        baseline[building_id] = tuple(sorted(by_polygon[r['polygon_id']] for r in original))
        pairs, best_category = [], 3
        eligible = [i for i in indices if records[i]['status']!='FAIL']
        for position, i in enumerate(eligible):
            a = records[i]
            for j in eligible[position+1:]:
                b = records[j]
                diff = angle_distance(a['azimuth_deg'],b['azimuth_deg'])
                if diff < 135:
                    continue
                category = 2-int(a['status']=='PASS')-int(b['status']=='PASS')
                if category < best_category:
                    pairs, best_category = [], category
                if category == best_category:
                    clearances = [clearance_value(a),clearance_value(b)]
                    quality = (a['unknown_count']+b['unknown_count'],
                        a['touching_buildings']+b['touching_buildings'],
                        -min(clearances),-max(clearances),-diff,
                        min(a['azimuth_deg'],b['azimuth_deg']),max(a['azimuth_deg'],b['azimuth_deg']))
                    pairs.append((quality,(i,j)))
        if pairs:
            groups[building_id] = [pair for _,pair in sorted(pairs)]
        else:
            best_status = min({'PASS':0,'REVIEW':1,'FAIL':2}[records[i]['status']] for i in indices)
            singles = [i for i in indices if {'PASS':0,'REVIEW':1,'FAIL':2}[records[i]['status']]==best_status]
            if best_status == 2:
                # 부적합 후보에서는 침범량·침범수 최소 조건까지 고정한 뒤 중첩을 비교한다.
                best_failure = min((records[i]['worst_penetration_m'],records[i]['penetration_count']) for i in singles)
                singles = [i for i in singles if (records[i]['worst_penetration_m'],records[i]['penetration_count'])==best_failure]
            groups[building_id] = [(i,) for i in sorted(singles,key=lambda i:single_key(records[i]))]
        assert baseline[building_id] in groups[building_id]
    return groups, baseline


def total_overlap(selection, matrix):
    indices = [i for bundle in selection.values() for i in bundle]
    return float(matrix[np.ix_(indices,indices)].sum()/2)


def optimize_overlap(groups, baseline, matrix, starts=16, max_sweeps=60, seed=20261016):
    """전체 조합을 유지한 다중 시작점 좌표 하강. 전역 최적해 보장은 하지 않는다."""
    rng = np.random.default_rng(seed)
    names = sorted(groups)
    arrays = {name: np.asarray(groups[name],dtype=int) for name in names}
    best = dict(baseline)
    best_score = total_overlap(best,matrix)
    history = []
    def quality(selection):
        return tuple(groups[name].index(selection[name]) for name in names)
    best_quality = quality(best)
    for start in range(starts):
        current = dict(baseline) if start==0 else {
            name:groups[name][int(rng.integers(len(groups[name])))] for name in names}
        history.append(dict(start=start,sweep=0,overlap_sum_m2=total_overlap(current,matrix),changed_buildings=0))
        converged = False
        for sweep in range(1,max_sweeps+1):
            changed = 0
            for name in (names if start==0 else list(rng.permutation(names))):
                others = [i for key,bundle in current.items() if key!=name for i in bundle]
                cost = matrix[:,others].sum(axis=1) if others else np.zeros(len(matrix))
                values = cost[arrays[name]].sum(axis=1)
                # 면적 최솟값을 고른다. 같은 면적이면 기존 품질순으로 앞선 조합 우선.
                chosen = groups[name][int(np.argmin(values))]
                if chosen != current[name]:
                    current[name] = chosen
                    changed += 1
            score = total_overlap(current,matrix)
            history.append(dict(start=start,sweep=sweep,overlap_sum_m2=score,changed_buildings=changed))
            current_quality = quality(current)
            if score < best_score or (score==best_score and current_quality < best_quality):
                best,best_score,best_quality = dict(current),score,current_quality
            if not changed:
                converged = True
                break
        print(f'중첩 공동선정 {start+1}/{starts}: {score:,.3f}m² / 최선 {best_score:,.3f}m² / 수렴 {converged}',flush=True)
    # 같은 총면적에서 각 건물의 기존 품질순으로 정렬된 첫 최소값을 고르는 결정적 최종 수렴.
    final_converged = False
    for sweep in range(max_sweeps):
        changed = 0
        for name in names:
            others = [i for key,bundle in best.items() if key!=name for i in bundle]
            cost = matrix[:,others].sum(axis=1) if others else np.zeros(len(matrix))
            chosen = groups[name][int(np.argmin(cost[arrays[name]].sum(axis=1)))]
            if chosen != best[name]:
                best[name] = chosen
                changed += 1
        if not changed:
            final_converged = True
            break
    final_score = total_overlap(best,matrix)
    assert final_score <= total_overlap(baseline,matrix)+AREA_TOL
    if final_converged:
        # 다른 건물을 고정한 상태에서 단일 건물의 쌍 전체를 바꿔도 개선할 수 없음.
        for name in names:
            others = [i for key,bundle in best.items() if key!=name for i in bundle]
            cost = matrix[:,others].sum(axis=1) if others else np.zeros(len(matrix))
            assert cost[list(best[name])].sum() <= cost[arrays[name]].sum(axis=1).min()+AREA_TOL
    return best,history,dict(starts=starts,max_sweeps=max_sweeps,seed=seed,
        local_optimum_verified=final_converged,baseline_overlap_sum_m2=total_overlap(baseline,matrix),
        optimized_overlap_sum_m2=final_score)


def surface_overlap_details(a,b):
    """두 선형 표면의 XY 교차 영역과 그 위의 부호 있는 높이 차 최솟값/최댓값."""
    intersection = shapely.force_2d(a['geometry']).intersection(shapely.force_2d(b['geometry']))
    if intersection.area <= AREA_TOL:
        return None
    coords = shapely.get_coordinates(intersection)
    oa,ca = surface_plane(a['geometry'])
    ob,cb = surface_plane(b['geometry'])
    delta = (coords-oa)@ca[:2]+ca[2]-((coords-ob)@cb[:2]+cb[2])
    low,high = float(delta.min()),float(delta.max())
    crossing = low<=TOL and high>=-TOL
    return dict(overlap_area_m2=float(intersection.area),z_a_minus_b_min_m=low,z_a_minus_b_max_m=high,
        min_vertical_gap_m=0. if crossing else min(abs(low),abs(high)),
        surfaces_cross_or_touch=crossing,
        geometry=MultiPolygon([intersection]) if intersection.geom_type=='Polygon' else intersection)


def overlap_audit(records):
    pairs = []
    for i,a in enumerate(records):
        for b in records[i+1:]:
            if a['building_id']==b['building_id']:
                continue
            details = surface_overlap_details(a,b)
            if details:
                pairs.append(dict(building_a=a['building_id'],polygon_a=a['polygon_id'],azimuth_a=a['azimuth_deg'],
                    building_b=b['building_id'],polygon_b=b['polygon_id'],azimuth_b=b['azimuth_deg'],**details))
    stats = dict(overlap_sum_m2=sum(r['overlap_area_m2'] for r in pairs),overlap_pair_count=len(pairs),
        crossing_pair_count=sum(r['surfaces_cross_or_touch'] for r in pairs),
        overlap_union_m2=float(shapely.union_all([r['geometry'] for r in pairs]).area) if pairs else 0.)
    return pairs,stats


def workflow_tests():
    self_test()
    def item(angle, status='PASS', clearance=10, unknown=0, worst=0, pid=None):
        return dict(azimuth_deg=angle, status=status, min_clearance_m=clearance,
                    unknown_count=unknown, touching_buildings=0, worst_penetration_m=worst,
                    penetration_count=int(worst>0), polygon_id=pid or str(angle))
    # 0도 최고 단일 후보에 고정하면 쌍을 잃는다. 90/270을 우선 선택해야 한다.
    selected, pair = choose_candidates([item(0, clearance=100), item(90), item(270)])
    assert pair=='SELECTED' and {r['azimuth_deg'] for r in selected}=={90,270}
    assert choose_candidates([item(10), item(350)])[1]=='NOT_FOUND'
    assert angle_distance(227,90)==137 and angle_distance(30,165)==135
    selected, _ = choose_candidates([item(0), item(180,'REVIEW',unknown=2)])
    assert len(selected)==2 and selected[0]['status']=='PASS'
    selected, _ = choose_candidates([item(0,'REVIEW',unknown=1), item(180,'REVIEW',unknown=1)])
    assert len(selected)==2
    selected, pair = choose_candidates([item(0,'FAIL',worst=20), item(180,'FAIL',worst=5)])
    assert len(selected)==1 and selected[0]['worst_penetration_m']==5 and pair=='NOT_FOUND'
    shape = Polygon([(0,0,10),(1,0,11),(0,1,11)])
    base = dict(building_id='A', center_x=0., center_y=0., z0_dem_agl_m=10.,
        slope_deg=8., length_m=1000., start_width_m=0., end_width_m=573.,
        surface_model_id='TEST', source_vp_id=1, source_file='test', source_row=0)
    a = dict(base, polygon_id='a', geometry=shape)
    b = dict(base, polygon_id='b', geometry=Polygon(list(shape.exterior.coords)[::-1]))
    c = dict(base, polygon_id='c', geometry=Polygon([(0,0,11),(1,0,12),(0,1,12)]))
    d = dict(a, polygon_id='d', building_id='B')
    unique, mapping = deduplicate([a,b,c,d])
    assert len(unique)==3 and mapping[1]['representative_polygon_id']=='a'
    assert geometry_key(a, False)==geometry_key(d, False)
    near = dict(a, polygon_id='near', geometry=Polygon([(0.001,0,10),(1.001,0,11),(0.001,1,11)]))
    mark_near_duplicates([a, near])
    assert a['near_duplicate_flag'] and near['near_duplicate_flag']
    # 높이가 다른 XY 중첩과 실제 표면 교차를 구분한다.
    low = dict(building_id='A',polygon_id='low',azimuth_deg=0,
        geometry=Polygon([(0,0,10),(10,0,10),(0,10,10)]))
    high = dict(building_id='B',polygon_id='high',azimuth_deg=180,
        geometry=Polygon([(0,0,20),(10,0,20),(0,10,20)]))
    audit=surface_overlap_details(low,high)
    assert np.isclose(audit['overlap_area_m2'],50) and audit['min_vertical_gap_m']==10
    assert not audit['surfaces_cross_or_touch']
    crossing=dict(high,geometry=Polygon([(0,0,5),(10,0,15),(0,10,5)]))
    assert surface_overlap_details(low,crossing)['surfaces_cross_or_touch']
    # 공동선정이 품질 등급/135도 조건을 지키며 실제 중첩을 줄이는지 검증.
    def opt_record(owner,pid,angle,geom,status='PASS'):
        return dict(item(angle,status,pid=pid),building_id=owner,geometry=geom)
    same=low['geometry']
    away=Polygon([(100,100,10),(110,100,10),(100,110,10)])
    optimization_records=[opt_record('A','a0',0,same),opt_record('A','a180',180,same),
        opt_record('B','b0',0,same),opt_record('B','b180',180,same),
        opt_record('B','b10',10,away),opt_record('B','b190',190,away)]
    groups,baseline=admissible_bundles(optimization_records,[{'building_id':'A'},{'building_id':'B'}])
    matrix=overlap_matrix(optimization_records)
    solution,_,info=optimize_overlap(groups,baseline,matrix,starts=2,max_sweeps=10)
    assert info['baseline_overlap_sum_m2']>0 and total_overlap(solution,matrix)==0
    for bundle in solution.values():
        assert len(bundle)==2 and angle_distance(*(optimization_records[i]['azimuth_deg'] for i in bundle))>=135
    optimization_records[4]['status']='REVIEW'
    optimization_records[5]['status']='REVIEW'
    groups,_=admissible_bundles(optimization_records,[{'building_id':'B'}])
    assert all(all(optimization_records[i]['status']=='PASS' for i in bundle) for bundle in groups['B'])
    # 세 표면의 중첩을 쌍 수(3쌍)가 아닌 실제 표면 수(3개)로 분류한다.
    count_records=[dict(building_id=name,geometry=shapely.box(x,0,x+3,2))
                   for name,x in [('A',0),('B',1),('C',2)]]
    regions=overlap_regions(count_records)
    assert np.isclose(sum(r['geometry'].area for r in regions if r['ols_count']==2),4.)
    assert np.isclose(sum(r['geometry'].area for r in regions if r['ols_count']==3),2.)
    assert not overlap_regions([dict(r,building_id='A') for r in count_records])
    assert not overlap_regions([])
    print('WORKFLOW TEST OK: 쌍 우선, 135도/순환각, 조건부/부적합, 링 방향, Z/소유관계 보존')
    print('OVERLAP TEST OK: XY/높이 교차 구분, 공동 면적 감소, 품질등급 보존')
    print('OVERLAP COUNT TEST OK: 2중/3중 영역 면적, 동일 건물 제외, 빈 입력')


def surface_altitudes(record):
    """속성 대신 실제 Polygon Z에서 시작·끝 고도를 읽고 일관성을 검증한다."""
    xyz=np.asarray(record['geometry'].exterior.coords)[:-1]
    if xyz.shape!=(3,3) or not np.isfinite(xyz).all():
        raise ValueError('고도 표시에는 유효한 삼각형 Polygon Z가 필요합니다.')
    start_index=int(np.argmin(xyz[:,2]))
    start=float(xyz[start_index,2])
    end_values=np.delete(xyz,start_index,axis=0)[:,2]
    if not np.allclose(end_values,end_values[0],rtol=0,atol=1e-7):
        raise ValueError('OLS 끝단의 고도가 일정하지 않습니다.')
    if not np.isclose(start,record['z0_dem_agl_m'],rtol=0,atol=1e-7):
        raise ValueError('Polygon Z와 시작고도 속성이 다릅니다.')
    return start,float(end_values[0])


def source_ids(record):
    """그림 식별자는 원본 헬리포트의 id 필드를 사용한다. 새 번호를 만들지 않는다."""
    values=record.get('source_vp_ids')
    if isinstance(values,str):
        values=json.loads(values)
    if values is None:
        values=[record['source_vp_id']]
    return sorted({int(value) for value in values})


def owner_label(record,separator=' / '):
    return separator.join(f'id_{value}' for value in source_ids(record))


def agl_labels(records,dem_path):
    """시작점과 끝단 두 꼭짓점의 현지 지표고를 빼 AGL 표시값을 구한다.

    끝단 범위는 두 꼭짓점의 값이며, 끝단 변 전체의 최솟값/최댓값이 아니다.
    Polygon Z는 그대로 두고 표시만 변환한다. NoData는 임의 대체하지 않는다.
    """
    import rasterio
    labels={}
    with rasterio.open(dem_path) as dem:
        transform=Transformer.from_crs(5186,dem.crs,always_xy=True)
        for record in records:
            surface_altitudes(record)
            xyz=np.asarray(record['geometry'].exterior.coords)[:-1]
            x,y=transform.transform(xyz[:,0],xyz[:,1])
            values=list(dem.sample(zip(x,y),indexes=1,masked=True))
            ground=np.array([np.nan if np.ma.is_masked(value[0]) else float(value[0]) for value in values])
            agl=xyz[:,2]-ground
            start_index=int(np.argmin(xyz[:,2]))
            start=agl[start_index]
            ends=np.delete(agl,start_index)
            start_text=f'{start:.2f}' if np.isfinite(start) else 'DEM 결측'
            if np.isfinite(ends).all():
                low,high=float(ends.min()),float(ends.max())
                end_text=f'{low:.2f}' if round(low,2)==round(high,2) else f'{low:.2f}~{high:.2f}'
            else:
                end_text='DEM 결측'
            labels[record['polygon_id']]=(start_text,end_text)
    return labels


def save_map(fig,out,stem,formats):
    from io import BytesIO
    for extension in formats:
        if extension not in ('png','pdf'):
            raise ValueError(f'지원하지 않는 지도 형식: {extension}')
        # 메모리에서 렌더링을 완료한 뒤 기록한다. Windows 기존 PNG의 w+b 열기 오류 방지.
        buffer=BytesIO()
        fig.savefig(buffer,format=extension,dpi=180)
        target=out/f'{stem}.{extension}'
        with target.open('r+b' if target.exists() else 'wb') as stream:
            stream.write(buffer.getbuffer())
            stream.truncate()


def draw_maps(selected, owners, obstacles, out, formats=('png','pdf'),dem_path=DATA/'송파구_DEM.tif'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False})
    heights=agl_labels(selected,dem_path)
    row_count=max(2,int(np.ceil(len(owners)/4)))
    fig,axes=plt.subplots(row_count,4,figsize=(24,7.5*row_count),squeeze=False)
    for ax,owner in zip(axes.flat,sorted(owners,key=lambda r:source_ids(r))):
        records=sorted([r for r in selected if r['building_id']==owner['building_id']],
                       key=lambda r:r['selected_rank'])
        cx,cy=owner['center_x'],owner['center_y']
        geometries=[shapely.force_2d(r['geometry']) for r in records]
        bounds=shapely.union_all(geometries).bounds if geometries else (cx-500,cy-500,cx+500,cy+500)
        xmin,ymin,xmax,ymax=bounds
        pad=max(xmax-xmin,ymax-ymin,200)*.08
        local=obstacles.cx[xmin-pad:xmax+pad,ymin-pad:ymax+pad]
        if len(local):
            local.plot(ax=ax,color='#e0e3e5',edgecolor='#c6cbce',linewidth=.2)
        for record in records:
            rank=record['selected_rank']
            color='#167a93' if rank==1 else '#d3861d'
            if record['status']=='FAIL':
                color='#bb2637'
            style='--' if record['status']=='REVIEW' else '-'
            x,y=shapely.force_2d(record['geometry']).exterior.xy
            start,end=heights[record['polygon_id']]
            ax.fill(x,y,color=color,alpha=.13)
            ax.plot(x,y,color=color,linestyle=style,linewidth=1.7,
                label=f"{rank}순위 {record['azimuth_deg']:g}° / {record['status']}\n시작 AGL {start} m\n끝단 양끝 AGL {end} m")
        ax.scatter([cx],[cy],marker='*',s=65,color='black',zorder=5)
        ax.set_xlim(xmin-pad,xmax+pad);ax.set_ylim(ymin-pad,ymax+pad)
        ax.set_aspect('equal');ax.ticklabel_format(style='plain',useOffset=False)
        ax.tick_params(labelsize=8);ax.grid(alpha=.12)
        ax.set_title(owner_label(owner),fontsize=13,pad=10)
        if records:
            ax.legend(fontsize=9,loc='best',framealpha=.95,labelspacing=.8)
    for ax in list(axes.flat)[len(owners):]:
        ax.set_visible(False)
    fig.suptitle('건물 id별 최종 OLS 방향과 지상고 AGL (m)\n1순위 청록 / 2순위 주황 · REVIEW 점선 · 별표: 시작점 · id_n: 원본 id\nAGL = Polygon Z - 해당 좌표의 DEM · 끝단 범위는 두 꼭짓점의 값',fontsize=17)
    fig.tight_layout(rect=(0,0,1,.92),h_pad=3,w_pad=2)
    save_map(fig,out,'ols_selected_map',formats)
    plt.close(fig)


def draw_overlap_comparison(baseline,selected,comparison,obstacles,out,formats=('png','pdf'),dem_path=DATA/'송파구_DEM.tif'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False})
    heights=agl_labels(selected,dem_path)
    table_rows=max(1,len({r['building_id'] for r in selected}))
    table_height=max(3.4,.34*(table_rows+1))
    fig=plt.figure(figsize=(20,11.6+table_height))
    grid=GridSpec(2,2,figure=fig,height_ratios=[9.6,table_height],hspace=.23,wspace=.17)
    axes=[fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1])]
    all_xy=[shapely.force_2d(r['geometry']) for r in baseline+selected]
    # 모든 건물이 NO_GEOMETRY여도 빈 결과 지도를 생성할 수 있게 한다.
    bounds=shapely.union_all(all_xy).bounds if all_xy else obstacles.total_bounds
    xmin,ymin,xmax,ymax=bounds
    pad=max(xmax-xmin,ymax-ymin)*.04
    local=obstacles.cx[xmin-pad:xmax+pad,ymin-pad:ymax+pad]
    ids=sorted({r['building_id'] for r in selected})
    labels={name:next(r for r in selected if r['building_id']==name) for name in ids}
    ids.sort(key=lambda name:source_ids(labels[name]))
    colors={name:plt.get_cmap('tab10')(i%10) for i,name in enumerate(ids)}
    for ax,records,stats,title in zip(axes,[baseline,selected],comparison,
            ['기존 독립 선정','건물 간 중첩 공동 최소화 / 최종 OLS']):
        if len(local):
            local.plot(ax=ax,color='#edf0f1',linewidth=0)
        for record in records:
            x,y=shapely.force_2d(record['geometry']).exterior.xy
            ax.fill(x,y,color=colors[record['building_id']],alpha=.12)
            ax.plot(x,y,color=colors[record['building_id']],linewidth=.85,
                linestyle='--' if record['status']=='REVIEW' else '-')
        overlap_pairs,_=overlap_audit(records)
        if overlap_pairs:
            gpd.GeoSeries([r['geometry'] for r in overlap_pairs],crs=5186).plot(
                ax=ax,color='#cf283b',alpha=.5,linewidth=0)
        for name in ids:
            record=next((r for r in records if r['building_id']==name),None)
            if record is not None:
                ax.annotate(owner_label(record,separator='\n'),(record['center_x'],record['center_y']),
                    xytext=(5,6),textcoords='offset points',fontsize=7.5,
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.75),zorder=10)
        ax.set_xlim(xmin-pad,xmax+pad);ax.set_ylim(ymin-pad,ymax+pad);ax.set_aspect('equal')
        ax.ticklabel_format(style='plain',useOffset=False);ax.grid(alpha=.15)
        ax.set_title(f"{title}\n쌍별 중첩 면적 합 {stats['overlap_sum_m2']:,.1f}m² / 중첩 {stats['overlap_pair_count']}쌍",fontsize=13)
        ax.set_xlabel('Easting (m)');ax.set_ylabel('Northing (m)')
    ax_table=fig.add_subplot(grid[1,:]);ax_table.axis('off')
    rows=[]
    for owner in ids:
        records=sorted([r for r in selected if r['building_id']==owner],key=lambda r:r['selected_rank'])
        text=[owner_label(records[0]),str(source_ids(records[0]))]
        for rank in (1,2):
            record=next((r for r in records if r['selected_rank']==rank),None)
            if record is None:
                text.extend(['—','—'])
            else:
                start,end=heights[record['polygon_id']]
                text.extend([f"{record['azimuth_deg']:g}° / {record['status']}",f'{start} → {end}'])
        rows.append(text)
    if rows:
        table=ax_table.table(cellText=rows,colLabels=['후보 건물 id','원본 id 값',
            '1순위 방위각 / 판정','1순위 시작 → 끝단 AGL (m)','2순위 방위각 / 판정','2순위 시작 → 끝단 AGL (m)'],
            cellLoc='center',colWidths=[.13,.13,.18,.19,.18,.19],loc='center',bbox=[0,.04,1,.90])
        table.auto_set_font_size(False);table.set_fontsize(10)
        for (r,c),cell in table.get_celld().items():
            cell.set_edgecolor('#d1d7dd');cell.set_linewidth(.5)
            if r==0:
                cell.set_facecolor('#e9f0f4');cell.set_text_props(weight='bold')
            elif r%2==0:
                cell.set_facecolor('#f7f9fa')
    ax_table.set_title('최종 OLS 지상고 AGL (m): 각 좌표의 Polygon Z - DEM · 끝단 범위는 두 꼭짓점의 값',
        fontsize=12,pad=8)
    fig.suptitle('건물 id별 최종 OLS 배치 및 지상고 AGL\n빨강: 다른 건물 표면과의 XY 중첩 / 점선: REVIEW / 같은 위치의 원본 id는 함께 표기',
        fontsize=17,y=.98)
    fig.subplots_adjust(top=.91,bottom=.04,left=.05,right=.98)
    save_map(fig,out,'ols_overlap_comparison',formats)
    plt.close(fig)


def overlap_regions(selected):
    """XY 경계를 분할하여 각 영역의 실제 OLS 수와 소유 건물 수를 계산한다."""
    geometries=[shapely.force_2d(r['geometry']) for r in selected]
    if not geometries:
        return []
    from shapely.ops import polygonize
    tree=shapely.STRtree(geometries)
    faces=polygonize(shapely.union_all([g.boundary for g in geometries]))
    regions=[]
    for face in faces:
        if face.area <= AREA_TOL:
            continue
        point=face.representative_point()
        indices=tree.query(point,predicate='intersects')
        owners={selected[int(i)]['building_id'] for i in indices}
        # 최적화 목적함수와 동일하게 서로 다른 소유 건물 간 중첩만 집계한다.
        if len(owners)<2:
            continue
        regions.append(dict(geometry=face,ols_count=len(indices),owner_count=len(owners)))
    return regions


def draw_overlap_counts(selected,obstacles,out):
    """최소화된 최종 배치와 영역별 OLS 중첩 개수만 PNG로 표시한다."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False})
    regions=overlap_regions(selected)
    counts=sorted({r['ols_count'] for r in regions})
    colors={n:plt.get_cmap('YlOrRd')(.35+.6*i/max(1,len(counts)-1)) for i,n in enumerate(counts)}
    fig,ax=plt.subplots(figsize=(13,14))
    geometries=[shapely.force_2d(r['geometry']) for r in selected]
    xmin,ymin,xmax,ymax=shapely.union_all(geometries).bounds if geometries else obstacles.total_bounds
    pad=max(xmax-xmin,ymax-ymin,200)*.06
    local=obstacles.cx[xmin-pad:xmax+pad,ymin-pad:ymax+pad]
    if len(local):
        local.plot(ax=ax,color='#e8ecef',linewidth=0,zorder=0)
    if geometries:
        gpd.GeoSeries(geometries,crs=5186).plot(ax=ax,facecolor='#e2f0f5',edgecolor='#5d8295',linewidth=.7,zorder=1)
    for n in counts:
        gpd.GeoSeries([r['geometry'] for r in regions if r['ols_count']==n],crs=5186).plot(
            ax=ax,color=colors[n],edgecolor='#a53b28',linewidth=.45,zorder=2)
    # 같은 개수의 인접 영역을 합쳐 표기 수를 줄인다. 분리된 영역은 각각 표시한다.
    for n in counts:
        merged=shapely.union_all([r['geometry'] for r in regions if r['ols_count']==n])
        for part in shapely.get_parts(merged):
            p=part.representative_point()
            ax.annotate(f'{n}개',(p.x,p.y),ha='center',va='center',fontsize=11,weight='bold',
                bbox=dict(boxstyle='round,pad=.15',facecolor='white',edgecolor='#9b4938',alpha=.9),zorder=4)
    owners={r['building_id']:r for r in selected}
    for record in owners.values():
        ax.scatter(record['center_x'],record['center_y'],s=30,c='#213e50',marker='*',zorder=5)
        ax.annotate(owner_label(record),(record['center_x'],record['center_y']),xytext=(5,7),
            textcoords='offset points',fontsize=8,bbox=dict(facecolor='white',edgecolor='none',alpha=.8),zorder=5)
    ax.set_xlim(xmin-pad,xmax+pad);ax.set_ylim(ymin-pad,ymax+pad);ax.set_aspect('equal')
    ax.ticklabel_format(style='plain',useOffset=False);ax.tick_params(labelsize=8);ax.grid(alpha=.15)
    ax.set_xlabel('Easting (m) · EPSG:5186');ax.set_ylabel('Northing (m)')
    ax.set_title('색칠된 영역: 서로 다른 건물의 OLS 중첩\n영역 안 숫자: 해당 위치에 겹치는 OLS 개수',fontsize=13,pad=14)
    if counts:
        ax.legend(handles=[Patch(facecolor=colors[n],edgecolor='#a53b28',label=f'OLS {n}개 중첩') for n in counts],
            loc='upper right',framealpha=.95,fontsize=11)
    else:
        ax.text(.5,.5,'면적이 있는 중첩 영역 없음',transform=ax.transAxes,ha='center',fontsize=13)
    fig.suptitle('중첩 면적 최소화 후 최종 OLS',fontsize=22,weight='bold',y=.975)
    fig.text(.5,.025,'서로 다른 건물 간 XY 중첩 기준 · 영역 위 숫자는 해당 위치에 겹치는 OLS 개수',
        ha='center',fontsize=11,color='#444444')
    fig.subplots_adjust(top=.90,bottom=.08,left=.09,right=.97)
    save_map(fig,out,'중첩',('png',))
    plt.close(fig)
    print(f'중첩 시각화 저장: {out / "중첩.png"}',flush=True)


def redraw_existing_pngs(out,building_path,dem_path=DATA/'송파구_DEM.tif'):
    """현재 선정 결과를 읽어서 PNG 세 개만 갱신한다. 최적화는 재실행하지 않는다."""
    source=out/'ols_selected_surfaces.gpkg'
    selected=gpd.read_file(source,layer='ols_selected_surfaces').to_dict('records')
    baseline=gpd.read_file(source,layer='ols_independent_baseline').to_dict('records')
    owners=gpd.read_file(out/'ols_candidates_unique.gpkg',layer='candidate_buildings')
    owners=owners.sort_values('host_source_row').to_dict('records')
    comparison=pd.read_csv(out/'ols_overlap_summary.csv').set_index('stage')
    comparison=comparison.loc[['independent_baseline','overlap_optimized']].to_dict('records')
    obstacles=gpd.read_file(building_path).to_crs(5186)
    draw_maps(selected,owners,obstacles,out,formats=('png',),dem_path=dem_path)
    draw_overlap_comparison(baseline,selected,comparison,obstacles,out,formats=('png',),dem_path=dem_path)
    draw_overlap_counts(selected,obstacles,out)
    print(f'PNG 3개 갱신 완료: {out}')


def main():
    import shutil
    import zipfile
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--png-only',action='store_true',help='현재 선정 결과로 고도·중첩 PNG 3개만 갱신')
    parser.add_argument('--overlap-only',action='store_true',help='현재 선정 결과로 중첩.png만 생성 (최적화 재실행 없음)')
    parser.add_argument('--buildings', type=Path, default=DATA/'송파구_건물_수정_update.shp')
    parser.add_argument('--dem', type=Path, default=DATA/'송파구_DEM.tif')
    parser.add_argument('--candidates', type=Path, default=ROOT/'자료/결과물/Polygon/8')
    parser.add_argument('--output', type=Path, default=DATA/'Separation_확보/최종_OLS')
    parser.add_argument('--optimization-starts',type=int,default=16)
    parser.add_argument('--optimization-max-sweeps',type=int,default=60)
    args = parser.parse_args()
    if args.self_test:
        workflow_tests(); return
    if args.overlap_only:
        out=args.output.resolve()
        selected=gpd.read_file(out/'ols_selected_surfaces.gpkg',layer='ols_selected_surfaces').to_dict('records')
        obstacles=gpd.read_file(args.buildings).to_crs(5186)
        draw_overlap_counts(selected,obstacles,out)
        return
    if args.png_only:
        redraw_existing_pngs(args.output.resolve(),args.buildings,args.dem)
        return
    if args.optimization_starts<1 or args.optimization_max_sweeps<1:
        parser.error('최적화 시작점 수와 반복 한도는 1 이상이어야 합니다.')
    import rasterio
    summary_path=args.candidates/'generation_summary.csv'
    from 장애물접근표면 import create_ols_triangle
    summary=pd.read_csv(summary_path).sort_values('id')
    check_inputs(summary,args.candidates)
    buildings=gpd.read_file(args.buildings).to_crs(5186)
    if buildings.geometry.isna().any() or buildings.is_empty.any() or not buildings.is_valid.all():
        raise ValueError('유효하지 않거나 비어 있는 장애물 형상을 정정해야 합니다.')
    if summary.host_index.isna().any() or not summary.host_index.between(0,len(buildings)-1).all():
        raise ValueError('후보 건물 원본 host_index 매핑을 확인해야 합니다.')
    footprints=np.array(shapely.force_2d(buildings.geometry.array),dtype=object)
    tree=shapely.STRtree(footprints)
    transform=Transformer.from_crs(4326,5186,always_xy=True)
    owners=[]
    for host, group in summary.groupby('host_index',sort=True):
        vp=next(group.itertuples())
        # 동일 물리 건물에 서로 다른 중심점이 있으면 중심점을 하나로 합치지 않는다.
        centers=group[['longitude','latitude','apex_z_m']].drop_duplicates()
        if len(centers)!=1:
            raise ValueError(f'host {host}: 서로 다른 중심/높이의 설계 그룹을 별도로 정의해야 합니다.')
        x,y=transform.transform(vp.longitude,vp.latitude)
        distance=Point(x,y).distance(footprints[int(host)])
        owners.append(dict(building_id=f'HOST_{int(host)}',host_source_row=int(host),
            source_vp_ids=json.dumps(group.id.astype(int).tolist()),center_x=x,center_y=y,
            z0_dem_agl_m=float(vp.apex_z_m),current_host_distance_m=distance,
            source_host_distance_m=float(vp.host_distance_m),
            owner_mapping_note='generation_summary host_index 소유관계 보존; 현 형상 거리 재검증 수치 참조',
            geometry=footprints[int(host)]))
    raw,source_files=[],[summary_path]
    for vp in summary.itertuples():
        path=args.candidates/str(int(vp.id))/f'Polygon_{int(vp.id)}.shp'
        if not path.exists():
            continue
        source_files.extend(p for p in path.parent.glob(path.stem+'.*') if p.is_file())
        frame=gpd.read_file(path).to_crs(5186)
        cx,cy=transform.transform(vp.longitude,vp.latitude)
        for source_row, feature in frame.iterrows():
            surface_plane(feature.geometry)
            raw.append(dict(building_id=f'HOST_{int(vp.host_index)}',source_vp_id=int(vp.id),
                polygon_id=f'vp{int(vp.id)}_row{source_row}',source_row=int(source_row),
                source_file=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                azimuth_deg=float(feature.rot_deg)%360,center_x=cx,center_y=cy,
                z0_dem_agl_m=float(feature.apex_z),slope_deg=float(feature.slope_deg),
                length_m=float(feature.height_m),start_width_m=0.,end_width_m=float(feature.width_m),
                surface_model_id='TRIANGLE_LINEAR_DEM_AGL',z_ref='DEM+AGL',generated=False,
                geometry=feature.geometry))
    unique,duplicate_map=deduplicate(raw)
    errors={}
    for owner in owners:
        if any(r['building_id']==owner['building_id'] for r in unique):
            continue
        vp=summary.loc[summary.host_index.eq(owner['host_source_row'])].iloc[0]
        try:
            generated=[]
            for angle in range(360):
                tri=create_ols_triangle(latitude=vp.latitude,longitude=vp.longitude,
                    apex_agl_m=vp.agl_m,ground_z_m=vp.ground_z_m,rotation_deg=angle,
                    horizontal_height_m=1000,base_width_m=573,slope_deg=8)
                xyz=np.asarray(tri.polygon.exterior.coords)
                x,y=transform.transform(xyz[:,0],xyz[:,1])
                generated.append(dict(building_id=owner['building_id'],source_vp_id=int(vp.id),
                    polygon_id=f"generated_{owner['building_id']}_{angle}",source_row=-1,source_file='GENERATED',
                    azimuth_deg=float(angle),center_x=owner['center_x'],center_y=owner['center_y'],
                    z0_dem_agl_m=float(vp.apex_z_m),slope_deg=8.,length_m=1000.,start_width_m=0.,
                    end_width_m=573.,surface_model_id='TRIANGLE_LINEAR_DEM_AGL',z_ref='DEM+AGL',
                    generated=True,geometry=Polygon(np.column_stack([x,y,xyz[:,2]]))))
            unique.extend(generated)
        except (ValueError,TypeError) as exc:
            errors[owner['building_id']]=str(exc)
    mark_near_duplicates(unique)
    # 모든 중복 제거를 마친 뒤 장애물 평가. 소유자가 달라도 동일한 3D 표면 계산은 캐시한다.
    agl=pd.to_numeric(buildings.A16,errors='coerce').to_numpy(float)
    with rasterio.open(args.dem) as dem:
        ground=sample_ground(footprints,5186,dem)
        start_ground=sample_ground([Point(r.longitude,r.latitude) for r in summary.itertuples()],4326,dem)
    if not np.allclose(start_ground,summary.ground_z_m,rtol=0,atol=1e-5):
        raise ValueError('기존 Polygon 기준 DEM과 현재 DEM 표본이 다릅니다.')
    roofs=np.where(np.isfinite(agl)&(agl>0)&np.isfinite(ground),ground+agl,np.nan)
    cache={}
    for i,record in enumerate(unique):
        key=geometry_key(record,include_owner=False)
        if key not in cache:
            cache[key]=evaluate(record['geometry'],footprints,roofs,tree)
        stats=cache[key]
        record.update(stats,status=stats['result'],penetration_count=stats['known_penetrations'],
            unknown_count=stats['unknown_buildings'],
            worst_penetration_m=max(0.,-stats['min_clearance_m']) if stats['min_clearance_m'] is not None and stats['known_penetrations'] else 0.)
        record['surface_id']=hashlib.sha256(repr(key[:-1]).encode()+key[-1]).hexdigest()[:20]
        if (i+1)%500==0:
            print(f'평가 {i+1}/{len(unique)}',flush=True)
    matrix=overlap_matrix(unique)
    bundles,baseline_selection=admissible_bundles(unique,owners)
    joint_selection,optimization_history,optimization_info=optimize_overlap(bundles,baseline_selection,matrix,
        starts=args.optimization_starts,max_sweeps=args.optimization_max_sweeps)
    baseline_records=[unique[i] for bundle in baseline_selection.values() for i in bundle]
    optimized_indices=[i for bundle in joint_selection.values() for i in bundle]
    for i,record in enumerate(unique):
        others=[j for j in optimized_indices if unique[j]['building_id']!=record['building_id']]
        record['overlap_with_other_selected_m2']=float(matrix[i,others].sum())
        record['overlap_other_selected_count']=int(np.count_nonzero(matrix[i,others]))
    selected,summaries,unresolved=[],[],[]
    for owner in owners:
        options=[r for r in unique if r['building_id']==owner['building_id']]
        chosen=sorted([unique[i] for i in joint_selection.get(owner['building_id'],())],key=single_key)
        pair_status='SELECTED' if len(chosen)==2 else 'NOT_FOUND'
        diff=angle_distance(chosen[0]['azimuth_deg'],chosen[1]['azimuth_deg']) if len(chosen)==2 else None
        reasons=[]
        if not options:
            reasons.append('NO_GEOMETRY: '+errors.get(owner['building_id'],'생성 가능한 후보 없음'))
        if options and all(r['status']=='FAIL' for r in options):
            reasons.append('ALL_CANDIDATES_FAIL')
        if any(r['status']=='REVIEW' for r in chosen):
            reasons.append('DATA_REVIEW')
        if any(r['unknown_count']>0 for r in chosen):
            reasons.append('UNKNOWN_OBSTACLE_HEIGHT')
        if any(r['touching_buildings']>0 for r in chosen):
            reasons.append('SURFACE_CONTACT')
        if pair_status=='NOT_FOUND':
            reasons.append('NO_ACCEPTABLE_PAIR_135_DEG')
        if abs(owner['current_host_distance_m']-owner['source_host_distance_m'])>0.01:
            reasons.append('OWNER_DISTANCE_CHANGED_CHECK_SOURCE_MAPPING')
        if sum(r['overlap_with_other_selected_m2'] for r in chosen)>AREA_TOL:
            reasons.append('REMAINING_INTER_BUILDING_XY_OVERLAP')
        for rank,item in enumerate(chosen,1):
            check=evaluate(item['geometry'],footprints,roofs,tree)
            brute=np.flatnonzero(shapely.intersects(footprints,shapely.force_2d(item['geometry'])))
            assert set(brute)==set(tree.query(shapely.force_2d(item['geometry']),predicate='intersects'))
            assert check['result']==item['status'] and check['min_clearance_m']==item['min_clearance_m']
            flag={'PASS':'PASS_SELECTED','REVIEW':'CONDITIONAL_REVIEW','FAIL':'UNSAFE_CANDIDATE'}[item['status']]
            selected.append(dict(item,selected_rank=rank,selection_flag=flag,
                pair_angle_diff_deg=diff,pair_status=pair_status,revalidated=True,
                source_vp_ids=owner['source_vp_ids']))
        counts=pd.Series([r['status'] for r in options],dtype=str).value_counts()
        summary_row=dict(building_id=owner['building_id'],source_vp_ids=owner['source_vp_ids'],
            host_source_row=owner['host_source_row'],raw_candidate_count=sum(r['building_id']==owner['building_id'] for r in raw),
            unique_candidate_count=len(options),generated_count=sum(r['generated'] for r in options),
            pass_count=int(counts.get('PASS',0)),review_count=int(counts.get('REVIEW',0)),fail_count=int(counts.get('FAIL',0)),
            selected_count=len(chosen),rank1_polygon_id=chosen[0]['polygon_id'] if chosen else None,
            rank1_azimuth_deg=chosen[0]['azimuth_deg'] if chosen else None,rank1_status=chosen[0]['status'] if chosen else 'NO_GEOMETRY',
            rank2_polygon_id=chosen[1]['polygon_id'] if len(chosen)==2 else None,
            rank2_azimuth_deg=chosen[1]['azimuth_deg'] if len(chosen)==2 else None,
            rank2_status=chosen[1]['status'] if len(chosen)==2 else None,pair_angle_diff_deg=diff,
            pair_status=pair_status,review_reasons='; '.join(reasons),
            current_host_distance_m=owner['current_host_distance_m'],source_host_distance_m=owner['source_host_distance_m'])
        old_indices=baseline_selection.get(owner['building_id'],())
        old_other=[j for name,bundle in baseline_selection.items() if name!=owner['building_id'] for j in bundle]
        summary_row.update(baseline_azimuths=json.dumps([unique[i]['azimuth_deg'] for i in old_indices]),
            baseline_overlap_with_other_buildings_m2=float(matrix[np.ix_(list(old_indices),old_other)].sum()),
            selected_overlap_with_other_buildings_m2=sum(r['overlap_with_other_selected_m2'] for r in chosen))
        summaries.append(summary_row)
        if reasons:
            unresolved.append(summary_row)
        print(f"{owner['building_id']}: {[(r['azimuth_deg'],r['status']) for r in chosen]} / 각도차 {diff}",flush=True)
    baseline_pairs,baseline_stats=overlap_audit(baseline_records)
    final_pairs,final_stats=overlap_audit(selected)
    assert np.isclose(final_stats['overlap_sum_m2'],optimization_info['optimized_overlap_sum_m2'],atol=AREA_TOL)
    assert np.isclose(baseline_stats['overlap_sum_m2'],optimization_info['baseline_overlap_sum_m2'],atol=AREA_TOL)
    assert final_stats['overlap_sum_m2']<=baseline_stats['overlap_sum_m2']+AREA_TOL
    reduction=baseline_stats['overlap_sum_m2']-final_stats['overlap_sum_m2']
    reduction_pct=100*reduction/baseline_stats['overlap_sum_m2'] if baseline_stats['overlap_sum_m2'] else 0.
    comparison=[dict(stage='independent_baseline',**baseline_stats),dict(stage='overlap_optimized',**final_stats)]
    print(f"중첩 면적: {baseline_stats['overlap_sum_m2']:,.3f} → {final_stats['overlap_sum_m2']:,.3f}m² / {reduction_pct:.2f}% 감소",flush=True)
    # 명시적 필수 검증: 소유자별 하나의 rank1, 최대 두 표면, 각도/중복/상태 보호.
    assert {r['building_id'] for r in summaries}=={r['building_id'] for r in owners}
    for owner in owners:
        chosen=[r for r in selected if r['building_id']==owner['building_id']]
        if not chosen:
            assert owner['building_id'] in errors
            continue
        assert sum(r['selected_rank']==1 for r in chosen)==1 and len(chosen) in (1,2)
        if len(chosen)==2:
            assert angle_distance(chosen[0]['azimuth_deg'],chosen[1]['azimuth_deg'])>=135
            assert all(r['status']!='FAIL' for r in chosen)
    assert len({geometry_key(r) for r in unique})==len(unique)
    assert len({geometry_key(r) for r in selected})==len(selected)
    assert len(duplicate_map)==len(raw)
    assert all(m['representative_polygon_id'] in {r['polygon_id'] for r in unique} for m in duplicate_map)
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    # 이전 버전 결과를 삭제하지 않고 하위 백업으로 이동하여 현재 결과와 혼동을 막는다.
    previous=[p for p in out.iterdir() if p.is_file() and p.name.startswith('ols_')]
    if previous:
        archive=out/'previous_runs'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        archive.mkdir(parents=True)
        for source in previous:
            if source.resolve().parent!=out:
                raise ValueError('백업 경로가 출력 폴더를 벗어났습니다.')
            shutil.move(str(source),str(archive/source.name))
    backup=out/'input_backup'
    backup.mkdir(exist_ok=True)
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    workflow=Path(__file__).with_name('Vertiport_Final_OLS_Selection_Workflow.md')
    with zipfile.ZipFile(backup/f'ols_candidates_raw_{stamp}.zip','w',compression=zipfile.ZIP_DEFLATED) as bundle:
        for path in sorted(set(source_files+[workflow,Path(__file__)])):
            bundle.write(path,str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else path.name)
    def flat(records):
        return pd.DataFrame([{k:v for k,v in r.items() if k!='geometry'} for r in records])
    flat(duplicate_map).to_csv(out/'ols_duplicate_map.csv',index=False,encoding='utf-8-sig')
    flat(unique).to_csv(out/'ols_candidate_scores.csv',index=False,encoding='utf-8-sig',float_format='%.9f')
    pd.DataFrame(summaries).to_csv(out/'ols_selection_summary.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(unresolved,columns=list(summaries[0])).to_csv(out/'ols_unresolved.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(comparison).to_csv(out/'ols_overlap_summary.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(optimization_history).to_csv(out/'ols_optimization_history.csv',index=False,encoding='utf-8-sig')
    pair_columns=['building_a','polygon_a','azimuth_a','building_b','polygon_b','azimuth_b',
        'overlap_area_m2','z_a_minus_b_min_m','z_a_minus_b_max_m','min_vertical_gap_m','surfaces_cross_or_touch']
    flat(final_pairs).reindex(columns=pair_columns).to_csv(out/'ols_surface_overlap_pairs.csv',index=False,encoding='utf-8-sig')
    unique_path=out/'ols_candidates_unique.gpkg'
    write_layer(unique,unique_path,'ols_candidates_unique','Polygon Z')
    write_layer(owners,unique_path,'candidate_buildings','Polygon')
    selected_path=out/'ols_selected_surfaces.gpkg'
    write_layer(selected,selected_path,'ols_selected_surfaces','Polygon Z')
    write_layer(baseline_records,selected_path,'ols_independent_baseline','Polygon Z')
    write_layer(final_pairs,selected_path,'ols_surface_overlaps','MultiPolygon')
    write_layer(baseline_pairs,selected_path,'ols_baseline_overlaps','MultiPolygon')
    write_layer([r for r in selected if r['selected_rank']==1],selected_path,'ols_selected_rank1','Polygon Z')
    write_layer([r for r in selected if r['selected_rank']==2],selected_path,'ols_selected_rank2','Polygon Z')
    # 서로 다른 소유자 간에도 완전히 같은 형상은 지도에서 한 번만 표시할 전용 레이어.
    map_features={}
    for record in selected:
        key=geometry_key(record,False)
        if key not in map_features:
            map_features[key]=dict(record,owner_building_ids=record['building_id'])
        else:
            map_features[key]['owner_building_ids']+=','+record['building_id']
    write_layer(list(map_features.values()),selected_path,'ols_selected_display_unique','Polygon Z')
    write_layer([{**r,'geometry':Point(o['center_x'],o['center_y'],o['z0_dem_agl_m'])}
                 for r,o in zip(summaries,owners)],selected_path,'candidate_building_decisions','Point Z')
    metadata=dict(workflow=str(workflow),building_input=str(args.buildings.resolve()),dem=str(args.dem.resolve()),
        owner_definition='generation_summary.csv의 기존 host_index로 8개 후보 건물 그룹; 원본 16 VP 관계 보존. 현재 footprint와 기록된 host 거리 차이를 별도 검토 표기',
        dedup='동일 소유자/중심/설계 및 normalize(Polygon Z) 완전일치. 반올림 없음. 1cm 이내 유사 형상은 near_duplicate 플래그만 표시. 서로 다른 소유자 관계 유지; 공유 형상 계산 캐시 및 지도 전용 중복제거',
        design='원본 길이 1000m/끝단573m/시작0m/경사8도 유지. 문서 8:1은 예시. 진북 기준 시계방향',
        height='DEM+AGL 동일 기준으로 비교. 입력 DEM의 AMSL 수직 datum은 독립 검증되지 않았음',
        obstacle_policy='A16 AGL; 결측/0/음수 및 DEM 결측은 UNKNOWN. A17 자동 대체/host 제외 없음',
        contact=f'침범 <-{TOL} FAIL; ±{TOL}m 접촉은 REVIEW; 나머지 결측 중첩 REVIEW; 모두 없으면 PASS',
        selection='건물별 최대 배정수 및 PASS/PASS > PASS/REVIEW > REVIEW/REVIEW 최상위 조합 고정, FAIL 포함 쌍 제외. 다음 건물 간 총 XY 중첩 면적 합 최소화. 면적 동률이면 기존 결측/접촉/여유고 순서. 단일 FAIL은 최소 침범량·침범수까지 먼저 고정',
        overlap_objective='서로 다른 건물의 선정 표면 쌍별 XY 교차면적 합(m²). 각 쌍 한 번 계산; 다중 중첩은 쌍별로 가산. 면적 없는 경계 접촉은 목적함수 0. union 면적은 별도 보조지표',
        overlap_optimizer=json.dumps(optimization_info,ensure_ascii=False),
        overlap_method='전체 허용 조합 보존, 고정 seed 16개 시작점(기본값)의 건물 단위 좌표 하강. 다른 건물 고정 시 해당 건물 모든 허용 쌍 탐색. 개선 불가능한 국소해 검증. 전역 최소 보장은 없음',
        overlap_height_audit='중첩 영역 모든 꼭짓점에서 두 선형 표면의 높이차 범위 및 최소 절대높이차 계산. 부호가 바뀌거나 0에 닿으면 3D 표면 교차/접촉. 높이차로 XY 면적을 할인하지 않음. 최소 운항 분리 기준 적용 없음',
        no_overlap='min_clearance NULL은 교차한 유효 건물 없음; PASS 후보 비교에서 제한 장애물 없음(+infinity)으로 취급',
        fallback='기존 후보가 없는 소유자만 360개 생성. 생성 실패는 NO_GEOMETRY. 기존 자료는 선행 필터된 후보이며 모든 360방향을 새로 검색한 결과가 아님',
        uncertainty='부분 교차 전체 꼭짓점 높이 최솟값 사용; 지형/운항/측면 이격은 선정 조건 아님. 건물 간 표면 중첩은 공동 목적함수에 포함',
        input_raw_count=str(len(raw)),owner_unique_count=str(len(unique)),evaluated_geometry_count=str(len(cache)),
        candidate_building_count=str(len(owners)),source_vp_count=str(len(summary)),
        created_utc=datetime.now(timezone.utc).isoformat())
    tables={'selection_summary':pd.DataFrame(summaries),'run_metadata':pd.DataFrame(metadata.items(),columns=['key','value']),
            'owner_vp_mapping':summary.assign(building_id=summary.host_index.map(lambda i:f'HOST_{int(i)}')),
            'duplicate_map':pd.DataFrame(duplicate_map),'overlap_summary':pd.DataFrame(comparison),
            'optimization_history':pd.DataFrame(optimization_history)}
    add_tables(unique_path,tables); add_tables(selected_path,tables)
    draw_maps(selected,owners,buildings,out,dem_path=args.dem)
    draw_overlap_comparison(baseline_records,selected,comparison,buildings,out,dem_path=args.dem)
    draw_overlap_counts(selected,buildings,out)
    notes=['# 후보 건물별 최종 OLS 선정',
        f'- 원본 {len(raw)}개 → 소유 건물별 고유 후보 {len(unique)}개. 후보 건물 {len(owners)}개, 원본 VP {len(summary)}개.',
        f'- 선정 표면 {len(selected)}개. 2개 선정 건물 {sum(r["selected_count"]==2 for r in summaries)}개.',
        f'- 건물 간 쌍별 중첩 면적 합: {baseline_stats["overlap_sum_m2"]:.3f} → {final_stats["overlap_sum_m2"]:.3f}m² ({reduction_pct:.2f}% 감소).',
        f'- 3D 표면 교차/접촉 쌍 수: {baseline_stats["crossing_pair_count"]} → {final_stats["crossing_pair_count"]}. XY 면적이 목적함수이며 3D 교차 수는 최소화하지 않았다.',
        '- 기존 독립 선정과 동일한 건물별 배정수 및 최상위 판정 조합을 유지. 중첩 우선으로 같은 판정 안의 결측 개수/여유고는 달라질 수 있음.',
        '- 건물별 원본 VP 관계, 원본 Polygon→대표 Polygon은 GPKG 표와 CSV에 보존.',
        '- 실선 PASS, 점선 REVIEW, 붉은색 FAIL. CONDITIONAL_REVIEW/UNSAFE_CANDIDATE는 승인된 안전 표면이 아님.',
        '- 기존 host_index에 따른 소유관계를 유지했으나 현재 footprint까지 거리와 기존 기록 거리가 다르므로 후보 건물 매핑은 추가 확인 대상.',
        '- 이전 실행 파일은 previous_runs, 원본 후보 백업은 input_backup에 보관.',
        '']
    notes.extend(f'- {k}: {v}' for k,v in metadata.items())
    (out/'ols_run_notes.md').write_text('\n'.join(notes),encoding='utf-8')
    loaded=gpd.read_file(selected_path,layer='ols_selected_surfaces')
    assert len(loaded)==len(selected) and loaded.geometry.has_z.all()
    assert set(loaded.building_id)|set(errors)=={o['building_id'] for o in owners}
    assert loaded.groupby('building_id').selected_rank.apply(lambda x:(x==1).sum()).eq(1).all()
    print(pd.DataFrame(summaries)[['building_id','source_vp_ids','selected_count','rank1_azimuth_deg','rank1_status','rank2_azimuth_deg','rank2_status','pair_angle_diff_deg']].to_string(index=False))
    print(f'완료: {out}')


if __name__ == '__main__':
    main()
