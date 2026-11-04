"""Inside-ring alternatives: hard XY no-crossing, or maximum Z gap at crossings.

Entry: python -B optimization_path.py --inside-new-options
No input or previous result is modified. No new altitude/separation threshold
is assumed: all previously recorded geometry and Landing_Area rules remain.
Vertical alternative retains ALL feasible opposite-FAF STAR/SID pairs and
maximizes the worst absolute Z difference over their actual XY intersections.
An XY intersection is not interpolated by Shapely's blended intersection Z:
each original 3D route is evaluated independently at the same XY location.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import math
from pathlib import Path
import time

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point

import optimization_path as core


def make_route(xyz):
    xyz = np.asarray(xyz,dtype=float)
    lengths = np.r_[0.,np.cumsum(np.linalg.norm(np.diff(xyz[:,:2],axis=0),axis=1))]
    line = LineString(xyz[:,:2])
    if not line.is_simple:
        raise ValueError('Crossing-height evaluation requires simple individual routes.')
    return dict(xyz=xyz,s=lengths,line=line)


def heights_at(route,xy):
    from shapely import points, line_locate_point
    distances = line_locate_point(route['line'],points(np.asarray(xy)))
    return np.interp(distances,route['s'],route['xyz'][:,2])


def crossing_height_profile(a,b,details=False):
    """Exact piecewise-linear heights, including overlapping line intervals.

    Break an overlap at BOTH routes' vertices. Difference in Z is linear
    between breaks, so a sign change yields a true zero-clearance point.
    Return inf for disjoint lines (no crossing to score).
    """
    hit = a['line'].intersection(b['line'])
    if hit.is_empty:
        return math.inf,[]
    minimum,records = math.inf,[]

    def evaluate(xy,kind,continuous=False):
        nonlocal minimum
        xy = np.asarray(xy,dtype=float)
        za,zb = heights_at(a,xy),heights_at(b,xy)
        signed = za-zb
        minimum = min(minimum,float(np.abs(signed).min()))
        if details:
            for point,z1,z2 in zip(xy,za,zb):
                records.append(dict(kind=kind,x=float(point[0]),y=float(point[1]),
                                    z_a=float(z1),z_b=float(z2),signed_dz=float(z1-z2),
                                    abs_dz=float(abs(z1-z2))))
        if continuous:
            for i in range(len(xy)-1):
                if signed[i]*signed[i+1]<0:
                    minimum = 0.
                    if details:
                        t = signed[i]/(signed[i]-signed[i+1])
                        point = xy[i]+t*(xy[i+1]-xy[i])
                        z1 = za[i]+t*(za[i+1]-za[i])
                        z2 = zb[i]+t*(zb[i+1]-zb[i])
                        records.append(dict(kind='overlap_zero',x=float(point[0]),y=float(point[1]),
                                            z_a=float(z1),z_b=float(z2),signed_dz=0.,abs_dz=0.))

    def visit(component):
        if component.geom_type=='Point':
            evaluate([component.coords[0]],'point')
        elif component.geom_type in ('LineString','LinearRing'):
            from shapely import points, distance
            coords = list(component.coords)
            for route in (a,b):
                xy = route['xyz'][:,:2]
                mask = distance(points(xy),component)<=1e-7
                coords.extend(xy[mask])
            coords = np.unique(np.asarray(coords),axis=0)
            order = np.argsort([component.project(Point(p)) for p in coords])
            evaluate(coords[order],'overlap_break',continuous=True)
        elif hasattr(component,'geoms'):
            for child in component.geoms:
                visit(child)
        else:
            raise ValueError(f'Unexpected route intersection: {component.geom_type}')
    visit(hit)
    return minimum,records


def build_vertical_matrix(pool,routes,crossings,output,key):
    """Compute only route pairs that could be co-selected, caching exact gaps."""
    cache = output/'crossing_z_gaps.npz'
    if cache.exists():
        with np.load(cache,allow_pickle=False) as saved:
            if saved['key'].item()==key:
                print('Reusing validated crossing-height matrix.',flush=True)
                return saved['gap'].copy()
    n = len(pool)
    result = np.full((n,n),np.inf)
    ii,jj = np.where(np.triu(crossings,1))
    pair,op,faf = pool.pair.to_numpy(),pool.operation.to_numpy(),pool.faf_id.to_numpy()
    valid = ~((pair[ii]==pair[jj])&((op[ii]==op[jj])|(faf[ii]==faf[jj])))
    edges = np.column_stack((ii[valid],jj[valid]))
    same_node = pool.iaf_id.to_numpy()
    started,completed,next_update = time.monotonic(),0,10000

    def block(first):
        values = []
        for i,j in edges[first:first+256]:
            # A common ring node has the same XYZ (350 m) in both routes.
            gap = 0. if same_node[i]==same_node[j] else crossing_height_profile(routes[i],routes[j])[0]
            values.append((int(i),int(j),gap))
        return values
    with ThreadPoolExecutor(max_workers=4) as executor:
        tasks = [executor.submit(block,i) for i in range(0,len(edges),256)]
        for task in as_completed(tasks):
            values = task.result()
            for i,j,gap in values:
                result[i,j] = result[j,i] = gap
            completed += len(values)
            if completed>=next_update:
                print(f'Crossing-height pairs {completed}/{len(edges)}; {time.monotonic()-started:.1f}s',flush=True)
                next_update += 10000
    np.savez_compressed(cache,key=key,gap=result)
    return result


def network_z_minimum(ids,crossings,z_gaps):
    block = np.ix_(ids,ids)
    mask = np.triu(crossings[block],1)
    values = z_gaps[block][mask]
    finite = values[np.isfinite(values)]
    return float(finite.min()) if len(finite) else math.inf


def maximize_crossing_height(pool,options,crossings,z_gaps,seed,config,fafs,seconds,tolerance):
    """Fixed full coverage, max-min Z gap. Never treat timeout as infeasibility."""
    count = sum(bool(v) for v in options.values())
    if len(seed)!=2*count:
        raise ValueError('Vertical alternative needs a full-coverage incumbent.')
    low = network_z_minimum(seed,crossings,z_gaps)
    if math.isinf(low):
        return seed,dict(status='full_coverage_has_no_planar_intersections',floor_m=None,probes=[])
    high = float(config.ring_z-min(f['xyz'][2] for f in fafs))
    high = max(high,low)
    best,probes = list(seed),[]
    started = time.monotonic()
    unknown = False
    while high-low>tolerance and time.monotonic()-started<seconds-1:
        target = (low+high)/2
        forbidden = crossings&(z_gaps<target-1e-7)
        remaining = seconds-(time.monotonic()-started)
        ids,status = core.traffic_noncrossing_subset(pool,options,forbidden,
            seconds=min(6.,remaining),required_pairs=count)
        probes.append(dict(target_z_m=target,**status))
        if len(ids)==2*count:
            actual = network_z_minimum(ids,crossings,z_gaps)
            if math.isinf(actual):
                return ids,dict(status='full_coverage_has_no_planar_intersections',floor_m=None,probes=probes)
            if actual<target-1e-6:
                raise RuntimeError('Vertical-gap MILP incumbent fails its gap constraint.')
            best,low = ids,max(low,actual)
        elif status.get('status_code')==2:
            high = target
        else:
            unknown = True
            break
        print(f'Vertical worst-gap bound: {low:.3f} .. {high:.3f} m',flush=True)
    return best,dict(status='fixed_coverage_maxmin_vertical_gap',floor_m=low,upper_bound_m=high,
                     tolerance_m=tolerance,optimal_to_tolerance_in_shortlist=bool(not unknown and high-low<=tolerance),
                     seconds=time.monotonic()-started,probes=probes)


def crossing_audit(result,output):
    """Report all exact intersections, each route Z, and min Z per route pair."""
    routes = [make_route(np.asarray(g.coords)) for g in result.geometry]
    detail,pairs,levels = [],[],[]
    points = []
    for i,a in enumerate(routes):
        for j in range(i+1,len(routes)):
            gap,items = crossing_height_profile(a,routes[j],details=True)
            if not items:
                continue
            ra,rb = result.iloc[i],result.iloc[j]
            common = dict(route_a=int(ra.cand_id),operation_a=ra.operation,pair_a=int(ra.pair),
                          route_b=int(rb.cand_id),operation_b=rb.operation,pair_b=int(rb.pair))
            pairs.append(dict(**common,min_abs_dz=gap,evaluation_locations=len(items)))
            for item in items:
                record = dict(**common,**item)
                detail.append(record)
                points.append(Point(item['x'],item['y'],(item['z_a']+item['z_b'])/2))
                for role,z in (('A',item['z_a']),('B',item['z_b'])):
                    levels.append(dict(route_id=common['route_a'] if role=='A' else common['route_b'],
                                       role=role,alt=z,geometry=Point(item['x'],item['y'],z)))
    pair_cols = ['route_a','operation_a','pair_a','route_b','operation_b','pair_b','min_abs_dz','evaluation_locations']
    detail_cols = pair_cols[:6]+['kind','x','y','z_a','z_b','signed_dz','abs_dz']
    pd.DataFrame(pairs,columns=pair_cols).to_csv(output/'crossing_pair_summary.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(detail,columns=detail_cols).to_csv(output/'crossing_details.csv',index=False,encoding='utf-8-sig')
    if detail:
        # Short field names are deliberate: Shapefile DBF allows 10 chars.
        attributes = [dict(route_a=r['route_a'],route_b=r['route_b'],pair_a=r['pair_a'],pair_b=r['pair_b'],
                           z_a=r['z_a'],z_b=r['z_b'],abs_dz=r['abs_dz'],kind=r['kind']) for r in detail]
        gpd.GeoDataFrame(attributes,geometry=points,crs=5186).to_file(
            output/'crossing_points_z.shp',engine='pyogrio',encoding='UTF-8',geometry_type='Point Z')
        gpd.GeoDataFrame(levels,crs=5186).to_file(output/'crossing_levels_z.shp',
            engine='pyogrio',encoding='UTF-8',geometry_type='Point Z')
        fig,ax = core.plt.subplots(figsize=(11,5))
        labels = [f"{p['operation_a']} P{p['pair_a']} / {p['operation_b']} P{p['pair_b']}" for p in pairs]
        ax.bar(np.arange(len(pairs)),[p['min_abs_dz'] for p in pairs],color='steelblue')
        ax.set(xticks=np.arange(len(pairs)),xticklabels=labels,ylabel='Minimum absolute Z difference (m)',
               title='Worst vertical spacing at each XY-intersecting route pair')
        ax.tick_params(axis='x',rotation=70,labelsize=7)
        ax.grid(axis='y',alpha=.2)
        fig.tight_layout()
        fig.savefig(output/'crossing_z_spacing.png',dpi=160)
        core.plt.close(fig)
        fig = core.plt.figure(figsize=(12,9))
        ax = fig.add_subplot(111,projection='3d')
        for _,row in result.iterrows():
            xyz = np.asarray(row.geometry.coords)
            ax.plot(*xyz.T,ls='-' if row.operation=='STAR' else '--',alpha=.7)
        for item in detail:
            ax.plot([item['x']]*2,[item['y']]*2,[item['z_a'],item['z_b']],color='red',lw=1.)
        ax.set(xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)',title='Vertical links at plan-view crossings')
        fig.savefig(output/'crossing_vertical_links.png',dpi=160)
        core.plt.close(fig)
    return dict(planar_crossing_pairs=len(pairs),evaluated_locations=len(detail),
                minimum_crossing_z_m=min((r['min_abs_dz'] for r in pairs),default=None))


def export_case(pool,ids,folder,context,scenario='inside_only'):
    nodes,headings,fafs,areas,audit,config,traffic = context
    result,spacing = core.traffic_export(pool,ids,folder,nodes,headings,fafs,areas,audit,config,traffic,scenario)
    height = crossing_audit(result,folder)
    info = dict(selected_pairs=int(result.pair.nunique()),approach_count=int((result.operation=='STAR').sum()),
                departure_count=int((result.operation=='SID').sum()),length_m=float(result.length_m.sum()),
                landing_area_crossings=int(result.cross_n.sum()),
                minimum_xy_m=min((r['min_xy_m'] for r in spacing),default=None),**height)
    (folder/'summary.json').write_text(json.dumps(info,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    return result,info


def write_report(output,summary):
    sections = []
    labels = dict(no_xy_crossing='우선안: XY 교차 금지',strict_no_xy_crossing='Landing_Area도 관통 금지',
                  vertical_clearance_alternative='대안: 전체 연결 + 교차점 수직 간격 최대화')
    for name,info in summary['results'].items():
        gap = info.get('minimum_crossing_z_m')
        gap_text = '교차점 없음' if gap is None else f'{gap:.3f}m'
        sections.append(f'''<h2>{labels[name]}</h2><p>FAF {info['selected_pairs']}쌍,
Approach {info['approach_count']}개 / Departure {info['departure_count']}개.
XY 교차 {info['planar_crossing_pairs']}쌍, 교차점 최소 수직 간격: {gap_text}.
Landing_Area 관통 {info['landing_area_crossings']}건.</p>
<p><a href="{name}/Approach_z.shp">Approach SHP</a> · <a href="{name}/Departure_z.shp">Departure SHP</a> ·
<a href="{name}/IAF_points_z.shp">IAF 점</a> · <a href="{name}/FAF_points_z.shp">FAF 점</a> ·
<a href="{name}/crossing_details.csv">교점별 두 항로 높이 CSV</a></p>
<img src="{name}/map_combined.png"><details><summary>Approach / Departure 개별 지도</summary>
<img src="{name}/map_star.png"><img src="{name}/map_sid.png"></details>''')
        if gap is not None:
            sections.append(f'<img src="{name}/crossing_z_spacing.png"><img src="{name}/crossing_vertical_links.png">')
    html = f'''<!doctype html><html lang="ko"><meta charset="utf-8"><title>Inner new options</title>
<style>body{{font-family:Malgun Gothic,sans-serif;max-width:1250px;margin:32px auto;line-height:1.7}}img{{max-width:100%}}.note{{padding:14px;background:#fff1d5}}</style>
<h1>링 내부: 평면 무교차 우선 / 교차점 수직 이격 대안</h1>
<p>사용자 확인에 따라 두 조건을 다른 시나리오로 구분했습니다. 교차 금지를 충족하는 안에서는
교차점의 Z 간격이 정의되지 않습니다. 교차를 허용하는 대안에서는 전체 연결 범위를 유지한 채
모든 교차점의 최소 |ZA−ZB|를 최대화합니다.</p>
<p class="note">수직 간격이 크다는 이유만으로 실제 운항 안전성이 입증되는 것은 아닙니다.
이전의 70m를 하드 기준으로 재적용하지 않았으며, 정해진 수직 분리값을 충족한다고 주장하지 않습니다.
시각·풍향·기체 운항 성능은 포함하지 않은 기하 분석입니다.</p>
<p>기존 조건 유지: 링 350m·내부 경계, FAF alt, 70↔60kt, CW RSR, 구배와 뱅크 제한,
총 선회 180° 형태 필터·자기교차 금지·반대편 FAF의 STAR/SID 배정·기존 Landing_Area 처리.
수직 간격을 늘리기 위해 원본 고도나 경로의 Z 보간을 임의로 변경하지 않았습니다.</p>
{''.join(sections)}
<h2>교점 높이 계산</h2><p>두 항로의 XY 교점을 구한 뒤, 각 항로의 원본 XYZ에서 독립적으로 높이를 보간합니다.
교점이 여러 개면 가장 작은 간격을 사용합니다. XY 중첩 구간은 두 항로의 모든 꺾임점을 포함해
높이 차가 0이 되는 위치까지 검사합니다. Shapely가 교점 기하에 넣는 혼합 Z는 사용하지 않습니다.</p>
<p>crossing_points_z.shp의 기하 Z는 표시용 중간 높이이며 z_a/z_b가 실제 항로 높이입니다.
crossing_levels_z.shp는 같은 XY에 각 항로의 실제 높이 점 두 개를 저장합니다.</p>
<p>최적성은 저장된 선별 후보 안의 의미입니다. 전체 연속 경로나 모든 원자료 후보 조합의 전역 최적성은
주장하지 않습니다. 시간 제한은 불가능성 증명으로 해석하지 않습니다.</p>
<h2>재실행</h2><pre>python -B 수정본_20261103/코드/optimization_path.py --inside-new-options</pre></html>'''
    (output/'report.html').write_text(html,encoding='utf-8')


def self_check():
    a = make_route([[0,0,350],[10,0,350]])
    b = make_route([[5,-5,300],[5,5,300]])
    gap,items = crossing_height_profile(a,b,True)
    assert gap==50 and len(items)==1 and items[0]['z_a']==350 and items[0]['z_b']==300
    a = make_route([[0,0,0],[10,0,10]])
    b = make_route([[0,0,10],[10,0,0]])
    gap,items = crossing_height_profile(a,b,True)
    assert gap==0 and any(p['kind']=='overlap_zero' for p in items)
    b = make_route([[0,0,70],[10,0,80]])
    assert crossing_height_profile(a,b)[0]==70
    b = make_route([[0,1,0],[10,1,10]])
    assert math.isinf(crossing_height_profile(a,b)[0])
    a = make_route([[0,0,0],[10,0,10]])
    b = make_route([[5,-5,5],[5,5,5]])
    assert crossing_height_profile(a,b)[0]==0
    print('Crossing heights passed: independent Z interpolation, overlap zero, constant offset, disjoint, equal height.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=core.BASE/'결과물/optimization_path_reviewed')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--qgis',type=Path,default=core.BASE/'QGIS')
    parser.add_argument('--search-seconds',type=float,default=45.)
    parser.add_argument('--z-tolerance',type=float,default=.5)
    parser.add_argument('--mode',choices=['both','planar-only','vertical-only'],default='both')
    parser.add_argument('--self-check',action='store_true')
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return
    if not all(math.isfinite(v) and v>0 for v in (args.search_seconds,args.z_tolerance)):
        raise ValueError('Search time and tolerance must be finite and positive.')
    source = args.source.resolve()
    output = (args.output or source/'inside_only/new_options').resolve()
    output.mkdir(parents=True,exist_ok=True)
    previous = json.loads((source/'summary.json').read_text(encoding='utf-8'))
    config,traffic = core.Config(**previous['config']),core.TrafficConfig(**previous['traffic_config'])
    nodes,headings,fafs,areas,audit,inputs = core.prepare(args.qgis.resolve(),config)
    hashes = core.fingerprint(inputs)
    if hashes!=previous['input_sha256']:
        raise ValueError('Source input files changed; regenerate the reviewed candidate inventory first.')
    frame = pd.read_csv(source/'review_candidate_metrics.csv')
    pool,options,policy = core.traffic_pool(frame,'inside_only',config,traffic,audit)
    key = hashlib.sha256(json.dumps(dict(version=1,config=core.asdict(config),inputs=hashes,
        candidate_ids=pool.cand_id.astype(int).tolist()),sort_keys=True).encode('utf-8')).hexdigest()
    _,distance,crossings,gaps = core.traffic_pair_matrices(pool,nodes,headings,fafs,config,traffic,
        source/'spacing_matrix_inside_only.npz',key)
    context = nodes,headings,fafs,areas,audit,config,traffic
    summary = dict(config=core.asdict(config),traffic_config=core.asdict(traffic),policy=policy,
                   input_sha256=hashes,interpretation='hard XY no-crossing first; crossing-Z alternative separately',
                   hard_vertical_threshold_m=None,results={})
    if args.mode!='vertical-only':
        ids,cardinality = core.traffic_noncrossing_subset(pool,options,crossings,args.search_seconds)
        if ids:
            ids,spread = core.traffic_maxmin_spacing(pool,options,crossings,distance,fafs,ids,seconds=args.search_seconds)
            active = set(pool.loc[ids,'pair'].astype(int))
            ids,_ = core.traffic_select(pool,{p:v for p,v in options.items() if p in active},distance,crossings,gaps,
                nodes,traffic,config,initial_indices=ids,forbid_crossings=True,
                minimum_spacing_m=max(0.,spread.get('achieved_minimum_m',0.)-1e-6))
        else:
            spread = dict(status='no_eligible_pairs')
        _,stats = export_case(pool,ids,output/'no_xy_crossing',context)
        if stats['planar_crossing_pairs']!=0:
            raise RuntimeError('Primary no-crossing scenario has XY intersections.')
        summary['results']['no_xy_crossing'] = dict(**stats,coverage_search=cardinality,spacing_search=spread)
        strict = core.traffic_strict_layout(pool,options,policy['strict_pairs'],distance,crossings,gaps,
            nodes,headings,fafs,areas,audit,config,traffic,output/'strict_no_xy_crossing','inside_only',seconds=args.search_seconds)
        result = gpd.read_file(output/'strict_no_xy_crossing/combined_routes_z.shp')
        height = crossing_audit(result,output/'strict_no_xy_crossing')
        summary['results']['strict_no_xy_crossing'] = dict(selected_pairs=strict['selected_pairs'],
            approach_count=int((result.operation=='STAR').sum()),departure_count=int((result.operation=='SID').sum()),
            landing_area_crossings=strict['landing_crossings'],minimum_xy_m=strict['exact_min_distance_m'],
            length_m=strict['length_m'],**height)
    if args.mode!='planar-only':
        routes = [make_route(core.traffic_xyz(row,nodes,headings,fafs,config)[0]) for _,row in pool.iterrows()]
        z_key = hashlib.sha256((key+':exact-crossing-height-v1').encode()).hexdigest()
        z_gaps = build_vertical_matrix(pool,routes,crossings,output,z_key)
        baseline_ids,_ = core.traffic_select(pool,options,distance,crossings,gaps,nodes,traffic,config)
        baseline_gap = network_z_minimum(baseline_ids,crossings,z_gaps)
        ids,search = maximize_crossing_height(pool,options,crossings,z_gaps,baseline_ids,config,fafs,
                                             args.search_seconds,args.z_tolerance)
        floor = search.get('floor_m')
        forbidden = crossings if floor is None else crossings&(z_gaps<floor-1e-6)
        ids,_ = core.traffic_select(pool,options,distance,crossings,gaps,nodes,traffic,config,
            initial_indices=ids,forbid_crossings=True,forbidden_pairs=forbidden)
        _,stats = export_case(pool,ids,output/'vertical_clearance_alternative',context)
        if floor is not None and stats['minimum_crossing_z_m'] is not None and stats['minimum_crossing_z_m']<floor-1e-5:
            raise RuntimeError('Exported vertical alternative violates optimized crossing-height floor.')
        summary['results']['vertical_clearance_alternative'] = dict(**stats,vertical_search=search,
            baseline_minimum_z_m=baseline_gap if math.isfinite(baseline_gap) else None)
    if core.fingerprint(inputs)!=hashes:
        raise RuntimeError('Original input files changed during analysis.')
    summary['original_inputs_unchanged'] = True
    (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    write_report(output,summary)
    print(json.dumps(summary['results'],ensure_ascii=False,indent=2),flush=True)
    print(f'Complete: {output}',flush=True)


if __name__=='__main__':
    main()
