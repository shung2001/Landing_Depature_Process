"""교차 위치에서 두 원본 선의 Z를 독립적으로 보간해 수직 간격을 측정한다.

기본 실행: python -B 수정본_20261103/코드/Vertical_Separation_gauge.py
검증:      위 명령에 --self-check

분석: Short_Final 내부 / Approach·Departure 내부 / 두 자료 사이.
항로 폴더에서는 Approach_z.shp와 Departure_z.shp만 읽는다. 이 둘과 내용이
중복되는 combined_routes_z.shp 및 Point SHP는 다시 읽지 않는다.
XY 좌표는 EPSG:5186으로 통일하고 Z는 변환하지 않는다. 두 입력 Z가 같은
수직 기준·미터 단위를 사용한다고 가정한다. 결과는 저장된 Polyline Z의
선분 기준이며, 바람·시각·기체 크기 또는 안전 분리 기준을 판정하지 않는다.

중요: Shapely가 2D 교점에 부여하는 혼합 Z는 사용하지 않는다. 교차하는
각 원본 XYZ 선분에서 따로 Z를 구한다. 중첩 선분에서는 양 끝과 Z 차의
부호가 바뀌는 위치를 평가한다. 공유 끝점도 삭제하지 않고 별도로 표시한다.
필요 패키지: geopandas, shapely >= 2, numpy, pandas, matplotlib, pyogrio.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import LineString, Point

BASE = Path(__file__).resolve().parents[1]
PROJECT = BASE.parent
DEFAULT_FINAL = PROJECT/'자료/최종정리/Short_Final.shp'
DEFAULT_ROUTES = BASE/'결과물/optimization_path_reviewed/inside_only/new_options/vertical_clearance_alternative'
DEFAULT_OUTPUT = BASE/'결과물/Separation'

EVENT_COLUMNS = [
    'cross_id','group','route_a','label_a','file_a','fid_a','part_a','source_id_a','operation_a',
    'route_b','label_b','file_b','fid_b','part_b','source_id_b','operation_b',
    'segment_a','segment_b','x_5186','y_5186','z_a_m','z_b_m','signed_dz_m','abs_dz_m',
    'xy_relation','event_type','endpoint_a','endpoint_b','shared_endpoint_only',
]


@dataclass
class Route:
    key: str
    label: str
    source: str
    fid: int
    part: int
    source_id: str
    operation: str
    xyz: np.ndarray
    segments: object = field(init=False,repr=False)
    tree: object = field(init=False,repr=False)
    line: object = field(init=False,repr=False)

    def __post_init__(self):
        self.xyz = np.asarray(self.xyz,dtype=float)
        if self.xyz.ndim!=2 or self.xyz.shape[1]<3 or len(self.xyz)<2:
            raise ValueError(f'{self.key}: finite XYZ vertices are required.')
        self.xyz = self.xyz[:,:3]
        if not np.isfinite(self.xyz).all():
            raise ValueError(f'{self.key}: nonfinite XYZ coordinate.')
        # Exact duplicate vertices are redundant; vertical XY-zero segments
        # have a Z interval instead of one height and must not be guessed.
        keep = np.r_[True,np.any(np.diff(self.xyz,axis=0)!=0,axis=1)]
        self.xyz = self.xyz[keep]
        if len(self.xyz)<2:
            raise ValueError(f'{self.key}: fewer than two distinct vertices.')
        if (np.linalg.norm(np.diff(self.xyz[:,:2],axis=0),axis=1)==0).any():
            raise ValueError(f'{self.key}: XY-zero vertical segment requires a separate interval model.')
        self.line = LineString(self.xyz[:,:2])
        self.segments = shapely.linestrings(np.stack((self.xyz[:-1,:2],self.xyz[1:,:2]),axis=1))
        self.tree = shapely.STRtree(self.segments)


def attribute(value,default=''):
    if value is None or pd.isna(value):
        return default
    if isinstance(value,(int,np.integer)):
        return str(int(value))
    if isinstance(value,(float,np.floating)) and float(value).is_integer():
        return str(int(value))
    return str(value)


def read_routes(path,tag,operation):
    path = Path(path).resolve()
    data = gpd.read_file(path,engine='pyogrio',fid_as_index=True)
    if data.crs is None:
        raise ValueError(f'Missing horizontal CRS: {path}')
    data = data.to_crs(5186)
    routes = []
    for feature_id,row in data.iterrows():
        fid = int(feature_id)
        geom = row.geometry
        if geom is None or geom.is_empty:
            raise ValueError(f'{path}, FID={fid}: empty geometry.')
        if geom.geom_type=='LineString':
            parts = [geom]
        elif geom.geom_type=='MultiLineString':
            parts = list(geom.geoms)
        else:
            raise ValueError(f'{path}, FID={fid}: expected LineString Z, got {geom.geom_type}.')
        identity = attribute(row.get('cand_id',row.get('id')),str(fid))
        op = attribute(row.get('operation'),operation)
        pair,faf = attribute(row.get('pair')),attribute(row.get('faf_id'))
        label = f'{operation} P{pair} F{faf}' if pair and faf else f'{operation} id={identity}'
        for part,geometry in enumerate(parts):
            if not geometry.has_z:
                raise ValueError(f'{path}, FID={fid}, part={part}: Z geometry is required.')
            routes.append(Route(key=f'{tag}:FID{fid}:part{part}',label=label,
                                source=str(path),fid=fid,part=part,source_id=identity,
                                operation=op,xyz=np.asarray(geometry.coords)))
    return routes


def route_input_paths(folder):
    """Choose one nonduplicated representation of the route network."""
    folder = Path(folder).resolve()
    for approach,departure in [('Approach_z.shp','Departure_z.shp'),('star_routes_z.shp','sid_routes_z.shp')]:
        if (folder/approach).exists() and (folder/departure).exists():
            return [(folder/approach,'AP','Approach'),(folder/departure,'DP','Departure')]
    combined = folder/'combined_routes_z.shp'
    if combined.exists():
        return [(combined,'RT','Route')]
    raise FileNotFoundError(f'No separate Approach/Departure or combined route lines found: {folder}')


def file_hashes(paths):
    result = {}
    for path in sorted(set(Path(p).resolve() for p in paths)):
        for suffix in ('.shp','.shx','.dbf','.prj','.cpg'):
            source = path.with_suffix(suffix)
            if source.exists():
                result[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
    return result


def segment_height(xyz,index,xy):
    a,b = xyz[index],xyz[index+1]
    delta = b[:2]-a[:2]
    t = float(np.clip(np.dot(np.asarray(xy)-a[:2],delta)/np.dot(delta,delta),0.,1.))
    return float(a[2]+t*(b[2]-a[2]))


def at_endpoint(route,xy,tolerance):
    return bool(min(np.linalg.norm(np.asarray(xy)-route.xyz[0,:2]),
                    np.linalg.norm(np.asarray(xy)-route.xyz[-1,:2]))<=tolerance)


def pair_intersections(a,b,tolerance=1e-6):
    """Analyze original segment pairs; preserve distinct Z at self-crossings."""
    if not a.line.intersects(b.line):
        return []
    indices = b.tree.query(a.segments,predicate='intersects')
    events = {}
    priority = dict(point=0,overlap_endpoint=1,overlap_zero=2)

    def add(xy,ia,ib,kind):
        x,y = map(float,xy)
        za,zb = segment_height(a.xyz,ia,(x,y)),segment_height(b.xyz,ib,(x,y))
        ea,eb = at_endpoint(a,(x,y),tolerance),at_endpoint(b,(x,y),tolerance)
        key = tuple(int(round(v/tolerance)) for v in (x,y,za,zb))
        record = dict(segment_a=int(ia),segment_b=int(ib),x_5186=x,y_5186=y,
                      z_a_m=za,z_b_m=zb,signed_dz_m=za-zb,abs_dz_m=abs(za-zb),
                      xy_relation='point' if kind=='point' else 'overlap',event_type=kind,
                      endpoint_a=ea,endpoint_b=eb,shared_endpoint_only=bool(kind=='point' and ea and eb))
        old = events.get(key)
        if old is None or priority[kind]>priority[old['event_type']]:
            events[key] = record

    def visit(geometry,ia,ib):
        if geometry.is_empty:
            return
        if geometry.geom_type=='Point':
            add(geometry.coords[0],ia,ib,'point')
        elif geometry.geom_type=='LineString':
            xy = np.asarray(geometry.coords)
            for point in (xy[0],xy[-1]):
                add(point,ia,ib,'overlap_endpoint')
            dz0 = segment_height(a.xyz,ia,xy[0])-segment_height(b.xyz,ib,xy[0])
            dz1 = segment_height(a.xyz,ia,xy[-1])-segment_height(b.xyz,ib,xy[-1])
            if dz0*dz1<0:
                t = dz0/(dz0-dz1)
                add(xy[0]+t*(xy[-1]-xy[0]),ia,ib,'overlap_zero')
        elif hasattr(geometry,'geoms'):
            for child in geometry.geoms:
                visit(child,ia,ib)
        else:
            raise ValueError(f'Unexpected line intersection: {geometry.geom_type}')
    for ia,ib in sorted(map(tuple,indices.T)):
        visit(a.segments[ia].intersection(b.segments[ib]),int(ia),int(ib))
    return sorted(events.values(),key=lambda r:(r['x_5186'],r['y_5186'],r['z_a_m'],r['z_b_m']))


def analyze_group(a_routes,b_routes,name,prefix,tolerance,internal=False):
    records = []
    for i,a in enumerate(a_routes):
        for j,b in enumerate(b_routes):
            if internal and j<=i:
                continue
            # Parts of one original feature are not independent route pairs.
            if a.source==b.source and a.fid==b.fid:
                continue
            events = pair_intersections(a,b,tolerance)
            common = dict(group=name,route_a=a.key,label_a=a.label,file_a=a.source,fid_a=a.fid,
                          part_a=a.part,source_id_a=a.source_id,operation_a=a.operation,
                          route_b=b.key,label_b=b.label,file_b=b.source,fid_b=b.fid,
                          part_b=b.part,source_id_b=b.source_id,operation_b=b.operation)
            for event in events:
                records.append(dict(cross_id=f'{prefix}{len(records)+1:04d}',**common,**event))
    result = pd.DataFrame(records,columns=EVENT_COLUMNS)
    for column in ('endpoint_a','endpoint_b','shared_endpoint_only'):
        result[column] = result[column].astype(bool)
    return result


def summary_for(data,name,a_count,b_count):
    meaningful = data[~data.shared_endpoint_only] if len(data) else data
    pairs = data[['route_a','route_b']].drop_duplicates() if len(data) else data
    return dict(group=name,a_route_parts=a_count,b_route_parts=b_count,
                intersection_evaluations=len(data),intersecting_route_pairs=len(pairs),
                shared_endpoint_evaluations=int(data.shared_endpoint_only.sum()) if len(data) else 0,
                non_shared_endpoint_evaluations=len(meaningful),
                overlap_evaluations=int(data.xy_relation.eq('overlap').sum()) if len(data) else 0,
                minimum_all_m=float(data.abs_dz_m.min()) if len(data) else None,
                minimum_non_shared_m=float(meaningful.abs_dz_m.min()) if len(meaningful) else None,
                maximum_all_m=float(data.abs_dz_m.max()) if len(data) else None)


def write_point_shapes(data,folder):
    attributes,points,levels,level_points = [],[],[],[]
    for row in data.itertuples():
        attributes.append(dict(cross_id=row.cross_id,route_a=row.route_a,route_b=row.route_b,
                               z_a_m=row.z_a_m,z_b_m=row.z_b_m,abs_dz_m=row.abs_dz_m,
                               relation=row.xy_relation,shared=row.shared_endpoint_only))
        points.append(Point(row.x_5186,row.y_5186,(row.z_a_m+row.z_b_m)/2))
        for role,route,z in [('A',row.route_a,row.z_a_m),('B',row.route_b,row.z_b_m)]:
            levels.append(dict(cross_id=row.cross_id,role=role,route_id=route,alt=z,abs_dz_m=row.abs_dz_m))
            level_points.append(Point(row.x_5186,row.y_5186,z))
    if not attributes:
        return
    gpd.GeoDataFrame(attributes,geometry=points,crs=5186).to_file(folder/'intersection_points_z.shp',
        engine='pyogrio',encoding='UTF-8',geometry_type='Point Z')
    gpd.GeoDataFrame(levels,geometry=level_points,crs=5186).to_file(folder/'intersection_levels_z.shp',
        engine='pyogrio',encoding='UTF-8',geometry_type='Point Z')


def draw_map(data,routes,folder,title):
    fig,ax = plt.subplots(figsize=(12,11))
    used = set()
    for route in routes:
        if route.key in used:
            continue
        used.add(route.key)
        color = '#777777' if route.operation=='Short_Final' else ('#1376bb' if route.operation=='STAR' else '#e68522')
        ax.plot(*route.xyz[:,:2].T,color=color,alpha=.65,lw=.9)
    if len(data):
        colors = ax.scatter(data.x_5186,data.y_5186,c=data.abs_dz_m,cmap='viridis',s=30,zorder=5,
                            vmin=0,vmax=max(1.,float(data.abs_dz_m.max())))
        fig.colorbar(colors,ax=ax,label='Absolute vertical separation (m)',shrink=.7)
        shared = data[data.shared_endpoint_only]
        if len(shared):
            ax.scatter(shared.x_5186,shared.y_5186,facecolors='none',edgecolors='orange',marker='s',s=65,
                       zorder=6,label='Shared XY endpoints')
            ax.legend(fontsize=8)
        if len(data)<=80:
            for r in data.itertuples():
                ax.annotate(r.cross_id,(r.x_5186,r.y_5186),xytext=(3,3),textcoords='offset points',fontsize=6)
    else:
        ax.text(.5,.5,'No XY intersections',transform=ax.transAxes,ha='center')
    ax.set(aspect='equal',title=title,xlabel='EPSG:5186 X (m)',ylabel='EPSG:5186 Y (m)')
    ax.ticklabel_format(style='plain',useOffset=False)
    ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(folder/'map_intersections.png',dpi=180)
    plt.close(fig)


def draw_separation(data,folder,title,page_size):
    sorted_data = data.sort_values(['abs_dz_m','cross_id'])
    fig,ax = plt.subplots(figsize=(max(11,min(20,len(data)*.20)),6))
    if len(data):
        xpos = np.arange(len(data))
        shared = sorted_data.shared_endpoint_only.to_numpy(dtype=bool)
        ax.bar(xpos,sorted_data.abs_dz_m,color=np.where(shared,'#b6b6b6','#1677bb'))
        if len(data)<=100:
            ax.set_xticks(xpos,sorted_data.cross_id,rotation=90,fontsize=6)
        ax.text(.02,.97,'Grey: shared XY endpoints; blue: other intersections / overlaps',
                transform=ax.transAxes,va='top',fontsize=9)
    else:
        ax.text(.5,.5,'No XY intersections',transform=ax.transAxes,ha='center')
    ax.set(title=title+' / absolute Z gap',xlabel='Intersection IDs, sorted by gap',ylabel='Vertical separation (m)')
    ax.grid(axis='y',alpha=.25)
    fig.tight_layout()
    fig.savefig(folder/'vertical_separation.png',dpi=180)
    plt.close(fig)
    pages = []
    for page,first in enumerate(range(0,len(data),page_size),1):
        chunk = data.iloc[first:first+page_size]
        fig,ax = plt.subplots(figsize=(max(10,len(chunk)*.65),7))
        x = np.arange(len(chunk))
        ax.vlines(x,chunk.z_a_m,chunk.z_b_m,color='#666666',linewidth=2)
        ax.scatter(x,chunk.z_a_m,color='#1677bb',marker='o',label='Route A height')
        ax.scatter(x,chunk.z_b_m,color='#df7f20',marker='^',label='Route B height')
        for k,r in enumerate(chunk.itertuples()):
            ax.annotate(f'{r.abs_dz_m:.3f} m',(k,(r.z_a_m+r.z_b_m)/2),xytext=(4,2),
                        textcoords='offset points',fontsize=7,rotation=90)
        ax.set_xticks(x,chunk.cross_id,rotation=60,ha='right',fontsize=8)
        ax.set(title=f'{title} / crossing heights page {page}',xlabel='CSV cross_id',ylabel='Original geometry Z (m)')
        ax.legend()
        ax.grid(alpha=.2)
        fig.tight_layout()
        filename = f'crossing_heights_{page:02d}.png'
        fig.savefig(folder/filename,dpi=180)
        plt.close(fig)
        pages.append(filename)
    return pages


def draw_3d(data,routes,folder,title):
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111,projection='3d')
    seen = set()
    for route in routes:
        if route.key not in seen:
            ax.plot(*route.xyz.T,alpha=.45,lw=.8)
            seen.add(route.key)
    for r in data.itertuples():
        ax.plot([r.x_5186]*2,[r.y_5186]*2,[r.z_a_m,r.z_b_m],color='crimson',lw=1.4)
        ax.scatter([r.x_5186]*2,[r.y_5186]*2,[r.z_a_m,r.z_b_m],s=5,color='crimson')
    ax.set(title=title+' / vertical links at XY intersections',xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)')
    fig.tight_layout()
    fig.savefig(folder/'vertical_links_3d.png',dpi=170)
    plt.close(fig)


def save_group(data,a_routes,b_routes,folder,title,page_size,write_shapes):
    folder.mkdir(parents=True,exist_ok=True)
    data.to_csv(folder/'intersections.csv',index=False,encoding='utf-8-sig',float_format='%.9f')
    data[~data.shared_endpoint_only].to_csv(folder/'non_shared_endpoint_intersections.csv',
        index=False,encoding='utf-8-sig',float_format='%.9f')
    data[data.shared_endpoint_only].to_csv(folder/'shared_endpoints.csv',index=False,
        encoding='utf-8-sig',float_format='%.9f')
    if len(data):
        pairs = data.groupby(['route_a','label_a','route_b','label_b'],as_index=False).agg(
            evaluations=('cross_id','count'),min_gap_m=('abs_dz_m','min'),max_gap_m=('abs_dz_m','max'))
    else:
        pairs = pd.DataFrame(columns=['route_a','label_a','route_b','label_b','evaluations','min_gap_m','max_gap_m'])
    pairs.to_csv(folder/'route_pair_summary.csv',index=False,encoding='utf-8-sig',float_format='%.9f')
    draw_map(data,a_routes+b_routes,folder,title)
    pages = draw_separation(data,folder,title,page_size)
    draw_3d(data,a_routes+b_routes,folder,title)
    if write_shapes:
        write_point_shapes(data,folder)
    stats = summary_for(data,folder.name,len(a_routes),len(b_routes))
    stats['height_png_pages'] = pages
    (folder/'summary.json').write_text(json.dumps(stats,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    return stats


def report(output,summary):
    sections = []
    labels = dict(Short_Final_internal='Short_Final 내부',Vertical_routes_internal='Approach·Departure 내부',
                  Short_Final_vs_vertical_routes='Short_Final과 Approach·Departure 사이')
    def metric(value):
        return '교점 없음' if value is None else f'{value:.3f}m'
    for name,stats in summary['groups'].items():
        pages = ''.join(f'<details><summary>{page}</summary><img src="{name}/{page}"></details>' for page in stats['height_png_pages'])
        sections.append(f'''<h2>{labels[name]}</h2><p>교차한 선 쌍 {stats['intersecting_route_pairs']}개,
교점·중첩 평가 위치 {stats['intersection_evaluations']}개.
공유 끝점 {stats['shared_endpoint_evaluations']}개. 전체 최소 수직 간격 {metric(stats['minimum_all_m'])},
공유 끝점만 제외한 최소 간격 {metric(stats['minimum_non_shared_m'])}.</p>
<p><a href="{name}/intersections.csv">모든 교점 수치 CSV</a> ·
<a href="{name}/non_shared_endpoint_intersections.csv">공유 끝점 제외 CSV</a> ·
<a href="{name}/shared_endpoints.csv">공유 끝점 CSV</a> ·
<a href="{name}/route_pair_summary.csv">선 쌍별 요약</a></p>
<img src="{name}/map_intersections.png"><img src="{name}/vertical_separation.png">
<details><summary>3D 수직 간격</summary><img src="{name}/vertical_links_3d.png"></details>{pages}''')
    html = f'''<!doctype html><html lang="ko"><meta charset="utf-8"><title>Vertical Separation Gauge</title>
<style>body{{font-family:Malgun Gothic,sans-serif;max-width:1250px;margin:32px auto;line-height:1.7}}img{{max-width:100%}}.note{{background:#fff1d7;padding:15px}}</style>
<h1>교차점 수직 Separation 측정</h1>
<p>Short_Final과 최적화 항로를 내부 교차·상호 교차로 구분해 측정했습니다.
Approach·Departure 파일만 읽으며 중복 combined 파일과 점 SHP는 다시 포함하지 않았습니다.</p>
<p class="note">|ZA−ZB|는 같은 XY 위치에서 각 원본 선분의 Z를 따로 보간한 값입니다.
교점 기하의 혼합 Z를 사용하지 않습니다. 입력의 Z 수직 기준·단위가 같다고 가정하며 Z 기준은 변환하지 않습니다.
두 선의 끝점이 같은 XY에서 만나는 위치는 따로 표시합니다. 이것이 의도된 연결인지는 코드가 단정하지 않습니다.
수치만 측정하며 70m 등의 안전 기준이나 실제 운항 시각을 적용하지 않습니다.</p>
<p>중첩 구간은 선분 경계와 높이 차가 0이 되는 위치를 모두 평가합니다. 따라서 평가 위치 수는
단순한 고유 XY 교점 수와 다를 수 있습니다. 자기교차한 선에 높이가 여러 개 있으면 각 선분의 높이를 보존합니다.
XY를 자동 스냅하지 않으며, 허용오차 {summary['xy_tolerance_m']:g}m는 끝점 분류·중복 제거에만 사용합니다.</p>
<p>점 SHP의 intersection_points_z는 표시용 중간 높이, intersection_levels_z는 각 선의 실제 높이 점 두 개입니다.
CSV와 그림의 ZA·ZB가 실제 비교값이며, 3D 그림은 축별 범위가 달라 실제 종횡비 도면은 아닙니다.</p>
{''.join(sections)}
<h2>재실행</h2><pre>python -B 수정본_20261103/코드/Vertical_Separation_gauge.py
python -B 수정본_20261103/코드/Vertical_Separation_gauge.py --self-check</pre></html>'''
    (output/'report.html').write_text(html,encoding='utf-8')
    (output/'README.md').write_text('''# Vertical Separation Gauge

report.html에서 전체 결과를 확인할 수 있습니다.

- Short_Final_internal: Short_Final 선들 사이
- Vertical_routes_internal: Approach/Departure 선들 사이
- Short_Final_vs_vertical_routes: 두 자료 사이
- all_intersections.csv: 모든 분석 수치 (UTF-8 BOM)
- summary.csv / summary.json: 분석별 건수와 최소 수직 간격
- 각 폴더 crossing_heights_XX.png: 각 교점의 두 Z와 절대 높이 차
- 각 폴더 map_intersections.png / vertical_separation.png / vertical_links_3d.png

shared_endpoint_only=True는 두 선(part)의 끝점이 같은 XY에서 만난다는 뜻입니다.
공유 끝점의 0m 간격을 일반적인 내부 교차와 구분하되 원자료에서 삭제하지 않습니다.
중첩 선분은 공유 끝점을 갖더라도 중첩 분석에 남깁니다.
안전 등급·충돌 시각·기체 성능을 판정하는 코드는 아닙니다.
''',encoding='utf-8')


def self_check():
    def route(key,xyz):
        return Route(key,key,key+'.shp',0,0,key,'test',np.asarray(xyz,dtype=float))
    a = route('A',[[0,0,350],[10,0,350]])
    b = route('B',[[5,-5,300],[5,5,300]])
    events = pair_intersections(a,b)
    assert len(events)==1 and events[0]['abs_dz_m']==50 and not events[0]['shared_endpoint_only']
    a = route('A',[[0,0,0],[10,0,10]])
    b = route('B',[[0,0,10],[10,0,0]])
    events = pair_intersections(a,b)
    assert len(events)==3 and min(e['abs_dz_m'] for e in events)==0
    assert any(e['event_type']=='overlap_zero' and e['x_5186']==5 for e in events)
    a = route('A',[[0,0,10],[1,0,20]])
    b = route('B',[[1,0,20],[1,1,30]])
    events = pair_intersections(a,b)
    assert len(events)==1 and events[0]['shared_endpoint_only'] and events[0]['abs_dz_m']==0
    a = route('A',[[-1,-1,0],[1,1,10],[-1,1,20],[1,-1,30]])
    b = route('B',[[-2,0,100],[2,0,100]])
    events = pair_intersections(a,b)
    assert len(events)==2 and sorted(e['abs_dz_m'] for e in events)==[75,95]
    a = route('A',[[0,0,0],[1,0,1]])
    b = route('B',[[0,1,0],[1,1,1]])
    assert not pair_intersections(a,b)
    empty = analyze_group([a],[b],'test','T',1e-6)
    assert empty.empty and list(empty[~empty.shared_endpoint_only].columns)==EVENT_COLUMNS
    print('Self-check passed: crossing, shared endpoint, overlap zero, multiple Z at one XY, disjoint.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--short-final',type=Path,default=DEFAULT_FINAL)
    parser.add_argument('--routes-dir',type=Path,default=DEFAULT_ROUTES)
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--xy-tolerance',type=float,default=1e-6)
    parser.add_argument('--page-size',type=int,default=20)
    parser.add_argument('--skip-shp',action='store_true')
    parser.add_argument('--self-check',action='store_true')
    args = parser.parse_args()
    if args.self_check:
        self_check()
        return
    if not math.isfinite(args.xy_tolerance) or args.xy_tolerance<=0 or args.page_size<1:
        raise ValueError('Tolerance must be positive/finite; page size must be >=1.')
    route_paths = route_input_paths(args.routes_dir)
    inputs = [args.short_final.resolve()]+[p for p,_,_ in route_paths]
    hashes = file_hashes(inputs)
    final = read_routes(args.short_final,'SF','Short_Final')
    traffic = [r for path,tag,op in route_paths for r in read_routes(path,tag,op)]
    output = args.output.resolve()
    output.mkdir(parents=True,exist_ok=True)
    groups = [('Short_Final_internal','SF',final,final,True,'Short Final: internal intersections'),
              ('Vertical_routes_internal','RT',traffic,traffic,True,'Approach / Departure: internal intersections'),
              ('Short_Final_vs_vertical_routes','MX',final,traffic,False,'Short Final vs Approach / Departure')]
    summary = dict(input_files=[str(p) for p in inputs],input_sha256=hashes,horizontal_crs='EPSG:5186',
                   z_datum_conversion=False,z_assumption='Both inputs use the same vertical datum and metres.',
                   xy_tolerance_m=args.xy_tolerance,short_final_route_parts=len(final),
                   approach_departure_route_parts=len(traffic),groups={})
    frames = []
    for name,prefix,a,b,internal,title in groups:
        data = analyze_group(a,b,name,prefix,args.xy_tolerance,internal)
        summary['groups'][name] = save_group(data,a,b,output/name,title,args.page_size,not args.skip_shp)
        frames.append(data)
        print(json.dumps(summary['groups'][name],ensure_ascii=False),flush=True)
    pd.concat(frames,ignore_index=True).to_csv(output/'all_intersections.csv',index=False,
        encoding='utf-8-sig',float_format='%.9f')
    pd.DataFrame([{k:v for k,v in row.items() if k!='height_png_pages'} for row in summary['groups'].values()]).to_csv(
        output/'summary.csv',index=False,encoding='utf-8-sig',float_format='%.9f')
    if file_hashes(inputs)!=hashes:
        raise RuntimeError('Input files changed during analysis.')
    summary['original_inputs_unchanged'] = True
    (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    report(output,summary)
    print(f'Complete: {output}',flush=True)


if __name__=='__main__':
    main()
