"""송파구 라인-건물 이격 분석. 실행: python 코드/지형/건물_높이_간격.py

필요 패키지: geopandas, shapely, numpy, pandas, matplotlib.
라인 Z(없으면 A16)와 건물 A16은 모두 AGL(m). DEM을 사용하지 않는다.
높이 간격은 AGL 값의 차이며, 서로 다른 지점의 절대고도 차이는 아니다.
"""
from pathlib import Path
import argparse
import json

import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely import force_2d, get_coordinates, make_valid
from shapely.geometry import Polygon, LineString
from shapely.ops import nearest_points

ROOT = Path(__file__).resolve().parents[2]


def load_ring(path, crs):
    frame = gpd.read_file(path)
    if frame.crs is None or len(frame) != 1:
        raise ValueError(f'{path}: CRS와 단일 닫힌 라인이 필요합니다.')
    frame = frame.to_crs(crs)
    geom = frame.geometry.iloc[0]
    if geom.geom_type != 'LineString' or not geom.is_ring:
        raise ValueError(f'{path}: 단순 폐합 LineString이 필요합니다.')
    z = get_coordinates(geom, include_z=True)[:, 2]
    if np.isfinite(z).all():
        if np.ptp(z) > 1e-6:
            raise ValueError(f'{path}: 가변 Z는 지원하지 않습니다.')
        height = float(z[0])
        if 'A16' in frame and abs(float(frame.A16.iloc[0]) - height) > 1e-6:
            raise ValueError(f'{path}: Z와 A16이 다릅니다.')
    else:
        height = float(frame.A16.iloc[0])
    if not np.isfinite(height):
        raise ValueError('라인 고도가 유효하지 않습니다.')
    line = force_2d(geom)
    return line, Polygon(line), height


def plot_target_buildings(buildings, rings, target_height, out):
    """지정 A16 건물 각각에 대해 두 라인까지의 정확한 2D 최단거리 도면 저장."""
    targets = buildings.loc[np.isclose(buildings.building_height_m, target_height, rtol=0, atol=1e-6)]
    rows = []
    for _, building in targets.iterrows():
        fig, axes = plt.subplots(1, 2, figsize=(15, 8))
        for ax, (name, (line, interior, altitude)) in zip(axes, rings.items()):
            a, b = nearest_points(building.geometry, line)
            distance = building.geometry.distance(line)
            gap = altitude - building.building_height_m
            rows.append(dict(source_row=int(building.source_row), A16=building.building_height_m,
                             line=name, line_height_m=altitude, vertical_gap_m=gap,
                             xy_distance_m=distance, fully_inside=interior.covers(building.geometry),
                             building_x=a.x, building_y=a.y, line_x=b.x, line_y=b.y,
                             geometry=LineString([a, b])))
            x, y = line.xy
            ax.plot(x, y, color='#2670ac', linewidth=1.8, label=f'{name} 라인')
            gpd.GeoSeries([building.geometry], crs=buildings.crs).plot(ax=ax, color='#df6239', edgecolor='#8b301b')
            ax.plot([a.x, b.x], [a.y, b.y], color='#c62728', linewidth=2, label='XY 최단 연결선')
            ax.scatter([a.x, b.x], [a.y, b.y], s=24, color='#c62728', zorder=5)
            ax.annotate(f'{distance:,.2f} m', ((a.x+b.x)/2, (a.y+b.y)/2),
                        xytext=(12, 10), textcoords='offset points', fontsize=12,
                        bbox=dict(facecolor='white', edgecolor='#c62728', alpha=.95))
            ax.annotate(f'A16 = {building.building_height_m:g} m', (a.x, a.y),
                        xytext=(10, 18), textcoords='offset points', fontsize=10)
            # 각 패널은 해당 원과 건물을 모두 포함해 거리의 위치 관계를 보여준다.
            xmin, ymin, xmax, ymax = line.union(building.geometry).bounds
            pad = max(xmax-xmin, ymax-ymin) * .08
            ax.set_xlim(xmin-pad, xmax+pad); ax.set_ylim(ymin-pad, ymax+pad)
            ax.set_aspect('equal'); ax.ticklabel_format(style='plain', useOffset=False)
            ax.set_xlabel('Easting (m)'); ax.set_ylabel('Northing (m)')
            ax.set_title(f'{name}: XY 최단거리 {distance:,.2f} m\n라인 AGL {altitude:g} m / AGL 높이 차 {gap:+.2f} m')
            ax.legend(loc='upper right'); ax.grid(alpha=.15)
        fig.suptitle(f'A16 {target_height:g}m 건물과 내곽·외곽의 평면 최단거리 (원본 행 {int(building.source_row)})\n건물 폴리곤과 라인 사이의 거리 / DEM 미사용', fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, .93))
        fig.savefig(out / f'A16_{target_height:g}m_행{int(building.source_row)}_XY최단거리.png', dpi=200)
        plt.close(fig)
    columns = ['source_row', 'A16', 'line', 'line_height_m', 'vertical_gap_m', 'xy_distance_m',
               'fully_inside', 'building_x', 'building_y', 'line_x', 'line_y']
    pd.DataFrame(rows, columns=columns).to_csv(out / f'A16_{target_height:g}m_XY최단거리.csv', index=False, encoding='utf-8-sig')
    if rows:
        gpd.GeoDataFrame(rows, crs=buildings.crs).to_file(out / f'A16_{target_height:g}m_최단연결선.geojson', driver='GeoJSON')
    else:
        print(f'주의: 분석 범위에 A16={target_height:g}m 건물이 없어 대상 PNG를 생성하지 않았습니다.')
    return len(targets)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=ROOT / '자료/안전성/송파구')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--buffer-m', type=float, default=100)
    parser.add_argument('--vertical-threshold-m', type=float, default=100)
    parser.add_argument('--building-file', default='송파구_건물_수정_update.shp', help='입력 폴더 기준 파일명 또는 절대 경로')
    parser.add_argument('--target-height-m', type=float, default=555, help='별도 XY 거리 PNG를 생성할 건물 A16')
    args = parser.parse_args()
    if args.buffer_m <= 0 or args.vertical_threshold_m < 0:
        parser.error('버퍼는 양수, 수직 임계값은 0 이상이어야 합니다.')
    out = args.output_dir or args.input_dir / 'Separation_확보'
    out.mkdir(parents=True, exist_ok=True)
    buildings = gpd.read_file(args.input_dir / args.building_file)
    crs = buildings.crs
    if crs is None or not crs.is_projected or any(abs(a.unit_conversion_factor - 1) > 1e-9 for a in crs.axis_info):
        raise ValueError('건물 CRS는 미터 단위 투영좌표계여야 합니다.')
    buildings['source_row'] = np.arange(len(buildings))
    buildings.geometry = buildings.geometry.map(lambda g: make_valid(force_2d(g)) if g is not None else None)
    buildings['building_height_m'] = pd.to_numeric(buildings.A16, errors='coerce')
    buildings.loc[~np.isfinite(buildings.building_height_m) | (buildings.building_height_m <= 0), 'building_height_m'] = np.nan
    rings = {name: load_ring(args.input_dir / f'{name}.shp', crs) for name in ['외곽', '내곽']}
    scope = rings['외곽'][1].buffer(args.buffer_m).union(rings['내곽'][1].buffer(args.buffer_m))
    buildings = buildings.loc[buildings.intersects(scope)].copy()
    buildings['height_reliable'] = buildings.building_height_m.notna()
    for name, (line, _, _) in rings.items():
        buildings[f'{name}_distance_m'] = buildings.distance(line)
    plt.rcParams['font.family'] = 'Malgun Gothic'
    plt.rcParams['axes.unicode_minus'] = False
    summaries = []
    for name, (line, interior, altitude) in rings.items():
        buffer = line.buffer(args.buffer_m, quad_segs=64)
        data = buildings.loc[buildings.intersects(interior) | buildings.intersects(buffer)].copy()
        data['in_buffer'] = data.intersects(buffer)
        data['intersects_interior'] = data.intersects(interior)
        data['fully_inside'] = data.geometry.map(interior.covers)
        data['line_height_m'] = altitude
        data['vertical_gap_m'] = altitude - data.building_height_m
        data['status'] = np.select(
            [~data.height_reliable, data.vertical_gap_m < 0,
             data.vertical_gap_m <= args.vertical_threshold_m],
            ['높이확인불가', '설정고도초과', '수직간격임계값이하'], default='수직간격확보')
        risk = data.loc[data.height_reliable & (data.vertical_gap_m <= args.vertical_threshold_m)].copy()
        columns = ['source_row', 'A0', 'A21', 'A16', 'building_height_m',
                   'height_reliable', 'line_height_m',
                   'vertical_gap_m', '외곽_distance_m', '내곽_distance_m', 'in_buffer',
                   'intersects_interior', 'fully_inside', 'status']
        columns = [c for c in columns if c in data]
        for suffix, frame in [('전체', data), (f'{args.buffer_m:g}m버퍼', data.loc[data.in_buffer]), ('주의건물', risk)]:
            frame[columns].to_csv(out / f'{name}_{suffix}.csv', index=False, encoding='utf-8-sig')
        data[columns + ['geometry']].to_file(out / f'{name}_분석.geojson', driver='GeoJSON')
        connectors = []
        for _, row in risk.iterrows():
            for target, (target_line, _, _) in rings.items():
                a, b = nearest_points(row.geometry, target_line)
                connectors.append({'source_row': row.source_row, 'target': target,
                                   'distance_m': a.distance(b), 'geometry': LineString([a, b])})
        if connectors:
            gpd.GeoDataFrame(connectors, crs=crs).to_file(out / f'{name}_주의건물_최단연결선.geojson', driver='GeoJSON')
        near = data.loc[data.in_buffer & data.height_reliable]
        summary = dict(line=name, line_height_m=altitude, scope_count=len(data),
                       buffer_count=int(data.in_buffer.sum()), buffer_valid_count=len(near),
                       buffer_unknown_count=int((data.in_buffer & ~data.height_reliable).sum()),
                       buffer_gap_min_m=float(near.vertical_gap_m.min()),
                       buffer_gap_median_m=float(near.vertical_gap_m.median()),
                       buffer_gap_max_m=float(near.vertical_gap_m.max()),
                       scope_exceed_count=int((data.status == '설정고도초과').sum()),
                       scope_within_threshold_count=int((data.status == '수직간격임계값이하').sum()),
                       scope_unknown_count=int((~data.height_reliable).sum()),
                       risk_min_xy_m=float(risk[f'{name}_distance_m'].min()) if len(risk) else None)
        summaries.append(summary)
        fig, (ax, hist) = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={'width_ratios': [1.5, 1]})
        data.plot(ax=ax, color='#c9cdd2', linewidth=0)
        gpd.GeoSeries([buffer], crs=crs).plot(ax=ax, color='#82c6e5', alpha=.35)
        if len(near):
            near.plot(ax=ax, column='vertical_gap_m', cmap='viridis', legend=True,
                      legend_kwds={'label': 'AGL 높이 차 (m)'}, linewidth=0)
        for other, (other_line, _, _) in rings.items():
            gpd.GeoSeries([other_line], crs=crs).plot(ax=ax, color='black', linewidth=1, linestyle='--' if other != name else '-')
        if len(risk):
            risk.plot(ax=ax, color='red', edgecolor='darkred')
            for _, row in risk.iterrows():
                a, b = nearest_points(row.geometry, line)
                ax.plot([a.x, b.x], [a.y, b.y], color='red', linewidth=.9)
                ax.annotate(f"행 {row.source_row}\nXY {a.distance(b):.1f}m\n수직 {row.vertical_gap_m:.1f}m", (a.x, a.y), fontsize=8)
        ax.set_title(f'{name} AGL {altitude:g}m / 버퍼 {args.buffer_m:g}m\n빨강: 범위 내 설정 높이 초과 또는 AGL 높이 차 {args.vertical_threshold_m:g}m 이하')
        ax.set_xlabel('Easting (m)'); ax.set_ylabel('Northing (m)'); ax.set_aspect('equal')
        hist.hist(near.vertical_gap_m, bins=35, color='#3679a8')
        hist.axvline(0, color='red', label='설정 고도')
        hist.axvline(args.vertical_threshold_m, color='orange', label='수직 간격 임계값')
        hist.set(title=f'버퍼 교차 건물 {int(data.in_buffer.sum()):,}개\n높이 유효 {len(near):,}개 / 미확인 {summary["buffer_unknown_count"]:,}개', xlabel='라인 AGL - 건물 A16 (AGL), m', ylabel='건물 수')
        hist.legend(); fig.tight_layout(); fig.savefig(out / f'{name}_높이간격.png', dpi=180); plt.close(fig)
    pd.DataFrame(summaries).to_csv(out / '분석요약.csv', index=False, encoding='utf-8-sig')
    target_count = plot_target_buildings(buildings, rings, args.target_height_m, out)
    notes = f'''# 송파구 건물 높이 간격 분석

- 입력 건물: {args.building_file}.
- 높이 간격(vertical_gap_m) = 라인 AGL − 건물 A16(AGL). DEM을 읽거나 더하지 않음.
- 사용자 확인에 따라 라인과 건물 모두 AGL(m). 지면 표고가 다른 두 지점의 절대고도 차이를 의미하지 않음.
- A16이 0 이하/결측이면 높이확인불가. 미확인 건물은 안전하다고 판정하지 않음.
- 별도 거리 도면: A16={args.target_height_m:g}m 건물 {target_count}개. 높이 일치 허용오차 0.000001m.
- 대상 건물별 PNG 및 CSV에 양쪽 라인의 XY 최단거리, AGL 높이 차, 최단 연결점 좌표를 기록.
- 버퍼 분석: 라인 양쪽 {args.buffer_m:g}m 버퍼와 건물 폴리곤이 2D 교차하는 건물.
- 전체 범위: 폐합 라인 내부와 교차하거나 버퍼와 교차하는 건물. 내부 판정은 원의 원판이며 내곽을 뺀 환형 영역이 아님.
- 주의건물: 전체 범위에서 수직 간격이 {args.vertical_threshold_m:g}m 이하인 유효 건물(고도 초과 포함).
- XY 거리: 건물 폴리곤과 라인 사이의 최단거리. 라인과 겹치면 0m. 모든 CSV에 내곽/외곽 거리를 함께 기록.
- source_row는 원본 SHP의 0부터 시작하는 행 번호. A0/A21도 함께 제공.
- GeoJSON은 원본 미터 좌표계 보존용 GIS 자료이며 일반 웹 지도의 WGS84 GeoJSON이 아님.
- CSV는 UTF-8 BOM. PNG의 회색 건물은 배경, 빨강은 주의건물, 색상 건물은 버퍼 내 유효 높이 분석 대상.
- 법적 안전 판정이 아니라 입력 자료와 위 가정에 따른 기하학적 이격 분석.
'''
    (out / '분석방법.md').write_text(notes, encoding='utf-8')
    print(json.dumps(summaries, ensure_ascii=False, indent=2))
    print(f'저장 완료: {out}')


if __name__ == '__main__':
    main()
