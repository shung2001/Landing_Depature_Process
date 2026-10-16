"""IAF–FAF 경로와 내곽_2000(1500ft)의 높이 간격을 CSV/PNG로 저장한다.

실행: python 코드/지형/FAF_외곽이동시_1500ft와의_간격.py
필요 패키지: geopandas, shapely, numpy, pandas, matplotlib.
원본 SHP는 변경하지 않는다. 실제 geometry Z(m)를 사용하며 DEM/AGL 변환은 하지 않는다.
교차점은 XY로 구하고 경로 Z를 누적 거리로 따로 보간한다.
Shapely의 3D 교차 결과 Z는 사용하지 않는다(서로 다른 높이를 평균할 수 있음).
"""

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from shapely import force_2d, get_coordinates
from shapely.geometry import Point


ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / '자료/안전성/송파구'
NAMES = ('IAF_FAF', 'IAF_FAF_Out')
OUTPUT_NAME = 'FAF_외곽이동시_1500ft와의_간격'


def coordinates(geometry, label):
    if geometry is None or geometry.is_empty or geometry.geom_type != 'LineString':
        raise ValueError(f'{label}: 비어 있지 않은 단일 LineStringZ가 필요합니다.')
    xyz = get_coordinates(geometry, include_z=True)
    if len(xyz) < 2 or not np.isfinite(xyz).all():
        raise ValueError(f'{label}: 유효한 XYZ 좌표가 필요합니다.')
    return xyz


def intersection_points(geometry):
    if geometry.is_empty:
        return []
    if geometry.geom_type == 'Point':
        return [geometry]
    if geometry.geom_type in ('MultiPoint', 'GeometryCollection'):
        return [p for part in geometry.geoms for p in intersection_points(part)]
    raise ValueError('경로와 내곽 선이 구간으로 겹칩니다. 점 교차 분석으로 처리할 수 없습니다.')


def analyze_route(row, layer, number, ring, reference_z, tolerance):
    xyz = coordinates(row.geometry, f'{layer} 행 {number}')
    line = force_2d(row.geometry)
    distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(xyz[:, :2], axis=0), axis=1))]
    if np.any(np.diff(distance) <= 0):
        raise ValueError(f'{layer} 행 {number}: 중복 XY 꼭짓점/수직 선분은 지원하지 않습니다.')
    candidates = [(line.project(p), 'intersection', p.distance(ring))
                  for p in intersection_points(line.intersection(ring))]
    # 원의 이산화 및 좌표 변환으로 수 mm 차이가 나는 끝점은 별도 유형으로 기록.
    for at, xy in ((0.0, xyz[0, :2]), (distance[-1], xyz[-1, :2])):
        separation = Point(xy).distance(ring)
        if separation <= tolerance:
            candidates.append((at, 'endpoint_within_tolerance', separation))
    # 가까운 교차점은 한 위치로 묶되 끝점의 정확한 Z를 우선한다.
    selected = []
    for candidate in sorted(candidates, key=lambda item: item[1] != 'endpoint_within_tolerance'):
        if not any(abs(candidate[0] - other[0]) <= tolerance for other in selected):
            selected.append(candidate)
    selected.sort()
    base = dict(layer=layer, feature_row=number, source_id=row.get('id', number))
    contacts = []
    for index, (at, kind, xy_distance) in enumerate(selected, 1):
        z = float(np.interp(at, distance, xyz[:, 2]))
        xy = line.interpolate(at)
        gap = z - reference_z
        contacts.append(dict(**base, contact_no=index, contact_type=kind,
                             distance_m=at, progress_pct=100 * at / distance[-1],
                             x=xy.x, y=xy.y, route_z_m=z, ring_z_m=reference_z,
                             signed_gap_m=gap, absolute_gap_m=abs(gap),
                             xy_distance_to_ring_m=xy_distance))
    sample_at = np.unique(np.r_[distance, np.linspace(0, distance[-1], 201),
                                [contact['distance_m'] for contact in contacts]])
    samples = []
    for at in sample_at:
        xy = line.interpolate(at)
        z = float(np.interp(at, distance, xyz[:, 2]))
        samples.append(dict(**base, distance_m=at, progress_pct=100 * at / distance[-1],
                            x=xy.x, y=xy.y, route_z_m=z, ring_reference_z_m=reference_z,
                            signed_height_difference_m=z - reference_z,
                            absolute_height_difference_m=abs(z - reference_z),
                            xy_distance_to_ring_m=xy.distance(ring)))
    summary = dict(**base, length_m=distance[-1], start_z_m=xyz[0, 2], end_z_m=xyz[-1, 2],
                   min_route_z_m=xyz[:, 2].min(), max_route_z_m=xyz[:, 2].max(),
                   ring_z_m=reference_z, min_signed_difference_m=xyz[:, 2].min() - reference_z,
                   max_signed_difference_m=xyz[:, 2].max() - reference_z,
                   contact_count=len(contacts),
                   min_contact_gap_m=min((c['absolute_gap_m'] for c in contacts), default=np.nan),
                   max_contact_gap_m=max((c['absolute_gap_m'] for c in contacts), default=np.nan))
    return dict(base=base, line=line, distance=distance, z=xyz[:, 2],
                contacts=contacts, samples=samples, summary=summary)


def configure_font():
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ('Malgun Gothic', 'NanumGothic', 'Noto Sans CJK KR'):
        if name in available:
            plt.rcParams['font.family'] = name
            break
    plt.rcParams['axes.unicode_minus'] = False


def plot_profiles(routes, layer, height, output):
    rows = (len(routes) + 3) // 4
    fig, axes = plt.subplots(rows, 4, figsize=(20, rows * 3.75), sharey=True, squeeze=False)
    for ax, route in zip(axes.flat, routes):
        at, z = route['distance'], route['z']
        ax.plot(at, z, color='#2166ac', lw=2, label='경로 Z')
        ax.axhline(height, color='#d6604d', ls='--', label=f'내곽 기준 {height:g}m')
        ax.fill_between(at, z, height, color='#92c5de', alpha=.25)
        for index, c in enumerate(route['contacts']):
            x, value = c['distance_m'], c['route_z_m']
            ax.plot(x, value, 'o', color='#7b3294', ms=4)
            ax.annotate('', xy=(x, value), xytext=(x, height),
                        arrowprops=dict(arrowstyle='<->', color='#7b3294', lw=1.4))
            # 위치별 값은 부호와 함께 표시하여 기준보다 위/아래를 구분.
            ax.annotate(f"ΔZ {c['signed_gap_m']:+.2f}m", xy=(x, (value + height) / 2),
                        xytext=(5 if x < at[-1] / 2 else -5, 12 + 12 * (index % 2)),
                        textcoords='offset points', fontsize=8,
                        ha='left' if x < at[-1] / 2 else 'right', color='#7b3294')
        s = route['summary']
        ax.set_title(f"행 {s['feature_row']} / id {s['source_id']} / 접점 {s['contact_count']}개", fontsize=10)
        ax.set_xlabel('저장된 시작점부터 누적 거리 (m)', fontsize=8)
        ax.set_ylabel('Z (m)', fontsize=9)
        ax.grid(alpha=.25)
        ax.margins(x=.14, y=.22)
    for ax in list(axes.flat)[len(routes):]:
        ax.set_visible(False)
    axes.flat[0].legend(fontsize=7, loc='lower left')
    fig.suptitle(f'{layer}: 내곽_2000과의 높이 간격\n'
                 '화살표 = XY 교차/허용거리 내 끝점의 높이 차이 (ΔZ = 경로 Z - 내곽 Z)', fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(output / f'{layer}_높이간격.png', dpi=180)
    plt.close(fig)


def plot_overview(routes, ring, height, output):
    fig, axes = plt.subplots(2, 2, figsize=(17, 13))
    for row, layer in enumerate(NAMES):
        subset = [r for r in routes if r['base']['layer'] == layer]
        ax, bars = axes[row]
        ax.plot(*ring.xy, color='#d6604d', lw=2, label=f'내곽_2000 Z={height:g}m')
        for route in subset:
            ax.plot(*route['line'].xy, lw=.8, color='#2166ac', alpha=.7)
            for contact in route['contacts']:
                ax.plot(contact['x'], contact['y'], 'o', color='#7b3294', ms=4)
        ax.set_aspect('equal')
        ax.set_title(f'{layer}: XY 경로와 접점')
        ax.set_xlabel('EPSG:5186 X (m)')
        ax.set_ylabel('EPSG:5186 Y (m)')
        ax.legend(fontsize=9)
        ax.ticklabel_format(style='plain', useOffset=False)
        positions = np.arange(len(subset))
        values = [r['summary']['min_contact_gap_m'] for r in subset]
        bars.bar(positions, values, color='#7b3294', alpha=.8)
        for x, value in zip(positions, values):
            bars.text(x, 0 if np.isnan(value) else value,
                      '없음' if np.isnan(value) else f'{value:.2f}',
                      ha='center', va='bottom', fontsize=8)
        bars.set_xticks(positions, [str(r['base']['feature_row']) for r in subset])
        bars.set_xlabel('피처 행 번호 (1부터, 중복 id와 구분)')
        bars.set_ylabel('접점의 최소 |ΔZ| (m)')
        bars.set_title('내곽과 만나는 위치의 최소 높이 간격')
        bars.set_ylim(0, max([v for v in values if np.isfinite(v)] + [1]) * 1.18)
        bars.grid(axis='y', alpha=.25)
    fig.suptitle('IAF–FAF와 1500ft 내곽의 높이 간격 (3차원 최단거리 아님)', fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, .96))
    fig.savefig(output / '전체_높이간격.png', dpi=180)
    plt.close(fig)


def closest_segment_to_ring(a0, a1, ring_xyz):
    """3D 선분 대 모든 내곽 선분의 최소거리. 내부 해와 네 경계 해를 비교한다."""
    u = a1 - a0
    start, v = ring_xyz[:-1], np.diff(ring_xyz, axis=0)
    w = a0 - start
    aa = float(u @ u)
    cc = np.einsum('ij,ij->i', v, v)
    if aa <= 0 or np.any(cc <= 0):
        raise ValueError('3D 거리 계산에 길이가 0인 선분이 있습니다.')
    bb, dd, ee = v @ u, w @ u, np.einsum('ij,ij->i', v, w)
    candidates = [(np.zeros(len(v)), np.clip(ee / cc, 0, 1)),
                  (np.ones(len(v)), np.clip((ee + bb) / cc, 0, 1)),
                  (np.clip(-dd / aa, 0, 1), np.zeros(len(v))),
                  (np.clip((bb - dd) / aa, 0, 1), np.ones(len(v)))]
    denominator = aa * cc - bb * bb
    nonparallel = denominator > 1e-12 * aa * cc
    s = np.divide(bb * ee - cc * dd, denominator,
                  out=np.zeros_like(cc), where=nonparallel)
    t = np.divide(aa * ee - bb * dd, denominator,
                  out=np.zeros_like(cc), where=nonparallel)
    valid = nonparallel & (s >= 0) & (s <= 1) & (t >= 0) & (t <= 1)
    candidates.append((s, t))
    best = None
    for index, (s, t) in enumerate(candidates):
        delta = w + s[:, None] * u - t[:, None] * v
        squared = np.einsum('ij,ij->i', delta, delta)
        if index == 4:
            squared[~valid] = np.inf
        at = int(np.argmin(squared))
        if best is None or squared[at] < best[0]:
            best = (float(squared[at]), a0 + s[at] * u, start[at] + t[at] * v[at],
                    float(s[at]), at + 1)
    return best


def analyze_buffer(route, ring_xyz, buffer_m):
    xyz = np.column_stack((np.asarray(route['line'].coords), route['z']))
    segments = []
    for index, (a, b) in enumerate(zip(xyz[:-1], xyz[1:])):
        squared, p, q, fraction, ring_segment = closest_segment_to_ring(a, b, ring_xyz)
        gap = float(np.sqrt(squared))
        at = route['distance'][index] + fraction * (route['distance'][index + 1] - route['distance'][index])
        segments.append(dict(**route['base'], segment_no=index + 1,
                             ring_segment_no=ring_segment, distance_along_route_m=at,
                             min_distance_3d_m=gap, horizontal_distance_m=float(np.linalg.norm(p[:2] - q[:2])),
                             signed_vertical_distance_m=float(p[2] - q[2]),
                             route_x=p[0], route_y=p[1], route_z=p[2],
                             ring_x=q[0], ring_y=q[1], ring_z=q[2],
                             combined_buffer_m=buffer_m, buffer_margin_m=gap - buffer_m,
                             buffer_overlap=gap <= buffer_m + 1e-8))
    closest = min(segments, key=lambda item: item['min_distance_3d_m'])
    route['buffer_segments'] = segments
    route['buffer_result'] = dict(closest, flagged_segment_count=sum(s['buffer_overlap'] for s in segments))


def plot_buffer(routes, ring, buffer_m, output):
    fig, axes = plt.subplots(2, 2, figsize=(18, 13))
    for row, layer in enumerate(NAMES):
        subset = [r for r in routes if r['base']['layer'] == layer]
        ax, bars = axes[row]
        ax.plot(*ring.xy, color='black', lw=2, label='내곽_2000')
        for route in subset:
            xy = np.asarray(route['line'].coords)
            ax.plot(xy[:, 0], xy[:, 1], color='#2166ac', lw=.8, alpha=.55)
            for index, segment in enumerate(route['buffer_segments']):
                if segment['buffer_overlap']:
                    ax.plot(xy[index:index + 2, 0], xy[index:index + 2, 1], color='#d73027', lw=2)
            c = route['buffer_result']
            ax.plot([c['route_x'], c['ring_x']], [c['route_y'], c['ring_y']], color='#762a83', lw=.8)
        ax.plot([], [], color='#d73027', lw=2, label='일부에 버퍼 겹침이 있는 선분')
        ax.set_title(f'{layer}: 겹침 선분의 XY 위치')
        ax.set_aspect('equal')
        ax.legend(fontsize=9)
        ax.ticklabel_format(style='plain', useOffset=False)
        values = [r['buffer_result']['min_distance_3d_m'] for r in subset]
        bars.bar(np.arange(len(subset)), values,
                 color=['#d73027' if r['buffer_result']['buffer_overlap'] else '#2166ac' for r in subset])
        bars.axhline(buffer_m, color='black', ls='--', label=f'합산 버퍼 {buffer_m:g}m')
        for i, value in enumerate(values):
            bars.text(i, value + 2, f'{value:.2f}', ha='center', fontsize=8)
        bars.set_xticks(np.arange(len(subset)), [str(r['base']['feature_row']) for r in subset])
        bars.set_xlabel('피처 행 번호 (1부터)')
        bars.set_ylabel('경로 전체와 내곽 경계의 3D 최단거리 (m)')
        bars.set_ylim(0, max(values + [buffer_m]) * 1.2)
        bars.set_title(f'{sum(r["buffer_result"]["buffer_overlap"] for r in subset)}/{len(subset)} 경로 버퍼 겹침')
        bars.grid(axis='y', alpha=.25)
        bars.legend()
    fig.suptitle(f'합산 {buffer_m:g}m 3D 거리 버퍼 판정 (각 선 {buffer_m / 2:g}m씩과 동등)\n'
                 '빨강: 최단거리 <= 합산 버퍼 / 파랑: 합산 버퍼 초과', fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, .94))
    fig.savefig(output / '3D_버퍼겹침_판정.png', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DIRECTORY)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--contact-tolerance', type=float, default=.01, help='끝점 접촉 허용 거리(m)')
    parser.add_argument('--combined-buffer', type=float, default=100, help='두 선의 버퍼 반경 합계(m)')
    args = parser.parse_args()
    if not np.isfinite(args.contact_tolerance) or args.contact_tolerance <= 0:
        raise ValueError('접촉 허용 거리는 0보다 큰 유한한 수여야 합니다.')
    if not np.isfinite(args.combined_buffer) or args.combined_buffer < 0:
        raise ValueError('합산 버퍼는 0 이상의 유한한 수여야 합니다.')
    base = args.input_dir.resolve()
    output = args.output_dir or base / '생성파일' / OUTPUT_NAME
    reference = gpd.read_file(base / '내곽_2000.shp')
    if len(reference) != 1 or reference.crs is None:
        raise ValueError('내곽_2000은 좌표계가 있는 단일 닫힌 선이어야 합니다.')
    reference = reference.to_crs(5186)
    xyz = coordinates(reference.geometry.iloc[0], '내곽_2000')
    ring = force_2d(reference.geometry.iloc[0])
    if not ring.is_ring or np.ptp(xyz[:, 2]) > 1e-6:
        raise ValueError('내곽_2000은 고도가 일정한 단순 폐합 선이어야 합니다.')
    height = float(xyz[0, 2])
    routes = []
    for layer in NAMES:
        frame = gpd.read_file(base / f'{layer}.shp')
        if frame.empty or frame.crs is None:
            raise ValueError(f'{layer}: 좌표계와 피처가 필요합니다.')
        frame = frame.to_crs(5186)
        for number, (_, row) in enumerate(frame.iterrows(), 1):
            routes.append(analyze_route(row, layer, number, ring, height, args.contact_tolerance))
    for route in routes:
        analyze_buffer(route, xyz, args.combined_buffer)
    # 계산을 모두 완료한 후 출력 폴더를 생성한다.
    output.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame([r['summary'] for r in routes])
    summary.to_csv(output / '경로별_간격요약.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame([c for r in routes for c in r['contacts']]).to_csv(
        output / '내곽_교차접점_높이간격.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame([s for r in routes for s in r['samples']]).to_csv(
        output / '경로별_높이프로파일.csv', index=False, encoding='utf-8-sig')
    configure_font()
    for layer in NAMES:
        plot_profiles([r for r in routes if r['base']['layer'] == layer], layer, height, output)
    plot_overview(routes, ring, height, output)
    buffer_summary = pd.DataFrame([r['buffer_result'] for r in routes])
    buffer_summary.to_csv(output / '3D_버퍼겹침_경로별판정.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame([s for r in routes for s in r['buffer_segments']]).to_csv(
        output / '3D_버퍼겹침_선분별판정.csv', index=False, encoding='utf-8-sig')
    plot_buffer(routes, ring, args.combined_buffer, output)
    text = (f'내곽_2000의 실제 Z: {height:.6f} m ({height / .3048:.3f} ft)\n'
            f'끝점 접촉 허용 거리: {args.contact_tolerance:g} m\n'
            'signed_gap_m = 경로 Z − 내곽 Z. 양수는 위, 음수는 아래.\n'
            'absolute_gap_m은 높이 간격의 크기이다.\n'
            '교차접점 CSV: XY 교차점 또는 허용 거리 내 끝점에서 각 자료의 Z를 독립 계산.\n'
            'endpoint_within_tolerance는 정확한 교차가 아니라 허용 거리 내 끝점이다.\n'
            '높이프로파일 CSV/음영: 내곽과 수평 위치가 다른 구간도 포함한 기준고도 차이.\n'
            '그 구간의 값은 내곽 선 자체와 같은 XY에서 측정한 간격이 아니다.\n'
            '전체 PNG의 막대는 각 경로의 여러 접점 중 최소 |ΔZ|를 표시한다.\n'
            '모든 접점의 간격은 개별 경로 PNG와 교차접점 CSV에 기록한다.\n'
            'feature_row는 1부터 시작하는 행 번호. 원본 id는 중복될 수 있다.\n'
            f'3D 버퍼 판정은 두 선의 합산 반경 {args.combined_buffer:g}m를 사용한다. 각 선에 100m씩 적용하지 않는다.\n'
            '경로와 내곽 경계를 구성하는 모든 3D 직선 선분 쌍의 최단거리를 계산한다.\n'
            'buffer_overlap=True는 3D 최단거리 <= 합산 반경 (수치 허용오차 1e-8m).\n'
            '연결을 의도한 끝점도 제외하지 않는다. 의도된 접속 여부는 별도 검토가 필요하다.\n'
            '선분별판정/빨간 선분은 선분 일부가 버퍼 내임을 의미하며 선분 전체의 침범을 뜻하지 않는다.\n'
            '이 기준은 등방성 3D 거리 버퍼이며 수평/수직 분리 기준이나 운항 안전성 판정이 아니다.\n'
            'Z는 원본 기준 그대로 사용. DEM/AGL 변환은 하지 않음.\n')
    (output / '분석설명.txt').write_text(text, encoding='utf-8-sig')
    print(f'내곽 기준 Z: {height:g}m / 경로 {len(routes)}개')
    print(summary.groupby('layer')[['contact_count', 'min_contact_gap_m', 'max_contact_gap_m']]
          .agg({'contact_count': 'sum', 'min_contact_gap_m': 'min', 'max_contact_gap_m': 'max'}).to_string())
    print(buffer_summary.groupby('layer').agg(overlap_routes=('buffer_overlap', 'sum'),
                                             min_3d_m=('min_distance_3d_m', 'min'),
                                             max_3d_m=('min_distance_3d_m', 'max')).to_string())
    print(f'CSV 5개, PNG 4개, 분석설명 저장 완료: {output.resolve()}')


if __name__ == '__main__':
    main()
