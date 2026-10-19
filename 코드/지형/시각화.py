"""최종정리 XYZ로 외곽 → 내곽 → FAF → OLS → 착륙을 그린다.

실행: python 코드/지형/시각화.py
필요 패키지: geopandas, shapely, numpy, matplotlib
고도는 원본 geometry Z(m)이며 AGL로 변환하지 않는다.
FAF_IAF_In.shp가 없으면 기존 이름 IAF_FAF.shp를 읽는다.
전체 altitude.png와 개별 altitude_route_01.png ~ altitude_route_16.png를 저장한다.
"""

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from shapely import get_coordinates
from shapely.geometry import LineString, Point


ROOT = Path(__file__).resolve().parents[2]
FT = 0.3048
OUTER = '#dc8b28'
INNER = '#2879ba'
INK = '#24384c'
TRANSFER = '#169b91'
LANDING = '#8952ae'


def descending(xyz):
    return xyz if xyz[0, 2] >= xyz[-1, 2] else xyz[::-1].copy()


def inner_arc(ring, start, end):
    """내곽 실측 선 위의 최단 원호. 운항 방향이 지정되지 않은 연결 예시."""
    line = LineString(ring[:, :2])
    a, b = (line.project(Point(x[:2])) for x in (start, end))
    forward = (b - a) % line.length
    delta = forward if forward <= line.length / 2 else forward - line.length
    stations = (a + np.linspace(0, delta, max(2, int(abs(delta) / 20) + 1))) % line.length
    xyz = np.array([[line.interpolate(s).x, line.interpolate(s).y, ring[0, 2]] for s in stations])
    # 수 mm인 원의 이산화 오차만 끝점에 흡수한다.
    xyz[0], xyz[-1] = start, end
    return xyz


def landing_surfaces(path, routes):
    frame = gpd.read_file(path).to_crs(5186)
    surfaces = []
    for geom in frame.geometry:
        if geom is None or geom.is_empty or geom.geom_type != 'Polygon':
            raise ValueError('최종OLS에는 유효한 PolygonZ가 필요합니다.')
        xyz = get_coordinates(geom.exterior, include_z=True)[:-1]
        if len(xyz) != 3 or not np.isfinite(xyz).all():
            raise ValueError('최종OLS의 삼각형 XYZ를 확인하세요.')
        touchdown = xyz[np.argmin(xyz[:, 2])]
        upper = xyz[~np.isclose(xyz[:, 2], touchdown[2])]
        if len(upper) != 2:
            raise ValueError('OLS 상단 모서리를 결정할 수 없습니다.')
        surfaces.append((xyz, np.vstack([upper.mean(axis=0), touchdown])))
    matched = {}
    for key, route in routes:
        distances = [np.linalg.norm(centerline[0] - route[-1]) for _, centerline in surfaces]
        index = int(np.argmin(distances))
        if distances[index] > .1:
            raise ValueError(f'경로 {key}: FAF와 OLS 상단 중앙의 불일치 {distances[index]:.3f}m')
        centerline = surfaces[index][1].copy()
        centerline[0] = route[-1]
        matched[key] = centerline
    return surfaces, matched


def load_xyz(path):
    frame = gpd.read_file(path)
    if frame.crs is None:
        raise ValueError(f'좌표계가 없습니다: {path}')
    frame = frame.to_crs(5186)
    result = []
    for index, row in frame.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty or geom.geom_type != 'LineString':
            raise ValueError(f'{path.name} 행 {index}: LineStringZ가 필요합니다.')
        xyz = get_coordinates(geom, include_z=True)
        if len(xyz) < 2 or not np.isfinite(xyz).all():
            raise ValueError(f'{path.name} 행 {index}: 유효한 XYZ가 필요합니다.')
        result.append((row.get('id', index + 1), xyz))
    if not result:
        raise ValueError(f'빈 레이어: {path}')
    return result


def configure_style():
    available = {font.name for font in font_manager.fontManager.ttflist}
    for name in ('Malgun Gothic', 'NanumGothic', 'Noto Sans CJK KR'):
        if name in available:
            plt.rcParams['font.family'] = name
            break
    plt.rcParams.update({'axes.unicode_minus': False, 'font.size': 10,
                         'text.color': INK, 'axes.labelcolor': INK,
                         'xtick.color': '#607080', 'ytick.color': '#607080',
                         'axes.edgecolor': '#ccd5df', 'savefig.facecolor': '#f7f9fc'})


def draw_individual(key, segments, rings, surfaces, local, output, dpi):
    """전체 그림과 동일한 경로 조합으로 개별 착륙 과정을 출력한다."""
    via_inner = len(segments) == 4
    mode = '1,500ft 내곽 경유' if via_inner else '외곽 직접 접근'
    landing = segments[-1][0]
    surface = min(surfaces, key=lambda item: np.linalg.norm(item[1][0] - landing[0]))[0]
    fig = plt.figure(figsize=(16, 11), facecolor='#f7f9fc')
    grid = fig.add_gridspec(2, 2, width_ratios=(1.35, 1), left=.04, right=.96,
                           top=.80, bottom=.20, wspace=.21, hspace=.40)
    ax3 = fig.add_subplot(grid[:, 0], projection='3d', facecolor='#f7f9fc')
    axmap = fig.add_subplot(grid[0, 1], facecolor='white')
    axprof = fig.add_subplot(grid[1, 1], facecolor='white')
    length = sum(np.linalg.norm(np.diff(xyz[:, :2], axis=0), axis=1).sum()
                 for xyz, _ in segments) / 1000
    fig.text(.055, .94, f'ROUTE {int(key):02d}  /  APPROACH → LANDING', fontsize=24, weight='bold')
    fig.text(.055, .896, f'{mode} → FAF → OLS 중심선 하강 → 착륙', fontsize=15)
    fig.text(.055, .855, f'수평 이동거리 {length:.2f} km    |    '
             f'출발 {segments[0][0][0, 2] / FT:,.0f} ft    →    '
             f'FAF {landing[0, 2] / FT:,.0f} ft    →    착륙 {landing[-1, 2] / FT:,.0f} ft',
             fontsize=11, color='#607080')
    for layer, color, label in rings:
        for _, xyz in layer:
            q = local(xyz)
            ax3.plot(*q.T, color=color, alpha=.20, lw=1.2)
            axmap.plot(q[:, 0], q[:, 1], color=color, alpha=.25, lw=1.2)
    polygon = local(surface)
    ax3.add_collection3d(Poly3DCollection([polygon], facecolors=LANDING,
                                        edgecolors=LANDING, alpha=.17, linewidths=.8))
    axmap.fill(polygon[:, 0], polygon[:, 1], color=LANDING, alpha=.15)
    at = 0.
    station_bounds = [0.]
    for xyz, color in segments:
        q = local(xyz)
        ax3.plot(*q.T, color=color, lw=2.6)
        axmap.plot(q[:, 0], q[:, 1], color=color, lw=2.2)
        distance = np.r_[0., np.cumsum(np.linalg.norm(np.diff(xyz[:, :2], axis=0), axis=1))] / 1000
        axprof.plot(at + distance, q[:, 2], color=color, lw=2.5)
        # 직선 OLS 구간도 포함해 평면도에 진행 방향을 표시한다.
        if distance[-1] > .01:
            a, b = distance[-1] * .45, distance[-1] * .62
            points = [np.array([np.interp(s, distance, q[:, i]) for i in (0, 1)]) for s in (a, b)]
            axmap.annotate('', xy=points[1], xytext=points[0],
                           arrowprops={'arrowstyle': '-|>', 'color': color, 'lw':1.8})
        at += distance[-1]
        station_bounds.append(at)

    marks = [(segments[0][0][0], 0., '출발', 'o', OUTER),
             (landing[0], station_bounds[-2], 'FAF', 'D', INK),
             (landing[-1], at, '착륙', 's', LANDING)]
    if via_inner:
        marks.insert(1, (segments[2][0][0], station_bounds[2], '내곽 IAF', 'o', INNER))
    for xyz, station, label, marker, color in marks:
        point = local(np.array([xyz]))[0]
        ax3.scatter(*point, marker=marker, color=color, s=45, edgecolors='white', depthshade=False)
        axmap.scatter(*point[:2], marker=marker, color=color, s=35, edgecolors='white', zorder=5)
        axmap.annotate(label, point[:2], xytext=(6, 6), textcoords='offset points', fontsize=9, color=color)
        axprof.scatter(station, point[2], marker=marker, color=color, s=35, zorder=5)
        axprof.annotate(label, (station, point[2]), xytext=(0, 10), textcoords='offset points',
                        ha='center', fontsize=9, color=color)
    ax3.set_title('01   개별 경로 · 3D', loc='left', fontsize=13, weight='bold')
    ax3.set_xlabel('동서 거리 (km)', labelpad=9)
    ax3.set_ylabel('남북 거리 (km)', labelpad=9)
    ax3.set_zlabel('원본 Z 고도 (ft)', labelpad=8)
    ax3.set_zlim(0, 2180)
    ax3.set_box_aspect((1, 1, .75))
    ax3.view_init(elev=26, azim=-58)
    for axis in (ax3.xaxis, ax3.yaxis, ax3.zaxis):
        axis.pane.fill = False
        axis._axinfo['grid']['color'] = (.82, .86, .9, .5)
    axmap.set_title('02   평면 경로 · 북쪽 ↑', loc='left', fontsize=13, weight='bold', pad=10)
    axmap.set_aspect('equal', adjustable='datalim')
    axmap.set_xlabel('동서 거리 (km)')
    axmap.set_ylabel('남북 거리 (km)')
    axprof.set_title('03   출발부터 착륙까지 고도 변화', loc='left', fontsize=13, weight='bold', pad=10)
    axprof.set_xlabel('외곽 출발점부터 누적 수평거리 (km)')
    axprof.set_ylabel('원본 Z 고도 (ft)')
    axprof.set_ylim(0, 2200)
    axprof.set_xlim(-.07 * length, length * 1.08)
    axprof.axhline(1500, color=INNER, ls='--', lw=.8, alpha=.4)
    for axis in (axmap, axprof):
        axis.spines[['top', 'right']].set_visible(False)
        axis.grid(color='#e6ebf1', lw=.65)
        axis.set_axisbelow(True)
    handles = [Line2D([0], [0], color=color, lw=2.5, label=label) for color, label in
               ((OUTER, '외곽 직접 접근'), (TRANSFER, '내곽 진입·경유'),
                (INNER, '내곽 IAF→FAF'), (LANDING, 'OLS→착륙'))]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .12), ncol=4, frameon=False)
    fig.text(.055, .09, '내곽 경유: 가장 가까운 진입부 + 실제 내곽의 최단 원호 연결 예시. 원본에 운항 방향은 지정되어 있지 않습니다.', fontsize=9, color='#607080')
    fig.text(.055, .063, 'OLS는 장애물 제한 표면입니다. FAF→OLS 저점 중심선은 착륙 시각화용 경로입니다.', fontsize=9, color='#607080')
    fig.text(.055, .036, '원본 geometry Z를 ft로 환산(AGL 아님) · EPSG:5186 · XY는 도심 기준 상대거리 · 3D 수직축 확대', fontsize=9, color='#607080')
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    print(f'Saved route {int(key):02d}: {output}')


def draw(input_dir, output, dpi=220, individual_only=False):
    configure_style()
    outer = load_xyz(input_dir / '외곽_2800.shp')
    inner = load_xyz(input_dir / '내곽_2000.shp')
    source = input_dir / 'FAF_IAF_In.shp'
    if not source.exists():
        source = input_dir / 'IAF_FAF.shp'
    raw_routes = load_xyz(source)
    # 낮은 끝점이 FAF. 저장된 선의 방향과 무관하게 하강 방향으로 정렬한다.
    routes = [(key, xyz if xyz[0, 2] >= xyz[-1, 2] else xyz[::-1].copy())
              for key, xyz in raw_routes]
    routes.sort(key=lambda item: item[0])
    transfers = [(key, descending(xyz)) for key, xyz in load_xyz(input_dir / 'Interchange_In.shp')]
    surfaces, landings = landing_surfaces(input_dir / '최종OLS.shp', routes)
    rings = [(outer, OUTER, '외곽'), (inner, INNER, '내곽')]
    heights = [float(np.median(np.concatenate([xyz[:, 2] for _, xyz in ring])))
               for ring in (outer, inner)]
    all_xyz = np.vstack([xyz for layer in (outer, inner) for _, xyz in layer])
    origin = (all_xyz[:, :2].min(axis=0) + all_xyz[:, :2].max(axis=0)) / 2

    def local(xyz):
        return np.column_stack(((xyz[:, :2] - origin) / 1000, xyz[:, 2] / FT))

    def route_color(xyz):
        return OUTER if abs(xyz[0, 2] - heights[0]) < abs(xyz[0, 2] - heights[1]) else INNER

    # 내곽 출발 경로는 가장 가까운 진입부에서 실제 내곽 선을 따라 연결한다.
    # 원호의 운항 방향/진입부 선택은 원본에 없어 최단거리 연결 예시로 명시한다.
    complete = []
    for key, xyz in routes:
        segments = []
        if abs(xyz[0, 2] - heights[1]) < .1:
            options = [(transfer, inner_arc(inner[0][1], transfer[-1], xyz[0]))
                       for _, transfer in transfers]
            transfer, arc = min(options, key=lambda pair: np.linalg.norm(np.diff(pair[1][:, :2], axis=0), axis=1).sum())
            segments += [(transfer, TRANSFER), (arc, TRANSFER)]
        segments += [(xyz, route_color(xyz)), (landings[key], LANDING)]
        complete.append((key, segments))

    fig = plt.figure(figsize=(18, 12), facecolor='#f7f9fc')
    grid = fig.add_gridspec(2, 2, width_ratios=(1.55, 1), height_ratios=(1, 1),
                           left=.035, right=.96, bottom=.16, top=.80,
                           wspace=.18, hspace=.36)
    ax3 = fig.add_subplot(grid[:, 0], projection='3d', facecolor='#f7f9fc')
    axmap = fig.add_subplot(grid[0, 1], facecolor='white')
    axprof = fig.add_subplot(grid[1, 1], facecolor='white')
    fig.text(.055, .94, 'APPROACH  /  IAF → FAF → OLS → LANDING', fontsize=23, weight='bold')
    fig.text(.055, .897, '외곽 직접 접근 + 1,500ft 내곽 경유 접근 → OLS 중심선 하강 → 착륙', fontsize=15)
    counts = [sum(abs(xyz[0, 2] - h) < .1 for _, xyz in routes) for h in heights]
    fig.text(.055, .858,
             f'총 {len(routes)}개 경로    |    외곽 {heights[0] / FT:,.0f} ft · {counts[0]}개'
             f'    |    내곽 {heights[1] / FT:,.0f} ft · {counts[1]}개', fontsize=11, color='#607080')

    for (layer, color, label), h in zip(rings, heights):
        for _, xyz in layer:
            q = local(xyz)
            ax3.plot(*q.T, color=color, lw=2.7, alpha=.9)
            ax3.plot(q[:, 0], q[:, 1], np.zeros(len(q)), color=color, lw=.8, alpha=.13)
            axmap.plot(q[:, 0], q[:, 1], color=color, lw=2, alpha=.6)
        axprof.axhline(h / FT, color=color, ls='--', lw=1, alpha=.65)

    # OLS는 장애물 제한 표면이다. 상단 중앙→저점 중심선은 착륙 시각화용이다.
    ax3.add_collection3d(Poly3DCollection([local(xyz) for xyz, _ in surfaces],
                         facecolors=LANDING, edgecolors=LANDING, alpha=.10, linewidths=.5))
    for xyz, _ in surfaces:
        q = local(xyz)
        axmap.fill(q[:, 0], q[:, 1], color=LANDING, alpha=.07)
    for _, xyz in transfers:
        q = local(xyz)
        ax3.plot(*q.T, color=TRANSFER, lw=1, alpha=.45)
        axmap.plot(q[:, 0], q[:, 1], color=TRANSFER, lw=.8, alpha=.45)
    for _, segments in complete:
        at = 0.
        for xyz, color in segments:
            q = local(xyz)
            distance = np.r_[0, np.cumsum(np.linalg.norm(np.diff(xyz[:, :2], axis=0), axis=1))] / 1000
            axprof.plot(at + distance, q[:, 2], color=color, lw=1.45, alpha=.8)
            at += distance[-1]
            if color in (TRANSFER, LANDING):
                ax3.plot(*q.T, color=color, lw=1.7, alpha=.8)
                axmap.plot(q[:, 0], q[:, 1], color=color, lw=1.2, alpha=.8)
        axprof.scatter(at, q[-1, 2], marker='s', s=18, color=LANDING, zorder=4)

    touchdown = np.unique(np.round([centerline[-1] for centerline in landings.values()], 5), axis=0)
    td = local(touchdown)
    ax3.scatter(*td.T, marker='s', color=LANDING, s=36, edgecolors='white', linewidths=.6, depthshade=False)
    axmap.scatter(td[:, 0], td[:, 1], marker='s', color=LANDING, s=28, edgecolors='white', linewidths=.5, zorder=6)

    endpoints = []
    for key, xyz in routes:
        q = local(xyz)
        color = route_color(xyz)
        ax3.plot(*q.T, color=color, lw=1.7, alpha=.9)
        ax3.scatter(*q[0], color=color, s=20, edgecolors='white', linewidths=.5)
        axmap.plot(q[:, 0], q[:, 1], color=color, lw=1.25, alpha=.8)
        axmap.scatter(*q[0, :2], color=color, s=12, zorder=4)
        # 각 경로의 중간에 실제 하강 방향 화살표.
        k = min(max(len(q) // 2, 1), len(q) - 1)
        j = max(0, k - max(1, len(q) // 12))
        axmap.annotate('', xy=q[k, :2], xytext=q[j, :2],
                       arrowprops={'arrowstyle': '-|>', 'color': color, 'lw':1.2})
        endpoints.append(q[-1])

    # 같은 FAF를 공유하는 접근 경로의 점을 중복 표시하지 않는다.
    endpoints = np.unique(np.round(endpoints, 6), axis=0)
    ax3.scatter(*endpoints.T, marker='D', color=INK, s=32,
                edgecolors='white', linewidths=.7, depthshade=False)
    axmap.scatter(endpoints[:, 0], endpoints[:, 1], marker='D', color=INK,
                  s=25, edgecolors='white', linewidths=.6, zorder=5)
    ax3.set_title('01   3D 접근 경로', loc='left', fontsize=13, weight='bold', pad=6)
    ax3.set_xlabel('동서 거리 (km)', labelpad=10)
    ax3.set_ylabel('남북 거리 (km)', labelpad=10)
    ax3.set_zlabel('원본 Z 고도 (ft)', labelpad=9)
    ax3.set_zlim(0, max(heights) / FT + 180)
    ax3.set_box_aspect((1, 1, .68))
    ax3.view_init(elev=26, azim=-58)
    ax3.tick_params(labelsize=9)
    for axis in (ax3.xaxis, ax3.yaxis, ax3.zaxis):
        axis.pane.fill = False
        axis._axinfo['grid']['color'] = (.82, .86, .9, .5)
    axmap.set_title('02   평면 경로 · 북쪽 ↑', loc='left', fontsize=13, weight='bold', pad=12)
    axmap.set_aspect('equal', adjustable='datalim')
    axmap.set_xlabel('동서 거리 (km)')
    axmap.set_ylabel('남북 거리 (km)')
    axprof.set_title('03   내곽 경유·직접 접근부터 착륙까지', loc='left', fontsize=13, weight='bold', pad=12)
    axprof.set_xlabel('외곽 출발점부터 누적 수평거리 (km)')
    axprof.set_ylabel('원본 Z 고도 (ft)')
    axprof.set_ylim(0, max(heights) / FT + 130)
    axprof.set_xlim(left=0)
    for axis in (axmap, axprof):
        axis.spines[['top', 'right']].set_visible(False)
        axis.grid(color='#e6ebf1', lw=.65)
        axis.set_axisbelow(True)

    handles = [Line2D([0], [0], color=c, lw=2.5, label=f'{label} / 해당 IAF 출발 경로')
               for _, c, label in rings]
    handles += [Line2D([0], [0], color=TRANSFER, lw=2.5, label='외곽→내곽 진입 / 내곽 경유'),
                Line2D([0], [0], color=LANDING, lw=2.5, label='OLS / 착륙 하강 중심선'),
                Line2D([0], [0], marker='D', color='none', markerfacecolor=INK,
                       label='FAF · OLS 진입'),
                Line2D([0], [0], marker='s', color='none', markerfacecolor=LANDING,
                       label='착륙 지점')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .085),
               ncol=3, frameon=False, fontsize=10)
    fig.text(.055, .062,
             '내곽 경유는 가장 가까운 Interchange_In 진입부와 실제 내곽의 최단 원호를 잇는 예시입니다. 원본에 운항 방향은 지정되어 있지 않습니다.',
             fontsize=9, color='#607080')
    fig.text(.055, .041,
             'OLS는 장애물 제한 표면이며, 보라색 중심선은 FAF에서 OLS 상단 중앙을 거쳐 저점(착륙 지점)으로 내려가는 시각화용 경로입니다.',
             fontsize=9, color='#607080')
    fig.text(.055, .017,
             '고도는 geometry Z(m)를 ft로 환산한 값이며 AGL이 아닙니다. '
             '실제 Z 사용 · EPSG:5186 · XY는 도심 기준 상대거리 · 3D 수직축 확대 · 자료: 최종정리 / IAF_FAF, Interchange_In, 최종OLS',
             fontsize=9, color='#607080')
    output.parent.mkdir(parents=True, exist_ok=True)
    if not individual_only:
        fig.savefig(output, dpi=dpi)
    plt.close(fig)
    if not individual_only:
        print(f'Saved: {output}')
    print(f'Routes: {len(routes)}, outer: {counts[0]}, inner: {counts[1]}, FAF endpoints: {len(endpoints)}')
    print(f'Interchanges: {len(transfers)}, OLS surfaces: {len(surfaces)}, matched landing paths: {len(landings)}, touchdown sites: {len(touchdown)}')
    for key, segments in complete:
        individual = output.with_name(f'{output.stem}_route_{int(key):02d}.png')
        draw_individual(key, segments, rings, surfaces, local, individual, dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=ROOT / '자료/최종정리')
    parser.add_argument('--output', type=Path, default=ROOT / '자료/결과물/altitude.png')
    parser.add_argument('--dpi', type=int, default=220)
    parser.add_argument('--individual-only', action='store_true', help='전체 그림을 덮어쓰지 않고 개별 경로만 생성')
    args = parser.parse_args()
    draw(args.input_dir, args.output, args.dpi, args.individual_only)


if __name__ == '__main__':
    main()
