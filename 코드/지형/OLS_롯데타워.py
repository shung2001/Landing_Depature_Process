"""Polygon 26과 롯데타워의 XY·종단·3D 이격 분석.

실행: python 코드/지형/OLS_롯데타워.py
원본 A16=555가 없으면 update 파일의 동일 형상을 확인한 뒤 555m 가정임을 명시한다.
필요 패키지: geopandas, shapely, numpy, pandas, rasterio, matplotlib.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from shapely import force_2d
from shapely.geometry import Polygon, LineString
from shapely.ops import nearest_points
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[2]
POLYGONS = ROOT / '자료/결과물/Polygon/8/안전지대/26/Polygon_26.shp'
BUILDINGS = ROOT / '자료/안전성/송파구/송파구_건물_수정.shp'
UPDATED = BUILDINGS.with_name('송파구_건물_수정_update.shp')
DEM = BUILDINGS.with_name('송파구_DEM.tif')
OUT = BUILDINGS.parent / 'Separation_확보/롯데월드_충돌'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'Malgun Gothic', 'axes.unicode_minus': False})
    buildings = gpd.read_file(BUILDINGS)
    targets = buildings.loc[pd.to_numeric(buildings.A16, errors='coerce').eq(555)].copy()
    assumption = targets.empty
    if assumption:
        updated = gpd.read_file(UPDATED).to_crs(buildings.crs)
        candidates = updated.loc[pd.to_numeric(updated.A16, errors='coerce').eq(555)]
        if candidates.empty:
            raise ValueError('원본과 update 모두 A16=555 대상이 없습니다.')
        indices = []
        for _, row in candidates.iterrows():
            matches = buildings.index[buildings.geometry.map(lambda x: x.equals(row.geometry))]
            if len(matches) != 1:
                raise ValueError('update 대상과 원본 형상의 유일한 일치를 확인할 수 없습니다.')
            indices.append(matches[0])
        targets = buildings.loc[indices].copy()
    targets['source_row'] = targets.index
    targets = targets.to_crs(5186)
    polygons = gpd.read_file(POLYGONS).to_crs(5186)
    if not polygons.z_ref.eq('DEM+AGL').all():
        raise ValueError('Polygon 높이 기준을 확인해야 합니다.')
    rows, links, surfaces = [], [], []
    with PdfPages(OUT / '전체_Polygon_분석도.pdf') as pdf, rasterio.open(DEM) as dem:
        for _, building in targets.iterrows():
            footprint = force_2d(building.geometry)
            if footprint.geom_type != 'Polygon' or not footprint.is_valid:
                raise ValueError('유효한 단일 건물 Polygon이 필요합니다.')
            point = gpd.GeoSeries([footprint.representative_point()], crs=5186).to_crs(dem.crs).iloc[0]
            sample = next(dem.sample([(point.x, point.y)], masked=True))[0]
            if np.ma.is_masked(sample) or not np.isfinite(sample):
                raise ValueError('건물 DEM 표본을 읽을 수 없습니다.')
            ground, roof = float(sample), float(sample) + 555
            origin = np.array(footprint.centroid.coords[0])
            ring = np.asarray(footprint.exterior.coords) - origin
            faces = [np.column_stack([ring[:-1], np.full(len(ring)-1, roof)])]
            for a, b in zip(ring[:-1], ring[1:]):
                faces.append(np.array([[*a, ground], [*b, ground], [*b, roof], [*a, roof]]))
            local_rows = []
            overview, ax_over = plt.subplots(figsize=(11, 10))
            ax_over.fill(*ring.T, color='#de7047', label='건물 footprint / 555m 가정' if assumption else '건물 footprint')
            scene = plt.figure(figsize=(12, 10))
            ax_scene = scene.add_subplot(111, projection='3d')
            ax_scene.add_collection3d(Poly3DCollection(faces, facecolor='#dc794a', alpha=.55, edgecolor='#8d4028', linewidth=.4))
            for source_row, feature in polygons.iterrows():
                xyz = np.asarray(feature.geometry.exterior.coords)[:-1]
                if xyz.shape != (3, 3) or not np.isfinite(xyz).all():
                    raise ValueError('유한한 3개 꼭짓점의 Polygon Z가 필요합니다.')
                xy = xyz[:, :2]
                surface = Polygon(xy)
                # 원본 Polygon Z 세 꼭짓점을 지나는 평면. 투영 XY는 미터 단위.
                coeff = np.linalg.solve(np.column_stack([xy - origin, np.ones(3)]), xyz[:, 2])
                height = lambda p: float(np.dot(np.r_[np.asarray(p)-origin, 1], coeff))
                assert np.allclose([height(p) for p in xy], xyz[:, 2], atol=1e-7)
                bp, sp = nearest_points(footprint, surface)
                z = height(sp.coords[0])
                horizontal = footprint.distance(surface)
                # 이 데이터에서는 표면 전체가 건물 바닥과 옥상 사이에 존재한다.
                # 따라서 XY 최단점과 동일 Z의 외벽점을 연결하면 3D 하한(XY 거리)을 달성한다.
                # 다른 높이 범위의 데이터에 잘못 적용하지 않도록 명시적으로 검증한다.
                if not (ground <= xyz[:, 2].min() and xyz[:, 2].max() <= roof):
                    raise ValueError('표면 높이가 건물 수직 구간 밖입니다. 일반 3D 최적화가 필요합니다.')
                connector = LineString([(bp.x, bp.y, z), (sp.x, sp.y, z)])
                assert np.isclose(np.linalg.norm(np.array(connector.coords[0])-connector.coords[1]), horizontal)
                apex_i = int(np.argmin(xyz[:, 2]))
                apex = xyz[apex_i]
                base = np.mean(np.delete(xyz, apex_i, axis=0), axis=0)
                axis = base[:2]-apex[:2]
                length = np.linalg.norm(axis)
                unit = axis/length
                centerline = LineString([apex[:2], base[:2]])
                projected_building = (np.asarray(footprint.exterior.coords)-apex[:2]) @ unit
                cross = (np.asarray(footprint.exterior.coords)-apex[:2]) @ np.array([-unit[1], unit[0]])
                overlap = footprint.intersection(surface)
                bearing = float(feature.rot_deg)
                record = dict(polygon_row=int(source_row), building_row=int(building.source_row), bearing_deg=bearing,
                    height_assumed=assumption, building_agl_m=555, building_ground_m=ground, building_roof_m=roof,
                    polygon_min_z_m=float(xyz[:, 2].min()), polygon_max_z_m=float(xyz[:, 2].max()),
                    xy_gap_m=horizontal, distance_3d_m=horizontal, centerline_xy_gap_m=footprint.distance(centerline),
                    xy_overlap_m2=overlap.area, geometric_contact=horizontal <= 1e-7,
                    surface_z_at_xy_nearest_m=z, roof_minus_surface_at_xy_nearest_m=roof-z,
                    vertical_interval_gap_m=0.0,
                    building_x=bp.x, building_y=bp.y, surface_x=sp.x, surface_y=sp.y,
                    nearest_z_m=z, building_cross_min_m=float(cross.min()), building_cross_max_m=float(cross.max()))
                rows.append(record); local_rows.append(record)
                links.append({**record, 'geometry': connector})
                surfaces.append({**record, 'geometry': feature.geometry})
                color = '#267eac' if bearing < 180 else '#8065af'
                local = xyz.copy(); local[:, :2] -= origin
                ax_over.plot(*np.vstack([local[:, :2], local[0, :2]]).T, color=color, alpha=.6, linewidth=.8)
                ax_over.text(*local[1:, :2].mean(axis=0), f'{bearing:g}°', fontsize=7, color=color)
                ax_scene.add_collection3d(Poly3DCollection([local], facecolor=color, alpha=.09, edgecolor=color, linewidth=.4))
                fig = plt.figure(figsize=(18, 6))
                ax = fig.add_subplot(131)
                ax.fill(*local[:, :2].T, color=color, alpha=.22, label='Polygon 표면의 XY 투영')
                ax.fill(*ring.T, color='#de7047', label='건물 footprint')
                line = np.array([bp.coords[0], sp.coords[0]])-origin
                ax.plot(*line.T, 'r.-', label=f'최단 간격 {horizontal:.2f}m')
                ax.plot(*np.array([apex[:2]-origin, base[:2]-origin]).T, '--', color=color, label='기하학적 중심선')
                ax.set(title=f'XY 평면 / 방위각 {bearing:g}°', xlabel='건물 중심 기준 동쪽 (m)', ylabel='북쪽 (m)', aspect='equal')
                ax.legend(fontsize=8); ax.grid(alpha=.2)
                ax = fig.add_subplot(132)
                ax.fill_between([projected_building.min(), projected_building.max()], ground, roof, color='#de7047', alpha=.5, label='건물의 종단면 투영')
                ax.plot([0, length], [apex[2], base[2]], color=color, linewidth=3, label='Polygon 중심선 높이')
                along = (np.array(sp.coords[0])-apex[:2]) @ unit
                ax.scatter([along], [z], color='red', label='XY 최근접 표면점')
                ax.axhline(roof, color='#de7047', linestyle=':', label=f'옥상 {roof:.2f}m')
                ax.set(title=f'종단 높이 투영 / 옥상 - 최근접 표면 {roof-z:.2f}m\n좌우 횡방향을 접은 투영이므로 겹침 ≠ 충돌', xlabel='시작점에서 중심선 방향 투영거리 (m)', ylabel='DEM + AGL 높이 (m)', ylim=(0, roof+40))
                ax.legend(fontsize=8); ax.grid(alpha=.2)
                ax = fig.add_subplot(133, projection='3d')
                ax.add_collection3d(Poly3DCollection(faces, facecolor='#de7047', alpha=.45, edgecolor='#8d4028', linewidth=.3))
                ax.add_collection3d(Poly3DCollection([local], facecolor=color, alpha=.4, edgecolor=color))
                ax.plot(line[:, 0], line[:, 1], [z, z], color='red', linewidth=2)
                bounds = np.vstack([local[:, :2], ring])
                ax.set(xlim=(bounds[:, 0].min()-30, bounds[:, 0].max()+30), ylim=(bounds[:, 1].min()-30, bounds[:, 1].max()+30), zlim=(0, roof+40), xlabel='동쪽 (m)', ylabel='북쪽 (m)', zlabel='DEM+AGL (m)', title=f'3D 표면↔건물 입체 최단거리 {horizontal:.2f}m')
                ax.set_box_aspect([np.ptp(bounds[:, 0])+60, np.ptp(bounds[:, 1])+60, roof+40])
                ax.view_init(elev=25, azim=-55)
                fig.suptitle(f'Polygon 26 / 원본 Polygon 행 {source_row} / 건물 행 {int(building.source_row)}'+(' / A16 결측: 555m 가정' if assumption else ''), fontsize=13)
                fig.tight_layout()
                fig.savefig(OUT / f'방위각_{bearing:03.0f}_건물{int(building.source_row)}_XY_높이_3D.png', dpi=160)
                pdf.savefig(fig); plt.close(fig)
            ax_over.set(aspect='equal', xlabel='건물 중심 기준 동쪽 (m)', ylabel='북쪽 (m)', title='20개 Polygon과 대상 건물의 XY 위치\n파랑: 88–97° / 보라: 218–227°')
            ax_over.legend(); ax_over.grid(alpha=.2); overview.tight_layout()
            overview.savefig(OUT / '전체_XY_배치.png', dpi=180); pdf.savefig(overview); plt.close(overview)
            coords = np.concatenate([np.asarray(g.exterior.coords)[:, :2]-origin for g in polygons.geometry]+[ring])
            ax_scene.set(xlim=(coords[:, 0].min()-30, coords[:, 0].max()+30), ylim=(coords[:, 1].min()-30, coords[:, 1].max()+30), zlim=(0, roof+40), xlabel='동쪽 (m)', ylabel='북쪽 (m)', zlabel='DEM+AGL (m)', title='전체 Polygon 3D / 건물은 footprint × 555m 수직기둥 모델')
            ax_scene.set_box_aspect([np.ptp(coords[:, 0])+60, np.ptp(coords[:, 1])+60, roof+40]); ax_scene.view_init(25, -55)
            scene.tight_layout(); scene.savefig(OUT / '전체_3D.png', dpi=180); pdf.savefig(scene); plt.close(scene)
    table = pd.DataFrame(rows)
    table.to_csv(OUT / 'Polygon별_간격.csv', index=False, encoding='utf-8-sig', float_format='%.6f')
    gpkg = OUT / '롯데월드_간격분석.gpkg'
    gpd.GeoDataFrame(surfaces, crs=5186).to_file(gpkg, layer='polygon_surfaces', driver='GPKG')
    gpd.GeoDataFrame(links, crs=5186).to_file(gpkg, layer='nearest_3d_connectors', driver='GPKG')
    targets.to_file(gpkg, layer='target_buildings', driver='GPKG')
    fig, ax = plt.subplots(figsize=(13, 6))
    labels = [f'{v:g}°' for v in table.bearing_deg]
    bars = ax.bar(labels, table.xy_gap_m, color=['#267eac' if v<180 else '#8065af' for v in table.bearing_deg])
    ax.bar_label(bars, fmt='%.2f', fontsize=8, padding=3)
    ax.set(ylabel='최단 간격 (m)', xlabel='Polygon 방위각', title='각 Polygon과 건물의 간격 / 이 데이터에서는 XY 거리 = 3D 입체 최단거리')
    ax.tick_params(axis='x', rotation=45); ax.grid(axis='y', alpha=.2); fig.tight_layout()
    fig.savefig(OUT / 'Polygon별_간격비교.png', dpi=180); plt.close(fig)
    minimum = table.loc[table.distance_3d_m.idxmin()]
    notes = f'''# 롯데월드 Polygon 26 간격 분석

## 입력과 높이 기준
- 요청 원본: {BUILDINGS.relative_to(ROOT)}
- Polygon: {POLYGONS.relative_to(ROOT)}
- 원본 A16=555 피처 수: {0 if assumption else len(targets)}.
- {'원본 A16이 결측이다. update의 A16=555 피처와 원본의 형상 완전 일치를 확인하고 원본 형상에 555m를 적용한 가정 분석이다. 원본 파일은 수정하지 않았다.' if assumption else '원본 A16=555 피처를 사용했다.'}
- source_row/building_row/polygon_row는 0부터 시작하는 원본 행 번호.
- EPSG:5186 미터 좌표계. Polygon Z의 DEM+AGL 높이를 보존했다.
- 건물 바닥은 Polygon 생성 코드와 동일하게 footprint 내부 대표점에서 송파구_DEM.tif의 셀 값을 추출했다. 건물은 그 높이부터 555m까지의 수직기둥이다. 실제 타워의 상층부 축소 형상은 포함하지 않는다.
- 수직 datum은 입력 DEM 기준이다. AGL 555m를 Polygon 절대 높이와 직접 빼지 않았다.

## 결과
- Polygon {len(polygons)}개 × 건물 {len(targets)}개 = {len(table)}개 비교.
- 최소 간격: {minimum.distance_3d_m:.3f}m, 방위각 {minimum.bearing_deg:g}°.
- 최대 간격: {table.distance_3d_m.max():.3f}m.
- XY 접촉/겹침: {int(table.geometric_contact.sum())}개.
- 건물 옥상 DEM+AGL: {table.building_roof_m.min():.3f}m.
- Polygon 높이: {table.polygon_min_z_m.min():.3f}~{table.polygon_max_z_m.max():.3f}m.

## 지표 해석
- xy_gap_m: 채워진 Polygon의 XY 투영과 건물 footprint의 최단거리. 꼭짓점끼리의 거리나 중심점 거리가 아니다.
- distance_3d_m: 경사진 삼각 표면과 건물 입체의 최단거리. 이번 입력에서는 삼각 표면 전체 높이가 건물 바닥~옥상 사이에 있으므로 XY 최단점과 같은 높이의 외벽점이 존재한다. 따라서 3D 거리의 하한인 XY 거리를 실제로 달성하며 두 값이 같다. 실행 시 이 조건을 모든 Polygon에 검증한다.
- centerline_xy_gap_m: 시작 꼭짓점~반대쪽 변 중점의 기하학적 중심선과 건물 사이 거리. 실제 비행 궤적이 제공된 것은 아니다.
- roof_minus_surface_at_xy_nearest_m: XY 최단 표면점에서 옥상까지 높이 차. 옥상이 더 높으면 양수. 수평으로 떨어진 점의 높이 비교이며 수직 안전 여유나 3D 최단거리가 아니다.
- vertical_interval_gap_m: 건물과 Polygon의 높이 구간 간 거리. 높이 구간이 겹치므로 0이다. XY가 떨어져 있으면 이것만으로 충돌하지 않는다.
- 종단 그림은 건물을 중심선 방향에 투영한다. 횡방향 거리를 생략하므로 그림에서 겹쳐 보이는 것은 실제 교차를 뜻하지 않는다. CSV의 building_cross_min/max_m에 횡방향 범위를 기록했다.
- Polygon은 두께 없는 표면, 건물은 수직기둥으로 모델링했다. 간격은 기하학적 결과이며 항공기 크기, 운항 오차, 요구 Separation 기준은 적용하지 않았다.

## 생성 자료
- Polygon별_간격.csv: 모든 비교, 최단점 좌표, 높이, 중심선 거리.
- 방위각_*_XY_높이_3D.png: Polygon별 XY/종단/3D 3패널.
- 전체_XY_배치.png, 전체_3D.png, Polygon별_간격비교.png.
- 전체_Polygon_분석도.pdf: 전체 및 개별 도면 모음.
- 롯데월드_간격분석.gpkg: EPSG:5186 Polygon Z, 3D 최근접 연결선, 원본 대상 건물 레이어.
- 코드: 코드/지형/OLS_롯데타워.py. 재실행 시 같은 이름의 생성 결과를 갱신한다.
'''
    (OUT / '분석보고서.md').write_text(notes, encoding='utf-8')
    hashes = {}
    for source in [POLYGONS, BUILDINGS, UPDATED, DEM]:
        parts = [source] if source.suffix == '.tif' else [source.with_suffix(ext) for ext in ['.shp', '.shx', '.dbf', '.prj']]
        for part in parts:
            hashes[str(part.relative_to(ROOT))] = hashlib.sha256(part.read_bytes()).hexdigest()
    (OUT / '입력_SHA256.json').write_text(json.dumps(hashes, ensure_ascii=False, indent=2), encoding='utf-8')
    print(table[['bearing_deg', 'xy_gap_m', 'distance_3d_m', 'roof_minus_surface_at_xy_nearest_m']].to_string(index=False))
    print(f'완료: {OUT}')


if __name__ == '__main__':
    main()
