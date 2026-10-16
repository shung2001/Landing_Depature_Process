"""최종 OLS의 끝단 선분 중심을 FAF MultiPoint Z로 저장한다.

실행: python -B 코드/polygon/FAF.py
출력: 자료/안전성/송파구/Separation_확보/최종_OLS/최종/FAF.shp

Polygon당 점 하나가 들어 있는 MultiPoint Z 피처 하나를 만든다.
id는 원본 id → fid → 1부터 시작하는 행 번호 순서로 선택한다.
x_5186/y_5186은 EPSG:5186 미터 좌표, alt는 원본 Polygon Z 기준(m)이다.
현재 입력의 Z는 DEM+AGL이며 alt를 현지 AGL로 임의 변환하지 않는다.
"""
from pathlib import Path
import argparse

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
from shapely.geometry import MultiPoint


ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / '자료/안전성/송파구/Separation_확보/최종_OLS/최종'


def terminal_midpoint(row):
    """시작 꼭짓점의 반대편 선분 양 끝의 X/Y/Z를 각각 평균한다."""
    polygon = row.geometry
    if polygon is None or polygon.is_empty or polygon.geom_type != 'Polygon':
        raise ValueError('비어 있지 않은 단일 Polygon이 필요합니다.')
    if not polygon.is_valid or len(polygon.interiors):
        raise ValueError('구멍 없는 유효한 삼각형 Polygon이 필요합니다.')
    xyz = np.asarray(polygon.exterior.coords, dtype=float)[:-1]
    if xyz.shape != (3, 3) or not np.isfinite(xyz).all():
        raise ValueError('유한한 XYZ 꼭짓점 세 개의 Polygon Z가 필요합니다.')
    if len(np.unique(xyz[:, :2], axis=0)) != 3:
        raise ValueError('삼각형 꼭짓점이 중복됩니다.')
    if 'center_x' in row.index and 'center_y' in row.index:
        center = np.array([row.center_x, row.center_y], dtype=float)
        if not np.isfinite(center).all():
            raise ValueError('시작 중심점 좌표가 유효하지 않습니다.')
        distances = np.linalg.norm(xyz[:, :2] - center, axis=1)
        apex = int(np.argmin(distances))
        if distances[apex] > 0.001:
            raise ValueError('center_x/y가 Polygon 시작 꼭짓점과 일치하지 않습니다.')
    else:
        minimum = np.min(xyz[:, 2])
        candidates = np.flatnonzero(np.isclose(xyz[:, 2], minimum, rtol=0, atol=1e-7))
        if len(candidates) != 1:
            raise ValueError('시작 좌표가 없고 최저 Z 꼭짓점이 유일하지 않아 끝단을 정할 수 없습니다.')
        apex = int(candidates[0])
    endpoints = np.delete(xyz, apex, axis=0)
    return endpoints.mean(axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DIRECTORY / '최종OLS.shp')
    parser.add_argument('--output', type=Path, default=DIRECTORY / 'FAF.shp')
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError('원본과 출력 경로가 같을 수 없습니다.')
    source = gpd.read_file(args.input)
    if source.empty or source.crs is None:
        raise ValueError('입력에 피처와 좌표계가 필요합니다.')
    source = source.to_crs(5186)
    id_field = next((field for field in ('id', 'fid') if field in source.columns), None)
    identifiers = source[id_field].copy() if id_field else pd.Series(np.arange(1, len(source) + 1), index=source.index)
    if identifiers.isna().any() or identifiers.duplicated().any():
        raise ValueError('원본 id/fid는 결측 없이 Polygon별로 고유해야 합니다.')
    # Shapefile의 실수형 fid가 정수라면 불필요한 .0 없이 정수로 저장한다.
    if pd.api.types.is_numeric_dtype(identifiers):
        values = identifiers.to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError('id가 유효하지 않습니다.')
        if np.equal(values, np.floor(values)).all():
            identifiers = identifiers.astype('int64')
    records = []
    for index, row in source.iterrows():
        try:
            x, y, z = terminal_midpoint(row)
        except ValueError as error:
            raise ValueError(f'입력 행 {index}, id={identifiers.loc[index]}: {error}') from error
        records.append(dict(id=identifiers.loc[index], x_5186=float(x), y_5186=float(y),
                            alt=float(z), geometry=MultiPoint([(x, y, z)])))
    result = gpd.GeoDataFrame(records, crs=5186)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pyogrio.write_dataframe(result, args.output, driver='ESRI Shapefile',
                            geometry_type='MultiPoint Z', encoding='UTF-8')
    saved = gpd.read_file(args.output)
    if len(saved) != len(source) or saved.crs.to_epsg() != 5186:
        raise RuntimeError('저장 후 피처 수 또는 좌표계가 다릅니다.')
    if list(saved.columns) != ['id', 'x_5186', 'y_5186', 'alt', 'geometry']:
        raise RuntimeError('저장 후 필드 구성이 다릅니다.')
    if not saved.geom_type.eq('MultiPoint').all() or not saved.has_z.all():
        raise RuntimeError('MultiPoint Z로 저장되지 않았습니다.')
    if saved.id.tolist() != result.id.tolist():
        raise RuntimeError('저장 후 id가 다릅니다.')
    for i, row in saved.iterrows():
        if len(row.geometry.geoms) != 1:
            raise RuntimeError('Polygon당 점 하나의 MultiPoint가 필요합니다.')
        point = row.geometry.geoms[0]
        if not np.allclose([point.x, point.y, point.z], [row.x_5186, row.y_5186, row.alt], rtol=0, atol=1e-7):
            raise RuntimeError('기하 좌표와 속성이 다릅니다.')
        if not np.allclose(point.coords[0], result.geometry.iloc[i].geoms[0].coords[0], rtol=0, atol=1e-7):
            raise RuntimeError('저장 후 중심점 좌표가 다릅니다.')
    print(f'생성 완료: {args.output}')
    print(f'MultiPoint Z {len(saved)}개 / EPSG:5186 / id 원본 필드: {id_field or "행 번호"}')
    print('alt: 원본 Polygon Z 기준(m), 현재 입력은 DEM+AGL')
    print(saved.drop(columns='geometry').to_string(index=False))


if __name__ == '__main__':
    main()
