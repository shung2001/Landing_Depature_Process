"""서로 다른 폴리곤 SHP가 겹치는 영역과 영역별 중복 파일 수를 계산한다.

각 SHP 내부의 모든 피처는 먼저 하나의 면 형상으로 dissolve한다. 이후에는
파일당 하나의 형상만 사용하고, 서로 다른 파일 쌍만 비교한다. 따라서 한 SHP
내부 피처끼리의 중복이나 동일 피처의 자기 자신 비교는 중복 수에 포함되지 않는다.

기본 입력:
    자료/결과물/Polygon/12.5/<하위 폴더>/*.shp

기본 출력:
    자료/결과물/Polygon/12.5/중복/중복.shp
    자료/결과물/Polygon/12.5/중복/중복.gpkg

Shapefile의 DBF 필드명은 최대 10바이트이므로 ``중복된_개수``를 온전히
저장할 수 없다. SHP에는 ``중복개수``를 사용하고, 함께 만드는 GeoPackage에는
요청한 원래 필드명 ``중복된_개수``를 사용한다.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import shapely
from pyproj import CRS
from shapely import Geometry, STRtree


EXACT_COUNT_FIELD = "중복된_개수"
SHP_COUNT_FIELD = "중복개수"
SHP_ENCODING = "CP949"


def project_root() -> Path:
    """이 스크립트의 위치를 기준으로 저장소 루트를 반환한다."""

    return Path(__file__).resolve().parents[2]


def default_input_dir() -> Path:
    return project_root() / "자료" / "결과물" / "Polygon" / "12.5"


def natural_path_key(path: Path) -> tuple[tuple[int, int | str], ...]:
    """숫자로 된 폴더를 1, 2, ..., 10 순서로 정렬한다."""

    key: list[tuple[int, int | str]] = []
    for part in path.parts:
        key.append((0, int(part)) if part.isdigit() else (1, part.casefold()))
    return tuple(key)


def discover_shapefiles(input_dir: Path, output_path: Path) -> list[Path]:
    """출력 폴더를 제외한 입력 SHP 목록을 찾는다."""

    if not input_dir.is_dir():
        raise FileNotFoundError(f"입력 폴더가 없습니다: {input_dir}")

    input_dir = input_dir.resolve()
    output_dir = output_path.resolve().parent
    paths = [
        path.resolve()
        for path in input_dir.rglob("*.shp")
        if not path.resolve().is_relative_to(output_dir)
    ]
    paths.sort(key=lambda path: natural_path_key(path.relative_to(input_dir)))
    if not paths:
        raise FileNotFoundError(f"입력 SHP를 찾지 못했습니다: {input_dir}")
    return paths


def polygon_parts(geometry: Geometry | None) -> Iterable[Geometry]:
    """GeometryCollection을 포함한 형상에서 면 형상만 재귀적으로 꺼낸다."""

    if geometry is None or geometry.is_empty:
        return
    if geometry.geom_type == "Polygon":
        yield geometry
        return
    if geometry.geom_type in {"MultiPolygon", "GeometryCollection"}:
        for part in geometry.geoms:
            yield from polygon_parts(part)


def clean_polygonal_geometry(geometry: Geometry | None) -> Geometry | None:
    """형상을 2차원 유효 면 형상으로 정리한다."""

    if geometry is None or geometry.is_empty:
        return None
    geometry = shapely.force_2d(geometry)
    if not geometry.is_valid:
        geometry = shapely.make_valid(geometry)
    parts = [part for part in polygon_parts(geometry) if part.area > 0]
    if not parts:
        return None
    result = shapely.union_all(parts)
    if not result.is_valid:
        result = shapely.make_valid(result)
        result = shapely.union_all(list(polygon_parts(result)))
    return result if not result.is_empty else None


def read_and_dissolve_shapefile(
    path: Path, expected_crs: CRS | None
) -> tuple[Geometry, CRS, int]:
    """SHP 내부 피처를 모두 합쳐 파일당 정확히 하나의 면 형상을 만든다."""

    frame = gpd.read_file(path, engine="pyogrio")
    if frame.crs is None:
        raise ValueError(f"좌표계가 정의되지 않은 입력입니다: {path}")

    current_crs = CRS.from_user_input(frame.crs)
    if expected_crs is not None and not current_crs.equals(expected_crs):
        raise ValueError(
            "입력 SHP의 좌표계가 서로 다릅니다. "
            f"기준={expected_crs.to_string()}, 입력={current_crs.to_string()}, 파일={path}"
        )

    geometries: list[Geometry] = []
    for geometry in frame.geometry:
        cleaned = clean_polygonal_geometry(geometry)
        if cleaned is not None:
            geometries.append(cleaned)
    if not geometries:
        raise ValueError(f"유효한 폴리곤이 없는 입력입니다: {path}")

    dissolved = clean_polygonal_geometry(shapely.union_all(geometries))
    if dissolved is None:
        raise ValueError(f"입력 피처를 합친 결과가 비어 있습니다: {path}")
    return dissolved, current_crs, len(geometries)


def load_dissolved_shapefiles(paths: list[Path]) -> tuple[list[Geometry], CRS]:
    """각 SHP를 파일당 하나의 형상으로 dissolve해 읽는다."""

    if len(paths) != len(set(paths)):
        raise ValueError("동일한 입력 SHP 경로가 두 번 이상 발견되었습니다.")

    dissolved_polygons: list[Geometry] = []
    common_crs: CRS | None = None
    for number, path in enumerate(paths, start=1):
        dissolved, current_crs, feature_count = read_and_dissolve_shapefile(
            path, common_crs
        )
        common_crs = current_crs if common_crs is None else common_crs
        dissolved_polygons.append(dissolved)
        print(
            f"[{number:>2}/{len(paths)}] {path}: "
            f"{feature_count:,}개 피처 -> 1개 형상"
        )

    if common_crs is None:
        raise RuntimeError("좌표계를 확인할 입력이 없습니다.")
    if len(dissolved_polygons) != len(paths):
        raise RuntimeError("파일당 하나의 형상으로 통합하지 못했습니다.")
    return dissolved_polygons, common_crs


def select_cross_file_overlaps(
    dissolved_polygons: list[Geometry],
) -> tuple[list[Geometry], int]:
    """자기 자신을 제외하고 실제 면적이 겹치는 서로 다른 SHP만 선택한다.

    공간 인덱스 결과에서 ``왼쪽 인덱스 < 오른쪽 인덱스``인 쌍만 남긴다.
    이 조건으로 자기 자신 쌍과 A-B/B-A 중복 비교를 동시에 제거한다.
    경계선이나 꼭짓점만 닿는 경우는 중복 면적으로 보지 않는다.
    """

    polygon_array = np.asarray(dissolved_polygons, dtype=object)
    tree = STRtree(dissolved_polygons)
    candidate_pairs = tree.query(polygon_array, predicate="intersects")

    cross_file = candidate_pairs[0] < candidate_pairs[1]
    left = candidate_pairs[0, cross_file]
    right = candidate_pairs[1, cross_file]
    if left.size == 0:
        raise RuntimeError("서로 다른 SHP끼리 겹치는 영역이 없습니다.")

    intersections = shapely.intersection(polygon_array[left], polygon_array[right])
    has_area = (~shapely.is_empty(intersections)) & (shapely.area(intersections) > 0)
    left = left[has_area]
    right = right[has_area]
    if left.size == 0:
        raise RuntimeError("서로 다른 SHP끼리 면적으로 겹치는 영역이 없습니다.")

    active_indices = np.unique(np.concatenate((left, right)))
    active_polygons = [dissolved_polygons[int(index)] for index in active_indices]
    return active_polygons, int(left.size)


def build_atomic_areas(dissolved_polygons: list[Geometry]) -> list[Geometry]:
    """모든 파일 경계를 이용해 중복 수가 일정한 최소 영역들로 분할한다."""

    boundaries = [geometry.boundary for geometry in dissolved_polygons]
    noded_lines = shapely.union_all(boundaries)
    faces = shapely.polygonize(shapely.get_parts(noded_lines))
    return [part for part in shapely.get_parts(faces) if part.area > 0]


def count_dissolved_shapefiles(
    faces: list[Geometry], dissolved_polygons: list[Geometry]
) -> np.ndarray:
    """각 영역을 덮는 dissolve된 SHP 형상의 개수를 센다.

    ``dissolved_polygons``에는 파일당 형상이 하나만 있으므로 같은 SHP는 어떤
    영역에서도 최대 1회만 집계된다.
    """

    points = np.asarray([face.representative_point() for face in faces], dtype=object)
    tree = STRtree(dissolved_polygons)
    matches = tree.query(points, predicate="within")
    if matches.size == 0:
        return np.zeros(len(faces), dtype=np.int32)
    return np.bincount(matches[0], minlength=len(faces)).astype(np.int32)


def dissolve_by_count(
    faces: list[Geometry], counts: np.ndarray, minimum_overlap: int
) -> tuple[list[Geometry], list[int]]:
    """같은 중복 수를 가지며 맞닿은 최소 영역들을 하나로 합친다."""

    geometries: list[Geometry] = []
    values: list[int] = []
    face_array = np.asarray(faces, dtype=object)
    for count in sorted(int(value) for value in np.unique(counts) if value >= minimum_overlap):
        merged = clean_polygonal_geometry(shapely.union_all(face_array[counts == count]))
        for part in polygon_parts(merged):
            if part.area > 0:
                geometries.append(part)
                values.append(count)
    return geometries, values


def calculate_overlaps(
    dissolved_polygons: list[Geometry], crs: CRS, minimum_overlap: int
) -> gpd.GeoDataFrame:
    """서로 다른 SHP가 minimum_overlap개 이상 겹친 연결 영역을 반환한다."""

    active_polygons, pair_count = select_cross_file_overlaps(dissolved_polygons)
    print(
        f"서로 다른 SHP 간 면 중복: {pair_count:,}쌍 "
        f"(자기 자신 비교 제외, 참여 SHP {len(active_polygons):,}개)"
    )

    print("경계 교차점을 계산하고 영역을 분할하는 중...")
    faces = build_atomic_areas(active_polygons)
    if not faces:
        raise RuntimeError("분할된 폴리곤 영역이 없습니다.")

    print(f"분할 영역 {len(faces):,}개의 중복 수를 계산하는 중...")
    counts = count_dissolved_shapefiles(faces, active_polygons)
    if counts.size and int(counts.max()) > len(active_polygons):
        raise RuntimeError("동일 SHP가 중복 집계되었습니다.")
    geometries, values = dissolve_by_count(faces, counts, minimum_overlap)
    if not geometries:
        raise RuntimeError(f"{minimum_overlap}개 이상 겹치는 영역이 없습니다.")

    return gpd.GeoDataFrame(
        {EXACT_COUNT_FIELD: pd.Series(values, dtype="int32")},
        geometry=geometries,
        crs=crs,
    )


def write_outputs(result: gpd.GeoDataFrame, output_path: Path) -> Path:
    """SHP와 정확한 필드명을 보존하는 보조 GeoPackage를 저장한다."""

    if output_path.suffix.casefold() != ".shp":
        raise ValueError(f"출력 파일 확장자는 .shp여야 합니다: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    shp_result = result.rename(columns={EXACT_COUNT_FIELD: SHP_COUNT_FIELD})
    pyogrio.write_dataframe(
        shp_result,
        output_path,
        driver="ESRI Shapefile",
        encoding=SHP_ENCODING,
        promote_to_multi=True,
    )

    gpkg_path = output_path.with_suffix(".gpkg")
    pyogrio.write_dataframe(
        result,
        gpkg_path,
        layer=output_path.stem,
        driver="GPKG",
        promote_to_multi=True,
    )
    return gpkg_path


def validate_outputs(
    output_path: Path, gpkg_path: Path, expected_rows: int, expected_crs: CRS
) -> None:
    """저장 파일을 다시 읽어 행 수, 필드, 자료형, 좌표계를 검증한다."""

    shp = gpd.read_file(output_path, engine="pyogrio", encoding=SHP_ENCODING)
    gpkg = gpd.read_file(gpkg_path, engine="pyogrio", layer=output_path.stem)

    if len(shp) != expected_rows or len(gpkg) != expected_rows:
        raise RuntimeError(
            f"저장 후 피처 수가 달라졌습니다: 계산={expected_rows}, "
            f"SHP={len(shp)}, GPKG={len(gpkg)}"
        )
    if SHP_COUNT_FIELD not in shp.columns or not pd.api.types.is_integer_dtype(
        shp[SHP_COUNT_FIELD]
    ):
        raise RuntimeError(f"SHP 정수 필드 검증 실패: {SHP_COUNT_FIELD}")
    if EXACT_COUNT_FIELD not in gpkg.columns or not pd.api.types.is_integer_dtype(
        gpkg[EXACT_COUNT_FIELD]
    ):
        raise RuntimeError(f"GPKG 정수 필드 검증 실패: {EXACT_COUNT_FIELD}")
    if shp.crs is None or not CRS.from_user_input(shp.crs).equals(expected_crs):
        raise RuntimeError("SHP 좌표계 검증에 실패했습니다.")
    if gpkg.crs is None or not CRS.from_user_input(gpkg.crs).equals(expected_crs):
        raise RuntimeError("GPKG 좌표계 검증에 실패했습니다.")
    if shp.geometry.is_empty.any() or shp.geometry.isna().any():
        raise RuntimeError("SHP에 비어 있거나 누락된 형상이 있습니다.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="여러 폴리곤 SHP가 겹치는 영역과 영역별 중복 SHP 개수를 계산합니다."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input_dir(),
        help="입력 SHP를 재귀 탐색할 폴더",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_input_dir() / "중복" / "중복.shp",
        help="출력 Shapefile 경로",
    )
    parser.add_argument(
        "--minimum-overlap",
        type=int,
        default=2,
        help="출력할 최소 중복 SHP 개수(기본값: 2)",
    )
    args = parser.parse_args()
    if args.minimum_overlap < 2:
        parser.error("--minimum-overlap은 2 이상이어야 합니다.")
    return args


def main() -> int:
    args = parse_args()
    try:
        input_dir = args.input_dir.resolve()
        output_path = args.output.resolve()
        paths = discover_shapefiles(input_dir, output_path)
        print(f"입력 SHP: {len(paths)}개")
        print("처리 기준: SHP별 전체 피처 dissolve 후 서로 다른 파일끼리만 비교")

        dissolved_polygons, crs = load_dissolved_shapefiles(paths)
        result = calculate_overlaps(dissolved_polygons, crs, args.minimum_overlap)
        gpkg_path = write_outputs(result, output_path)
        validate_outputs(output_path, gpkg_path, len(result), crs)

        minimum = int(result[EXACT_COUNT_FIELD].min())
        maximum = int(result[EXACT_COUNT_FIELD].max())
        print(f"완료: {len(result):,}개 영역, 중복 수 {minimum}~{maximum}")
        print(f"SHP:  {output_path} (정수 필드: {SHP_COUNT_FIELD})")
        print(f"GPKG: {gpkg_path} (정수 필드: {EXACT_COUNT_FIELD})")
        return 0
    except Exception as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
