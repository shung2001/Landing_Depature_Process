"""안전지대에 선택된 각 폴리곤의 중첩 개수를 계산한다.

기본 입력::

    자료/결과물/Polygon/12.5/안전지대/<id>/Polygon_<id>.shp

기본 출력::

    자료/결과물/Polygon/12.5/안전지대/안전지대_중첩.csv

각 입력 피처를 독립된 폴리곤으로 취급하고, 다른 버티포트 ID의 폴리곤과
실제 면적이 겹치는 횟수를 ``중첩개수``에 기록한다. 같은 버티포트 안에서
서로 다른 방위각끼리 겹치는 것은 기본적으로 제외한다. 경계선이나 꼭짓점만
닿는 경우도 면 중첩으로 집계하지 않는다.

``중첩개수``는 겹치는 상대 폴리곤 피처의 수이고, ``중첩ID개수``는 그 상대
폴리곤들이 속한 고유 버티포트 ID의 수이다. ``중첩된 id``에는 해당 ID들을
JSON 배열 문자열로 저장한다.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Iterable
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from pyproj import CRS
from shapely import Geometry, STRtree


SOURCE_ID_FIELD = "heli_id"
BEARING_FIELD = "방위각"
SOURCE_BEARING_FIELD = "rot_deg"
OVERLAP_COUNT_FIELD = "중첩개수"
OVERLAP_ID_COUNT_FIELD = "중첩ID개수"
OVERLAP_IDS_FIELD = "중첩된 id"
OUTPUT_NAME = "안전지대_중첩.csv"
SHP_ENCODING = "CP949"
INTERSECTION_CHUNK_SIZE = 20_000

Identifier = int | float | str


def project_root() -> Path:
    """이 스크립트의 위치를 기준으로 저장소 루트를 반환한다."""

    return Path(__file__).resolve().parents[2]


def default_input_dir() -> Path:
    return (
        project_root()
        / "자료"
        / "결과물"
        / "Polygon"
        / "12.5"
        / "안전지대"
    )


def default_output_path() -> Path:
    return default_input_dir() / OUTPUT_NAME


def normalize_identifier(value: object) -> Identifier:
    """ID를 비교와 JSON 저장에 사용할 수 있는 스칼라로 정규화한다."""

    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return int(number) if number.is_integer() else number

    text = str(value).strip()
    if not text:
        raise ValueError("빈 ID는 사용할 수 없습니다.")
    if text.lstrip("+-").isdigit():
        return int(text)
    return text


def identifier_sort_key(identifier: Identifier) -> tuple[int, float | str]:
    """숫자 ID를 먼저 숫자 순으로, 나머지를 문자열 순으로 정렬한다."""

    if isinstance(identifier, (int, float)):
        return 0, float(identifier)
    return 1, identifier.casefold()


def natural_path_key(path: Path) -> tuple[tuple[int, int | str], ...]:
    """숫자 폴더가 1, 2, ..., 10 순서가 되도록 경로 정렬 키를 만든다."""

    return tuple(
        (0, int(part)) if part.isdigit() else (1, part.casefold())
        for part in path.parts
    )


def discover_shapefiles(input_dir: Path) -> list[Path]:
    """안전지대의 ID별 하위 폴더에서 Shapefile을 찾는다."""

    if not input_dir.is_dir():
        raise FileNotFoundError(f"안전지대 입력 폴더가 없습니다: {input_dir}")

    paths = [path.resolve() for path in input_dir.glob("*/*.shp")]
    paths.sort(key=lambda path: natural_path_key(path.relative_to(input_dir.resolve())))
    if not paths:
        raise FileNotFoundError(f"안전지대 Shapefile을 찾지 못했습니다: {input_dir}")
    if len(paths) != len(set(paths)):
        raise ValueError("동일한 안전지대 Shapefile이 두 번 이상 발견되었습니다.")
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
    """입력 형상을 유효한 2차원 Polygon 또는 MultiPolygon으로 정리한다."""

    if geometry is None or geometry.is_empty:
        return None
    cleaned = shapely.force_2d(geometry)
    if not cleaned.is_valid:
        cleaned = shapely.make_valid(cleaned)
    parts = [part for part in polygon_parts(cleaned) if part.area > 0]
    if not parts:
        return None
    result = shapely.union_all(parts)
    return result if not result.is_empty and result.area > 0 else None


def read_safety_polygons(
    paths: list[Path],
) -> tuple[pd.DataFrame, np.ndarray, CRS]:
    """모든 안전지대 피처의 식별 정보와 정리된 형상을 읽는다."""

    records: list[dict[str, object]] = []
    geometries: list[Geometry] = []
    common_crs: CRS | None = None

    for file_number, path in enumerate(paths, start=1):
        frame = gpd.read_file(path, engine="pyogrio", encoding=SHP_ENCODING)
        if frame.empty:
            raise ValueError(f"비어 있는 안전지대 Shapefile입니다: {path}")
        if frame.crs is None:
            raise ValueError(f"좌표계가 없는 안전지대 Shapefile입니다: {path}")

        current_crs = CRS.from_user_input(frame.crs)
        if common_crs is None:
            common_crs = current_crs
        elif not current_crs.equals(common_crs):
            raise ValueError(
                "안전지대 Shapefile의 좌표계가 서로 다릅니다. "
                f"기준={common_crs.to_string()}, "
                f"입력={current_crs.to_string()}, 파일={path}"
            )

        if SOURCE_ID_FIELD not in frame.columns:
            raise ValueError(f"{SOURCE_ID_FIELD!r} 필드가 없습니다: {path}")
        bearing_source = (
            BEARING_FIELD if BEARING_FIELD in frame.columns else SOURCE_BEARING_FIELD
        )
        if bearing_source not in frame.columns:
            raise ValueError(
                f"{BEARING_FIELD!r} 또는 {SOURCE_BEARING_FIELD!r} 필드가 없습니다: {path}"
            )

        identifiers = frame[SOURCE_ID_FIELD].map(normalize_identifier)
        if identifiers.nunique(dropna=False) != 1:
            values = sorted(set(identifiers), key=identifier_sort_key)
            raise ValueError(
                f"한 파일에는 {SOURCE_ID_FIELD!r}가 하나만 있어야 합니다: "
                f"값={values}, 파일={path}"
            )
        identifier = identifiers.iloc[0]

        bearings = pd.to_numeric(frame[bearing_source], errors="coerce")
        if bearings.isna().any() or not np.isfinite(
            bearings.to_numpy(dtype=float)
        ).all():
            raise ValueError(f"유효하지 않은 방위각이 있습니다: {path}")

        valid_in_file = 0
        for source_row, (bearing, geometry) in enumerate(
            zip(bearings, frame.geometry, strict=True), start=1
        ):
            cleaned = clean_polygonal_geometry(geometry)
            if cleaned is None:
                raise ValueError(
                    f"유효한 면 형상이 아닌 피처가 있습니다: "
                    f"행={source_row}, 파일={path}"
                )
            records.append(
                {
                    "id": identifier,
                    BEARING_FIELD: float(bearing),
                    "입력 파일": str(path),
                    "입력 행번호": source_row,
                }
            )
            geometries.append(cleaned)
            valid_in_file += 1

        print(
            f"[{file_number:>2}/{len(paths)}] id={identifier}: "
            f"안전지대 polygon {valid_in_file:,}개"
        )

    if common_crs is None or not records:
        raise RuntimeError("중첩을 계산할 안전지대 polygon이 없습니다.")

    metadata = pd.DataFrame.from_records(records)
    duplicated = metadata.duplicated(subset=["id", BEARING_FIELD], keep=False)
    if duplicated.any():
        values = metadata.loc[duplicated, ["id", BEARING_FIELD]].to_dict(
            orient="records"
        )
        raise ValueError(f"같은 ID에 중복된 방위각이 있습니다: {values[:10]}")

    return metadata, np.asarray(geometries, dtype=object), common_crs


def calculate_overlap_counts(
    metadata: pd.DataFrame,
    geometries: np.ndarray,
    include_same_id: bool = False,
) -> tuple[np.ndarray, list[set[Identifier]], int]:
    """각 폴리곤과 면적으로 겹치는 상대 폴리곤 수와 ID 집합을 계산한다."""

    if len(metadata) != len(geometries):
        raise ValueError("폴리곤 메타데이터와 형상 개수가 서로 다릅니다.")

    tree = STRtree(geometries)
    pairs = tree.query(geometries, predicate="intersects")
    unique_pairs = pairs[0] < pairs[1]
    left = pairs[0, unique_pairs]
    right = pairs[1, unique_pairs]

    identifiers = metadata["id"].to_numpy(dtype=object)
    if not include_same_id:
        different_id = identifiers[left] != identifiers[right]
        left = left[different_id]
        right = right[different_id]

    counts = np.zeros(len(metadata), dtype=np.int32)
    overlap_ids: list[set[Identifier]] = [set() for _ in range(len(metadata))]
    area_pair_count = 0

    for start in range(0, len(left), INTERSECTION_CHUNK_SIZE):
        stop = min(start + INTERSECTION_CHUNK_SIZE, len(left))
        chunk_left = left[start:stop]
        chunk_right = right[start:stop]
        intersections = shapely.intersection(
            geometries[chunk_left], geometries[chunk_right]
        )
        has_area = (~shapely.is_empty(intersections)) & (
            shapely.area(intersections) > 0
        )
        active_left = chunk_left[has_area]
        active_right = chunk_right[has_area]
        if active_left.size == 0:
            continue

        np.add.at(counts, active_left, 1)
        np.add.at(counts, active_right, 1)
        area_pair_count += int(active_left.size)
        for left_index, right_index in zip(active_left, active_right, strict=True):
            left_number = int(left_index)
            right_number = int(right_index)
            overlap_ids[left_number].add(identifiers[right_number])
            overlap_ids[right_number].add(identifiers[left_number])

    return counts, overlap_ids, area_pair_count


def build_result(
    metadata: pd.DataFrame,
    counts: np.ndarray,
    overlap_ids: list[set[Identifier]],
) -> pd.DataFrame:
    """계산 결과를 CSV용 테이블로 만든다."""

    if len(metadata) != len(counts) or len(metadata) != len(overlap_ids):
        raise ValueError("중첩 결과의 길이가 입력 polygon 개수와 다릅니다.")

    result = metadata.copy()
    result[OVERLAP_COUNT_FIELD] = counts.astype(np.int32)
    result[OVERLAP_ID_COUNT_FIELD] = np.fromiter(
        (len(values) for values in overlap_ids),
        dtype=np.int32,
        count=len(overlap_ids),
    )
    result[OVERLAP_IDS_FIELD] = [
        json.dumps(
            sorted(values, key=identifier_sort_key),
            ensure_ascii=False,
            separators=(",", ":"),
        )
        for values in overlap_ids
    ]
    return result[
        [
            "id",
            BEARING_FIELD,
            OVERLAP_COUNT_FIELD,
            OVERLAP_ID_COUNT_FIELD,
            OVERLAP_IDS_FIELD,
            "입력 파일",
            "입력 행번호",
        ]
    ]


def write_and_validate_result(result: pd.DataFrame, output_path: Path) -> None:
    """중첩 결과를 UTF-8 BOM CSV로 저장하고 핵심 필드와 값을 검증한다."""

    if result.empty:
        raise RuntimeError("CSV에 저장할 중첩 결과가 없습니다.")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False, encoding="utf-8-sig")

    saved = pd.read_csv(output_path, encoding="utf-8-sig")
    required = {
        "id",
        BEARING_FIELD,
        OVERLAP_COUNT_FIELD,
        OVERLAP_ID_COUNT_FIELD,
        OVERLAP_IDS_FIELD,
    }
    missing = required - set(saved.columns)
    if missing or len(saved) != len(result):
        raise RuntimeError(
            f"출력 CSV 검증 실패: 누락={sorted(missing)}, "
            f"계산={len(result)}, 저장={len(saved)}, 파일={output_path}"
        )
    if not np.array_equal(
        saved[OVERLAP_COUNT_FIELD].to_numpy(dtype=np.int32),
        result[OVERLAP_COUNT_FIELD].to_numpy(dtype=np.int32),
    ):
        raise RuntimeError(f"저장된 중첩개수가 계산 결과와 다릅니다: {output_path}")
    if not np.array_equal(
        saved[OVERLAP_ID_COUNT_FIELD].to_numpy(dtype=np.int32),
        result[OVERLAP_ID_COUNT_FIELD].to_numpy(dtype=np.int32),
    ):
        raise RuntimeError(f"저장된 중첩ID개수가 계산 결과와 다릅니다: {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="선택된 안전지대의 각 polygon이 다른 안전지대 polygon과 겹치는 개수를 계산합니다."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input_dir(),
        help="ID별 안전지대 Shapefile 하위 폴더가 있는 경로",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_output_path(),
        help="중첩 결과 CSV 경로",
    )
    parser.add_argument(
        "--include-same-id",
        action="store_true",
        help="같은 버티포트 ID의 서로 다른 방위각 polygon끼리도 중첩으로 집계",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        input_dir = args.input_dir.resolve()
        output_path = args.output.resolve()
        paths = discover_shapefiles(input_dir)

        print(f"입력 안전지대 Shapefile: {len(paths)}개")
        print(
            "중첩 기준: 자기 자신 제외, 면적이 0보다 큰 교차만 집계, "
            + (
                "같은 ID 포함"
                if args.include_same_id
                else "같은 ID의 다른 방위각 polygon 제외"
            )
        )
        metadata, geometries, crs = read_safety_polygons(paths)
        counts, overlap_ids, pair_count = calculate_overlap_counts(
            metadata,
            geometries,
            include_same_id=args.include_same_id,
        )
        result = build_result(metadata, counts, overlap_ids)
        write_and_validate_result(result, output_path)

        print(f"좌표계: {crs.to_string()}")
        print(f"입력 polygon: {len(result):,}개")
        print(f"면적으로 겹치는 polygon 쌍: {pair_count:,}쌍")
        print(
            f"polygon별 중첩개수: 최소 {int(counts.min()):,}, "
            f"최대 {int(counts.max()):,}, 평균 {float(counts.mean()):,.2f}"
        )
        print(f"출력 CSV: {output_path}")
        return 0
    except Exception as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
