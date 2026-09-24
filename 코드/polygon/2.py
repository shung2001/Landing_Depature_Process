"""버티포트별로 다른 ID와의 중첩이 가장 적은 방위각 폴리곤을 선별한다.

입력
----
* ``자료/결과물/Polygon/12.5/<id>/Polygon_<id>.shp``
* ``자료/결과물/Polygon/12.5/중복/중복.gpkg``

``중복.gpkg``의 ``중복된_id``를 이용해 각 방위각 폴리곤이 면적으로 겹치는
다른 버티포트 ID의 종류 수를 계산한다. 같은 다른 ID가 여러 중복 영역에
나타나더라도 한 번만 집계한다. 각 ID에서 중첩 수가 적은 순서, 같은 경우에는
방위각이 작은 순서로 상위 10개를 선택하고, 중첩된 개수가 5개 이하인
폴리곤은 상위 10개 포함 여부와 관계없이 모두 선택한다.

출력
----
* ``자료/결과물/Polygon/안전지대/<id>/Polygon_<id>.shp``
* ``자료/결과물/Polygon/안전지대/안전지대_요약.csv``

Shapefile DBF의 필드명은 최대 10바이트이므로 요청한 ``중첩된 개수`` 대신
SHP에는 ``중첩개수``를 사용한다. 전체 요약 CSV에는 요청한 정확한 필드명
``중첩된 개수``를 사용한다. ``방위각``은 SHP와 CSV에서 동일하게 사용한다.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
import shapely
from pyproj import CRS
from shapely import STRtree


SOURCE_ID_FIELD = "heli_id"
SOURCE_BEARING_FIELD = "rot_deg"
OVERLAP_IDS_FIELD = "중복된_id"
OTHER_IDS_FIELD = "중첩된 id"
SHP_OVERLAP_COUNT_FIELD = "중첩개수"
EXACT_OVERLAP_COUNT_FIELD = "중첩된 개수"
BEARING_FIELD = "방위각"
RANK_FIELD = "순위"
SUMMARY_NAME = "안전지대_요약.csv"
SHP_ENCODING = "CP949"
ALWAYS_INCLUDE_MAX_OVERLAP = 5

Identifier = int | float | str


def project_root() -> Path:
    """이 스크립트의 위치를 기준으로 저장소 루트를 반환한다."""

    return Path(__file__).resolve().parents[2]


def default_input_dir() -> Path:
    return project_root() / "자료" / "결과물" / "Polygon" / "12.5"


def default_overlap_path() -> Path:
    return default_input_dir() / "중복" / "중복.gpkg"


def default_output_dir() -> Path:
    return project_root() / "자료" / "결과물" / "Polygon" / "12.5" / "안전지대"


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


def identifier_folder_name(identifier: Identifier) -> str:
    """ID를 안전한 단일 폴더명으로 변환한다."""

    text = str(identifier)
    if text in {"", ".", ".."} or any(character in text for character in '<>:"/\\|?*'):
        raise ValueError(f"폴더명으로 사용할 수 없는 ID입니다: {identifier!r}")
    return text


def discover_candidate_shapefiles(input_dir: Path) -> list[Path]:
    """중복 결과 폴더를 제외하고 ID별 입력 Shapefile을 찾는다."""

    if not input_dir.is_dir():
        raise FileNotFoundError(f"입력 폴더가 없습니다: {input_dir}")

    paths = [
        path.resolve()
        for path in input_dir.glob("*/*.shp")
        if path.parent.name != "중복"
    ]
    if not paths:
        raise FileNotFoundError(f"ID별 입력 Shapefile을 찾지 못했습니다: {input_dir}")

    def path_key(path: Path) -> tuple[int, float | str, str]:
        folder = path.parent.name
        identifier: Identifier = normalize_identifier(folder)
        group, value = identifier_sort_key(identifier)
        return group, value, path.name.casefold()

    paths.sort(key=path_key)
    return paths


def read_overlap_areas(
    overlap_path: Path,
) -> tuple[gpd.GeoDataFrame, list[frozenset[Identifier]], CRS]:
    """중복 영역과 영역별 ID 집합을 읽고 검증한다."""

    if not overlap_path.is_file():
        raise FileNotFoundError(f"중복 결과 파일이 없습니다: {overlap_path}")

    read_kwargs: dict[str, object] = {"engine": "pyogrio"}
    if overlap_path.suffix.casefold() == ".gpkg":
        read_kwargs["layer"] = overlap_path.stem
    else:
        read_kwargs["encoding"] = SHP_ENCODING
    frame = gpd.read_file(overlap_path, **read_kwargs)

    if frame.crs is None:
        raise ValueError(f"중복 결과에 좌표계가 없습니다: {overlap_path}")
    if OVERLAP_IDS_FIELD not in frame.columns:
        raise ValueError(
            f"중복 결과에 {OVERLAP_IDS_FIELD!r} 필드가 없습니다: {overlap_path}"
        )
    if frame.empty:
        raise ValueError(f"중복 결과가 비어 있습니다: {overlap_path}")
    if frame.geometry.isna().any() or frame.geometry.is_empty.any():
        raise ValueError(f"중복 결과에 비어 있거나 누락된 형상이 있습니다: {overlap_path}")

    identifier_sets: list[frozenset[Identifier]] = []
    for row_number, value in enumerate(frame[OVERLAP_IDS_FIELD], start=1):
        try:
            parsed = json.loads(value)
        except (TypeError, json.JSONDecodeError) as exc:
            raise ValueError(
                f"중복 결과 {row_number}행의 ID 목록을 해석할 수 없습니다: {value!r}"
            ) from exc
        if not isinstance(parsed, list) or len(parsed) < 2:
            raise ValueError(
                f"중복 결과 {row_number}행의 ID 목록은 2개 이상이어야 합니다: {value!r}"
            )
        identifiers = frozenset(normalize_identifier(item) for item in parsed)
        if len(identifiers) != len(parsed):
            raise ValueError(
                f"중복 결과 {row_number}행의 ID 목록에 중복이 있습니다: {value!r}"
            )
        identifier_sets.append(identifiers)

    return frame, identifier_sets, CRS.from_user_input(frame.crs)


def read_candidate_shapefile(
    path: Path, expected_crs: CRS
) -> tuple[gpd.GeoDataFrame, Identifier]:
    """방위각별 후보 Shapefile을 읽고 ID, 방위각, 좌표계를 검증한다."""

    frame = gpd.read_file(path, engine="pyogrio")
    if frame.crs is None:
        raise ValueError(f"입력 Shapefile에 좌표계가 없습니다: {path}")
    current_crs = CRS.from_user_input(frame.crs)
    if not current_crs.equals(expected_crs):
        raise ValueError(
            f"입력과 중복 결과의 좌표계가 다릅니다: "
            f"중복={expected_crs.to_string()}, 입력={current_crs.to_string()}, 파일={path}"
        )
    if frame.empty:
        raise ValueError(f"입력 Shapefile이 비어 있습니다: {path}")
    if frame.geometry.isna().any() or frame.geometry.is_empty.any():
        raise ValueError(f"입력에 비어 있거나 누락된 형상이 있습니다: {path}")

    missing = {SOURCE_ID_FIELD, SOURCE_BEARING_FIELD} - set(frame.columns)
    if missing:
        raise ValueError(f"입력 필드가 누락되었습니다: {sorted(missing)}, 파일={path}")

    identifiers = {
        normalize_identifier(value)
        for value in frame[SOURCE_ID_FIELD]
        if not pd.isna(value)
    }
    if len(identifiers) != 1:
        raise ValueError(
            f"입력 파일 하나에는 {SOURCE_ID_FIELD!r}가 정확히 하나여야 합니다: "
            f"값={sorted(identifiers, key=str)}, 파일={path}"
        )
    identifier = next(iter(identifiers))

    bearings = pd.to_numeric(frame[SOURCE_BEARING_FIELD], errors="coerce")
    if bearings.isna().any() or not np.isfinite(bearings.to_numpy(dtype=float)).all():
        raise ValueError(f"유효하지 않은 방위각이 있습니다: {path}")
    if bearings.duplicated().any():
        duplicates = sorted(bearings[bearings.duplicated(keep=False)].unique())
        raise ValueError(f"중복된 방위각이 있습니다: {duplicates}, 파일={path}")

    return frame, identifier


def calculate_overlap_scores(
    candidates: gpd.GeoDataFrame,
    current_identifier: Identifier,
    overlap_geometries: np.ndarray,
    overlap_identifier_sets: list[frozenset[Identifier]],
    overlap_tree: STRtree,
) -> tuple[np.ndarray, list[str]]:
    """후보별로 면적이 겹치는 다른 ID의 고유 개수와 목록을 계산한다."""

    candidate_geometries = np.asarray(candidates.geometry, dtype=object)
    pairs = overlap_tree.query(candidate_geometries, predicate="intersects")
    other_identifiers: list[set[Identifier]] = [set() for _ in candidates.index]

    if pairs.size:
        intersections = shapely.intersection(
            candidate_geometries[pairs[0]], overlap_geometries[pairs[1]]
        )
        has_area = (~shapely.is_empty(intersections)) & (shapely.area(intersections) > 0)
        for candidate_index, overlap_index in zip(
            pairs[0, has_area], pairs[1, has_area], strict=True
        ):
            identifiers = overlap_identifier_sets[int(overlap_index)]
            if current_identifier in identifiers:
                other_identifiers[int(candidate_index)].update(
                    identifiers - {current_identifier}
                )

    counts = np.fromiter(
        (len(identifiers) for identifiers in other_identifiers),
        dtype=np.int32,
        count=len(other_identifiers),
    )
    identifier_texts = [
        json.dumps(
            sorted(identifiers, key=identifier_sort_key),
            ensure_ascii=False,
            separators=(",", ":"),
        )
        for identifiers in other_identifiers
    ]
    return counts, identifier_texts


def select_safest_candidates(
    candidates: gpd.GeoDataFrame,
    counts: np.ndarray,
    other_identifier_texts: list[str],
    top_n: int,
) -> tuple[gpd.GeoDataFrame, pd.DataFrame]:
    """상위 N개와 중첩 수 5 이하인 안전한 후보를 선택한다."""

    working = candidates.copy()
    working[SHP_OVERLAP_COUNT_FIELD] = pd.Series(counts, index=working.index, dtype="int32")
    working[BEARING_FIELD] = pd.to_numeric(
        working[SOURCE_BEARING_FIELD], errors="raise"
    ).astype("float64")
    working["__other_ids"] = other_identifier_texts
    working["__source_order"] = np.arange(len(working), dtype=np.int32)

    ranked = working.sort_values(
        [SHP_OVERLAP_COUNT_FIELD, BEARING_FIELD, "__source_order"],
        kind="stable",
    )
    top_indices = ranked.head(top_n).index
    selected = ranked.loc[
        (ranked[SHP_OVERLAP_COUNT_FIELD] <= ALWAYS_INCLUDE_MAX_OVERLAP)
        | ranked.index.isin(top_indices)
    ].copy()
    details = pd.DataFrame(
        {
            RANK_FIELD: np.arange(1, len(selected) + 1, dtype=np.int32),
            EXACT_OVERLAP_COUNT_FIELD: selected[SHP_OVERLAP_COUNT_FIELD].to_numpy(
                dtype=np.int32
            ),
            BEARING_FIELD: selected[BEARING_FIELD].to_numpy(dtype=float),
            OTHER_IDS_FIELD: selected["__other_ids"].to_list(),
        }
    )
    selected = selected.drop(columns=["__other_ids", "__source_order"])
    return selected, details


def write_selected_shapefile(
    selected: gpd.GeoDataFrame, identifier: Identifier, output_dir: Path
) -> Path:
    """선택된 후보를 ID별 폴더의 Shapefile로 저장한다."""

    folder_name = identifier_folder_name(identifier)
    id_dir = output_dir / folder_name
    id_dir.mkdir(parents=True, exist_ok=True)
    output_path = id_dir / f"Polygon_{folder_name}.shp"
    pyogrio.write_dataframe(
        selected,
        output_path,
        driver="ESRI Shapefile",
        encoding=SHP_ENCODING,
        promote_to_multi=True,
    )
    return output_path


def validate_selected_output(
    output_path: Path,
    expected_rows: int,
    expected_crs: CRS,
    expected_counts: np.ndarray,
    expected_bearings: np.ndarray,
) -> None:
    """저장한 Shapefile을 다시 읽어 행, 필드, 값과 좌표계를 검증한다."""

    saved = gpd.read_file(output_path, engine="pyogrio", encoding=SHP_ENCODING)
    if len(saved) != expected_rows:
        raise RuntimeError(
            f"저장 후 폴리곤 수가 다릅니다: 예상={expected_rows}, "
            f"저장={len(saved)}, 파일={output_path}"
        )
    missing = {SHP_OVERLAP_COUNT_FIELD, BEARING_FIELD} - set(saved.columns)
    if missing:
        raise RuntimeError(f"저장 필드가 누락되었습니다: {sorted(missing)}, 파일={output_path}")
    if not np.array_equal(
        saved[SHP_OVERLAP_COUNT_FIELD].to_numpy(dtype=np.int32), expected_counts
    ):
        raise RuntimeError(f"저장된 중첩 개수가 계산 결과와 다릅니다: {output_path}")
    if not np.allclose(
        saved[BEARING_FIELD].to_numpy(dtype=float), expected_bearings, rtol=0, atol=1e-9
    ):
        raise RuntimeError(f"저장된 방위각이 계산 결과와 다릅니다: {output_path}")
    if saved.geometry.isna().any() or saved.geometry.is_empty.any():
        raise RuntimeError(f"저장 결과에 비어 있거나 누락된 형상이 있습니다: {output_path}")
    if saved.crs is None or not CRS.from_user_input(saved.crs).equals(expected_crs):
        raise RuntimeError(f"저장 결과의 좌표계가 다릅니다: {output_path}")


def build_summary_rows(
    identifier: Identifier,
    details: pd.DataFrame,
    input_count: int,
    input_path: Path,
    output_path: Path,
) -> list[dict[str, object]]:
    """선택된 폴리곤 한 개당 하나의 전체 요약 행을 만든다."""

    rows: list[dict[str, object]] = []
    for record in details.to_dict(orient="records"):
        rows.append(
            {
                "id": identifier,
                RANK_FIELD: record[RANK_FIELD],
                EXACT_OVERLAP_COUNT_FIELD: record[EXACT_OVERLAP_COUNT_FIELD],
                BEARING_FIELD: record[BEARING_FIELD],
                OTHER_IDS_FIELD: record[OTHER_IDS_FIELD],
                "전체 후보 개수": input_count,
                "입력 파일": str(input_path),
                "출력 파일": str(output_path),
            }
        )
    return rows


def write_and_validate_summary(rows: list[dict[str, object]], output_dir: Path) -> Path:
    """전체 선택 결과를 CSV로 저장한 뒤 핵심 필드를 검증한다."""

    if not rows:
        raise RuntimeError("요약 CSV에 저장할 결과가 없습니다.")
    summary = pd.DataFrame(rows)
    summary_path = output_dir / SUMMARY_NAME
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    saved = pd.read_csv(summary_path, encoding="utf-8-sig")
    required = {
        "id",
        RANK_FIELD,
        EXACT_OVERLAP_COUNT_FIELD,
        BEARING_FIELD,
        OTHER_IDS_FIELD,
    }
    missing = required - set(saved.columns)
    if missing or len(saved) != len(summary):
        raise RuntimeError(
            f"요약 CSV 검증에 실패했습니다: 누락={sorted(missing)}, "
            f"계산={len(summary)}, 저장={len(saved)}"
        )
    return summary_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="ID별로 다른 ID와의 중첩이 가장 적은 방위각 폴리곤을 선택합니다."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input_dir(),
        help="ID별 Polygon Shapefile이 있는 경사 폴더",
    )
    parser.add_argument(
        "--overlap",
        type=Path,
        default=default_overlap_path(),
        help="중복된_id 필드가 있는 중복 GPKG 또는 SHP",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir(),
        help="ID별 안전지대 폴더와 전체 요약 CSV를 저장할 폴더",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="ID별로 선택할 폴리곤 개수(기본값: 10)",
    )
    args = parser.parse_args()
    if args.top_n < 1:
        parser.error("--top-n은 1 이상이어야 합니다.")
    return args


def main() -> int:
    args = parse_args()
    try:
        input_dir = args.input_dir.resolve()
        overlap_path = args.overlap.resolve()
        output_dir = args.output_dir.resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        paths = discover_candidate_shapefiles(input_dir)
        overlap_frame, overlap_identifier_sets, overlap_crs = read_overlap_areas(
            overlap_path
        )
        overlap_geometries = np.asarray(overlap_frame.geometry, dtype=object)
        overlap_tree = STRtree(overlap_geometries)

        print(f"입력 ID Shapefile: {len(paths)}개")
        print(f"중복 영역: {len(overlap_frame):,}개")
        print(
            "선정 기준: 다른 ID의 고유 중첩 수 오름차순, "
            f"방위각 오름차순, ID별 상위 {args.top_n}개 + "
            f"중첩 {ALWAYS_INCLUDE_MAX_OVERLAP}개 이하 전체"
        )

        seen_identifiers: set[Identifier] = set()
        summary_rows: list[dict[str, object]] = []
        for number, path in enumerate(paths, start=1):
            candidates, identifier = read_candidate_shapefile(path, overlap_crs)
            if identifier in seen_identifiers:
                raise ValueError(f"동일한 ID의 입력 파일이 두 개 이상입니다: id={identifier}")
            seen_identifiers.add(identifier)

            counts, other_identifier_texts = calculate_overlap_scores(
                candidates,
                identifier,
                overlap_geometries,
                overlap_identifier_sets,
                overlap_tree,
            )
            selected, details = select_safest_candidates(
                candidates, counts, other_identifier_texts, args.top_n
            )
            output_path = write_selected_shapefile(selected, identifier, output_dir)
            validate_selected_output(
                output_path,
                len(selected),
                overlap_crs,
                selected[SHP_OVERLAP_COUNT_FIELD].to_numpy(dtype=np.int32),
                selected[BEARING_FIELD].to_numpy(dtype=float),
            )
            summary_rows.extend(
                build_summary_rows(
                    identifier,
                    details,
                    len(candidates),
                    path,
                    output_path,
                )
            )

            minimum = int(selected[SHP_OVERLAP_COUNT_FIELD].min())
            maximum = int(selected[SHP_OVERLAP_COUNT_FIELD].max())
            bearings = ", ".join(f"{value:g}" for value in selected[BEARING_FIELD])
            print(
                f"[{number:>2}/{len(paths)}] id={identifier}: "
                f"{len(candidates)}개 중 {len(selected)}개 선택, "
                f"중첩 {minimum}~{maximum}, 방위각 [{bearings}]"
            )

        summary_path = write_and_validate_summary(summary_rows, output_dir)
        print(f"완료: ID {len(seen_identifiers)}개, 선택 폴리곤 {len(summary_rows):,}개")
        print(f"출력 폴더: {output_dir}")
        print(f"요약 CSV: {summary_path}")
        return 0
    except Exception as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())