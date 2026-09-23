"""송파구 옥상 헬리포트별 회전 삼각형 장애물 제한 표면(OLS)을 생성한다.

헬리포트 Shapefile의 ``id``, ``x_long``, ``y_lati``, ``A16``을 읽어 각 id를
기준 꼭짓점으로 사용한다. 꼭짓점은 고정하고 진북 기준 시계방향 0~359도로
회전한 360개 Polygon Z를 id별 Shapefile 하나에 저장한다. 건물 Shapefile은
각 헬리포트의 최근접 host 건물과 거리를 검증하는 데 사용한다. 0~359도 후보를
메모리에서 먼저 생성한 뒤, host를 제외한 건물 또는 자기 자신을 제외한
헬리포트가 OLS 표면고에 닿거나 위로 관통하는 후보만 제거하고 안전한 Polygon만
최종 Shapefile로 저장한다. XY 평면에서 footprint가 겹치는 것만으로는 충돌로
판정하지 않으며, 생성 후보끼리도 intersection 검사하지 않는다. 검증 전 후보
Shapefile은 디스크에 저장하지 않는다.

기본 형상
---------
* 수평 중심선 길이(삼각형 높이): 1,000 m
* 밑변 폭: 573 m
* OLS 종단 경사각: 12.5도
* 수평 좌표계: WGS84 (EPSG:4326, 좌표 순서 경도/위도)
* 기준점 AGL: 헬리포트 ``A16``
* Polygon Z: 기준점 DEM 지표고 + ``A16`` + OLS 상승량
* 충돌: 장애물의 ``DEM + A16`` 상단고가 중첩 위치의 OLS 표면고 이상인 경우
* 회전 범위: 0도 이상 360도 미만, 1도 간격(총 360개 Polygon)
* 출력: ``자료/결과물/Polygon/<OLS 경사>/<id>/Polygon_<id>.shp``
* 충돌 기록: ``자료/결과물/Polygon/<OLS 경사>/generation_summary.csv``

안전 방위각이 하나 이상 남은 id만 Shapefile과 id 폴더를 출력한다. 모든
방위각이 충돌한 id는 Shapefile을 만들지 않지만 요약 CSV에는 포함한다.

``A16``이 결측 또는 0 이하이면 동일 좌표·동일 Geometry의 다른 헬리포트
레코드에 있는 유효한 ``A16``으로 보완한다. 보완 여부는 ``agl_src`` 필드와
``generation_summary.csv``에 기록한다. 건물 ``A16``은 host 검증 속성으로만
저장하고 헬리포트 AGL을 임의로 대체하지 않는다.

EPSG:4326은 수평 2차원 좌표계이므로 Polygon Z는 EPSG:4326 타원체고가 아니다.
기준점 DEM 지표고에 헬리포트 A16을 더해 꼭짓점 표면고를 구한 뒤, OLS
상승량을 적용한다. ``z_ref`` 속성에는 ``DEM+AGL``임을 명시한다.

실행 예시::

    python -B "코드/polygon/장애물접근표면.py"

위 명령은 기본 건물·헬리포트·DEM을 읽고 OLS 경사/id별 폴더와 Shapefile을
만든다. 각 Shapefile에는 0~359도 중 3차원 높이 충돌이 없는 ``rot_deg``만
기록된다. 건물 ``A16``이 0 또는 결측이면 높이 미상으로 보고, XY 중첩만으로
해당 방위각을 제거하지 않는다.
충돌각은 ``[0,10] & [20,30]`` 형식으로 요약 CSV에 저장한다. 0도와 360도는
같은 방향이므로 기본 출력에서는 360도를 중복 생성하지 않는다. 입력 경로와
필드명은 아래 ``SCRIPT_CONFIG`` 또는 명령행 인수에서 바꿀 수 있다.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
import os
from pathlib import Path
import re
from typing import Sequence

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Geod, Transformer
import rasterio
from rasterio.windows import Window
import shapely
from shapely.geometry import Polygon


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SONGPA_DATA_DIR = PROJECT_ROOT / "자료" / "안전성" / "송파구"
DEFAULT_BUILDING_PATH = SONGPA_DATA_DIR / "송파구_건물_4326.shp"
DEFAULT_HELIPAD_PATH = SONGPA_DATA_DIR / "송파구_옥상헬리포트_4326.shp"
DEFAULT_DEM_PATH = SONGPA_DATA_DIR / "송파구_DEM.tif"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "자료" / "결과물" / "Polygon"


@dataclass(frozen=True)
class ScriptConfig:
    """명령행 인수 없이 실행할 때 사용할 스크립트 내부 설정."""

    building_path: Path = DEFAULT_BUILDING_PATH
    helipad_path: Path = DEFAULT_HELIPAD_PATH
    dem_path: Path = DEFAULT_DEM_PATH
    output_root: Path = DEFAULT_OUTPUT_ROOT
    id_field: str = "id"
    longitude_field: str = "x_long"
    latitude_field: str = "y_lati"
    agl_field: str = "A16"
    host_match_crs: str = "EPSG:5186"
    start_rotation_deg: float = 0.0
    end_rotation_deg: float = 360.0
    rotation_step_deg: float = 1.0
    include_end_angle: bool = False
    triangle_height_m: float = 1_000.0
    base_width_m: float = 573.0
    ols_slope_deg: float = 8.0
    dem_band: int = 1


# =========================== 스크립트 설정 ===========================
# 아래 입력 경로와 필드명을 바꾸면 parser 인수 없이 일괄 실행할 수 있다.
# 명령행 인수를 함께 지정하면 해당 명령행 값이 이 설정보다 우선한다.
SCRIPT_CONFIG = ScriptConfig(
    building_path=DEFAULT_BUILDING_PATH,
    helipad_path=DEFAULT_HELIPAD_PATH,
    dem_path=DEFAULT_DEM_PATH,
    output_root=DEFAULT_OUTPUT_ROOT,
    id_field="id",
    longitude_field="x_long",
    latitude_field="y_lati",
    agl_field="A16",
    host_match_crs="EPSG:5186",
    start_rotation_deg=0.0,
    end_rotation_deg=360.0,
    rotation_step_deg=1.0,
    include_end_angle=False,
    triangle_height_m=1_000.0,
    base_width_m=573.0,
    ols_slope_deg=8.0,
    dem_band=1,
)
# ===================================================================


WGS84_GEOD = Geod(ellps="WGS84")


@dataclass(frozen=True)
class Vertex:
    """삼각형 꼭짓점 하나의 WGS84 수평좌표와 절대 표면고."""

    name: str
    latitude: float
    longitude: float
    z_m: float

    def xyz_coordinate(self) -> tuple[float, float, float]:
        """Shapefile Polygon Z 좌표 순서인 (경도, 위도, 표면고)을 반환한다."""

        return self.longitude, self.latitude, self.z_m


@dataclass(frozen=True)
class DemSample:
    """기준점 위치에서 읽은 DEM 셀 정보."""

    elevation_m: float
    path: Path
    crs: str
    band: int
    row: int
    column: int
    x: float
    y: float


@dataclass(frozen=True)
class HelipadRecord:
    """Shapefile 한 행에서 검증해 얻은 헬리포트 생성 정보."""

    identifier: str
    folder_name: str
    source_index: int
    latitude: float
    longitude: float
    agl_m: float
    agl_source: str
    host_index: int
    host_distance_m: float
    host_a16_m: float
    host_a17_m: float
    excluded_helipad_indices: tuple[int, ...]


@dataclass(frozen=True)
class GenerationResult:
    """헬리포트 하나의 출력 결과 요약."""

    helipad: HelipadRecord
    dem_sample: DemSample
    output_path: Path | None
    polygon_count: int
    collision_angles: tuple[float, ...]


@dataclass(frozen=True)
class OLSTriangle:
    """꼭짓점 세 개와 생성에 사용한 OLS 형상 정보."""

    apex: Vertex
    base_right: Vertex
    base_left: Vertex
    base_center: Vertex
    rotation_deg: float
    horizontal_height_m: float
    base_width_m: float
    slope_deg: float
    apex_agl_m: float
    ground_z_m: float

    @property
    def vertices(self) -> tuple[Vertex, Vertex, Vertex]:
        """외곽링 반시계 방향 순서의 꼭짓점 세 개를 반환한다."""

        return self.apex, self.base_right, self.base_left

    @property
    def polygon(self) -> Polygon:
        """경도/위도/DEM 기반 표면고의 3차원 Shapely Polygon을 반환한다."""

        return Polygon([vertex.xyz_coordinate() for vertex in self.vertices])

    @property
    def vertical_rise_m(self) -> float:
        """기준 꼭짓점에서 밑변까지의 수직 상승량."""

        return self.base_center.z_m - self.apex.z_m

    @property
    def apex_z_m(self) -> float:
        """DEM 지표고와 입력 AGL을 합한 기준 꼭짓점 표면고."""

        return self.apex.z_m

    @property
    def base_z_m(self) -> float:
        """OLS 경사를 적용한 밑변 표면고."""

        return self.base_center.z_m


@dataclass(frozen=True)
class CollisionLayers:
    """3차원 충돌 검사에 사용하는 투영 도형과 장애물 상단고."""

    analysis_crs: str
    buildings: gpd.GeoDataFrame
    helipads: gpd.GeoDataFrame
    building_top_z_m: np.ndarray
    helipad_top_z_m: np.ndarray


def _finite_number(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name}은(는) 유한한 숫자여야 합니다: {value!r}")
    return value


def _identifier_text(value: object) -> str:
    """Shapefile id 값을 폴더와 속성에 사용할 안정적인 문자열로 바꾼다."""

    if pd.isna(value):
        raise ValueError("헬리포트 id에 결측값이 있습니다.")
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    text = str(value).strip()
    if not text:
        raise ValueError("헬리포트 id에 빈 문자열이 있습니다.")
    return text


def _safe_folder_name(identifier: str) -> str:
    """id를 Windows에서도 안전한 단일 폴더명으로 검증·정리한다."""

    folder = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", identifier).strip().rstrip(".")
    if not folder or folder in {".", ".."}:
        raise ValueError(f"폴더명으로 사용할 수 없는 id입니다: {identifier!r}")
    reserved = {
        "CON", "PRN", "AUX", "NUL",
        *(f"COM{number}" for number in range(1, 10)),
        *(f"LPT{number}" for number in range(1, 10)),
    }
    if folder.upper() in reserved:
        folder = f"id_{folder}"
    return folder


def load_source_layers(
    building_path: str | Path = SCRIPT_CONFIG.building_path,
    helipad_path: str | Path = SCRIPT_CONFIG.helipad_path,
    id_field: str = SCRIPT_CONFIG.id_field,
    longitude_field: str = SCRIPT_CONFIG.longitude_field,
    latitude_field: str = SCRIPT_CONFIG.latitude_field,
    agl_field: str = SCRIPT_CONFIG.agl_field,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """건물과 헬리포트 Shapefile을 읽고 기본 구조를 검증한다."""

    building_file = Path(building_path).expanduser().resolve()
    helipad_file = Path(helipad_path).expanduser().resolve()
    for label, path in (("건물", building_file), ("헬리포트", helipad_file)):
        if not path.is_file():
            raise ValueError(f"{label} Shapefile을 찾을 수 없습니다: {path}")

    try:
        buildings = gpd.read_file(building_file)
        helipads = gpd.read_file(helipad_file)
    except Exception as error:
        raise ValueError(f"입력 Shapefile을 읽을 수 없습니다: {error}") from error

    if buildings.crs is None:
        raise ValueError(f"건물 Shapefile에 CRS가 없습니다: {building_file}")
    if helipads.crs is None:
        raise ValueError(f"헬리포트 Shapefile에 CRS가 없습니다: {helipad_file}")
    if buildings.empty:
        raise ValueError(f"건물 Shapefile에 Feature가 없습니다: {building_file}")
    if helipads.empty:
        raise ValueError(f"헬리포트 Shapefile에 Feature가 없습니다: {helipad_file}")

    required = {id_field, longitude_field, latitude_field, agl_field}
    missing = sorted(required - set(helipads.columns))
    if missing:
        raise ValueError(
            f"헬리포트 Shapefile에 필수 필드가 없습니다: {', '.join(missing)}"
        )
    if "A16" not in buildings.columns or "A17" not in buildings.columns:
        raise ValueError("건물 Shapefile에 host 검증용 A16/A17 필드가 없습니다.")

    for label, frame in (("건물", buildings), ("헬리포트", helipads)):
        invalid_geometry = frame.geometry.isna() | frame.geometry.is_empty
        if invalid_geometry.any():
            indices = frame.index[invalid_geometry].tolist()
            raise ValueError(f"{label} Geometry가 비어 있는 행이 있습니다: {indices}")
        invalid_shape = ~frame.geometry.is_valid
        if invalid_shape.any():
            indices = frame.index[invalid_shape].tolist()
            raise ValueError(f"{label} Geometry가 유효하지 않은 행이 있습니다: {indices}")

    return (
        buildings.to_crs("EPSG:4326").reset_index(drop=True),
        helipads.to_crs("EPSG:4326").reset_index(drop=True),
    )


def prepare_helipad_records(
    buildings: gpd.GeoDataFrame,
    helipads: gpd.GeoDataFrame,
    id_field: str = SCRIPT_CONFIG.id_field,
    longitude_field: str = SCRIPT_CONFIG.longitude_field,
    latitude_field: str = SCRIPT_CONFIG.latitude_field,
    agl_field: str = SCRIPT_CONFIG.agl_field,
    host_match_crs: str = SCRIPT_CONFIG.host_match_crs,
) -> list[HelipadRecord]:
    """헬리포트 필드를 검증하고 최근접 host 건물 및 유효 AGL을 연결한다.

    A16이 결측 또는 0 이하이면 동일 좌표·동일 Geometry를 가진 다른 헬리포트
    레코드의 양수 A16을 사용한다. 건물 A16은 검증 정보로만 기록하며 헬리포트
    AGL을 대체하지 않는다.
    """

    prepared = helipads.copy().reset_index(drop=False, names="__source_index")
    prepared["__id"] = prepared[id_field].map(_identifier_text)
    if prepared["__id"].str.casefold().duplicated(keep=False).any():
        duplicates = prepared.loc[
            prepared["__id"].str.casefold().duplicated(keep=False), "__id"
        ].tolist()
        raise ValueError(f"중복된 헬리포트 id가 있습니다: {duplicates}")

    prepared["__lon"] = pd.to_numeric(prepared[longitude_field], errors="coerce")
    prepared["__lat"] = pd.to_numeric(prepared[latitude_field], errors="coerce")
    prepared["__agl"] = pd.to_numeric(prepared[agl_field], errors="coerce")
    invalid_coordinates = (
        ~np.isfinite(prepared["__lon"])
        | ~np.isfinite(prepared["__lat"])
        | ~prepared["__lon"].between(-180.0, 180.0)
        | ~prepared["__lat"].between(-90.0, 90.0)
    )
    if invalid_coordinates.any():
        bad_ids = prepared.loc[invalid_coordinates, "__id"].tolist()
        raise ValueError(f"경도/위도가 유효하지 않은 헬리포트 id: {bad_ids}")

    coordinate_points = gpd.GeoSeries(
        gpd.points_from_xy(prepared["__lon"], prepared["__lat"]),
        index=prepared.index,
        crs="EPSG:4326",
    )
    outside_geometry = [
        not geometry.covers(point)
        for geometry, point in zip(prepared.geometry, coordinate_points)
    ]
    if any(outside_geometry):
        bad_ids = prepared.loc[outside_geometry, "__id"].tolist()
        raise ValueError(
            "x_long/y_lati 좌표가 헬리포트 Geometry 밖에 있는 id: "
            f"{bad_ids}"
        )

    prepared["__peer_key"] = [
        (round(lon, 12), round(lat, 12), geometry.wkb_hex)
        for lon, lat, geometry in zip(
            prepared["__lon"], prepared["__lat"], prepared.geometry
        )
    ]
    peer_values: dict[object, list[float]] = {}
    peer_indices: dict[object, tuple[int, ...]] = {}
    for key, group in prepared.groupby("__peer_key", sort=False):
        values = sorted(
            {
                float(value)
                for value in group["__agl"]
                if np.isfinite(value) and float(value) > 0.0
            }
        )
        peer_values[key] = values
        peer_indices[key] = tuple(int(value) for value in group["__source_index"])

    agl_values: list[float] = []
    agl_sources: list[str] = []
    unresolved: list[str] = []
    for row in prepared.to_dict("records"):
        if np.isfinite(row["__agl"]) and float(row["__agl"]) > 0.0:
            agl_values.append(float(row["__agl"]))
            agl_sources.append(agl_field)
            continue
        candidates = peer_values[row["__peer_key"]]
        if len(candidates) == 1:
            agl_values.append(candidates[0])
            agl_sources.append(f"peer_{agl_field}")
        else:
            unresolved.append(row["__id"])
            agl_values.append(math.nan)
            agl_sources.append("unresolved")
    if unresolved:
        raise ValueError(
            f"{agl_field}을 결정할 수 없는 헬리포트 id: {unresolved}"
        )
    prepared["__resolved_agl"] = agl_values
    prepared["__agl_source"] = agl_sources

    building_match = buildings[["A16", "A17", "geometry"]].copy()
    building_match["__host_index"] = buildings.index.to_numpy()
    left = prepared[["__source_index", "geometry"]].to_crs(host_match_crs)
    right = building_match.to_crs(host_match_crs)
    nearest = gpd.sjoin_nearest(
        left,
        right,
        how="left",
        distance_col="__host_distance",
    )
    nearest = (
        nearest.sort_values(["__source_index", "__host_distance", "__host_index"])
        .drop_duplicates("__source_index", keep="first")
        .set_index("__source_index")
    )

    records: list[HelipadRecord] = []
    used_folders: dict[str, str] = {}
    for row in prepared.to_dict("records"):
        folder = _safe_folder_name(row["__id"])
        collision_key = folder.casefold()
        if collision_key in used_folders:
            raise ValueError(
                f"id 폴더명이 충돌합니다: "
                f"{used_folders[collision_key]!r}, {row['__id']!r}"
            )
        used_folders[collision_key] = row["__id"]
        host = nearest.loc[row["__source_index"]]
        records.append(
            HelipadRecord(
                identifier=row["__id"],
                folder_name=folder,
                source_index=int(row["__source_index"]),
                latitude=float(row["__lat"]),
                longitude=float(row["__lon"]),
                agl_m=float(row["__resolved_agl"]),
                agl_source=row["__agl_source"],
                host_index=int(host["__host_index"]),
                host_distance_m=float(host["__host_distance"]),
                host_a16_m=float(host["A16"]),
                host_a17_m=float(host["A17"]),
                excluded_helipad_indices=peer_indices[row["__peer_key"]],
            )
        )
    return records


def sample_dem_at_reference_point(
    latitude: float,
    longitude: float,
    dem_path: str | Path = SCRIPT_CONFIG.dem_path,
    band: int = SCRIPT_CONFIG.dem_band,
) -> DemSample:
    """WGS84 기준점이 포함된 DEM 셀의 지표고를 읽는다.

    기준점은 DEM CRS로 변환되며, 해당 위치를 포함하는 셀의 값을 사용한다.
    DEM 범위 밖 또는 NoData 셀이면 임의의 대체값을 사용하지 않고 실패한다.
    """

    lat = _finite_number("latitude", latitude)
    lon = _finite_number("longitude", longitude)
    if not -90.0 <= lat <= 90.0:
        raise ValueError(f"위도는 -90~90도여야 합니다: {lat}")
    if not -180.0 <= lon <= 180.0:
        raise ValueError(f"경도는 -180~180도여야 합니다: {lon}")
    if not isinstance(band, int) or isinstance(band, bool) or band < 1:
        raise ValueError(f"DEM 밴드는 1 이상의 정수여야 합니다: {band!r}")

    path = Path(dem_path).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"DEM 파일을 찾을 수 없습니다: {path}")

    try:
        with rasterio.open(path) as dataset:
            if dataset.crs is None:
                raise ValueError(f"DEM에 CRS가 없습니다: {path}")
            if band > dataset.count:
                raise ValueError(
                    f"DEM 밴드 {band}가 없습니다. 사용 가능한 밴드: 1~{dataset.count}"
                )

            transformer = Transformer.from_crs(
                "EPSG:4326", dataset.crs, always_xy=True
            )
            x, y = transformer.transform(lon, lat)
            row, column = dataset.index(x, y)
            if not (0 <= row < dataset.height and 0 <= column < dataset.width):
                raise ValueError(
                    "기준점이 DEM 범위 밖입니다: "
                    f"lat={lat}, lon={lon}, DEM CRS 좌표=({x:.3f}, {y:.3f}), "
                    f"DEM bounds={dataset.bounds}"
                )

            cell = dataset.read(
                band,
                window=Window(column, row, 1, 1),
                masked=True,
            )
            if bool(np.ma.getmaskarray(cell)[0, 0]):
                raise ValueError(
                    "기준점이 DEM NoData 셀에 있습니다: "
                    f"lat={lat}, lon={lon}, row={row}, column={column}"
                )
            elevation = float(cell[0, 0])
            if not math.isfinite(elevation):
                raise ValueError(
                    f"기준점의 DEM 고도가 유한하지 않습니다: {elevation!r}"
                )

            return DemSample(
                elevation_m=elevation,
                path=path,
                crs=dataset.crs.to_string(),
                band=band,
                row=row,
                column=column,
                x=float(x),
                y=float(y),
            )
    except rasterio.errors.RasterioError as error:
        raise ValueError(f"DEM을 읽을 수 없습니다: {path}: {error}") from error


def _validate_parameters(
    latitude: float,
    longitude: float,
    apex_agl_m: float,
    ground_z_m: float,
    horizontal_height_m: float,
    base_width_m: float,
    slope_deg: float,
    rotation_deg: float,
) -> tuple[float, float, float, float, float, float, float, float]:
    values = (
        _finite_number("latitude", latitude),
        _finite_number("longitude", longitude),
        _finite_number("apex_agl_m", apex_agl_m),
        _finite_number("ground_z_m", ground_z_m),
        _finite_number("horizontal_height_m", horizontal_height_m),
        _finite_number("base_width_m", base_width_m),
        _finite_number("slope_deg", slope_deg),
        _finite_number("rotation_deg", rotation_deg),
    )
    lat, lon, agl, ground, height, width, slope, rotation = values
    if not -90.0 <= lat <= 90.0:
        raise ValueError(f"위도는 -90~90도여야 합니다: {lat}")
    if not -180.0 <= lon <= 180.0:
        raise ValueError(f"경도는 -180~180도여야 합니다: {lon}")
    if agl < 0.0:
        raise ValueError(f"AGL은 0 m 이상이어야 합니다: {agl}")
    if height <= 0.0:
        raise ValueError(f"삼각형 높이는 0 m보다 커야 합니다: {height}")
    if width <= 0.0:
        raise ValueError(f"밑변 폭은 0 m보다 커야 합니다: {width}")
    if not 0.0 <= slope < 90.0:
        raise ValueError(f"OLS 경사각은 0도 이상 90도 미만이어야 합니다: {slope}")
    return lat, lon, agl, ground, height, width, slope, rotation


def create_ols_triangle(
    latitude: float,
    longitude: float,
    apex_agl_m: float,
    ground_z_m: float,
    rotation_deg: float = SCRIPT_CONFIG.start_rotation_deg,
    horizontal_height_m: float = SCRIPT_CONFIG.triangle_height_m,
    base_width_m: float = SCRIPT_CONFIG.base_width_m,
    slope_deg: float = SCRIPT_CONFIG.ols_slope_deg,
) -> OLSTriangle:
    """기준 꼭짓점을 중심으로 회전 가능한 삼각형 OLS를 생성한다.

    ``rotation_deg``는 꼭짓점에서 밑변 중심을 보는 초기 방위각이다. WGS84
    타원체의 측지선을 따라 ``horizontal_height_m``만큼 이동해 밑변 중심을
    구한 다음, 그 지점의 중심선에 직교하도록 밑변 양 끝을 배치한다.

    기준 꼭짓점의 표면고는 DEM 지표고와 입력 AGL을 더해 구한다. OLS는
    기준점에서 바깥쪽으로 상승하며 밑변의 표면고는 다음과 같다.

    ``base_z = ground_z + apex_agl + horizontal_height * tan(slope_deg)``
    """

    lat, lon, agl, ground, height, width, slope, rotation = _validate_parameters(
        latitude,
        longitude,
        apex_agl_m,
        ground_z_m,
        horizontal_height_m,
        base_width_m,
        slope_deg,
        rotation_deg,
    )

    bearing = rotation % 360.0
    base_center_lon, base_center_lat, back_azimuth = WGS84_GEOD.fwd(
        lon, lat, bearing, height
    )
    # fwd가 반환하는 것은 도착점에서 출발점을 향하는 역방위각이다.
    # 180도를 더해 도착점의 진행방향을 얻고, 그 좌우 90도에 밑변을 놓는다.
    terminal_azimuth = (back_azimuth + 180.0) % 360.0
    half_width = width / 2.0
    right_lon, right_lat, _ = WGS84_GEOD.fwd(
        base_center_lon, base_center_lat, terminal_azimuth + 90.0, half_width
    )
    left_lon, left_lat, _ = WGS84_GEOD.fwd(
        base_center_lon, base_center_lat, terminal_azimuth - 90.0, half_width
    )

    vertical_rise = height * math.tan(math.radians(slope))
    apex_z = ground + agl
    base_z = apex_z + vertical_rise

    return OLSTriangle(
        apex=Vertex("apex", lat, lon, apex_z),
        base_right=Vertex("base_right", right_lat, right_lon, base_z),
        base_left=Vertex("base_left", left_lat, left_lon, base_z),
        base_center=Vertex(
            "base_center", base_center_lat, base_center_lon, base_z
        ),
        rotation_deg=rotation,
        horizontal_height_m=height,
        base_width_m=width,
        slope_deg=slope,
        apex_agl_m=agl,
        ground_z_m=ground,
    )


def rotation_angles(
    start_deg: float = SCRIPT_CONFIG.start_rotation_deg,
    end_deg: float = SCRIPT_CONFIG.end_rotation_deg,
    step_deg: float = SCRIPT_CONFIG.rotation_step_deg,
    include_end: bool = SCRIPT_CONFIG.include_end_angle,
) -> tuple[float, ...]:
    """시작각부터 종료각까지 일정 간격의 회전각을 만든다.

    기본값은 0도 이상 360도 미만이므로 0~359도의 360개 값이다. 360도는
    0도와 같은 형상이지만 ``include_end=True``이면 360도 레코드도 추가한다.
    """

    start = _finite_number("start_deg", start_deg)
    end = _finite_number("end_deg", end_deg)
    step = _finite_number("step_deg", step_deg)
    if end < start:
        raise ValueError(f"종료각은 시작각 이상이어야 합니다: {start} > {end}")
    if step <= 0.0:
        raise ValueError(f"회전 간격은 0도보다 커야 합니다: {step}")

    tolerance = max(1.0, abs(start), abs(end)) * 1e-12
    angles: list[float] = []
    index = 0
    while True:
        angle = start + index * step
        inside = angle <= end + tolerance if include_end else angle < end - tolerance
        if not inside:
            break
        angles.append(
            end
            if include_end and math.isclose(angle, end, abs_tol=tolerance)
            else angle
        )
        index += 1
        if index > 1_000_000:
            raise ValueError("생성할 회전각이 1,000,000개를 초과합니다.")

    if not angles:
        raise ValueError("지정한 범위에서 생성되는 회전각이 없습니다.")
    return tuple(angles)


def create_rotated_ols_triangles(
    latitude: float,
    longitude: float,
    apex_agl_m: float,
    ground_z_m: float,
    angles_deg: Sequence[float],
    horizontal_height_m: float = SCRIPT_CONFIG.triangle_height_m,
    base_width_m: float = SCRIPT_CONFIG.base_width_m,
    slope_deg: float = SCRIPT_CONFIG.ols_slope_deg,
) -> list[OLSTriangle]:
    """동일한 기준점을 중심으로 지정된 각도별 OLS Polygon을 만든다."""

    if len(angles_deg) == 0:
        raise ValueError("angles_deg에 회전각을 하나 이상 지정해야 합니다.")
    return [
        create_ols_triangle(
            latitude=latitude,
            longitude=longitude,
            apex_agl_m=apex_agl_m,
            ground_z_m=ground_z_m,
            rotation_deg=angle,
            horizontal_height_m=horizontal_height_m,
            base_width_m=base_width_m,
            slope_deg=slope_deg,
        )
        for angle in angles_deg
    ]


def _sample_dem_at_geometry_points(
    frame: gpd.GeoDataFrame,
    dem_path: str | Path,
    band: int,
) -> np.ndarray:
    """각 Geometry 내부 대표점에서 DEM 지표고를 읽는다.

    대표점은 건물 footprint 내부가 보장되므로 centroid가 오목한 Polygon 밖으로
    벗어나는 문제를 피한다. NoData 또는 범위 밖인 위치는 ``nan``으로 남긴다.
    """

    path = Path(dem_path).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"DEM 파일을 찾을 수 없습니다: {path}")
    if frame.crs is None:
        raise ValueError("DEM 표본을 읽을 입력 레이어에 CRS가 없습니다.")
    try:
        with rasterio.open(path) as dataset:
            if dataset.crs is None:
                raise ValueError(f"DEM에 CRS가 없습니다: {path}")
            if band < 1 or band > dataset.count:
                raise ValueError(
                    f"DEM 밴드 {band}가 없습니다. 사용 가능한 밴드: 1~{dataset.count}"
                )
            points = frame.geometry.representative_point()
            transformer = Transformer.from_crs(
                frame.crs, dataset.crs, always_xy=True
            )
            coordinates = [
                transformer.transform(point.x, point.y) for point in points
            ]
            elevations = np.full(len(frame), np.nan, dtype=float)
            for index, value in enumerate(
                dataset.sample(coordinates, indexes=band, masked=True)
            ):
                cell = value[0]
                if not np.ma.is_masked(cell) and math.isfinite(float(cell)):
                    elevations[index] = float(cell)
            return elevations
    except rasterio.errors.RasterioError as error:
        raise ValueError(f"DEM을 읽을 수 없습니다: {path}: {error}") from error


def prepare_collision_layers(
    buildings: gpd.GeoDataFrame,
    helipads: gpd.GeoDataFrame,
    helipad_records: Sequence[HelipadRecord],
    helipad_dem_samples: dict[str, DemSample],
    dem_path: str | Path,
    dem_band: int,
    analysis_crs: str,
) -> CollisionLayers:
    """XY 후보 검색용 도형과 ``DEM + A16`` 장애물 상단고를 준비한다.

    건물 ``A16``이 0·결측이거나 DEM을 읽을 수 없는 경우 높이 미상(``nan``)으로
    둔다. 높이 미상 장애물은 XY footprint가 겹치더라도 충돌로 단정하지 않는다.
    헬리포트 높이는 peer 보완까지 끝난 ``HelipadRecord.agl_m``을 사용한다.
    """

    if "A16" not in buildings.columns:
        raise ValueError("건물 Shapefile에 높이 필드 A16이 없습니다.")

    building_ground = _sample_dem_at_geometry_points(
        buildings, dem_path=dem_path, band=dem_band
    )
    building_agl = pd.to_numeric(buildings["A16"], errors="coerce").to_numpy(
        dtype=float
    )
    valid_building_height = (
        np.isfinite(building_ground)
        & np.isfinite(building_agl)
        & (building_agl > 0.0)
    )
    building_top = np.full(len(buildings), np.nan, dtype=float)
    building_top[valid_building_height] = (
        building_ground[valid_building_height]
        + building_agl[valid_building_height]
    )

    helipad_top = np.full(len(helipads), np.nan, dtype=float)
    for record in helipad_records:
        sample = helipad_dem_samples[record.identifier]
        helipad_top[record.source_index] = sample.elevation_m + record.agl_m

    buildings_metric = buildings.copy()
    buildings_metric.geometry = buildings_metric.geometry.map(shapely.force_2d)
    buildings_metric = buildings_metric.to_crs(analysis_crs)
    helipads_metric = helipads.copy()
    helipads_metric.geometry = helipads_metric.geometry.map(shapely.force_2d)
    helipads_metric = helipads_metric.to_crs(analysis_crs)
    return CollisionLayers(
        analysis_crs=analysis_crs,
        buildings=buildings_metric,
        helipads=helipads_metric,
        building_top_z_m=building_top,
        helipad_top_z_m=helipad_top,
    )


def _obstacle_penetrates_triangle(
    triangle: OLSTriangle,
    footprint: Polygon,
    apex_xy: np.ndarray,
    direction: np.ndarray,
    obstacle_geometry: Polygon,
    obstacle_top_z_m: float,
    tolerance_m: float = 1e-6,
) -> bool:
    """장애물 상단이 footprint 중첩 위치의 OLS 평면에 닿는지 판정한다."""

    if not math.isfinite(obstacle_top_z_m):
        return False
    overlap = footprint.intersection(obstacle_geometry)
    if overlap.is_empty:
        return False
    coordinates = shapely.get_coordinates(overlap)
    if len(coordinates) == 0:
        point = overlap.representative_point()
        coordinates = np.array([[point.x, point.y]], dtype=float)
    along_distances = np.clip(
        (coordinates[:, :2] - apex_xy) @ direction,
        0.0,
        triangle.horizontal_height_m,
    )
    rise_per_metre = triangle.vertical_rise_m / triangle.horizontal_height_m
    minimum_surface_z = (
        triangle.apex_z_m + float(np.min(along_distances)) * rise_per_metre
    )
    return obstacle_top_z_m >= minimum_surface_z - tolerance_m


def validate_candidate_triangles(
    triangles: Sequence[OLSTriangle],
    collision_layers: CollisionLayers,
    helipad: HelipadRecord,
) -> tuple[list[OLSTriangle], tuple[float, ...]]:
    """장애물 상단고와 OLS 표면고를 비교해 안전/충돌 후보로 분리한다.

    XY intersection은 장애물 후보 검색에만 사용하며 그 자체로 충돌이 아니다.
    후보 Polygon끼리는 비교하지 않는다. 최근접 host 건물과 자기 헬리포트의
    동일 Geometry 레코드는 검사 대상에서 제외한다. 높이 미상 장애물 역시 XY
    중첩만으로 방위각을 제거하지 않는다.
    """

    excluded_helipads = set(helipad.excluded_helipad_indices)
    transformer = Transformer.from_crs(
        "EPSG:4326", collision_layers.analysis_crs, always_xy=True
    )
    safe: list[OLSTriangle] = []
    conflicts: list[float] = []
    for triangle in triangles:
        projected_vertices = [
            transformer.transform(vertex.longitude, vertex.latitude)
            for vertex in triangle.vertices
        ]
        footprint = Polygon(projected_vertices)
        apex_xy = np.asarray(projected_vertices[0], dtype=float)
        base_center_xy = np.asarray(
            transformer.transform(
                triangle.base_center.longitude, triangle.base_center.latitude
            ),
            dtype=float,
        )
        direction = base_center_xy - apex_xy
        direction /= np.linalg.norm(direction)
        building_hits = set(
            collision_layers.buildings.sindex.query(
                footprint, predicate="intersects"
            ).tolist()
        )
        building_hits.discard(helipad.host_index)
        helipad_hits = set(
            collision_layers.helipads.sindex.query(
                footprint, predicate="intersects"
            ).tolist()
        )
        helipad_hits.difference_update(excluded_helipads)

        building_collision = any(
            _obstacle_penetrates_triangle(
                triangle=triangle,
                footprint=footprint,
                apex_xy=apex_xy,
                direction=direction,
                obstacle_geometry=collision_layers.buildings.geometry.iloc[index],
                obstacle_top_z_m=collision_layers.building_top_z_m[index],
            )
            for index in building_hits
        )
        helipad_collision = any(
            _obstacle_penetrates_triangle(
                triangle=triangle,
                footprint=footprint,
                apex_xy=apex_xy,
                direction=direction,
                obstacle_geometry=collision_layers.helipads.geometry.iloc[index],
                obstacle_top_z_m=collision_layers.helipad_top_z_m[index],
            )
            for index in helipad_hits
        )
        if building_collision or helipad_collision:
            conflicts.append(triangle.rotation_deg)
        else:
            safe.append(triangle)
    return safe, tuple(conflicts)


def format_angle_ranges(
    angles_deg: Sequence[float], step_deg: float = 1.0
) -> str:
    """연속 충돌각을 ``[0,10] & [20,30]`` 형식으로 압축한다."""

    if len(angles_deg) == 0:
        return ""
    step = _finite_number("step_deg", step_deg)
    if step <= 0.0:
        raise ValueError(f"각도 간격은 0보다 커야 합니다: {step}")
    angles = sorted(set(float(angle) for angle in angles_deg))
    ranges: list[tuple[float, float]] = []
    start = previous = angles[0]
    tolerance = max(1.0, abs(step)) * 1e-9
    for angle in angles[1:]:
        if math.isclose(angle - previous, step, abs_tol=tolerance):
            previous = angle
            continue
        ranges.append((start, previous))
        start = previous = angle
    ranges.append((start, previous))
    return " & ".join(f"[{start:g},{end:g}]" for start, end in ranges)


def _triangle_record(
    triangle: OLSTriangle,
    dem_sample: DemSample,
    helipad: HelipadRecord,
) -> dict[str, float | int | str | Polygon]:
    """10자 이하 DBF 필드명으로 Shapefile 레코드를 만든다."""

    return {
        "heli_id": helipad.identifier,
        "rot_deg": triangle.rotation_deg,
        "src_fid": helipad.source_index,
        "agl_src": helipad.agl_source,
        "apex_lat": triangle.apex.latitude,
        "apex_lon": triangle.apex.longitude,
        "apex_agl": triangle.apex_agl_m,
        "ground_z": triangle.ground_z_m,
        "apex_z": triangle.apex_z_m,
        "right_lat": triangle.base_right.latitude,
        "right_lon": triangle.base_right.longitude,
        "left_lat": triangle.base_left.latitude,
        "left_lon": triangle.base_left.longitude,
        "base_z": triangle.base_z_m,
        "rise_m": triangle.vertical_rise_m,
        "height_m": triangle.horizontal_height_m,
        "width_m": triangle.base_width_m,
        "slope_deg": triangle.slope_deg,
        "h_crs": "EPSG:4326",
        "z_ref": "DEM+AGL",
        "dem_file": dem_sample.path.name,
        "dem_crs": dem_sample.crs,
        "dem_band": dem_sample.band,
        "dem_row": dem_sample.row,
        "dem_col": dem_sample.column,
        "host_fid": helipad.host_index,
        "host_dist": helipad.host_distance_m,
        "host_a16": helipad.host_a16_m,
        "host_a17": helipad.host_a17_m,
        "geometry": triangle.polygon,
    }


def save_shapefile(
    triangles: Sequence[OLSTriangle],
    output_path: str | Path,
    dem_sample: DemSample,
    helipad: HelipadRecord,
    template_triangle: OLSTriangle | None = None,
) -> Path:
    """충돌하지 않는 Polygon Z를 EPSG:4326 Shapefile로 저장한다.

    안전각이 하나도 없으면 ``template_triangle``의 속성 스키마를 이용해 행이
    0개인 Polygon Z Shapefile을 만든다.
    """

    output = Path(output_path).expanduser().resolve()
    if output.suffix.lower() != ".shp":
        raise ValueError(f"Shapefile 출력 경로는 .shp여야 합니다: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    if len(triangles) > 0:
        records = [
            _triangle_record(triangle, dem_sample, helipad)
            for triangle in triangles
        ]
        frame = gpd.GeoDataFrame(
            records,
            geometry="geometry",
            crs="EPSG:4326",
        )
        frame.to_file(
            output,
            driver="ESRI Shapefile",
            engine="pyogrio",
            encoding="UTF-8",
            index=False,
        )
    else:
        if template_triangle is None:
            raise ValueError("빈 Shapefile 스키마를 만들 template_triangle이 없습니다.")
        template = gpd.GeoDataFrame(
            [_triangle_record(template_triangle, dem_sample, helipad)],
            geometry="geometry",
            crs="EPSG:4326",
        )
        empty = template.iloc[0:0].copy()
        empty.to_file(
            output,
            driver="ESRI Shapefile",
            engine="pyogrio",
            encoding="UTF-8",
            index=False,
            geometry_type="Polygon Z",
        )
    return output


def remove_generated_shapefile(output_path: str | Path, allowed_root: str | Path) -> None:
    """안전각이 없을 때 이전 실행에서 남은 동일 이름 Shapefile만 제거한다."""

    output = Path(output_path).expanduser().resolve()
    root = Path(allowed_root).expanduser().resolve()
    if not output.is_relative_to(root):
        raise ValueError(f"출력 루트 밖의 파일은 제거할 수 없습니다: {output}")
    sidecars = (
        ".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix", ".fix", ".sbn", ".sbx"
    )
    for suffix in sidecars:
        candidate = output.with_suffix(suffix)
        if candidate.is_file():
            candidate.unlink()
    metadata = Path(f"{output}.xml")
    if metadata.is_file():
        metadata.unlink()
    try:
        output.parent.rmdir()
    except OSError:
        # 생성물 외의 파일이 있으면 폴더는 보존한다.
        pass


def assert_shapefile_unlocked(output_path: str | Path, allowed_root: str | Path) -> None:
    """기존 출력 세트가 QGIS 등 다른 프로세스에 잠겨 있지 않은지 확인한다."""

    output = Path(output_path).expanduser().resolve()
    root = Path(allowed_root).expanduser().resolve()
    if not output.is_relative_to(root):
        raise ValueError(f"출력 루트 밖의 파일은 검사할 수 없습니다: {output}")
    sidecars = (
        ".shp", ".shx", ".dbf", ".prj", ".cpg", ".qix", ".fix", ".sbn", ".sbx"
    )
    for suffix in sidecars:
        candidate = output.with_suffix(suffix)
        if not candidate.is_file():
            continue
        try:
            os.rename(candidate, candidate)
        except PermissionError as error:
            raise ValueError(
                "기존 출력 파일이 다른 프로그램에서 사용 중입니다. "
                f"QGIS 등에서 레이어를 닫은 뒤 다시 실행하세요: {candidate}"
            ) from error


def generate_all_helipad_outputs(
    helipad_records: Sequence[HelipadRecord],
    buildings: gpd.GeoDataFrame,
    helipad_frame: gpd.GeoDataFrame,
    angles_deg: Sequence[float],
    output_root: str | Path = SCRIPT_CONFIG.output_root,
    dem_path: str | Path = SCRIPT_CONFIG.dem_path,
    dem_band: int = SCRIPT_CONFIG.dem_band,
    horizontal_height_m: float = SCRIPT_CONFIG.triangle_height_m,
    base_width_m: float = SCRIPT_CONFIG.base_width_m,
    slope_deg: float = SCRIPT_CONFIG.ols_slope_deg,
    angle_step_deg: float = SCRIPT_CONFIG.rotation_step_deg,
) -> list[GenerationResult]:
    """OLS 경사/id 폴더에 충돌각을 제외한 Polygon Shapefile을 만든다."""

    if len(helipad_records) == 0:
        raise ValueError("생성할 헬리포트가 없습니다.")
    if len(angles_deg) == 0:
        raise ValueError("생성할 회전각이 없습니다.")

    # 출력 폴더를 만들기 전에 모든 기준점의 DEM을 먼저 검증한다.
    samples = {
        helipad.identifier: sample_dem_at_reference_point(
            latitude=helipad.latitude,
            longitude=helipad.longitude,
            dem_path=dem_path,
            band=dem_band,
        )
        for helipad in helipad_records
    }
    collision_layers = prepare_collision_layers(
        buildings=buildings,
        helipads=helipad_frame,
        helipad_records=helipad_records,
        helipad_dem_samples=samples,
        dem_path=dem_path,
        dem_band=dem_band,
        analysis_crs=SCRIPT_CONFIG.host_match_crs,
    )

    root = Path(output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    slope_value = _finite_number("slope_deg", slope_deg)
    slope_directory = root / f"{slope_value:g}"
    slope_directory.mkdir(parents=True, exist_ok=True)
    validated: list[
        tuple[
            HelipadRecord,
            DemSample,
            list[OLSTriangle],
            list[OLSTriangle],
            tuple[float, ...],
            Path,
        ]
    ] = []
    for helipad in helipad_records:
        sample = samples[helipad.identifier]
        # 1) 해당 id의 모든 방위각 후보를 메모리에서 생성한다.
        candidate_triangles = create_rotated_ols_triangles(
            latitude=helipad.latitude,
            longitude=helipad.longitude,
            apex_agl_m=helipad.agl_m,
            ground_z_m=sample.elevation_m,
            angles_deg=angles_deg,
            horizontal_height_m=horizontal_height_m,
            base_width_m=base_width_m,
            slope_deg=slope_value,
        )
        # 2) 저장 전에 후보를 검증하고 3차원 높이 충돌 방위각을 제거한다.
        # XY intersection은 장애물 후보 검색에만 쓰며 그 자체로 제거하지 않는다.
        # candidate_triangles 상호 간에는 intersection 검사를 수행하지 않는다.
        safe_triangles, collision_angles = validate_candidate_triangles(
            triangles=candidate_triangles,
            collision_layers=collision_layers,
            helipad=helipad,
        )
        id_directory = slope_directory / helipad.folder_name
        output = id_directory / f"Polygon_{helipad.folder_name}.shp"
        validated.append(
            (
                helipad,
                sample,
                candidate_triangles,
                safe_triangles,
                collision_angles,
                output,
            )
        )

    # 실제 파일을 하나라도 변경하기 전에 전체 기존 출력의 잠금 상태를 검사한다.
    for _, _, _, _, _, output in validated:
        assert_shapefile_unlocked(output, slope_directory)

    results: list[GenerationResult] = []
    for (
        helipad,
        sample,
        candidate_triangles,
        safe_triangles,
        collision_angles,
        output,
    ) in validated:
        # 3) 검증을 통과한 Polygon만 최종 출력 경로에 저장한다.
        final_output: Path | None
        if safe_triangles:
            final_output = save_shapefile(
                safe_triangles,
                output,
                sample,
                helipad,
                template_triangle=candidate_triangles[0],
            )
        else:
            remove_generated_shapefile(output, slope_directory)
            final_output = None
        results.append(
            GenerationResult(
                helipad=helipad,
                dem_sample=sample,
                output_path=final_output,
                polygon_count=len(safe_triangles),
                collision_angles=collision_angles,
            )
        )

    summary = pd.DataFrame(
        [
            {
                "id": result.helipad.identifier,
                "source_index": result.helipad.source_index,
                "longitude": result.helipad.longitude,
                "latitude": result.helipad.latitude,
                "agl_m": result.helipad.agl_m,
                "agl_source": result.helipad.agl_source,
                "ground_z_m": result.dem_sample.elevation_m,
                "apex_z_m": result.dem_sample.elevation_m + result.helipad.agl_m,
                "host_index": result.helipad.host_index,
                "host_distance_m": result.helipad.host_distance_m,
                "host_a16_m": result.helipad.host_a16_m,
                "host_a17_m": result.helipad.host_a17_m,
                "ols_slope_deg": slope_value,
                "polygon_count": result.polygon_count,
                "collision_count": len(result.collision_angles),
                "collision": format_angle_ranges(
                    result.collision_angles, step_deg=angle_step_deg
                ),
                "output_path": str(result.output_path) if result.output_path else "",
            }
            for result in results
        ]
    )
    summary_path = slope_directory / "generation_summary.csv"
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    legacy_summary = root / "generation_summary.csv"
    if legacy_summary.is_file() and legacy_summary != summary_path:
        legacy_summary.unlink()
    return results


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "송파구 옥상헬리포트의 id/x_long/y_lati/A16을 읽어 OLS 경사/id별 "
            "폴더에 건물·다른 헬리포트가 OLS 표면고에 닿거나 관통하지 않는 "
            "삼각형 Polygon Z만 저장합니다."
        )
    )
    parser.add_argument(
        "--buildings",
        type=Path,
        default=SCRIPT_CONFIG.building_path,
        help=f"건물 Shapefile 경로(기본: {SCRIPT_CONFIG.building_path})",
    )
    parser.add_argument(
        "--helipads",
        type=Path,
        default=SCRIPT_CONFIG.helipad_path,
        help=f"옥상 헬리포트 Shapefile 경로(기본: {SCRIPT_CONFIG.helipad_path})",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=SCRIPT_CONFIG.output_root,
        help=f"id별 출력 폴더의 상위 경로(기본: {SCRIPT_CONFIG.output_root})",
    )
    parser.add_argument(
        "--dem",
        type=Path,
        default=SCRIPT_CONFIG.dem_path,
        help=f"기준점 지표고를 읽을 DEM 경로(기본: {SCRIPT_CONFIG.dem_path})",
    )
    parser.add_argument(
        "--dem-band",
        type=int,
        default=SCRIPT_CONFIG.dem_band,
        help=f"DEM 고도 밴드 번호(기본: {SCRIPT_CONFIG.dem_band})",
    )
    parser.add_argument(
        "--id-field",
        default=SCRIPT_CONFIG.id_field,
        help=f"헬리포트 id 필드(기본: {SCRIPT_CONFIG.id_field})",
    )
    parser.add_argument(
        "--longitude-field",
        default=SCRIPT_CONFIG.longitude_field,
        help=f"헬리포트 경도 필드(기본: {SCRIPT_CONFIG.longitude_field})",
    )
    parser.add_argument(
        "--latitude-field",
        default=SCRIPT_CONFIG.latitude_field,
        help=f"헬리포트 위도 필드(기본: {SCRIPT_CONFIG.latitude_field})",
    )
    parser.add_argument(
        "--agl-field",
        default=SCRIPT_CONFIG.agl_field,
        help=f"헬리포트 AGL 필드(기본: {SCRIPT_CONFIG.agl_field})",
    )
    parser.add_argument(
        "--host-match-crs",
        default=SCRIPT_CONFIG.host_match_crs,
        help=f"최근접 건물 거리 계산 CRS(기본: {SCRIPT_CONFIG.host_match_crs})",
    )
    parser.add_argument(
        "--start-angle",
        type=float,
        default=SCRIPT_CONFIG.start_rotation_deg,
        help=f"첫 회전각(도, 포함, 기본: {SCRIPT_CONFIG.start_rotation_deg:g})",
    )
    parser.add_argument(
        "--end-angle",
        type=float,
        default=SCRIPT_CONFIG.end_rotation_deg,
        help=(
            "마지막 회전각 경계(도, 기본적으로 미포함, 기본: "
            f"{SCRIPT_CONFIG.end_rotation_deg:g})"
        ),
    )
    parser.add_argument(
        "--angle-step",
        type=float,
        default=SCRIPT_CONFIG.rotation_step_deg,
        help=f"회전각 간격(도, 기본: {SCRIPT_CONFIG.rotation_step_deg:g})",
    )
    parser.add_argument(
        "--include-end-angle",
        action=argparse.BooleanOptionalAction,
        default=SCRIPT_CONFIG.include_end_angle,
        help="종료각 포함 여부(360도는 0도와 같은 중복 형상임)",
    )
    parser.add_argument(
        "--height",
        type=float,
        default=SCRIPT_CONFIG.triangle_height_m,
        help=f"삼각형의 수평 중심선 길이(m, 기본: {SCRIPT_CONFIG.triangle_height_m:g})",
    )
    parser.add_argument(
        "--base-width",
        type=float,
        default=SCRIPT_CONFIG.base_width_m,
        help=f"밑변 폭(m, 기본: {SCRIPT_CONFIG.base_width_m:g})",
    )
    parser.add_argument(
        "--ols-angle",
        type=float,
        default=SCRIPT_CONFIG.ols_slope_deg,
        help=f"바깥쪽으로 상승하는 OLS 종단 경사각(도, 기본: {SCRIPT_CONFIG.ols_slope_deg:g})",
    )
    return parser


def _format_batch_summary(
    results: Sequence[GenerationResult],
    angles: Sequence[float],
    output_root: Path,
    slope_deg: float,
) -> str:
    fallback_ids = [
        result.helipad.identifier
        for result in results
        if result.helipad.agl_source.startswith("peer_")
    ]
    host_distances = [result.helipad.host_distance_m for result in results]
    slope_directory = output_root / f"{float(slope_deg):g}"
    safe_count = sum(result.polygon_count for result in results)
    collision_count = sum(len(result.collision_angles) for result in results)
    rows = [
        "헬리포트 id별 삼각형 OLS Shapefile 생성 완료",
        "  수평 좌표계: EPSG:4326 (WGS84)",
        "  Polygon Z: 기준점 DEM 지표고 + 헬리포트 A16 + OLS 상승량",
        f"  헬리포트/id 폴더 수: {len(results)}",
        f"  Shapefile이 생성된 id 수: {sum(result.output_path is not None for result in results)}",
        f"  id별 검사 방위각 수: {len(angles)}",
        f"  안전 Polygon 수: {safe_count}",
        f"  3차원 높이 충돌 제외 방위각 수: {collision_count}",
        f"  회전각 범위: {angles[0]:g}° ~ {angles[-1]:g}°",
        f"  peer A16 보완 id: {fallback_ids or '없음'}",
        f"  최근접 건물 거리 범위: {min(host_distances):.3f}~"
        f"{max(host_distances):.3f} m",
        f"  출력 루트: {output_root}",
        f"  OLS 경사 출력 폴더: {slope_directory}",
        f"  요약 CSV: {slope_directory / 'generation_summary.csv'}",
    ]
    return "\n".join(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        buildings, helipad_frame = load_source_layers(
            building_path=args.buildings,
            helipad_path=args.helipads,
            id_field=args.id_field,
            longitude_field=args.longitude_field,
            latitude_field=args.latitude_field,
            agl_field=args.agl_field,
        )
        helipads = prepare_helipad_records(
            buildings=buildings,
            helipads=helipad_frame,
            id_field=args.id_field,
            longitude_field=args.longitude_field,
            latitude_field=args.latitude_field,
            agl_field=args.agl_field,
            host_match_crs=args.host_match_crs,
        )
        angles = rotation_angles(
            start_deg=args.start_angle,
            end_deg=args.end_angle,
            step_deg=args.angle_step,
            include_end=args.include_end_angle,
        )
        results = generate_all_helipad_outputs(
            helipad_records=helipads,
            buildings=buildings,
            helipad_frame=helipad_frame,
            angles_deg=angles,
            output_root=args.output_root,
            dem_path=args.dem,
            dem_band=args.dem_band,
            horizontal_height_m=args.height,
            base_width_m=args.base_width,
            slope_deg=args.ols_angle,
            angle_step_deg=args.angle_step,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))

    output_root = Path(args.output_root).expanduser().resolve()
    print(_format_batch_summary(results, angles, output_root, args.ols_angle))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
