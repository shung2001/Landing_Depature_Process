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

처리 흐름(큰 그림)
-----------------
1. 건물/헬리포트 Shapefile을 읽고 좌표계와 필수 필드를 검사한다.
2. 각 헬리포트의 기준 좌표, AGL(지면으로부터의 높이), host 건물을 정리한다.
3. DEM에서 헬리포트와 건물 위치의 지표고를 읽어 절대 높이를 계산한다.
4. 헬리포트마다 지정한 방위각 수만큼 삼각형 OLS 후보를 메모리에 만든다.
5. 각 후보와 건물/다른 헬리포트의 수평 중첩 및 높이를 비교한다.
6. 충돌하지 않는 후보만 Shapefile로 저장하고 결과를 CSV로 요약한다.

처음 읽을 때 알아둘 용어
-----------------------
* OLS(Obstacle Limitation Surface): 장애물이 침범하면 안 되는 가상의 제한 표면
* AGL(Above Ground Level): 해당 위치의 지면으로부터 잰 높이
* DEM(Digital Elevation Model): 각 격자 셀에 지표고가 저장된 래스터 자료
* CRS(Coordinate Reference System): 좌표 숫자가 지구상의 위치를 뜻하는 방식
* Geometry/Feature: 각각 공간 도형/도형과 속성을 합친 한 개의 공간 객체
* footprint: 3차원 물체를 위에서 내려다본 2차원 바닥 도형
* host 건물: 해당 옥상 헬리포트가 올라가 있다고 간주하는 가장 가까운 건물

기본 형상
---------
* 수평 중심선 길이(삼각형 높이): 1,000 m
* 밑변 폭: 573 m
* OLS 종단 경사각: 8도(``SCRIPT_CONFIG.ols_slope_deg``에서 변경 가능)
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

# 타입 힌트에서 아직 정의되지 않은 클래스를 문자열로 감싸지 않고 쓸 수 있게 한다.
# 실행 결과에는 영향을 주지 않고, 아래 함수 선언을 더 읽기 쉽게 해 주는 설정이다.
from __future__ import annotations

# 표준 라이브러리: Python 설치 시 기본으로 제공된다.
import argparse
from dataclasses import dataclass
import math
import os
from pathlib import Path
import re
from typing import Sequence

# 외부 라이브러리: 표/공간 자료/좌표 변환/래스터 처리를 담당한다.
import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Geod, Transformer
import rasterio
from rasterio.windows import Window
import shapely
from shapely.geometry import Polygon


# __file__은 '현재 실행 중인 이 파일'의 경로이다. parents[2]는 파일 위치에서
# 두 단계 위로 올라가 프로젝트 루트를 찾는다. 절대 경로로 바꾸므로 실행한
# 터미널의 현재 폴더가 달라도 항상 같은 입력 자료를 찾을 수 있다.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SONGPA_DATA_DIR = PROJECT_ROOT / "자료" / "안전성" / "송파구"

# 사용자가 명령행 옵션을 주지 않았을 때 사용할 기본 입출력 경로이다.
DEFAULT_BUILDING_PATH = SONGPA_DATA_DIR / "송파구_건물_4326.shp"
DEFAULT_HELIPAD_PATH = SONGPA_DATA_DIR / "송파구_옥상헬리포트_4326.shp"
DEFAULT_DEM_PATH = SONGPA_DATA_DIR / "송파구_DEM.tif"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "자료" / "결과물" / "Polygon"


@dataclass(frozen=True)
class ScriptConfig:
    """명령행 인수 없이 실행할 때 사용할 값을 한곳에 모아 둔 설정 객체.

    ``@dataclass``는 아래 필드를 받는 생성자 등을 자동으로 만들어 준다.
    ``frozen=True``이므로 한 번 만든 설정 객체의 값을 실수로 바꿀 수 없다.
    경로/필드명/수치 설정을 묶어 두면 함수 곳곳에 같은 값을 반복하지 않아도 된다.
    """

    # 입력 자료와 결과가 저장될 위치
    building_path: Path = DEFAULT_BUILDING_PATH  # 건물 폴리곤 Shapefile
    helipad_path: Path = DEFAULT_HELIPAD_PATH  # 옥상 헬리포트 Shapefile
    dem_path: Path = DEFAULT_DEM_PATH  # 지표고가 들어 있는 DEM 래스터
    output_root: Path = DEFAULT_OUTPUT_ROOT  # 결과 폴더의 최상위 위치

    # 헬리포트 Shapefile에서 읽을 열(column) 이름
    id_field: str = "id"  # 각 헬리포트를 구분하는 고유값
    longitude_field: str = "x_long"  # 경도(X)
    latitude_field: str = "y_lati"  # 위도(Y)
    agl_field: str = "A16"  # 지면으로부터 헬리포트까지의 높이(m)

    # EPSG:4326은 각도 단위라 거리 계산에 부적합하다. 따라서 최근접 건물을
    # 찾을 때는 미터 단위의 한국 중부원점 좌표계 EPSG:5186으로 잠시 변환한다.
    host_match_crs: str = "EPSG:5186"

    # 회전각은 진북 0도에서 시작해 시계 방향으로 증가한다.
    start_rotation_deg: float = 0.0
    end_rotation_deg: float = 360.0
    rotation_step_deg: float = 1.0
    include_end_angle: bool = False  # False면 360도는 제외(0도와 같은 방향)

    # 삼각형 OLS의 기본 크기와 기울기
    triangle_height_m: float = 1_000.0  # 꼭짓점부터 밑변 중심까지의 수평거리
    base_width_m: float = 573.0  # 밑변 전체 폭
    ols_slope_deg: float = 8.0  # 바깥쪽으로 올라가는 종단 경사각
    dem_band: int = 1  # 여러 밴드 중 고도가 들어 있는 밴드 번호(1부터 시작)


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


# 위/경도 좌표에서 '방위각을 따라 몇 m 이동한 지점'을 계산하는 측지 객체이다.
# 단순한 평면 삼각법 대신 WGS84 타원체를 사용하므로 지구 곡률이 반영된다.
WGS84_GEOD = Geod(ellps="WGS84")


@dataclass(frozen=True)
class Vertex:
    """삼각형 꼭짓점 하나를 나타내는 불변(immutable) 데이터 묶음."""

    name: str  # 사람이 알아보기 위한 이름: apex, base_right 등
    latitude: float  # WGS84 위도(Y), 단위는 도(degree)
    longitude: float  # WGS84 경도(X), 단위는 도(degree)
    z_m: float  # DEM 기준의 절대 표면고, 단위는 m

    def xyz_coordinate(self) -> tuple[float, float, float]:
        """Shapefile Polygon Z 좌표 순서인 (경도, 위도, 표면고)을 반환한다."""

        # GIS 좌표는 보통 (X, Y, Z), 즉 (경도, 위도, 높이) 순서를 사용한다.
        return self.longitude, self.latitude, self.z_m


@dataclass(frozen=True)
class DemSample:
    """기준점 위치에서 읽은 DEM 값과 그 값의 출처를 함께 보관한다."""

    elevation_m: float  # DEM 셀에 저장된 지표고(m)
    path: Path  # 읽은 DEM 파일 경로
    crs: str  # DEM 자체의 좌표계
    band: int  # 값을 읽은 래스터 밴드
    row: int  # DEM 격자의 행 번호
    column: int  # DEM 격자의 열 번호
    x: float  # 기준점을 DEM 좌표계로 바꾼 X
    y: float  # 기준점을 DEM 좌표계로 바꾼 Y


@dataclass(frozen=True)
class HelipadRecord:
    """입력 한 행을 검증한 뒤 후속 계산에 필요한 값만 정리한 객체."""

    identifier: str  # 원본의 헬리포트 id를 문자열로 정규화한 값
    folder_name: str  # Windows 폴더명으로 안전하게 바꾼 id
    source_index: int  # 원본 GeoDataFrame에서의 행 번호
    latitude: float  # OLS 꼭짓점의 위도
    longitude: float  # OLS 꼭짓점의 경도
    agl_m: float  # 검증 또는 peer 보완을 마친 헬리포트 AGL(m)
    agl_source: str  # AGL이 원본인지 peer 보완값인지 나타내는 문자열
    host_index: int  # 가장 가까운 host 건물의 행 번호
    host_distance_m: float  # 헬리포트와 host 건물 사이의 최단거리(m)
    host_a16_m: float  # 결과 확인용 host 건물 A16
    host_a17_m: float  # 결과 확인용 host 건물 A17
    # 같은 위치/Geometry를 공유하는 행들은 '자기 자신'으로 보아 충돌 검사에서 제외
    excluded_helipad_indices: tuple[int, ...]


@dataclass(frozen=True)
class GenerationResult:
    """헬리포트 하나의 처리가 끝난 뒤 CSV와 화면 출력에 쓸 결과."""

    helipad: HelipadRecord  # 어떤 헬리포트의 결과인지 식별하는 정보
    dem_sample: DemSample  # 기준점에서 읽은 DEM 정보
    output_path: Path | None  # 안전각이 없으면 파일이 없으므로 None
    polygon_count: int  # 최종 저장된 안전 OLS 수
    collision_angles: tuple[float, ...]  # 충돌로 제외된 방위각 모음


@dataclass(frozen=True)
class OLSTriangle:
    """삼각형 도형과 그 도형을 재현·설명하는 계산 조건을 함께 보관한다."""

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

        # base_center는 밑변 가운데를 계산하기 위한 보조점이라 외곽 꼭짓점이 아니다.
        return self.apex, self.base_right, self.base_left

    @property
    def polygon(self) -> Polygon:
        """경도/위도/DEM 기반 표면고의 3차원 Shapely Polygon을 반환한다."""

        # 리스트 컴프리헨션은 꼭짓점 3개를 각각 (X, Y, Z) 튜플로 바꾼다.
        # Shapely의 Polygon은 마지막 점을 처음 점과 자동으로 연결해 닫힌 면을 만든다.
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
    """충돌 검사를 반복할 때 재사용할 도형과 높이 배열을 묶은 객체."""

    analysis_crs: str  # 교차/거리 계산에 사용하는 미터 단위 CRS
    buildings: gpd.GeoDataFrame  # 2차원으로 만들고 투영한 건물 도형
    helipads: gpd.GeoDataFrame  # 2차원으로 만들고 투영한 헬리포트 도형
    building_top_z_m: np.ndarray  # 각 건물 상단의 절대고(DEM + AGL)
    helipad_top_z_m: np.ndarray  # 각 헬리포트 상단의 절대고(DEM + AGL)


def _finite_number(name: str, value: float) -> float:
    """입력값을 ``float``로 바꾸고 NaN/무한대가 아닌지 확인한다.

    함수명 앞의 밑줄(``_``)은 이 파일 내부에서만 쓰는 보조 함수라는 관례이다.
    ``NaN``은 '숫자가 아님', ``inf``는 무한대를 뜻하며 거리 계산에 쓸 수 없다.
    """

    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name}은(는) 유한한 숫자여야 합니다: {value!r}")
    return value


def _identifier_text(value: object) -> str:
    """여러 자료형으로 들어온 id를 손실 없는 문자열로 정규화한다.

    예를 들어 정수 ``3``과 실수 ``3.0``은 모두 ``"3"``으로 만든다. 그대로
    ``str(3.0)``을 사용하면 폴더명이 불필요하게 ``3.0``이 되는 것을 막기 위함이다.
    """

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
    """id를 Windows에서도 안전한 단일 폴더명으로 검증·정리한다.

    Windows가 파일명에서 금지하는 문자와 제어 문자는 ``_``로 바꾼다. 또한
    ``CON``이나 ``COM1``처럼 운영체제가 장치명으로 예약한 이름도 피한다.
    """

    # 정규표현식 대괄호 안의 문자 중 하나라도 만나면 '_'로 치환한다.
    folder = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", identifier).strip().rstrip(".")
    if not folder or folder in {".", ".."}:
        raise ValueError(f"폴더명으로 사용할 수 없는 id입니다: {identifier!r}")
    # * 앞의 표현식은 COM1~COM9, LPT1~LPT9 문자열을 집합 안에 펼친다.
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
    """건물과 헬리포트 Shapefile을 읽고 기본 구조를 검증한다.

    Parameters
    ----------
    building_path, helipad_path
        입력 Shapefile 경로. ``str``과 ``Path``를 모두 받을 수 있다.
    id_field, longitude_field, latitude_field, agl_field
        헬리포트 파일에서 반드시 존재해야 하는 열 이름.

    Returns
    -------
    tuple[GeoDataFrame, GeoDataFrame]
        좌표계를 EPSG:4326으로 통일한 ``(건물, 헬리포트)`` 두 레이어.

    Notes
    -----
    함수 시작에서 입력을 엄격하게 확인하면 뒤쪽의 복잡한 공간 계산 도중에
    원인을 알기 어려운 오류가 나는 일을 줄일 수 있다. 이런 방식을 '빠르게
    실패하기(fail fast)'라고 한다.
    """

    # expanduser()는 경로의 '~'를 사용자 폴더로, resolve()는 절대 경로로 바꾼다.
    building_file = Path(building_path).expanduser().resolve()
    helipad_file = Path(helipad_path).expanduser().resolve()

    # 두 파일에 같은 검사를 적용하므로 (설명, 경로) 쌍을 순회한다.
    for label, path in (("건물", building_file), ("헬리포트", helipad_file)):
        if not path.is_file():
            raise ValueError(f"{label} Shapefile을 찾을 수 없습니다: {path}")

    # GeoPandas는 .shp뿐 아니라 함께 있는 .dbf/.shx/.prj 등도 같이 읽는다.
    # 라이브러리의 긴 원본 예외 대신 사용자가 이해하기 쉬운 메시지로 감싼다.
    try:
        buildings = gpd.read_file(building_file)
        helipads = gpd.read_file(helipad_file)
    except Exception as error:
        raise ValueError(f"입력 Shapefile을 읽을 수 없습니다: {error}") from error

    # CRS가 없으면 같은 숫자가 경도/위도인지 미터 좌표인지 알 수 없다.
    if buildings.crs is None:
        raise ValueError(f"건물 Shapefile에 CRS가 없습니다: {building_file}")
    if helipads.crs is None:
        raise ValueError(f"헬리포트 Shapefile에 CRS가 없습니다: {helipad_file}")
    if buildings.empty:
        raise ValueError(f"건물 Shapefile에 Feature가 없습니다: {building_file}")
    if helipads.empty:
        raise ValueError(f"헬리포트 Shapefile에 Feature가 없습니다: {helipad_file}")

    # set의 차집합을 이용해 '필수 열 - 실제 열 = 빠진 열'을 한 번에 구한다.
    required = {id_field, longitude_field, latitude_field, agl_field}
    missing = sorted(required - set(helipads.columns))
    if missing:
        raise ValueError(
            f"헬리포트 Shapefile에 필수 필드가 없습니다: {', '.join(missing)}"
        )
    if "A16" not in buildings.columns or "A17" not in buildings.columns:
        raise ValueError("건물 Shapefile에 host 검증용 A16/A17 필드가 없습니다.")

    # 결측/빈 도형과 self-intersection 등 유효하지 않은 도형을 미리 거른다.
    for label, frame in (("건물", buildings), ("헬리포트", helipads)):
        invalid_geometry = frame.geometry.isna() | frame.geometry.is_empty
        if invalid_geometry.any():
            indices = frame.index[invalid_geometry].tolist()
            raise ValueError(f"{label} Geometry가 비어 있는 행이 있습니다: {indices}")
        invalid_shape = ~frame.geometry.is_valid
        if invalid_shape.any():
            indices = frame.index[invalid_shape].tolist()
            raise ValueError(f"{label} Geometry가 유효하지 않은 행이 있습니다: {indices}")

    # 이후 함수가 입력 파일의 원래 CRS를 신경 쓰지 않도록 둘 다 WGS84로 통일한다.
    # reset_index(drop=True)는 행 번호를 0, 1, 2, ...로 다시 붙인다. 이 연속된
    # 행 번호는 뒤에서 NumPy 높이 배열의 인덱스로도 사용된다.
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

    처리 단계는 다음과 같다.

    1. 원본 행 번호를 보존하고 id/숫자 필드를 정규화한다.
    2. 속성의 기준 좌표가 실제 헬리포트 Geometry 안에 있는지 검사한다.
    3. 같은 위치의 다른 레코드를 이용해 비어 있는 AGL을 보완한다.
    4. 미터 좌표계에서 가장 가까운 host 건물을 연결한다.
    5. 이후 함수가 쓰기 편한 ``HelipadRecord`` 목록으로 변환한다.
    """

    # copy()를 쓰므로 호출자가 넘긴 원본 GeoDataFrame은 수정되지 않는다.
    # reset_index(..., names=...)로 기존 행 번호를 __source_index 열에 보존한다.
    # '__' 접두사는 원본 열과 구분하려고 이 함수 안에서 임시로 쓰는 이름이다.
    prepared = helipads.copy().reset_index(drop=False, names="__source_index")
    prepared["__id"] = prepared[id_field].map(_identifier_text)

    # casefold()는 대소문자 비교용 변환이다. 즉 'ABC'와 'abc'도 같은 id로 본다.
    # duplicated(keep=False)는 중복 그룹의 첫 행까지 모두 True로 표시한다.
    if prepared["__id"].str.casefold().duplicated(keep=False).any():
        duplicates = prepared.loc[
            prepared["__id"].str.casefold().duplicated(keep=False), "__id"
        ].tolist()
        raise ValueError(f"중복된 헬리포트 id가 있습니다: {duplicates}")

    # 숫자로 바꿀 수 없는 문자열은 errors='coerce' 때문에 NaN이 된다. 다음의
    # isfinite 검사에서 NaN도 잡히므로 모든 잘못된 값을 같은 방식으로 보고한다.
    prepared["__lon"] = pd.to_numeric(prepared[longitude_field], errors="coerce")
    prepared["__lat"] = pd.to_numeric(prepared[latitude_field], errors="coerce")
    prepared["__agl"] = pd.to_numeric(prepared[agl_field], errors="coerce")
    # | 는 Series끼리의 논리 OR이다. 아래 조건 중 하나라도 참이면 잘못된 좌표다.
    invalid_coordinates = (
        ~np.isfinite(prepared["__lon"])
        | ~np.isfinite(prepared["__lat"])
        | ~prepared["__lon"].between(-180.0, 180.0)
        | ~prepared["__lat"].between(-90.0, 90.0)
    )
    if invalid_coordinates.any():
        bad_ids = prepared.loc[invalid_coordinates, "__id"].tolist()
        raise ValueError(f"경도/위도가 유효하지 않은 헬리포트 id: {bad_ids}")

    # 속성 열의 경도/위도로 실제 Point Geometry를 만든다. 이 점이 입력 Polygon의
    # 내부 또는 경계에 있어야 OLS 기준점으로 신뢰할 수 있다.
    coordinate_points = gpd.GeoSeries(
        gpd.points_from_xy(prepared["__lon"], prepared["__lat"]),
        index=prepared.index,
        crs="EPSG:4326",
    )
    # covers는 내부뿐 아니라 경계 위의 점도 포함한다는 점에서 contains와 다르다.
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

    # 동일한 경도/위도와 동일 Geometry를 가진 행을 peer(동일 대상의 다른 기록)로
    # 묶는다. 부동소수점의 미세한 표현 차이는 소수 12자리 반올림으로 줄이고,
    # Geometry는 도형을 바이트로 직렬화한 WKB 문자열로 정확히 비교한다.
    prepared["__peer_key"] = [
        (round(lon, 12), round(lat, 12), geometry.wkb_hex)
        for lon, lat, geometry in zip(
            prepared["__lon"], prepared["__lat"], prepared.geometry
        )
    ]
    # peer_values: 그룹별로 사용할 수 있는 양수 AGL 후보
    # peer_indices: 그룹별 원본 행 번호(나중에 자기 자신 충돌 제외에 사용)
    peer_values: dict[object, list[float]] = {}
    peer_indices: dict[object, tuple[int, ...]] = {}
    for key, group in prepared.groupby("__peer_key", sort=False):
        # set comprehension으로 같은 높이가 여러 번 있어도 한 번만 남긴다.
        values = sorted(
            {
                float(value)
                for value in group["__agl"]
                if np.isfinite(value) and float(value) > 0.0
            }
        )
        peer_values[key] = values
        peer_indices[key] = tuple(int(value) for value in group["__source_index"])

    # 원본 AGL이 유효하면 그대로 쓰고, 아니면 peer 그룹의 유일한 양수값을 쓴다.
    # 후보가 없거나 서로 다른 후보가 둘 이상이면 임의로 고르지 않고 실패한다.
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

    # 최근접 조인에는 필요한 열만 복사한다. 원본 건물 행 번호를 별도 열에 넣어
    # 좌표 변환/공간 조인 뒤에도 어떤 건물인지 추적할 수 있게 한다.
    building_match = buildings[["A16", "A17", "geometry"]].copy()
    building_match["__host_index"] = buildings.index.to_numpy()
    # 위/경도는 '도' 단위이므로 그대로 거리를 재면 m가 아니다. 두 레이어를 같은
    # 미터 단위 CRS로 바꾼 후 sjoin_nearest로 가장 가까운 건물을 찾는다.
    left = prepared[["__source_index", "geometry"]].to_crs(host_match_crs)
    right = building_match.to_crs(host_match_crs)
    nearest = gpd.sjoin_nearest(
        left,
        right,
        how="left",
        distance_col="__host_distance",
    )
    # 거리가 정확히 같은 건물이 여럿이면 host_index가 작은 것을 택해 실행할 때마다
    # 같은 결과가 나오게 한다. drop_duplicates는 헬리포트당 한 행만 남긴다.
    nearest = (
        nearest.sort_values(["__source_index", "__host_distance", "__host_index"])
        .drop_duplicates("__source_index", keep="first")
        .set_index("__source_index")
    )

    # 검증/보완이 끝난 DataFrame 행을 명시적인 데이터 클래스로 옮긴다.
    records: list[HelipadRecord] = []
    used_folders: dict[str, str] = {}
    for row in prepared.to_dict("records"):
        folder = _safe_folder_name(row["__id"])
        collision_key = folder.casefold()

        # 서로 다른 id가 금지문자 치환 후 같은 폴더명으로 합쳐지는 경우를 막는다.
        # 예: 'A:B'와 'A?B'는 둘 다 'A_B'가 된다.
        if collision_key in used_folders:
            raise ValueError(
                f"id 폴더명이 충돌합니다: "
                f"{used_folders[collision_key]!r}, {row['__id']!r}"
            )
        used_folders[collision_key] = row["__id"]
        # set_index를 해 두었으므로 원본 행 번호로 최근접 건물 결과를 바로 찾는다.
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

    반환되는 ``DemSample``에는 고도뿐 아니라 파일/밴드/행/열도 들어 있다.
    덕분에 나중에 '어떤 셀에서 이 값이 나왔는가'를 추적할 수 있다.
    """

    # 숫자와 범위를 먼저 확인해 좌표 변환 라이브러리에 잘못된 값을 넘기지 않는다.
    lat = _finite_number("latitude", latitude)
    lon = _finite_number("longitude", longitude)
    if not -90.0 <= lat <= 90.0:
        raise ValueError(f"위도는 -90~90도여야 합니다: {lat}")
    if not -180.0 <= lon <= 180.0:
        raise ValueError(f"경도는 -180~180도여야 합니다: {lon}")
    # Python에서 bool은 int의 하위 타입이라 isinstance(True, int)가 True이다.
    # 그래서 True/False가 밴드 번호로 통과하지 않도록 별도로 제외한다.
    if not isinstance(band, int) or isinstance(band, bool) or band < 1:
        raise ValueError(f"DEM 밴드는 1 이상의 정수여야 합니다: {band!r}")

    path = Path(dem_path).expanduser().resolve()
    if not path.is_file():
        raise ValueError(f"DEM 파일을 찾을 수 없습니다: {path}")

    try:
        # with 블록을 벗어나면 rasterio가 파일 핸들을 자동으로 닫는다.
        with rasterio.open(path) as dataset:
            if dataset.crs is None:
                raise ValueError(f"DEM에 CRS가 없습니다: {path}")
            if band > dataset.count:
                raise ValueError(
                    f"DEM 밴드 {band}가 없습니다. 사용 가능한 밴드: 1~{dataset.count}"
                )

            # always_xy=True는 좌표 순서를 언제나 (X=경도, Y=위도)로 고정한다.
            # 일부 CRS가 위도/경도 순서를 정의하더라도 코드의 순서가 뒤집히지 않는다.
            transformer = Transformer.from_crs(
                "EPSG:4326", dataset.crs, always_xy=True
            )
            x, y = transformer.transform(lon, lat)
            # 실세계 좌표 (x, y)가 래스터 배열의 몇 행/몇 열인지 계산한다.
            row, column = dataset.index(x, y)
            if not (0 <= row < dataset.height and 0 <= column < dataset.width):
                raise ValueError(
                    "기준점이 DEM 범위 밖입니다: "
                    f"lat={lat}, lon={lon}, DEM CRS 좌표=({x:.3f}, {y:.3f}), "
                    f"DEM bounds={dataset.bounds}"
                )

            # DEM 전체를 읽지 않고 필요한 1×1 셀만 읽어 메모리와 시간을 절약한다.
            cell = dataset.read(
                band,
                window=Window(column, row, 1, 1),
                masked=True,
            )
            # masked=True이면 NoData 셀이 마스크된다. 마스크가 True인 셀은 사용할 수 없다.
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
    # from error는 원래 예외를 연결해 두어 디버깅할 때 상세 원인도 볼 수 있게 한다.
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
    """삼각형 생성 인수를 float로 통일하고 허용 범위를 검사한다.

    검사를 ``create_ols_triangle``과 분리해 두면 본 계산 부분은 좌표 기하 로직에
    집중할 수 있다. 반환 순서는 입력과 같으며 각 값은 유한한 ``float``이다.
    """

    # 튜플을 먼저 만들고 구조 분해 할당하여 변환/검사를 한곳에 모은다.
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

    Parameters
    ----------
    latitude, longitude
        움직이지 않는 삼각형 꼭짓점(apex)의 WGS84 좌표.
    apex_agl_m, ground_z_m
        꼭짓점의 지상고와 해당 위치의 DEM 지표고. 둘을 더하면 절대 표면고이다.
    rotation_deg
        진북 0도, 동쪽 90도, 남쪽 180도, 서쪽 270도인 시계방향 방위각.
    horizontal_height_m, base_width_m
        삼각형 중심선 길이와 밑변 폭(m).
    slope_deg
        꼭짓점에서 밑변으로 갈수록 높아지는 OLS 경사각(도).

    Returns
    -------
    OLSTriangle
        세 외곽 꼭짓점, 밑변 중심, 높이 계산 정보를 모두 담은 객체.
    """

    # 이후 계산은 검증을 통과한 float만 사용한다.
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

    # 나머지 연산으로 -10도는 350도, 370도는 10도처럼 0~360도 안으로 맞춘다.
    # 단, 결과 속성 rotation_deg에는 사용자가 넣은 원래 각도를 보존한다.
    bearing = rotation % 360.0

    # Geod.fwd(시작 경도, 시작 위도, 방위각, 거리)는 타원체 표면을 따라 이동한
    # (도착 경도, 도착 위도, 도착점의 역방위각)을 반환한다.
    base_center_lon, base_center_lat, back_azimuth = WGS84_GEOD.fwd(
        lon, lat, bearing, height
    )
    # fwd가 반환하는 것은 도착점에서 출발점을 향하는 역방위각이다.
    # 180도를 더해 도착점의 진행방향을 얻고, 그 좌우 90도에 밑변을 놓는다.
    terminal_azimuth = (back_azimuth + 180.0) % 360.0
    # 밑변 중심에서 진행방향 기준 오른쪽/왼쪽으로 절반 폭만큼 이동한다.
    # 두 점 사이의 전체 거리는 width가 된다.
    half_width = width / 2.0
    right_lon, right_lat, _ = WGS84_GEOD.fwd(
        base_center_lon, base_center_lat, terminal_azimuth + 90.0, half_width
    )
    left_lon, left_lat, _ = WGS84_GEOD.fwd(
        base_center_lon, base_center_lat, terminal_azimuth - 90.0, half_width
    )

    # math.tan은 라디안 단위를 받으므로 degrees → radians 변환이 필요하다.
    # 직각삼각형에서 tan(경사각) = 수직상승량 / 수평거리이므로
    # 수직상승량 = 수평거리 × tan(경사각)이다.
    vertical_rise = height * math.tan(math.radians(slope))

    # 절대고 = 지표고(DEM) + 지면 위 높이(AGL). 밑변은 여기에 경사 상승량을 더한다.
    apex_z = ground + agl
    base_z = apex_z + vertical_rise

    # 위치 인수는 Vertex 선언 순서(name, latitude, longitude, z_m)에 맞춘다.
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

    ``range``는 정수 간격만 처리하므로, 0.5도 같은 실수 간격도 지원하기 위해
    직접 반복문으로 값을 만든다. 반환값은 수정 불가능한 튜플이다.
    """

    start = _finite_number("start_deg", start_deg)
    end = _finite_number("end_deg", end_deg)
    step = _finite_number("step_deg", step_deg)
    if end < start:
        raise ValueError(f"종료각은 시작각 이상이어야 합니다: {start} > {end}")
    if step <= 0.0:
        raise ValueError(f"회전 간격은 0도보다 커야 합니다: {step}")

    # 부동소수점에서는 0.1을 반복해서 더한 결과가 0.3과 정확히 같지 않을 수 있다.
    # tolerance는 그 아주 작은 오차 때문에 종료각이 빠지거나 추가되는 것을 막는다.
    tolerance = max(1.0, abs(start), abs(end)) * 1e-12
    angles: list[float] = []
    index = 0
    while True:
        # 이전 angle에 step을 누적하지 않고 매번 start + index*step으로 계산하면
        # 반복할수록 부동소수점 오차가 쌓이는 현상을 조금 줄일 수 있다.
        angle = start + index * step
        inside = angle <= end + tolerance if include_end else angle < end - tolerance
        if not inside:
            break
        # 계산값이 종료각과 사실상 같다면 사용자가 입력한 end 자체를 저장한다.
        angles.append(
            end
            if include_end and math.isclose(angle, end, abs_tol=tolerance)
            else angle
        )
        index += 1
        # 지나치게 작은 step이 메모리를 고갈시키는 것을 막는 안전장치이다.
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
    """같은 기준점에서 방위각만 바꾼 OLS 삼각형 목록을 만든다.

    아래 리스트 컴프리헨션은 일반적인 ``for`` 반복문을 짧게 쓴 것이다. 각
    ``angle``에 대해 ``create_ols_triangle(...)``을 한 번 호출하고 그 반환값을
    새 리스트에 차례로 넣는다. 이 단계에서는 파일을 저장하지 않는다.
    """

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

    기준 헬리포트와 달리 장애물은 여러 개이므로, 값을 읽지 못한 한 객체 때문에
    전체 작업을 중단하지 않는다. 그 객체의 높이만 ``nan``(미상)으로 표시한다.
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
            # centroid(무게중심)는 ㄷ자 같은 오목한 도형의 바깥에 놓일 수 있다.
            # representative_point는 반드시 각 도형 내부에 있는 점을 돌려준다.
            points = frame.geometry.representative_point()
            transformer = Transformer.from_crs(
                frame.crs, dataset.crs, always_xy=True
            )
            # 모든 대표점을 DEM과 같은 CRS로 변환한다.
            coordinates = [
                transformer.transform(point.x, point.y) for point in points
            ]
            # 먼저 모든 값을 '높이 미상'인 nan으로 채우고, 정상 표본만 덮어쓴다.
            elevations = np.full(len(frame), np.nan, dtype=float)

            # dataset.sample은 좌표별 래스터 값을 순서대로 생성하는 iterator이다.
            # enumerate가 0부터 인덱스를 붙여 GeoDataFrame 행과 같은 위치에 기록한다.
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

    이 함수에서 비싼 전처리를 한 번만 수행해 ``CollisionLayers``로 묶는다. 이후
    360개 회전 후보가 같은 좌표 변환 도형과 높이 배열을 반복해서 재사용한다.
    """

    if "A16" not in buildings.columns:
        raise ValueError("건물 Shapefile에 높이 필드 A16이 없습니다.")

    # 건물별 지표고와 AGL을 각각 1차원 NumPy 배열로 만든다.
    building_ground = _sample_dem_at_geometry_points(
        buildings, dem_path=dem_path, band=dem_band
    )
    building_agl = pd.to_numeric(buildings["A16"], errors="coerce").to_numpy(
        dtype=float
    )
    # 벡터화 연산: Python for문 없이 모든 건물에 세 조건을 한 번에 적용한다.
    # DEM과 AGL이 둘 다 유한하고 AGL이 양수인 건물만 높이를 아는 것으로 본다.
    valid_building_height = (
        np.isfinite(building_ground)
        & np.isfinite(building_agl)
        & (building_agl > 0.0)
    )
    building_top = np.full(len(buildings), np.nan, dtype=float)

    # Boolean 배열이 True인 위치만 골라 '지표고 + 건물높이'를 저장한다.
    building_top[valid_building_height] = (
        building_ground[valid_building_height]
        + building_agl[valid_building_height]
    )

    # 헬리포트는 앞 단계에서 AGL과 기준점 DEM을 엄격히 검증했으므로 각 원본
    # 행 번호 위치에 확정된 상단고를 넣을 수 있다.
    helipad_top = np.full(len(helipads), np.nan, dtype=float)
    for record in helipad_records:
        sample = helipad_dem_samples[record.identifier]
        helipad_top[record.source_index] = sample.elevation_m + record.agl_m

    # 충돌의 수평 부분은 2차원 footprint로 판단한다. 입력 Geometry에 Z가 있더라도
    # force_2d로 제거한다. 높이는 위의 별도 top_z 배열과 비교한다.
    buildings_metric = buildings.copy()
    buildings_metric.geometry = buildings_metric.geometry.map(shapely.force_2d)
    # EPSG:4326의 도 단위가 아니라 미터 좌표계에서 벡터 거리/투영을 계산한다.
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
    """장애물 상단이 겹치는 지점의 OLS 높이에 닿거나 뚫는지 판정한다.

    ``footprint``와 ``obstacle_geometry``는 모두 미터 단위 2차원 도형이다.
    ``apex_xy``는 삼각형 꼭짓점 벡터, ``direction``은 꼭짓점에서 밑변으로 향하는
    길이 1의 단위벡터이다. 반환값이 ``True``이면 해당 방위각을 제거한다.

    OLS는 중심선 진행거리 ``d``에 따라 선형으로 높아진다.

    ``surface_z(d) = apex_z + d * (vertical_rise / horizontal_height)``

    장애물과 겹치는 좌표 중 꼭짓점에 가장 가까운 곳의 OLS가 가장 낮다. 장애물
    상단이 그 최저 높이 이상이면 겹친 영역 어딘가에서 표면에 닿는다고 판단한다.
    """

    # 높이를 모르는 장애물은 수평 도형만으로 충돌이라고 추측하지 않는다.
    if not math.isfinite(obstacle_top_z_m):
        return False

    # intersection은 두 footprint가 실제로 공통으로 차지하는 도형을 계산한다.
    overlap = footprint.intersection(obstacle_geometry)
    if overlap.is_empty:
        return False
    # 겹친 선/면/점의 좌표를 N×2 NumPy 배열로 꺼낸다.
    coordinates = shapely.get_coordinates(overlap)
    if len(coordinates) == 0:
        point = overlap.representative_point()
        coordinates = np.array([[point.x, point.y]], dtype=float)
    # (중첩점 - 꼭짓점) @ direction은 벡터 내적(dot product)이다. 즉 각 중첩점을
    # 중심선 방향에 투영한 거리 d를 구한다. clip은 수치오차나 경계 조건 때문에
    # d가 삼각형 범위 0~height 밖으로 나가지 않도록 제한한다.
    along_distances = np.clip(
        (coordinates[:, :2] - apex_xy) @ direction,
        0.0,
        triangle.horizontal_height_m,
    )
    # 수평으로 1 m 갈 때 OLS 표면이 몇 m 올라가는지를 계산한다.
    rise_per_metre = triangle.vertical_rise_m / triangle.horizontal_height_m
    minimum_surface_z = (
        triangle.apex_z_m + float(np.min(along_distances)) * rise_per_metre
    )
    # tolerance_m만큼 여유를 두어 0.0000001 m 같은 계산 오차로 판정이 뒤집히지
    # 않게 한다. '닿는 것'도 안전하지 않다고 보므로 >= 연산자를 쓴다.
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

    Returns
    -------
    tuple[list[OLSTriangle], tuple[float, ...]]
        첫 번째 값은 안전한 삼각형 목록, 두 번째 값은 충돌한 회전각 튜플이다.
    """

    # list보다 set은 '이 인덱스가 제외 대상인가?'를 빠르게 검사할 수 있다.
    excluded_helipads = set(helipad.excluded_helipad_indices)

    # 후보 삼각형은 저장용 EPSG:4326 좌표이므로 충돌 계산용 미터 CRS로 바꾼다.
    transformer = Transformer.from_crs(
        "EPSG:4326", collision_layers.analysis_crs, always_xy=True
    )
    safe: list[OLSTriangle] = []
    conflicts: list[float] = []
    for triangle in triangles:
        # 각 외곽 꼭짓점의 (경도, 위도)를 미터 단위 (X, Y)로 변환한다. Z는 아래의
        # 별도 평면식에서 계산하므로 여기서는 전달하지 않는다.
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
        # 방향벡터 = 도착점 - 시작점. 그 길이(norm)로 나누면 길이가 1인 단위벡터다.
        direction = base_center_xy - apex_xy
        direction /= np.linalg.norm(direction)

        # 공간 인덱스(sindex)는 모든 건물을 하나씩 검사하지 않고, 이 삼각형과
        # 수평으로 교차할 가능성이 있는 행 번호만 빠르게 찾아 준다.
        building_hits = set(
            collision_layers.buildings.sindex.query(
                footprint, predicate="intersects"
            ).tolist()
        )
        # discard는 값이 없어도 오류가 나지 않는다. 현재 헬리포트가 놓인 host
        # 건물은 의도된 중첩이므로 장애물 후보에서 제외한다.
        building_hits.discard(helipad.host_index)
        helipad_hits = set(
            collision_layers.helipads.sindex.query(
                footprint, predicate="intersects"
            ).tolist()
        )
        # 자기 자신 및 같은 Geometry의 peer 헬리포트를 한꺼번에 제외한다.
        helipad_hits.difference_update(excluded_helipads)

        # any(...)는 후보 중 하나라도 True이면 즉시 평가를 멈추고 True가 된다.
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
        # 건물 또는 다른 헬리포트 중 하나라도 표면을 침범하면 충돌각으로 분류한다.
        if building_collision or helipad_collision:
            conflicts.append(triangle.rotation_deg)
        else:
            safe.append(triangle)
    return safe, tuple(conflicts)


def format_angle_ranges(
    angles_deg: Sequence[float], step_deg: float = 1.0
) -> str:
    """연속된 충돌각을 사람이 읽기 쉬운 범위 문자열로 압축한다.

    예를 들어 ``[0, 1, 2, 8, 9]``와 간격 1을 받으면
    ``"[0,2] & [8,9]"``를 반환한다. 실제 계산용 자료가 아니라 요약 CSV의
    표시용 문자열을 만드는 함수이다.
    """

    if len(angles_deg) == 0:
        return ""
    step = _finite_number("step_deg", step_deg)
    if step <= 0.0:
        raise ValueError(f"각도 간격은 0보다 커야 합니다: {step}")
    # set으로 중복을 없애고 sorted로 작은 각도부터 정렬한다.
    angles = sorted(set(float(angle) for angle in angles_deg))
    ranges: list[tuple[float, float]] = []
    # start는 현재 연속 구간의 시작, previous는 바로 전에 본 각도이다.
    start = previous = angles[0]
    tolerance = max(1.0, abs(step)) * 1e-9
    for angle in angles[1:]:
        # 현재값 - 이전값이 step과 같으면 같은 연속 구간을 계속 확장한다.
        if math.isclose(angle - previous, step, abs_tol=tolerance):
            previous = angle
            continue
        # 간격이 끊겼다면 지금까지의 구간을 저장하고 새 구간을 시작한다.
        ranges.append((start, previous))
        start = previous = angle
    # 반복문 안에서는 '다음 구간이 시작될 때'만 저장하므로 마지막 구간을 추가한다.
    ranges.append((start, previous))
    return " & ".join(f"[{start:g},{end:g}]" for start, end in ranges)


def _triangle_record(
    triangle: OLSTriangle,
    dem_sample: DemSample,
    helipad: HelipadRecord,
) -> dict[str, float | int | str | Polygon]:
    """삼각형 객체 하나를 Shapefile의 속성 한 행(dict)으로 펼친다.

    ESRI Shapefile의 속성은 DBF 파일에 저장되며 필드명이 10자로 제한된다.
    따라서 ``horizontal_height_m`` 같은 긴 Python 이름 대신 ``height_m``처럼
    짧은 이름을 사용한다. 마지막 ``geometry`` 값은 공간 도형 자체이다.
    """

    return {
        "heli_id": helipad.identifier,  # 헬리포트 식별자
        "rot_deg": triangle.rotation_deg,  # 진북 기준 시계방향 회전각
        "src_fid": helipad.source_index,  # 헬리포트 원본 행 번호
        "agl_src": helipad.agl_source,  # AGL 원본/peer 보완 여부
        "apex_lat": triangle.apex.latitude,  # 기준 꼭짓점 위도
        "apex_lon": triangle.apex.longitude,  # 기준 꼭짓점 경도
        "apex_agl": triangle.apex_agl_m,  # 기준 꼭짓점의 지상고
        "ground_z": triangle.ground_z_m,  # 기준 꼭짓점 DEM 지표고
        "apex_z": triangle.apex_z_m,  # 지표고 + AGL
        "right_lat": triangle.base_right.latitude,  # 밑변 오른쪽 위도
        "right_lon": triangle.base_right.longitude,  # 밑변 오른쪽 경도
        "left_lat": triangle.base_left.latitude,  # 밑변 왼쪽 위도
        "left_lon": triangle.base_left.longitude,  # 밑변 왼쪽 경도
        "base_z": triangle.base_z_m,  # 밑변의 절대 표면고
        "rise_m": triangle.vertical_rise_m,  # 꼭짓점→밑변 수직 상승량
        "height_m": triangle.horizontal_height_m,  # 수평 중심선 길이
        "width_m": triangle.base_width_m,  # 밑변 폭
        "slope_deg": triangle.slope_deg,  # OLS 종단 경사각
        "h_crs": "EPSG:4326",  # X/Y 좌표 기준
        "z_ref": "DEM+AGL",  # Z가 타원체고가 아님을 명시
        "dem_file": dem_sample.path.name,  # DEM 파일명
        "dem_crs": dem_sample.crs,  # DEM 좌표계
        "dem_band": dem_sample.band,  # 사용한 DEM 밴드
        "dem_row": dem_sample.row,  # 기준점 DEM 행
        "dem_col": dem_sample.column,  # 기준점 DEM 열
        "host_fid": helipad.host_index,  # 제외한 host 건물 행 번호
        "host_dist": helipad.host_distance_m,  # host까지 거리(m)
        "host_a16": helipad.host_a16_m,  # host의 원본 A16
        "host_a17": helipad.host_a17_m,  # host의 원본 A17
        "geometry": triangle.polygon,  # 최종 Polygon Z
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

    현재 전체 실행 흐름은 안전각이 없을 때 이 함수를 호출하지 않지만, 빈 결과를
    명시적으로 저장하려는 다른 호출자도 사용할 수 있도록 빈 파일 기능을 둔다.
    반환값은 정규화된 출력 ``Path``이다.
    """

    # 확장자를 검사해 사용자가 실수로 CSV 등 다른 파일에 쓰는 것을 막는다.
    output = Path(output_path).expanduser().resolve()
    if output.suffix.lower() != ".shp":
        raise ValueError(f"Shapefile 출력 경로는 .shp여야 합니다: {output}")
    # parents=True는 상위 폴더까지 만들고, exist_ok=True는 이미 있어도 허용한다.
    output.parent.mkdir(parents=True, exist_ok=True)
    if len(triangles) > 0:
        # 삼각형 객체를 DBF 속성 + Geometry 형태의 dict로 하나씩 바꾼다.
        records = [
            _triangle_record(triangle, dem_sample, helipad)
            for triangle in triangles
        ]
        # 일반 dict 목록에 geometry 열과 CRS를 지정해 공간 DataFrame으로 만든다.
        frame = gpd.GeoDataFrame(
            records,
            geometry="geometry",
            crs="EPSG:4326",
        )
        # Shapefile은 실제로 .shp/.shx/.dbf/.prj 등 여러 파일 묶음으로 저장된다.
        # index=False는 GeoDataFrame의 행 번호를 불필요한 속성 열로 쓰지 않는다.
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
        # 행이 0개이면 GeoPandas가 열 자료형을 추론할 수 없으므로 먼저 정상 행
        # 하나로 스키마를 만든 뒤 iloc[0:0]으로 행만 제거한다.
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
    """안전각이 없을 때 이전 실행의 동일 Shapefile 구성 파일만 제거한다.

    Shapefile은 여러 sidecar 파일의 묶음이므로 ``.shp``만 지우면 불완전한 찌꺼기가
    남는다. 삭제 전에 경로가 허용된 출력 루트 아래인지 확인하여 입력 자료나 다른
    파일을 잘못 지우는 일을 막는다.
    """

    output = Path(output_path).expanduser().resolve()
    root = Path(allowed_root).expanduser().resolve()
    # resolve된 절대 경로끼리 비교하므로 '..'를 이용한 루트 탈출도 차단된다.
    if not output.is_relative_to(root):
        raise ValueError(f"출력 루트 밖의 파일은 제거할 수 없습니다: {output}")
    # 드라이버/프로그램에 따라 생길 수 있는 대표적인 Shapefile 부속 확장자들이다.
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
    # 폴더가 완전히 비었을 때만 지운다. 다른 파일이 있으면 rmdir가 실패하지만
    # 의도된 안전 동작이므로 OSError를 무시하고 폴더를 보존한다.
    try:
        output.parent.rmdir()
    except OSError:
        # 생성물 외의 파일이 있으면 폴더는 보존한다.
        pass


def assert_shapefile_unlocked(output_path: str | Path, allowed_root: str | Path) -> None:
    """기존 출력 세트가 QGIS 등 다른 프로세스에 잠겼는지 미리 확인한다.

    Windows에서는 QGIS가 열어 둔 파일을 덮어쓰려 하면 작업 중간에 실패할 수 있다.
    자기 자신과 같은 이름으로 ``os.rename``을 시도하면 내용은 바꾸지 않으면서
    쓰기 가능한 상태인지 확인할 수 있다.
    """

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
            # 출발/도착 경로가 같으므로 실제 파일 위치나 내용은 변하지 않는다.
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
    """전체 계산을 지휘하고 헬리포트별 결과 파일과 요약 CSV를 만든다.

    이 함수는 세부 수학을 직접 처리하기보다 앞에서 정의한 작은 함수들을 순서대로
    호출하는 '오케스트레이터(orchestrator)' 역할을 한다. 중요한 설계 원칙은
    가능한 검증을 모두 끝낸 다음 파일을 변경하는 것이다.

    처리 순서
    ---------
    1. 모든 헬리포트 기준점의 DEM과 공통 충돌 레이어를 준비한다.
    2. 헬리포트마다 모든 방위각 후보를 만들고 안전/충돌로 나눈다.
    3. 모든 기존 출력 파일이 잠기지 않았는지 확인한다.
    4. 안전 후보만 Shapefile로 쓰고 전체 결과를 CSV로 저장한다.

    Returns
    -------
    list[GenerationResult]
        입력 헬리포트와 같은 순서의 처리 결과 목록.
    """

    if len(helipad_records) == 0:
        raise ValueError("생성할 헬리포트가 없습니다.")
    if len(angles_deg) == 0:
        raise ValueError("생성할 회전각이 없습니다.")

    # 출력 폴더를 만들기 전에 모든 기준점의 DEM을 먼저 검증한다.
    # dict comprehension으로 {헬리포트 id: 해당 기준점 DEM 정보} 사전을 만든다.
    samples = {
        helipad.identifier: sample_dem_at_reference_point(
            latitude=helipad.latitude,
            longitude=helipad.longitude,
            dem_path=dem_path,
            band=dem_band,
        )
        for helipad in helipad_records
    }
    # 건물/헬리포트의 투영 도형과 상단고는 모든 회전각에서 공통으로 사용된다.
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
    # :g 형식은 8.0을 '8'처럼 불필요한 소수점 없이 폴더명으로 만든다.
    slope_directory = root / f"{slope_value:g}"
    slope_directory.mkdir(parents=True, exist_ok=True)
    # 아직 파일을 쓰지 않고 검증 결과를 메모리에 모은다. 아래 긴 타입은 튜플에
    # (헬리포트, DEM, 전체 후보, 안전 후보, 충돌각, 출력경로)가 든다는 뜻이다.
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
        # 예: .../Polygon/8/헬리포트ID/Polygon_헬리포트ID.shp
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
    # 밑줄 변수는 '이 값은 여기서 사용하지 않는다'는 Python 관례이다.
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
        # 빈 list는 False로 평가된다. 안전 삼각형이 하나 이상일 때만 저장한다.
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
        # 화면 출력과 CSV 작성에 필요한 최소 결과를 데이터 클래스로 보관한다.
        results.append(
            GenerationResult(
                helipad=helipad,
                dem_sample=sample,
                output_path=final_output,
                polygon_count=len(safe_triangles),
                collision_angles=collision_angles,
            )
        )

    # 리스트 컴프리헨션으로 헬리포트마다 dict 한 개를 만들고 표로 변환한다.
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
                # 조건 표현식: 출력이 없으면 None 대신 빈 문자열을 CSV에 쓴다.
                "output_path": str(result.output_path) if result.output_path else "",
            }
            for result in results
        ]
    )
    summary_path = slope_directory / "generation_summary.csv"
    # utf-8-sig는 Excel에서 한글 CSV를 열 때 글자가 깨질 가능성을 줄인다.
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    legacy_summary = root / "generation_summary.csv"
    # 과거 버전이 상위 폴더에 만든 요약 파일은 새 결과와 혼동되지 않게 제거한다.
    if legacy_summary.is_file() and legacy_summary != summary_path:
        legacy_summary.unlink()
    return results


def _build_parser() -> argparse.ArgumentParser:
    """터미널에서 받을 명령행 옵션과 도움말을 정의한다.

    ``argparse``는 문자열로 들어온 옵션을 ``Path``/``float``/``int`` 등 지정한
    자료형으로 변환한다. 옵션을 생략하면 ``SCRIPT_CONFIG``의 값이 사용된다.
    이 함수는 옵션을 '정의'만 하고, 실제 해석은 ``main``의 ``parse_args``가 한다.
    parser는 원하는 변수를 터미널에서 입력하는 방식으로 input과는 조금 다르다.
    만일 방위각을 하나만 설정하길 희망한다면, python  코드\polygon\장애물접근표면.py --start-angle 0의 형태로 진행하면 됨.
    """

    parser = argparse.ArgumentParser(
        description=(
            "송파구 옥상헬리포트의 id/x_long/y_lati/A16을 읽어 OLS 경사/id별 "
            "폴더에 건물·다른 헬리포트가 OLS 표면고에 닿거나 관통하지 않는 "
            "삼각형 Polygon Z만 저장합니다."
        )
    )
    # 입력/출력 경로 옵션 -------------------------------------------------
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
    # 입력 Shapefile의 필드명 옵션 ---------------------------------------
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
    # 좌표계와 회전 범위 옵션 --------------------------------------------
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
    # action=BooleanOptionalAction은 --include-end-angle과
    # --no-include-end-angle 두 형태를 자동으로 지원한다.
    parser.add_argument(
        "--include-end-angle",
        action=argparse.BooleanOptionalAction,
        default=SCRIPT_CONFIG.include_end_angle,
        help="종료각 포함 여부(360도는 0도와 같은 중복 형상임)",
    )
    # OLS 형상 크기/경사 옵션 ---------------------------------------------
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
    """전체 실행 결과를 터미널에 출력할 여러 줄 문자열로 만든다."""

    # peer AGL을 쓴 id만 골라 사용자가 보완 발생 여부를 바로 알 수 있게 한다.
    fallback_ids = [
        result.helipad.identifier
        for result in results
        if result.helipad.agl_source.startswith("peer_")
    ]
    host_distances = [result.helipad.host_distance_m for result in results]
    slope_directory = output_root / f"{float(slope_deg):g}"
    # generator expression을 sum에 넘겨 모든 헬리포트의 개수를 합산한다.
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
    # 문자열 리스트 사이에 줄바꿈을 끼워 하나의 출력 문자열로 합친다.
    return "\n".join(rows)


def main(argv: Sequence[str] | None = None) -> int:
    """명령행 실행의 시작점이며 전체 함수를 올바른 순서로 호출한다.

    ``argv=None``이면 실제 터미널 인수를 사용한다. 테스트에서는 문자열 목록을
    직접 넘겨 원하는 옵션을 재현할 수 있다. 성공하면 운영체제 종료 코드 0을,
    입력 오류가 있으면 ``argparse``를 통해 오류 메시지와 종료 코드 2를 낸다.
    """

    parser = _build_parser()

    # 여기서 사용자가 입력한 문자열 옵션이 args.buildings 같은 속성으로 바뀐다.
    args = parser.parse_args(argv)

    try:
        # 1) 원본 공간 파일 읽기 및 공통 검증
        buildings, helipad_frame = load_source_layers(
            building_path=args.buildings,
            helipad_path=args.helipads,
            id_field=args.id_field,
            longitude_field=args.longitude_field,
            latitude_field=args.latitude_field,
            agl_field=args.agl_field,
        )
        # 2) 헬리포트별 AGL/host/기준 좌표 정리
        helipads = prepare_helipad_records(
            buildings=buildings,
            helipads=helipad_frame,
            id_field=args.id_field,
            longitude_field=args.longitude_field,
            latitude_field=args.latitude_field,
            agl_field=args.agl_field,
            host_match_crs=args.host_match_crs,
        )
        # 3) 검사할 회전각 목록 생성
        angles = rotation_angles(
            start_deg=args.start_angle,
            end_deg=args.end_angle,
            step_deg=args.angle_step,
            include_end=args.include_end_angle,
        )
        # 4) OLS 생성 → 충돌 판정 → 파일 저장 → CSV 작성
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
    # 예상 가능한 입력/파일 오류는 긴 traceback 대신 사용법과 함께 보여 준다.
    except (OSError, ValueError) as error:
        parser.error(str(error))

    # 5) 사용자가 확인할 수 있도록 최종 통계를 표준 출력에 표시한다.
    output_root = Path(args.output_root).expanduser().resolve()
    print(_format_batch_summary(results, angles, output_root, args.ols_angle))
    return 0


# 이 파일을 직접 실행할 때만 main()을 호출한다. 다른 Python 파일이 import하면
# 함수/클래스 정의만 불러오고 전체 생성 작업은 자동으로 시작하지 않는다.
if __name__ == "__main__":
    # SystemExit에 0을 전달해 운영체제/배치 스크립트에 성공 상태를 알린다.
    raise SystemExit(main())
