"""헬리포트별 양방향 선형 접근면과 건물의 3차원 충돌을 검사한다.

입력은 읽기 전용이며, 모든 산출물과 임시 파일은 자료/결과물에만 저장한다.
접근면 수치는 제공되지 않았으므로 CONFIG 또는 명령행에서 명시해야 한다.
수치 미설정 시에도 입력 검증 결과와 판정 불가(빈 값)인 각도표를 저장하고
종료 코드 2를 반환한다. 빈 값은 '충돌 없음'을 뜻하지 않는다.

실행 환경: requirements.txt를 설치한 Python 또는 QGIS Python.
pip install gdal은 필요 없다. Rasterio → 이미 설치된 osgeo → QGIS GDAL 실행 파일
순서로 DEM 읽기 환경을 자동 선택하므로 일반 python 명령으로 실행할 수 있다.
    python -B 코드/지형/collide_detection.py --help
    python -B 코드/지형/collide_detection.py --inspect-only

가정: 건물은 일정 높이의 수직 기둥, 헬리포트는 host 건물 지붕 높이,
접근면은 일정한 비음수 경사의 사다리꼴이다. 각도는 격자북 0°, 시계 방향.
DEM과 건물 높이는 m 단위이며 DEM의 수직 기준을 공통으로 사용한다.
법정 장애물제한표면의 종류/수치나 수직 기준의 적합성은 자동으로 결정하지 않는다.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict, dataclass
import hashlib
import json
import logging
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Iterator

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import CRS, Transformer
import pyogrio
import shapely
from shapely.affinity import affine_transform
from shapely.geometry import Point, Polygon
from shapely.ops import transform, unary_union


# ============================== CONFIG ==============================
PROJECT_ROOT = Path(r"C:\Users\choih\Desktop\연구\Landing_Depature_Process")
SCRIPT_DIR = PROJECT_ROOT / "코드" / "지형"
SCRIPT_PATH = SCRIPT_DIR / "collide_detection.py"
INPUT_DIR = PROJECT_ROOT / "자료" / "안전성" / "송파구"
BUILDING_PATH = INPUT_DIR / "송파구_건물_수정.shp"
HELIPAD_PATH = INPUT_DIR / "송파구_옥상헬리포트_수정.shp"
DEM_PATH = INPUT_DIR / "송파구_DEM.tif"
# --host-source reference를 명시한 경우에만 읽는다.
HOST_REFERENCE_PATH = INPUT_DIR / "송파구_건물.shp"
OUTPUT_DIR = PROJECT_ROOT / "자료" / "결과물"
TEMP_DIR = OUTPUT_DIR / "_temp"
COLLISION_DETAIL_PATH = OUTPUT_DIR / "collision_detail.csv"
BUILDING_SUMMARY_PATH = OUTPUT_DIR / "collision_building_summary.csv"
ANGLE_SUMMARY_PATH = OUTPUT_DIR / "helipad_angle_summary.csv"
COLLISION_GPKG_PATH = OUTPUT_DIR / "collision_buildings.gpkg"
UNMATCHED_HELIPADS_PATH = OUTPUT_DIR / "unmatched_helipads.csv"
LOG_PATH = OUTPUT_DIR / "analysis.log"
VALIDATION_PATH = OUTPUT_DIR / "validation_report.txt"
HEIGHT_DIAGNOSTICS_PATH = OUTPUT_DIR / "building_height_diagnostics.csv"
HOST_MATCHES_PATH = OUTPUT_DIR / "helipad_host_matches.csv"
ANGLE_DIAGNOSTICS_PATH = OUTPUT_DIR / "angle_diagnostics.csv"

DEM_BACKEND = "auto"  # auto / rasterio / gdal / gdal-cli
GDAL_BIN_DIR: Path | None = None  # 별도 설치 경로가 있으면 bin 폴더 지정
GDAL_SEARCH_ROOTS = (
    Path(os.environ.get("ProgramFiles", r"C:\Program Files")),
    Path(r"C:\OSGeo4W"),
    Path(r"C:\OSGeo4W64"),
)

ANALYSIS_CRS = "EPSG:5186"
SHP_ENCODING = "UTF-8"
BUILDING_HEIGHT_FIELD = "A16"  # 지상 높이(m), 데이터 제공자의 필드 정의 확인 필요
BUILDING_FLOORS_FIELD = "A26"
BUILDING_ID_FIELD: str | None = None  # None: 원본 SHP FID로 B000000 형식 생성
HELIPAD_ID_FIELD: str | None = None   # None: 원본 SHP FID로 H000000 형식 생성
GROUND_STATISTIC = "median"          # 건물과 면적으로 겹치는 DEM 셀의 통계
MIN_DEM_COVERAGE = 0.999999
MIN_HOST_OVERLAP_RATIO = 0.5
FLOOR_HEIGHT_M: float | None = None   # 명시한 경우에만 층수로 결측 높이 추정
HOST_SOURCE = "buildings"            # buildings / reference / helipad

# 아래 다섯 값은 의도적으로 미설정. 임의의 항공 기준을 대입하지 않는다.
APPROACH_LENGTH_M: float | None = None
APPROACH_INNER_WIDTH_M: float | None = None
APPROACH_OUTER_WIDTH_M: float | None = None
APPROACH_SLOPE_PERCENT: float | None = None
APPROACH_START_OFFSET_M: float | None = None
HELIPAD_ROOF_OFFSET_M = 0.0
PENETRATION_TOLERANCE_M = 1e-6
AXIS_ANGLES = tuple(range(180))
# ====================================================================

DETAIL_COLUMNS = [
    "helipad_id", "host_building_id", "building_id", "building_agl_m",
    "building_ground_z_m", "building_top_z_m", "axis_angle_deg",
    "side_bearing_deg", "min_surface_z_m", "penetration_m", "collision",
]
BUILDING_SUMMARY_COLUMNS = [
    "helipad_id", "building_id", "building_top_z_m", "collision_angle_count",
    "collision_axis_angles_deg", "collision_axis_ranges_deg", "max_penetration_m",
]
ANGLE_COLUMNS = [
    "helipad_id", "axis_angle_deg", "collision_count", "colliding_building_ids",
    "max_penetration_m", "is_collision",
]
MATCH_COLUMNS = [
    "helipad_id", "source_fid", "host_building_id", "host_source",
    "overlap_ratio", "helipad_z_m", "center_x_m", "center_y_m", "status",
]
UNMATCHED_COLUMNS = ["helipad_id", "source_fid", "reason", "candidate_host_ids"]
DIAGNOSTIC_COLUMNS = [
    "helipad_id", "axis_angle_deg", "status", "known_collision_count",
    "unknown_building_count", "unknown_building_ids",
]
INTERNAL_PREFIX = "__cd_"
LOG = logging.getLogger("collide_detection")


@dataclass(frozen=True)
class SurfaceConfig:
    length_m: float
    inner_width_m: float
    outer_width_m: float
    slope_percent: float
    start_offset_m: float
    roof_offset_m: float = HELIPAD_ROOF_OFFSET_M

    def __post_init__(self) -> None:
        if not all(math.isfinite(v) for v in asdict(self).values()):
            raise ValueError("접근면 수치는 모두 유한한 값이어야 합니다.")
        if min(self.length_m, self.inner_width_m, self.outer_width_m) <= 0:
            raise ValueError("접근면 길이와 폭은 0보다 커야 합니다.")
        if self.outer_width_m < self.inner_width_m:
            raise ValueError("끝 폭은 시작 폭 이상이어야 합니다.")
        if min(self.slope_percent, self.start_offset_m, self.roof_offset_m) < 0:
            raise ValueError("경사, 시작 거리, 지붕 이격 높이는 음수일 수 없습니다.")


@dataclass
class DemData:
    values: np.ndarray
    valid: np.ndarray
    geotransform: tuple[float, ...]
    crs: CRS
    backend: str


@dataclass(frozen=True)
class DemStat:
    elevation_m: float
    coverage_ratio: float
    valid_cell_count: int
    status: str


@dataclass
class Inputs:
    buildings: gpd.GeoDataFrame
    helipads: gpd.GeoDataFrame
    host_reference: gpd.GeoDataFrame | None
    dem: DemData


def validate_paths() -> None:
    """입출력 위치와 심볼릭 링크/정션의 실제 목적지를 검사한다."""
    root = PROJECT_ROOT.resolve()
    if Path(__file__).resolve() != SCRIPT_PATH.resolve():
        raise ValueError(f"프로그램 위치는 {SCRIPT_PATH}이어야 합니다.")
    for path, relative in (
        (SCRIPT_DIR, Path("코드") / "지형"),
        (OUTPUT_DIR, Path("자료") / "결과물"),
    ):
        if path.resolve() != root / relative:
            raise ValueError(f"허용되지 않은 실제 경로: {path.resolve()}")
    output_paths = [
        COLLISION_DETAIL_PATH, BUILDING_SUMMARY_PATH, ANGLE_SUMMARY_PATH,
        COLLISION_GPKG_PATH, UNMATCHED_HELIPADS_PATH, LOG_PATH, VALIDATION_PATH,
        HEIGHT_DIAGNOSTICS_PATH, HOST_MATCHES_PATH, ANGLE_DIAGNOSTICS_PATH,
    ]
    for path in output_paths:
        if path.resolve().parent != OUTPUT_DIR.resolve():
            raise ValueError(f"결과 경로가 OUTPUT_DIR를 벗어났습니다: {path}")
    if TEMP_DIR.resolve() != OUTPUT_DIR.resolve() / "_temp":
        raise ValueError("임시 경로가 OUTPUT_DIR/_temp를 벗어났습니다.")
    for path in (BUILDING_PATH, HELIPAD_PATH, DEM_PATH):
        if not path.is_file():
            raise FileNotFoundError(path)
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def input_fingerprints() -> dict[str, str]:
    """원본 데이터 폴더의 모든 일반 파일을 해시하여 쓰기 여부를 검증한다."""
    fingerprints = {}
    for path in sorted(INPUT_DIR.iterdir()):
        if path.is_file():
            with path.open("rb") as stream:
                fingerprints[path.name] = hashlib.file_digest(stream, "sha256").hexdigest()
    return fingerprints


@contextmanager
def output_workspace() -> Iterator[Path]:
    """이번 실행의 임시 파일만 정리하고, 외부 파일/기존 임시 파일은 보존한다."""
    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="collision_", dir=TEMP_DIR) as name:
        stage = Path(name)
        keys = ("TMP", "TEMP", "TMPDIR", "CPL_TMPDIR", "GDAL_PAM_ENABLED")
        previous = {key: os.environ.get(key) for key in keys}
        old_tempdir = tempfile.tempdir
        for key in keys[:-1]:
            os.environ[key] = str(stage)
        os.environ["GDAL_PAM_ENABLED"] = "NO"
        tempfile.tempdir = str(stage)
        try:
            yield stage
        finally:
            tempfile.tempdir = old_tempdir
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value
    try:
        TEMP_DIR.rmdir()  # 빈 폴더일 때만 제거
    except OSError:
        pass


def _read_vector(path: Path, prefix: str, id_field: str | None) -> gpd.GeoDataFrame:
    frame = gpd.read_file(path, engine="pyogrio", encoding=SHP_ENCODING, fid_as_index=True)
    if frame.crs is None:
        raise ValueError(f"좌표계가 없는 입력: {path}")
    if any(str(c).startswith(INTERNAL_PREFIX) for c in frame.columns):
        raise ValueError(f"예약된 내부 필드 접두사 {INTERNAL_PREFIX}: {path}")
    frame["__cd_fid"] = frame.index.astype("int64")
    if id_field is None:
        ids = pd.Series([f"{prefix}{fid:06d}" for fid in frame.index], index=frame.index)
    else:
        if id_field not in frame:
            raise ValueError(f"ID 필드 {id_field} 없음: {path}")
        if frame[id_field].isna().any():
            raise ValueError(f"ID 필드 {id_field}에 결측값이 있습니다.")
        ids = frame[id_field].astype(str).str.strip()
        if ids.eq("").any() or ids.duplicated().any():
            raise ValueError(f"ID 필드 {id_field}는 비어 있지 않은 고유값이어야 합니다.")
    frame["__cd_id"] = ids
    frame = frame.to_crs(ANALYSIS_CRS).reset_index(drop=True)
    return frame


class DemBackendUnavailable(RuntimeError):
    """데이터 오류와 구분하여 다른 DEM 읽기 환경으로 전환한다."""


def find_gdal_cli() -> tuple[Path, Path]:
    """같은 설치의 gdalinfo와 gdal_translate를 찾는다. pip GDAL 설치 불필요."""
    suffix = ".exe" if os.name == "nt" else ""
    directories = []
    if GDAL_BIN_DIR is not None:
        directories.append(GDAL_BIN_DIR)
    else:
        for name in ("gdalinfo", "gdal_translate"):
            found = shutil.which(name)
            if found:
                directories.append(Path(found).parent)
        qgis_bins = [path for root in GDAL_SEARCH_ROOTS for path in root.glob("QGIS */bin")]
        directories.extend(sorted(qgis_bins, key=lambda p: tuple(
            int(part) for part in re.findall(r"\d+", p.parent.name)), reverse=True))
        directories.extend(root / "bin" for root in GDAL_SEARCH_ROOTS)
    for directory in directories:
        info = directory / ("gdalinfo" + suffix)
        translate = directory / ("gdal_translate" + suffix)
        if info.is_file() and translate.is_file():
            return info.resolve(), translate.resolve()
    raise DemBackendUnavailable("QGIS/OSGeo4W의 gdalinfo, gdal_translate 실행 파일을 찾지 못했습니다.")


def _read_envi_plane(path: Path, shape: tuple[int, int], data_type: int) -> np.ndarray:
    """GDAL이 임시 생성한 단일 밴드 ENVI를 헤더의 byte order에 맞춰 읽는다."""
    header = path.with_suffix(".hdr").read_text(encoding="utf-8")
    fields = {key.lower(): int(value) for key, value in re.findall(
        r"^\s*(samples|lines|bands|header offset|data type|byte order)\s*=\s*(\d+)\s*$",
        header, re.MULTILINE | re.IGNORECASE,
    )}
    if (fields.get("lines"), fields.get("samples")) != shape or fields.get("bands") != 1:
        raise ValueError("GDAL 임시 래스터 크기/밴드 수가 원본과 다릅니다.")
    if fields.get("data type") != data_type or fields.get("byte order") not in (0, 1):
        raise ValueError("지원하지 않는 GDAL 임시 래스터 데이터 형식입니다.")
    endian = "<" if fields["byte order"] == 0 else ">"
    dtype = np.dtype(endian + {1: "u1", 5: "f8"}[data_type])
    offset = fields.get("header offset", 0)
    if path.stat().st_size != offset + math.prod(shape) * dtype.itemsize:
        raise ValueError("GDAL 임시 래스터의 파일 크기가 올바르지 않습니다.")
    return np.fromfile(path, dtype=dtype, offset=offset).reshape(shape)


def _load_dem_cli(path: Path) -> DemData:
    """Python ABI와 무관한 QGIS/GDAL 실행 파일로 원본을 읽기만 한다."""
    info_exe, translate_exe = find_gdal_cli()
    with output_workspace() as stage:
        env = os.environ.copy()
        # 외부 실행 파일의 GDAL/PROJ는 같은 설치의 데이터 디렉터리를 사용한다.
        install_root = info_exe.parent.parent
        for key, relative in (("GDAL_DATA", "share/gdal"), ("PROJ_DATA", "share/proj"),
                              ("PROJ_LIB", "share/proj")):
            location = install_root / relative
            if location.is_dir():
                env[key] = str(location)

        def run(arguments: list[str]) -> str:
            result = subprocess.run(
                arguments, capture_output=True, text=True, encoding="utf-8", errors="replace",
                env=env, cwd=stage, shell=False,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            if result.returncode:
                raise RuntimeError(f"GDAL 실행 실패 ({Path(arguments[0]).name}): {result.stderr.strip()}")
            return result.stdout

        source = str(path.resolve())
        metadata = json.loads(run([str(info_exe), "-json", source]))
        wkt = metadata.get("coordinateSystem", {}).get("wkt")
        if not wkt or "geoTransform" not in metadata:
            raise ValueError("DEM 좌표계 또는 affine geotransform이 없습니다.")
        shape = (metadata["size"][1], metadata["size"][0])
        band = metadata["bands"][0]
        arrays = []
        for selection, datatype, filename, envi_type in (
            ("1", "Float64", "dem_values.bsq", 5),
            ("mask,1", "Byte", "dem_mask.bsq", 1),
        ):
            destination = stage / filename
            run([str(translate_exe), "-q", "-of", "ENVI", "-ot", datatype,
                 "-b", selection, "-mask", "none", "-co", "INTERLEAVE=BSQ",
                 source, str(destination)])
            arrays.append(_read_envi_plane(destination, shape, envi_type))
        values, mask = arrays
        valid = (mask != 0) & np.isfinite(values)
        nodata = band.get("noDataValue")
        if nodata is not None:
            valid &= values != float(nodata)
        values = values * band.get("scale", 1.0) + band.get("offset", 0.0)
        return DemData(values, valid & np.isfinite(values), tuple(metadata["geoTransform"]),
                       CRS.from_wkt(wkt), "gdal-cli")


def _load_dem_rasterio(path: Path) -> DemData:
    try:
        import rasterio
    except (ImportError, OSError) as error:
        raise DemBackendUnavailable(f"Rasterio를 불러올 수 없습니다: {error}") from error
    with rasterio.Env(GDAL_PAM_ENABLED="NO"):
        with rasterio.open(path, "r") as dataset:
            if dataset.crs is None:
                raise ValueError("DEM 좌표계가 없습니다.")
            raw = dataset.read(1, masked=True)
            values = np.asarray(raw.data, dtype=float)
            valid = ~np.ma.getmaskarray(raw) & np.isfinite(values)
            if dataset.nodata is not None:
                valid &= values != dataset.nodata
            values = values * dataset.scales[0] + dataset.offsets[0]
            return DemData(values, valid & np.isfinite(values),
                           dataset.transform.to_gdal(), CRS(dataset.crs), "rasterio")


def _load_dem_gdal(path: Path) -> DemData:
    try:
        from osgeo import gdal, gdal_array  # NumPy bridge의 DLL도 사용 가능해야 한다.
    except (ImportError, OSError) as error:
        raise DemBackendUnavailable(f"기설치 GDAL Python 모듈을 불러올 수 없습니다: {error}") from error
    gdal.UseExceptions()
    gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
    dataset = gdal.OpenEx(str(path), gdal.OF_RASTER | gdal.OF_READONLY)
    if dataset is None or not dataset.GetProjection():
        raise ValueError("DEM을 열 수 없거나 좌표계가 없습니다.")
    band = dataset.GetRasterBand(1)
    values = band.ReadAsArray().astype(float)
    valid = (band.GetMaskBand().ReadAsArray() != 0) & np.isfinite(values)
    nodata = band.GetNoDataValue()
    if nodata is not None:
        valid &= values != nodata
    scale, offset = band.GetScale(), band.GetOffset()
    values = values * (1.0 if scale is None else scale) + (0.0 if offset is None else offset)
    result = DemData(values, valid & np.isfinite(values), dataset.GetGeoTransform(),
                     CRS.from_wkt(dataset.GetProjection()), "gdal")
    band = None
    dataset = None
    return result


def _load_dem(path: Path, backend: str = DEM_BACKEND) -> DemData:
    """미설치/DLL import 오류에만 대체 경로 사용. 손상된 데이터 오류는 그대로 알린다."""
    loaders = {"rasterio": _load_dem_rasterio, "gdal": _load_dem_gdal, "gdal-cli": _load_dem_cli}
    if backend != "auto" and backend not in loaders:
        raise ValueError(f"지원하지 않는 DEM backend: {backend}")
    failures = []
    for name, loader in loaders.items():
        if backend not in ("auto", name):
            continue
        try:
            result = loader(path)
        except DemBackendUnavailable as error:
            failures.append(str(error))
            LOG.info("DEM backend %s unavailable: %s", name, error)
        else:
            LOG.info("DEM backend: %s", result.backend)
            return result
    raise RuntimeError(
        "DEM 읽기 환경이 없습니다. pip install gdal은 필요하지 않습니다. "
        "python -m pip install --only-binary=:all: rasterio 로 설치하거나 "
        "QGIS를 설치하고 --dem-backend auto로 실행하세요. "
        "별도 설치 위치는 CONFIG의 GDAL_BIN_DIR에 지정할 수 있습니다.\n"
        + "\n".join(failures)
    )


def load_inputs(host_source: str = HOST_SOURCE, dem_backend: str = DEM_BACKEND) -> Inputs:
    crs = CRS(ANALYSIS_CRS)
    if not crs.is_projected or not all(
        math.isclose(axis.unit_conversion_factor, 1.0) for axis in crs.axis_info[:2]
    ):
        raise ValueError("분석 좌표계는 m 단위의 투영좌표계여야 합니다.")
    buildings = _read_vector(BUILDING_PATH, "B", BUILDING_ID_FIELD)
    helipads = _read_vector(HELIPAD_PATH, "H", HELIPAD_ID_FIELD)
    reference = None
    if host_source == "reference":
        reference = _read_vector(HOST_REFERENCE_PATH, "R", None)
    return Inputs(buildings, helipads, reference, _load_dem(DEM_PATH, dem_backend))


def inspect_inputs(inputs: Inputs) -> dict:
    report: dict = {"analysis_crs": ANALYSIS_CRS, "shp_encoding": SHP_ENCODING}
    for label, frame in (("buildings", inputs.buildings), ("helipads", inputs.helipads)):
        height = pd.to_numeric(frame.get(BUILDING_HEIGHT_FIELD, pd.Series(dtype=float)),
                               errors="coerce")
        text_columns = frame.select_dtypes(include=["object", "string"])
        corrupt = sum(int(text_columns[c].astype("string").str.contains(
            "\ufffd", regex=False, na=False).sum()) for c in text_columns)
        report[label] = {
            "active_feature_count": len(frame),
            "fields": [c for c in frame.columns if not c.startswith(INTERNAL_PREFIX)],
            "bounds": frame.total_bounds.tolist(),
            "geometry_types": frame.geom_type.value_counts().to_dict(),
            "invalid_geometry_count": int((~frame.geometry.is_valid).sum()),
            "duplicate_geometry_count": int(frame.geometry.to_wkb().duplicated().sum()),
            "missing_or_nonpositive_height_count": int((height.isna() | height.le(0)).sum()),
            "text_cells_with_existing_replacement_characters": corrupt,
        }
        LOG.info("%s: %d active features, duplicate geometries=%d", label, len(frame),
                 report[label]["duplicate_geometry_count"])
    report["dem"] = {
        "backend": inputs.dem.backend, "crs": inputs.dem.crs.to_string(),
        "shape": list(inputs.dem.values.shape),
        "geotransform": list(inputs.dem.geotransform),
        "valid_cell_count": int(inputs.dem.valid.sum()),
        "statistic": GROUND_STATISTIC, "minimum_coverage": MIN_DEM_COVERAGE,
    }
    return report


def _usable_polygon(geometry):
    if geometry is None or geometry.is_empty:
        return None
    geometry = shapely.force_2d(geometry)
    if not geometry.is_valid:
        geometry = shapely.make_valid(geometry)
    if geometry.geom_type in ("Polygon", "MultiPolygon"):
        return geometry if geometry.area > 0 else None
    if geometry.geom_type == "GeometryCollection":
        parts = [_usable_polygon(part) for part in geometry.geoms]
        parts = [part for part in parts if part is not None]
        return unary_union(parts) if parts else None
    return None


def get_dem_stat_for_geometry(
    geometry, dem: DemData, geometry_crs=ANALYSIS_CRS, statistic=GROUND_STATISTIC,
) -> DemStat:
    """겹치는 셀만 사용. NoData/범위 밖을 0으로 대체하거나 주변으로 보간하지 않는다.

    polygon의 면적과 실제 유효 셀의 교집합 면적으로 coverage를 구한다.
    평균/중앙값은 양의 면적으로 겹치는 셀에 대한 비가중 통계이다.
    회전 geotransform도 역변환하여 처리한다.
    """
    if statistic not in {"median", "mean", "min", "max"}:
        raise ValueError(f"지원하지 않는 DEM 통계: {statistic}")
    if geometry is None or geometry.is_empty:
        return DemStat(math.nan, 0.0, 0, "invalid_geometry")
    if CRS(geometry_crs) != dem.crs:
        transformer = Transformer.from_crs(geometry_crs, dem.crs, always_xy=True)
        geometry = transform(transformer.transform, geometry)
    gt = dem.geotransform
    inverse = np.linalg.inv(np.array([[gt[1], gt[2]], [gt[4], gt[5]]]))
    origin = -inverse @ np.array([gt[0], gt[3]])
    pixel_geometry = affine_transform(geometry, [*inverse[0], *inverse[1], *origin])
    rows, cols = dem.values.shape
    if pixel_geometry.geom_type == "Point":
        col, row = math.floor(pixel_geometry.x), math.floor(pixel_geometry.y)
        if 0 <= col < cols and 0 <= row < rows and dem.valid[row, col]:
            return DemStat(float(dem.values[row, col]), 1.0, 1, "ok")
        return DemStat(math.nan, 0.0, 0, "dem_nodata_or_outside")
    if pixel_geometry.area <= 0:
        return DemStat(math.nan, 0.0, 0, "invalid_geometry")
    xmin, ymin, xmax, ymax = pixel_geometry.bounds
    c0, c1 = max(0, math.floor(xmin)), min(cols, math.ceil(xmax))
    r0, r1 = max(0, math.floor(ymin)), min(rows, math.ceil(ymax))
    if c0 >= c1 or r0 >= r1:
        return DemStat(math.nan, 0.0, 0, "dem_outside")
    rr, cc = np.mgrid[r0:r1, c0:c1]
    cells = shapely.box(cc.ravel(), rr.ravel(), cc.ravel() + 1, rr.ravel() + 1)
    areas = shapely.area(shapely.intersection(cells, pixel_geometry))
    selected = (areas > 0) & dem.valid[r0:r1, c0:c1].ravel()
    coverage = min(1.0, float(areas[selected].sum() / pixel_geometry.area))
    count = int(selected.sum())
    if count == 0 or coverage < MIN_DEM_COVERAGE:
        return DemStat(math.nan, coverage, count, "dem_incomplete_coverage")
    values = dem.values[r0:r1, c0:c1].ravel()[selected]
    elevation = float(getattr(np, statistic)(values))
    return DemStat(elevation, coverage, count, "ok")


def prepare_building_heights(
    buildings: gpd.GeoDataFrame, dem: DemData, floor_height_m: float | None = FLOOR_HEIGHT_M,
) -> tuple[gpd.GeoDataFrame, pd.DataFrame]:
    if BUILDING_HEIGHT_FIELD not in buildings:
        raise ValueError(f"건물 높이 필드가 없습니다: {BUILDING_HEIGHT_FIELD}")
    if floor_height_m is not None and (not math.isfinite(floor_height_m) or floor_height_m <= 0):
        raise ValueError("층고는 0보다 큰 유한한 값이어야 합니다.")
    prepared = buildings.copy()
    prepared.geometry = prepared.geometry.map(_usable_polygon)
    if prepared.geometry.isna().any():
        ids = prepared.loc[prepared.geometry.isna(), "__cd_id"].tolist()
        raise ValueError(f"위치를 판정할 수 없는 건물 도형이 있습니다: {ids[:20]}")
    raw = pd.to_numeric(prepared[BUILDING_HEIGHT_FIELD], errors="coerce")
    agl = raw.where(np.isfinite(raw) & raw.gt(0))
    sources = pd.Series("attribute", index=prepared.index)
    sources.loc[agl.isna()] = "unknown"
    if floor_height_m is not None:
        if BUILDING_FLOORS_FIELD not in prepared:
            raise ValueError(f"층수 필드가 없습니다: {BUILDING_FLOORS_FIELD}")
        floors = pd.to_numeric(prepared[BUILDING_FLOORS_FIELD], errors="coerce")
        estimate = agl.isna() & np.isfinite(floors) & floors.gt(0)
        agl.loc[estimate] = floors.loc[estimate] * floor_height_m
        sources.loc[estimate] = "floor_estimate"
    dem_stats = []
    for number, geometry in enumerate(prepared.geometry, 1):
        dem_stats.append(get_dem_stat_for_geometry(geometry, dem))
        if number % 5000 == 0:
            LOG.info("DEM sampling: %d / %d", number, len(prepared))
    prepared["__cd_agl"] = agl
    prepared["__cd_ground"] = [s.elevation_m for s in dem_stats]
    prepared["__cd_top"] = prepared["__cd_ground"] + agl
    diagnostics = pd.DataFrame({
        "building_id": prepared["__cd_id"], "source_fid": prepared["__cd_fid"],
        "raw_height_m": raw, "height_source": sources,
        "building_agl_m": agl, "building_ground_z_m": prepared["__cd_ground"],
        "building_top_z_m": prepared["__cd_top"],
        "dem_coverage_ratio": [s.coverage_ratio for s in dem_stats],
        "dem_status": [s.status for s in dem_stats],
    })
    LOG.info("Height preparation: %d buildings, %d unknown tops", len(prepared),
             int(prepared["__cd_top"].isna().sum()))
    return prepared, diagnostics


def match_helipads_to_host_buildings(
    helipads: gpd.GeoDataFrame, buildings: gpd.GeoDataFrame, dem: DemData,
    host_source: str = HOST_SOURCE, reference: gpd.GeoDataFrame | None = None,
    floor_height_m: float | None = FLOOR_HEIGHT_M,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """면적 중첩으로 host 선택. 근처의 별개 건물을 거리만으로 매칭하지 않는다."""
    hosts = buildings
    if host_source == "reference":
        if reference is None:
            raise ValueError("host reference 입력이 없습니다.")
        # 27,000개 전체의 DEM을 다시 계산할 필요 없이 공간 후보만 준비한다.
        good = helipads.geometry.notna() & ~helipads.geometry.is_empty
        queries = reference.sindex.query(helipads.loc[good].geometry, predicate="intersects")
        hosts = reference.iloc[np.unique(queries[1])].copy().reset_index(drop=True)
        hosts, _ = prepare_building_heights(hosts, dem, floor_height_m)
    elif host_source == "helipad":
        hosts, _ = prepare_building_heights(helipads, dem, floor_height_m)
    elif host_source != "buildings":
        raise ValueError(f"지원하지 않는 host source: {host_source}")
    matches, unmatched = [], []
    for position, helipad in helipads.iterrows():
        geometry = helipad.geometry
        reason, candidates = "no_spatial_host", []
        if geometry is None or geometry.is_empty:
            reason = "invalid_helipad_geometry"
        else:
            geometry = shapely.force_2d(geometry)
            if not geometry.is_valid:
                geometry = shapely.make_valid(geometry)
            if host_source == "helipad":
                candidates = [(1.0, position)]
            else:
                for index in hosts.sindex.query(geometry, predicate="intersects"):
                    host_geom = hosts.iloc[index].geometry
                    if geometry.geom_type == "Point":
                        ratio = float(host_geom.covers(geometry))
                    elif geometry.area > 0:
                        ratio = host_geom.intersection(geometry).area / geometry.area
                    else:
                        continue
                    if ratio >= MIN_HOST_OVERLAP_RATIO:
                        candidates.append((ratio, int(index)))
            candidates.sort(key=lambda item: (-item[0], str(hosts.iloc[item[1]]["__cd_id"])))
        if candidates and len(candidates) > 1 and math.isclose(
            candidates[0][0], candidates[1][0], abs_tol=1e-9, rel_tol=0,
        ):
            reason = "ambiguous_spatial_hosts"
            selected = None
        else:
            selected = candidates[0] if candidates else None
        center = geometry.centroid if geometry is not None and not geometry.is_empty else None
        row = {
            "helipad_id": helipad["__cd_id"], "source_fid": helipad["__cd_fid"],
            "host_building_id": None, "host_source": host_source, "overlap_ratio": math.nan,
            "helipad_z_m": math.nan,
            "center_x_m": center.x if center is not None else math.nan,
            "center_y_m": center.y if center is not None else math.nan, "status": reason,
        }
        if selected is None:
            unmatched.append({
                "helipad_id": row["helipad_id"], "source_fid": row["source_fid"],
                "reason": reason,
                "candidate_host_ids": ",".join(str(hosts.iloc[i]["__cd_id"]) for _, i in candidates),
            })
        else:
            ratio, index = selected
            host = hosts.iloc[index]
            host_id = str(host["__cd_id"])
            # 다른 host 입력을 선택해도 동일한 도형의 장애물 행을 제외할 수 있게 ID 연결.
            if host_source != "buildings":
                same = [i for i in buildings.sindex.query(host.geometry, predicate="intersects")
                        if buildings.iloc[i].geometry.equals(host.geometry)]
                if len(same) == 1:
                    host_id = str(buildings.iloc[same[0]]["__cd_id"])
                elif host_source == "helipad":
                    host_id = "HOST_" + host_id
            row.update(host_building_id=host_id, overlap_ratio=ratio,
                       helipad_z_m=host["__cd_top"],
                       status="matched" if np.isfinite(host["__cd_top"]) else "host_height_unknown")
        matches.append(row)
    return pd.DataFrame(matches, columns=MATCH_COLUMNS), pd.DataFrame(unmatched, columns=UNMATCHED_COLUMNS)


def _bearing_vectors(bearing_deg: float) -> tuple[np.ndarray, np.ndarray]:
    radians = math.radians(bearing_deg % 360)
    direction = np.array([math.sin(radians), math.cos(radians)])
    right = np.array([math.cos(radians), -math.sin(radians)])
    return direction, right


def create_approach_polygon(center: Point, side_bearing_deg: float, config: SurfaceConfig) -> Polygon:
    direction, right = _bearing_vectors(side_bearing_deg)
    start = np.array([center.x, center.y]) + config.start_offset_m * direction
    end = start + config.length_m * direction
    return Polygon([
        start - right * config.inner_width_m / 2,
        start + right * config.inner_width_m / 2,
        end + right * config.outer_width_m / 2,
        end - right * config.outer_width_m / 2,
    ])


def extract_geometry_vertices(geometry) -> np.ndarray:
    """홀, MultiPolygon, 선, 점, GeometryCollection의 모든 XY 꼭짓점."""
    return shapely.get_coordinates(geometry, include_z=False)


def calculate_min_surface_z(
    intersection_geometry, center: Point, side_bearing_deg: float,
    helipad_z_m: float, config: SurfaceConfig,
) -> float:
    """z = z_roof + offset_z + (s - start_offset) * slope / 100.

    선형 면의 최소값은 실제 교집합 꼭짓점에서 발생한다. 건물 중심점이나
    원래 건물 꼭짓점만 검사하면 놓치는, 면 경계의 새 교차점도 포함한다.
    """
    vertices = extract_geometry_vertices(intersection_geometry)
    if len(vertices) == 0:
        return math.nan
    direction, _ = _bearing_vectors(side_bearing_deg)
    along = (vertices - np.array([center.x, center.y])) @ direction - config.start_offset_m
    distances = np.clip(along, 0.0, config.length_m)
    return float(helipad_z_m + config.roof_offset_m + distances.min() * config.slope_percent / 100)


def check_surface_collisions(
    buildings: gpd.GeoDataFrame, matches: pd.DataFrame,
    config: SurfaceConfig | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """이벤트는 (helipad, axis, side, building), 각도 건물 수는 양쪽을 합친 고유 수."""
    details, angles, diagnostics = [], [], []
    spatial_index = buildings.sindex
    building_ids = buildings["__cd_id"].to_numpy(dtype=str)
    tops = buildings["__cd_top"].to_numpy(dtype=float)
    agls = buildings["__cd_agl"].to_numpy(dtype=float)
    grounds = buildings["__cd_ground"].to_numpy(dtype=float)
    geometries = buildings.geometry.to_numpy()
    for number, match in enumerate(matches.to_dict("records"), 1):
        helipad_id = match["helipad_id"]
        unavailable = match["status"] if match["status"] != "matched" else None
        if config is None:
            unavailable = "surface_parameters_missing" if unavailable is None else unavailable + ";surface_parameters_missing"
        for axis in AXIS_ANGLES:
            collisions: dict[str, float] = {}
            unknown: set[str] = set()
            if unavailable is None:
                center = Point(match["center_x_m"], match["center_y_m"])
                for bearing in (axis, axis + 180):
                    surface = create_approach_polygon(center, bearing, config)
                    candidates = np.sort(spatial_index.query(surface, predicate="intersects"))
                    candidates = candidates[building_ids[candidates] != match["host_building_id"]]
                    known = np.isfinite(tops[candidates])
                    unknown.update(building_ids[candidates[~known]].tolist())
                    candidates = candidates[known]
                    # 비음수 경사 면의 시작 높이 이하 건물은 벡터 연산으로 제외한다.
                    candidates = candidates[tops[candidates] > (
                        match["helipad_z_m"] + config.roof_offset_m + PENETRATION_TOLERANCE_M)]
                    for index in candidates:
                        building_id = building_ids[index]
                        intersection = geometries[index].intersection(surface)
                        if intersection.is_empty:
                            continue
                        minimum = calculate_min_surface_z(intersection, center, bearing,
                                                          match["helipad_z_m"], config)
                        penetration = float(tops[index] - minimum)
                        if penetration > PENETRATION_TOLERANCE_M:
                            details.append({
                                "helipad_id": helipad_id, "host_building_id": match["host_building_id"],
                                "building_id": building_id, "building_agl_m": agls[index],
                                "building_ground_z_m": grounds[index],
                                "building_top_z_m": tops[index], "axis_angle_deg": axis,
                                "side_bearing_deg": bearing, "min_surface_z_m": minimum,
                                "penetration_m": penetration, "collision": True,
                            })
                            collisions[building_id] = max(collisions.get(building_id, 0.0), penetration)
            incomplete = unavailable is not None or bool(unknown)
            status = unavailable or ("unknown_obstacle_heights" if unknown else "evaluated")
            angles.append({
                "helipad_id": helipad_id, "axis_angle_deg": axis,
                "collision_count": pd.NA if incomplete else len(collisions),
                "colliding_building_ids": ",".join(sorted(collisions)),
                "max_penetration_m": max(collisions.values()) if collisions else (math.nan if incomplete else 0.0),
                "is_collision": True if collisions else (pd.NA if incomplete else False),
            })
            diagnostics.append({
                "helipad_id": helipad_id, "axis_angle_deg": axis, "status": status,
                "known_collision_count": len(collisions), "unknown_building_count": len(unknown),
                "unknown_building_ids": ",".join(sorted(unknown)),
            })
        LOG.info("Helipad %d/%d %s: %s", number, len(matches), helipad_id, unavailable or "checked")
    angle_frame = pd.DataFrame(angles, columns=ANGLE_COLUMNS)
    angle_frame["collision_count"] = angle_frame["collision_count"].astype("Int64")
    angle_frame["is_collision"] = angle_frame["is_collision"].astype("boolean")
    return (pd.DataFrame(details, columns=DETAIL_COLUMNS), angle_frame,
            pd.DataFrame(diagnostics, columns=DIAGNOSTIC_COLUMNS))


def merge_consecutive_angles(angles) -> str:
    """연속 각도를 오름차순으로 병합한다. 179와 0은 별도 구간으로 표시한다."""
    values = sorted({int(angle) for angle in angles})
    if not values:
        return ""
    spans, start, previous = [], values[0], values[0]
    for value in values[1:]:
        if value != previous + 1:
            spans.append(str(start) if start == previous else f"{start}-{previous}")
            start = value
        previous = value
    spans.append(str(start) if start == previous else f"{start}-{previous}")
    return "; ".join(spans)


def _building_summary(details: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (helipad_id, building_id), group in details.groupby(["helipad_id", "building_id"], sort=True):
        axes = sorted(set(int(value) for value in group.axis_angle_deg))
        rows.append({
            "helipad_id": helipad_id, "building_id": building_id,
            "building_top_z_m": group.building_top_z_m.iloc[0],
            "collision_angle_count": len(axes),
            "collision_axis_angles_deg": ",".join(map(str, axes)),
            "collision_axis_ranges_deg": merge_consecutive_angles(axes),
            "max_penetration_m": group.penetration_m.max(),
        })
    return pd.DataFrame(rows, columns=BUILDING_SUMMARY_COLUMNS)


def _collision_geodataframe(original: gpd.GeoDataFrame, details: pd.DataFrame) -> gpd.GeoDataFrame:
    selected = original.loc[original["__cd_id"].isin(details.building_id)].copy()
    source_ids = selected["__cd_id"].copy()
    source_fids = selected["__cd_fid"].copy()
    selected = selected.drop(columns=[c for c in selected if c.startswith(INTERNAL_PREFIX)])
    additions = ["building_id", "source_fid", "collision_count", "helipad_count", "max_pen_m",
                 "collision_axis_angles_deg", "collision_axis_ranges_deg", "helipad_axis_ranges"]
    for column in additions:
        if column in selected:
            renamed = "source_" + column
            while renamed in selected or renamed in additions:
                renamed = "source_" + renamed
            selected = selected.rename(columns={column: renamed})
    selected["building_id"] = source_ids
    selected["source_fid"] = source_fids.astype("int64")
    grouped = {str(key): group for key, group in details.groupby("building_id")}
    selected["collision_count"] = pd.Series(
        [len(grouped[key]) for key in source_ids], index=selected.index, dtype="int64")
    selected["helipad_count"] = pd.Series(
        [grouped[key].helipad_id.nunique() for key in source_ids], index=selected.index, dtype="int64")
    selected["max_pen_m"] = pd.Series(
        [grouped[key].penetration_m.max() for key in source_ids], index=selected.index, dtype="float64")
    selected["collision_axis_angles_deg"] = pd.Series([
        ",".join(map(str, sorted(set(grouped[key].axis_angle_deg)))) for key in source_ids
    ], index=selected.index, dtype="object")
    selected["collision_axis_ranges_deg"] = pd.Series([
        merge_consecutive_angles(grouped[key].axis_angle_deg) for key in source_ids
    ], index=selected.index, dtype="object")
    selected["helipad_axis_ranges"] = pd.Series([
        json.dumps({str(helipad): merge_consecutive_angles(group.axis_angle_deg)
                    for helipad, group in grouped[key].groupby("helipad_id")}, ensure_ascii=False)
        for key in source_ids
    ], index=selected.index, dtype="object")
    # 원본 도형/속성 보존. 입력 FID는 source_fid 필드로 저장하며 GPKG FID와 분리한다.
    return selected


def save_outputs(
    stage: Path, inputs: Inputs, details: pd.DataFrame, angles: pd.DataFrame,
    unmatched: pd.DataFrame, matches: pd.DataFrame, height_diagnostics: pd.DataFrame,
    angle_diagnostics: pd.DataFrame, report: dict, before_hashes: dict[str, str],
) -> None:
    summary = _building_summary(details)
    tables = {
        COLLISION_DETAIL_PATH: details, BUILDING_SUMMARY_PATH: summary,
        ANGLE_SUMMARY_PATH: angles, UNMATCHED_HELIPADS_PATH: unmatched,
        HEIGHT_DIAGNOSTICS_PATH: height_diagnostics, HOST_MATCHES_PATH: matches,
        ANGLE_DIAGNOSTICS_PATH: angle_diagnostics,
    }
    expected = len(inputs.helipads) * len(AXIS_ANGLES)
    if len(angles) != expected or angles.duplicated(["helipad_id", "axis_angle_deg"]).any():
        raise RuntimeError("각도표의 개수/고유성 검증 실패")
    for _, group in angles.groupby("helipad_id"):
        if set(group.axis_angle_deg) != set(AXIS_ANGLES):
            raise RuntimeError("누락된 축 각도가 있습니다.")
    if details.duplicated(["helipad_id", "axis_angle_deg", "side_bearing_deg", "building_id"]).any():
        raise RuntimeError("중복 collision event")
    for path, frame in tables.items():
        frame.to_csv(stage / path.name, index=False, encoding="utf-8-sig")
    colliding = _collision_geodataframe(inputs.buildings, details)
    gpkg_stage = stage / COLLISION_GPKG_PATH.name
    pyogrio.write_dataframe(colliding, gpkg_stage, layer="collision_buildings", driver="GPKG",
                            geometry_type="MultiPolygon", promote_to_multi=True)
    reread = pyogrio.read_dataframe(gpkg_stage, layer="collision_buildings")
    if len(reread) != details.building_id.nunique():
        raise RuntimeError("GPKG 재읽기/피처 수 검증 실패")
    after_hashes = input_fingerprints()
    if before_hashes != after_hashes:
        raise RuntimeError("실행 중 입력 폴더 파일이 변경되었습니다. 결과는 게시하지 않았습니다.")
    report["validation"] = {
        "input_files_unchanged": True, "input_sha256": before_hashes,
        "angle_rows": len(angles), "expected_angle_rows": expected,
        "collision_events": len(details), "colliding_buildings": len(colliding),
        "unmatched_helipads": len(unmatched),
        "unevaluated_or_partial_angles": int(angle_diagnostics.status.ne("evaluated").sum()),
        "gpkg_reopen_verified": True,
        "output_paths": [str(path) for path in [*tables, COLLISION_GPKG_PATH, VALIDATION_PATH, LOG_PATH]],
    }
    report["definitions"] = {
        "axis": "0..179, grid north=0, clockwise, bearings=(axis, axis+180)",
        "helipad_center": "input geometry centroid in analysis CRS; Point remains unchanged",
        "building_top": f"DEM {GROUND_STATISTIC} + positive {BUILDING_HEIGHT_FIELD} (metres)",
        "helipad_z": "host building top + explicitly configured roof offset",
        "collision": f"building_top - minimum_surface_z > {PENETRATION_TOLERANCE_M} m; XY boundary contact included",
        "detail": "positive collision events only; one row per helipad/axis/side/building",
        "angle_collision_count": "unique buildings across both sides; blank if count is incomplete",
        "is_collision": "True=known collision; False=fully evaluated with none; blank=undetermined",
        "max_penetration": "maximum known penetration; a lower bound if some heights are unknown",
        "gpkg_collision_count": "number of detail events; two sides may count separately",
        "ids": "original SHP FID-based unless CONFIG specifies unique non-null ID fields",
        "duplicates": "input helipad features retained separately, including identical geometries",
        "no_data": "zero/missing AGL and insufficient DEM coverage remain unknown",
        "text_encoding": "CSV UTF-8-SIG; already-corrupted source strings are preserved",
        "scope": "linear analytical surface only; no implicit statutory standard or vertical datum conversion",
    }
    with (stage / VALIDATION_PATH.name).open("w", encoding="utf-8") as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2, default=str, allow_nan=False)
        stream.write("\n")
    # 지정된 결과 파일만 개별 교체한다. 입력 폴더에는 쓰지 않는다.
    for target in [*tables, COLLISION_GPKG_PATH, VALIDATION_PATH]:
        (stage / target.name).replace(target)
    LOG.info("Outputs verified and saved: %s", OUTPUT_DIR)


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--length-m", type=float, default=APPROACH_LENGTH_M)
    parser.add_argument("--inner-width-m", type=float, default=APPROACH_INNER_WIDTH_M)
    parser.add_argument("--outer-width-m", type=float, default=APPROACH_OUTER_WIDTH_M)
    parser.add_argument("--slope-percent", type=float, default=APPROACH_SLOPE_PERCENT,
                        help="5는 경사 5%%(0.05)를 의미")
    parser.add_argument("--start-offset-m", type=float, default=APPROACH_START_OFFSET_M,
                        help="헬리포트 도형 중심에서 접근면 시작점까지 거리")
    parser.add_argument("--roof-offset-m", type=float, default=HELIPAD_ROOF_OFFSET_M)
    parser.add_argument("--host-source", choices=("buildings", "reference", "helipad"), default=HOST_SOURCE)
    parser.add_argument("--floor-height-m", type=float, default=FLOOR_HEIGHT_M,
                        help="명시한 경우에만 0/결측 AGL을 A26 × 층고로 추정")
    parser.add_argument("--dem-backend", choices=("auto", "rasterio", "gdal", "gdal-cli"),
                        default=DEM_BACKEND,
                        help="auto: Rasterio → 기설치 osgeo → QGIS GDAL 실행 파일. pip gdal 설치 불필요")
    parser.add_argument("--inspect-only", action="store_true",
                        help="충돌 계산을 생략하고 판정 불가 각도표/입력 검증 산출물을 생성")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    validate_paths()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8"),
                                  logging.StreamHandler(sys.stdout)], force=True)
    values = (args.length_m, args.inner_width_m, args.outer_width_m, args.slope_percent, args.start_offset_m)
    surface = None
    if all(value is not None for value in values):
        surface = SurfaceConfig(*values, roof_offset_m=args.roof_offset_m)
    missing = [name for name, value in zip(
        ("length_m", "inner_width_m", "outer_width_m", "slope_percent", "start_offset_m"), values,
    ) if value is None]
    before_hashes = input_fingerprints()
    with output_workspace() as stage:
        inputs = load_inputs(args.host_source, args.dem_backend)
        report = inspect_inputs(inputs)
        report["surface_config"] = asdict(surface) if surface is not None else None
        report["missing_surface_parameters"] = missing
        report["host_source"] = args.host_source
        report["floor_height_m"] = args.floor_height_m
        report["inspect_only"] = args.inspect_only
        buildings, height_diagnostics = prepare_building_heights(inputs.buildings, inputs.dem, args.floor_height_m)
        report["height_validation"] = {
            "unknown_top_count": int(height_diagnostics.building_top_z_m.isna().sum()),
            "unknown_agl_count": int(height_diagnostics.building_agl_m.isna().sum()),
            "unknown_ground_count": int(height_diagnostics.building_ground_z_m.isna().sum()),
            "estimated_agl_count": int(height_diagnostics.height_source.eq("floor_estimate").sum()),
        }
        matches, unmatched = match_helipads_to_host_buildings(
            inputs.helipads, buildings, inputs.dem, args.host_source, inputs.host_reference, args.floor_height_m)
        details, angles, diagnostics = check_surface_collisions(
            buildings, matches, None if args.inspect_only else surface)
        if args.inspect_only:
            diagnostics["status"] = "inspection_only;" + diagnostics["status"]
        incomplete = bool(missing or args.inspect_only or diagnostics.status.ne("evaluated").any())
        report["status"] = "INCOMPLETE" if incomplete else "COMPLETE"
        save_outputs(stage, inputs, details, angles, unmatched, matches, height_diagnostics,
                     diagnostics, report, before_hashes)
    if incomplete:
        LOG.warning("INCOMPLETE: blank results mean undetermined, not collision-free. See %s", VALIDATION_PATH)
        if missing:
            LOG.warning("Missing surface parameters: %s", ", ".join(missing))
        return 2
    LOG.info("COMPLETE")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        LOG.exception("분석 실패. 입력 원본은 읽기 전용으로 취급됩니다.")
        raise SystemExit(1)
