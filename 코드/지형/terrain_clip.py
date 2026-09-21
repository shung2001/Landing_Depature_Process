"""송파구 행정경계로 한반도 90 m DEM을 잘라 GeoTIFF로 저장한다.

기본 실행:
    python 코드/지형/terrain_clip.py

다른 파일을 사용할 때:
    python 코드/지형/terrain_clip.py --dem INPUT.img --cutline BOUNDARY.shp \
        --output OUTPUT.tif
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DEM = Path(
    r"C:\Users\choih\Desktop\연구\전체자료\한반도_DEM\한반도90m_GRS80.img"
)
DEFAULT_CUTLINE = PROJECT_ROOT / "자료" / "안전성" / "송파구" / "송파구_지도.shp"
DEFAULT_OUTPUT = PROJECT_ROOT / "자료" / "안전성" / "송파구" / "송파구_DEM.tif"


def qgis_version(path: Path) -> tuple[int, ...]:
    """QGIS 설치 경로에서 버전 숫자를 추출한다."""
    match = re.search(r"QGIS\s+([0-9.]+)", str(path), flags=re.IGNORECASE)
    return tuple(int(part) for part in match.group(1).split(".")) if match else ()


def find_gdalwarp(explicit_path: Path | None = None) -> Path:
    """PATH와 일반적인 QGIS/OSGeo4W 설치 위치에서 gdalwarp를 찾는다."""
    if explicit_path is not None:
        path = explicit_path.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"gdalwarp 실행 파일이 없습니다: {path}")
        return path

    from_path = shutil.which("gdalwarp")
    if from_path:
        return Path(from_path).resolve()

    candidates: list[Path] = []
    program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
    candidates.extend(program_files.glob("QGIS */bin/gdalwarp.exe"))
    candidates.extend(Path(r"C:\OSGeo4W").glob("bin/gdalwarp.exe"))
    candidates.extend(Path(r"C:\OSGeo4W64").glob("bin/gdalwarp.exe"))

    existing = [path.resolve() for path in candidates if path.is_file()]
    if existing:
        return max(existing, key=lambda path: (qgis_version(path), path.stat().st_mtime))

    raise FileNotFoundError(
        "gdalwarp를 찾을 수 없습니다. QGIS/GDAL을 설치하거나 "
        "--gdalwarp 옵션으로 실행 파일 경로를 지정하세요."
    )


def build_command(
    gdalwarp: Path,
    dem: Path,
    cutline: Path,
    output: Path,
    target_crs: str,
    resolution: float,
    resampling: str,
    nodata: float,
) -> list[str]:
    """출력 형식에 맞는 gdalwarp 명령을 구성한다."""
    common_options = [
        str(gdalwarp),
        "-overwrite",
        "-cutline",
        str(cutline),
        "-crop_to_cutline",
        "-t_srs",
        target_crs,
        "-tr",
        str(resolution),
        str(resolution),
        "-tap",
        "-r",
        resampling,
        "-dstnodata",
        str(nodata),
        "-multi",
        "-wo",
        "NUM_THREADS=ALL_CPUS",
    ]

    suffix = output.suffix.lower()
    if suffix in {".tif", ".tiff"}:
        format_options = [
            "-of",
            "GTiff",
            "-co",
            "COMPRESS=DEFLATE",
            "-co",
            "PREDICTOR=3",
            "-co",
            "TILED=YES",
            "-co",
            "BIGTIFF=IF_SAFER",
            "-co",
            "NUM_THREADS=ALL_CPUS",
        ]
    elif suffix == ".img":
        format_options = ["-of", "HFA", "-co", "COMPRESSED=YES"]
    else:
        raise ValueError("출력 확장자는 .tif, .tiff 또는 .img여야 합니다.")

    return common_options + format_options + [str(dem), str(output)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="행정경계 Shapefile로 DEM을 잘라 GeoTIFF 또는 IMG로 저장합니다."
    )
    parser.add_argument("--dem", type=Path, default=DEFAULT_DEM, help="입력 DEM 경로")
    parser.add_argument(
        "--cutline", type=Path, default=DEFAULT_CUTLINE, help="클립 경계 Shapefile 경로"
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT, help="출력 .tif 또는 .img 경로"
    )
    parser.add_argument(
        "--target-crs",
        default="EPSG:5186",
        help="출력 좌표계(기본값: 송파구 경계와 같은 EPSG:5186)",
    )
    parser.add_argument(
        "--resampling",
        choices=("near", "bilinear", "cubic", "cubicspline", "lanczos"),
        default="bilinear",
        help="DEM 재표본화 방식(기본값: bilinear)",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=90.0,
        help="출력 픽셀 크기(m, 기본값: 90)",
    )
    parser.add_argument(
        "--nodata", type=float, default=-9999.0, help="경계 밖 NoData 값(기본값: -9999)"
    )
    parser.add_argument(
        "--gdalwarp", type=Path, help="gdalwarp.exe 경로(생략하면 자동 탐색)"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    dem = args.dem.expanduser().resolve()
    cutline = args.cutline.expanduser().resolve()
    output = args.output.expanduser().resolve()

    for label, path in (("입력 DEM", dem), ("클립 경계", cutline)):
        if not path.is_file():
            raise FileNotFoundError(f"{label} 파일이 없습니다: {path}")
    if args.resolution <= 0:
        raise ValueError("출력 픽셀 크기는 0보다 커야 합니다.")

    output.parent.mkdir(parents=True, exist_ok=True)
    gdalwarp = find_gdalwarp(args.gdalwarp)
    command = build_command(
        gdalwarp=gdalwarp,
        dem=dem,
        cutline=cutline,
        output=output,
        target_crs=args.target_crs,
        resolution=args.resolution,
        resampling=args.resampling,
        nodata=args.nodata,
    )

    print(f"GDAL: {gdalwarp}")
    print(f"입력: {dem}")
    print(f"경계: {cutline}")
    print(f"출력: {output}")
    subprocess.run(command, check=True)

    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"출력 파일이 정상적으로 생성되지 않았습니다: {output}")

    print(f"완료: {output} ({output.stat().st_size / 1024:.1f} KiB)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError, RuntimeError) as error:
        print(f"오류: {error}", file=sys.stderr)
        raise SystemExit(1) from error
    except subprocess.CalledProcessError as error:
        print(f"오류: gdalwarp 실행 실패(종료 코드 {error.returncode})", file=sys.stderr)
        raise SystemExit(error.returncode) from error
