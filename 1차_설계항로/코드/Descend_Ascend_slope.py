"""Approach_z.shp의 항로별 고도 프로파일과 구간별 Z축 경사를 그린다.

경사(도) = degrees(atan2(다음 정점 Z - 현재 정점 Z, 두 정점의 XY 수평거리)).
음수는 하강, 양수는 상승이다. 실행: python Descend_Ascend_slope.py
필요 패키지: geopandas, numpy, matplotlib.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "결과물" / "Outer" / "vertical_clearance_alternative"
INPUT_SHP = DATA_DIR / "Approach_z.shp"
OUTPUT_PNG = DATA_DIR / "각_항로별_slope.png"


def configure_font() -> None:
    # Windows의 한글 표시를 우선하고 다른 환경에서는 설치된 한글 글꼴을 사용한다.
    available = {font.name for font in matplotlib.font_manager.fontManager.ttflist}
    for name in ("Malgun Gothic", "NanumGothic", "Noto Sans CJK KR", "AppleGothic"):
        if name in available:
            plt.rcParams["font.family"] = name
            break
    plt.rcParams["axes.unicode_minus"] = False


def route_profile(geometry):
    if geometry is None or geometry.is_empty or geometry.geom_type != "LineString":
        raise ValueError("모든 피처는 비어 있지 않은 3차원 LineString이어야 합니다.")
    if not geometry.has_z:
        raise ValueError("LineString에 Z 좌표가 없습니다.")

    xyz = np.asarray(geometry.coords, dtype=float)
    if len(xyz) < 2 or not np.isfinite(xyz).all():
        raise ValueError("정점이 부족하거나 좌표에 유효하지 않은 값이 있습니다.")

    xy_length_m = np.hypot(np.diff(xyz[:, 0]), np.diff(xyz[:, 1]))
    if np.any(xy_length_m <= 0):
        raise ValueError("수평거리가 0인 구간이 있어 경사를 계산할 수 없습니다.")

    distance_km = np.r_[0, np.cumsum(xy_length_m)] / 1000
    slope_deg = np.degrees(np.arctan2(np.diff(xyz[:, 2]), xy_length_m))
    return distance_km, xyz[:, 2], slope_deg


def main() -> None:
    if not INPUT_SHP.is_file():
        raise FileNotFoundError(INPUT_SHP)
    routes = gpd.read_file(INPUT_SHP)
    if routes.empty:
        raise ValueError("입력 SHP에 항로가 없습니다.")
    if routes.crs is None or not routes.crs.is_projected:
        raise ValueError("수평거리 계산을 위해 투영 좌표계가 필요합니다.")
    if not np.isclose(routes.crs.axis_info[0].unit_conversion_factor, 1.0):
        raise ValueError("XY 좌표 단위가 미터인 투영 좌표계가 필요합니다.")

    configure_font()
    fig, axes = plt.subplots(len(routes), 2, figsize=(15, max(3.2 * len(routes), 5)),
                             squeeze=False, constrained_layout=True)
    fig.suptitle("접근 항로별 고도 및 Z축 경사", fontsize=17, fontweight="bold")

    for row_no, (_, route) in enumerate(routes.iterrows()):
        distance_km, z_m, slope_deg = route_profile(route.geometry)
        pair = route.get("pair", row_no + 1)
        cand_id = route.get("cand_id", "-")
        label = f"항로 {pair} (cand_id {cand_id})"
        altitude_ax, slope_ax = axes[row_no]

        altitude_ax.plot(distance_km, z_m, color="#176b9a", linewidth=2)
        altitude_ax.scatter([distance_km[0], distance_km[-1]], [z_m[0], z_m[-1]],
                            color="#176b9a", s=20, zorder=3)
        altitude_ax.set_title(f"{label} · 고도 프로파일", loc="left", fontsize=11)
        altitude_ax.set_ylabel("고도 Z (m)")
        altitude_ax.set_xlim(0, distance_km[-1])
        altitude_ax.grid(alpha=0.25)

        # 각 값은 해당 XY 정점 구간 전체의 경사다. 계단형으로 구간 범위를 보여준다.
        slope_ax.stairs(slope_deg, distance_km, color="#c45236", linewidth=1.5)
        slope_ax.axhline(0, color="#555555", linewidth=0.9)
        slope_ax.set_title(
            f"구간별 경사 · 최저 {slope_deg.min():+.2f}° / 최고 {slope_deg.max():+.2f}°",
            loc="left", fontsize=11,
        )
        slope_ax.set_ylabel("경사 (°)  -하강 / +상승")
        slope_ax.set_xlim(0, distance_km[-1])
        slope_ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        slope_ax.grid(alpha=0.25)

        for ax in (altitude_ax, slope_ax):
            ax.set_xlabel("시작점부터 누적 수평거리 (km)")

        print(f"{label}: {distance_km[-1]:.3f} km, "
              f"Z {z_m[0]:.2f} → {z_m[-1]:.2f} m, "
              f"구간 경사 {slope_deg.min():+.3f}° ~ {slope_deg.max():+.3f}°")

    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=180, facecolor="white")
    plt.close(fig)
    print(f"저장 완료: {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
