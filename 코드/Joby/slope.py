"""Joby S4 비행 궤적 CSV에 구간별 경사각을 추가한다.

각 행의 ``slope``는 이전 행에서 현재 행까지의 비행 구간에 대한
부호 있는 경사각(degree)이다. 위도·경도로 수평 대권거리를 구하고,
피트 단위인 고도 차를 미터로 변환한 뒤 ``atan2(고도 차, 수평거리)``를
적용한다. 따라서 상승 구간은 양수, 하강 구간은 음수로 기록된다.

첫 행은 이전 위치가 없어 경사각을 계산할 수 없으므로 빈 값(NaN)으로
저장한다.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE_CSV = Path(
    r"C:\Users\c\Desktop\대학생활\학연생\새로운_경로방식"
    r"\자료\심화\교수님_자료\새로운_자료\자료"
    r"\260501_Joby_S4\260501_Joby_S4\3f7e9f84_kepler.csv"
)
PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_CSV = PROJECT_ROOT / "기초자료" / "Joby_S4_Slope_5.csv"

LATITUDE_COLUMN = "latitude"
LONGITUDE_COLUMN = "longitude"
ALTITUDE_COLUMN = "Altitude"
SLOPE_COLUMN = "slope"

EARTH_RADIUS_M = 6_371_008.8
FEET_TO_METRES = 0.3048


def calculate_slope_degrees(data: pd.DataFrame) -> pd.Series:
    """이전 행부터 현재 행까지의 부호 있는 경사각을 계산한다.

    수평거리는 WGS84 위·경도에 Haversine 공식을 적용한 대권거리로
    근사한다. 입력값이 비어 있거나 숫자로 변환되지 않는 구간은 NaN으로
    남긴다. 수평 이동이 0인 구간은 ``atan2`` 정의에 따라 고도 변화가
    있으면 ±90도, 고도 변화도 없으면 0도가 된다.
    """

    required_columns = {
        LATITUDE_COLUMN,
        LONGITUDE_COLUMN,
        ALTITUDE_COLUMN,
    }
    missing_columns = required_columns.difference(data.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"원본 CSV에 필수 컬럼이 없습니다: {missing}")

    latitude = pd.to_numeric(data[LATITUDE_COLUMN], errors="coerce") #pd.to_numeric(): 자료형을 실수형으로 변환하는 것
    longitude = pd.to_numeric(data[LONGITUDE_COLUMN], errors="coerce")
    altitude_ft = pd.to_numeric(data[ALTITUDE_COLUMN], errors="coerce")

    invalid_latitude = latitude.notna() & ~latitude.between(-90.0, 90.0) #notna(): 결측값의 경우 'False'로 출력. between(): 범위 이외의 수치들에 대해서 'True'를 출력
    invalid_longitude = longitude.notna() & ~longitude.between(-180.0, 180.0)
    if invalid_latitude.any() or invalid_longitude.any():
        raise ValueError("위도 또는 경도 값이 유효 범위를 벗어났습니다.")

    latitude_rad = np.radians(latitude)
    longitude_rad = np.radians(longitude)

    delta_latitude = latitude_rad - latitude_rad.shift(1)
    delta_longitude = longitude_rad - longitude_rad.shift(1) # 한 행 아래로 이동

    haversine_a = (
        np.sin(delta_latitude / 2.0) ** 2
        + np.cos(latitude_rad.shift(1))
        * np.cos(latitude_rad)
        * np.sin(delta_longitude / 2.0) ** 2
    )
    # 부동소수점 반올림으로 a가 [0, 1]을 아주 조금 벗어나는 경우를 방지한다.
    haversine_a = haversine_a.clip(lower=0.0, upper=1.0)
    central_angle = 2.0 * np.arctan2(
        np.sqrt(haversine_a), np.sqrt(1.0 - haversine_a)
    )
    horizontal_distance_m = EARTH_RADIUS_M * central_angle

    altitude_change_m = altitude_ft.diff() * FEET_TO_METRES
    slope = np.degrees(np.arctan2(altitude_change_m, horizontal_distance_m))

    # 원본 DataFrame의 인덱스를 그대로 보존하고 첫 행은 계산 불가로 표시한다.
    result = pd.Series(slope, index=data.index, dtype="float64", name=SLOPE_COLUMN)
    if not result.empty:
        result.iloc[0] = math.nan
    return result


def add_slope_column(source: Path = SOURCE_CSV, output: Path = OUTPUT_CSV) -> Path:
    """원본 CSV를 읽어 slope 컬럼을 추가하고 UTF-8 BOM CSV로 저장한다."""

    if not source.is_file():
        raise FileNotFoundError(f"원본 CSV를 찾을 수 없습니다: {source}")

    data = pd.read_csv(source, encoding="utf-8-sig")
    if data.empty:
        raise ValueError(f"원본 CSV에 데이터 행이 없습니다: {source}")

    data[SLOPE_COLUMN] = calculate_slope_degrees(data).round(6)

    output.parent.mkdir(parents=True, exist_ok=True) # 상위 디렉터리 생성용
    data.to_csv(output, index=False, encoding="utf-8-sig")

    # 저장 직후 기본 무결성을 확인해 잘못된 파일을 정상 결과로 안내하지 않는다.
    saved = pd.read_csv(output, encoding="utf-8-sig")
    if len(saved) != len(data):
        raise RuntimeError(
            f"저장 전후 행 수가 다릅니다: {len(data)}행 -> {len(saved)}행"
        )
    if SLOPE_COLUMN not in saved.columns:
        raise RuntimeError(f"저장 결과에 {SLOPE_COLUMN!r} 컬럼이 없습니다.")

    return output


def main() -> None:
    output = add_slope_column()
    print(f"저장 완료: {output}")


if __name__ == "__main__":
    main()
