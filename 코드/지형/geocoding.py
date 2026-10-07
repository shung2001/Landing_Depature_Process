"""카카오 로컬 API로 주소를 찾고 EPSG:5186 좌표 열을 추가한다.

카카오 API의 WGS84 경위도를 pyproj로 변환하고, 입력 CSV의 ``주 소`` 열
바로 뒤에 ``x_5186``(Easting), ``y_5186``(Northing) 열을 추가한다.
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

try:
    from pyproj import Transformer
except ImportError:
    Transformer = None  # type: ignore[assignment]


# ============================ 사용자 설정 ============================
# 카카오 디벨로퍼스의 'REST API 키'를 따옴표 안에 입력하세요.
# 실제 키가 들어간 파일은 외부에 공유하거나 Git에 커밋하지 마세요.
KAKAO_REST_API_KEY = "f1161b2a855cf332983a74de75da4167"

INPUT_CSV = Path("기초자료/2020년 헬기 이착륙장 현황_착륙가능.csv")
OUTPUT_CSV = Path("기초자료/2020년 헬기 이착륙장 현황_착륙가능_5186.csv")
REQUEST_TIMEOUT = 10.0
RETRIES = 3
REQUEST_INTERVAL = 0.05
COORDINATE_PRECISION = 3
# ===================================================================


ADDRESS_URL = "https://dapi.kakao.com/v2/local/search/address.json"
KEYWORD_URL = "https://dapi.kakao.com/v2/local/search/keyword.json"
ADDRESS_COLUMN = "주 소"
BUILDING_COLUMN = "건물명"
X_COLUMN = "x_5186"
Y_COLUMN = "y_5186"


class KakaoAPIError(RuntimeError):
    """카카오 API 요청 자체가 실패한 경우."""


@dataclass(frozen=True)
class GeocodeResult:
    longitude: float
    latitude: float


def normalize_text(value: str | None) -> str:
    """CSV 셀의 앞뒤 및 연속 공백을 정리한다."""

    return re.sub(r"\s+", " ", value or "").strip()


class KakaoGeocoder:
    """주소 검색 후 필요할 때 건물명 기반 장소 검색으로 보완한다."""

    def __init__(
        self,
        api_key: str,
        *,
        timeout: float = 10.0,
        retries: int = 3,
        request_interval: float = 0.05,
        opener: Callable[..., Any] = urlopen,
    ) -> None:
        self.api_key = api_key
        self.timeout = timeout
        self.retries = retries
        self.request_interval = request_interval
        self.opener = opener
        self._cache: dict[tuple[str, str], GeocodeResult | None] = {}

    def _request(self, url: str, params: dict[str, str | int]) -> dict[str, Any]:
        request = Request(
            f"{url}?{urlencode(params)}",
            headers={
                "Authorization": f"KakaoAK {self.api_key}",
                "User-Agent": "heliport-epsg5186-geocoder/1.0",
            },
        )

        for attempt in range(self.retries + 1):
            try:
                with self.opener(request, timeout=self.timeout) as response:
                    payload = response.read().decode("utf-8")
                if self.request_interval > 0:
                    time.sleep(self.request_interval)
                return json.loads(payload)
            except HTTPError as exc:
                # 인증 오류 등은 재시도해도 해결되지 않는다.
                if exc.code < 500 and exc.code != 429:
                    detail = exc.read().decode("utf-8", errors="replace")
                    raise KakaoAPIError(
                        f"카카오 API HTTP {exc.code}: {detail}"
                    ) from exc
                last_error: Exception = exc
            except (URLError, TimeoutError, json.JSONDecodeError) as exc:
                last_error = exc

            if attempt < self.retries:
                time.sleep(min(2**attempt, 8))

        raise KakaoAPIError(
            f"카카오 API 요청이 {self.retries + 1}회 실패했습니다: {last_error}"
        )

    def _address_search(self, query: str) -> GeocodeResult | None:
        data = self._request(
            ADDRESS_URL,
            {"query": query, "analyze_type": "similar", "size": 1},
        )
        documents = data.get("documents") or []
        if not documents:
            return None
        return GeocodeResult(
            longitude=float(documents[0]["x"]),
            latitude=float(documents[0]["y"]),
        )

    def _keyword_search(self, query: str) -> GeocodeResult | None:
        data = self._request(KEYWORD_URL, {"query": query, "size": 1})
        documents = data.get("documents") or []
        if not documents:
            return None
        return GeocodeResult(
            longitude=float(documents[0]["x"]),
            latitude=float(documents[0]["y"]),
        )

    def geocode(self, address: str, building_name: str) -> GeocodeResult | None:
        """주소를 우선 검색하고, 실패하면 주소와 건물명으로 장소를 검색한다."""

        address = normalize_text(address)
        building_name = normalize_text(building_name)
        cache_key = (address, building_name)
        if cache_key in self._cache:
            return self._cache[cache_key]

        result = self._address_search(address) if address else None
        if result is None and building_name:
            # 주소를 함께 넣어 동명 건물의 다른 지역 결과가 선택될 가능성을 낮춘다.
            query = " ".join(part for part in (address, building_name) if part)
            result = self._keyword_search(query)

        self._cache[cache_key] = result
        return result


def insert_coordinate_columns(fieldnames: Iterable[str]) -> list[str]:
    """주소 열 바로 뒤에 좌표 열을 배치하고 기존 좌표 열의 중복은 제거한다."""

    columns = [name for name in fieldnames if name not in {X_COLUMN, Y_COLUMN}]
    if ADDRESS_COLUMN not in columns:
        raise ValueError(f"입력 CSV에 '{ADDRESS_COLUMN}' 열이 없습니다.")

    address_index = columns.index(ADDRESS_COLUMN) + 1
    columns[address_index:address_index] = [X_COLUMN, Y_COLUMN]
    return columns


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    # utf-8-sig는 BOM이 있거나 없는 UTF-8 CSV를 모두 읽는다.
    with path.open("r", encoding="utf-8-sig", newline="") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames is None:
            raise ValueError("입력 CSV에 헤더가 없습니다.")
        return reader.fieldnames, list(reader)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Excel에서 한글이 깨지지 않도록 UTF-8 BOM으로 저장한다.
    with path.open("w", encoding="utf-8-sig", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def geocode_csv(
    input_path: Path,
    output_path: Path,
    geocoder: KakaoGeocoder,
    *,
    precision: int = 3,
) -> tuple[int, int, int]:
    if Transformer is None:
        raise RuntimeError("pyproj가 없습니다. 'pip install pyproj'로 설치해 주세요.")

    original_columns, rows = read_csv(input_path)
    output_columns = insert_coordinate_columns(original_columns)
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:5186", always_xy=True)

    success = 0
    not_found = 0
    skipped = 0

    for index, row in enumerate(rows, start=1):
        address = normalize_text(row.get(ADDRESS_COLUMN))
        building = normalize_text(row.get(BUILDING_COLUMN))
        if not address and not building:
            row[X_COLUMN] = ""
            row[Y_COLUMN] = ""
            skipped += 1
            continue

        result = geocoder.geocode(address, building)
        if result is None:
            row[X_COLUMN] = ""
            row[Y_COLUMN] = ""
            not_found += 1
            print(
                f"[{index}/{len(rows)}] 검색 결과 없음: {address or '-'} / {building or '-'}",
                file=sys.stderr,
            )
            continue

        x_5186, y_5186 = transformer.transform(result.longitude, result.latitude)
        row[X_COLUMN] = f"{x_5186:.{precision}f}"
        row[Y_COLUMN] = f"{y_5186:.{precision}f}"
        success += 1

        if index % 20 == 0 or index == len(rows):
            print(f"진행: {index}/{len(rows)}행 (성공 {success}, 미검색 {not_found})")

    write_csv(output_path, output_columns, rows)
    return success, not_found, skipped


def main() -> int:
    api_key = KAKAO_REST_API_KEY.strip()
    if not api_key:
        print(
            "오류: 파일 상단의 KAKAO_REST_API_KEY에 카카오 REST API 키를 입력하세요.",
            file=sys.stderr,
        )
        return 2
    if RETRIES < 0:
        print("오류: RETRIES는 0 이상이어야 합니다.", file=sys.stderr)
        return 2
    if COORDINATE_PRECISION < 0:
        print("오류: COORDINATE_PRECISION은 0 이상이어야 합니다.", file=sys.stderr)
        return 2

    geocoder = KakaoGeocoder(
        api_key,
        timeout=REQUEST_TIMEOUT,
        retries=RETRIES,
        request_interval=REQUEST_INTERVAL,
    )

    try:
        success, not_found, skipped = geocode_csv(
            INPUT_CSV,
            OUTPUT_CSV,
            geocoder,
            precision=COORDINATE_PRECISION,
        )
    except (OSError, ValueError, RuntimeError, KakaoAPIError) as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 1

    print(
        f"완료: {OUTPUT_CSV} "
        f"(좌표 {success}행, 미검색 {not_found}행, 빈 행 {skipped}행)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
