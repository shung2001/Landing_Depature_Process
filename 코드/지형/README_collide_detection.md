# 헬리포트 접근면 충돌 분석

프로그램: `코드/지형/collide_detection.py`

모든 출력과 로그는 `자료/결과물`, 임시 파일은 그 아래 `_temp`에 저장한다.
입력 GIS 파일은 읽기 전용으로 열고, 실행 전후 입력 폴더의 SHA-256을 비교한다.
정상 종료하면 이번 실행의 임시 파일을 정리한다. `terrain_clip.py`는 호출하지 않는다.

## 현재 입력으로 확인한 사항

- 수정 건물 SHP: 활성 피처 27,261개. DBF의 물리적 레코드 수와 다를 수 있다.
- 헬리포트 SHP: 36개. 서로 동일한 도형 18쌍을 입력 그대로 유지한다.
- 수정 건물 SHP에서 이 헬리포트의 host 건물이 제외되어 있어 기본 매칭 결과는 36개 모두 미매칭이다.
- 기본 높이 필드는 A16, 층수 필드는 A26이다. A16이 0/결측인 건물 11,978개를 0m로 간주하지 않는다.
- DEM은 90m 해상도이며, 건물 면적의 유효 DEM 피복률이 기준에 미달하면 지표고도 미확정으로 처리한다.
- 일부 원본 한글 속성은 이미 U+FFFD(�) 문자를 포함한다. UTF-8-SIG 출력은 추가적인 인코딩 문제를 방지하며 원본에 이미 손실된 글자를 복구하지는 않는다.
- 접근면 길이, 시작 폭, 끝 폭, 경사, 시작점 거리가 제공되지 않아 초기 CONFIG는 `None`이다.

현재 생성된 결과는 입력 검증 결과이다. 각도표의 빈 판정은 **판정 불가**이며,
충돌 건물 GPKG가 비어 있는 것을 '충돌 건물 없음'으로 해석하면 안 된다.

## 실행

Python 3.11 이상을 사용한다. **`pip install gdal`은 필요하지 않다.**
일반 Python에서 Rasterio를 먼저 사용하며, 없으면 이미 설치된 GDAL Python 모듈,
그것도 없으면 QGIS/OSGeo4W의 `gdalinfo`와 `gdal_translate` 실행 파일을 자동으로 찾는다.
QGIS의 실행 파일은 Python 버전과 독립적으로 사용하므로 일반 Python에서도 실행 가능하다.

프로젝트 루트에서 PowerShell:

```powershell
python -B '코드\지형\collide_detection.py' --inspect-only
```

QGIS의 실행 파일을 직접 선택하여 확인할 수도 있다.

```powershell
python -B '코드\지형\collide_detection.py' --inspect-only --dem-backend gdal-cli
```

설치 위치가 일반적인 위치와 다르면 코드 상단 `GDAL_BIN_DIR`에 QGIS의 `bin` 경로를 지정한다.
프로그램이 패키지를 자동 설치하거나 시스템의 Python/GDAL 환경을 변경하지는 않는다.
외부 GDAL 방식에서 필요한 임시 래스터와 마스크는 `자료/결과물/_temp`에만 생성하고 정리한다.

새 Python 환경에 의존성을 설치할 때는 다음 명령을 사용한다. requirements는 바이너리 wheel만
허용하므로 C++ 소스 빌드를 시도하지 않는다. 호환 wheel이 없는 환경에서는 QGIS 환경을 사용할 수 있다.

```powershell
python -m pip install -r '코드\지형\requirements.txt'
```

GDAL 설치 로그의 `Microsoft Visual C++ 14.0 or greater is required`는 소스 빌드 중 발생한 오류이다.
GDAL Python 바인딩을 직접 빌드하려면 컴파일러뿐 아니라 GDAL 라이브러리와 개발 헤더도 필요하다
([GDAL 설치 문서](https://gdal.org/en/stable/api/python/python_bindings.html)).
이 프로그램은 위의 자동 읽기 경로로 실행할 수 있다.

실제 분석 시 코드 상단 CONFIG의 아래 값을 입력하고 `--inspect-only` 없이 실행한다.
각 값은 동일한 이름의 명령행 옵션으로도 지정할 수 있다.

| CONFIG | 명령행 옵션 | 의미 |
| --- | --- | --- |
| APPROACH_LENGTH_M | --length-m | 접근면 길이(m) |
| APPROACH_INNER_WIDTH_M | --inner-width-m | 시작 폭 전체(m) |
| APPROACH_OUTER_WIDTH_M | --outer-width-m | 끝 폭 전체(m), 시작 폭 이상 |
| APPROACH_SLOPE_PERCENT | --slope-percent | 경사(%), 5는 0.05 |
| APPROACH_START_OFFSET_M | --start-offset-m | 헬리포트 중심에서 면 시작점까지 거리(m) |
| HELIPAD_ROOF_OFFSET_M | --roof-offset-m | host 지붕보다 높은 헬리포트 추가 높이(m), 기본 0 |

좌우 각 변의 확폭률이 `d`라면 `끝 폭 = 시작 폭 + 2 × 길이 × d`로 바꿔 입력한다.
구간별 경사가 다른 면은 현재의 단일 선형 접근면과 다르므로 별도 구현이 필요하다.
면의 방향은 격자북 0°, 동쪽 90°이며 진북 보정은 하지 않는다.

host 높이 자료는 다음 중 명시적으로 선택한다.

- `--host-source buildings`: 지정 수정 건물 SHP와 공간 매칭. 기본값.
- `--host-source reference`: CONFIG의 `HOST_REFERENCE_PATH`에 지정한 원본 건물 SHP를 host 확인용으로만 추가 읽기. 장애물 후보는 계속 수정 건물 SHP이다.
- `--host-source helipad`: 헬리포트 SHP 자체가 host 건물의 도형/높이 속성을 갖고 있을 때 사용. 일반적인 헬리포트 점/작은 패드 도형에는 적용하지 않는다.

높이는 기본적으로 추정하지 않는다. `--floor-height-m`를 지정한 경우에만 A16의 0/결측값을
양의 A26 × 지정 층고로 보완한다. 층수도 0/결측이면 계속 판정 불가이다.
높이 통계 방법과 필드명 등은 코드 상단 CONFIG에서 관리한다.

## 계산과 결과 해석

각 헬리포트의 중심은 입력 도형의 평면 중심점(centroid)이다. 건물 지상 높이와 DEM 중앙값을
더하여 건물 상단 절대고도를 구하고, host 건물 상단을 헬리포트 기준 높이로 사용한다.
DEM 값의 수직 기준/단위와 건물 높이 필드의 의미는 입력 제공 기준에 맞아야 한다.

축 각도 0~179°마다 `axis`와 `axis + 180°` 두 접근면을 생성한다. 건물과 접근면의 실제
교집합 꼭짓점에서 선형 면의 최소 높이를 계산한다. XY 경계 접촉도 검사에 포함한다.
`건물 상단 - 최소 면 높이 > 0.000001m`인 경우만 충돌 이벤트이다.
맞닿는 높이는 충돌로 세지 않으며, host 건물은 장애물 후보에서 제외한다.

| 파일 | 내용 |
| --- | --- |
| collision_detail.csv | 충돌 이벤트만 저장. 헬리포트·축·방향·건물별 1행 |
| collision_building_summary.csv | 헬리포트·건물별 충돌 축 개수, 정렬한 각도와 연속 구간 |
| helipad_angle_summary.csv | 모든 헬리포트의 180개 축 각도. 현재 입력은 6,480행 |
| collision_buildings.gpkg | 충돌 건물 원본 속성과 도형, 집계 필드. 레이어 이름 collision_buildings |
| unmatched_helipads.csv | 공간 host 미매칭/모호한 매칭 사유 |
| building_height_diagnostics.csv | 모든 건물의 ID, 원본 FID, 높이 출처, DEM 피복률, 계산 높이 |
| helipad_host_matches.csv | 전체 헬리포트의 host 매칭 상태와 중심, 기준 고도 |
| angle_diagnostics.csv | 각도별 판정 상태와 알려진 충돌 개수, 높이 불명 건물 ID |
| validation_report.txt | JSON 형식의 설정·입력 점검·출력 경로·원본 해시·결과 의미 |
| analysis.log | 실행 로그 |

CSV는 모두 UTF-8-SIG이다. 건물/헬리포트 ID는 기본적으로 입력 SHP의 실제 FID에서 생성하여
빈 원본 ID나 중복 A0 때문에 피처가 합쳐지지 않게 한다. GPKG의 `source_fid`는 원본 FID이며,
원본 속성 A0/A1 등도 보존한다. GPKG 자체 FID와 `source_fid`는 다른 필드이다.

각도별 `collision_count`는 양방향에서 충돌한 고유 건물 수이다. GPKG의 `collision_count`는
상세 이벤트 행 수로, 동일 축의 양방향 충돌은 2회가 될 수 있다. 범위는 0~179°에서 병합하며
179°와 0°를 연결하지 않는다. GPKG의 `helipad_axis_ranges`는 헬리포트별 범위를 JSON으로 저장한다.

각도표의 `is_collision=True`는 확인된 충돌이 있다는 의미이고, `False`는 검사가 완료되었으며
충돌이 없다는 의미이다. host/면 설정/후보 건물 높이가 불명인 경우 판정은 빈 값이다.
확인된 충돌과 높이 불명 후보가 동시에 있으면 `is_collision=True`, `collision_count`는 빈 값,
`max_penetration_m`는 확인된 충돌 중 최대값(전체 최대값의 하한)이 된다. 사유는 진단 CSV에 기록한다.

종료 코드: 0=전체 각도 판정 완료, 2=설정/데이터 미확정 또는 검증 모드, 1=실행 오류.
정상 저장 시 같은 이름의 기존 결과 파일을 갱신한다. 검증 모드도 결과 파일을 갱신한다.

## 회귀 검증

```powershell
python -B '코드\지형\test_collide_detection.py'
```

합성 도형을 사용하여 방향/거리, 교집합의 최소 높이, 양방향 중복 집계, 결측 판정,
DEM NoData/피복률/회전/좌표변환, host 매칭, 빈 GPKG, 한글 속성과 CSV BOM을 검증한다.
