# Kepler 이동 흐름 사용법

## 먼저 볼 파일

1. **kepler_ring_flow.csv**: 계속 수요를 생성해 채운 뒤 실제로 순환하는 링 전용 실험. 링별 마지막 15분 최소 동시 이동 대수 합을 최대로 하는 시험 조건을 선정했다.
2. **kepler_best_throughput_flow.csv**: 접근선→버티포트→출항 전체 OD. 추가 대기 바퀴가 없는 최단시간 경로 실험 중 착륙+이륙 처리량이 가장 높은 조건이다.
3. **kepler_OD_flow.csv**: 추가 순환을 포함한 전체 OD 실험 중 동시 순환 지속을 우선한 조건이다.

- 처리량 우선 조건: OD_laps0_R15_O14_I14, seed=20261020, 입력 수요 15대/h.
- 착륙+이륙 25.2회/h, 왕복 완료 출항 10.8대/h, 안정 판정 False.
- 시험한 조건에서 목적별로 선정했다. 전역 최적해나 안정성이 입증된 최대 수용량을 뜻하지 않는다.

## CSV로 재생

1. Kepler에서 Add Data로 원하는 CSV를 불러온다.
2. Point 레이어의 Latitude=latitude, Longitude=longitude를 지정한다.
3. 색상은 aircraft_id 또는 ring, 툴팁은 aircraft_id·phase·time_s·altitude_m를 지정한다.
4. Filters에서 timestamp를 선택한다. 시간창을 링 실험은 5~10초, 전체 OD는 10~20초로 좁히고 재생한다.
5. kepler_ring_paths.csv를 함께 추가하면 geometry 열의 GeoJSON 선으로 원본 링을 표시할 수 있다.

CSV는 시간별 포인트이다. 과거 점을 모두 켜두면 한 기체가 여러 점으로 보이므로 현재 동시 대수를 보려면 좁은 시간창을 사용한다. 누적 표시 대수를 동시 체공 대수로 해석하지 않는다.

## 움직이는 궤적과 꼬리 표시

동일한 이름의 .geojson 파일도 제공한다. Kepler Trip 레이어용 [경도, 위도, 고도, Unix 초] 좌표이다. .geojson을 불러와 Trip 레이어로 표시하고 Trail Length를 약 10~30초로 두면 기체별 이동 방향을 볼 수 있다. 도착과 출발은 다른 세그먼트로 분리했다.

## 시각과 좌표

- 좌표: EPSG:4326, longitude/latitude는 도, altitude_m는 원본 SHP 고도(m).
- timestamp: UTC ISO 8601. timestamp_ms: Unix 밀리초.
- 시뮬레이션 0초=2026-10-20 00:00:00 UTC로 정한 가상 시각이며 실측 운항 기록이 아니다.
- 링 전용은 5초 표본과 진입 시각, 전체 OD는 10초 표본이다.

공식 문서: [CSV 형식](https://docs.kepler.gl/docs/user-guides/b-kepler-gl-workflow/a-add-data-to-the-map), [시간 재생](https://docs.kepler.gl/docs/user-guides/h-playback), [Trip 레이어](https://docs.kepler.gl/docs/user-guides/c-types-of-layers/k-trip).

재생성: `python -B -X utf8 Network/코드/export_kepler_flows.py`
