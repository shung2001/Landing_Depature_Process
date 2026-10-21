# Network 용량 분석

`코드/해야할 과제.md`의 SHP 연결 관계와 운항 규칙을 적용하는 3차원 방향성 네트워크 실험이다. 원본은 프로젝트의 `자료/최종정리`에서 읽고 생성물은 `Network/결과물`에 저장한다.

## 추가된 연속 생성·순환 실험과 Kepler

**기존 안정 왕복 실험의 최대 관측 체공 대수는 링 자체의 최대 순환 용량을 뜻하지 않는다.**
`결과물/순환확장/순환분석보고서.md`에서 두 지표를 분리한 추가 실험을 확인한다.

```powershell
python -B -X utf8 Network/코드/validate_circulation.py
python -B -X utf8 Network/코드/circulation_capacity.py --workers 4
python -B -X utf8 Network/코드/export_kepler_flows.py
```

- 원본 SHP 위 사전 배치 순환: 균등·무작위 간격, 분리거리 4단계, 상한 초과 배치까지 900초씩 검증한다.
- 연속 생성 순환: 링별 120/360/720대/h를 1시간 공급한다. 단일·다중 진입점, 14/28/기하학 상한, 두 시드를 비교한다. 진입한 기체는 끝까지 이동한다. 착륙이나 링 간 환승은 포함하지 않는 독립 링 시험이다.
- 전체 OD 과부하: 기존 접근·착륙·출항 정책과 추가 순환 0/2바퀴를 비교하고, 점유 대수와 실제 동시 이동 대수를 따로 기록한다.
- Kepler: `kepler_ring_flow.csv`는 최대 지속 순환, `kepler_best_throughput_flow.csv`는 추가 바퀴 없는 최단시간 경로 중 처리량 최대 조건, `kepler_OD_flow.csv`는 전체 OD 중 동시 이동 우선 조건이다. 각 파일에 시간별 위도·경도·고도·기체 ID가 들어 있다.
- 동일 이름의 `.geojson`은 시간 좌표를 포함한 Trip 레이어용이다. 자세한 설정은 `결과물/순환확장/Kepler_사용법.md`를 참고한다.

## 실행

접근·착륙·지상 체류·이륙·출항을 포함한 10·20·30·40·50·60대 실험은 `python -B -X utf8 Network/코드/full_operation_fleet.py`로 실행한다. 결과와 Kepler 궤적은 [전체 운항 비교](결과물/고정대수_20_60/실험과_학습자료.md)에 있다. 총 투입 대수이며 동시 체공 대수는 별도로 기록한다.

프로젝트 최상위 폴더에서 실행한다.

```powershell
python -m pip install -r Network/코드/requirements.txt
python -B -X utf8 Network/코드/validate_capacity.py
python -B -X utf8 Network/코드/network_capacity.py --workers 6
python -B -X utf8 Network/코드/report_capacity.py
python -B -X utf8 Network/코드/sensitivity_capacity.py
```

짧은 동작 확인은 별도 폴더에 저장한다. 전체 분석 결과와 혼용하지 않는다.

```powershell
python -B -X utf8 Network/코드/network_capacity.py --quick --output Network/결과물/검증
python -B -X utf8 Network/코드/report_capacity.py --output Network/결과물/검증
```

## 적용 내용

- 외곽은 반시계, 내곽은 시계 방향이다. XYZ가 일치하는 끝점으로만 연결하며 다른 고도의 평면 교차는 연결하지 않는다.
- 이동시간을 가중치로 Dijkstra를 적용한다. SHP의 Speed 계열 필드를 우선하고, 누락 시 `capacity_config.json`에 명시한 kt 속도를 사용한다. 현재 Short_Final 가정은 30kt이다.
- 출입점은 `Name`, 버티포트는 `fid`로 식별한다. 도착 후 20초 호버링과 지상 체류를 거쳐 같은 기체가 출발한다. 버티포트 용량은 각각 1이다.
- 강원 10%, 나머지 3개 출발지 각각 30%를 10대 블록마다 섞는다. 마지막 불완전 블록에서는 비율이 달라질 수 있다. 목적지는 fid 1과의 거리에 따른 확률로 생성한다.
- 각 기체의 50m 버퍼는 중심 간 최소 100m로 검사한다. 링 위에서는 설정한 100/250/500/1,000m 분리거리도 매 시간 단계의 연속 상대운동으로 검사한다.
- 교차구간 선진입 예약을 유지한다. 새 경합에는 착륙 최종접근과 링 진입에 우선권을 준다. 최종접근용 교차구간을 예약하기 전 착륙장을 함께 확보하여 착륙장과 해당 교차구간을 서로 기다리는 상황을 방지한다.
- 각 링의 기본 용량은 14대다. 설정 분리거리가 허용하는 기하학적 상한이 더 작으면 상한을 적용한다. 한 번에 한쪽 링을 1대씩 늘리는 난수 과정을 저장하고 지정한 체크포인트와 상한에서 비교한다.
- `capacity_scenarios.csv`에는 평균 원호 간격과 균등 배치 시 현 길이를 함께 기록한다. `selected_spacing.csv`는 선정 조건의 실제 최근접 간격을 10초 단위로 기록한다.

## 결과 확인

우선 `결과물/분석보고서.md`와 `결과물/summary.json`을 확인한다.

- `selection_status=stable_tested`: 시험한 조건 중 모든 반복이 안정 기준을 만족하는 조건을 선정했다.
- `selection_status=diagnostic_only_no_stable_condition`: 안정 조건이 없다. 선정 조건과 체공 최댓값은 혼잡 진단용이며 최대 운영 용량으로 쓰지 않는다.
- `max_tested_stable_arrivals_per_h`: 시험한 수요 중 안정 기준을 통과한 가장 큰 실제 생성 수요다. 요청 수요는 정수 대수로 절삭될 수 있다.
- `max_stable_observed_airborne`: 안정 조건에서 관측한 최대 동시 체공 대수다. 링의 설정 용량 합이나 전역 최대 수용량과는 다른 지표다.
- `OD.csv`: 조건·시드·기체·방향별 출발지/목적지, 계획 시각, 실제 운항 시각, 완료/진행/대기 상태다. 도착과 출발이 각각 한 행이다.
- `capacity_expansion_history.csv`: 재현 가능한 1대 단위 증가 순서다. 모든 중간 단계를 시뮬레이션한 것은 아니다.
- `capacity_results.csv`: 반복별 상세 결과이며 `mature_aircraft`와 `completed_mature_aircraft`로 안정 판정 표본 수를 확인할 수 있다.

기본 수요 격자에는 낮은 수요 5·10대/h도 포함한다. 혼잡한 수요만 시험해 안정 구간을 놓치는 것을 줄이기 위한 설정이다. 90분·두 시드에서 얻은 값이므로 최대 용량을 정밀하게 추정하려면 안정 경계 주변 수요를 세분하고 운항시간과 반복 수를 늘려야 한다.

## 남아 있는 모델 한계

Holding_point의 Capacity는 유한 대기열에 적용하지만 실제 racetrack 위치와 진입·이탈 궤적은 계산하지 않는다. 따라서 Holding 기체의 실제 공간 충돌 검증은 포함하지 않는다. 현재 회피 제어는 공중 정지를 허용하고 가감속·선회반경·상승하강률을 반영하지 않는다. 교차 예약은 보수적이며 다른 구간의 혼잡 또는 순환 대기가 처리량을 제한할 수 있다. 이번 착륙장 예약 순서 수정이 모든 형태의 교착을 없앤다는 의미는 아니다.

원본 경로에 없는 속도, 지상 체류시간, 목적지 확률의 거리 척도는 가정이다. 보고서의 처리량과 동시 체공 대수는 이 가정 아래의 유한 실험 결과이며, 실제 항공기의 인증 용량으로 해석하지 않는다.
