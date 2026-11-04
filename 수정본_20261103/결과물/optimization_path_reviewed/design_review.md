# 초기 조건 대조 및 설계 재검토

기존 결과는 수학적으로 허용된 RSR 후보를 골랐지만, 큰 원형 선회를 제한하지 않아
구배를 맞추기 위해 거의 한 바퀴 도는 항로가 선택되었다. 또 일부 FAF 쌍에서
Landing_Area 회피가 불가능하다는 이유로 다른 쌍까지 관통 금지를 완화했다.
두 부분을 수정했다. 이 문서는 실제 비행절차 승인 문서가 아니다.

| 요구사항 | 이번 처리 |
|---|---|
| 이번 범위는 내곽 링과 FAF | 외곽 링·Interchange·Outer→IF는 계산하지 않음 |
| IF 제외 | 새로운 IF나 중간 항법 지점을 추가하지 않음 |
| 내곽·외곽 높이 350m | 이번 내곽 시작/합류 높이 350m 고정. 불리해도 임의로 낮추지 않음 |
| FAF 고도는 alt | 기하 Z나 평균 259.3m로 대체하지 않음. 입력 수직 기준과 동일하다고 가정 |
| 링 후보 1000개 | 공급된 내곽 꼭짓점으로 확인한 반경 2000m 링에 1000개 균등 배치 |
| 내곽 속도 70kt, FAF 60kt | Approach 70→60kt, Departure 60→70kt. 원본의 80kt 속성은 사용하지 않음 |
| CSC, 내곽 시계방향 RSR | 양방향 각각 새 RSR 생성. Approach 역순을 Departure로 사용하지 않음 |
| FAF 진입·이탈 방향 | Approach는 FAF→해당 OLS 꼭짓점 방향으로 끝나고, Departure는 그 반대 방향으로 출발 |
| 링 접선 | Approach 시작과 Departure 끝이 시계방향 링 접선에 일치하는지 해석적으로 검증 |
| 최대 하강 2.36°·250ft/NM | 더 엄격한 250ft/NM 적용. 동일 원리 요청에 따라 Departure 상승에도 같은 수치 사용 |
| 뱅크 최대 30°, 최소 선회반경 | 70kt 최대 속도에서 약 229m 이상. 229/350/500/750m 후보 비교 |
| 낮은 뱅크 선호 | 길이와 함께 연성 비용에 반영 |
| FAA/NASA 기준 | 사용한 반경은 협조선회 관계식의 모델값. 기종별 승인 반경이라고 주장하지 않음 |
| Landing_Area 우선 회피 | 각 FAF 쌍에 대해 회피 가능한 경우 hard 유지. 불가능한 쌍만 관통 개수·길이·수직 높이 차 penalty |
| 한 쌍의 FAF 중 하나를 STAR에 사용 | 다른 FAF를 반드시 SID에 사용. 한 FAF를 양쪽에 중복 배정하지 않음 |
| 항로 간 간격 최대한 확보·교차 최소화 | 모든 Approach/Departure 조합의 XY 교차와 최소 간격을 함께 평가. 100m는 변경 가능한 목표이며 hard 보장값 아님 |
| 이전 수직 70m 조건 제외 | Outer→IF와 Inner Ring에 관한 조건이므로 이번 항로 쌍에 재적용하지 않음 |
| 노드 집중 억제 | 링 노드와 주변 노드의 근접성에 penalty |
| 링 내부/외부 비교 | 내부는 원호의 해석적 최대 반경을 검사. 외부 허용은 경계 제약만 해제. 외부 필수는 별도 보조 비교 |
| Approach/Departure SHP 분리 | 각 시나리오에 Approach_z.shp와 Departure_z.shp 별도 저장, 각각 Z 검증 |
| 조건을 만족하지 못한 경우 | 누락 FAF 쌍을 표시. 고도·구배·회전 방향을 임의로 완화하거나 원형 우회로 채우지 않음 |

## 이번 피드백에 따른 추가 형태 필터

총 시계방향 선회량 180° 이하, 각 원호 180° 이하, 자기교차 금지를 기본으로 적용한다.
이는 한 바퀴에 가까운 우회 해를 제거하기 위한 **모델 선택**이며 초기 사용자 수치나
FAA/NASA의 공통 규정이라고 주장하지 않는다. --max-turn-deg와 --max-arc-deg로 조정할 수 있다.
350m와 FAF alt 차이를 이 필터 안에서 해소할 길이가 없으면 그 조합은 제외한다.

## 범위와 남은 한계

FAF 이후 Landing Area→버티포트의 -1m/s² 감속과 35kt 역천이는 처음에 정한
IAF→FAF 및 반대편 FAF→링 범위 밖이므로 이번 경로에 끼워 넣지 않았다.
별도 IF나 마지막 직선 구간을 덧붙이면 요청한 CSC 형태 자체가 바뀌므로 임의 추가하지 않았다.
현재 FAF에서 방향 접선은 맞추지만, 실제 롤 전이·2.36°에서 8°로의 종단 전이·기체 성능은 모델링하지 않는다.
무풍 TAS 가정이며 KIAS/TAS 변환, 실제 운항 시각, 지형·건물 장애물 분석도 포함하지 않는다.
따라서 이 결과는 조건별 기하 후보 분석이며 실제 운항 가능한 절차의 승인이나 안전성 보증이 아니다.

불관용 기준 비교안(strict_noncrossing)은 Landing_Area 관통과 항로 교차를 모두 금지한 채
연결 가능한 FAF 쌍 수를 먼저 최대화하고, 그 수를 유지하며 최소 간격을 개선한다.
Landing_Area 페널티 허용 부분 구성(noncrossing_subset)과 전체 연결 비교안은 별도 파일로 남긴다.
서로 다른 조건의 결과를 하나의 최종 항로처럼 합치지 않는다.

개별 IAF/FAF/반경 후보를 평가한 뒤 공간적으로 분산된 후보 집합에서 조합을 최적화한다.
전역 최적성은 증명하지 않는다. 교차 0건 MILP 탐색이 시간 제한에 도달하면 이를 불가능성 증명으로 쓰지 않는다.
교차가 남으면 최종 CSV와 지도에 그대로 표시한다.

## 실제 실행 요약

```json
{
  "geometry_review": {
    "total_gradient_bank_feasible": 111466,
    "shape_feasible": 10175,
    "excessive_turn_rejected": 101291,
    "self_intersection_rejected": 0,
    "max_arc_deg": 180.0,
    "max_total_turn_deg": 180.0,
    "angle_limit_is_model_assumption_not_regulatory_standard": true
  },
  "scenarios": {
    "inside_only": {
      "landing_mode": "pairwise_strict_then_penalized",
      "eligible_candidates": 7816,
      "shortlisted_candidates": 975,
      "strict_pairs": [
        3,
        6,
        7
      ],
      "penalized_pairs": [
        1,
        2,
        4,
        5,
        8
      ],
      "strict_missing_pairs": [
        1,
        2,
        4,
        5,
        8
      ],
      "missing_pairs": [],
      "zero_crossing_search": {
        "status": "time_limit_no_zero_crossing_incumbent",
        "attempts": 35,
        "seconds": 20.295077099999617,
        "infeasibility_proven": false
      },
      "noncrossing_subset": {
        "score": 145285.79227990762,
        "selected_pairs": 4,
        "selected_routes": 8,
        "route_crossing_pairs": 0,
        "minimum_distance_m": 165.10191183523528,
        "below_target_pairs": 0,
        "global_optimality_proven": false,
        "multistarts": 80,
        "solver": {
          "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
          "selected_pairs": 4,
          "status_code": 0,
          "minimum_spacing_m": 0.0,
          "coverage_optimal_in_shortlist": true,
          "not_a_full_candidate_optimality_proof": true,
          "seconds": 1.9913004000009096,
          "crossing_constraints": 103035
        },
        "spacing_search": {
          "status": "bounded_maxmin_search",
          "retained_pairs": 4,
          "achieved_minimum_m": 165.10191183523528,
          "upper_bound_m": 168.3342753071908,
          "endpoint_upper_bound_m": 540.2929859734319,
          "tolerance_m": 2.0,
          "maxmin_optimal_to_tolerance_in_shortlist": false,
          "seconds": 19.88166689999889,
          "probes": [
            {
              "goal_m": 100.0,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 4,
              "status_code": 0,
              "minimum_spacing_m": 100.0,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.6535017999995034,
              "crossing_constraints": 118888
            },
            {
              "goal_m": 324.64414691549314,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 324.64414691549314,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.3942270999996254,
              "crossing_constraints": 145426
            },
            {
              "goal_m": 216.81972738652377,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 216.81972738652377,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.4744190000001254,
              "crossing_constraints": 134139
            },
            {
              "goal_m": 162.90751762203908,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 4,
              "status_code": 0,
              "minimum_spacing_m": 162.90751762203908,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.396774299999379,
              "crossing_constraints": 126962
            },
            {
              "goal_m": 190.96081961087953,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 190.96081961087953,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.5189159999990807,
              "crossing_constraints": 131011
            },
            {
              "goal_m": 178.0313657230574,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 178.0313657230574,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.444930800000293,
              "crossing_constraints": 129364
            },
            {
              "goal_m": 171.56663877914633,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 171.56663877914633,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.389258199998949,
              "crossing_constraints": 128400
            },
            {
              "goal_m": 168.3342753071908,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 168.3342753071908,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.437375699999393,
              "crossing_constraints": 127910
            }
          ]
        },
        "missing_pairs": [
          1,
          2,
          3,
          8
        ],
        "length_m": 13714.243547072307,
        "landing_crossings": 6,
        "exact_crossing_pairs": 0,
        "exact_min_distance_m": 165.10191183522005
      },
      "score": 21295057.257534653,
      "selected_pairs": 8,
      "selected_routes": 16,
      "route_crossing_pairs": 21,
      "minimum_distance_m": 0.0,
      "below_target_pairs": 29,
      "global_optimality_proven": false,
      "multistarts": 80,
      "star_count": 8,
      "sid_count": 8,
      "length_m": 38800.615385644945,
      "landing_crossings": 8,
      "outside_routes": 0,
      "exact_crossing_pairs": 21,
      "exact_min_distance_m": 0.0,
      "median_pair_distance_m": 763.9590982642288,
      "below_spacing_target_pairs": 29,
      "roundtrip_and_analytic_boundary_validated": true,
      "strict_noncrossing": {
        "score": 16.28622196519872,
        "selected_pairs": 3,
        "selected_routes": 6,
        "route_crossing_pairs": 0,
        "minimum_distance_m": 162.20582113343752,
        "below_target_pairs": 0,
        "global_optimality_proven": false,
        "multistarts": 80,
        "solver": {
          "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
          "selected_pairs": 3,
          "status_code": 0,
          "minimum_spacing_m": 0.0,
          "coverage_optimal_in_shortlist": true,
          "not_a_full_candidate_optimality_proof": true,
          "seconds": 0.11733290000120178,
          "crossing_constraints": 11599
        },
        "spacing_search": {
          "status": "bounded_maxmin_search",
          "retained_pairs": 3,
          "achieved_minimum_m": 162.20582113343752,
          "upper_bound_m": 163.42857788506208,
          "endpoint_upper_bound_m": 512.9629477284523,
          "tolerance_m": 2.0,
          "maxmin_optimal_to_tolerance_in_shortlist": true,
          "seconds": 0.9336428999995405,
          "probes": [
            {
              "goal_m": 100.0,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 3,
              "status_code": 0,
              "minimum_spacing_m": 100.0,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.08976080000138609,
              "crossing_constraints": 13303
            },
            {
              "goal_m": 320.2166588244754,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 320.2166588244754,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.10357350000049337,
              "crossing_constraints": 18908
            },
            {
              "goal_m": 223.84351437248702,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 223.84351437248702,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09500300000036077,
              "crossing_constraints": 17374
            },
            {
              "goal_m": 175.65694214649284,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 175.65694214649284,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.08440499999960593,
              "crossing_constraints": 15861
            },
            {
              "goal_m": 151.56365603349573,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 3,
              "status_code": 0,
              "minimum_spacing_m": 151.56365603349573,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09052109999902314,
              "crossing_constraints": 14812
            },
            {
              "goal_m": 164.65133463668664,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 164.65133463668664,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.081419599999208,
              "crossing_constraints": 15384
            },
            {
              "goal_m": 159.14853088178353,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 3,
              "status_code": 0,
              "minimum_spacing_m": 159.14853088178353,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.10690040000008594,
              "crossing_constraints": 15027
            },
            {
              "goal_m": 162.1349924818548,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 3,
              "status_code": 0,
              "minimum_spacing_m": 162.1349924818548,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09547339999880933,
              "crossing_constraints": 15183
            },
            {
              "goal_m": 163.42857788506208,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 163.42857788506208,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.08104629999979807,
              "crossing_constraints": 15254
            }
          ]
        },
        "missing_pairs": [
          1,
          2,
          4,
          5,
          8
        ],
        "length_m": 11867.920552874915,
        "landing_crossings": 0,
        "exact_crossing_pairs": 0,
        "exact_min_distance_m": 162.20582113343752
      }
    },
    "outside_allowed": {
      "landing_mode": "pairwise_strict_then_penalized",
      "eligible_candidates": 9267,
      "shortlisted_candidates": 1011,
      "strict_pairs": [
        3,
        6,
        7,
        8
      ],
      "penalized_pairs": [
        1,
        2,
        4,
        5
      ],
      "strict_missing_pairs": [
        1,
        2,
        4,
        5
      ],
      "missing_pairs": [],
      "zero_crossing_search": {
        "status": "infeasible_in_shortlist",
        "attempts": 12,
        "seconds": 2.4883458000003884,
        "not_a_full_candidate_infeasibility_proof": true
      },
      "noncrossing_subset": {
        "score": 145285.79227990762,
        "selected_pairs": 4,
        "selected_routes": 8,
        "route_crossing_pairs": 0,
        "minimum_distance_m": 165.10191183523528,
        "below_target_pairs": 0,
        "global_optimality_proven": false,
        "multistarts": 80,
        "solver": {
          "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
          "selected_pairs": 4,
          "status_code": 0,
          "minimum_spacing_m": 0.0,
          "coverage_optimal_in_shortlist": true,
          "not_a_full_candidate_optimality_proof": true,
          "seconds": 1.7623924999988958,
          "crossing_constraints": 97772
        },
        "spacing_search": {
          "status": "bounded_maxmin_search",
          "retained_pairs": 4,
          "achieved_minimum_m": 165.10191183523528,
          "upper_bound_m": 165.3382045919854,
          "endpoint_upper_bound_m": 540.2929859734319,
          "tolerance_m": 2.0,
          "maxmin_optimal_to_tolerance_in_shortlist": true,
          "seconds": 19.165500499999325,
          "probes": [
            {
              "goal_m": 100.0,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 4,
              "status_code": 0,
              "minimum_spacing_m": 100.0,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.2383454000009806,
              "crossing_constraints": 116475
            },
            {
              "goal_m": 321.5851925422938,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 321.5851925422938,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.1665502000014385,
              "crossing_constraints": 147395
            },
            {
              "goal_m": 212.23129582672473,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 212.23129582672473,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.0540658999998413,
              "crossing_constraints": 132752
            },
            {
              "goal_m": 157.5543474689402,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 4,
              "status_code": 0,
              "minimum_spacing_m": 157.5543474689402,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.0603716999994504,
              "crossing_constraints": 124374
            },
            {
              "goal_m": 187.22164716819708,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 187.22164716819708,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.095273900000393,
              "crossing_constraints": 128927
            },
            {
              "goal_m": 174.71682283893327,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 174.71682283893327,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.08777259999988,
              "crossing_constraints": 127010
            },
            {
              "goal_m": 168.46441067430135,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 168.46441067430135,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.0755069000006188,
              "crossing_constraints": 126162
            },
            {
              "goal_m": 165.3382045919854,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 165.3382045919854,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.0893450999992638,
              "crossing_constraints": 125690
            },
            {
              "goal_m": 163.77510155082743,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 4,
              "status_code": 0,
              "minimum_spacing_m": 163.77510155082743,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 2.105206799998996,
              "crossing_constraints": 125549
            }
          ]
        },
        "missing_pairs": [
          1,
          2,
          3,
          8
        ],
        "length_m": 13714.243547072307,
        "landing_crossings": 6,
        "exact_crossing_pairs": 0,
        "exact_min_distance_m": 165.10191183522005
      },
      "score": 16235464.682588423,
      "selected_pairs": 8,
      "selected_routes": 16,
      "route_crossing_pairs": 16,
      "minimum_distance_m": 0.0,
      "below_target_pairs": 27,
      "global_optimality_proven": false,
      "multistarts": 80,
      "star_count": 8,
      "sid_count": 8,
      "length_m": 38682.79941261235,
      "landing_crossings": 6,
      "outside_routes": 5,
      "exact_crossing_pairs": 16,
      "exact_min_distance_m": 0.0,
      "median_pair_distance_m": 748.3995369428945,
      "below_spacing_target_pairs": 27,
      "roundtrip_and_analytic_boundary_validated": true,
      "strict_noncrossing": {
        "score": 16.674839897968965,
        "selected_pairs": 3,
        "selected_routes": 6,
        "route_crossing_pairs": 0,
        "minimum_distance_m": 198.18789584368287,
        "below_target_pairs": 0,
        "global_optimality_proven": false,
        "multistarts": 80,
        "solver": {
          "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
          "selected_pairs": 3,
          "status_code": 0,
          "minimum_spacing_m": 0.0,
          "coverage_optimal_in_shortlist": true,
          "not_a_full_candidate_optimality_proof": true,
          "seconds": 0.18780050000168558,
          "crossing_constraints": 19183
        },
        "spacing_search": {
          "status": "bounded_maxmin_search",
          "retained_pairs": 3,
          "achieved_minimum_m": 197.80266011301563,
          "upper_bound_m": 198.90970707494714,
          "endpoint_upper_bound_m": 770.1682957590357,
          "tolerance_m": 2.0,
          "maxmin_optimal_to_tolerance_in_shortlist": true,
          "seconds": 1.7015360000004875,
          "probes": [
            {
              "goal_m": 453.223790470398,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 453.223790470398,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.18664650000027905,
              "crossing_constraints": 35255
            },
            {
              "goal_m": 294.7515378260792,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 294.7515378260792,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.1716377000011562,
              "crossing_constraints": 31820
            },
            {
              "goal_m": 215.51541150391978,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 215.51541150391978,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.16452170000047772,
              "crossing_constraints": 28970
            },
            {
              "goal_m": 175.89734834284008,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 3,
              "status_code": 0,
              "minimum_spacing_m": 175.89734834284008,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.18200100000103703,
              "crossing_constraints": 27115
            },
            {
              "goal_m": 195.9346425934936,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 3,
              "status_code": 0,
              "minimum_spacing_m": 195.9346425934936,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.19028240000079677,
              "crossing_constraints": 28103
            },
            {
              "goal_m": 206.6590358084677,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 206.6590358084677,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.17064170000048762,
              "crossing_constraints": 28628
            },
            {
              "goal_m": 202.23084796074167,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 202.23084796074167,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.1717090999991342,
              "crossing_constraints": 28417
            },
            {
              "goal_m": 200.01675403687864,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 200.01675403687864,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.16526789999988978,
              "crossing_constraints": 28282
            },
            {
              "goal_m": 198.90970707494714,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 198.90970707494714,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.164861799999926,
              "crossing_constraints": 28254
            }
          ]
        },
        "missing_pairs": [
          1,
          2,
          4,
          5,
          8
        ],
        "length_m": 11552.958285678707,
        "landing_crossings": 0,
        "exact_crossing_pairs": 0,
        "exact_min_distance_m": 198.1878958436852
      }
    },
    "outside_detour": {
      "landing_mode": "pairwise_strict_then_penalized",
      "eligible_candidates": 1608,
      "shortlisted_candidates": 578,
      "strict_pairs": [
        6,
        7,
        8
      ],
      "penalized_pairs": [
        1,
        2,
        3,
        4,
        5
      ],
      "strict_missing_pairs": [
        1,
        2,
        3,
        4,
        5
      ],
      "missing_pairs": [
        1,
        2,
        3,
        4,
        5
      ],
      "zero_crossing_search": {
        "status": "infeasible_in_shortlist",
        "attempts": 12,
        "seconds": 0.2464030000010098,
        "not_a_full_candidate_infeasibility_proof": true
      },
      "noncrossing_subset": {
        "score": 5.91773285599195,
        "selected_pairs": 1,
        "selected_routes": 2,
        "route_crossing_pairs": 0,
        "minimum_distance_m": 1236.0679774998018,
        "below_target_pairs": 0,
        "global_optimality_proven": false,
        "multistarts": 80,
        "solver": {
          "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
          "selected_pairs": 1,
          "status_code": 0,
          "minimum_spacing_m": 0.0,
          "coverage_optimal_in_shortlist": true,
          "not_a_full_candidate_optimality_proof": true,
          "seconds": 0.06781140000020969,
          "crossing_constraints": 11600
        },
        "spacing_search": {
          "status": "bounded_maxmin_search",
          "retained_pairs": 1,
          "achieved_minimum_m": 1236.0679774998018,
          "upper_bound_m": 1238.0284730348721,
          "endpoint_upper_bound_m": 1992.3914666437208,
          "tolerance_m": 2.0,
          "maxmin_optimal_to_tolerance_in_shortlist": true,
          "seconds": 0.8473823000003904,
          "probes": [
            {
              "goal_m": 1529.9997119626119,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1529.9997119626119,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09603729999980715,
              "crossing_constraints": 32628
            },
            {
              "goal_m": 1298.8038346220574,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1298.8038346220574,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09442529999978433,
              "crossing_constraints": 32068
            },
            {
              "goal_m": 1183.2058959517801,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 1,
              "status_code": 0,
              "minimum_spacing_m": 1183.2058959517801,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.1056365000004007,
              "crossing_constraints": 31970
            },
            {
              "goal_m": 1267.4359060609295,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1267.4359060609295,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.0948645000007673,
              "crossing_constraints": 31988
            },
            {
              "goal_m": 1251.7519417803655,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1251.7519417803655,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09645790000104171,
              "crossing_constraints": 31988
            },
            {
              "goal_m": 1243.9099596400838,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1243.9099596400838,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09476200000062818,
              "crossing_constraints": 31988
            },
            {
              "goal_m": 1239.9889685699427,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1239.9889685699427,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.10431450000032783,
              "crossing_constraints": 31988
            },
            {
              "goal_m": 1238.0284730348721,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1238.0284730348721,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.0958864000003814,
              "crossing_constraints": 31988
            }
          ]
        },
        "missing_pairs": [
          1,
          2,
          3,
          4,
          5,
          7,
          8
        ],
        "length_m": 5251.107800504865,
        "landing_crossings": 0,
        "exact_crossing_pairs": 0,
        "exact_min_distance_m": 1236.0679774998018
      },
      "score": 5006021.206894411,
      "selected_pairs": 3,
      "selected_routes": 6,
      "route_crossing_pairs": 5,
      "minimum_distance_m": 0.0,
      "below_target_pairs": 5,
      "global_optimality_proven": false,
      "multistarts": 80,
      "star_count": 3,
      "sid_count": 3,
      "length_m": 14857.945233109227,
      "landing_crossings": 0,
      "outside_routes": 6,
      "exact_crossing_pairs": 5,
      "exact_min_distance_m": 0.0,
      "median_pair_distance_m": 625.2449353832104,
      "below_spacing_target_pairs": 5,
      "roundtrip_and_analytic_boundary_validated": true,
      "strict_noncrossing": {
        "score": 5.91773285599195,
        "selected_pairs": 1,
        "selected_routes": 2,
        "route_crossing_pairs": 0,
        "minimum_distance_m": 1236.0679774998018,
        "below_target_pairs": 0,
        "global_optimality_proven": false,
        "multistarts": 80,
        "solver": {
          "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
          "selected_pairs": 1,
          "status_code": 0,
          "minimum_spacing_m": 0.0,
          "coverage_optimal_in_shortlist": true,
          "not_a_full_candidate_optimality_proof": true,
          "seconds": 0.06697859999985667,
          "crossing_constraints": 11600
        },
        "spacing_search": {
          "status": "bounded_maxmin_search",
          "retained_pairs": 1,
          "achieved_minimum_m": 1236.0679774998018,
          "upper_bound_m": 1238.0284730348721,
          "endpoint_upper_bound_m": 1992.3914666437208,
          "tolerance_m": 2.0,
          "maxmin_optimal_to_tolerance_in_shortlist": true,
          "seconds": 0.8469746000009764,
          "probes": [
            {
              "goal_m": 1529.9997119626119,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1529.9997119626119,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.0966472000000067,
              "crossing_constraints": 32628
            },
            {
              "goal_m": 1298.8038346220574,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1298.8038346220574,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09584220000033383,
              "crossing_constraints": 32068
            },
            {
              "goal_m": 1183.2058959517801,
              "status": "Optimization terminated successfully. (HiGHS Status 7: Optimal)",
              "selected_pairs": 1,
              "status_code": 0,
              "minimum_spacing_m": 1183.2058959517801,
              "coverage_optimal_in_shortlist": true,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.10793659999944794,
              "crossing_constraints": 31970
            },
            {
              "goal_m": 1267.4359060609295,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1267.4359060609295,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09523789999911969,
              "crossing_constraints": 31988
            },
            {
              "goal_m": 1251.7519417803655,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1251.7519417803655,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.0968261000016355,
              "crossing_constraints": 31988
            },
            {
              "goal_m": 1243.9099596400838,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1243.9099596400838,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.0952028000010614,
              "crossing_constraints": 31988
            },
            {
              "goal_m": 1239.9889685699427,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1239.9889685699427,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09697940000114613,
              "crossing_constraints": 31988
            },
            {
              "goal_m": 1238.0284730348721,
              "status": "The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)",
              "selected_pairs": 0,
              "status_code": 2,
              "minimum_spacing_m": 1238.0284730348721,
              "coverage_optimal_in_shortlist": false,
              "not_a_full_candidate_optimality_proof": true,
              "seconds": 0.09591329999966547,
              "crossing_constraints": 31988
            }
          ]
        },
        "missing_pairs": [
          1,
          2,
          3,
          4,
          5,
          7,
          8
        ],
        "length_m": 5251.107800504865,
        "landing_crossings": 0,
        "exact_crossing_pairs": 0,
        "exact_min_distance_m": 1236.0679774998018
      }
    }
  }
}
```
