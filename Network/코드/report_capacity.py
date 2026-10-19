"""Produce the Korean analysis report and audit/visual artifacts from completed runs."""
import argparse
import sys
sys.dont_write_bytecode = True
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point
from shapely import get_coordinates
from network_capacity import Network, ROOT, OUTPUT, configure_plot, csv, plots, capacity_expansion_history


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    out = parser.parse_args().output
    cfg = json.loads((out / 'effective_config.json').read_text(encoding='utf-8'))
    summary = json.loads((out / 'summary.json').read_text(encoding='utf-8'))
    runs = pd.read_csv(out / 'capacity_results.csv')
    aggregate = pd.read_csv(out / 'capacity_summary.csv')
    scenarios = pd.read_csv(out / 'capacity_scenarios.csv')
    od = pd.read_csv(out / 'OD.csv')
    net = Network(cfg)
    configure_plot()
    plots(net, runs, scenarios.to_dict('records'), pd.read_csv(out / 'selected_timeline.csv'), out)

    # Explicit planned and actual times for each direction of the same aircraft.
    inbound = od.movement == 'arrival'
    od['movement_ready_s'] = np.where(inbound, od.requested_s, od.landed_s + cfg['turnaround_s'])
    od['movement_started_s'] = np.where(inbound, od.admitted_s, od.departure_s)
    od['movement_finished_s'] = np.where(inbound, od.landed_s, od.exited_s)
    od['movement_status'] = np.where(od.movement_finished_s.notna(), 'completed',
                             np.where(od.movement_started_s.notna(), 'in_progress',
                               np.where(od.movement_ready_s.notna(), 'waiting', 'awaiting_landing')))
    od.to_csv(out / 'OD.csv', index=False, encoding='utf-8-sig')
    source_manifest = []
    for name in net.frames:
        for suffix in ('.shp', '.shx', '.dbf', '.prj', '.cpg'):
            path = net.source / (name + suffix)
            if path.exists():
                source_manifest.append(dict(path=str(path.relative_to(ROOT)), bytes=path.stat().st_size,
                                            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    csv(source_manifest, out / 'source_manifest.csv')
    csv(capacity_expansion_history(net, cfg), out / 'capacity_expansion_history.csv')
    holding_rows = []
    for name, hid in net.holding_for.items():
        h = net.holdings[hid]
        holding_rows.append(dict(terminal=name, holding_id=hid, **{k: v for k, v in h.items() if k != 'node'},
                                nominal_lap_s=(2 * math.pi * h['radius'] + 2 * h['straight']) / (h['speed'] * 1852 / 3600)))
    csv(holding_rows, out / 'holding_parameters.csv')
    fids = sorted(net.ports)
    anchor = net.nodes[net.ports[1]][:2]
    distances = np.array([np.linalg.norm(net.nodes[net.ports[f]][:2] - anchor) for f in fids])
    probs = np.exp(-distances / cfg['destination_distance_scale_m']); probs /= probs.sum()
    csv([dict(fid=f, distance_from_fid1_m=d, destination_probability=p)
         for f, d, p in zip(fids, distances, probs)], out / 'destination_priorities.csv')
    resource_load = np.zeros(len(net.conflicts))
    ring_load = np.zeros(2)
    pad_load = {fid: 0. for fid in fids}
    for name in net.terminals:
        origin_probability = .1 if name == '강원' else .3
        for fid, destination_probability in zip(fids, probs):
            weight = origin_probability * destination_probability
            for kind in ('arrival', 'departure'):
                route = net.paths[(kind, name, fid)]
                for claim in route['claims'][:-1]:
                    for resource in claim:
                        resource_load[resource] += weight * cfg['step_s']
                for ri in (0, 1):
                    ring_load[ri] += weight * np.count_nonzero(route['ring'][:-1] == ri) * cfg['step_s']
            arrival = net.paths[('arrival', name, fid)]
            # Under the implemented policy the pad is claimed before the final
            # corridor. These are minimum occupancies without any traffic delay.
            pad_load[fid] += weight * ((len(arrival['xyz']) - 1 - arrival['pad_gate']) * cfg['step_s']
                                       + cfg['turnaround_s'] + 100 / (cfg['speed_defaults_kts']['Short_Final'] * 1852 / 3600))
    resource_bounds = [dict(resource=f'crossing_{i}', type='crossing', expected_occupancy_s_per_arrival=x,
                            capacity=1, arrival_rate_upper_bound_per_h=3600 / x)
                       for i, x in enumerate(resource_load) if x > 0]
    resource_bounds += [dict(resource=f'pad_{fid}', type='pad', expected_occupancy_s_per_arrival=x,
                             capacity=1, arrival_rate_upper_bound_per_h=3600 / x) for fid, x in pad_load.items()]
    for scenario in scenarios.itertuples():
        for ri, cap in enumerate((scenario.outer_capacity, scenario.inner_capacity)):
            if ring_load[ri] > 0:
                resource_bounds.append(dict(resource=f'{scenario.scenario}_ring_{ri}', type='ring',
                    expected_occupancy_s_per_arrival=ring_load[ri], capacity=cap,
                    arrival_rate_upper_bound_per_h=3600 * cap / ring_load[ri]))
    resource_bounds = pd.DataFrame(resource_bounds).sort_values('arrival_rate_upper_bound_per_h')
    resource_bounds.to_csv(out / 'resource_upper_bounds.csv', index=False, encoding='utf-8-sig')
    selected_bound = resource_bounds[(resource_bounds.type != 'ring') |
        resource_bounds.resource.str.startswith(summary['selected_aggregate']['scenario'] + '_')].iloc[0]

    bounds = []
    for sep, group in aggregate.groupby('separation_m'):
        stable = group[group.stable_all_replicates]
        if len(stable):
            row = stable.sort_values('arrival_demand_per_h').iloc[-1]
            next_rates = group[group.arrival_demand_per_h > row.arrival_demand_per_h]
            next_test = next_rates.arrival_demand_per_h.min() if len(next_rates) else np.nan
            bounds.append(dict(separation_m=sep, largest_tested_stable_arrivals_per_h=row.arrival_demand_per_h,
                movements_mean_per_h=row.movements_mean_per_h, scenario=row.scenario,
                airborne_peak_in_stable_runs=int(row.peak_airborne_max),
                next_tested_arrivals_per_h=next_test,
                next_rate_any_stable=bool(next_rates[next_rates.arrival_demand_per_h == next_test].stable_all_replicates.any())))
        else:
            bounds.append(dict(separation_m=sep, largest_tested_stable_arrivals_per_h=np.nan,
                               movements_mean_per_h=np.nan, scenario='none', airborne_peak_in_stable_runs=np.nan,
                               next_tested_arrivals_per_h=np.nan, next_rate_any_stable=False))
    bounds = pd.DataFrame(bounds)
    bounds.to_csv(out / 'capacity_bounds.csv', index=False, encoding='utf-8-sig')
    selection = summary['selected_aggregate']
    selected = runs[(runs.scenario == selection['scenario']) &
                    np.isclose(runs.arrival_demand_per_h, selection['arrival_demand_per_h'])]
    tracks = pd.read_csv(out / 'selected_trajectories.csv')
    spacing_rows = []
    for (time_s, ring), group in tracks[tracks.ring >= 0].groupby(['time_s', 'ring']):
        if len(group) < 2:
            continue
        ring = int(ring)
        name = ('외곽_2800', '내곽_2000')[ring]
        curve = LineString(get_coordinates(net.frames[name].geometry.iloc[0])[:, :2])
        xyz = group[['x', 'y', 'z']].to_numpy()
        ids = group.aircraft_id.to_numpy()
        stations = np.array([curve.project(Point(p[:2])) for p in xyz])
        distances = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
        np.fill_diagonal(distances, np.inf)
        for i in range(len(group)):
            j = int(distances[i].argmin())
            delta = abs(stations[i] - stations[j])
            spacing_rows.append(dict(time_s=time_s, ring=name, occupancy=len(group), aircraft_id=int(ids[i]),
                nearest_aircraft_id=int(ids[j]), chord_distance_m=distances[i, j],
                ring_arc_distance_m=min(delta, curve.length - delta),
                required_distance_m=selection['separation_m']))
    spacing = pd.DataFrame(spacing_rows, columns=['time_s', 'ring', 'occupancy', 'aircraft_id',
        'nearest_aircraft_id', 'chord_distance_m', 'ring_arc_distance_m', 'required_distance_m'])
    spacing.to_csv(out / 'selected_spacing.csv', index=False, encoding='utf-8-sig')
    if len(spacing):
        assert (spacing.chord_distance_m >= spacing.required_distance_m - 1e-6).all()
        fig, ax = plt.subplots(figsize=(12, 5))
        for name, data in spacing.groupby('ring'):
            grouped_spacing = data.groupby('time_s').chord_distance_m.min()
            ax.plot(grouped_spacing.index / 60, grouped_spacing.values, '.-', ms=3, label=name)
        ax.axhline(selection['separation_m'], color='red', ls='--', label='요구 분리거리')
        ax.set_xlabel('시간 (분)'); ax.set_ylabel('링 내 가장 가까운 기체 중심 간 거리 (m)')
        ax.set_title('선정 조건의 실제 링 간격 · 10초 표본'); ax.grid(alpha=.2); ax.legend()
        fig.tight_layout(); fig.savefig(out / 'selected_spacing.png', dpi=180); plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.text(.5, .5, '선정 조건의 10초 표본에 같은 링의 기체 쌍이 없습니다.\n실제 링 간격 통계를 산출하지 않습니다.',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_axis_off()
        fig.tight_layout(); fig.savefig(out / 'selected_spacing.png', dpi=180); plt.close(fig)
    stage = od[(od.scenario == selection['scenario']) &
               np.isclose(od.demand_per_h, round(selection['arrival_demand_per_h'])) &
               (od.seed == summary['selected']['seed']) & inbound]
    # Target demand is an integer configuration value; generated demand can be
    # slightly smaller because the horizon contains a whole number of aircraft.
    if not len(stage):
        target_rate = min(cfg['arrival_rates_per_hour'], key=lambda r: abs(r - selection['arrival_demand_per_h']))
        stage = od[(od.scenario == selection['scenario']) & (od.demand_per_h == target_rate) &
                   (od.seed == summary['selected']['seed']) & inbound]
    matrix = stage.pivot_table(index='origin_name', columns='destination_fid', values='aircraft_id', aggfunc='count', fill_value=0)
    matrix.to_csv(out / 'selected_OD_matrix.csv', encoding='utf-8-sig')
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    image = axes[0].imshow(matrix.values, cmap='Blues')
    axes[0].set_xticks(range(len(matrix.columns)), matrix.columns)
    axes[0].set_yticks(range(len(matrix.index)), matrix.index)
    axes[0].set_xlabel('버티포트 fid'); axes[0].set_title('선정 조건의 도착 OD')
    for i in range(len(matrix)):
        for j in range(len(matrix.columns)):
            axes[0].text(j, i, str(matrix.iloc[i, j]), ha='center', va='center')
    axes[1].bar([str(f) for f in fids], probs, color='#287eb5')
    axes[1].set_title('fid 1 기준 거리별 도착 확률'); axes[1].set_xlabel('버티포트 fid'); axes[1].set_ylabel('확률')
    fig.tight_layout(); fig.savefig(out / 'OD_distribution.png', dpi=180); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for sep, data in scenarios.groupby('separation_m'):
        axes[0].plot(data.outer_capacity, data.inner_capacity, 'o-', label=f'{sep} m')
    axes[0].set_xlabel('외곽 설정 용량 (대)'); axes[0].set_ylabel('내곽 설정 용량 (대)')
    axes[0].set_title('비교한 용량 조합'); axes[0].legend(); axes[0].grid(alpha=.2)
    wait_cols = ['blocked_pad_steps', 'blocked_reservation_steps', 'blocked_ring_steps', 'blocked_separation_steps']
    values = selected[wait_cols].mean() * cfg['step_s']
    axes[1].bar(['착륙장', '교차구간 예약', '링 용량', '기체 분리'], values, color=['#8952ae', '#d78c28', '#287eb5', '#169b91'])
    axes[1].set_title('선정 조건의 기내 대기 원인 (반복 평균)'); axes[1].set_ylabel('누적 aircraft·s')
    fig.tight_layout(); fig.savefig(out / 'capacity_and_bottlenecks.png', dpi=180); plt.close(fig)
    stable_all = aggregate[aggregate.stable_all_replicates]
    stable_peak = int(stable_all.peak_airborne_max.max()) if len(stable_all) else None
    proof = dict(all_od_reachable=len(net.paths) == 64,
                 buffer_violation_count=int(runs.buffer_violations.sum()),
                 minimum_observed_swept_distance_m=float(runs.min_swept_separation_m.min()),
                 all_outer_occupancy_within_capacity=bool((runs.peak_outer <= runs.outer_capacity).all()),
                 all_inner_occupancy_within_capacity=bool((runs.peak_inner <= runs.inner_capacity).all()),
                 all_swept_ring_spacing_valid=bool((runs.min_ring_swept_separation_m.isna() |
                     (runs.min_ring_swept_separation_m >= runs.separation_m - 1e-6)).all()),
                 all_replicates_stable_scenario_count=int(aggregate.stable_all_replicates.sum()),
                 stable_peak_airborne=stable_peak)
    (out / 'validation_summary.json').write_text(json.dumps(proof, ensure_ascii=False, indent=2), encoding='utf-8')
    assert proof['buffer_violation_count'] == 0
    assert proof['all_outer_occupancy_within_capacity'] and proof['all_inner_occupancy_within_capacity']
    assert proof['all_swept_ring_spacing_valid']
    lines = [
        '# 네트워크 용량 분석', '',
        '## 결과 해석', '',
        ('- **선정 상태: 시험 범위 내 안정 조건 확인.**' if len(stable_all) else
         '- **선정 상태: 안정 조건 없음. 아래 선정 조건은 진단용이며 최대 운영 용량으로 해석할 수 없다.**'),
        f'- 분석 네트워크: 노드 {len(net.nodes)}개, 방향성 링크 {len(net.edges)}개, 도착·출발 OD {len(net.paths)}개 모두 연결됨.',
        f'- 실험: {len(scenarios)}개 분리거리·링 용량 조합 × {len(cfg["arrival_rates_per_hour"])}개 수요 × {cfg["replicates"]}개 난수 시드 = {len(runs)}개 조건.',
        f'- 1회 {cfg["duration_s"]/60:.0f}분, 앞 {cfg["warmup_s"]/60:.0f}분 제외. 시간 간격 {cfg["step_s"]}초.',
        f'- 모든 반복이 안정 기준을 만족한 조건 수: {len(stable_all)}개.',
        f'- 선정 조건: {selection["scenario"]}, 도착 수요 {selection["arrival_demand_per_h"]:.2f}대/h.',
        f'- 추적 대상으로 선정한 단일 반복의 완료율 판정 표본은 {summary["selected"]["mature_aircraft"]}대이다. '
        '작은 표본의 안정 판정은 장시간·다수 반복 실험으로 보강해야 한다.',
        f'- 선정 조건의 평균 착륙+이륙 처리량: **{selection["movements_mean_per_h"]:.1f}회/h** '
        f'(반복 범위 {selection["movements_min_per_h"]:.1f}~{selection["movements_max_per_h"]:.1f}회/h).',
        f'- 선정 조건의 동시 체공 관측 최댓값: **{selection["peak_airborne_max"]:.0f}대**. '
        + (f'모든 안정 조건에서의 관측 최댓값: {stable_peak}대.' if stable_peak is not None else '안정 조건의 체공 최댓값은 산출하지 못했다.'),
        f'- 전체 조건의 동시 체공 최댓값 {summary["max_observed_airborne"]}대는 혼잡·정체 조건을 포함하므로 운영 가능한 용량으로 해석하면 안 된다.',
        f'- 고정 OD 확률과 자유비행 점유시간에서 계산한 평균 자원부하 상한: '
        f'{selected_bound.arrival_rate_upper_bound_per_h:.1f}도착/h '
        f'(선정 조건의 병목 자원 {selected_bound.resource}). 대기 때문에 실제 처리량은 이보다 낮을 수 있다. 이 값은 확률 혼합의 장기 평균 상한이다.',
        '- 이 결과는 아래 속도·점유·교차 예약·정지 허용 가정 아래의 유한 실험값이다. 전역 최적 용량이나 인증된 안전 용량이 아니다.', '',
        '### 분리거리별 가장 높은 안정 수요', '',
        '| 링 분리거리(m) | 가장 높은 검증 수요(대/h) | 평균 처리량(회/h) | 해당 조건 체공 최고(대) | 다음 시험 수요(대/h) |',
        '|---:|---:|---:|---:|---:|']
    for _, row in bounds.iterrows():
        if pd.isna(row.largest_tested_stable_arrivals_per_h):
            lines.append(f'| {row.separation_m:g} | 안정 조건 없음 | — | — | — |')
            continue
        lines.append(f'| {row.separation_m:g} | {row.largest_tested_stable_arrivals_per_h:.2f} | '
                     f'{row.movements_mean_per_h:.1f} | {row.airborne_peak_in_stable_runs:g} | {row.next_tested_arrivals_per_h:.2f} |')
    lines += ['', '수요는 정수 대수로 생성하므로 15대/h 입력이 90분 동안 22대(14.67대/h)처럼 소폭 절삭될 수 있다. '
              '최대값은 시험한 수요 격자 안의 값이며 인접한 수요 구간은 추가 세분 실험이 필요하다. '
              '난수·교차 우선순위 때문에 결과가 단조 증가하지 않을 수 있다.', '',
              '## 적용 규칙', '',
              '1. 모든 자료를 EPSG:5186으로 변환해 XYZ 0.1m 이내 끝점을 연결했다. XY에서 교차해도 고도가 다르면 접속점을 만들지 않았다.',
              '2. 외곽 반시계, 내곽 시계. In은 하강/진입, Out은 상승/출발 방향. Short_Final은 양방향으로 사용하되 동일 기체의 착륙 뒤에만 이륙한다.',
              '3. 3D 길이/속도의 시간 가중치로 Dijkstra 최단시간 경로를 구했다. 혼잡 시 재경로 탐색은 하지 않고 진입을 기다린다.',
              '4. 각 기체를 반경 50m의 3D 구로 모델링해 중심 간 최소 100m를 적용했다. 링 위 두 기체에는 추가로 100/250/500/1,000m 조건을 적용했다.',
              '5. 한 시간 간격 내부의 상대 선형 운동 최소거리를 계산한다. 기록 시점에서만 떨어져 있는 교차도 탐지한다. 이동 궤적은 1초 단위 선형 근사이다.',
              '6. 20m 간격으로 샘플링한 서로 다른 비링 링크 사이 120m 이내를 잠재 교차로 지정하고, 실제 이동 중 상대 링크까지 150m 이내인 주변 구간을 독점 예약한다. '
              '같은 교차의 선진입 예약은 선점하지 않는다. 새 경합에는 착륙 최종접근, 링 진입, 기존 항로 순으로 우선권을 주고 동률은 네트워크 진입 시각을 적용한다.',
              '7. 교차 예약은 충돌 회피를 위한 보수적 제어 가정이다. 안전거리가 충족되어도 동시에 진입하지 못하는 경우가 있어 처리량을 낮출 수 있다. '
              '예약·분리 조건이 막히면 기체는 현재 위치에서 정지한다. 이러한 공중 정지는 속도/가속도/최소선회반경을 갖춘 비행역학 모델로 재검증해야 한다.',
              '8. 각 버티포트 용량 1. 최종접근에 필요한 교차구간의 첫 예약 전에 패드를 함께 확보해 착륙장-교차구간 간 순환 대기를 방지한다. 도착 후 20초 호버링을 완료한 뒤 착륙한다. '
              f'착륙 후 {cfg["turnaround_s"]}초 점유를 가정하고 이륙한다. 출발 기체가 패드에서 100m 이상 벗어나면 다음 패드 예약을 허용한다.',
              '9. 지상에 있는 기체는 패드 용량으로 관리하며 공중 50m 구 충돌 검사에는 포함하지 않는다. 같은 고도의 가까운 패드 간 동시 비행은 공중 거리 검사 대상이다.',
              '10. Holding은 가장 가까운 접근 시작점에 연결한 유한 대기열(각 Capacity=2, 총 8대)로 반영했다. 초과 수요는 외부 대기열에 남는다. '
              'Racetrack의 Radius/straight/Speed_kts로 계산한 이론 회전시간은 holding_parameters.csv에 기록했으나, 실제 궤도·진입 연결선·시간별 회전 위치는 모델링하지 않았다. '
              '따라서 체공 대수는 경로 위 기체만 세고 Holding은 별도로 기록한다.',
              '11. 출발지 비율은 강원 10%, 나머지 각각 30%로 10대 블록마다 섞는다. 마지막 불완전 블록은 비율이 다를 수 있다. '
              '목적지는 fid 1에서의 평면거리를 사용한 exp(-거리/500m) 확률이며 fid 1이 가장 높다. 귀환 목적지도 같은 출입점 비율로 생성한다.',
              '12. 링 초기 용량은 각각 14, 난수로 한쪽씩 1대 증가시킨 과정은 capacity_expansion_history.csv에 기록했다. '
              '모든 중간 1대 증분을 시뮬레이션한 것은 아니며 14/28/56 및 기하학적 한계 지점에서 비교했다.',
              '13. 동시 점유 한계는 원 위 중심 간 현 길이로 계산한 floor(π / asin(분리거리/(2r)))를 넘지 않도록 했다. '
              '1,000m에서 내곽 14대의 평균 원호 간격은 약 898m이므로 내곽은 12대로 제한된다.', '',
              '## 데이터에 없어 적용한 가정', '',
              '- 경로 SHP에는 Speed 필드가 없고 Holding_point에만 Speed_kts=90이 있다. '
              '접근·Expert_Off·링 90kt, Interchange·IAF↔FAF 60kt, Short_Final 전체 구간 30kt를 가정했다. '
              '실제 Speed/Speed_kts 필드가 추가되면 코드가 우선 읽는다(단위 kt 가정).',
              '- Expert_Off는 XY만 있어서 알려진 끝점의 609.6m 고도를 보간했다. Holding_point의 Z는 Altitude 필드를 썼다.',
              '- 착륙장 지상 점유 60초, 목적지 거리 척도 500m, 우선순위의 동률 처리, 교차 예약 범위는 모델 가정이며 설정/코드에 명시했다.',
              '- 짧은 최종접근의 가감속, 상승·하강률 한계, 바람, 에너지, 패드 이륙 준비 절차는 미포함이다.', '',
              '## 안정 판단과 검증', '',
              '- 준비 구간 이후부터 종료 1,200초 전까지 생성된 기체 중 왕복 출항 완료 비율 95% 이상.',
              '- 같은 코호트에서 완료 기체의 왕복 추가 지연 95백분위 300초 이하. 미완료 기체는 지연 평균에서 제외하지만 완료율 기준에서 실패 처리한다.',
              '- 준비 구간 이후 전체 미완료 왕복 대기열의 선형 증가율 2대/h 이하.',
              '- 모든 반복이 위 조건을 만족해야 안정 조건으로 선정한다. 안정 여부는 정해진 운항 시간 내의 실험 판정이며 무한 시간 안정성 증명은 아니다.',
              '- 기본 링 용량이 어느 순간에도 진입을 막지 않은 실행은 더 큰 용량에서도 시간 단계별 상태가 동일하다. '
              '이 경우 동일 결과를 재사용하고 capacity_results.csv의 computation에 identical_nonbinding_capacity로 표기했다.',
              f'- 전 조건의 50m 버퍼 위반: {proof["buffer_violation_count"]}건. 관측 최소 연속 분리거리: {proof["minimum_observed_swept_distance_m"]:.3f}m.',
              '- validate_capacity.py: 64개 OD, 링 방향, XYZ 접속, 시간 사이 교차, 수직 분리, OD 비율, 착륙→이륙·패드 점유 검증.', '',
              '## 재실행', '',
              '```powershell', 'python -B -X utf8 Network/코드/validate_capacity.py',
              'python -B -X utf8 Network/코드/network_capacity.py',
              'python -B -X utf8 Network/코드/report_capacity.py', '```', '',
              '속도·점유시간·수요·분리거리·시드는 Network/코드/capacity_config.json에서 변경한다. 원본 SHP는 수정하지 않는다.', '',
              '## 생성물', '',
              '- OD.csv: 모든 조건의 기체별 왕복 OD, 계획/실제 시각, 상태. 출입점은 Name, 버티포트는 fid.',
              '- OD_routes.csv, network_edges.csv, network_nodes.csv: 경로·시간·속도 출처·XYZ 그래프.',
              '- capacity_results.csv, capacity_summary.csv, capacity_bounds.csv: 개별 반복·반복 통합·분리거리별 결과.',
              '- selected_timeline.csv, selected_trajectories.csv(10초 간격), selected_events.csv: 선정 조건의 추적 자료.',
              '- selected_spacing.csv, selected_spacing.png: 기체별 최근접 기체와 실제 3D 현 거리·링 원호 간격(10초 간격).',
              '- network_3d.png, capacity_comparison.png, selected_timeline.png, OD_distribution.png, capacity_and_bottlenecks.png: 시각화.',
              '- source_manifest.csv, input_audit.csv, conflict_corridors.csv, holding_parameters.csv, destination_priorities.csv: 데이터·모델 근거.',
              '- resource_upper_bounds.csv: 고정 OD 혼합·최소 점유시간으로 계산한 자원별 장기 평균 도착률 상한.',
              '- 검증/ 폴더: 짧은 기능 시험 결과로, 본 분석 결과와 구분한다.']
    if (out / '순환확장/순환분석보고서.md').exists():
        lines[2:2] = ['**추가 순환 실험:** 이 문서의 체공 관측값은 안정 왕복 조건의 통계이며 '
                      '링 자체의 최대 순환 대수가 아니다. 계속 생성·진입시키는 확장 실험은 '
                      '[순환분석보고서](순환확장/순환분석보고서.md)에서 확인한다.', '']
    (out / '분석보고서.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(json.dumps(proof, ensure_ascii=False, indent=2))
    print(bounds.to_string(index=False))


if __name__ == '__main__':
    main()
