"""Export full arrival/landing/turnaround/departure experiments, 10 to 60 aircraft."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from network_capacity import Network, OUTPUT, demand, simulate, csv
from circulation_capacity import kepler_csv
from export_kepler_flows import trips


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT / '고정대수_20_60')
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(Path(__file__).with_name('capacity_config.json').read_text(encoding='utf-8'))
    cfg.update(duration_s=5400, warmup_s=0, arrival_extra_ring_laps=0)
    net = Network(cfg)
    # The capacity is a ceiling, not a forced occupancy. Keep it identical in
    # all six runs so the comparison changes only the requested fleet size.
    caps, separation = (60, 60), 100
    results = []
    for count in range(10, 61, 10):
        name = f'full_OD_{count}aircraft'
        request_cfg = dict(cfg, duration_s=3600)
        requests = demand(net, request_cfg, count, cfg['seed'], count=count)
        metrics, history, states, positions, events = simulate(
            net, cfg, requests, caps, separation, name, cfg['seed'], detail=True)
        state = pd.DataFrame(states)
        assert len(state) == count and state.aircraft_id.nunique() == count
        # Verify observed phase ordering for every aircraft that reached it.
        for earlier, later in [('requested_s', 'admitted_s'), ('admitted_s', 'landed_s'),
                               ('landed_s', 'departure_s'), ('departure_s', 'exited_s')]:
            reached = state[later].notna()
            assert state.loc[reached, earlier].notna().all()
            assert (state.loc[reached, later] >= state.loc[reached, earlier]).all()
        assert metrics['buffer_violations'] == 0
        state['status'] = np.select(
            [state.exited_s.notna(), state.departure_s.notna(), state.landed_s.notna(),
             state.admitted_s.notna()],
            ['completed', 'departing', 'grounded', 'arriving'], default='waiting')
        metrics.update(total_input_aircraft=count, injection_duration_s=3600,
                       observation_duration_s=5400,
                       admitted_total=int(state.admitted_s.notna().sum()),
                       landed_total=int(state.landed_s.notna().sum()),
                       departed_total=int(state.departure_s.notna().sum()),
                       completed_total=int(state.exited_s.notna().sum()),
                       incomplete_total=int(state.exited_s.isna().sum()))
        results.append(metrics)
        state.to_csv(out / f'{name}_aircraft.csv', index=False, encoding='utf-8-sig')
        history.to_csv(out / f'{name}_timeline.csv', index=False, encoding='utf-8-sig')
        csv(events, out / f'{name}_events.csv')
        frame = kepler_csv(positions, out / f'{name}_flow.csv', name, cfg['seed'])
        trips(frame, out / f'{name}_flow.geojson')
        csv(results, out / 'full_operation_summary.csv')
        print(f'{name}: landed={metrics["landed_total"]}, completed={metrics["completed_total"]}, '
              f'incomplete={metrics["incomplete_total"]}, peak_airborne={metrics["peak_airborne"]}', flush=True)
    csv(net.route_rows, out / 'OD_routes.csv')
    (out / 'effective_config.json').write_text(json.dumps(dict(
        config=cfg, counts=list(range(10, 61, 10)), injection_duration_s=3600,
        drain_duration_s=1800, ring_capacities=caps, separation_m=separation,
        experiment='full_OD'), ensure_ascii=False, indent=2), encoding='utf-8')
    # Write the report without an optional markdown-table dependency.
    lines = ['# 전체 운항 10~60대 비교', '',
             '기존 독립 링 순환 결과를 접근→링→착륙→지상 체류→이륙→출항 전체 OD 결과로 교체했다.', '',
             '- 총 투입 기체 수: 10·20·30·40·50·60대. 각 기체가 한 번 도착하고 같은 기체가 출발한다.',
             '- 0~3600초에 분산 투입하고 5400초까지 관찰한다. 이후 1800초는 추가 투입이 없다.',
             '- 전 조건 동일 시드, 링별 한도 60대, 분리거리 100m, 추가 순환 0바퀴. 한도는 실제 점유 대수가 아니다.',
             '- 원본 접근선·최단시간 OD·패드·교차 예약·분리 검사·20초 호버링·60초 지상 체류를 적용했다.',
             '- 전체 운항 절차를 모델링하지만 혼잡하면 끝까지 완료하지 못할 수 있다. 정체 기체를 강제로 완료시키지 않았다.', '',
             '| 총 투입(대) | 착륙(대) | 이륙(대) | 출항 완료(대) | 미완료(대) | 최대 동시 체공(대) |',
             '|---:|---:|---:|---:|---:|---:|']
    for r in results:
        lines.append('| ' + ' | '.join(str(r[k]) for k in (
            'total_input_aircraft', 'landed_total', 'departed_total', 'completed_total',
            'incomplete_total', 'peak_airborne')) + ' |')
    lines += ['', '## 파일과 재생', '',
              '- full_OD_Naircraft_flow.geojson: Kepler Trip 레이어용. N이 총 투입 대수다. 원하는 조건 하나를 켜고 시간 재생한다.',
              '- full_OD_Naircraft_flow.csv: 10초 간격 공중 경로 표본. 지상 체류·외부 대기·Holding은 공간 궤적으로 표시하지 않는다.',
              '- full_OD_Naircraft_aircraft.csv: 기체별 요청·진입·착륙·이륙·출항 시각과 종료 상태. 빈 시각은 미도달이다.',
              '- full_OD_Naircraft_timeline.csv: 시간별 체공·링 이동·지상·Holding·외부 대기·누적 착륙·출항 대수.',
              '- full_operation_summary.csv: 여섯 조건의 비교 결과. landed_total 등은 90분 전체 누적이며 per_h는 90분 관찰 시간 기준이다.',
              '- 시간 원점은 2026-10-20 00:00:00 UTC로 설정한 가상 시각이다.', '',
              '총 60대 투입은 60대 동시 비행과 다르다. 기존 공중 정지 허용, 보수적 교차 예약, Holding 공간 궤적 미포함 등의 한계가 남는다. 단일 시드 실험이며 실제 운항 용량 인증이나 최적해가 아니다.', '',
              '검증: 조건별 요청 대수와 고유 ID, 요청→진입→착륙→이륙→출항 시각 순서, 기존 시뮬레이터의 연속 분리 검사.', '',
              '## 재실행', '', '```powershell',
              'python -B -X utf8 Network/코드/full_operation_fleet.py', '```', '',
              '## 학습 자료', '',
              '- [SimPy 공유 자원](https://simpy.readthedocs.io/en/stable/simpy_intro/shared_resources.html)',
              '- [NASA 용량·처리량 분석](https://ntrs.nasa.gov/api/citations/20205001421/downloads/AIAA_Vertiport_Analysis_final_20200507.pdf?attachment=true)',
              '- [OR-Tools 스케줄링](https://developers.google.com/optimization/scheduling/job_shop)']
    (out / '실험과_학습자료.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
