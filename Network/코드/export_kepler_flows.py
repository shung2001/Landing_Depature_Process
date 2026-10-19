"""Export the best tested productive OD run and native Kepler Trip animations."""
from pathlib import Path
import json
import pandas as pd

from network_capacity import Network, OUTPUT, demand, simulate, csv, json_text
from circulation_capacity import kepler_csv


def trips(frame, path):
    features = []
    # Keep arrival, ground dwell, and departure distinct; never draw a straight
    # artificial line across a missing phase or between different aircraft.
    for (aid, phase), group in frame.groupby(['aircraft_id', 'phase'], sort=True):
        group = group.sort_values('time_s').drop_duplicates('time_s')
        if len(group) < 2:
            continue
        coordinates = group[['longitude', 'latitude', 'altitude_m', 'timestamp_ms']].to_numpy()
        coordinates[:, 3] /= 1000
        features.append(dict(type='Feature', properties=dict(aircraft_id=str(aid), phase=phase,
            scenario=group.scenario.iloc[0]), geometry=dict(type='LineString', coordinates=coordinates.tolist())))
    path.write_text(json.dumps(dict(type='FeatureCollection', features=features),
                               ensure_ascii=False, allow_nan=False), encoding='utf-8')
    return len(features)


def main():
    out = OUTPUT / '순환확장'
    cfg = json.loads((out / 'effective_config.json').read_text(encoding='utf-8'))
    runs = pd.read_csv(out / 'full_OD_stress_results.csv')
    cfg['arrival_extra_ring_laps'] = 0
    net = Network(cfg)
    # Include balanced demand as well as saturated demand: higher input can
    # reduce output when the current intersection controller gridlocks.
    balanced = []
    for rate in (5, 10, 15, 30, 60):
        for seed in (cfg['seed'], cfg['seed'] + 1):
            name = f'OD_laps0_R{rate}_O14_I14'
            metrics, *_ = simulate(net, cfg, demand(net, cfg, rate, seed), (14, 14), 100, name, seed)
            metrics.update(extra_laps=0, target_arrivals_per_h=rate)
            balanced.append(metrics)
            print(f'Throughput search: {name}, seed={seed}, movements={metrics["movements_per_h"]:.1f}/h', flush=True)
    candidates = pd.concat([runs[runs.extra_laps == 0], pd.DataFrame(balanced)], ignore_index=True)
    candidates.to_csv(out / 'throughput_search_results.csv', index=False, encoding='utf-8-sig')
    # The operational objective here is completed landing+takeoff throughput,
    # then more moving aircraft and fewer waiting aircraft. Restrict to ordinary
    # shortest-time paths so additional diagnostic laps are not called optimal.
    selected = candidates.sort_values(
        ['movements_per_h', 'external_exits_per_h', 'end_backlog', 'seed', 'scenario'],
        ascending=[False, False, True, True, True]).iloc[0]
    rate, seed = int(selected.target_arrivals_per_h), int(selected.seed)
    metrics, history, states, positions, events = simulate(net, cfg, demand(net, cfg, rate, seed),
        (int(selected.outer_capacity), int(selected.inner_capacity)), selected.separation_m,
        selected.scenario, seed, detail=True)
    kepler_csv(positions, out / 'kepler_best_throughput_flow.csv', selected.scenario, seed)
    csv(states, out / 'best_throughput_OD.csv')
    csv(events, out / 'best_throughput_events.csv')
    history.to_csv(out / 'best_throughput_timeline.csv', index=False, encoding='utf-8-sig')
    csv(net.route_rows, out / 'best_throughput_routes.csv')
    exported = {}
    for stem in ('kepler_ring_flow', 'kepler_OD_flow', 'kepler_best_throughput_flow'):
        frame = pd.read_csv(out / (stem + '.csv'))
        exported[stem] = dict(rows=len(frame), aircraft=int(frame.aircraft_id.nunique()),
                             trip_segments=trips(frame, out / (stem + '.geojson')))
    (out / 'kepler_export_summary.json').write_text(json_text(dict(
        selection_objective='highest_tested_shortest_route_landing_takeoff_throughput',
        selected=metrics, exports=exported)), encoding='utf-8')
    guide = [
        '# Kepler 이동 흐름 사용법', '',
        '## 먼저 볼 파일', '',
        '1. **kepler_ring_flow.csv**: 계속 수요를 생성해 채운 뒤 실제로 순환하는 링 전용 실험. '
        '링별 마지막 15분 최소 동시 이동 대수 합을 최대로 하는 시험 조건을 선정했다.',
        '2. **kepler_best_throughput_flow.csv**: 접근선→버티포트→출항 전체 OD. '
        '추가 대기 바퀴가 없는 최단시간 경로 실험 중 착륙+이륙 처리량이 가장 높은 조건이다.',
        '3. **kepler_OD_flow.csv**: 추가 순환을 포함한 전체 OD 실험 중 동시 순환 지속을 우선한 조건이다.', '',
        f'- 처리량 우선 조건: {selected.scenario}, seed={seed}, 입력 수요 {rate}대/h.',
        f'- 착륙+이륙 {metrics["movements_per_h"]:.1f}회/h, '
        f'왕복 완료 출항 {metrics["external_exits_per_h"]:.1f}대/h, 안정 판정 {metrics["stable"]}.',
        '- 시험한 조건에서 목적별로 선정했다. 전역 최적해나 안정성이 입증된 최대 수용량을 뜻하지 않는다.', '',
        '## CSV로 재생', '',
        '1. Kepler에서 Add Data로 원하는 CSV를 불러온다.',
        '2. Point 레이어의 Latitude=latitude, Longitude=longitude를 지정한다.',
        '3. 색상은 aircraft_id 또는 ring, 툴팁은 aircraft_id·phase·time_s·altitude_m를 지정한다.',
        '4. Filters에서 timestamp를 선택한다. 시간창을 링 실험은 5~10초, 전체 OD는 10~20초로 좁히고 재생한다.',
        '5. kepler_ring_paths.csv를 함께 추가하면 geometry 열의 GeoJSON 선으로 원본 링을 표시할 수 있다.', '',
        'CSV는 시간별 포인트이다. 과거 점을 모두 켜두면 한 기체가 여러 점으로 보이므로 '
        '현재 동시 대수를 보려면 좁은 시간창을 사용한다. 누적 표시 대수를 동시 체공 대수로 해석하지 않는다.', '',
        '## 움직이는 궤적과 꼬리 표시', '',
        '동일한 이름의 .geojson 파일도 제공한다. Kepler Trip 레이어용 '
        '[경도, 위도, 고도, Unix 초] 좌표이다. .geojson을 불러와 Trip 레이어로 표시하고 '
        'Trail Length를 약 10~30초로 두면 기체별 이동 방향을 볼 수 있다. '
        '도착과 출발은 다른 세그먼트로 분리했다.', '',
        '## 시각과 좌표', '',
        '- 좌표: EPSG:4326, longitude/latitude는 도, altitude_m는 원본 SHP 고도(m).',
        '- timestamp: UTC ISO 8601. timestamp_ms: Unix 밀리초.',
        '- 시뮬레이션 0초=2026-10-20 00:00:00 UTC로 정한 가상 시각이며 실측 운항 기록이 아니다.',
        '- 링 전용은 5초 표본과 진입 시각, 전체 OD는 10초 표본이다.', '',
        '공식 문서: [CSV 형식](https://docs.kepler.gl/docs/user-guides/b-kepler-gl-workflow/a-add-data-to-the-map), '
        '[시간 재생](https://docs.kepler.gl/docs/user-guides/h-playback), '
        '[Trip 레이어](https://docs.kepler.gl/docs/user-guides/c-types-of-layers/k-trip).', '',
        '재생성: `python -B -X utf8 Network/코드/export_kepler_flows.py`', '']
    (out / 'Kepler_사용법.md').write_text('\n'.join(guide), encoding='utf-8')
    report = out / '순환분석보고서.md'
    text = report.read_text(encoding='utf-8').split('\n## 처리량 우선 추가 탐색')[0]
    text += ('\n## 처리량 우선 추가 탐색\n\n'
             '과부하 조건만으로 최적 조건을 정하지 않도록 추가 바퀴 없는 기본 용량에서 '
             '5/10/15/30/60대/h도 두 시드로 시험했다. 같은 1시간·준비 10분 조건으로 비교했다.\n\n'
             f'- 처리량 최대: {selected.scenario}, seed={seed}, {metrics["movements_per_h"]:.1f}회/h.\n'
             f'- 안정 기준 만족: {metrics["stable"]}. 처리량 최대와 안정 운영 가능 여부는 다른 판단이다.\n'
             '- 상세 비교: throughput_search_results.csv.\n'
             '- Kepler 파일 선택·재생 방법: [Kepler_사용법](Kepler_사용법.md).\n')
    report.write_text(text, encoding='utf-8')
    print(json_text(dict(selected=metrics, exports=exported)))


if __name__ == '__main__':
    main()
