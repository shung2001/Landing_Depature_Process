"""Separate SHP ring circulation, sustained admission, and full-OD stress tests.

python -B -X utf8 Network/코드/circulation_capacity.py --workers 4
Ring-only experiments isolate circulation from pads and interchange reservations.
They do not certify the same capacity for connected arrival/departure operations.
"""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer
import matplotlib.pyplot as plt

from network_capacity import (Network, OUTPUT, KT, RINGS, length, demand, simulate,
                              csv, json_text, configure_plot)


class RingTrack:
    """Directed original SHP polyline, parameterised by arc length in metres."""
    def __init__(self, net, ring):
        edges = [e for e in net.edges if e.ring == ring]
        following = {e.u: e for e in edges}
        first, edge = edges[0], edges[0]
        points, visited = [], set()
        while edge.id not in visited:
            visited.add(edge.id)
            points.extend(edge.xyz[:-1])
            edge = following[edge.v]
        if edge.id != first.id or len(visited) != len(edges):
            raise ValueError('Ring must be one closed directed cycle')
        points.append(first.xyz[0])
        self.xyz = np.asarray(points)
        self.arc = length(self.xyz)
        self.length = float(self.arc[-1])
        self.radius = self.length / (2 * math.pi)
        self.speed = min(e.speed for e in edges) * KT
        self.ring = ring
        self.name = RINGS[ring]
        self.entries = {}
        # External approach endpoints for outer; actual Interchange_In endpoints
        # nearest each external direction for inner. This is a ring-entry test;
        # the connector flight itself is not simulated here.
        inner = [net.nodes[e.v] for e in net.edges if e.layer == 'Interchange_In']
        for name in net.terminals:
            path = net.paths[('arrival', name, 1)]['path']
            outer = next(net.nodes[net.edges[e].v] for e in path if net.edges[e].layer == '접근')
            point = outer if ring == 0 else min(inner, key=lambda p: np.linalg.norm(p[:2] - outer[:2]))
            self.entries[name] = float(self.arc[np.linalg.norm(self.xyz - point, axis=1).argmin()])

    def at(self, station):
        s = np.asarray(station) % self.length
        return np.column_stack([np.interp(s, self.arc, self.xyz[:, k]) for k in range(3)])

    def geometric_max(self, separation):
        return int(math.floor(math.pi / math.asin(separation / (2 * self.radius))))


def pair_minimum(old, new):
    if len(old) < 2:
        return math.inf
    a, b = np.triu_indices(len(old), 1)
    relative = old[a] - old[b]
    velocity = (new[a] - old[a]) - (new[b] - old[b])
    norm = np.einsum('ij,ij->i', velocity, velocity)
    at = np.clip(-np.einsum('ij,ij->i', relative, velocity) / np.maximum(norm, 1e-20), 0, 1)
    return float(np.linalg.norm(relative + at[:, None] * velocity, axis=1).min())


def preload(track, separation, count, placement, seed, duration=900, step=1):
    """Fly an existing fleet for >2 outer laps; no admission or pad bottleneck."""
    rng = np.random.default_rng(seed)
    required_arc = 2 * track.radius * math.asin(separation / (2 * track.radius))
    if placement == 'uniform':
        stations = np.arange(count) * track.length / count
    else:
        spare = track.length - count * required_arc
        if spare < 0:
            return dict(ring=track.name, separation_m=separation, count=count, placement=placement,
                        seed=seed, valid=False, reason='geometric_limit', min_swept_m=None,
                        laps_per_aircraft=0., moving_aircraft=0)
        gaps = required_arc + rng.dirichlet(np.ones(count) * 2) * spare
        stations = np.r_[0., np.cumsum(gaps[:-1])] + rng.uniform(0, track.length)
    minimum = pair_minimum(track.at(stations), track.at(stations))
    valid = minimum >= separation - 1e-6
    if valid:
        for t in np.arange(0., duration, step):
            old = track.at(stations + track.speed * t)
            new = track.at(stations + track.speed * min(t + step, duration))
            minimum = min(minimum, pair_minimum(old, new))
            if minimum < separation - 1e-6:
                valid = False
                break
    return dict(ring=track.name, separation_m=separation, count=count, placement=placement,
                seed=seed, valid=valid, reason='passed' if valid else 'spacing_violation',
                min_swept_m=minimum, laps_per_aircraft=track.speed * duration / track.length if valid else 0.,
                moving_aircraft=count if valid else 0)


def feed(track, net, cfg, separation, capacity, rate, seed, mode, detail=False):
    """Keep generating demand while admitted aircraft continuously circulate.

    All aircraft have the same SHP/default speed on a given ring. Admission needs
    a gap both ahead and behind; nobody is teleported, redistributed, or stopped.
    This experiment deliberately has no sink and tests sustained occupancy.
    """
    duration, step = cfg['duration_s'], cfg['step_s']
    requests = demand(net, cfg, rate, seed)
    queues = {name: deque() for name in net.terminals}
    phases, admitted, request_ids = [], [], []
    records, history, states = [], [], []
    cursor, minimum, ten_seconds = 0, math.inf, 0.
    required_arc = 2 * track.radius * math.asin(separation / (2 * track.radius)) + 2.
    fixed = track.entries[next(iter(net.terminals))]
    scenario = f'{mode}_S{separation}_C{capacity}_R{rate}_seed{seed}'

    def record(t):
        if not detail or not phases:
            return
        xyz = track.at(np.asarray(phases) + track.speed * t)
        for i, p in enumerate(xyz):
            records.append(dict(time_s=t, aircraft_id=f'{track.ring}_{request_ids[i]}',
                x=p[0], y=p[1], z=p[2], ring=track.ring, phase='circulation',
                distance_m=track.speed * (t - admitted[i]), scenario=scenario))

    for tick in range(round(duration / step) + 1):
        t = tick * step
        while cursor < len(requests) and requests[cursor]['requested_s'] <= t:
            request = requests[cursor]
            queues[request['origin']].append(request)
            cursor += 1
        # Oldest eligible head among access queues; a blocked access need not
        # prevent an independent access point from admitting an aircraft.
        for name in sorted(queues, key=lambda n: queues[n][0]['requested_s'] if queues[n] else math.inf):
            if not queues[name] or len(phases) >= capacity:
                continue
            station = fixed if mode == 'single_entry' else track.entries[name]
            if phases:
                delta = (np.asarray(phases) + track.speed * t - station) % track.length
                if np.minimum(delta, track.length - delta).min() < required_arc:
                    continue
            request = queues[name].popleft()
            phases.append(station - track.speed * t)
            admitted.append(t)
            request_ids.append(request['aircraft_id'])
            states.append(dict(**request, admitted_s=t, ring=track.name, scenario=scenario))
            if detail:
                record(t)
        if tick % max(1, round(5 / step)) == 0 or t == duration:
            history.append(dict(time_s=t, ring=track.name, generated=cursor, moving=len(phases),
                                external_queue=cursor - len(phases), scenario=scenario))
            record(t)
        if t >= duration:
            break
        old = track.at(np.asarray(phases) + track.speed * t)
        new = track.at(np.asarray(phases) + track.speed * min(duration, t + step))
        minimum = min(minimum, pair_minimum(old, new))
        if minimum < separation - 1e-6:
            raise AssertionError(f'{track.name}/{scenario}: unsafe admission or circulation')
        if len(phases) >= 10:
            ten_seconds += step
    tail = [h['moving'] for h in history if h['time_s'] >= .75 * duration]
    metrics = dict(ring=track.name, scenario=scenario, separation_m=separation, capacity=capacity,
        rate_per_h=rate, seed=seed, mode=mode, generated=len(requests), admitted=len(phases),
        peak_moving=len(phases), tail_min_moving=min(tail), seconds_at_least_10_moving=ten_seconds,
        completed_laps=sum(int(track.speed * (duration - start) / track.length) for start in admitted),
        end_queue=len(requests) - len(phases), min_swept_m=minimum if math.isfinite(minimum) else None,
        buffer_violations=0)
    return metrics, history, records, states


def kepler_csv(rows, path, experiment, seed):
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError('No trajectory samples to export')
    if 'scenario' not in frame:
        frame['scenario'] = experiment
    # Admission events can coincide with scheduled samples.
    frame = frame.drop_duplicates(['scenario', 'time_s', 'aircraft_id'], keep='last')
    lon, lat = Transformer.from_crs(5186, 4326, always_xy=True).transform(frame.x, frame.y)
    frame['longitude'], frame['latitude'], frame['altitude_m'] = lon, lat, frame.z
    frame['timestamp'] = (pd.Timestamp('2026-10-20T00:00:00Z') +
                          pd.to_timedelta(frame.time_s, unit='s')).dt.strftime('%Y-%m-%dT%H:%M:%SZ')
    frame['timestamp_ms'] = (pd.to_datetime(frame.timestamp, utc=True).astype('int64') // 1_000_000)
    frame['experiment'], frame['seed'] = experiment, seed
    frame['aircraft_id'] = experiment + '_' + frame.aircraft_id.astype(str)
    frame.sort_values(['aircraft_id', 'time_s']).to_csv(path, index=False, encoding='utf-8-sig')
    return frame


def init_worker(cfg):
    global NET, CFG, TRACKS
    CFG, NET = cfg, Network(cfg)
    TRACKS = [RingTrack(NET, ri) for ri in (0, 1)]


def run_job(job):
    kind, options = job
    if kind == 'preload':
        ri, sep, count, placement, seed = options
        return kind, preload(TRACKS[ri], sep, count, placement, seed)
    if kind == 'feed':
        ri, sep, cap, rate, seed, mode = options
        metrics, history, _, _ = feed(TRACKS[ri], NET, CFG, sep, cap, rate, seed, mode)
        return kind, (metrics, history)
    raise ValueError(kind)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--output', type=Path, default=OUTPUT / '순환확장')
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(Path(__file__).with_name('capacity_config.json').read_text(encoding='utf-8'))
    cfg.update(duration_s=3600, warmup_s=600)
    net = Network(cfg)
    tracks = [RingTrack(net, ri) for ri in (0, 1)]
    jobs = []
    for track in tracks:
        for sep in (100, 250, 500, 1000):
            maximum = track.geometric_max(sep)
            counts = sorted(set([10, 14, 20, 28, 56, maximum, maximum + 1]))
            for count in counts:
                for placement in ('uniform', 'random_gaps'):
                    jobs.append(('preload', (track.ring, sep, count, placement, cfg['seed'])))
            for cap in sorted(set([min(14, maximum), min(28, maximum), maximum])):
                for rate in (120, 360, 720):
                    for mode in ('single_entry', 'weighted_entries'):
                        for seed in (cfg['seed'], cfg['seed'] + 1):
                            jobs.append(('feed', (track.ring, sep, cap, rate, seed, mode)))
    preloads, feeds, timelines = [], [], []
    with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=(cfg,)) as pool:
        futures = [pool.submit(run_job, job) for job in jobs]
        for i, future in enumerate(as_completed(futures), 1):
            kind, result = future.result()
            if kind == 'preload':
                preloads.append(result)
            else:
                metrics, history = result
                feeds.append(metrics)
                timelines.extend(history)
            if i % 20 == 0 or i == len(jobs):
                print(f'Ring experiments {i}/{len(jobs)}', flush=True)
    csv(preloads, out / 'preloaded_circulation.csv')
    csv(feeds, out / 'continuous_feed_results.csv')
    csv(timelines, out / 'continuous_feed_timeline.csv')
    preloads, feeds = pd.DataFrame(preloads), pd.DataFrame(feeds)

    # Select a common rate/seed/mode/spacing for both rings, with capacity at
    # the geometric ceiling. Both tracks are independent ring-only experiments.
    ceiling = feeds[feeds.apply(lambda r: r.capacity == tracks[RINGS.index(r.ring)].geometric_max(r.separation_m), axis=1)]
    choice = ceiling.groupby(['separation_m', 'rate_per_h', 'seed', 'mode']).agg(
        sustained=('tail_min_moving', 'sum'), peak=('peak_moving', 'sum')).reset_index().sort_values(
            ['sustained', 'peak', 'separation_m', 'rate_per_h', 'seed', 'mode'],
            ascending=[False, False, True, True, True, True]).iloc[0]
    ring_tracks, selected_timeline, selected_states, selected_metrics = [], [], [], []
    for track in tracks:
        metrics, history, records, states = feed(track, net, cfg, float(choice.separation_m),
            track.geometric_max(choice.separation_m), int(choice.rate_per_h), int(choice.seed), choice['mode'], True)
        ring_tracks.extend(records)
        selected_timeline.extend(history)
        selected_states.extend(states)
        selected_metrics.append(metrics)
    kepler_csv(ring_tracks, out / 'kepler_ring_flow.csv', 'ring_only', int(choice.seed))
    csv(selected_states, out / 'ring_admissions_OD.csv')
    csv(selected_timeline, out / 'selected_ring_timeline.csv')

    # Compare the existing full-OD policy with actual extra laps in each visited
    # ring. This exposes freezes instead of calling parked aircraft circulation.
    od_metrics, od_histories = [], []
    best_key, best_detail, best_cfg = None, None, None
    for extra in (0, 2):
        options = dict(cfg, arrival_extra_ring_laps=extra)
        model = Network(options)
        for rate in (120, 360):
            for caps in ((14, 14), (28, 28)):
                for seed in (cfg['seed'], cfg['seed'] + 1):
                    name = f'OD_laps{extra}_R{rate}_O{caps[0]}_I{caps[1]}'
                    metrics, history, states, _, _ = simulate(model, options, demand(model, options, rate, seed),
                        caps, 100, name, seed)
                    metrics.update(extra_laps=extra, target_arrivals_per_h=rate)
                    od_metrics.append(metrics)
                    history['scenario'], history['seed'] = name, seed
                    od_histories.extend(history.to_dict('records'))
                    # Prefer sustained simultaneous movement over stopped peaks,
                    # then completed movements. Deterministic ties retain first.
                    score = (metrics['outer_seconds_at_least_10_moving'] + metrics['inner_seconds_at_least_10_moving'],
                             metrics['peak_moving_outer'] + metrics['peak_moving_inner'], metrics['movements_per_h'])
                    if best_key is None or score > best_key:
                        best_key = score
                        best_detail = (name, seed, rate, caps)
                        best_cfg = options
                    print(f'{name} seed={seed}: occupied={metrics["peak_outer"]}/{metrics["peak_inner"]} '
                          f'moving={metrics["peak_moving_outer"]}/{metrics["peak_moving_inner"]}', flush=True)
                    csv(od_metrics, out / 'full_OD_stress_results.csv')
    csv(od_metrics, out / 'full_OD_stress_results.csv')
    csv(od_histories, out / 'full_OD_stress_timeline.csv')
    name, seed, rate, caps = best_detail
    model = Network(best_cfg)
    metrics, history, states, positions, events = simulate(model, best_cfg,
        demand(model, best_cfg, rate, seed), caps, 100, name, seed, detail=True)
    kepler_csv(positions, out / 'kepler_OD_flow.csv', name, seed)
    csv(states, out / 'selected_full_OD.csv')
    csv(events, out / 'selected_full_OD_events.csv')
    history.to_csv(out / 'selected_full_OD_timeline.csv', index=False, encoding='utf-8-sig')
    csv(model.route_rows, out / 'selected_full_OD_routes.csv')
    export_ring_paths(tracks, out)
    summary = dict(ring_experiments=len(jobs), full_OD_experiments=len(od_metrics),
        selected_feed=choice.to_dict(), selected_ring_metrics=selected_metrics,
        selected_full_OD=dict(scenario=name, seed=seed, rate=rate, caps=caps,
                             extra_laps=best_cfg['arrival_extra_ring_laps'], metrics=metrics),
        preload_max=preloads[preloads.valid].groupby(['ring', 'separation_m'])['count'].max().reset_index().to_dict('records'))
    (out / 'summary.json').write_text(json_text(summary), encoding='utf-8')
    (out / 'effective_config.json').write_text(json_text(cfg), encoding='utf-8')
    report(preloads, feeds, pd.DataFrame(od_metrics), pd.DataFrame(selected_timeline), summary, out)
    print(json_text(summary), flush=True)


def export_ring_paths(tracks, out):
    transform = Transformer.from_crs(5186, 4326, always_xy=True)
    rows = []
    for track in tracks:
        lon, lat = transform.transform(track.xyz[:, 0], track.xyz[:, 1])
        # Kepler GeoJSON column: native line layer, without connecting aircraft
        # or asking Kepler to infer point ordering.
        geometry = {'type': 'Feature', 'properties': {'ring': track.name},
                    'geometry': {'type': 'LineString', 'coordinates':
                                 [[float(x), float(y), float(z)] for x, y, z in zip(lon, lat, track.xyz[:, 2])]}}
        rows.append(dict(ring=track.name, geometry=json.dumps(geometry, ensure_ascii=False)))
    csv(rows, out / 'kepler_ring_paths.csv')


def report(preloads, feeds, ods, timeline, summary, out):
    configure_plot()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ring, data in timeline.groupby('ring'):
        axes[0].plot(data.time_s / 60, data.moving, label=ring)
        axes[1].plot(data.time_s / 60, data.external_queue, label=ring)
    for ax in axes:
        ax.set_xlabel('시간 (분)'); ax.grid(alpha=.2); ax.legend()
    axes[0].set_title('계속 생성하며 순환 중인 실제 기체 수'); axes[0].set_ylabel('순환 기체 (대)')
    axes[1].set_title('진입 간격 부족으로 남은 외부 수요'); axes[1].set_ylabel('대기 기체 (대)')
    fig.tight_layout(); fig.savefig(out / 'continuous_circulation.png', dpi=170); plt.close(fig)
    lines = ['# 순환 용량 확장 실험', '',
        '기존 분석의 안정 왕복 처리량과 링 자체의 순환 가능 대수를 분리했다. '
        '아래 링 실험에서는 기체가 실제 SHP 선형을 따라 계속 이동하며 한 단계 안의 모든 기체 쌍 최소거리도 검사한다.', '',
        f'- 링 실험 {summary["ring_experiments"]}개, 전체 OD 실험 {summary["full_OD_experiments"]}개.',
        '- 링별 진입 수요 120/360/720대/h를 1시간 동안 계속 생성한다. '
        '링 자체 실험에는 착륙·출항이 없으므로 진입한 기체는 종료까지 순환한다.',
        '- 단일 진입점과 출발지 비율을 반영한 여러 진입점을 비교한다. '
        '다중 진입은 링 접속점에서 시작하며 접근선·Interchange 이동은 생략한다.',
        '- 외곽 반시계/내곽 시계, 반경 50m 버퍼, 분리거리 100/250/500/1000m를 적용한다. '
        '고정된 SHP/default 속도로 움직이며 진입 후 재배치나 정지는 하지 않는다.',
        '- 기존 14대, 확대 28대, 기하학 상한으로 진입 한도를 비교한다. '
        '공간이 부족하면 새 기체가 외부에서 기다린다.', '',
        '## 미리 배치한 기체의 순환', '',
        '균등 배치와 최소 간격을 둔 무작위 간격 배치를 각각 900초 동안 시험했다. '
        '이는 처음부터 배치한 순환 가능 대수이며, 접근선으로 그 대수까지 채울 수 있다는 뜻은 아니다. '
        '상한 초과 조건도 의도적으로 시험하고 실패로 기록했다.', '',
        '| 링 | 분리거리(m) | 검증된 최대 사전배치 순환 대수 |', '|---|---:|---:|']
    for row in summary['preload_max']:
        lines.append(f'| {row["ring"]} | {row["separation_m"]} | {row["count"]} |')
    lines += ['', '## 계속 생성하여 실제 진입시킨 결과', '',
        '같은 분리거리·수요·시드·진입 방식을 두 링에 적용한 조건 중 마지막 15분의 '
        '최소 동시 순환 대수 합이 가장 큰 조건을 Kepler용으로 선정했다. '
        '두 링은 독립 시험이며 링 간 환승이나 착륙 처리량의 최적해는 아니다.', '',
        '| 링 | 생성 대수 | 실제 진입/계속 순환 | 마지막 15분 최소 순환 | 10대 이상 동시 이동 시간(s) | 최소 거리(m) |',
        '|---|---:|---:|---:|---:|---:|']
    for row in summary['selected_ring_metrics']:
        lines.append(f'| {row["ring"]} | {row["generated"]} | {row["admitted"]} | {row["tail_min_moving"]} | '
                     f'{row["seconds_at_least_10_moving"]:.0f} | {row["min_swept_m"]:.2f} |')
    lines += ['', '## 전체 OD의 과부하와 추가 순환 비교', '',
        '접근선→버티포트→출항의 기존 교차 예약 정책을 유지하며 수요 120/360대/h, '
        '링 한도 14/28대, 두 시드, 방문 링마다 추가 0/2바퀴를 비교했다. '
        '점유 대수와 실제 이동 중인 대수를 따로 기록한다. '
        '추가 바퀴는 최단시간 경로에 의도적으로 더한 순환 대기 시험으로, 최단경로라고 부르지 않는다.', '',
        f'- 전체 OD에서 외곽 최대 점유 {ods.peak_outer.max():.0f}대 / 최대 동시 이동 {ods.peak_moving_outer.max():.0f}대.',
        f'- 전체 OD에서 내곽 최대 점유 {ods.peak_inner.max():.0f}대 / 최대 동시 이동 {ods.peak_moving_inner.max():.0f}대.',
        '- Kepler OD 조건은 10대 이상 이동 지속시간, 동시 이동 최대 대수 합, 착륙+이륙 처리량 순으로 선정했다. '
        '이는 시험한 조건의 비교 결과이며 전역 최적 운항 해는 아니다.', '',
        '## Kepler CSV', '',
        '- `kepler_ring_flow.csv`: 계속 생성한 링 순환 실험에서 선정한 시간별 점 궤적(5초 간격과 진입 시점).',
        '- `kepler_OD_flow.csv`: 전체 OD에서 선정한 실제 이동·대기 궤적(10초 간격).',
        '- `kepler_ring_paths.csv`: GeoJSON geometry 열로 외곽·내곽의 원본 경로를 표시한다.',
        '- 좌표는 EPSG:5186에서 EPSG:4326으로 변환했다. longitude/latitude는 도, altitude_m는 원본 Z(m)다.',
        '- timestamp는 UTC ISO 8601이며 2026-10-20 00:00:00 UTC를 시뮬레이션 0초로 정한 가상 시각이다. timestamp_ms도 함께 제공한다.',
        '- Kepler에서 CSV를 추가하고 Point 레이어의 Lat/Lng를 latitude/longitude로 지정한다. '
        '색상은 aircraft_id 또는 ring을 사용하고 timestamp에 시간 필터를 적용한다. '
        '시간창을 약 5~15초로 좁혀 재생하면 현재 위치, 길게 잡으면 이동 흔적을 볼 수 있다.',
        '- 점 CSV는 시간별 포인트 데이터다. 각 기체의 점을 다른 기체와 연결한 선으로 해석하지 않는다. '
        '기체별 정렬 순서는 aircraft_id, time_s다.', '',
        '## 재실행', '', '```powershell',
        'python -B -X utf8 Network/코드/circulation_capacity.py --workers 4', '```', '']
    lines += ['Kepler 처리량 우선 CSV와 Trip 애니메이션 추가 생성:', '', '```powershell',
              'python -B -X utf8 Network/코드/export_kepler_flows.py', '```', '']
    (out / '순환분석보고서.md').write_text('\n'.join(lines), encoding='utf-8')


if __name__ == '__main__':
    main()
