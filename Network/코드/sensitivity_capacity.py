"""Check selected-condition sensitivity to missing operational inputs and time step."""
import sys
sys.dont_write_bytecode = True
import copy
import json
from network_capacity import Network, OUTPUT, demand, simulate, csv


def main():
    cfg = json.loads((OUTPUT / 'effective_config.json').read_text(encoding='utf-8'))
    summary = json.loads((OUTPUT / 'summary.json').read_text(encoding='utf-8'))
    selected = summary['selected_aggregate']
    variants = [dict(name='base'), dict(name='time_step_0.5s', step_s=.5),
                dict(name='turnaround_0s', turnaround_s=0),
                dict(name='turnaround_300s', turnaround_s=300),
                dict(name='other_speeds_minus20pct', speed_factor=.8)]
    rows = []
    for variant in variants:
        options = copy.deepcopy(cfg)
        for field in ('step_s', 'turnaround_s'):
            if field in variant:
                options[field] = variant[field]
        if 'speed_factor' in variant:
            for layer in options['speed_defaults_kts']:
                if layer != 'Short_Final':
                    options['speed_defaults_kts'][layer] *= variant['speed_factor']
        net = Network(options)
        for rep in range(options['replicates']):
            seed = options['seed'] + rep
            requests = demand(net, options, selected['arrival_demand_per_h'], seed)
            metrics, *_ = simulate(net, options, requests,
                (int(selected['outer_capacity']), int(selected['inner_capacity'])),
                selected['separation_m'], selected['scenario'], seed)
            metrics.update(variant=variant['name'], step_s=options['step_s'], turnaround_s=options['turnaround_s'])
            rows.append(metrics)
            print(f'{variant["name"]}: seed={seed}, movements={metrics["movements_per_h"]:.1f}, '
                  f'p95={metrics["p95_delay_s"]:.1f}, stable={metrics["stable"]}', flush=True)
    csv(rows, OUTPUT / 'sensitivity_results.csv')
    lines = ['\n## 가정과 시간 간격 민감도', '',
             '선정된 용량·분리거리·수요를 고정하고 각 조건을 같은 두 시드로 다시 계산했다. '
             '각 변형의 최대 용량을 다시 탐색한 것은 아니다.', '',
             '| 변형 | 평균 처리량(회/h) | 안정 반복/전체 | 최악 95백분위 지연(s) |',
             '|---|---:|---:|---:|']
    for variant in variants:
        group = [r for r in rows if r['variant'] == variant['name']]
        lines.append(f'| {variant["name"]} | {sum(r["movements_per_h"] for r in group)/len(group):.1f} | '
                     f'{sum(r["stable"] for r in group)}/{len(group)} | {max(r["p95_delay_s"] for r in group):.1f} |')
    lines += ['', 'Short_Final 속도는 모든 변형에서 30kt로 유지했다. step 0.5s는 시간 이산화의 민감도 확인이며 '
              '비행역학의 유효성 검증은 아니다. 속도/점유시간이 달라지면 기준 분석의 용량값을 그대로 적용하지 않는다.',
              '', '재실행: `python -B -X utf8 Network/코드/sensitivity_capacity.py`', '']
    report = OUTPUT / '분석보고서.md'
    text = report.read_text(encoding='utf-8').split('\n## 가정과 시간 간격 민감도')[0]
    report.write_text(text + '\n'.join(lines), encoding='utf-8')


if __name__ == '__main__':
    main()
