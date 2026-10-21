"""Validate and export fixed fleets on each independent SHP ring."""
import argparse
import json
from pathlib import Path

import numpy as np

from network_capacity import Network, OUTPUT, csv
from circulation_capacity import RingTrack, preload, kepler_csv
from export_kepler_flows import trips


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=int, nargs='+', default=[20, 60])
    parser.add_argument('--output', type=Path, default=OUTPUT / '독립링_20_60')
    args = parser.parse_args()
    if any(n < 2 for n in args.counts):
        parser.error('counts must be at least 2')
    cfg = json.loads(Path(__file__).with_name('capacity_config.json').read_text(encoding='utf-8'))
    net = Network(cfg)
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    for ri in (0, 1):
        track = RingTrack(net, ri)
        for count in args.counts:
            for separation in cfg['separations_m']:
                result = preload(track, separation, count, 'uniform', cfg['seed'])
                results.append(result)
                print(result, flush=True)
                if separation != 100 or not result['valid']:
                    continue
                rows = []
                stations = np.arange(count) * track.length / count
                for t in range(0, 901, 5):
                    for i, p in enumerate(track.at(stations + track.speed * t)):
                        rows.append(dict(time_s=t, aircraft_id=str(i), x=p[0], y=p[1],
                                         z=p[2], ring=ri, phase='circulation'))
                name = f'ring{ri}_{count}aircraft'
                frame = kepler_csv(rows, args.output / f'{name}.csv', name, cfg['seed'])
                assert frame.groupby('time_s').aircraft_id.nunique().eq(count).all()
                assert trips(frame, args.output / f'{name}.geojson') == count
    csv(results, args.output / 'validation_results.csv')
    (args.output / 'effective_config.json').write_text(
        json.dumps(dict(source_config=cfg, counts=args.counts, duration_s=900,
                        validation_step_s=1, export_step_s=5, placement='uniform',
                        experiment='independent_ring_preload'), ensure_ascii=False, indent=2),
        encoding='utf-8')


if __name__ == '__main__':
    main()
