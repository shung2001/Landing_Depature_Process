"""IAF_FAF/IAF_FAF_Out의 실제 끝점에 대응하는 점 Z를 연결해 원본을 갱신한다.

실행: python 코드/polygon/connection_1.py
검토: python 코드/polygon/connection_1.py --dry-run
필요 패키지: pyshp, pyproj. 같은 폴더의 connection.py도 필요하다.

현재 각 선 파일의 8개 곡선은 내곽_IAF(457.2m), 8개는 외곽_IAF_Interchange
(609.6m)에 연결된다. FAF는 각 점의 실제 Z를 사용한다. ID/행 순서는 사용하지 않는다.
XY 경로, 점 파일, DBF 속성, 좌표계는 유지한다. 고도는 곡선의 누적 평면 길이에
비례해 보간하며 수직 기준/단위는 변환하지 않는다. 모든 끝점을 검증한 후
원본을 백업하고 SHP/SHX 내용만 갱신하여 Windows 접근 권한을 보존한다.
"""

import argparse
from datetime import datetime
import math
from pathlib import Path
import shutil
import tempfile

import shapefile
from pyproj import Transformer

from connection import DIRECTORY, finite, open_reader, read_crs, stage_lines


LINE_NAMES = ('IAF_FAF', 'IAF_FAF_Out')
FAF_RELATIVE = Path('Separation_확보/최종_OLS/최종/FAF.shp')


def load_z_points(path, target_crs):
    """PointZ와 하나의 점을 담은 MultiPointZ를 읽는다. 실제 Z만 사용한다."""
    transform = Transformer.from_crs(read_crs(path), target_crs, always_xy=True)
    points = []
    with open_reader(path) as reader:
        if reader.shapeType not in (shapefile.POINTZ, shapefile.MULTIPOINTZ):
            raise ValueError(f'{path}: PointZ 또는 MultiPointZ가 필요합니다.')
        for index, shape in enumerate(reader.iterShapes()):
            label = f'{path.name} 점 {index + 1}'
            if len(shape.points) != 1 or len(shape.z) != 1:
                raise ValueError(f'{label}: 피처당 점이 정확히 하나여야 합니다.')
            x, y = transform.transform(*shape.points[0])
            points.append(tuple(finite(v, label) for v in (x, y, shape.z[0])))
    if not points:
        raise ValueError(f'{path}: 점 자료가 비어 있습니다.')
    return points


def report_matches(path, groups, tolerance):
    """FAF 반대편 끝점의 출처를 검증하고 집계한다."""
    counts = {name: 0 for name in groups if name != 'FAF'}
    with open_reader(path) as reader:
        for index, shape in enumerate(reader.iterShapes()):
            if len(shape.points) < 2:
                raise ValueError(f'{path.name} 선 {index + 1}: 끝점이 부족합니다.')
            faf_ends = [any(math.dist(xy, p[:2]) <= tolerance for p in groups['FAF'])
                        for xy in (shape.points[0], shape.points[-1])]
            if sum(faf_ends) != 1:
                raise ValueError(f'{path.name} 선 {index + 1}: FAF 쪽 끝점이 유일하지 않습니다.')
            xy = shape.points[-1] if faf_ends[0] else shape.points[0]
            matches = [name for name, points in groups.items() if name != 'FAF'
                       for p in points if math.dist(xy, p[:2]) <= tolerance]
            if len(matches) != 1:
                nearest = min(math.dist(xy, p[:2]) for name, points in groups.items()
                              if name != 'FAF' for p in points)
                raise ValueError(f'{path.name} 선 {index + 1}: IAF 대응점 {len(matches)}개, '
                                 f'최근접 거리 {nearest:g}m. XY를 자동 이동하지 않습니다.')
            counts[matches[0]] += 1
    print(f'{path.name} 끝점 매칭: {counts}')


def update(directory=DIRECTORY, tolerance=0.01, dry_run=False):
    directory = Path(directory).resolve()
    tolerance = finite(tolerance, 'tolerance')
    if tolerance <= 0:
        raise ValueError('허용 거리는 0보다 커야 합니다.')
    paths = [directory / f'{name}.shp' for name in LINE_NAMES]
    crs = read_crs(paths[0])
    if not crs.is_projected or any(a.unit_conversion_factor != 1 for a in crs.axis_info[:2]):
        raise ValueError('연결선은 미터 단위 투영 좌표계여야 합니다.')
    if read_crs(paths[1]) != crs:
        raise ValueError('두 연결선의 좌표계가 다릅니다.')
    groups = {
        'FAF': load_z_points(directory / FAF_RELATIVE, crs),
        '외곽_IAF_Interchange': load_z_points(directory / '외곽_IAF_Interchange.shp', crs),
        '내곽_IAF': load_z_points(directory / '내곽_IAF.shp', crs),
    }
    iaf = groups['외곽_IAF_Interchange'] + groups['내곽_IAF']
    with tempfile.TemporaryDirectory(prefix='.connection_1_', dir=directory) as temporary:
        stage_dir = Path(temporary)
        for path in paths:
            report_matches(path, groups, tolerance)
            stage_lines(path, stage_dir / path.name, groups['FAF'], iaf, tolerance)
        if dry_run:
            print('검토 완료: 원본은 변경하지 않았습니다.')
            return None
        backup = directory / 'connection_1_backup' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        backup.mkdir(parents=True)
        for path in paths:
            for suffix in ('.shp', '.shx', '.dbf', '.prj', '.cpg'):
                source = path.with_suffix(suffix)
                if source.exists():
                    shutil.copy2(source, backup / source.name)
        changed = []
        try:
            for path in paths:
                for suffix in ('.shp', '.shx'):
                    target = path.with_suffix(suffix)
                    changed.append(target)
                    shutil.copyfile(stage_dir / target.name, target)
        except OSError:
            for target in reversed(changed):
                shutil.copyfile(backup / target.name, target)
            raise
        print(f'갱신 완료. 원본 백업: {backup}')
        return backup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=DIRECTORY)
    parser.add_argument('--tolerance', type=float, default=0.01, help='끝점 매칭 허용 거리(m)')
    parser.add_argument('--dry-run', action='store_true')
    try:
        update(**vars(parser.parse_args()))
    except (ValueError, OSError, shapefile.ShapefileException) as error:
        parser.exit(1, f'오류: {error}\n')


if __name__ == '__main__':
    main()
