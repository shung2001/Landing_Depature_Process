"""내곽/외곽 점을 PointZ로, Interchange 곡선을 PolyLineZ로 함께 갱신.

필요 패키지: pip install pyshp pyproj
실행: python 코드/polygon/connection.py
고도 변경: python 코드/polygon/connection.py --inner-z 457.2 --outer-z 609.6
검토만 수행: 위 명령에 --dry-run 추가.
최종정리 Short_Final 갱신: python 코드/polygon/connection.py --short-final
Short_Final 모드는 FAF/버티포트의 XYZ와 최종OLS 대응 관계로 기존 파일을
갱신한다. 누락 연결을 보완하며 별도 출력/백업 파일을 생성하지 않는다.

고도는 같은 수직 기준의 미터 값이어야 한다. DEM/AGL/해발고도 변환은 하지 않는다.
기본 고도는 사용자가 지정한 내곽 457.2m, 외곽 609.6m이다.
--inner-z/--outer-z 또는 --inner-z-field/--outer-z-field로 변경할 수 있다.
네 자료 모두 기존 XY/속성/좌표계와 파일 접근 권한을 보존하고, SHP/SHX를 백업 후 갱신한다.
중간 꼭짓점 Z는 평면 곡선의 누적 길이 비율로 선형 보간한다.
"""

import argparse
import io
from contextlib import ExitStack, contextmanager
from datetime import datetime
import math
from pathlib import Path
import shutil
import tempfile

import shapefile
from pyproj import CRS, Transformer


ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / '자료/안전성/송파구'
LINE_NAMES = ('Interchange_In', 'Interchange_out')
POINT_NAMES = ('내곽_Interchange', '외곽_IAF_Interchange')


def finite(value, label):
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f'{label}: 유효한 숫자가 필요합니다: {value!r}') from error
    if not math.isfinite(number):
        raise ValueError(f'{label}: 유한한 숫자가 필요합니다.')
    return number


def read_crs(path):
    return CRS.from_wkt(path.with_suffix('.prj').read_text(encoding='utf-8-sig'))


def encoding_for(path):
    cpg = path.with_suffix('.cpg')
    return cpg.read_text(encoding='utf-8-sig').strip() if cpg.exists() else 'utf-8'


@contextmanager
def open_reader(path, geometry_only=False):
    """필수 파일을 직접 열어 PyShp의 암묵적인 파일 탐색을 피한다."""
    path = Path(path).resolve()
    try:
        with ExitStack() as stack:
            files = {}
            for suffix in ('.shp', '.shx') if geometry_only else ('.shp', '.shx', '.dbf'):
                source = path.with_suffix(suffix)
                try:
                    files[suffix[1:]] = stack.enter_context(source.open('rb'))
                except OSError as error:
                    raise OSError(f'파일을 열 수 없습니다: {source} ({error})') from error
            reader = stack.enter_context(shapefile.Reader(
                **files, encoding=encoding_for(path)))
            yield reader
    except shapefile.ShapefileException as error:
        raise shapefile.ShapefileException(f'{path}: {error}') from error


def load_points(path, target_crs, override=None, field=None, allow_multipoint=False):
    """점 순서/id에 의존하지 않고 매칭할 수 있도록 동일한 평면 좌표계로 변환."""
    transform = Transformer.from_crs(read_crs(path), target_crs, always_xy=True)
    result = []
    with open_reader(path) as reader:
        accepted = (shapefile.POINT, shapefile.POINTZ)
        if allow_multipoint:
            accepted += (shapefile.MULTIPOINT, shapefile.MULTIPOINTZ)
        if reader.shapeType not in accepted:
            raise ValueError(f'{path.name}: Point 또는 PointZ 자료가 필요합니다.')
        if field and field not in [f[0] for f in reader.fields[1:]]:
            raise ValueError(f'{path.name}: 고도 필드 {field!r}가 없습니다.')
        for index, item in enumerate(reader.iterShapeRecords()):
            shape = item.shape
            label = f'{path.name} 점 {index}'
            if len(shape.points) != 1:
                raise ValueError(f'{label}: 단일 점이 필요합니다.')
            x, y = transform.transform(*shape.points[0])
            x, y = finite(x, label), finite(y, label)
            if override is not None:
                z = override
            elif field:
                z = item.record[field]
            elif shape.shapeType in (shapefile.POINTZ, shapefile.MULTIPOINTZ):
                z = shape.z[0]
            else:
                raise ValueError(f'{label}: Z가 없습니다. 고도 값 또는 고도 필드를 지정하세요.')
            result.append((x, y, finite(z, label)))
    if not result:
        raise ValueError(f'{path.name}: 점 자료가 비어 있습니다.')
    return result


def match_endpoint(xy, points, tolerance, label):
    matches = [p for p in points if math.dist(xy, p[:2]) <= tolerance]
    if len(matches) != 1:
        raise ValueError(f'{label}: 허용 거리 {tolerance}m 내 대응점 {len(matches)}개 (1개 필요).')
    return matches[0][2]


def lift_line(shape, inner, outer, tolerance, label):
    if shape.shapeType not in (shapefile.POLYLINE, shapefile.POLYLINEZ):
        raise ValueError(f'{label}: PolyLine 또는 PolyLineZ가 필요합니다.')
    if list(shape.parts) != [0] or len(shape.points) < 2:
        raise ValueError(f'{label}: 꼭짓점 2개 이상의 단일 곡선이 필요합니다.')
    xy = [tuple(finite(v, label) for v in p) for p in shape.points]
    # In/Out 모두 실제 선의 저장 방향을 보존한다. 반대 방향으로 저장된 선도 처리.
    start_outer = any(math.dist(xy[0], p[:2]) <= tolerance for p in outer)
    start_group, end_group = (outer, inner) if start_outer else (inner, outer)
    z0 = match_endpoint(xy[0], start_group, tolerance, label + ' 시작점')
    z1 = match_endpoint(xy[-1], end_group, tolerance, label + ' 끝점')
    lengths = [0.0]
    for a, b in zip(xy, xy[1:]):
        lengths.append(lengths[-1] + math.dist(a, b))
    if lengths[-1] <= 0:
        raise ValueError(f'{label}: 곡선의 평면 길이가 0입니다.')
    xyz = [(x, y, z0 + (z1 - z0) * length / lengths[-1])
           for (x, y), length in zip(xy, lengths)]
    xyz[0], xyz[-1] = (*xy[0], z0), (*xy[-1], z1)
    return xyz


def stage_lines(path, stage_path, inner, outer, tolerance):
    """DBF/PRJ/CPG는 그대로 두고 교체할 SHP/SHX만 생성하고 검증한다."""
    with open_reader(path) as reader:
        lines = [lift_line(s, inner, outer, tolerance, f'{path.name} 선 {i}')
                 for i, s in enumerate(reader.iterShapes())]
        if not lines or len(lines) != len(reader):
            raise ValueError(f'{path.name}: 곡선/속성 레코드 개수가 유효하지 않습니다.')
    with shapefile.Writer(shp=str(stage_path.with_suffix('.shp')),
                          shx=str(stage_path.with_suffix('.shx')),
                          shapeType=shapefile.POLYLINEZ) as writer:
        for xyz in lines:
            writer.linez([xyz])
    # 원본 DBF를 복사하지 않고 기하 파일만 다시 열어 저장된 XYZ를 확인한다.
    with open_reader(stage_path, geometry_only=True) as saved:
        shapes = list(saved.iterShapes())
        if saved.shapeType != shapefile.POLYLINEZ or len(shapes) != len(lines):
            raise RuntimeError(f'{path.name}: PolyLineZ 저장 검증 실패.')
        for shape, expected in zip(shapes, lines):
            actual = [(*p, z) for p, z in zip(shape.points, shape.z)]
            if actual != expected:
                raise RuntimeError(f'{path.name}: XYZ 저장 검증 실패.')
    zs = [p[2] for line in lines for p in line]
    print(f'{path.name}: {len(lines)}개 곡선, {sum(map(len, lines))}개 꼭짓점, '
          f'Z {min(zs):g} ~ {max(zs):g} m')


def stage_points(path, stage_path, points):
    """원래 좌표계의 XY와 M을 보존하고 연결선과 같은 Z를 저장한다."""
    with open_reader(path) as reader:
        shapes = list(reader.iterShapes())
        if len(shapes) != len(points) or len(shapes) != len(reader):
            raise ValueError(f'{path.name}: 점/속성 개수가 일치하지 않습니다.')
    expected = [(*shape.points[0], point[2]) for shape, point in zip(shapes, points)]
    with shapefile.Writer(shp=str(stage_path.with_suffix('.shp')),
                          shx=str(stage_path.with_suffix('.shx')),
                          shapeType=shapefile.POINTZ) as writer:
        for shape, xyz in zip(shapes, expected):
            measures = getattr(shape, 'm', None)
            writer.pointz(*xyz, m=measures[0] if measures else None)
    with open_reader(stage_path, geometry_only=True) as saved:
        actual = [(*s.points[0], s.z[0]) for s in saved.iterShapes()]
        if saved.shapeType != shapefile.POINTZ or actual != expected:
            raise RuntimeError(f'{path.name}: PointZ 저장 검증 실패.')
    zs = [point[2] for point in points]
    print(f'{path.name}: {len(points)}개 점, Z {min(zs):g} ~ {max(zs):g} m')


def update(directory, inner_z=None, outer_z=None, inner_z_field=None,
           outer_z_field=None, tolerance=0.01, dry_run=False):
    directory = Path(directory).resolve()
    tolerance = finite(tolerance, 'tolerance')
    if tolerance <= 0:
        raise ValueError('허용 거리는 0보다 커야 합니다.')
    paths = [directory / f'{name}.shp' for name in LINE_NAMES]
    target_crs = read_crs(paths[0])
    if not target_crs.is_projected or any(
            axis.unit_conversion_factor != 1 for axis in target_crs.axis_info[:2]):
        raise ValueError('연결선은 미터 단위의 투영 좌표계여야 합니다.')
    if read_crs(paths[1]) != target_crs:
        raise ValueError('두 연결선의 좌표계가 다릅니다.')
    if inner_z is None and inner_z_field is None:
        inner_z = 457.2
    if outer_z is None and outer_z_field is None:
        outer_z = 609.6
    point_paths = [directory / f'{name}.shp' for name in POINT_NAMES]
    inner = load_points(point_paths[0], target_crs, inner_z, inner_z_field)
    outer = load_points(point_paths[1], target_crs, outer_z, outer_z_field)
    all_paths = paths + point_paths
    with tempfile.TemporaryDirectory(prefix='.connection_', dir=directory) as temporary:
        stage_dir = Path(temporary)
        for path in paths:
            stage_lines(path, stage_dir / path.name, inner, outer, tolerance)
        for path, points in zip(point_paths, (inner, outer)):
            stage_points(path, stage_dir / path.name, points)
        if dry_run:
            print('검토 완료: 원본 파일은 갱신하지 않았습니다.')
            return None
        backup = directory / 'connection_backup' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        backup.mkdir(parents=True)
        for path in all_paths:
            for suffix in ('.shp', '.shx', '.dbf', '.prj', '.cpg'):
                source = path.with_suffix(suffix)
                if source.exists():
                    shutil.copy2(source, backup / source.name)
        replaced = []
        try:
            for path in all_paths:
                for suffix in ('.shp', '.shx'):
                    destination = path.with_suffix(suffix)
                    # 임시 디렉터리의 제한된 Windows ACL이 원본에 전파되지 않도록
                    # 기존 파일을 유지하고 내용만 복사한다. 부분 쓰기 실패도 복원한다.
                    replaced.append(destination)
                    shutil.copyfile(stage_dir / destination.name, destination)
        except OSError:
            for destination in reversed(replaced):
                shutil.copy2(backup / destination.name, destination)
            raise
        print(f'갱신 완료. 원본 백업: {backup}')
        return backup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=None)
    parser.add_argument('--short-final', action='store_true', help='최종정리의 Short_Final을 3D 연결선으로 제자리 갱신')
    for group in ('inner', 'outer'):
        selection = parser.add_mutually_exclusive_group()
        selection.add_argument(f'--{group}-z', type=float, help='모든 점에 적용할 고도(m)')
        selection.add_argument(f'--{group}-z-field', help='고도로 사용할 점 속성 필드')
    parser.add_argument('--tolerance', type=float, default=0.01, help='끝점 매칭 허용 거리(m)')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    try:
        options = vars(args).copy()
        short_final = options.pop('short_final')
        options['directory'] = options['directory'] or (
            ROOT / '자료/최종정리' if short_final else DIRECTORY)
        if short_final:
            if any(options[name] is not None for name in
                   ('inner_z', 'outer_z', 'inner_z_field', 'outer_z_field')):
                raise ValueError('Short_Final은 원본 점의 Z를 사용합니다. 내곽/외곽 고도 옵션을 제거하세요.')
            update_short_final(options['directory'], options['tolerance'], options['dry_run'])
        else:
            update(**options)
    except (ValueError, OSError, shapefile.ShapefileException) as error:
        parser.exit(1, f'오류: {error}\n')


def update_short_final(directory, tolerance=0.01, dry_run=False):
    """원본 점 XYZ를 정확히 잇고, 메모리 검증/복원으로 기존 파일만 갱신."""
    directory = Path(directory).resolve()
    path = directory / 'Short_Final.shp'
    tolerance = finite(tolerance, 'tolerance')
    if tolerance <= 0:
        raise ValueError('허용 거리는 0보다 커야 합니다.')
    crs = read_crs(path)
    if not crs.is_projected or any(a.unit_conversion_factor != 1 for a in crs.axis_info[:2]):
        raise ValueError('Short_Final은 미터 단위 투영 좌표계여야 합니다.')
    faf = load_points(directory / 'FAF.shp', crs, allow_multipoint=True)
    ports = load_points(directory / '각_건물의_최종_버티포트.shp', crs)

    def unique_match(point, candidates, label, xyz=False):
        matches = [i for i, candidate in enumerate(candidates)
                   if math.dist(point if xyz else point[:2],
                                candidate if xyz else candidate[:2]) <= tolerance]
        if len(matches) != 1:
            raise ValueError(f'{label}: 대응점 {len(matches)}개 (1개 필요).')
        return matches[0]

    # 가장 가까운 건물로 추정하지 않고 각 FAF를 생성한 OLS의 저점으로 대응.
    ols_path = directory / '최종OLS.shp'
    transform = Transformer.from_crs(read_crs(ols_path), crs, always_xy=True)
    targets = {}
    with open_reader(ols_path) as reader:
        for i, shape in enumerate(reader.iterShapes()):
            if shape.shapeType != shapefile.POLYGONZ or list(shape.parts) != [0]:
                raise ValueError(f'OLS {i}: 단일 PolygonZ가 필요합니다.')
            vertices = [(*transform.transform(*xy), finite(z, f'OLS {i} Z'))
                        for xy, z in zip(shape.points, shape.z)]
            if vertices[-1] == vertices[0]:
                vertices.pop()
            if len(vertices) != 3:
                raise ValueError(f'OLS {i}: 삼각형이 필요합니다.')
            vertices.sort(key=lambda p: p[2])
            low, upper1, upper2 = vertices
            midpoint = tuple((a + b) / 2 for a, b in zip(upper1, upper2))
            fi = unique_match(midpoint, faf, f'OLS {i} FAF', xyz=True)
            vi = unique_match(low, ports, f'OLS {i} 버티포트', xyz=True)
            if fi in targets:
                raise ValueError(f'FAF {fi}: 중복 OLS 대응입니다.')
            targets[fi] = vi
    if set(targets) != set(range(len(faf))):
        raise ValueError('모든 FAF에 대응하는 OLS가 필요합니다.')

    with open_reader(path) as reader:
        fields = reader.fields[1:]
        records = [list(record) for record in reader.records()]
        shapes = list(reader.iterShapes())
    if len(shapes) != len(records):
        raise ValueError('Short_Final의 기하/속성 개수가 다릅니다.')
    lines, used, corrected = [], set(), 0
    for i, shape in enumerate(shapes):
        if shape.shapeType not in (shapefile.POLYLINE, shapefile.POLYLINEZ) or list(shape.parts) != [0] or len(shape.points) != 2:
            raise ValueError(f'Short_Final {i}: 두 꼭짓점의 단일 연결선이 필요합니다.')
        xy = shape.points
        matches = [(end, fi) for end in (0, 1) for fi, point in enumerate(faf)
                   if math.dist(xy[end], point[:2]) <= tolerance]
        if len(matches) != 1:
            raise ValueError(f'Short_Final {i}: FAF 끝점이 유일하지 않습니다.')
        end, fi = matches[0]
        if fi in used:
            raise ValueError(f'Short_Final {i}: 중복 FAF 연결입니다.')
        old_vi = unique_match(xy[1 - end], ports, f'Short_Final {i} 버티포트')
        vi = targets[fi]
        corrected += old_vi != vi
        endpoints = [None, None]
        endpoints[end], endpoints[1 - end] = faf[fi], ports[vi]
        lines.append(endpoints)
        used.add(fi)
    id_column = next((i for i, field in enumerate(fields) if field[0].lower() == 'id'), None)
    ids = {record[id_column] for record in records} if id_column is not None else set()
    missing = sorted(set(targets) - used)
    for fi in missing:
        lines.append([faf[fi], ports[targets[fi]]])
        record = [None] * len(fields)
        if id_column is not None:
            new_id = 1
            while new_id in ids:
                new_id += 1
            record[id_column] = new_id
            ids.add(new_id)
        records.append(record)

    buffers = {name: io.BytesIO() for name in ('shp', 'shx', 'dbf')}
    with shapefile.Writer(**buffers, shapeType=shapefile.POLYLINEZ,
                          encoding=encoding_for(path)) as writer:
        writer.fields = fields
        for line, record in zip(lines, records):
            writer.linez([line])
            writer.record(*record)
    content = {name: buffer.getvalue() for name, buffer in buffers.items()}

    def validate(reader):
        actual = [[(*xy, z) for xy, z in zip(s.points, s.z)] for s in reader.iterShapes()]
        if reader.shapeType != shapefile.POLYLINEZ or actual != lines:
            raise ValueError('Short_Final XYZ 저장 검증 실패.')
        if [list(record) for record in reader.records()] != records:
            raise ValueError('Short_Final 속성 저장 검증 실패.')

    with shapefile.Reader(**{name: io.BytesIO(data) for name, data in content.items()},
                          encoding=encoding_for(path)) as reader:
        validate(reader)
    print(f'Short_Final: {len(shapes)} → {len(lines)}개 PolyLineZ, '
          f'누락 보완 {len(missing)}개, 버티포트 대응 보정 {corrected}개')
    if dry_run:
        print('검토 완료: 원본 파일은 갱신하지 않았습니다.')
        return
    destinations = {name: path.with_suffix('.' + name) for name in content}
    originals = {name: dest.read_bytes() for name, dest in destinations.items()}
    changed = []
    try:
        for name, dest in destinations.items():
            changed.append(name)
            with dest.open('r+b') as stream:
                stream.write(content[name])
                stream.truncate()
        with open_reader(path) as reader:
            validate(reader)
    except Exception:
        for name in reversed(changed):
            with destinations[name].open('r+b') as stream:
                stream.write(originals[name])
                stream.truncate()
        raise
    print('기존 Short_Final.shp/.shx/.dbf 갱신 완료. 모든 연결 끝점 XYZ 일치 확인.')


if __name__ == '__main__':
    main()
