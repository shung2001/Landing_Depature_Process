"""IAF_FAF/IAF_FAF_Out의 실제 끝점에 대응하는 점 Z를 연결해 원본을 갱신한다.

실행: python 코드/polygon/connection_1.py
검토: python 코드/polygon/connection_1.py --dry-run
특정 행 검토: python 코드/polygon/connection_1.py --out-rows 18 19 --dry-run
행 선택 옵션이 없으면 두 파일의 활성 레코드 전체를 처리한다. 행 번호는 1부터,
QGIS FID는 0부터 시작한다. 속성 id와는 별개이다. 삭제 표시 레코드는 보존만 한다.
필요 패키지: pyshp, pyproj. 같은 폴더의 connection.py도 필요하다.

연결 대상은 내곽_IAF(457.2m), 외곽_IAF_Interchange(609.6m)이다.
실제 곡선 수는 편집 상태에 따라 달라진다. FAF는 각 점의 실제 Z를 사용한다.
끝점 대응에는 ID/행 순서 대신 XY 위치를 사용한다.
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

# connection.py의 공통 기능을 재사용한다.
# DIRECTORY: 기본 자료 폴더 / finite: NaN·무한대 등을 거르는 숫자 검사
# open_reader: SHP·SHX·DBF를 명시적으로 열고 자동으로 닫는 판독기
# read_crs: PRJ에서 좌표계 읽기 / lift_line: 끝점 Z 매칭과 누적 XY 거리 비례 보간
from connection import DIRECTORY, finite, open_reader, read_crs, lift_line


# 갱신할 선 파일 이름(확장자 제외). In/Out 구분은 저장된 선의 방향을 뜻하지 않는다.
LINE_NAMES = ('IAF_FAF', 'IAF_FAF_Out')
# 기본 자료 폴더를 기준으로 한 FAF 점 파일의 상대 경로.
FAF_RELATIVE = Path('Separation_확보/최종_OLS/최종/FAF.shp')


def load_z_points(path, target_crs):
    """점 파일의 XY를 선의 좌표계로 맞추고 (x, y, z) 튜플 목록을 반환한다.

    path: PointZ 또는 점 하나짜리 MultiPointZ SHP 경로.
    target_crs: 연결선의 좌표계. XY만 변환하며 원본 Z(m)는 그대로 사용한다.
    2D 점, 여러 점이 든 피처, 비유한 좌표, 빈 자료는 오류로 처리한다.
    """
    # 축 순서가 위도·경도로 바뀌지 않도록 항상 X, Y 순서로 변환한다.
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
            # Z는 속성값으로 추정하지 않고 기하에 저장된 실제 높이를 읽는다.
            points.append(tuple(finite(v, label) for v in (x, y, shape.z[0])))
    if not points:
        raise ValueError(f'{path}: 점 자료가 비어 있습니다.')
    return points


def report_matches(path, groups, tolerance, rows=None):
    """처리할 선의 끝점 대응을 검증하고 0부터 시작하는 레코드 인덱스를 반환한다.

    groups: 'FAF' 및 IAF 그룹 이름을 키로 하는 XYZ 목록 사전.
    tolerance: 끝점과 대응점 사이에 허용할 평면 거리(m).
    rows: 1부터 시작하는 처리 행 목록. None이면 모든 활성 행을 검사한다.
    출력에는 IAF 그룹별 매칭 개수를 표시하며 원본 파일을 수정하지 않는다.
    """
    counts = {name: 0 for name in groups if name != 'FAF'}
    selected = []
    with open_reader(path) as reader:
        if rows is not None:
            if not rows or any(row < 1 or row > len(reader) for row in rows):
                raise ValueError(f'{path.name}: 행 번호는 1~{len(reader)} 범위여야 합니다.')
        for index, shape in enumerate(reader.iterShapes()):
            # 내부 인덱스/FID는 0부터, 사용자가 입력하는 행 번호는 1부터 시작한다.
            if rows is not None and index + 1 not in rows:
                continue
            record = reader.record(index)
            # DBF 삭제 표시는 SHP 기하를 물리적으로 지우지 않는다.
            # 전체 처리에서는 건너뛰고, 사용자가 직접 지정한 삭제 행은 오류로 알린다.
            if record is None:
                if rows is not None:
                    raise ValueError(f'{path.name} 행 {index + 1}: DBF에서 삭제 표시된 레코드입니다.')
                continue
            label = f'{path.name} 행 {index + 1} (FID={index}, id={record.as_dict().get("id", "없음")})'
            if len(shape.points) < 2:
                raise ValueError(f'{label}: 끝점이 부족합니다.')
            faf_ends = [any(math.dist(xy, p[:2]) <= tolerance for p in groups['FAF'])
                        for xy in (shape.points[0], shape.points[-1])]
            if sum(faf_ends) != 1:
                # FAF에 닿는 끝점이 0개/2개이면 선의 어느 쪽을 FAF로 볼지 결정할 수 없다.
                # 이 단계는 FAF에 닿는 '끝점 수'를 검사하며 개별 대응점의 유일성은
                # 뒤의 lift_line에서 추가로 검증한다.
                distances = [min(math.dist(xy, p[:2]) for p in groups['FAF'])
                             for xy in (shape.points[0], shape.points[-1])]
                reason = '두 끝점 모두 FAF와 불일치' if not any(faf_ends) else '두 끝점 모두 FAF와 일치'
                raise ValueError(f'{label}: {reason}. 허용 거리 {tolerance:g}m, '
                                 f'시작/끝점의 최근접 FAF 거리 {distances[0]:.6f}/{distances[1]:.6f}m. '
                                 '행 미지정 시 모든 활성 행을 검사합니다. XY를 자동 이동하지 않습니다.')
            # FAF가 아닌 반대편 끝점이 내곽/외곽 IAF 중 어느 점과 연결되는지 찾는다.
            xy = shape.points[-1] if faf_ends[0] else shape.points[0]
            matches = [name for name, points in groups.items() if name != 'FAF'
                       for p in points if math.dist(xy, p[:2]) <= tolerance]
            if len(matches) != 1:
                # 일치점이 없거나 여러 개면 임의 선택·자동 스냅 없이 중단한다.
                nearest = min(math.dist(xy, p[:2]) for name, points in groups.items()
                              if name != 'FAF' for p in points)
                raise ValueError(f'{label}: IAF 대응점 {len(matches)}개, '
                                 f'최근접 거리 {nearest:g}m. XY를 자동 이동하지 않습니다.')
            counts[matches[0]] += 1
            selected.append(index)
    print(f'{path.name} 끝점 매칭: {counts}')
    return selected


def stage_selected_lines(path, destination, faf, iaf, tolerance, selected):
    """선택 행의 Z를 보간해 임시 SHP/SHX를 만들고 다시 읽어 검증한다.

    path: 원본 선 파일 / destination: 임시 선 파일 경로.
    faf, iaf: 끝점 Z를 가져올 XYZ 목록 / tolerance: XY 매칭 허용 거리(m).
    selected: report_matches가 반환한 0 기반 인덱스 목록.
    원본 DBF를 재작성하지 않으며 삭제/미선택 행도 같은 순서로 보존한다.
    """
    # 반복문에서 선택 여부를 빠르게 확인하고 중복 선택을 한 번으로 처리한다.
    selected = set(selected)
    with open_reader(path) as reader:
        shapes = list(reader.iterShapes())
        if reader.shapeType != shapefile.POLYLINEZ or len(shapes) != reader.numRecords:
            raise ValueError(f'{path.name}: SHP/DBF 개수가 일치하는 PolyLineZ가 필요합니다.')
    # 저장 직후 좌표와 비교할 보간 결과를 행 인덱스별로 보관한다.
    expected = {}
    # DBF를 새로 쓰지 않고 기하(SHP)와 위치 인덱스(SHX)만 생성한다.
    with shapefile.Writer(shp=str(destination.with_suffix('.shp')),
                          shx=str(destination.with_suffix('.shx')),
                          shapeType=shapefile.POLYLINEZ) as writer:
        for index, shape in enumerate(shapes):
            if index in selected:
                # 끝점의 실제 Z를 고정하고 중간 꼭짓점은 누적 평면 거리 비율로 보간한다.
                xyz = lift_line(shape, faf, iaf, tolerance, f'{path.name} 행 {index + 1}')
                writer.linez([xyz])
                expected[index] = xyz
            else:
                # 삭제 행도 없애지 않아야 기존 DBF와 SHP의 레코드 대응이 유지된다.
                writer.shape(shape)
    # 임시 파일에는 DBF가 없으므로 기하만 열어 레코드 수와 저장 좌표를 확인한다.
    with open_reader(destination, geometry_only=True) as saved:
        actual = list(saved.iterShapes())
        if len(actual) != len(shapes):
            raise RuntimeError(f'{path.name}: 저장 후 레코드 수 불일치.')
        for index, shape in enumerate(actual):
            if index in expected:
                if [(*xy, z) for xy, z in zip(shape.points, shape.z)] != expected[index]:
                    raise RuntimeError(f'{path.name} 행 {index + 1}: XYZ 저장 검증 실패.')
            else:
                before = shapes[index]
                # 미선택 행은 XY·파트 구분·Z·M(측정값)이 모두 같은지 확인한다.
                if (list(shape.points), list(shape.parts), list(shape.z), list(shape.m)) != (
                        list(before.points), list(before.parts), list(before.z), list(before.m)):
                    raise RuntimeError(f'{path.name} 행 {index + 1}: 미선택 기하 보존 실패.')
    print(f'{path.name}: {len(selected)}개 활성 행 갱신 대상, {len(shapes) - len(selected)}개 행 보존')


def update(directory=DIRECTORY, tolerance=0.01, dry_run=False, in_rows=None, out_rows=None):
    """검증→임시 저장→백업→원본 갱신을 순서대로 수행하는 전체 처리 함수.

    directory: 입력 및 백업 폴더의 기준 경로 / tolerance: 끝점 매칭 거리(m).
    dry_run=True: 임시 결과까지 검증하고 원본 갱신 없이 None을 반환한다.
    in_rows/out_rows: 각 파일의 1 기반 행 목록. 둘 다 None이면 양쪽 전체 활성 행,
    한쪽만 지정하면 해당 파일의 지정 행만 처리한다.
    실제 갱신 시에는 백업 폴더 Path를 반환한다.
    """
    # 경로를 절대 경로로 고정하고 거리 기준이 유한한 양수인지 검사한다.
    directory = Path(directory).resolve()
    tolerance = finite(tolerance, 'tolerance')
    if tolerance <= 0:
        raise ValueError('허용 거리는 0보다 커야 합니다.')
    paths = [directory / f'{name}.shp' for name in LINE_NAMES]
    crs = read_crs(paths[0])
    # 거리 비교와 보간에 사용하는 XY가 같은 미터 단위 투영 좌표계여야 한다.
    if not crs.is_projected or any(a.unit_conversion_factor != 1 for a in crs.axis_info[:2]):
        raise ValueError('연결선은 미터 단위 투영 좌표계여야 합니다.')
    if read_crs(paths[1]) != crs:
        raise ValueError('두 연결선의 좌표계가 다릅니다.')
    # 세 점 자료의 XY를 연결선 좌표계로 통일한다. Z 기준/단위는 변환하지 않는다.
    groups = {
        'FAF': load_z_points(directory / FAF_RELATIVE, crs),
        '외곽_IAF_Interchange': load_z_points(directory / '외곽_IAF_Interchange.shp', crs),
        '내곽_IAF': load_z_points(directory / '내곽_IAF.shp', crs),
    }
    iaf = groups['외곽_IAF_Interchange'] + groups['내곽_IAF']
    # 행 선택이 하나라도 있으면 선택 옵션을 준 파일만 처리 대상으로 남긴다.
    selections = dict(zip(LINE_NAMES, (in_rows, out_rows)))
    if in_rows is not None or out_rows is not None:
        paths = [path for path in paths if selections[path.stem] is not None]
    # 원본을 건드리기 전에 모든 대상 파일의 임시 결과를 완성한다.
    # 오류가 발생해도 임시 폴더는 문맥 관리자에 의해 정리된다.
    with tempfile.TemporaryDirectory(prefix='.connection_1_', dir=directory) as temporary:
        stage_dir = Path(temporary)
        for path in paths:
            selected = report_matches(path, groups, tolerance, selections[path.stem])
            stage_selected_lines(path, stage_dir / path.name, groups['FAF'], iaf, tolerance, selected)
        if dry_run:
            # 검토 모드에서는 아래의 백업 생성·원본 쓰기를 수행하지 않는다.
            print('검토 완료: 원본은 변경하지 않았습니다.')
            return None
        # 실행별 고유 폴더에 원본 기하·속성·좌표계·인코딩 파일을 함께 복사한다.
        backup = directory / 'connection_1_backup' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        backup.mkdir(parents=True)
        for path in paths:
            for suffix in ('.shp', '.shx', '.dbf', '.prj', '.cpg'):
                source = path.with_suffix(suffix)
                if source.exists():
                    shutil.copy2(source, backup / source.name)
        # 쓰기 도중 실패할 경우를 위해 변경을 시도한 파일을 순서대로 기록한다.
        changed = []
        try:
            for path in paths:
                for suffix in ('.shp', '.shx'):
                    target = path.with_suffix(suffix)
                    changed.append(target)
                    # 파일 자체를 교체하지 않고 내용만 덮어 기존 Windows 접근 권한을 유지한다.
                    # DBF·PRJ·CPG와 원본 점 파일은 갱신하지 않는다.
                    shutil.copyfile(stage_dir / target.name, target)
        except OSError:
            # 부분적으로 쓰인 파일도 포함해 백업 내용으로 되돌린 뒤 원래 오류를 전달한다.
            for target in reversed(changed):
                shutil.copyfile(backup / target.name, target)
            raise
        print(f'갱신 완료. 원본 백업: {backup}')
        return backup


def main():
    """명령행 옵션을 해석해 update를 호출하고 예상 가능한 오류를 간결하게 출력한다."""
    # 문서 문자열을 도움말로 사용한다. --help로 실행법과 처리 범위를 확인할 수 있다.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=DIRECTORY)
    parser.add_argument('--tolerance', type=float, default=0.01, help='끝점 매칭 허용 거리(m)')
    parser.add_argument('--dry-run', action='store_true')
    # nargs='+'는 하나 이상의 행 번호를 받는다. 예: --out-rows 18 19 20
    parser.add_argument('--in-rows', type=int, nargs='+', help='IAF_FAF의 처리할 행 번호(1부터). 지정 시 해당 행만 처리')
    parser.add_argument('--out-rows', type=int, nargs='+', help='IAF_FAF_Out의 처리할 행 번호(1부터). 지정 시 해당 행만 처리')
    try:
        # 옵션 이름을 update의 매개변수 이름에 맞춰 전달한다.
        update(**vars(parser.parse_args()))
    except (ValueError, OSError, shapefile.ShapefileException) as error:
        # 입력/파일 오류는 종료 코드 1과 함께 출력한다. 예상 밖의 코드 오류는 숨기지 않는다.
        parser.exit(1, f'오류: {error}\n')


if __name__ == '__main__':
    # 직접 실행할 때만 시작한다. 다른 코드에서 import할 때에는 원본을 갱신하지 않는다.
    main()
