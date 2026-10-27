"""Estimate still-air bank-limited speeds from the supplied PolyLineZ geometry.

No GIS package required: read SHP/DBF using struct, fit circles using numpy.
These are geometric estimates, not certified aircraft or procedure limits.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import struct
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / '자료' / '최종정리'
OUTPUT = ROOT / '자료' / '결과물' / 'RF_Turn_속도분석'


def read_lines(path):
    raw = path.read_bytes()
    dbf = path.with_suffix('.dbf').read_bytes()
    count, header, reclen = struct.unpack_from('<IHH', dbf, 4)
    fields = []
    pos = 32
    while dbf[pos] != 13:
        desc = dbf[pos:pos + 32]
        fields.append((desc[:11].split(b'\0')[0].decode('ascii'), desc[16]))
        pos += 32
    assert 'UNIT["Meter",1.0]' in path.with_suffix('.prj').read_text()
    offset = 100
    row = 0
    while offset < len(raw):
        number, words = struct.unpack_from('>ii', raw, offset)
        content = raw[offset + 8:offset + 8 + words * 2]
        record = dbf[header + row * reclen:header + (row + 1) * reclen]
        attributes, col = {}, 1
        for name, size in fields:
            attributes[name] = record[col:col + size].decode('ascii').strip()
            col += size
        row += 1
        offset += 8 + words * 2
        if record[:1] == b'*':
            continue
        kind, = struct.unpack_from('<i', content)
        assert kind == 13, kind
        nparts, npoints = struct.unpack_from('<ii', content, 36)
        starts = list(struct.unpack_from('<' + 'i' * nparts, content, 44)) + [npoints]
        point_start = 44 + 4 * nparts
        xy = np.frombuffer(content, dtype='<f8', count=npoints * 2, offset=point_start).reshape(-1, 2).copy()
        z = np.frombuffer(content, dtype='<f8', count=npoints, offset=point_start + npoints * 16 + 16).copy()
        yield number, attributes, [(xy[a:b], z[a:b]) for a, b in zip(starts[:-1], starts[1:])]
    assert row == count


def circle_fit(xy):
    origin = xy.mean(axis=0)
    p = xy - origin
    a = np.column_stack((2 * p, np.ones(len(p))))
    b = np.sum(p * p, axis=1)
    coef, _, rank, _ = np.linalg.lstsq(a, b, rcond=None)
    if rank < 3:
        return math.inf, math.inf
    radius = math.sqrt(max(0, coef[2] + np.dot(coef[:2], coef[:2])))
    residual = np.linalg.norm(p - coef[:2], axis=1) - radius
    return radius, float(np.max(np.abs(residual)))


def local_radii(xy, stride=1):
    # Circumcircle at every triple; stride sensitivity diagnoses digitizing noise.
    a, b, c = xy[:-2 * stride], xy[stride:-stride], xy[2 * stride:]
    ab, ac, bc = b - a, c - a, c - b
    cross = np.abs(ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0])
    product = np.linalg.norm(ab, axis=1) * np.linalg.norm(ac, axis=1) * np.linalg.norm(bc, axis=1)
    radii = np.full(len(cross), np.inf)
    np.divide(product, 2 * cross, out=radii, where=cross > 1e-9)
    return radii


def speed(radius, bank):
    return math.sqrt(9.80665 * radius * math.tan(math.radians(bank))) / (1852 / 3600)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = []
    hashes = {}
    for label in ('In', 'Out'):
        path = SOURCE / f'IAF_FAF_{label}_Update.shp'
        for suffix in ('.shp', '.shx', '.dbf', '.prj'):
            file = path.with_suffix(suffix)
            hashes[file.name] = hashlib.sha256(file.read_bytes()).hexdigest()
        for row, attr, parts in read_lines(path):
            for part, (xy, z) in enumerate(parts, 1):
                keep = np.r_[True, np.linalg.norm(np.diff(xy, axis=0), axis=1) > 1e-7]
                xy, z = xy[keep], z[keep]
                steps = np.linalg.norm(np.diff(xy, axis=0), axis=1)
                s = np.r_[0, steps.cumsum()]
                rlocal = local_radii(xy)
                i = int(np.argmin(rlocal)) + 1
                rmin = float(rlocal.min())
                rfit, residual = circle_fit(xy)
                result = dict(layer=label, id=attr['id'], row_1based=row, part=part,
                              vertices=len(xy), length_xy_m=float(s[-1]),
                              start_z_m=float(z[0]), end_z_m=float(z[-1]),
                              average_abs_slope_deg=math.degrees(math.atan(abs(float(z[-1] - z[0])) / s[-1])),
                              max_abs_slope_deg=float(np.degrees(np.arctan(np.max(np.abs(np.diff(z) / steps))))),
                              circle_fit_radius_m=rfit, circle_fit_max_residual_m=residual,
                              min_local_radius_m=rmin, median_local_radius_m=float(np.median(rlocal)),
                              min_radius_station_from_stored_start_m=float(s[i]),
                              min_radius_x=float(xy[i, 0]), min_radius_y=float(xy[i, 1]))
                for stride in (2, 3, 5):
                    result[f'min_radius_stride_{stride}_m'] = float(local_radii(xy, stride).min()) if len(xy) > 2 * stride else None
                for bank in (10, 15, 20, 25, 30):
                    result[f'v_bank_{bank}_kt'] = speed(rmin, bank)
                for v in (60, 70, 80, 90, 100, 120):
                    result[f'bank_at_{v}kt_deg'] = math.degrees(math.atan((v * 1852 / 3600) ** 2 / (9.80665 * rmin)))
                result['max_abs_vertical_speed_at_60kt_horizontal_fpm'] = (
                    60 * 1852 / 3600 * float(np.max(np.abs(np.diff(z) / steps))) * 60 / 0.3048)
                results.append(result)
    target = OUTPUT / 'RF_Turn_ID별_속도.csv'
    with target.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    (OUTPUT / 'input_sha256.json').write_text(json.dumps(hashes, indent=2), encoding='utf-8')
    report = [
        '# RF Turn 곡률 기반 속도 추정', '',
        '입력: 자료/최종정리/IAF_FAF_In_Update.shp 및 IAF_FAF_Out_Update.shp. 각 11개, 총 22개 활성 피처.',
        '좌표계: PRJ의 KGD2002 Central Belt 2010, 평면 거리 단위 m. 속성 id를 사용하며 중복은 원본 1-based 행으로 구분.', '',
        '## 방법과 가정', '',
        '- XY 전체 점에 원을 최소제곱 적합하고 인접 세 꼭짓점 외접원 반경으로 교차 검증했다. 2/3/5 간격 꼭짓점에서도 반경을 확인했다.',
        '- 22개 모두 단일 일정 반경 원호다. 전체 원 적합의 최대 반경 잔차는 4.2e-10m 미만이며, 꼭짓점 간격을 바꾼 반경 차이는 상대값 1.4e-9 미만이다. 이 값은 기하 적합 오차이며 실측 정확도가 아니다.',
        '- V = sqrt(g R tan(phi)), g=9.80665m/s², 1kt=1852/3600m/s. 선회반경은 수평 반경이다.',
        '- 무풍, 정상 협조 수평선회 근사. 바람/경사/롤 전이/기체 성능/항법오차/장애물 보호구역을 포함한 운항 허용속도는 아니다.',
        '- 무풍에서 표의 속도는 TAS와 GS가 같다고 놓은 수평 속도다. IAS/CAS 제한으로 직접 사용할 수 없다.',
        '- 15/20/25°는 민감도 분석 가정이다. UAM의 일률적인 법적 뱅크각 제한으로 해석하지 않는다.',
        '- 참고: FAA Order 8260.58D, 1-2-5, Formula 1-2-11 및 RF construction. RF 절차 설계의 반경/속도/뱅크각 관계와 바람 고려를 참고했다.',
        '  https://www.faa.gov/documentLibrary/media/Order/Order_8260.58D.pdf',
        '  https://www.faa.gov/regulations_policies/orders_notices/index.cfm/go/document.information/documentID/1043458', '',
        '## ID별 결과', '',
        '| 파일 | id | 원본 행(1부터) | R(m) | 15°(kt) | 20°(kt) | 25°(kt) | 최대 절대경사(°) |',
        '|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in results:
        report.append(f"| {r['layer']} | {r['id']} | {r['row_1based']} | {r['min_local_radius_m']:.1f} | {r['v_bank_15_kt']:.1f} | {r['v_bank_20_kt']:.1f} | {r['v_bank_25_kt']:.1f} | {r['max_abs_slope_deg']:.2f} |")
    report += ['', '## 해석', '',
               '- 전체 공통 속도를 제약하는 경로는 In id=1, R=350.367m. 뱅크각 15/20/25°에서 각각 58.98/68.74/77.81kt.',
               '- Out만 비교하면 id=3의 원본 4행, R=628.759m가 제약한다. 15/20/25°에서 79.01/92.09/104.23kt.',
               '- 연구용 무풍 공통속도 후보: 15° 가정 55kt, 20° 가정 60~65kt, 25° 가정 70~75kt. 이는 기하 상한보다 낮게 정한 비교 후보이며 실제 운항 승인을 의미하지 않는다.',
               '- 순풍 10kt를 TAS에 가산하는 단순 보수 설계 예에서는 20°의 TAS 상한이 68.74-10=58.74kt가 된다. 실제 RF 바람 검토는 선회 중 방위와 풍향에 따른 요구 뱅크각을 계산해야 한다.',
               '- Z는 저장값을 사용했다. In의 대다수 IAF는 457.2m, 일부는 609.6m로, 이전 대화의 고도 하향안은 이 파일에 반영되어 있지 않다.',
               '- In id=1 최대 경사 14.12°, In id=5는 18.01°. 수평속도 60kt에서 각각 약 1,528/1,976ft/min의 수직속도가 필요하다. Out은 저장 방향에서 대체로 상승하므로 경사는 절대값으로 표기했다.',
               '- Z값에 따른 수직 비행성능과 접속부 접선/롤 진입·이탈을 별도로 검토해야 한다. 본 계산은 접속부를 포함한 RF 인증 적합성 검토가 아니다.',
               '- 외곽·내곽 고도만 낮추고 XY를 유지하면 표의 수평 선회반경 상한은 동일하지만 수직성능 요구는 변한다.', '']
    (OUTPUT / 'RF_Turn_속도_분석보고서.md').write_text('\n'.join(report), encoding='utf-8')
    print('layer id row n length Rmin Rmedian Rfit residual V15 V20 V25 slope')
    for r in results:
        print(r['layer'], r['id'], r['row_1based'], r['vertices'], *[round(r[k], 3) for k in (
            'length_xy_m', 'min_local_radius_m', 'median_local_radius_m', 'circle_fit_radius_m',
            'circle_fit_max_residual_m', 'v_bank_15_kt', 'v_bank_20_kt', 'v_bank_25_kt', 'max_abs_slope_deg')])


if __name__ == '__main__':
    main()
