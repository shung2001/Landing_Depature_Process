# Vertical Separation Gauge

report.html에서 전체 결과를 확인할 수 있습니다.

- Short_Final_internal: Short_Final 선들 사이
- Vertical_routes_internal: Approach/Departure 선들 사이
- Short_Final_vs_vertical_routes: 두 자료 사이
- all_intersections.csv: 모든 분석 수치 (UTF-8 BOM)
- summary.csv / summary.json: 분석별 건수와 최소 수직 간격
- 각 폴더 crossing_heights_XX.png: 각 교점의 두 Z와 절대 높이 차
- 각 폴더 map_intersections.png / vertical_separation.png / vertical_links_3d.png

shared_endpoint_only=True는 두 선(part)의 끝점이 같은 XY에서 만난다는 뜻입니다.
공유 끝점의 0m 간격을 일반적인 내부 교차와 구분하되 원자료에서 삭제하지 않습니다.
중첩 선분은 공유 끝점을 갖더라도 중첩 분석에 남깁니다.
안전 등급·충돌 시각·기체 성능을 판정하는 코드는 아닙니다.
