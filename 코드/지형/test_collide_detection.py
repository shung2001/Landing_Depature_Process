"""합성 도형으로 수치 판정과 판정 불가 처리의 회귀를 검증한다.

python -B 코드/지형/test_collide_detection.py
임시 검증 파일도 자료/결과물/_temp 아래에서만 만들고 종료 시 제거한다.
"""

import math
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
from shapely.geometry import GeometryCollection, LineString, MultiPolygon, Point, box

import collide_detection as cd


def buildings(geometries, tops):
    size = len(geometries)
    return gpd.GeoDataFrame({
        "__cd_id": [f"B{i}" for i in range(size)], "__cd_fid": list(range(size)),
        "__cd_top": tops, "__cd_ground": np.zeros(size), "__cd_agl": tops,
        "A16": tops, "A26": [2] * size, "name": ["한글 건물"] * size,
    }, geometry=geometries, crs=cd.ANALYSIS_CRS)


def helipads(geometries):
    return gpd.GeoDataFrame({
        "__cd_id": [f"H{i}" for i in range(len(geometries))],
        "__cd_fid": list(range(len(geometries))),
    }, geometry=geometries, crs=cd.ANALYSIS_CRS)


def match(status="matched"):
    return pd.DataFrame([{
        "helipad_id": "H0", "source_fid": 0, "host_building_id": "HOST",
        "host_source": "buildings", "overlap_ratio": 1.0, "helipad_z_m": 10.0,
        "center_x_m": 0.0, "center_y_m": 0.0, "status": status,
    }])


def dem(values, valid=None, geotransform=(0, 1, 0, 2, 0, -1)):
    values = np.asarray(values, dtype=float)
    return cd.DemData(values, np.isfinite(values) if valid is None else np.asarray(valid),
                      geotransform, cd.CRS(cd.ANALYSIS_CRS), "synthetic")


class CollisionTests(unittest.TestCase):
    def setUp(self):
        self.config = cd.SurfaceConfig(100, 20, 20, 10, 0)

    def test_bearings_and_start_offset(self):
        config = cd.SurfaceConfig(100, 20, 40, 10, 5)
        north = cd.create_approach_polygon(Point(0, 0), 0, config)
        east = cd.create_approach_polygon(Point(0, 0), 90, config)
        np.testing.assert_allclose(north.bounds, [-20, 5, 20, 105])
        np.testing.assert_allclose(east.bounds, [5, -20, 105, 20], atol=1e-12)
        self.assertAlmostEqual(cd.calculate_min_surface_z(
            Point(0, 15), Point(0, 0), 0, 50, config), 51)

    def test_actual_intersection_vertices_not_building_centroid(self):
        # 건물의 원래 꼭짓점은 모두 면 밖. 새 교점으로 최소 높이를 구해야 한다.
        footprint = box(-50, 20, 50, 30)
        intersection = footprint.intersection(cd.create_approach_polygon(Point(0, 0), 0, self.config))
        self.assertAlmostEqual(cd.calculate_min_surface_z(
            intersection, Point(0, 0), 0, 10, self.config), 12)

    def test_vertices_of_multipart_holes_lines_points(self):
        geometry = GeometryCollection([
            MultiPolygon([box(-1, 20, 1, 30), box(-1, 40, 1, 50)]),
            LineString([(0, 15), (0, 16)]), Point(0, 12),
        ])
        self.assertAlmostEqual(cd.calculate_min_surface_z(
            geometry, Point(0, 0), 0, 10, self.config), 11.2)
        self.assertTrue(math.isnan(cd.calculate_min_surface_z(
            GeometryCollection(), Point(0, 0), 0, 10, self.config)))

    def test_all_angles_clear_collision_and_ground_reference(self):
        frame = buildings([box(-2, 20, 2, 25), box(1000, 1000, 1001, 1001)], [13.0, 1000.0])
        details, angles, diagnostics = cd.check_surface_collisions(frame, match(), self.config)
        self.assertEqual(len(angles), 180)
        self.assertEqual(set(angles.axis_angle_deg), set(range(180)))
        axis0 = angles.loc[angles.axis_angle_deg.eq(0)].iloc[0]
        axis90 = angles.loc[angles.axis_angle_deg.eq(90)].iloc[0]
        self.assertEqual(axis0.collision_count, 1)
        self.assertTrue(axis0.is_collision)
        self.assertFalse(axis90.is_collision)
        self.assertEqual(axis90.collision_count, 0)
        self.assertEqual(axis90.max_penetration_m, 0)
        event = details.loc[details.axis_angle_deg.eq(0)].iloc[0]
        self.assertAlmostEqual(event.penetration_m, 1)
        self.assertTrue(diagnostics.status.eq("evaluated").all())

    def test_host_excluded_two_sides_count_unique_building(self):
        frame = buildings([box(-2, -25, 2, 25), box(-1, -1, 1, 1)], [13.0, 500])
        frame.loc[1, "__cd_id"] = "HOST"
        details, angles, _ = cd.check_surface_collisions(frame, match(), self.config)
        self.assertEqual(len(details.loc[details.axis_angle_deg.eq(0)]), 2)
        self.assertEqual(angles.iloc[0].collision_count, 1)
        self.assertNotIn("HOST", set(details.building_id))
        summary = cd._building_summary(details)
        self.assertEqual(summary.iloc[0].collision_angle_count, 180)
        self.assertEqual(summary.iloc[0].collision_axis_ranges_deg, "0-179")

    def test_unknown_heights_are_not_collision_free(self):
        frame = buildings([box(-1, 20, 1, 25)], [math.nan])
        details, angles, diagnostics = cd.check_surface_collisions(frame, match(), self.config)
        self.assertTrue(details.empty)
        self.assertTrue(pd.isna(angles.iloc[0].is_collision))
        self.assertTrue(pd.isna(angles.iloc[0].collision_count))
        self.assertEqual(diagnostics.iloc[0].unknown_building_count, 1)
        self.assertFalse(angles.loc[angles.axis_angle_deg.eq(90)].iloc[0].is_collision)

    def test_known_collision_with_unknown_candidate(self):
        frame = buildings([box(-1, 20, 1, 25), box(-1, 30, 1, 35)], [15, math.nan])
        _, angles, diagnostics = cd.check_surface_collisions(frame, match(), self.config)
        self.assertTrue(angles.iloc[0].is_collision)
        self.assertTrue(pd.isna(angles.iloc[0].collision_count))
        self.assertEqual(diagnostics.iloc[0].known_collision_count, 1)

    def test_unmatched_and_missing_config_have_180_unknown_rows(self):
        frame = buildings([box(-1, 20, 1, 25)], [50])
        for matches, config in [(match("no_spatial_host"), self.config), (match(), None)]:
            details, angles, _ = cd.check_surface_collisions(frame, matches, config)
            self.assertTrue(details.empty)
            self.assertEqual(len(angles), 180)
            self.assertTrue(angles.is_collision.isna().all())

    def test_equality_and_vertical_tolerance_not_penetration(self):
        frame = buildings([box(-1, 20, 1, 25)], [12])
        details, angles, _ = cd.check_surface_collisions(frame, match(), self.config)
        self.assertFalse(angles.iloc[0].is_collision)
        self.assertTrue(details.loc[details.axis_angle_deg.eq(0)].empty)

    def test_ranges_sorted_deduplicated_without_wrap(self):
        self.assertEqual(cd.merge_consecutive_angles([37, 32, 31, 33, 34, 35, 36, 112, 113, 113]),
                         "31-37; 112-113")
        self.assertEqual(cd.merge_consecutive_angles([179, 0, 1]), "0-1; 179")
        self.assertEqual(cd.merge_consecutive_angles([]), "")

    def test_invalid_surface_parameters(self):
        for values in [(0, 20, 20, 5, 0), (10, 20, 10, 5, 0), (10, 20, 20, -5, 0),
                       (10, 20, 20, 5, -1), (10, 20, 20, math.nan, 0)]:
            with self.assertRaises(ValueError):
                cd.SurfaceConfig(*values)


class DemAndMatchingTests(unittest.TestCase):
    def test_dem_touch_boundary_excluded_and_nodata_not_zero(self):
        raster = dem([[10, 20], [30, 40]])
        self.assertEqual(cd.get_dem_stat_for_geometry(box(0.1, 1.1, 0.9, 1.9), raster).elevation_m, 10)
        self.assertEqual(cd.get_dem_stat_for_geometry(box(0, 1, 1, 2), raster).valid_cell_count, 1)
        self.assertEqual(cd.get_dem_stat_for_geometry(box(0, 1, 2, 2), raster).elevation_m, 15)
        self.assertEqual(cd.get_dem_stat_for_geometry(Point(0.5, 0.5), raster).elevation_m, 30)
        self.assertTrue(math.isnan(cd.get_dem_stat_for_geometry(box(-1, 1, 1, 2), raster).elevation_m))
        raster.valid[0, 0] = False
        stat = cd.get_dem_stat_for_geometry(box(0, 1, 2, 2), raster)
        self.assertTrue(math.isnan(stat.elevation_m))
        self.assertAlmostEqual(stat.coverage_ratio, 0.5)

    def test_rotated_geotransform(self):
        raster = dem([[10, 20], [30, 40]], geotransform=(10, 0, 1, 20, 1, 0))
        stat = cd.get_dem_stat_for_geometry(box(10.1, 20.1, 10.9, 20.9), raster)
        self.assertEqual(stat.elevation_m, 10)

    def test_agl_plus_dem_and_opt_in_floor_estimate(self):
        frame = buildings([box(0.1, 1.1, 0.9, 1.9), box(1.1, 1.1, 1.9, 1.9)], [5, 0])
        raster = dem([[10, 20], [30, 40]])
        prepared, _ = cd.prepare_building_heights(frame, raster)
        self.assertEqual(prepared.iloc[0]["__cd_top"], 15)
        self.assertTrue(math.isnan(prepared.iloc[1]["__cd_top"]))
        prepared, diagnostic = cd.prepare_building_heights(frame, raster, floor_height_m=3)
        self.assertEqual(prepared.iloc[1]["__cd_top"], 26)
        self.assertEqual(diagnostic.iloc[1].height_source, "floor_estimate")

    def test_matching_overlap_ambiguous_and_no_nearest_fallback(self):
        frame = buildings([box(0, 0, 10, 10)], [20])
        pads = helipads([box(2, 2, 4, 4), box(10.01, 2, 11, 3)])
        matches, unmatched = cd.match_helipads_to_host_buildings(pads, frame, dem([[0]]))
        self.assertEqual(matches.iloc[0].host_building_id, "B0")
        self.assertEqual(matches.iloc[0].helipad_z_m, 20)
        self.assertEqual(unmatched.iloc[0].helipad_id, "H1")
        ambiguous = buildings([box(0, 0, 10, 10), box(0, 0, 10, 10)], [20, 30])
        matches, unmatched = cd.match_helipads_to_host_buildings(pads.iloc[:1], ambiguous, dem([[0]]))
        self.assertEqual(unmatched.iloc[0].reason, "ambiguous_spatial_hosts")

    def test_opt_in_reference_and_helipad_host_sources(self):
        raster = dem([[10, 20], [30, 40]])
        obstacles = buildings([box(1.1, 0.1, 1.9, 0.9)], [10])
        reference = buildings([box(0.1, 1.1, 0.9, 1.9)], [30])
        reference.loc[0, "__cd_id"] = "R0"
        pads = helipads([box(0.2, 1.2, 0.8, 1.8)])
        matches, unmatched = cd.match_helipads_to_host_buildings(
            pads, obstacles, raster, host_source="reference", reference=reference)
        self.assertTrue(unmatched.empty)
        self.assertEqual(matches.iloc[0].host_building_id, "R0")
        self.assertEqual(matches.iloc[0].helipad_z_m, 40)
        pads["A16"] = 50
        matches, unmatched = cd.match_helipads_to_host_buildings(
            pads, obstacles, raster, host_source="helipad")
        self.assertTrue(unmatched.empty)
        self.assertEqual(matches.iloc[0].host_building_id, "HOST_H0")
        self.assertEqual(matches.iloc[0].helipad_z_m, 60)

    def test_point_host_matching_and_crs_transform(self):
        frame = buildings([box(0, 0, 10, 10)], [20])
        pads = helipads([Point(5, 5)])
        matches, unmatched = cd.match_helipads_to_host_buildings(pads, frame, dem([[0]]))
        self.assertTrue(unmatched.empty)
        self.assertEqual(matches.iloc[0].host_building_id, "B0")
        raster = dem([[42]], geotransform=(200000, 100, 0, 550000, 0, -100))
        to_wgs84 = cd.Transformer.from_crs(cd.ANALYSIS_CRS, "EPSG:4326", always_xy=True)
        lon, lat = to_wgs84.transform(200050, 549950)
        stat = cd.get_dem_stat_for_geometry(Point(lon, lat), raster, geometry_crs="EPSG:4326")
        self.assertEqual(stat.elevation_m, 42)


class OutputTests(unittest.TestCase):
    def test_nonempty_and_empty_gpkg_and_csv_bom(self):
        frame = buildings([box(-1, 20, 1, 25)], [15])
        details, _, _ = cd.check_surface_collisions(frame, match(), cd.SurfaceConfig(100, 20, 20, 10, 0))
        cd.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        cd.TEMP_DIR.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="test_", dir=cd.TEMP_DIR) as directory:
            for label, events in [("nonempty", details), ("empty", details.iloc[:0])]:
                result = cd._collision_geodataframe(frame, events)
                target = Path(directory) / (label + ".gpkg")
                pyogrio.write_dataframe(result, target, layer="collision_buildings", driver="GPKG",
                                        geometry_type="MultiPolygon", promote_to_multi=True)
                reread = pyogrio.read_dataframe(target)
                self.assertEqual(len(reread), events.building_id.nunique())
                self.assertIn("collision_count", reread.columns)
                if not reread.empty:
                    self.assertEqual(reread.iloc[0]["name"], "한글 건물")
                    self.assertEqual(reread.iloc[0].collision_count, len(events))
                csv_path = Path(directory) / (label + ".csv")
                events.to_csv(csv_path, index=False, encoding="utf-8-sig")
                self.assertEqual(csv_path.read_bytes()[:3], b"\xef\xbb\xbf")
        try:
            cd.TEMP_DIR.rmdir()
        except OSError:
            pass


class DemBackendTests(unittest.TestCase):
    def require_cli(self):
        try:
            cd.find_gdal_cli()
        except cd.DemBackendUnavailable as error:
            self.skipTest(str(error))

    def test_auto_fallback_without_rasterio_or_osgeo(self):
        self.require_cli()
        if not cd.DEM_PATH.is_file():
            self.skipTest("실제 DEM 입력이 없는 환경")
        before = cd.input_fingerprints()
        expected = cd._load_dem(cd.DEM_PATH)
        # Python 패키지 설치 실패 상황을 재현하되 패키지를 삭제/변경하지 않는다.
        with patch.dict(cd.sys.modules, {"rasterio": None, "osgeo": None}):
            actual = cd._load_dem(cd.DEM_PATH)
        self.assertEqual(actual.backend, "gdal-cli")
        self.assertEqual(expected.crs, actual.crs)
        self.assertEqual(expected.geotransform, actual.geotransform)
        np.testing.assert_array_equal(expected.values, actual.values)
        np.testing.assert_array_equal(expected.valid, actual.valid)
        self.assertEqual(before, cd.input_fingerprints())

    def test_cli_preserves_float64_mask_nodata_scale_and_rotation(self):
        self.require_cli()
        try:
            import rasterio
        except ImportError:
            self.skipTest("합성 GeoTIFF 작성용 Rasterio 없음")
        with cd.output_workspace() as stage:
            source = stage / "한글 공간 test.tif"
            raw = np.array([[1.123456789012345, -9999, 3], [4, 5, 6]], dtype="float64")
            affine = rasterio.Affine(2, 0.2, 205920, 0.2, -2, 549360)
            with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
                with rasterio.open(source, "w", driver="GTiff", height=2, width=3,
                                   count=1, dtype="float64", crs=cd.ANALYSIS_CRS,
                                   transform=affine, nodata=-9999) as dataset:
                    dataset.write(raw, 1)
                    dataset.write_mask(np.array([[255, 255, 255], [0, 255, 255]], dtype="uint8"))
                    dataset.scales = (2.5,)
                    dataset.offsets = (7.0,)
            expected = cd._load_dem(source, "rasterio")
            actual = cd._load_dem(source, "gdal-cli")
            np.testing.assert_array_equal(expected.values, actual.values)
            np.testing.assert_array_equal(actual.values, raw * 2.5 + 7)
            np.testing.assert_array_equal(expected.valid, actual.valid)
            np.testing.assert_array_equal(actual.valid, [[True, False, True], [False, True, True]])
            self.assertEqual(expected.geotransform, actual.geotransform)
            self.assertEqual(expected.crs, actual.crs)

    def test_backend_unavailable_gives_binary_install_guidance(self):
        unavailable = cd.DemBackendUnavailable("not installed")
        with patch.object(cd, "_load_dem_rasterio", side_effect=unavailable), \
             patch.object(cd, "_load_dem_gdal", side_effect=unavailable), \
             patch.object(cd, "_load_dem_cli", side_effect=unavailable):
            with self.assertRaisesRegex(RuntimeError, "--only-binary=:all: rasterio"):
                cd._load_dem(cd.DEM_PATH)

    def test_bad_dem_does_not_silently_change_backend(self):
        with patch.object(cd, "_load_dem_rasterio", side_effect=ValueError("DEM 좌표계가 없습니다")), \
             patch.object(cd, "_load_dem_cli") as cli:
            with self.assertRaisesRegex(ValueError, "좌표계"):
                cd._load_dem(cd.DEM_PATH)
            cli.assert_not_called()


if __name__ == "__main__":
    unittest.main(verbosity=2)
