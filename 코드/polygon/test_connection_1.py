"""선택 행, 삭제 레코드와 원본 보존 검사."""
from pathlib import Path
import struct
import tempfile
import unittest

import shapefile

import connection_1 as connection


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'lines.shp'
        with shapefile.Writer(str(self.path), shapeType=shapefile.POLYLINEZ) as writer:
            writer.field('id', 'N')
            writer.linez([[(0, 0, 0), (2, 0, 0)]])
            writer.record(9)
            writer.linez([[(100, 0, 1), (102, 0, 2)]])
            writer.record(9)
        self.groups = {'FAF': [(0, 0, 10)], 'IAF': [(2, 0, 20)]}

    def test_explicit_selection_ignores_other_active_rows(self):
        self.assertEqual(connection.report_matches(self.path, self.groups, .01, [1]), [0])
        with self.assertRaisesRegex(ValueError, 'FID=1, id=9'):
            connection.report_matches(self.path, self.groups, .01)

    def test_deleted_geometry_is_preserved_without_matching(self):
        dbf = self.path.with_suffix('.dbf')
        data = bytearray(dbf.read_bytes())
        header, record_size = struct.unpack_from('<HH', data, 8)
        data[header + record_size] = ord('*')
        dbf.write_bytes(data)
        selected = connection.report_matches(self.path, self.groups, .01)
        self.assertEqual(selected, [0])
        output = self.path.with_name('staged.shp')
        connection.stage_selected_lines(self.path, output, self.groups['FAF'],
                                        self.groups['IAF'], .01, selected)
        self.assertEqual(dbf.read_bytes(), data)
        with connection.open_reader(output, geometry_only=True) as reader:
            shapes = list(reader.iterShapes())
        self.assertEqual(list(shapes[0].z), [10, 20])
        self.assertEqual(list(shapes[1].points), [(100, 0), (102, 0)])
        self.assertEqual(list(shapes[1].z), [1, 2])
        with self.assertRaisesRegex(ValueError, '삭제'):
            connection.report_matches(self.path, self.groups, .01, [2])

    def test_selection_range_validation(self):
        with self.assertRaises(ValueError):
            connection.report_matches(self.path, self.groups, .01, [0])
        with self.assertRaises(ValueError):
            connection.report_matches(self.path, self.groups, .01, [3])


if __name__ == '__main__':
    unittest.main()
