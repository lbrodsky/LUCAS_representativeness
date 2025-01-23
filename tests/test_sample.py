import sys
import tempfile
import pytest
from pathlib import Path

from osgeo import gdal, ogr
gdal.DontUseExceptions()

sys.path.insert(0, str(Path(__file__).parent.parent))
from region_grow.process_tiles import define_parser, process_tiles

@pytest.fixture(scope="class")
def fixture():
    output_path = tempfile.mkdtemp()
    version = 1
    gpkg_file = "netherlands_osm_clcplus_14580_sample.gpkg"

    return {
        'output_path': output_path,
        'version' : version,
        'gpkg_path': Path(output_path, f"v{version}", gpkg_file)
    }


class TestSample:
    def test_001(self, fixture):
        """Run computation. Check output directory content."""
        parser = define_parser()
        args = parser.parse_args([
            '--tiles_dir', '/data/tiles',
            '--tables_dir', '/opt/region_grow/tables',
            '--lucas_points', '/data/sample_points.gpkg',
            '--lucas_thr_points', '/data/sample_thr_points.gpkg',
            '--dst_dir', fixture['output_path'],
            '--shp_thr', '0.5',
            '--region_max_size', '100',
            '--version', str(fixture['version'])
        ])
        process_tiles(args)

        assert Path(fixture['output_path'], f"v{fixture['version']}").is_dir()
        assert Path(fixture['gpkg_path']).is_file()

    def test_002(self, fixture):
        """Open GPKG file by GDAL and check number of layers, layer names and geometry types."""
        ds = ogr.Open(fixture['gpkg_path'])
        assert ds is not None
        assert ds.GetDriver().GetName() == "GPKG"

        assert ds.GetLayerCount() == 4

        ref = {
            "lucas_points_buffer": "Polygon",
            "lucas_region_grow": "Polygon",
            "lucas_urban_grow": "Polygon",
            "lucas_original_points": "Point"
        }
        for lyr in ds:
            assert lyr.GetName() in ref.keys()
            assert ogr.GeometryTypeToName(lyr.GetGeomType()) == ref[lyr.GetName()]

        ds.Close()

    def test_003(self, fixture):
        """Check RG output layer."""
        with ogr.Open(fixture['gpkg_path']) as ds:
            lyr = ds.GetLayerByName("lucas_region_grow")

            assert lyr.GetFeatureCount() == 1

            feature = lyr.GetNextFeature()
            assert feature is not None

            # check attributes
            ref = {
                "point_id": 40283224,
                "tile_id": "netherlands_osm_clcplus_14580_sample.tif",
                "lc1_h": "C10",
                "lc1_name": "Broadleaved woodland",
                "osm_codes": "40",
                "osm_names": "Woodland",
                "osm_lc": None,
                "multiclass": 1,
                "geom_type": 1,
                "obs_type": "1", # why string?
                "obs_dist": 3,
                "gps_prec": 4,
                "gprec_src": "gps",
                "lc": 4,
                "similarity": 2,
                "point_update": 0,
                "rectangularity": 0.513,
                "length": 54,
                "width": 30,
                "area": 1164,
                "ratio": 0.555555555555556,
                "shape_gen": 1
            }
            for field_name in feature.keys():
                if field_name not in ref:
                    continue
                field_value = feature.GetField(field_name)
                if isinstance(ref[field_name], float):
                    assert ref[field_name] - field_value < 1e-12
                else:
                    assert ref[field_name] == field_value

            # check geometry
            assert feature.GetGeometryRef().GetArea() == ref["area"]
