import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

class TestSample:
    def test_001(self):
        """Run computation. Check output directory content."""
        from region_grow.process_tiles import define_parser, process_tiles

        parser = define_parser()
        args = parser.parse_args([
            '--tiles_dir', '/data/tiles',
            '--tables_dir', '/opt/region_grow/tables',
            '--lucas_points', '/data/sample_points.gpkg',
            '--lucas_thr_points', '/data/sample_thr_points.gpkg',
            '--dst_dir', '/tmp/output',
            '--shp_thr', '0.5',
            '--region_max_size', '100',
            '--version', '1'
        ])
        process_tiles(args)
