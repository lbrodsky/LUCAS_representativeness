#!/usr/bin/env python3

# python3 ./utils/crop_tile.py -t tests/data/tiles/netherlands_osm_clcplus_14580.tif -l tests/data/sample_points.gpkg -lt tests/data/sample_thr_points.gpkg

import argparse
from pathlib import Path

from osgeo import gdal, ogr
gdal.UseExceptions()

def open_ds(filename):
    ds = ogr.Open(filename)

    return ds, ds.GetLayer()

def crop_tiles(tile_path, lucas_points, lucas_thr_points,
               selected_points=None, grow=500):
    ds_lp, lyr_lp = open_ds(lucas_points)
    ds_thr_lp, lyr_thr_lp = open_ds(lucas_thr_points)

    if selected_points is not None:
        lyr_lp.SetAttributeFilter("where {}".format(
            ' and '.join([f"point_id = {x}" for x in selected_points])
        ))
  
    while True:
        feat = lyr_lp.GetNextFeature()
        if feat is None:
            break
        point_id = feat.GetField("point_id")
        geom_lp = feat.GetGeometryRef()

        # get thr geometry
        lyr_thr_lp.ResetReading()
        lyr_thr_lp.SetAttributeFilter(f"point_id = {point_id}")
        feat_thr = lyr_thr_lp.GetNextFeature()
        geom_thr_lp = feat_thr.GetGeometryRef()

        # ulx, uly, lrx, lry
        extent = (min(geom_lp.GetX(), geom_thr_lp.GetX()) - grow,
                  max(geom_lp.GetY(), geom_thr_lp.GetY()) + grow,
                  max(geom_lp.GetX(), geom_thr_lp.GetX()) + grow,                  
                  min(geom_lp.GetY(), geom_thr_lp.GetY()) - grow)

        crop_tile(Path(tile_path), extent)

    ds_lp = None
    ds_lp_thr = None

def crop_tile(filepath, extent):
    ds = gdal.Open(str(filepath))

    output_filename = filepath.parent / f"{filepath.stem}_sample{filepath.suffix}"
    gdal.Translate(str(output_filename), ds, projWin=extent)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tile',
                        required=True, help='Tile filename.')
    parser.add_argument('-l', '--lucas_points',
                        required=True, help='LUCAS points filename.')
    parser.add_argument('-lt', '--lucas_thr_points', 
                        required=True, help='LUCAS theoretical points.')
    parser.add_argument('-s', '--selected_points', nargs='?',
                        help='List of selected points to be processed.')
    
    args = parser.parse_args()

    selected_points = None    
    if args.selected_points is not None:
        selected_points = args.selected_points.split(',')

    crop_tiles(args.tile, args.lucas_points, args.lucas_thr_points, selected_points)
