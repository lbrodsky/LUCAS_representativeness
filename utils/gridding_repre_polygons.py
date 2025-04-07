#!/usr/bin/env python3

# Gridding representativness polygons to sensor spatial resolution, e.g. Sentinel-2 

import os
import sys
import argparse
import logging
from pathlib import Path

import numpy as np
from osgeo import ogr
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.wkt import loads

from utils import latest_version

ogr.UseExceptions()

def grid_bounds(geom, size):
    """Get gridding bounds. 
    """ 
    minx, miny, maxx, maxy = geom.bounds
    minx_r = round(minx, -1)
    miny_r = round(miny, -1)
    maxx_r = round(maxx, -1)
    maxy_r = round(maxy, -1)
    nx = int((maxx_r - minx_r)/size) + 1
    ny = int((maxy_r - miny_r)/size) + 1
    gx, gy = np.linspace(minx_r, maxx_r, nx), np.linspace(miny_r, maxy_r, ny)
    grid = []
    for i in range(len(gx)-1):
        for j in range(len(gy)-1):
            poly_ij = Polygon([[gx[i],gy[j]],[gx[i],gy[j+1]],[gx[i+1],gy[j+1]],[gx[i+1],gy[j]]])
            grid.append(poly_ij)
    return grid


def partition(geom, size):
    """Partition spatial extent to grid.
    """
    grid = []
    grid_ = grid_bounds(geom, size)
    for g in grid_:
        isect = shapely.intersection(g, geom)
        grid.append(isect)
    return grid

def remove_small(grid, grid_size):
    """Remove small pixel polygons
    """
    remove_list = []
    i = 0
    for i in range(len(grid)):
        if grid[i].area < (grid_size * grid_size):
            # remove polygon
            remove_list.append(grid[i])

    grid_s2 = [i for i in grid if i not in remove_list]

    return grid_s2


def main(tile, grid_size=10):
    """Gridding repre polygons main process.
    """
    layer_name = "lucas_region_grow"
    ds = ogr.Open(tile, update=True)
    source_layer = ds.GetLayerByName(layer_name)
    if source_layer is None:
        return # skip empty tiles

    # Check if the layer exists and delete it
    target_layer_name = f"{layer_name}_gridded_s2"
    if ds.GetLayerByName(target_layer_name):
        ds.DeleteLayer(target_layer_name)
    target_layer = ds.CreateLayer(target_layer_name, source_layer.GetSpatialRef(), source_layer.GetGeomType())

    field_blacklist = [
        "rectangularity",
        "length",
        "width",
        "area",
        "ratio",
        "shape_gen"
    ]

    # Copy field definitions (schema) from source to target
    source_layer_def = source_layer.GetLayerDefn()
    for i in range(source_layer_def.GetFieldCount()):
        field_def = source_layer_def.GetFieldDefn(i)
        if field_def.GetName() not in field_blacklist:
            target_layer.CreateField(field_def)

    feat_idx = 1
    count = source_layer.GetFeatureCount()
    for feature in source_layer:
        logging.debug(f'Point ID: {feature.GetField("point_id")} | {feat_idx}/{count}')
        feat_idx += 1

        # get selected feature geom
        repre_geom = loads(feature.GetGeometryRef().ExportToWkt())
        # create grid witin repre polygon
        grid = partition(repre_geom, grid_size)
        # cleaned grid for Sentinel-2
        grid_s2 = remove_small(grid, grid_size)
        # merge Sentinel-2 polygon pixels
        merged_grid_s2 = unary_union(grid_s2)
        # replace original repre polygon with S2 polygon
        if merged_grid_s2.geom_type != 'Polygon':
            nparts = len(merged_grid_s2.geoms)
            if nparts == 0:
                continue # skip empty features

            # take polygon with largest area
            sel_idx = -1
            area_max = -1
            for idx, p in enumerate(merged_grid_s2.geoms):
                if p.area > area_max:
                    sel_idx = idx

            merged_grid_s2 = merged_grid_s2.geoms[sel_idx]

        new_feature = ogr.Feature(target_layer.GetLayerDefn())
        new_feature.SetGeometry(ogr.CreateGeometryFromWkt(merged_grid_s2.wkt))
        for i in range(source_layer_def.GetFieldCount()):
            if source_layer_def.GetFieldDefn(i).GetName() not in field_blacklist:
                new_feature.SetField(i, feature.GetField(i))

        target_layer.CreateFeature(new_feature)
        new_feature = None

    logging.info(f"Number of skipped features: {source_layer.GetFeatureCount()-target_layer.GetFeatureCount()}")

def init_logging(dst_dir, version):
    """Initialize the processing log.
    """
    log_file = os.path.join(
        dst_dir, f'log_gridding_repre_polygons_v{version}.txt'
    )
    log_level = logging.INFO

    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',
                        level=log_level,
                        handlers=[logging.FileHandler(log_file, mode='w'),
                                  logging.StreamHandler()])
    logging.info(f'Log file: {log_file}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Performs gridding of representative polygons to Sentinel-2.')
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to directory with source files')
    parser.add_argument('--version', type=int,
                        help='Version to process (default: latest)')
    parser.add_argument('--grid_size', default=10, type=int,
                        help='Grid/pixel size of target sensor')

    args = parser.parse_args()

    if args.version is None:
        version = latest_version(Path(args.data_path).glob("*"))
    else:
        version = args.version
    init_logging(args.data_path, version)

    try:
        gpkg_fn = list(Path(args.data_path).rglob(f'v{version}/*.gpkg'))
    except:
        logging.error('No original region grow files available.')

    count = len(gpkg_fn)
    i = 0
    for fn in sorted(gpkg_fn):
        i += 1
        logging.info(f"Processing tile {i} of {count}...\n{fn}")
        main(fn, args.grid_size)
