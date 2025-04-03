#!/usr/bin/env python3

# Gridding representativness polygons to sensor spatial resolution, e.g. Sentinel-2 

import os
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
pd.options.mode.copy_on_write = True
import geopandas as gpd
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
from pyogrio.errors import DataSourceError

from utils import latest_version

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


def main(tile, dst_path, grid_size=10):
    """Gridding repre polygons main process.
    """
    try:
        repre_poly = gpd.read_file(tile, layer="lucas_region_grow")
    except DataSourceError:
        return # skip empty files

    repre_poly_gridded = repre_poly.head(0)

    # loop over features
    for feature_ix in range(len(repre_poly)):
        print(f'Point ID: {repre_poly.loc[feature_ix, "point_id"]} | {feature_ix+1}/{repre_poly.shape[0]}', file=sys.stderr)
        # get selected feature geom
        repre_geom = repre_poly.loc[feature_ix, 'geometry']
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

        new_repre_poly = repre_poly.loc[feature_ix].copy()
        new_repre_poly["geometry"] = merged_grid_s2
        repre_poly_gridded.loc[len(repre_poly_gridded)] = new_repre_poly.dropna(axis=1, how="all")

    print(f"Number of skipped features: {len(repre_poly)-len(repre_poly_gridded)}",
          file=sys.stderr)

    # save result
    repre_poly_gridded.set_crs(repre_poly.crs)
    repre_poly_gridded.set_crs(repre_poly.crs).to_file(tile, layer="lucas_region_grow_gridded", driver="GPKG")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='LC nomenclature without class 13 (others)')
    parser.add_argument('--src_path', type=str, required=True,
                        help='Path to directory with source files')
    parser.add_argument('--version', type=int,
                        help='Version to process (default: latest)')
    parser.add_argument('--dst_path', type=str, required=True,
                        help='Path to target directory')
    parser.add_argument('--grid_size', default=10, type=int,
                        help='Grid/pixel size of target sensor')

    args = parser.parse_args()

    if args.version is None:
        version = latest_version(Path(args.src_path).glob("*"))
    else:
        version = args.version
    
    try:
        gpkg_fn = list(Path(args.src_path).rglob(f'v{version}/*.gpkg'))
    except:
        print('No original region grow files available.')

    count = len(gpkg_fn)
    i = 0
    for fn in sorted(gpkg_fn):
        i += 1
        print(f"Processing tile {i} of {count}...\n{fn}", file=sys.stderr)
        main(fn, args.dst_path, args.grid_size)
        if i == 3:
            break
