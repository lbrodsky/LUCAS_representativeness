#!/usr/bin/env python3

# Gridding representativness polygons to sensor spatial resolution, e.g. Sentinel-2 

import os 
import glob
import argparse
from shapely.geometry import Polygon
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
# from shapely import intersection
import shapely
from shapely.ops import unary_union

"""
# Example: Sentinel-2 gridding  

python3 ./utils/gridding_repre_polygons.py \
    -s /lucas/lucas_RG_2018/ \
    -d /lucas/lucas_RG_2018/ \
    -g 10   
    
"""


parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='LC nomenclature without class 13 (others)')
parser.add_argument('-s', '--src_path', metavar='src_path', type=str, help='Path to directory with source files')
parser.add_argument('-d', '--dst_path', metavar='dst_path', type=str, help='Path to directory with source files')
parser.add_argument('-g', '--grid_size', metavar='grid_size', default=10, type=int, help='Grid/pixel size of target sensor')
args = parser.parse_args()


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


def main(repre, dst_path, grid_size=10):
    """Gridding repre polygons main process.
    """

    repre_poly = gpd.read_file(repre)
    repre_poly_s2 = repre_poly.copy()

    # loop over features
    i = 0
    for feature_ix in range(len(repre_poly)):
        i += 1
        print(f'Point ID: {repre_poly.loc[feature_ix, "point_id"]} | {i}/{repre_poly.shape[0]}')
        # get selected feature geom
        repre_geom = repre_poly.loc[feature_ix, 'geometry']
        # create grid witin repre polygon
        grid = partition(repre_geom, grid_size)
        # cleaned grid for Sentinel-2
        grid_s2 = remove_small(grid, grid_size)
        # merge Sentinel-2 polygon pixels
        merged_grid_s2 = unary_union(grid_s2)
        # replace original repre polygon with S2 polygon
        repre_poly_s2.loc[feature_ix, "geometry"] = merged_grid_s2


    # save result
    if os.path.basename(repre).split('.')[-1] == 'shp':
        repre_poly_s2.to_file(repre.replace('.shp', '_s2.gpkg'), driver="GPKG")
    elif os.path.basename(repre).split('.')[-1] == 'gpkg':
        repre_poly_s2.to_file(repre.replace('.gpkg', '_s2.gpkg'), driver="GPKG")
    else:
        print('Differernt file format!')


if __name__ == "__main__":

    args = parser.parse_args()
    cfg = {}
    cfg['src_path'] = args.src_path
    cfg['dst_path'] = args.dst_path

    # PARAMS
    cfg['grid_size'] = args.grid_size
    # grid_size =  # 10.0
    shp_fn = None
    try:
        shp_fn = glob.glob(os.path.join(cfg['src_path'], '*region_grow_recoded.gpkg'))
    except:
        shp_fn = glob.glob(os.path.join(cfg['src_path'], '*region_grow.shp'))

    if not shp_fn:
        print('RG not available, use COP')
        shp_fn = glob.glob(os.path.join(cfg['src_path'], '*COP_2018.gpkg'))
    shp_fn.sort()
    for repre in shp_fn:
        print(repre)
        main(repre, cfg['dst_path'], cfg['grid_size'])
        
