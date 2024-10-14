#!/usr/bin/env python

import os
import shutil
import argparse
from glob import glob
from pathlib import Path

from osgeo import gdal, ogr

def feature_count(fn):
    ds = ogr.Open(fn)
    lyr = ds.GetLayer()
    count = lyr.GetFeatureCount()
    ds = None

    return count

def merge_geometries(data_dir, dst_dir, first):
    """ogr2ogr merge over countries
    """
    code = os.path.basename(os.path.dirname(data_dir)).split('_')[2]
    basename = "lucas_representativeness"
    dst_fn_eu = os.path.join(dst_dir, f"eu_{basename}.gpkg")
    dst_fn = os.path.join(dst_dir, f"{code}_{basename}.gpkg")

    print(f'Processing:\n {dst_fn_eu}\n {dst_fn}')

    vector_format = 'GPKG'
    cnt = 0
    no_points = 0

    for pfn in Path(data_dir).glob("*.gpkg"):
        ds = ogr.Open(str(pfn))
        if ds is None:
            # empty datasource -> no LUCAS points in tile
            continue
        for lyr in ds:
            count = lyr.GetFeatureCount()
            if count < 1:
                continue

            gdal.VectorTranslate(dst_fn, ds.GetName(), format=vector_format, layerName=lyr.GetName(),
                                 accessMode='append' if cnt > 0 else 'overwrite')
            gdal.VectorTranslate(dst_fn_eu, ds.GetName(), format=vector_format, layerName=lyr.GetName(),
                                 accessMode='append' if not first else 'overwrite')
            cnt += 1
            no_points += count
        ds = None

def main(dirs, dst_dir):
    """Merge representative areas and other products for all EU countries
    """
    print('Starting merge procedure...')

    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    first = True
    for cntr in dirs:
        versions = cntr.glob('v*')
        v_max = 0
        for v in versions:
            v_num = int(v.name[1:])
            if v_num > v_max: v_max = v_num
        if v_max > 0:
            merge_geometries(v.parent / f"v{v_max}", dst_dir, first)
            first = False
    print('Done')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory',
                        required=True, help="Data directory")
    parser.add_argument('--target',
                        required=True, help="Target directory name")
    args = parser.parse_args()
    
    main(Path(args.directory).glob("*"), os.path.join(args.directory, args.target))
