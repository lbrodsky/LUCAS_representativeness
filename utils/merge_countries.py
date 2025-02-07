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

def create_ds(filename, vector_format="GPKG"):
    """Create GDAL datasource."""
    driver = ogr.GetDriverByName(vector_format)
    if Path(filename).exists():
        driver.DeleteDataSource(filename)
    return driver.CreateDataSource(filename)

def merge_geometries(dst_fn, dst_ds_eu, data_dir, vector_format="GPKG"):
    """ogr2ogr merge over countries
    """
    print(f'Generating: {dst_fn}')

    # create country-based datasource
    dst_ds = create_ds(dst_fn)
    country_code = Path(dst_fn).name.split('_')[0]

    # open input data sources
    for pfn in Path(data_dir).glob("*.gpkg"):
        ds = ogr.Open(str(pfn))
        if ds is None:
            # empty datasource -> no LUCAS points in tile
            continue
        for lyr in ds:
            count = lyr.GetFeatureCount()
            if count < 1:
                continue

            layer_name = lyr.GetName()
            gdal.VectorTranslate(dst_ds.GetName(), ds.GetName(), format=vector_format,
                                 layers=[layer_name], accessMode='append',
                                 layerName=f'{country_code}_{layer_name}')
            gdal.VectorTranslate(dst_ds_eu.GetName(), ds.GetName(), format=vector_format,
                                 layers=[layer_name], accessMode='append',
                                 layerName=f'eu_{layer_name}')
        ds = None

    dst_ds = None

def main(dirs, dst_dir):
    """Merge representative areas and other products for all EU countries
    """
    basename = "lucas_representativeness"
    print('Starting merge procedure...')

    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    dst_fn_eu = os.path.join(dst_dir, f"eu_{basename}.gpkg")
    dst_ds_eu = create_ds(dst_fn_eu)
    print(f'Generating: {dst_fn_eu}')

    for cntr in dirs:
        versions = cntr.glob('v*')
        v_max = 0
        for v in versions:
            v_num = int(v.name[1:])
            if v_num > v_max: v_max = v_num
        if v_max > 0:
            src_dir = v.parent / f"v{v_max}"
            code = os.path.basename(os.path.dirname(src_dir)).split('_')[2]
            dst_fn = os.path.join(dst_dir, f"{code}_{basename}.gpkg")
            print(f"Processing {src_dir}...")
            merge_geometries(dst_fn, dst_ds_eu, src_dir)

    dst_ds_eu = None

    print('Done')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory',
                        required=True, help="Data directory")
    parser.add_argument('--target',
                        required=True, help="Target directory name")
    args = parser.parse_args()
    
    main(Path(args.directory).glob("*"), os.path.join(args.directory, args.target))
