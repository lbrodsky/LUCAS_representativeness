#!/usr/bin/env python

import os
import sys
import shutil
import argparse
from glob import glob
from pathlib import Path

from osgeo import gdal, ogr

from utils import latest_version

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

def main(dirs, dst_dir, version=None):
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
        if version > 0:
            src_dir = cntr / f"v{version}"
            if not src_dir.exists():
                print(f"WARNING: {src_dir} doesn't exists. Skipped.", file=sys.stderr)
                continue
            code = src_dir.parent.name.split('_')[0]
            dst_fn = os.path.join(dst_dir, f"{code}_{basename}.gpkg")
            print(f"Processing {src_dir}...")
            merge_geometries(dst_fn, dst_ds_eu, src_dir)

    dst_ds_eu = None

    print('Done')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--src_path', type=str, required=True,
                        help='Path to directory with source files')
    parser.add_argument('--version', type=int,
                        help='Version to process (default: latest)')
    parser.add_argument('--dst_path', type=str, required=True,
                        help='Path to target directory')
    args = parser.parse_args()

    if args.version is None:
        version = latest_version(Path(args.src_path).glob("*"))
    else:
        version = args.version

    if Path(args.dst_path).exists():
        shutil.rmtree(args.dst_path)

    main(Path(args.src_path).glob("*"), args.dst_path, version)
