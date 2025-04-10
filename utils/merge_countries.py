#!/usr/bin/env python

import os
import sys
import shutil
import argparse
import logging

from glob import glob
from pathlib import Path

from osgeo import gdal, ogr
gdal.DontUseExceptions()

from utils import latest_version, country_codes

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
    logging.info(f'Generating: {dst_fn}')

    # create country-based datasource
    dst_ds = create_ds(dst_fn)
    country_name = Path(dst_fn).name.split('_')[0]

    # open input data sources
    for pfn in Path(data_dir).glob("*.gpkg"):
        ds = ogr.Open(pfn)
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
                                 layerName=f'{country_codes[country_name].lower()}_{layer_name}')
            gdal.VectorTranslate(dst_ds_eu.GetName(), ds.GetName(), format=vector_format,
                                 layers=[layer_name], accessMode='append',
                                 layerName=f'eu_{layer_name}')
        ds = None

    dst_ds = None

def main(dirs, dst_dir, version=None):
    """Merge representative areas and other products for all EU countries
    """
    basename = "lucas_representativeness"
    logging.info('Starting merge procedure...')

    dst_fn_eu = os.path.join(dst_dir, f"european_union_{basename}.gpkg")
    dst_ds_eu = create_ds(dst_fn_eu)
    logging.info(f'Generating: {dst_fn_eu}')

    for cntr in dirs:
        if cntr.name == Path(dst_dir).name:
            # skip output dir
            continue
        if not cntr.is_dir():
            continue

        if version is None:
            version = latest_version(Path(cntr).parent.glob('*'))

        if version > 0:
            src_dir = cntr / f"v{version}"
            if not src_dir.exists():
                logging.warning(f"WARNING: {src_dir} doesn't exists. Skipped.", file=sys.stderr)
                continue
            code = src_dir.parent.name.split('_')[0]
            dst_fn = os.path.join(dst_dir, f"{code}_{basename}.gpkg")
            logging.info(f"Processing {src_dir}...")
            merge_geometries(dst_fn, dst_ds_eu, src_dir)

    dst_ds_eu = None

    logging.info('Done')

def init_logging(dst_dir, basename, version=None):
    """Initialize the processing log.
    """
    if version is not None:
        log_file = os.path.join(
            dst_dir, f'log_{basename}_v{version}.txt'
        )
    else:
        log_file = os.path.join(
            dst_dir, f'log_{basename}.txt'
        )

    log_level = logging.INFO

    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',
                        level=log_level,
                        handlers=[logging.FileHandler(log_file, mode='w'),
                                  logging.StreamHandler()])
    logging.info(f'Log file: {log_file}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--src_path', type=str, required=True,
                        help='Path to directory with source files')
    parser.add_argument('--version', type=int,
                        help='Version to process (default: latest)')
    parser.add_argument('--dst_path', type=str, required=True,
                        help='Path to target directory')
    args = parser.parse_args()

    if Path(args.dst_path).exists():
        shutil.rmtree(args.dst_path)
    if not os.path.exists(args.dst_path):
        os.makedirs(args.dst_path)

    init_logging(args.dst_path, Path(__file__).stem)

    main(Path(args.src_path).glob("*"), args.dst_path, args.version)
