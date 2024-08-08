#!/usr/bin/env

import os
from glob import glob
import shutil

from osgeo import gdal, ogr

def feature_count(fn):
    ds = ogr.Open(fn)
    lyr = ds.GetLayer()
    count = lyr.GetFeatureCount()
    ds = None

    return count

def merge_geometries(data_dir, dst_dir, prod, first):
    """ogr2ogr merge over countries
    """
    code = os.path.basename(os.path.dirname(data_dir)).split('_')[2]
    tokens = prod.split('_')
    dst_fn_eu = os.path.join(
        dst_dir,
        '_'.join(['eu', tokens[1], tokens[2], tokens[3]]))
    dst_fn = os.path.join(
        dst_dir,
        '_'.join([code, tokens[1], tokens[2], tokens[3]]))

    print(f'Processing:\n {dst_fn_eu}\n {dst_fn}')

    flag = ''
    vector_format = 'ESRI Shapefile'
    cnt = 0
    no_points = 0
    prod_files = glob(os.path.join(data_dir, prod))
    for pfn in prod_files:
        no_points += feature_count(pfn)
        if feature_count(pfn) < 1:
            continue

        gdal.VectorTranslate(dst_fn, pfn, format=vector_format,
                             accessMode='append' if cnt > 0 else 'overwrite')
        gdal.VectorTranslate(dst_fn_eu, pfn, format=vector_format,
                             accessMode='append' if not first else 'overwrite')
        cnt += 1

    print(f'Numbner of points in layer: {no_points}')
    if no_points == 0:
        print(f'Creating empty layer from {os.path.basename(pfn)}')
        # TODO: copy
        shutil.copy(pfn, dst_fn)
        shutil.copy(pfn.replace('.shp', '.dbf'), dst_fn.replace('.shp', '.dbf'))
        shutil.copy(pfn.replace('.shp', '.shx'), dst_fn.replace('.shp', '.shx'))
        shutil.copy(pfn.replace('.shp', '.prj'), dst_fn.replace('.shp', '.prj'))

def main(dirs, products, dst_dir):
    """Merge representative areas and other products for all EU countries
    """
    print('Starting merge')

    if not os.path.exists(dst_dir):
        os.makedirs(dst_dir)

    for prod in products:
        print('---')
        print(prod)
        print('---')
        prod_files = []
        first = True
        for cntr in dirs:
            versions = glob(os.path.join(cntr, 'v*'), recursive=True)
            v_max = 0
            for v in versions:
                v_num = int(v.split('/')[-1][1:])
                if v_num > v_max: v_max = v_num
            if v_max > 0:
                data_dir = versions[0] if len(versions) == 1 else versions[0][:versions[0].rfind('/')+2] + str(v_max)
                merge_geometries(data_dir, dst_dir, prod, first)
                first = False

if __name__ == "__main__":

    
    path = '/Users/lukas/Work/prfuk/ownCloud/tmp/lucas/RegionGrow/region_grow_es/v3'
    dst_dir = 'merge'
    dirs = glob(os.path.join(path, '*'), recursive=True)
    print(dirs)

    products = ['*_lucas_region_grow.shp', '*_lucas_urban_grow.shp',
                '*_lucas_updated_points.shp', '*_lucas_points_buffer.shp',
                '*_lucas_original_points.shp', '*_lucas_nomatch_points.shp']

    main(dirs, products, os.path.join(path, dst_dir))

