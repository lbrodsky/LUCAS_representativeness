#!/usr/bin/env python3

import os
import logging
import yaml
import csv
import pandas as pd
from osgeo import gdal
from osgeo import ogr
from representativeness_exceptions import ConfigError, IllegalArgumentError
gdal.UseExceptions()

def open_LUCAS_points(lucas_points_fn):
    """Open LUCAS vector file.
    """
    driver = ogr.GetDriverByName('GPKG')
    vec_ds = driver.Open(lucas_points_fn, 0)
    lucas_points_basename = os.path.basename(lucas_points_fn)
    if vec_ds is None:
        raise ConfigError(f'Could not open {vec_ds}')
    logging.debug(f'Opened {lucas_points_fn}')

    layer = vec_ds.GetLayer()
    featureCount = layer.GetFeatureCount()
    logging.info(
        f"Number of ALL features in {lucas_points_basename}: {featureCount}"
    )

    return vec_ds, layer


def read_yaml(fn):
    """Read yaml configuration.
    """
    with open(fn, 'r') as stream:
        tab = yaml.load(stream, Loader=yaml.FullLoader)

    return tab


def csv2dict(csv_fn, k='code', v='name', delim=';'):
    """Read csv translation table into dictionary.
    """
    d = {}
    with open(csv_fn) as f:
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            d[row[k]] = row[v]

    return d


def read_semantic_similarity(similarity_fn):
    """Read land cover semantic similarity table from csv into Pandas DataFrame.
    """
    df = pd.read_csv(similarity_fn, sep=';', header=0)

    return df
