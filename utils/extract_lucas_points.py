#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys
import argparse
from pathlib import Path

from osgeo import gdal

from utils import country_codes

def create_gpkg(output_file, layer_name, sql, dbname, dbuser):
    gdal.VectorTranslate(
        destNameOrDestDS=str(output_file),
        srcDS=f'PG:dbname={dbname} user={dbuser}',
        format='GPKG',
        accessMode='overwrite',
        datasetCreationOptions=[],
        layerName=layer_name,
        SQLStatement=sql,
    )

def main(dst_dir, year):
    if not dst_dir.exists():
        dst_dir.mkdir(parents=True, exist_ok=True)

    dbname = "lucas"
    dbuser = "postgres"   
    filebase = "lucas_points"

    for country_name, country in country_codes(year).items():
        print(f"Processing {country_name}...", file=sys.stderr)

        output_file = (
            dst_dir / f"{country_name}_{filebase}_gps_{year}.gpkg"
        )
        sql = f"""
        select
        point_id,
        geom,
        survey_lc1_h lc1_h,
        survey_gps_prec gps_prec,
        survey_obs_dist obs_dist,
        survey_obs_type obs_type,
        survey_obs_direct obs_direct
        from data.lucas_points
        where survey_year = {year}
        AND point_nuts0 = '{country}'
        """
        create_gpkg(output_file, 'eu_lucas_points_gps', sql, dbname, dbuser)

    # theoretical points
    thr_output = dst_dir / f"eu_{filebase}_thr_{year}.gpkg"
    thr_sql = f"""
    select
    point_id,
    geom_thr as geom,
    survey_lc1_h lc1_h,
    survey_gps_prec gps_prec,
    survey_obs_dist obs_dist,
    survey_obs_type obs_type,
    survey_obs_direct obs_direct
    from data.lucas_points
    where survey_year = {year}
    """
    create_gpkg(output_file, 'eu_lucas_points_thr', sql, dbname, dbuser)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dst_dir', help='Destination directory.')
    parser.add_argument('--year', help='Year to be downloaded.', type=int, default=2018)

    args = parser.parse_args()

    main(Path(args.dst_dir), args.year)
