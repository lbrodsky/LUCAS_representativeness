#!/usr/bin/env python3

# Compute representativeness stats per country

import os
import argparse
from pathlib import Path

import pandas as pd
import geopandas as gpd
from joblib import Parallel, delayed

from utils import latest_version, country_codes

def main(data_dir, version, products):
    """Merge representative areas and other products for all EU countries
    """
    print('Starting stat analysis')
    df = pd.DataFrame(columns=['Country', 'Version', 'No_points',
                               'RG', 'RG_pct', 'No_match', 'No_match_pct',
                               'Multiclass', 'Multiclass_pct',
                               'Geom_update', 'Geom_update_pct', 'Sim_1', 'Sim_1_pct'])

    for cntr in sorted(country_codes.keys()):
        cntr_code = country_codes[cntr]
        # original points
        cntr_fn = os.path.join(data_dir, f'{cntr}_lucas_representativeness.gpkg')
        df_original_points = gpd.read_file(cntr_fn, layer=f'{cntr_code}_lucas_original_points')

        # region grow polygons
        df_region_grow = gpd.read_file(cntr_fn, layer=f'{cntr_code}_lucas_region_grow')

        # no match points
        df_nomatch_points = gpd.read_file(cntr_fn, layer=f'{cntr_code}_lucas_nomatch_points')

        # analysis
        # rg / no match points
        number_points = df_original_points.shape[0]
        number_procesed = df_region_grow.shape[0]
        no_processed = df_nomatch_points.shape[0]
        print(f'There are {number_procesed} points processed, but {no_processed} points did not match out of {number_points}.')

        # multi-class cases
        number_multiclass = (df_region_grow['multiclass'] > 1).sum()
        print(f'There are multiclass {number_multiclass} points out of {number_points}.')
        print(f'It is {round(number_multiclass / number_points * 100, 2)} %')

        # Updated points
        number_updated_points = (df_region_grow['point_update'] > 0).sum()
        print(f'There are {number_updated_points} updated points out of {number_points}.')
        print(f'It is {round(number_updated_points / number_points * 100, 2)} %')

        # Class similarity
        sim_1 = (df_region_grow['similarity'] == 1).sum()
        print(f'Detected similarity ~ 1 = {sim_1} out of {number_points}.')

        cntr_name = cntr.replace('-', ' ').capitalize()
        if number_points > 0:
            row = {'Country': cntr_name, 'Version': version, 'No_points': number_points,
                   'RG': number_procesed, 'RG_pct': round(number_procesed / number_points * 100, 1),
                   'No_match': no_processed, 'No_match_pct': round(no_processed / number_points * 100, 1),
                   'Multiclass': number_multiclass, 'Multiclass_pct': round(number_multiclass / number_points * 100, 1),
                   'Geom_update': number_updated_points, 'Geom_update_pct': round(number_updated_points / number_points * 100, 1),
                   'Sim_1': sim_1, 'Sim_1_pct': round(sim_1 / number_points * 100, 1)}
        else:
            row = {'Country': cntr_name, 'Version': version, 'No_points': number_points,
                   'RG': number_procesed, 'RG_pct': 0,
                   'No_match': no_processed, 'No_match_pct': 0,
                   'Multiclass': number_multiclass, 'Multiclass_pct': 0,
                   'Geom_update': number_updated_points, 'Geom_update_pct': 0,
                   'Sim_1': 0, 'Sim_1_pct': 0}

        df = pd.concat([df, pd.DataFrame.from_dict([row])], ignore_index=True)

    eu_lucas_repre_stat_fn = os.path.join(data_dir, 'european-union_lucas_representativeness_stat.csv')
    print(f"Stats saved in {eu_lucas_repre_stat_fn}")
    df.to_csv(eu_lucas_repre_stat_fn, sep=';')
    print('Done.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', help='Directory with region grow representative polygons for countries.')
    parser.add_argument('--workers', type=int, default=8,
                        help='Number of workers for parallel CPU processing.')

    args = parser.parse_args()

    # determine version from log file
    version = -1
    for fn in Path(args.data_path).glob('log_merge_countries_v*.txt'):
        version = fn.stem.split('_')[-1].lstrip('v')

    main(args.data_path, version, ['*_lucas_representativess.gpkg'])
