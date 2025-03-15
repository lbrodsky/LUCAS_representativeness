#!/usr/bin/env python3

# Compute representativeness stats per country
# Example usage:
"""
python3 ./utils/country_repre_stats.py \
    -s ./repre_v1 \
    -d ./repre_v1/stats
"""

import os
import pandas as pd
import geopandas as gpd
import argparse
from joblib import Parallel, delayed


def main(dir, dst_path, version, products):
    """Merge representative areas and other products for all EU countries
    dir = '/Users/lukas/Work/prfuk/ownCloud/tmp/lucas/RegionGrow_EU/merge_eu_rect05_gen6_v40'
    dst_path = '/Users/lukas/Work/prfuk/ownCloud/tmp/lucas/RegionGrow_EU/merge_eu_rect05_gen6_v40/stats'
    """
    print('Starting stat analysis')
    df = pd.DataFrame(columns=['Country', 'Version', 'No_points',
                               'RG', 'RG_pct', 'No_match', 'No_match_pct',
                               'Multiclass', 'Multiclass_pct',
                               'Geom_update', 'Geom_update_pct', 'Sim_1', 'Sim_1_pct'])

    countries = ['at', 'be', 'bg', 'cy', 'cz', 'de', 'dk', 'ee', 'es', 'fi', 'fr', 'gb', 'gr', 'hr', 'hu', 'ie',
                 'it', 'lt', 'lu', 'lv', 'nl', 'mt', 'pl', 'pt', 'ro', 'se', 'si', 'sk']
    # countries = ['eu']
    v_max = version # int(os.path.basename(dir).split('_')[-1][1:])

    gpkg_products = ['lucas_region_grow',
                     'lucas_updated_points', 'lucas_points_buffer',
                     'lucas_original_points', 'lucas_nomatch_points']

    for cntr in countries:
        print(cntr)
        # original points
        # original_points_fn = os.path.join(dir, f'{cntr}_lucas_original_points.gpkg')
        cntr_fn = os.path.join(dir, f'{cntr}_lucas_representativeness.gpkg')
        # df_original_points = gpd.read_file(original_points_fn, layer='lucas_original_points')
        df_original_points = gpd.read_file(cntr_fn, layer='lucas_original_points')

        # region grow polygons
        # region_grow_fn = os.path.join(dir, f'{cntr}_lucas_region_grow.gpkg')
        # df_region_grow = gpd.read_file(region_grow_fn)
        df_region_grow = gpd.read_file(cntr_fn, layer='lucas_region_grow')

        # no match points
        # nomatch_points_fn = os.path.join(dir, f'{cntr}_lucas_nomatch_points.gpkg')
        # df_nomatch_points = gpd.read_file(nomatch_points_fn)
        df_nomatch_points = gpd.read_file(cntr_fn, layer='lucas_nomatch_points')

        # ANALYSIS
        # rg / no match points
        number_points = df_original_points.shape[0]
        number_procesed = df_region_grow.shape[0]
        no_processed = df_nomatch_points.shape[0]
        print(f'There are {number_procesed} points processed, but {no_processed} points did not match out of {number_points}.')

        # Multi-class cases
        number_multiclass = (df_region_grow['multiclass'] > 1).sum() # .unique() # .size
        print(f'There are multiclass {number_multiclass} points out of {number_points}.')
        print(f'It is {round(number_multiclass / number_points * 100, 2)} %')

        # Updated points
        number_updated_points = (df_region_grow['point_update'] > 0).sum()
        print(f'There are {number_updated_points} updated points out of {number_points}.')
        print(f'It is {round(number_updated_points / number_points * 100, 2)} %')

        # Class similarity
        sim_1 = (df_region_grow['similarity'] == 1).sum()
        print(f'Detected similarity ~ 1 = {sim_1} out of {number_points}.')

        # Update table
        # df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
        # df = df.append({
        # {
        #     'Country': cntr, 'Version': v_max, 'No_points': number_points,
        #     'RG': number_procesed,
        #     'RG_pct': round(number_procesed / number_points * 100, 1),
        #     'No_match': no_processed,
        #     'No_match_pct': round(no_processed / number_points * 100, 1),
        #     'Multiclass': number_multiclass,
        #     'Multiclass_pct': round(number_multiclass / number_points * 100, 1),
        #     'Geom_update': number_updated_points,
        #     'Geom_update_pct': round(number_updated_points / number_points * 100, 1),
        #     'Sim_1': sim_1, 'Sim_1_pct': round(sim_1 / number_points * 100, 1)
        # }

        if number_points > 0:
            row = {'Country': cntr, 'Version': v_max, 'No_points': number_points,
                   'RG': number_procesed, 'RG_pct': round(number_procesed / number_points * 100, 1),
                   'No_match': no_processed, 'No_match_pct': round(no_processed / number_points * 100, 1),
                   'Multiclass': number_multiclass, 'Multiclass_pct': round(number_multiclass / number_points * 100, 1),
                   'Geom_update': number_updated_points, 'Geom_update_pct': round(number_updated_points / number_points * 100, 1),
                   'Sim_1': sim_1, 'Sim_1_pct': round(sim_1 / number_points * 100, 1)}
        else:
            row = {'Country': cntr, 'Version': v_max, 'No_points': number_points,
                   'RG': number_procesed, 'RG_pct': 0,
                   'No_match': no_processed, 'No_match_pct': 0,
                   'Multiclass': number_multiclass, 'Multiclass_pct': 0,
                   'Geom_update': number_updated_points, 'Geom_update_pct': 0,
                   'Sim_1': 0, 'Sim_1_pct': 0}

        df = pd.concat([df, pd.DataFrame.from_dict([row])], ignore_index=True)

    eu_lucas_repre_stat_fn = os.path.join(dst_path, 'eu_lucas_representativeness_stat.csv')
    df.to_csv(eu_lucas_repre_stat_fn, sep=';')
    print('End.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--src_dir', type=str, help='Directory with region grow representatrive polygons for countries.')
    parser.add_argument('-d', '--dst_dir', type=str, help='Destination directory for results.')
    parser.add_argument('-v', '--version', type=str, help='Processed version.')
    parser.add_argument('-w', '--workers', metavar='workers', type=int, default=12,
                        help='Number of workers for parallel CPU processing.')

    args = parser.parse_args()

    products = ['*_lucas_region_grow.gpkg']

    main(args.src_dir, args.dst_dir, args.version, products)
