#!/usr/bin/env python3

import os
import sys
import time
import joblib 
from joblib import Parallel, delayed 
import argparse
import glob
import math
import logging
import numpy as np
from osgeo import ogr
from osgeo import osr
from shapely.geometry import Polygon
from osgeo import gdal
gdal.UseExceptions() 
import multiprocessing

from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from region_grow.read_data import *
from region_grow.data_utils import *
from region_grow.region_grow import Point, RegionGrow
from region_grow.exceptions import ConfigError, IllegalArgumentError

# constants
RAD2DEGREE = 180 / math.pi
DEGREE2RAD = math.pi / 180

__version__ = "1.0"

def init_logging(dst_dir, args):
    """Initialize the processing log.
    """
    log_file = os.path.join(
        dst_dir, 'log_rg_process_thr_' + str(args.shp_thr) + '_tiles.txt'
    )
    log_level = {'info': logging.INFO, 'debug': logging.DEBUG}[args.log_level]

    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',
                        level=log_level,
                        handlers=[logging.FileHandler(log_file, mode='w'),
                                  logging.StreamHandler()])

    # log parameter settings
    logging.info('---')
    logging.info('Parameter settings')
    logging.info('---')
    logging.info(f'Process will be logged into: {log_file}')
    for arg, value in vars(args).items():
        logging.info(f"{arg}: {value}")

def prepare_inputs(args):
    """Prepare variables for input tiles, point filter and translation tables.
    """
    if args.selected_tiles is not None:
        tiles = []
        tiles_to_process = args.selected_tiles.split(',')

        for tile in tiles_to_process:
            tile_ = glob.glob(args.tiles_dir + f'/*{tile}.tif')
            if tile_ != []:
                tiles.append(tile_[0])
    else:
        tiles = glob.glob(args.tiles_dir + '/*.tif')

    if len(tiles) == 0:
        raise ConfigError(f'No available tiles to process in directory {args.tiles_dir}')

    if args.selected_points is None:
        selected_points = []
    else:
        selected_points = args.selected_points.split(',')

    tables = {
        'similarity_fn': os.path.join(args.tables_dir, 'semsim_l3.csv'),
        'lc1_2_lc_fn': os.path.join(args.tables_dir, 'lucas_lc1_to_lc.yaml'),
        'lc1_2_osm_fn': os.path.join(args.tables_dir, 'lucas_lc1_to_osm.yaml'),
        'lc1_codes_2_names_fn': os.path.join(args.tables_dir, 'lc1_h.csv'),
        'osm_2_lc_fn': os.path.join(args.tables_dir, 'osm_to_lc.yaml'),
        'osm_names_fn': os.path.join(args.tables_dir, 'osm_names.yaml'),
        'lc_codes_names_fn': os.path.join(args.tables_dir, 'lc_codes_names.yaml')
    }

    return tiles, selected_points, tables


def check_inputs(args):
    """Check if OSM tile and LUCAS points are valid input files
    """
    # OSM land cover tiles
    if not os.path.isdir(args.tiles_dir):
        raise ConfigError(f'Tiles directory {args.tiles_dir} does not exist')

    # LUCAS points
    if not os.path.isfile(args.lucas_points):
        raise ConfigError(f'LUCAS points file {args.lucas_points} not available.')


def select_LUCAS_points(img_ds, layer):
    """Apply spatial filter on LUCAS points.
    """
    # tile extent (xmin, ymin, xmax, ymax)
    extent = get_extent(img_ds)
    logging.debug(f'Extent is: {extent}')
    # spatial filter to current tile
    tile_ext_geom = create_polygon(extent)
    layer.SetSpatialFilter(tile_ext_geom)

    return layer.GetFeatureCount()

def polar_to_cartesian(centre, angle, radius):
    """Convert polar coords to cartesian.
    """
    x = math.cos(angle) * radius + centre[0]
    y = math.sin(angle) * radius + centre[1]
    return (x,y)

def make_arc(centre, radius, rad_from_azimuth, rad_to_azimuth, rad_angle_resolution):
    """Create arc for wedge buffer.
    """
    cartesian = []
    rad_az = rad_from_azimuth
    if rad_from_azimuth < rad_to_azimuth:
        while rad_az < rad_to_azimuth:
            cartesian.append(polar_to_cartesian(centre, rad_az, radius))
            rad_az = rad_az + rad_angle_resolution
    else:
        while rad_az > rad_to_azimuth:
            cartesian.append(polar_to_cartesian(centre, rad_az, radius))
            rad_az = rad_az - rad_angle_resolution
    cartesian.append(polar_to_cartesian(centre, rad_to_azimuth, radius))
    return cartesian

def wedge_buffer(centre, radius, azimuth, opening_angle, inner_radius=0, angle_resolution=10):
    """Create wedge buffer.
    """
    # make azimuth 0 north
    azimuth += 90
    rad_from_azimuth = (azimuth - opening_angle * 0.5) * DEGREE2RAD
    rad_to_azimuth = (azimuth + opening_angle * 0.5) * DEGREE2RAD
    rad_angle_res = angle_resolution * DEGREE2RAD

    cartesian_coords = make_arc(centre, radius, rad_from_azimuth, rad_to_azimuth, rad_angle_res)

    if inner_radius <= 0:
        cartesian_coords.append(centre)
    else:
        # reverse arc at inner radius
        cartesian_coords += make_arc(centre, inner_radius, rad_to_azimuth, rad_from_azimuth, rad_angle_res)

    # Close ring
    cartesian_coords.append(cartesian_coords[0])
    return cartesian_coords

def calculate_azimuth(ax, ay, bx, by):
    """Computes the bearing in degrees

    dx=Math.abs(x1-x2) and dy=Math.abs(y1-y2)
    tan(theta) = dx/dy

    :param ax (int): The x-coordinate of the first point defining a line.
    :param ay (int): The y-coordinate of the first point defining a line.
    :param bx (int): The x-coordinate of the second point defining a line.
    :param by (int): The y-coordinate of the second point defining a line.

    :return float: bearing in degrees
    """
    two_pi = math.pi * 2
    theta = math.atan2(bx - ax, ay - by)
    if theta < 0.0:
        theta += two_pi
    return math.degrees(theta)


def create_buffer_mask(geometry, geometry2, gps_prec_val, out_buffer_layer, output_fields,
                       tile_shape, srs, geo_transform):
    """Create buffer around LUCAS point geometry and convert it into raster / NumPy array.
    """
    x1 = geometry.GetX()
    y1 = geometry.GetY()
    x2 = geometry2.GetX()
    y2 = geometry2.GetY()

    dist = math.dist([x1, y1], [x2, y2])

    if dist >= gps_prec_val:
        # create wedge vector
        azimuth = calculate_azimuth(x2, y2, x1, y1)
        opening_angle = 120
        poly = wedge_buffer((x1, y1), gps_prec_val, azimuth, opening_angle)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for p in poly:
            ring.AddPoint_2D(p[0], p[1])
        # create polygon
        point_buffer = ogr.Geometry(ogr.wkbPolygon)
        point_buffer.AddGeometry(ring)
    else:
        # create point buffer
        point_buffer = geometry.Buffer(gps_prec_val)

    buf_feature = ogr.Feature(out_buffer_layer.GetLayerDefn())
    buf_layer_defn = out_buffer_layer.GetLayerDefn()
    for field, value in output_fields.items():
        if buf_layer_defn.GetFieldIndex(field) > -1: # skip exluded fields
            buf_feature.SetField(field, value)
    buf_feature.SetGeometry(point_buffer)
    out_buffer_layer.CreateFeature(buf_feature)

    # save to memory
    drv = ogr.GetDriverByName('Memory')
    dst_ds = drv.CreateDataSource('')
    dst_layer = dst_ds.CreateLayer('', srs=srs, geom_type=ogr.wkbPolygon)
    label_field = ogr.FieldDefn("label", ogr.OFTInteger)
    dst_layer.CreateField(label_field)
    feature_defn = dst_layer.GetLayerDefn()
    buffer_feature = ogr.Feature(feature_defn)
    buffer_feature.SetGeometry(point_buffer)
    buffer_feature.SetField("label", 1)
    dst_layer.CreateFeature(buffer_feature)
    buffer_feature = None

    # create array mask
    ycount, xcount = tile_shape
    target_ds = gdal.GetDriverByName('MEM').Create('', xcount, ycount, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform(geo_transform)
    target_ds.SetProjection(srs.ExportToWkt())
    gdal.RasterizeLayer(target_ds, [1], dst_layer, burn_values=[1])
    band_mask = target_ds.GetRasterBand(1)
    dst_ds.Close()

    return band_mask.ReadAsArray()


def get_patches(tile_img, buffer_mask):
    """Mask input OSM tile based on buffer mask.
    """
    return tile_img * buffer_mask


def get_osm_lc_codes_from_buffer(lc_masked):
    """Get unique land cover classes, except 0, no data.
    """
    unique_values = np.unique(lc_masked)
    return list(np.delete(unique_values, 0))


def get_similarity(osm_lc_code, lc1, lc_mappings):
    """Compare OSM and LUCAS land cover similarities.
    """
    row_select = lc_mappings['sim_tab'].loc[lc_mappings['sim_tab']['LUCAS_L2_codes'] == lc1] 
    sim_val = row_select[str(osm_lc_code)].item()
    logging.debug(f'Similarity value: {sim_val}')

    return sim_val


def get_patch_area(patches, osm_lc_code, lc_mappings):
    """Measure the area of the masked patch.
    """
    osm_code = 0
    unique_patches_codes = list(np.unique(patches))
    for k, v in lc_mappings['osm_2_lc'].items():
        if v == osm_lc_code:
            if k in unique_patches_codes:
                osm_code = k
    patch = patches == osm_code
    return np.sum(patch)


def update_seed_point(patches, osm_code_seed):
    """Update geometry of the LUCAS point, the seed point for region grow.
    """
    lc_masked_correct_px = (patches == osm_code_seed) * 1
    count = int(np.sum(lc_masked_correct_px))
    logging.debug(f'Sum of mask points: {count}')
    # https://stackoverflow.com/questions/38933566/quickly-calculating-centroid-of-binary-numpy-array
    # get the coordinates of the center point
    yx = np.argwhere(lc_masked_correct_px == 1).sum(0) / count
    y_center, x_center = None, None

    if yx[0] > 0 and yx[1] > 0:
        y_center, x_center = list(map(lambda x: int(x + 0.5), yx))
        logging.debug(f'Center point: {y_center} {x_center}')

    seed_value = patches[y_center, x_center]
    if seed_value == osm_code_seed:
        logging.debug('Centroid correct')
    else:
        logging.debug('Seed value mis-match')
        neighbours_ne = (
            (-1, 0), (-1, 1), (0, 1), (1, 1)
        )
        neighbours_sw = (
            (1, 0), (1, -1), (0, -1), (-1, -1)
        )
        directions = [neighbours_ne, neighbours_sw]
        for neighbours in directions:
            multiplier = 1
            max_multiplier = 4
            current_point_x = x_center
            current_point_y = y_center

            while multiplier <= max_multiplier and (seed_value != osm_code_seed):
                for y, x in neighbours:
                    x_center = current_point_x + (x * multiplier)
                    y_center = current_point_y + (y * multiplier)

                    # check position within tile patch
                    if x_center < 0 or y_center < 0 or x_center >= patches.shape[1] or y_center >= patches.shape[0]:
                        continue

                    seed_value = patches[y_center, x_center]
                    if seed_value == osm_code_seed:
                        multiplier = max_multiplier
                        logging.debug(f'Center point: {y_center} {x_center}')
                        break

                multiplier += 1
            if seed_value == osm_code_seed:
                break

        if seed_value != osm_code_seed:
            logging.debug('No centroid found')
            x_center, y_center = None, None


    return x_center, y_center


def get_tile_patch(full_img, lucas_point, tile_size, pixel_size):
    """Extract patch from full image array.
    """
    y_up = lucas_point.y - round(tile_size/2/pixel_size)
    y_down = lucas_point.y + round(tile_size/2/pixel_size)
    x_left = lucas_point.x - round(tile_size/2/pixel_size)
    x_right = lucas_point.x + round(tile_size/2/pixel_size)

    if y_up < 0:
        # if index out of tile, shift lower index to get symetrical square patch
        y_up_diff = abs(y_up)
        y_down = y_down - (y_up_diff -1)
        y_up = 0
        logging.warning('Upper Y is out of the image - using 0 instead')
    if y_down >= full_img.shape[0]:
        y_down = full_img.shape[0] - 1
        logging.warning(
            'Lower Y is out of the image - using the lowest coordinate instead'
        )
    if x_left < 0:
        # treat the index out of tile
        x_left_diff = abs(x_left)
        x_right = x_right + (x_left_diff -1)
        x_left = 0
        logging.warning('Left X is out of the image - using 0 instead')
    if x_right >= full_img.shape[1]:
        x_right = full_img.shape[1] - 1
        logging.warning(
            'Right X is out of the image - using the rightmost coordinate  instead'
        )

    # get the patch
    return full_img[y_up:y_down, x_left:x_right]


def get_grown_region(x_array, y_array, tile_img, region_max_size, pixel_size, shp_thr, threshold, connectivity):
    """Run region grow process for given x, y coordinates and return grown area as NumPy array.
    """
    seed_point = Point(x_array, y_array)
    patch = get_tile_patch(tile_img, seed_point, region_max_size, pixel_size)
    tile_center_y = int(patch.shape[0] / 2)
    tile_center_x = int(patch.shape[1] / 2)
    tile_center_point = [Point(tile_center_x, tile_center_y)]
    rg = RegionGrow(patch, tile_center_point, shp_thr, threshold, connectivity)
    grown_point, rect = rg.grow()

    return patch, grown_point, rect


def get_repre_data(t, geometry, geometry2, gps_prec_val, outputs, output_fields, tile_img, srs, geo_transform, point_id,
                   lc_mappings, lc1, region_max_size, pixel_size, shp_thr, threshold, connectivity,
                   obs_direct, max_multiplier):
    """Calculate representativess and related data.
    """
    patch, grown_point, rectangularity = None, None, None
    x_array, y_array = None, None
    similarity = None
    osm_code_seed = None
    x_lucas = geometry.GetX()
    y_lucas = geometry.GetY()
    # create buffer around LUCAS point based on GPS presision
    buffer_mask = create_buffer_mask(geometry, geometry2, gps_prec_val, outputs['points_buffer']['ds_layer'], output_fields,
                                     tile_img.shape, srs, geo_transform)

    # read OSM codes (10, 20, ..., 100) within the buffer area
    patches = get_patches(tile_img, buffer_mask)
    osm_codes = get_osm_lc_codes_from_buffer(patches)
    if osm_codes == []:
        logging.info(f'No OSM land cover over LUCAS point ID: {point_id}!')
        save_point(geometry, outputs['original_points']['ds_layer'], output_fields)
        return None

    # translate OSM to LUCAS lc coding
    lucas_osm_code = lc_mappings['lucas_lc1_2_osm'][lc1]

    # setting: coding for LUCAS lc and OSM land cover comparison - similarity
    # 0 .. no similarity -> different classes
    # 1 .. partial similarity
    # 2 .. full match

    # init
    grown_point = None
    # pt_update = 0 pt_update is solved in data_utils.py
    prec_multiplier = 1
    max_similarity = 0
    osm_lc_l1 = None
    rectangularity = 0

    if len(osm_codes) == 1 and osm_codes[0] == lucas_osm_code:
        # one code within buffer area -> run directly RG on original LUCAS geometry
        max_similarity = 2
        # get array index of seed point based on LUCAS geometry
        x_array, y_array = xy_to_array(t, geometry.GetX(), geometry.GetY())
        # patch is a symmetrical subset (e.g. 100m) of tile around LUCAS point
        # grown_point is a result of RG algorithm, Numpy array with coding 1..representative region, 0..not region
        patch, grown_point, rectangularity = get_grown_region(x_array, y_array, tile_img, region_max_size, pixel_size,
                                                              shp_thr, threshold, connectivity)

        osm_lc_l1 = lc_mappings['osm_2_lc'][osm_codes[0]]

    else:
        # multiple codes within buffer area
        # -> run RG on correct patch of LUCAS buffer or increse the buffer size
        while True:
            max_similarity = 0
            max_area = 0
            osm_code_seed = None
            i = 0

            for osm_code in osm_codes:
                similarity = get_similarity(osm_code, lc1, lc_mappings)
                if similarity == 2:
                    osm_code_seed = osm_codes[i]
                    max_similarity = similarity
                    max_area = None
                    break
                elif similarity == 1:
                    # returns number of pixels
                    area = get_patch_area(patches, osm_code, lc_mappings)
                    if area > max_area:
                        max_area = area
                        osm_code_seed = osm_codes[i]
                        max_similarity = similarity
                i += 1

            # if first round and not full match and obs direct N or E
            # then point falls on the edge of two classes (linear
            # feature, etc.)  see LUCAS 2018 doc:
            # https://ec.europa.eu/eurostat/documents/205002/8072634/LUCAS-2018-C1-Instructions.pdf
            # page 25
            if prec_multiplier == 1 and similarity != 2 and obs_direct != 1:
                osm_code_seed = None

            if osm_code_seed is not None:
                # if 1 or 2, update seed point
                # update point geometry cannot be decided here!
                # pt_update = 1
                x_array, y_array = update_seed_point(patches, osm_code_seed)

                if x_array is not None:
                    try:
                        patch, grown_point, rectangularity = get_grown_region(x_array, y_array, tile_img,
                                                                              region_max_size, pixel_size,
                                                                              shp_thr, threshold, connectivity)

                    except IndexError:
                        raise IndexError(f'RG error at point {point_id}')

                    # update the x, y geo coordinates based on shifted patch
                    x_lucas, y_lucas = array_to_xy(geo_transform, x_array, y_array, pixel_size)
                    osm_lc_l1 = lc_mappings['osm_2_lc'][osm_code_seed]
                    break

            # stop enlarging buffer when max_multiplier reached
            if prec_multiplier == max_multiplier:
                break

            # continue with updated buffer size
            prec_multiplier += 1
            buffer_mask = create_buffer_mask(geometry, geometry2, prec_multiplier * gps_prec_val,
                                             outputs['points_buffer']['ds_layer'],
                                             output_fields,
                                             tile_img.shape, srs, geo_transform)
            patches = get_patches(tile_img, buffer_mask)
            osm_codes = get_osm_lc_codes_from_buffer(patches)

    return {'osm_codes': osm_codes,
            'osm_lc_l1': osm_lc_l1,
            'max_similarity': max_similarity,
            'osm_lc': osm_code_seed,
            'rectangularity': rectangularity,
            'grown_point': grown_point,
            'patch': patch,
            'x_array': x_array, 'y_array': y_array,
            'prec_multiplier': prec_multiplier,
            'x_lucas': x_lucas, 'y_lucas': y_lucas
            }


def process_single_tile(config):
    """Main region grow process at level of single tile.
    """
    # Unpack configuration parameters
    t = config[0]
    lucas_fn = config[1]
    lucas_thr_fn = config[2]
    output_dir = config[3]
    selected_points = config[4]
    lc_mappings = config[5]
    region_max_size = config[6].region_max_size
    threshold = 1 # OSM land cover is coded with integers
    shp_thr = config[6].shp_thr
    connectivity = config[6].connectivity
    max_multiplier = config[6].multiplier
    shp_generalize_dist = config[6].shp_generalize_dist

    # read OSM raster data
    img_ds = gdal.Open(t)
    pixel_size = img_ds.GetGeoTransform()[1]
    geo_transform = img_ds.GetGeoTransform()
    tile_img = img_ds.ReadAsArray()
    # open ST_LUCAS point layer
    lucas_ds, lucas_layer = open_LUCAS_points(lucas_fn) 
    feature_count = select_LUCAS_points(img_ds, lucas_layer)
    # open ST_LUCAS theoretical point layer
    lucas_thr_ds, lucas_thr_layer = open_LUCAS_points(lucas_thr_fn)
    logging.info(f"Selected features from {os.path.basename(t)}: {feature_count}")
    logging.info('---')

    # define outputs
    srs = osr.SpatialReference()
    srs.ImportFromWkt(img_ds.GetProjectionRef())
    tile_id = os.path.basename(t).split('.')[0]
    outputs = define_outputs(output_dir, tile_id, ['region_grow', 'nomatch_points', 'points_buffer',
                                                   'original_points', 'updated_points'],
                             srs, define_fields())

    # cycle over selected LUCAS points and run Region Grow
    feature = lucas_layer.GetNextFeature()
    while feature:
        point_id = str(int(feature.GetField('point_id')))
        lc1 = feature.GetField('lc1_h')
        # check if lc is a real class
        if lc1 == '-1':
            feature = lucas_layer.GetNextFeature()
            continue
        # dev filter
        if selected_points and point_id not in selected_points:
            feature = lucas_layer.GetNextFeature()
            continue
        logging.debug(f'Processing PID/tile: {point_id, t}') 
        # get GPS accuracey from LUCAS in-situ measurements -> buffer
        gps_prec_val = feature.GetField('gps_prec')
        if gps_prec_val in (-1, 0):
            # parameter calculated as median of LUCAS GPS precision
            gps_prec_val = 5.0

        # get obs_dist: difference between GPS and theoretical point
        obs_dist = feature.GetField('obs_dist')
        logging.debug(f'GPS prec: {gps_prec_val} | obs_dist: {obs_dist}')
        obs_type = feature.GetField('obs_type')
        logging.debug(f'obs_type: {obs_type}')
        obs_direct = feature.GetField('obs_direct')

        # output vector attributes
        output_fields = {
            'point_id': int(point_id),
            'lc1_h': lc1,
            'lc1_name': lc_mappings['lucas_lc1_codes_2_names'][lc1],
            'gps_prec': gps_prec_val,
            'tile_id': os.path.basename(t)
        }

        # get LUCAS  geometry: 0 .. none, 1 .. measured, 2 .. theoretical
        geometry_type = 0
        # type insitu, office
        gps_prec_val_updated = max(gps_prec_val, geo_transform[1])
        # obs_type: 1 In situ < 100 m | obs_dist in meters,
        # obs_type: 7 In Office PI
        # out: and obs_dist < 2 * gps_prec_val
        repre_data = None
        # measured geometry
        geometry = feature.GetGeometryRef()
        # theoretical geometry
        lucas_thr_layer.ResetReading()
        lucas_thr_layer.SetAttributeFilter(f"point_id = {point_id}")
        thr_feature = lucas_thr_layer.GetNextFeature()
        lucas_thr_layer.SetAttributeFilter("")
        geometry_thr = thr_feature.GetGeometryRef()

        if obs_type == 1 and (obs_dist is not None and obs_dist < 100) and geometry is not None:
            # take GPS measured
            geometry_type = 1
            # geometry = feature.GetGeometryRef()
            repre_data = get_repre_data(t, geometry, geometry_thr, gps_prec_val_updated, outputs, output_fields, tile_img, srs, geo_transform, point_id,
                   lc_mappings, lc1, region_max_size, pixel_size, shp_thr, threshold, connectivity,
                   obs_direct, max_multiplier)

        # obs type insitu but within GPS distances 1, large distance 2, insitu PI 3
        elif (obs_dist is not None and obs_dist <= 2000.0) or (obs_type == 7): 
            # take theoretical point geometry
            geometry_type = 2
            logging.debug(f'Reading geometry from theoretical points, ID: {thr_feature.GetField("point_id")} ')
            # theoretical coords fall on edge of pixels -> increase gps prec val
            repre_data = get_repre_data(t, geometry_thr, geometry, gps_prec_val_updated, outputs, output_fields, tile_img, srs, geo_transform, point_id, lc_mappings, lc1, region_max_size, pixel_size, shp_thr, threshold, connectivity,obs_direct, max_multiplier)
            # obs_type 3: In situ PI
            if obs_type == 3 and geometry is not None and (repre_data is not None and repre_data['max_similarity'] != 2):
                repre_data_tmp = repre_data
                repre_data = get_repre_data(t, geometry, geometry_thr, gps_prec_val_updated, outputs, output_fields, tile_img, srs,
                                            geo_transform, point_id,
                                            lc_mappings, lc1, region_max_size, pixel_size, shp_thr, threshold,
                                            connectivity,
                                            obs_direct, max_multiplier)
                geometry_type = 1
                if (repre_data is None) or (repre_data is not None and repre_data_tmp is not None and repre_data['max_similarity'] <= repre_data_tmp['max_similarity']):
                    geometry_type = 2
                    repre_data = repre_data_tmp
        else: 
            logging.info('No match on measured point. Thr point too far.') 

        if repre_data is None:
            geometry_type = 0
            repre_data = {
                    'osm_codes': [], 
                    'osm_names': [], 
                    'osm_lc': '', 
                    'osm_lc_l1': lc1,
                    'max_similarity': 0,
                    'rectangularity': 0, 
                    'grown_point': None 
                    }

        output_fields["osm_codes"] = ', '.join(str(e) for e in repre_data['osm_codes'])
        output_fields["osm_names"] = ', '.join(str(e) for e in [lc_mappings['osm_names'][i] for i in repre_data['osm_codes']])
        output_fields["multiclass"] = len(repre_data['osm_codes'])
        output_fields["obs_type"] = obs_type
        output_fields["obs_dist"] = obs_dist
        output_fields["lc"] = lc_mappings['lucas_lc1_2_lc'][lc1[0:-1]]
        output_fields["similarity"] = repre_data['max_similarity']
        output_fields["rectangularity"] = round(repre_data['rectangularity'], 3)

        # save outputs
        if repre_data is not None and  repre_data['grown_point'] is not None:
            patch_y_size, patch_x_size = repre_data['patch'].shape
            patch_geo_transform = get_patch_geo_transform(geo_transform, img_ds.RasterYSize,
                                                          repre_data['x_array'], img_ds.RasterYSize - repre_data['y_array'],
                                                          patch_x_size, patch_y_size)
            # save vectors based on grown_point numpy array and return the RG geometry
            vectorize_grown_point(repre_data['grown_point'], outputs['region_grow']['ds_layer'], patch_geo_transform,
                                  img_ds.GetProjectionRef(), output_fields, geometry, shp_generalize_dist)

            logging.info(f'POINT ID: {point_id} | multiplier: {repre_data["prec_multiplier"]} | buffer: {repre_data["prec_multiplier"] * gps_prec_val}'
                         f' | tile: {os.path.basename(t)}')
        else:
            output_fields["point_update"] = 0
            save_point(geometry, outputs['nomatch_points']['ds_layer'], output_fields)
            logging.info(f'POINT ID: {point_id} | nomatch')

        save_point(geometry, outputs['original_points']['ds_layer'], output_fields)
        # save LUCAS points with updated geo coordinates
        if output_fields["point_update"]:
            updated_geometry = ogr.Geometry(ogr.wkbPoint)
            updated_geometry.AddPoint_2D(repre_data['x_lucas'], repre_data['y_lucas'])
            save_point(updated_geometry, outputs['updated_points']['ds_layer'], output_fields)

        logging.debug(f'Fields: {output_fields}')

        feature = lucas_layer.GetNextFeature()

    # close and save outputs
    save_outputs(outputs)
    lucas_ds = None
    lucas_thr_ds = None
    time.sleep(0.1) 


def process_tiles(args):
    """Main process of spatial representativness
    with region Grow algorithm for all given tiles.
    """
    start = time.time()
    if (args.tables_dir is None):
        args.tables_dir = str(Path(__file__).parent / "tables")

    dst_dir = os.path.join(args.dst_dir, 'v' + str(args.version))
    if not os.path.isdir(dst_dir):
        os.makedirs(dst_dir)

    init_logging(dst_dir, args)
    check_inputs(args)
    tiles, selected_points, tables = prepare_inputs(args)

    logging.info('---')
    logging.info('Region Grow Process')
    logging.info('---')

    # translation tables
    lc_mappings = {
        'sim_tab': read_semantic_similarity(tables['similarity_fn']),
        'lucas_lc1_2_osm': read_yaml(tables['lc1_2_osm_fn']),
        'lucas_lc1_2_lc': read_yaml(tables['lc1_2_lc_fn']),
        'lucas_lc1_codes_2_names': csv2dict(tables['lc1_codes_2_names_fn'], k='code', v='name', delim=';'),
        'osm_2_lc': read_yaml(tables['osm_2_lc_fn']),
        'osm_names': read_yaml(tables['osm_names_fn']),
        'lc_codes_names': read_yaml(tables['lc_codes_names_fn'])
    }
    # process info
    logging.info(f'Number of tiles to process: {len(tiles)}')
    logging.info(f'Number of CPUs: {multiprocessing.cpu_count()}')
    logging.info(f'Number of workers to run: {args.workers}')
    logging.info('---')

    if args.workers > 1:
        # joblib parallelization
        logging.info(f'Number of tiles to process: {len(tiles)}')
        Parallel(n_jobs=args.workers, backend="threading", verbose=10)(
            delayed(process_single_tile)([t, args.lucas_points, args.lucas_thr_points, dst_dir, selected_points, lc_mappings, args]) for t in tiles)
    else:
        # process using one core
        for t in tiles:
            process_single_tile([t, args.lucas_points, args.lucas_thr_points, dst_dir, selected_points, lc_mappings, args])

    # closing the input vector data
    end = time.time()
    logging.info(f"Process duration: {round((end - start) / 60., 2)} minutes.")


def define_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tiles_dir',
                        type=str, required=True,
                        help='Directory with OSM rasterized tiles.')
    parser.add_argument('--lucas_points',
                        type=str, required=True,
                        help='LUCAS points filename.')
    parser.add_argument('--lucas_thr_points',
                        type=str, required=True,
                        help='LUCAS theoretical points.')
    parser.add_argument('--selected_points',
                        type=str,
                        help='List of selected points (separated by comma) to be processed.')
    parser.add_argument('--selected_tiles',
                        type=str,
                        help='List of selected tiles (separated by comma) to be processed.')
    parser.add_argument('--shp_thr',
                        type=float, default=0.5,
                        help='Region grow shape threshold (default: %(default)s).')
    parser.add_argument('--region_max_size',
                        type=int, default=100,
                        help='Maximum size of region (default: %(default)s).')
    parser.add_argument('--tables_dir',
                        type=str,
                        help='Directory of translation and similarity tables (default: tables).')
    parser.add_argument('--version',
                        type=int, required=True,
                        help='Version of the RG run.')
    parser.add_argument('--dst_dir',
                        type=str, required=True,
                        help='Target directory where the results are stored.')
    parser.add_argument('--workers',
                        type=int, default=8,
                        help='Number of workers for parallel CPU processing (default: %(default)s).')
    parser.add_argument('--connectivity',
                        type=int, default=8,
                        help='Select 4 or 8 neigbouring pixels (default: %(default)s).')
    parser.add_argument('--multiplier',
                        type=int, default=3,
                        help='Set max multiplier for extending buffer zone around point (default: %(default)s).')
    parser.add_argument('--shp_generalize_dist',
                        type=int, default=6,
                        help='Set distance for shape generalization (<= 0 to disable) (default: %(default)s).')
    parser.add_argument('--log_level',
                        type=str, choices=('info', 'debug'),
                        default='info', help='Logging level (default: %(default)s).')

    return parser

if __name__ == "__main__":
    parser = define_parser()

    args = parser.parse_args()

    process_tiles(args)
