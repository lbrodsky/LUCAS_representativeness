#!/usr/bin/env python3

import os
import copy
import numpy as np
from osgeo import ogr
from osgeo import gdal
from osgeo import osr
import shapely
from shapely.geometry import Point
if int(shapely.__version__[0]) < 2:
    from shapely.wkb import loads as loads_wkb
else:
    from shapely import from_wkb as loads_wkb
import logging


def get_extent(ds):
    """Return list of corner coordinates from a gdal Dataset.
    """
    xmin, xpixel, _, ymax, _, ypixel = ds.GetGeoTransform()
    width, height = ds.RasterXSize, ds.RasterYSize
    xmax = xmin + width * xpixel
    ymin = ymax + height * ypixel

    return (xmin, ymin, xmax, ymax)


def create_polygon(extent):
    """Create polygon geometry from coordinates.
    """
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(extent[0], extent[1])
    ring.AddPoint(extent[0], extent[3])
    ring.AddPoint(extent[2], extent[3])
    ring.AddPoint(extent[2], extent[1])
    ring.AddPoint(extent[0], extent[1])
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    return poly


def define_fields():
    """Define set of vector attribues
       ['point_id', 'lc1_h', 'lc1_name', 'osm_code', 'lc',
       'gps_precision', ' gprec_src', 'lc1_codes',  ... ]
       to be written in the resulting vector files.
    """

    return {
        "point_id": ogr.FieldDefn('point_id', ogr.OFTInteger),
        'tile_id': ogr.FieldDefn('tile_id', ogr.OFTString),
        "lc1_h": ogr.FieldDefn('lc1_h', ogr.OFTString),
        "lc1_name": ogr.FieldDefn('lc1_name', ogr.OFTString),
        "osm_codes": ogr.FieldDefn('osm_codes', ogr.OFTString),
        "osm_names": ogr.FieldDefn('osm_names', ogr.OFTString),
        "osm_lc": ogr.FieldDefn('osm_lc', ogr.OFTString),
        "multiclass": ogr.FieldDefn('multiclass', ogr.OFTInteger),
        "geom_type": ogr.FieldDefn('geom_type', ogr.OFTInteger),
        "obs_type": ogr.FieldDefn('obs_type', ogr.OFTString),
        "obs_dist": ogr.FieldDefn('obs_dist', ogr.OFTInteger),
        "gps_prec":ogr.FieldDefn('gps_prec', ogr.OFTReal),
        "gprec_src": ogr.FieldDefn('gprec_src', ogr.OFTString),
        "lc": ogr.FieldDefn('lc', ogr.OFTInteger),
        "similarity": ogr.FieldDefn('similarity', ogr.OFTInteger),
        "pt_update": ogr.FieldDefn('pt_update', ogr.OFTInteger),
        "rect": ogr.FieldDefn('rect', ogr.OFTReal),
        "length": ogr.FieldDefn('length', ogr.OFTReal),
        "width": ogr.FieldDefn('width', ogr.OFTReal),
        "area": ogr.FieldDefn('area', ogr.OFTReal),
        "ratio": ogr.FieldDefn('ratio', ogr.OFTReal)
    }


def create_layer(fn, srs, field_defs):
    """Create new vector layer with attribute definitions.
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    vector_path = fn
    out_data = driver.CreateDataSource(vector_path)
    layer_name = os.path.basename(vector_path).split('.')[0]
    out_layer = out_data.CreateLayer(layer_name, srs)

    field_exclude = []
    # '_'.join(layer_name.split('_')[3:])
    if '_'.join(layer_name.split('_')[-3:]) not in ['lucas_region_grow', 'sentinel2_region_grow', 'lucas_urban_grow']:
        field_exclude = ['lc_update_def', 'gprec_src', 'unc', 'pup']

    for k, v in field_defs.items():
        if k not in field_exclude:
            out_layer.CreateField(v)

    return out_data, out_layer

def define_outputs(output_dir, tile_id, layers_kw, srs, field_defs):
    """Define output vector layers.
    """
    outputs = {}
    for layer in layers_kw:
        vector_path = os.path.join(output_dir,
                                   tile_id + '_lucas_' + layer + '.shp')
        out_rg_ds, out_rg_layer = create_layer(vector_path, srs, field_defs)
        outputs[layer] = {'ds': out_rg_ds, 'ds_layer': out_rg_layer}

    return outputs


def save_outputs(outputs):
    """Save created layers into output vectors.
    """
    for v in outputs.values():
        driver = name = None
        if v['ds_layer'].GetFeatureCount() < 1:
            driver = v['ds'].GetDriver()
            name = v['ds'].GetName()
        v['ds'] = None
        if driver:
            driver.DeleteDataSource(name)


def xy_to_array(raster_fn, x_lucas, y_lucas):
    """Find NumPy array indexes based on spatial coordinates.
    """
    ds = gdal.Open(raster_fn)
    ulx, xres, xskew, uly, yskew, yres = ds.GetGeoTransform()
    x_cell = round(((x_lucas - ulx) / xres))
    y_cell = round(((uly - y_lucas) / abs(yres)))

    return x_cell, y_cell


def get_patch_geo_transform(geo_transform, y_size, x_cell, y_cell, patch_x_size, patch_y_size):
    """Extract patch geo metadata from full image array.
    """
    ulx, xres, xskew, uly, yskew, yres = geo_transform
    lry = uly + (y_size * yres) # tile coord

    if patch_x_size % 2 > 0:
        # -1  ...
        ulx_patch = ulx + (x_cell * xres) - ((patch_x_size - 1) * xres / 2)
    else:
        ulx_patch = ulx + (x_cell * xres) - ((patch_x_size) * xres / 2)

    if patch_y_size % 2 > 0:
        uly_patch = lry + (y_cell * abs(yres)) + ((patch_y_size - 1) * abs(yres) / 2)
    else:
        uly_patch = lry + (y_cell * abs(yres)) + ((patch_y_size) * abs(yres) / 2)
    # y_cell_upper = round(((uly - y_geo) / abs(yres)))
    # uly_patch = uly - (y_cell_upper * abs(yres)) + ((patch_y_size - 1) * abs(yres) / 2)

    tile_geo_transform = ulx_patch, geo_transform[1], geo_transform[2], uly_patch, geo_transform[4], geo_transform[5]

    return tile_geo_transform


def convert_polygons2multi(layer):
    """Convert multiple single polygons into multipolygon object.
    """
    m_feat = None
    m_poly = ogr.Geometry(ogr.wkbMultiPolygon)

    while True:
        feat = layer.GetNextFeature()
        if feat is None:
            break
        if m_feat is None:
            m_feat = feat

        geometry = feat.GetGeometryRef()
        # update multipolygon geometry
        logging.debug(f'Geometry: {geometry}')
        m_poly.AddGeometry(geometry)
        # delete polygon feature
        layer.DeleteFeature(feat.GetFID())

    # add multipolygon feature
    m_feat.SetGeometry(m_poly)
    layer.CreateFeature(m_feat)

    return m_feat


def vectorize_grown_point(grown_point, out_rg_layer, geo_transform, geo_proj, output_fields, lucas_geometry, urban=False):
    """Convert the NumPy grown region into vector.
    """
    grown_point = grown_point.astype(np.byte)
    logging.debug(f"Size of RG: {np.sum(grown_point)}")

    # Save to temp  memory
    ycount, xcount = grown_point.shape
    temp_raster = gdal.GetDriverByName('MEM').Create('', xcount, ycount, 1, gdal.GDT_Byte)
    temp_raster.SetGeoTransform(geo_transform)
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(geo_proj)
    temp_raster.SetProjection(raster_srs.ExportToWkt())
    # polygonize only region = 1, not 2
    grown_point_1 = np.copy(grown_point)
    np.place(grown_point_1, (grown_point_1 == 2), [0])
    temp_raster.GetRasterBand(1).WriteArray(grown_point_1)

    temp_band = temp_raster.GetRasterBand(1)
    # vectorize temp raster
    gdal.Polygonize(temp_band, temp_band, out_rg_layer, 0, ["8CONNECTED=8"], callback=None)

    # update fields point_id, lc1_h, uncertain, pt_update
    logging.debug(
        f"Region Grow multi-polygon feature count: {out_rg_layer.GetFeatureCount()}"
    )

    ### is applied due to 8CONNECTED=8 ###
    out_rg_layer.ResetReading()
    if not urban:
        # non-urban: remove small polygonized features
        out_rg_layer.SetAttributeFilter('point_id = 1')

        rg_feat_area = {}
        if out_rg_layer.GetFeatureCount() > 1:
            rg_area_max = -1
            rg_fid_max = -1
            while True:
                rg_feat = out_rg_layer.GetNextFeature()
                if rg_feat is None:
                    break
                rg_fid = rg_feat.GetFID()
                rg_area = rg_feat.GetGeometryRef().GetArea()
                if rg_area > rg_area_max:
                    rg_area_max = rg_area
                    rg_fid_max = rg_fid
                rg_feat_area[rg_fid] = rg_area

            rg_feat = out_rg_layer.GetFeature(rg_fid_max)  # process feature with max area
        else:
            rg_feat = out_rg_layer.GetNextFeature()

        # delete small polygonized features if exists and unset filter
        for fid in rg_feat_area.keys():
            if fid != rg_fid_max:
                out_rg_layer.DeleteFeature(fid)
        out_rg_layer.SetAttributeFilter('')

    else:
        # urban: convert polygonized features into multipolygon feature
        # rg_feat = out_rg_layer.GetNextFeature()
        rg_feat = convert_polygons2multi(out_rg_layer)

    rg_geometry = rg_feat.GetGeometryRef()

    # overlay rg_geom polygon with geometry of LUCAS point
    if not lucas_geometry.Within(rg_geometry):
        output_fields["pt_update"] = 1
    else:
        output_fields["pt_update"] = 0

    # calculate polygon width and length
    output_fields["length"], output_fields["width"] = get_length_width(rg_geometry)
    output_fields["area"] = rg_geometry.GetArea()
    output_fields["ratio"] = output_fields["width"] / output_fields["length"]

    for field, value in output_fields.items():
        rg_feat.SetField(field, value)

    out_rg_layer.SetFeature(rg_feat)

    # perform shape generalization
    # rg_geometry = rg_geometry.Buffer(-6, options=["JOIN_STYLE=MITRE"]).Buffer(6)
    rg_geometry = rg_geometry.Buffer(-6).Buffer(5.9)
    rg_feat.SetGeometry(rg_geometry)
    out_rg_layer.SetFeature(rg_feat)

    point_id = rg_feat.GetField("point_id")
    temp_raster.AddBand(gdal.GDT_Byte)
    temp_raster.GetRasterBand(2).WriteArray(np.zeros(grown_point.shape))
    arr = temp_raster.GetRasterBand(2).ReadAsArray()
    gdal.RasterizeLayer(temp_raster, [2], out_rg_layer, burn_values=[1],
                        options=[f"where='point_id = {point_id}'", "ALL_TOUCHED=TRUE"])
    arr = temp_raster.GetRasterBand(2).ReadAsArray()
    temp_driver = ogr.GetDriverByName('Memory')
    temp_ds = temp_driver.CreateDataSource('memory')
    temp_layer = temp_ds.CreateLayer('temp', geom_type=ogr.wkbPolygon)
    temp_layer.CreateField(ogr.FieldDefn('Value', ogr.OFTInteger))
    temp_band2 = temp_raster.GetRasterBand(2)
    gdal.Polygonize(temp_band2, temp_band2, temp_layer, 0, [], callback=None)
    temp_feat = temp_layer.GetNextFeature()
    rg_feat.SetGeometry(temp_feat.GetGeometryRef())
    out_rg_layer.SetFeature(rg_feat)

def get_length_width(geometry):
    """Calculate length and with from region grow
    """
    # logging.debug(f'Geometry export: {bytes(geometry.ExportToIsoWkb())}')
    poly = loads_wkb(bytes(geometry.ExportToIsoWkb())) 
    # get minimum bounding box around polygon
    box = poly.minimum_rotated_rectangle
    # get coordinates of polygon vertices
    x, y = box.exterior.coords.xy
    # get length of bounding box edges
    edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))
    # get length of polygon as the longest edge of the bounding box
    length = max(edge_length)
    # get width of polygon as the shortest edge of the bounding box
    width = min(edge_length)

    return length, width


def array_to_xy(geo_transform, x_idx, y_idx, pixel_size):
    """Find geo-coordinates based on array indexes.
    """
    ulx, xres, xskew, uly, yskew, yres = geo_transform
    x_coord = ulx + (x_idx * xres) + (pixel_size / 2.)
    y_coord = uly - (y_idx * abs(yres)) - (pixel_size / 2.)

    return x_coord, y_coord


def save_point(lucas_geometry, out_point_layer, output_fields):
    """Save vector points geometry and attributes into layer.
    """
    # save LUCAS original geometry and LC attribute
    layer_defn = out_point_layer.GetLayerDefn()
    point_feature = ogr.Feature(layer_defn)
    for field, value in output_fields.items():
        if layer_defn.GetFieldIndex(field) > -1: # skip exluded fields
            point_feature.SetField(field, value)
    point_feature.SetGeometry(lucas_geometry)
    # out_point_layer.SetFeature(point_feature)
    out_point_layer.CreateFeature(point_feature)
