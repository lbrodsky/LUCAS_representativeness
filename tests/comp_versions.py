import argparse
from pathlib import Path

from osgeo import ogr

def main(dir_left, dir_right):
    count = 0
    mismatch = 0
    for l_tile in Path(dir_left).glob("*.gpkg"):
        print(f"Processing {l_tile}...")
        count += 1
        
        # search for right tile
        r_tile = Path(dir_right, l_tile.name.replace('lucas_representativeness', 'osm_clcplus'))
        print(r_tile)
        if not Path(r_tile).exists():
            print("ERROR: Right tile not found!")
            mismatch +=1
            continue

        if comp_tiles(l_tile, r_tile) is False:
            mismatch += 1

    print(f"Tiles processed: {count}")
    print(f"Mismatch: {mismatch}")

def comp_tiles(l_tile, r_tile, epsilon=1e-12):
    l_ds, l_layer = open_layer(str(l_tile))
    r_ds, r_layer = open_layer(str(r_tile))

    success = True

    if l_layer is None and r_layer is None:
        return success

    if (l_layer is None and r_layer is not None) or (l_layer is not None and r_layer is None):
        print(f"ERROR: Layer mismatch (left:{l_layer}, {r_ds.GetName()}:{r_layer}")
        success = False
        
    if l_layer.GetFeatureCount() != r_layer.GetFeatureCount():
        print(f"ERROR: Feature count mismatch ({l_layer.GetFeatureCount()}x{r_layer.GetFeatureCount()}: {r_ds.GetName()})!")
        success = False
    
    for l_feat in l_layer:
        point_id = l_feat.GetField('point_id')
        r_layer.SetAttributeFilter(f"point_id = {point_id}")
        r_layer.ResetReading()
        r_feat = r_layer.GetNextFeature()

        # check attributes
        for i in range(l_feat.GetFieldCount()):
            l_field_name = l_feat.GetFieldDefnRef(i).GetName()
            l_field_value = l_feat.GetField(i)
            r_field_value = r_feat.GetField(l_field_name)

            l_field_defn = l_layer.GetLayerDefn().GetFieldDefn(i)
            l_field_type = l_field_defn.GetType()
            r_field_defn = r_layer.GetLayerDefn().GetFieldDefn(r_layer.GetLayerDefn().GetFieldIndex(l_field_name))
            r_field_type = r_field_defn.GetType()

            if l_field_type != r_field_type: # avoid field type mismatch
                mismatch = str(l_field_value) != str(r_field_value)
            else:
                if l_field_type == ogr.OFTReal:
                    mismatch = abs(l_field_value - r_field_value) > epsilon
                else:
                    mismatch = l_field_value != r_field_value

            if mismatch is True:
                print(f"ERROR: Field ({l_field_name}) mismatch: {l_field_value} vs {r_field_value} [{r_ds.GetName()}: {point_id}]!")
                success = False

            # check geometry
            l_geom = l_feat.GetGeometryRef()
            r_geom = r_feat.GetGeometryRef()
            if abs(l_geom.GetArea() - r_geom.GetArea()) > 1e-12:
                print(f"ERROR: Point id {point_id} area mismatch!")
                success = False

            if l_geom.Equal(r_geom) is False:
                print(f"ERROR: Point id {point_id} geom mismatch!")
                success = False

        # if success is False:
        #     break

    l_ds = r_ds = None

    return success

def open_layer(tile):
    ds = ogr.Open(tile)
    return ds, ds.GetLayer("lucas_region_grow") if ds is not None else None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir1',
                        type=str, required=True,
                        help='Directory with RG areas (left).')
    parser.add_argument('--dir2',
                        type=str, required=True,
                        help='Directory with RG areas (right).')

    args = parser.parse_args()

    main(args.dir1, args.dir2)
