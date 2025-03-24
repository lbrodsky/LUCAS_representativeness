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
        r_tile = Path(dir_right, l_tile.name)
        if not Path(r_tile).exists():
            print("Right tile not found!")
            mismatch +=1
            continue

        if comp_tiles(l_tile, r_tile) is False:
            mismatch += 1

    print(f"Tiles processed: {count}")
    print(f"Mismatch: {mismatch}")

def comp_tiles(l_tile, r_tile):
    l_ds, l_layer = open_layer(str(l_tile))
    r_ds, r_layer = open_layer(str(r_tile))

    if l_layer is None and r_layer is None:
        return True

    if (l_layer is None and r_layer is not None) or (l_layer is not None and r_layer is None):
        print(f"Layer mismatch (left:{l_layer}, {r_ds.GetName()}:{r_layer}")
        return False
        
    if l_layer.GetFeatureCount() != r_layer.GetFeatureCount():
        print(f"Feature count mismatch ({l_layer.GetFeatureCount()}x{r_layer.GetFeatureCount()}: {r_ds.GetName()})!")
        return False
    
    for l_feat in l_layer:
        point_id = l_feat.GetField('point_id')
        r_layer.SetAttributeFilter(f"point_id = {point_id}")
        r_layer.ResetReading()
        r_feat = r_layer.GetNextFeature()

        for i in range(l_feat.GetFieldCount()):
            l_field_name = l_feat.GetFieldDefnRef(i).GetName()
            l_field_value = l_feat.GetField(i)
            r_field_value = r_feat.GetField(l_field_name)
            if str(l_field_value) != str(r_field_value): # avoid field type mismatch
                print(f"Field ({l_field_name}) mismatch: {l_field_value} vs {r_field_value} [{r_ds.GetName()}: {point_id}]!")
                return False

            if abs(l_feat.GetGeometryRef().GetArea() - r_feat.GetGeometryRef().GetArea()) > 1e-12:
                print(f"Point id {point_id} area mismatch!")
                return False

    l_ds = r_ds = None

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
