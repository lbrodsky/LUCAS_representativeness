#!/usr/bin/bash -e

# Run region grow process

script_dir="$(dirname "$(realpath "$0")")"

if test -z $1; then
    echo "Usage: $0 data_dir year version"
    exit 1
fi
if test -z $2; then
    echo "Usage: $0 data_dir year version"
fi
if test -z $3; then
    echo "Usage: $0 data_dir year version"
fi

data_dir="$1"
year="$2"
version="$3"
for dir in `ls ${data_dir}/osm_clcplus/${year}/*osm_clcplus*/ -d`; do
    country=`basename $dir | cut -d'_' -f1`
    echo "Processing ${country} v${version}..."

    docker run --rm  \
           --user $UID \
           -v ${script_dir}/../region_grow/:/opt/region_grow \
           -v ${data_dir}:/data \
           lucas_representativeness:latest \
           python3 /opt/region_grow/process_tiles.py \
           --tiles_dir /data/osm_clcplus/${year}/${country}_osm_clcplus_${year}/ \
           --lucas_points /data/lucas_points/${country}_lucas_points_gps_${year}.gpkg \
           --lucas_thr_points /data/lucas_points/eu_lucas_points_thr_${year}.gpkg \
           --dst_dir /data/lucas_representativeness/${year}/${country}_lucas_representativeness \
           --version ${version} \
           --log_level debug \
	   --workers 20
done

exit 0
