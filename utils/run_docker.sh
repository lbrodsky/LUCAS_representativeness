#!/usr/bin/bash -e

# Run region grow process

script_dir="$(dirname "$(realpath "$0")")"

# settings
version=1
data_dir=/data/geoharmonizer/lucas/representativeness/2018
workers=10
shp_thr=0.5
region_max_size=100

# run all countries at once, updated cntr versions
declare -A cntrv

countries=('at' 'be' 'bg' 'cy' 'cz' 'de' 'dk' 'ee' 'es' 'fi' 'fr' 'gb' 'gr' 'hr' 'hu' 'ie' 'it' 'lt' 'lu' 'lv' 'nl' 'mt' 'pl' 'pt' 'ro' 'se' 'si' 'sk')
for country in ${countries[@]}; do
    echo "Processing ${country} v${version}..."

    # define paths
    osm=osm_${country}
    rg=region_grow_${country} 

    # run docker container
    docker run --rm  \
           --user $UID \
           -v ${script_dir}/../region_grow/:/opt/region_grow \
           -v ${data_dir}:/data \
           lucas_representativeness:latest \
           python3 /opt/region_grow/process_tiles.py \
           --tables_dir /opt/region_grow/tables/ \
           --tiles_dir /data/${osm}/ \
           --lucas_points /data/lucas_points/${country}_lucas_points.gpkg \
           --lucas_thr_points /data/lucas_points/eu_lucas_points_thr.gpkg \
           --shp_thr $shp_thr --region_max_size $region_max_size \
           --workers $workers \
           --dst_dir /data/region_grow/$rg \
           --log_level info \
           --version ${version}
done

exit 0
