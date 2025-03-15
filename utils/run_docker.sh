#!/usr/bin/bash

# Run Region Grow process

# Run all countries at once, updated cntr versions
declare -A cntrv

countries=('at' 'be' 'bg' 'cy' 'cz' 'de' 'dk' 'ee' 'es' 'fi' 'fr' 'gb' 'gr' 'hr' 'hu' 'ie' 'it' 'lt' 'lu' 'lv' 'nl' 'mt' 'pl' 'pt' 'ro' 'se' 'si' 'sk')

# countries=('cz')

for country in ${countries[@]}; do
    # version=${cntrv[$country]}
    version=44
    echo $country $version


# DOCKER MANAGEMENT
# docker container ls
# docker container stop ID
# docker container remove ID

# PATH
osm=osm_${country} # osm_czech-republic
rg=region_grow_${country} 


# Docker ctu-geoharmonizer v1.3 with installed shapely 

docker run --rm  \
  --user $UID \
  -v /home/brodsky/src/LUCAS_representativeness/region_grow/:/home/region_grow \
  -v /data/geoharmonizer/lucas/representativeness/osm2018/$osm/merged_osm_clcplus:/data/representativeness/osm2018/$osm/merged_osm_clcplus/ \
  -v /data/geoharmonizer/lucas/representativeness/LUCAS_points/:/data/representativeness/LUCAS_points \
  -v /data/geoharmonizer/lucas/representativeness/Region_Grow/$rg:/data/representativeness/Region_Grow/$rg \
  lucas_representativeness:latest \
  python3 /home/region_grow/process_tiles.py \
  --tiles_dir /data/representativeness/osm2018/$osm/merged_osm_clcplus \
  --lucas_points /data/representativeness/LUCAS_points/${country}_lucas_points_2018.gpkg \
  --lucas_thr_points /data/representativeness/LUCAS_points/eu_lucas_points_thr.gpkg \
  --shp_thr 0.5 --region_max_size 100 \
  --workers 22 \
  -tdir /home/region_grow/tables/ \
  --dst_dir /data/representativeness/Region_Grow/$rg \
  --log_level debug \
  --version ${version} 
  
done # country / version loop 

