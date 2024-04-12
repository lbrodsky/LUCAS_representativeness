#!/usr/bin/bash

# RERGION GROW software versions:
# RG: version 0.9.0 - updated similarity table, multi-class translation bug fixed
# RG: version 0.9.1 - sim tab level-3 
# RG: version 0.9.2 - paralelisation with joblib
# RG: version 0.9.3 - numpy y,x correction
# RG: version 0.9.4 - coords correction 
# RG: version 0.9.5 - asymetrical parch coords correction
# ---
# LOGING:
# at: version 3 - RG v.0.9.2
# be: version 1 - RG v.0.9.2   
# bg: version 1 - RG v.0.9.2 
# cy: version 1 - RG v.0.9.2 (Cyprus)  
# cz: version 9 - RG v.0.9.2  
# de: version 1 - RG v.0.9.2 
# dk: version 1 - RG v.0.9.2 
# ee: version 4 - RG v.0.9.2 (Estonia)
# es: version 1 - RG v.0.9.2 
# fi: version 1 - RG v.0.9.2 
# fr: version 1 - RG v.0.9.2  
# gr: version 2 - RG v.0.9.2 
# hr: version 1 - RG v.0.9.2 (Croatia) 
# hu: version 3 - RG v.0.9.2
# ie: version 5 - RG v.0.9.2
# it: version 4 - RG v.0.9.5
# lt: version 1 - RG v.0.9.2 (Lithuania)   
# lu: version 1 - RG v.0.9.2 (Luxemburg) 
# lv: version 1 - RG v.0.9.2 (Latvia) 
# nl: version 2 - RG v.0.9.3 old info - No OSM land cover over LUCAS point ID: 40003246!  
# mt: version 1 - RG v.0.9.2 
# pl: version 1 - RG v.0.9.2
# pt: version 1 - RG v.0.9.2 
# ro: version 1 - RG v.0.9.2
# se: version 1 - RG v.0.9.2 (Sweden)
# sk: version 4 - RG v.0.9.2 
# si: version 1 - RG v.0.9.2
# gb: version 1 - RG v.0.9.2 (UK) 

# Run all countries at once, updated cntr versions
declare -A cntrv

countries=('at' 'be' 'bg' 'cy' 'cz' 'de' 'dk' 'ee' 'es' 'fi' 'fr' 'gb' 'gr' 'hr' 'hu' 'ie' 'it' 'lt' 'lu' 'lv' 'nl' 'mt' 'pl' 'pt' 'ro' 'se' 'si' 'sk')

# countries=('si')   

for country in ${countries[@]}; do
    # version=${cntrv[$country]}
    version=28     
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
  -v /home/brodsky/src/ctu-geoharmonizer/land-cover/lucas/representativeness/region_grow/:/home/region_grow \
  -v /data/geoharmonizer/lucas/representativeness/osm2018/$osm/clcplus:/data/representativeness/osm2018/$osm/clcplus/ \
  -v /data/geoharmonizer/lucas/representativeness/LUCAS_points/:/data/representativeness/LUCAS_points \
  -v /data/geoharmonizer/lucas/representativeness/Region_Grow/$rg:/data/representativeness/Region_Grow/$rg \
  ctu-geoharmonizer:1.3 \
  python3 /home/region_grow/process_tiles.py \
  --tiles_dir /data/representativeness/osm2018/$osm/clcplus \
  --lucas_points /data/representativeness/LUCAS_points/${country}_lucas_points_2018.gpkg \
  --lucas_thr_points /data/representativeness/LUCAS_points/eu_lucas_points_thr.gpkg \
  --shp_thr 0.7 --region_max_size 100 \
  --workers 20 \
  -tdir /home/region_grow/tabels/ \
  --dst_dir /data/representativeness/Region_Grow/$rg \
  --log_level debug \
  --version ${version} 
  
  # &> ./ctu-geoharmonizer/land-cover/lucas/representativeness/region_grow/runs/nohup_eu_rect.out
  # nohup_${country}_v${version}.out & 
  
done # country / version loop 

