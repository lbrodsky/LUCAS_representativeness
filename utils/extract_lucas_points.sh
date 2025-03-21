#!/bin/bash -e

if test -z $1; then
    echo "year not specified"
    exit 1
fi

YEAR=$1
DATA_DIR=/data/lucas_points
DB=lucas

mkdir -p ${DATA_DIR}
createdb -U postgres $DB
7z e -so ${DATA_DIR}/db_st_lucas_dump.sql.7z | sed '/^GRANT /d; /^ALTER DEFAULT PRIVILEGES /d' | psql -U postgres -d $DB

filename=eu_lucas_points
rm -f ${DATA_DIR}/${filename}_*.gpkg
ogr2ogr -f GPKG ${DATA_DIR}/${filename}_gps_${YEAR}.gpkg -overwrite PG:"dbname=lucas user=postgres" -nln eu_lucas_points_gps \
	-sql "select point_id,geom_gps as geom,lc1_h,gps_prec,obs_dist,obs_type,obs_direct from data.lucas_points where geom_gps is not NULL and survey_year = ${YEAR}"
ogr2ogr -f GPKG ${DATA_DIR}/${filename}_thr_${YEAR}.gpkg PG:"dbname=lucas user=postgres" -nln eu_lucas_points_thr \
	-sql "select point_id,geom_thr as geom,lc1_h,gps_prec,obs_dist,obs_type,obs_direct from data.lucas_points where geom_thr is not NULL and survey_year = ${YEAR}"

exit 0
