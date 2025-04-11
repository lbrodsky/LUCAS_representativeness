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

filebase=lucas_points
rm -f ${DATA_DIR}/${filename}_*.gpkg

declare -A country_codes
country_codes=(
  ["austria"]="AT"
  ["belgium"]="BE"
  ["bulgaria"]="BG"
  ["croatia"]="HR"
  ["cyprus"]="CY"
  ["czech-republic"]="CZ"
  ["denmark"]="DK"
  ["estonia"]="EE"
  ["finland"]="FI"
  ["france"]="FR"
  ["germany"]="DE"
  ["great-britain"]="UK"
  ["greece"]="EL"
  ["hungary"]="HU"
  ["ireland"]="IE"
  ["italy"]="IT"
  ["latvia"]="LV"
  ["lithuania"]="LT"
  ["luxembourg"]="LU"
  ["malta"]="MT"
  ["netherlands"]="NL"
  ["poland"]="PL"
  ["portugal"]="PT"
  ["romania"]="RO"
  ["slovakia"]="SK"
  ["slovenia"]="SI"
  ["spain"]="ES"
  ["sweden"]="SE"
)

for country_name in "${!country_codes[@]}"; do
    country=${country_codes[$country_name]}
    echo "Processing ${country}..."
    ogr2ogr -f GPKG ${DATA_DIR}/${country_name}_${filebase}_gps_${YEAR}.gpkg -overwrite PG:"dbname=lucas user=postgres" -nln eu_lucas_points_gps \
	    -sql "select point_id,geom,lc1_h,gps_prec,obs_dist,obs_type,obs_direct from data.lucas_points where survey_year = ${YEAR} AND nuts0 = '${country}'"
done

ogr2ogr -f GPKG ${DATA_DIR}/eu_${filebase}_thr_${YEAR}.gpkg PG:"dbname=lucas user=postgres" -nln eu_lucas_points_thr \
 	-sql "select point_id,geom_thr as geom,lc1_h,gps_prec,obs_dist,obs_type,obs_direct from data.lucas_points where survey_year = ${YEAR}"

exit 0
