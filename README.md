# LUCAS representativeness

Code to compute LUCAS representativeness regions.

## Build Docker image and run on sample data

Build docker image:

```
docker build --tag lucas_representativeness:latest docker/
```

Run Docker container on sample data:

```
docker run --rm  --user `id -u` \
 -v `pwd`:/opt \
 -v ./tests/data:/data \
 lucas_representativeness:latest \
 python3 /opt/region_grow/process_tiles.py \
 --tiles_dir /data/tiles \
 --tables_dir /opt/region_grow/tables \
 --lucas_points /data/sample_points.gpkg \
 --lucas_thr_points /data/sample_thr_points.gpkg \
 --dst_dir /data/output \
 --shp_thr 0.5 \
 --region_max_size 100 \
 --version 1
```

Run tests on sample data:

```
docker run --rm  --user `id -u` \
 -v `pwd`:/opt \
 -v ./tests/data:/data \
 lucas_representativeness:latest \
 python3 -m pytest /opt/tests/test_sample.py
```

## Full EU coverage computation

### Data preparation

Create data dir:

```
mkdir data
```

Extract LUCAS points for 2018:

```
wget https://geoforall.fsv.cvut.cz/extracts/st_lucas/db_st_lucas_dump.sql.7z -P data
docker pull postgis/postgis:16-3.4
docker run --name lucas-postgis -v ./data:/data -v `pwd`/utils:/opt/utils -e POSTGRES_PASSWORD=lucas -d postgis/postgis:16-3.4
docker exec -ti lucas-postgis bash -c "apt-get update; apt-get -y install p7zip gdal-bin"
docker exec -ti lucas-postgis bash /opt/utils/extract_lucas_points.sh 2018
docker stop lucas-postgis
docker rm lucas-postgis
```
