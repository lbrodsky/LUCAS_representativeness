# LUCAS representativeness

Code to compute spatially representative areas. Original
shape-constrained region-growing algorithm was specifically designed
to enhance the accuracy and consistency of machine learning approaches
for land cover and land use mapping by ensuring that the delineated
areas accurately represent the surrounding landscape. The algorithm
takes as input the OSM/CLC+ land cover product (ZENODO TBD), which
provides a harmonized and detailed land cover dataset suitable for
large-scale analysis.

The tool is demonstrated on LUCAS points in 2018 retrieved from
[ST_LUCAS](https://geoforall.fsv.cvut.cz/st_lucas/) system.

## Build Docker image

```
docker build --tag lucas_representativeness:latest docker/
```

## Run on sample data

Run Docker container on sample data:

```
docker run --rm  --user `id -u` \
 -v `pwd`:/opt \
 -v ./tests/data:/data \
 lucas_representativeness:latest \
 python3 /opt/region_grow/process_tiles.py \
 --tiles_dir /data/tiles \
 --lucas_points /data/sample_points.gpkg \
 --lucas_thr_points /data/sample_thr_points.gpkg \
 --dst_dir /data/output \
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

## Run complete EU coverage computation

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
docker exec -ti lucas-postgis bash -c "apt-get update; apt-get -y install p7zip-full gdal-bin"
docker exec -ti lucas-postgis bash /opt/utils/extract_lucas_points.sh 2018
docker stop lucas-postgis
docker rm lucas-postgis
```

Download OSM/CLC+ product (TBD):

```
mkdir -p data/osm_clcplus/2018/
wget ...
(cd data/osm_clcplus/2018/; for f in *.7z; do 7z x $f; done)
```

### Perform RG area computation

```
./utils/run_docker.sh ./data
```

Output data is stored in `lucas_representativeness` directory.

### Postprocessing

```
docker run --rm --user `id -u` \
 -v `pwd`:/opt -v ./data:/data lucas_representativeness:latest \
 python3 /opt/utils/gridding_repre_polygons.py \
 --data_path ./data/lucas_representativeness/2018/
```

```
docker run --rm --user 1001 -v /home/landamar/git/st_lucas/LUCAS_representativeness:/opt -v ./data:/data lucas_representativeness:latest  python3 /opt/utils/merge_countries.py --src_path ./data/lucas_representativeness/2018/ --dst_path ./data/lucas_representativeness/2018/merged
```