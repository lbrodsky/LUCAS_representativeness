# LUCAS representativeness

Code to compute LUCAS representativeness regions.

## Docker

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