FROM ghcr.io/osgeo/gdal:alpine-normal-3.9.2

RUN apk add py-pip gcc python3-dev linux-headers musl-dev geos-dev py3-pytest
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m pip install -r /tmp/requirements.txt --break-system-packages
