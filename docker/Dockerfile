FROM ubuntu:22.04

RUN apt-get update && apt-get install --yes software-properties-common gettext-base
RUN add-apt-repository --yes ppa:ubuntugis/ubuntugis-unstable
RUN apt-get update

RUN apt-get install --yes grass-core grass-doc grass-dev
RUN apt-get install --yes python3-pip

COPY requirements.txt /tmp/requirements.txt
RUN python3 -m pip install -r /tmp/requirements.txt

RUN sed -i '/^#\sdeb-src /s/^# *//' /etc/apt/sources.list
RUN apt-get update
RUN apt-get -y build-dep grass-core
RUN mkdir /opt/src
WORKDIR /opt/src
RUN git clone https://github.com/osgeo/grass.git
WORKDIR /opt/src/grass
RUN git checkout releasebranch_8_3
COPY configure.sh /tmp/configure.sh
RUN /tmp/configure.sh
RUN make -j8
RUN make install
