# Docker for CONCOCT (http://github.com/BinPro/CONCOCT) v1.1.0
# VERSION 1.1.0
#
# This docker creates and sets up an Ubuntu environment with all
# dependencies for CONCOCT v1.1.0 installed.
#
# To login to the docker with a shared directory from the host do:
#
# docker run -v /my/host/shared/directory:/my/docker/location -i -t alneberg/concoct_1.1.0 /bin/bash
#

FROM ubuntu:18.04
COPY . /opt/CONCOCT

# Get basic ubuntu packages needed
RUN apt-get update -qq
RUN apt-get install -qq wget build-essential libgsl0-dev git zip unzip bedtools python-pip samtools

RUN pip install --upgrade pip

RUN wget --no-check-certificate https://github.com/BinPro/integration_test_data/archive/v1.1.tar.gz
RUN mkdir /opt/CONCOCT/tests/test_data/integration_test_data
RUN tar -xvzf v1.1.tar.gz -C /opt/CONCOCT/tests/test_data/integration_test_data --strip-components=1

# Install python dependencies and fetch and install CONCOCT 1.1.0
RUN cd /opt/CONCOCT;\
    pip install -r requirements.txt;\

RUN cd /opt/CONCOCT/;\
    python setup.py install

RUN cd /opt/CONCOCT/;\
    nosetests
