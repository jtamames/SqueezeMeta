FROM ubuntu:20.04

LABEL container.base.image="ubuntu:20.04"
LABEL software.name="DAS_Tool"
LABEL software.description="DAS Tool is an automated method that integrates the results of a flexible number of binning algorithms to calculate an optimized, non-redundant set of bins from a single assembly."
LABEL software.website="https://github.com/cmks/DAS_Tool"
LABEL software.license="BSD"
LABEL software.citation="Sieber et al., 2018, Nature Microbiology (https://doi.org/10.1038/s41564-018-0171-1)"

RUN DEBIAN_FRONTEND="noninteractive" apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -yq \
    autoconf \
    automake \
    cmake \
    git \
    libpcre3 \
    libpcre3-dev \
    libgsl0-dev \
    libgomp1 \
    lzma \
    ncbi-blast+ \
    r-base \
    ruby-full \
    wget

# Prodigal
RUN cd /tmp && \
    wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux && \
    mv prodigal.linux /bin/prodigal && \
    chmod +x /bin/prodigal

# Pullseq
RUN cd /tmp && \
    git clone https://github.com/bcthomas/pullseq.git && \
    cd pullseq && \
    ./bootstrap && \
    ./configure --prefix=/ && \
    make && \
    make install && \
    rm -rf /tmp/pullseq.zip /tmp/pullseq-master
# RUN cd /tmp && \
#     wget https://github.com/bcthomas/pullseq/releases/download/1.0.2/pullseq_v1.0.2_linux64.zip && \
#     unzip pullseq_v1.0.2_linux64.zip && \
#     mv pullseq seqdiff /bin && \
#     chmod +x /bin/pullseq && \
#     chmod +x /bin/seqdiff && \
#     rm -rf /tmp/pullseq*

# Diamond
RUN cd /tmp && \
  wget https://github.com/bbuchfink/diamond/releases/download/v2.0.14/diamond-linux64.tar.gz && \
  tar xfvz diamond-linux64.tar.gz && \
  mv diamond /bin/diamond && \
  chmod +x /bin/diamond && \
  rm -rf /tmp/diamond*

# Usearch
RUN cd /tmp && \
  wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz && \
  gunzip usearch11.0.667_i86linux32.gz && \
  mv usearch11.0.667_i86linux32 /bin/usearch && \
  chmod +x /bin/usearch

# DAS Tool
RUN mkdir -p /opt/DAS_Tool
ADD ./DAS_Tool /opt/DAS_Tool/DAS_Tool
ADD ./src /opt/DAS_Tool/src
ADD ./db.zip /opt/DAS_Tool/db.zip
# ADD ./sample_data /opt/DAS_Tool/sample_data
# ADD ./sample_output /opt/DAS_Tool/sample_output
RUN ln -s /opt/DAS_Tool/DAS_Tool /bin/DAS_Tool
RUN unzip -o /opt/DAS_Tool/db.zip -d /opt/DAS_Tool/db && \
  R -e "install.packages(c('data.table','magrittr','docopt'), repos='http://cran.us.r-project.org')" && \
  rm /opt/DAS_Tool/db.zip && \
  chmod +x /opt/DAS_Tool/DAS_Tool
