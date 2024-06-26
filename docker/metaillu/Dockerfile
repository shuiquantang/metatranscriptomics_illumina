#FROM brunneis/python:3.7.0-ubuntu-18.04
FROM brunneis/python:3.8.3-ubuntu-20.04
# stang/metaillu:v2

RUN apt-get update

RUN apt-get install -y wget zip

#install others
RUN apt-get update
RUN apt-get install -y cpanminus build-essential
RUN cpanm Parallel::Loops
RUN apt-get install -y python3-pip


#isntall sourmash
RUN pip3 install importlib-metadata
RUN apt-get install -y zlib1g-dev libblas-dev liblapack-dev
RUN apt-get install -y libjpeg8-dev libtiff-dev libfreetype6 libfreetype6-dev libwebp-dev libopenjp2-7-dev libopenjp2-7-dev
Run apt-get install -y libatlas-base-dev gfortran
RUN pip install numpy
RUN pip install pillow
RUN pip install scipy
RUN apt-get install -y libffi-dev libssl-dev
RUN pip install cffi
RUN pip install importlib-metadata
RUN pip install sourmash

RUN apt-get install -y default-jre
RUN cpanm Text::CSV_PP
RUN apt-get install -y libbz2-dev


ARG BUILD_DATE
ENV KRAKEN_VERSION=2.1.3 \
    TOOLS_DIR=/usr/local/NGS

LABEL name="Kraken2 ${KRAKEN_VERSION}" \
	  maintainer="Thomas J. Taylor <thomas@infotechsoft.com>" \
	  build-date="${BUILD_DATE}"

# Dependencies
RUN apt-get update && \
    apt-get upgrade -y perl && \
    apt-get install -y git make sed wget g++ rsync locales && \
    locale-gen --purge en_US.UTF-8 && \
    DEBIAN_FRONTEND="noninteractive" dpkg-reconfigure locales && \
    update-locale LANG=en_US.UTF-8 LANGUAGE=en_US.UTF-8 LC_ALL=en_US.UTF-8
    
ENV KRAKEN_SRC=https://github.com/DerrickWood/kraken2/archive/refs/tags/v${KRAKEN_VERSION}.tar.gz 

RUN apt-get install -y libz-dev

# Kraken2
RUN mkdir -p $TOOLS_DIR && \
    cd $TOOLS_DIR && \
    wget -qO- $KRAKEN_SRC | tar xz && \
    cd kraken2-*/ && \
    mkdir bin && \
    ./install_kraken2.sh ./bin/

ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}" \
    PATH="${PATH}:$TOOLS_DIR/kraken2-${KRAKEN_VERSION}/bin"

# Clean
RUN apt-get remove -y git make locales && apt-get autoclean -y


ENTRYPOINT [""]

