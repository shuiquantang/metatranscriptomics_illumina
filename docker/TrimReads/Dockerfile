FROM ubuntu:20.04
# stang/trimreads:v2
#To use docker with the access to GPU, you need to use '--gpus all'

#sudo docker run --gpus all stang/trimreads:v1 /bin/bash

ENV DEBIAN_FRONTEND=noninteractive

# install sortmerna
RUN apt-get update && apt install -y software-properties-common && apt install -y python3-pip
RUN pip3 install ribodetector

RUN apt-get install -y default-jre
RUN apt-get install -y cpanminus build-essential
RUN cpanm Parallel::Loops
RUN cpanm Text::CSV_PP
RUN apt-get install unzip
