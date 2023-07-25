FROM ubuntu:18.04
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    mkdir /root/.conda &&\
    bash Miniconda3-latest-Linux-x86_64.sh -b &&\
    rm -f Miniconda3-latest-Linux-x86_64.sh

RUN conda install python=3.7 &&\
    conda config --add channels defaults &&\
    conda config --add channels bioconda


RUN conda install sibeliaz=1.2.2 numpy=1.17.3 pandas=1.0.1 biopython=1.74 scipy=1.3.1 scikit-learn=0.23.1 treelib=1.6.1
RUN apt-get update
RUN apt-get install -y python3-pip
RUN pip3 install bidict==0.21.3

RUN apt-get update
RUN apt-get install -y git

RUN git clone https://github.com/liaoherui/StrainScan.git /strainscan

RUN chmod 755 /strainscan/library/jellyfish-linux
RUN chmod 755 /strainscan/library/dashing_s128

RUN echo "conda activate strainscan" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

RUN apt-get update
RUN apt-get install -y cpanminus build-essential
RUN cpanm Parallel::Loops

RUN apt-get update
RUN DEBIAN_FRONTEND="noninteractive" apt-get install -y awscli

RUN mkdir /diamond_db

RUN wget -O /diamond_db/NCBIATB.dmnd https://zymo-files.s3.amazonaws.com/WGS_Pipeline/shotgun_database/20200809_ref/diamond/NCBIATB.dmnd

RUN wget -O /diamond_db/VFDBS.dmnd https://zymo-files.s3.amazonaws.com/WGS_Pipeline/shotgun_database/20200809_ref/diamond/VFDBS.dmnd

RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.15/diamond-linux64.tar.gz \
&& tar xzf diamond-linux64.tar.gz \
&& mv diamond /usr/bin
RUN chmod +x /usr/bin/diamond

