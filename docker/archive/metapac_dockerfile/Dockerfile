FROM ubuntu:18.04

# stang/metapac:v1

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2

# Download and install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# Set up Miniconda environment variables
ENV PATH /opt/conda/bin:$PATH

# Configure Miniconda channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# install Minimap2
RUN apt-get install --yes wget tar bzip2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 
RUN tar -jxvf minimap2-2.24_x64-linux.tar.bz2
RUN mv minimap2-2.24_x64-linux /minimap2 && rm minimap2-2.24_x64-linux.tar.bz2
ENV PATH "$PATH:/minimap2/"

# Install hifiasm
RUN conda install -c bioconda hifiasm_meta

# Install fastqc
RUN conda install -c bioconda fastqc

# Install GTDB-Tk
RUN conda create -n gtdbtk -c conda-forge -c bioconda gtdbtk=2.1.1 numpy=1.23.1

# parallel::loops
RUN conda install -c bioconda perl-parallel-loops

# Install MetaBAT2
RUN conda install metabat2 -y

# Install SemiBin2
#RUN conda create -n SemiBin
#RUN conda activate SemiBin
RUN conda install -c conda-forge -c bioconda semibin


# Install DAS_Tool
RUN conda install -c bioconda das_tool -y

# download and unzip pb-metagenome-tools
RUN apt-get install -y zip unzip
RUN wget https://github.com/PacificBiosciences/pb-metagenomics-tools/archive/refs/tags/v2.0.0.zip
RUN unzip v2.0.0.zip
RUN mv pb-metagenomics-tools-2.0.0 pb-metagenomics-tools && rm v2.0.0.zip
ENV PATH "$PATH:/pb-metagenomics-tools/HiFi-MAG-Pipeline/scripts/"

# install seaborn for using plot scripts
RUN apt install -y python3-pip
RUN pip3 install seaborn

# install checkM2
RUN wget https://github.com/chklovski/CheckM2/archive/refs/tags/1.0.1.zip

RUN unzip 1.0.1.zip && rm 1.0.1.zip
RUN conda env create -n checkm2 -f CheckM2-1.0.1/checkm2.yml
ENV PATH /opt/conda/envs/checkm2/bin:$PATH
SHELL ["/bin/bash", "-c"]
RUN echo "source activate checkm2" >> ~/.bashrc
RUN pip install CheckM2
RUN conda install -n checkm2 seaborn

RUN pip3 install sourmash

RUN apt-get -y install cpanminus
RUN perl -MCPAN -e "install Text::CSV_PP"


ENV DEBIAN_FRONTEND=noninteractive LANG=en_US.UTF-8 LC_ALL=C.UTF-8 LANGUAGE=en_US.UTF-8

## GTDB-Tk
RUN python3 -m pip install gtdbtk
ENV GTDBTK_DATA_PATH=/home/ubuntu/metapac_database/gtdb_r207_v2
RUN apt-get update
RUN apt-get install -y hmmer mash prodigal

RUN wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip -q && \
    unzip pplacer-linux-v1.1.alpha19.zip && \
    mv pplacer-Linux-v1.1.alpha19/* /usr/bin && \
    rm pplacer-linux-v1.1.alpha19.zip && \
    rm -rf pplacer-Linux-v1.1.alpha19

RUN wget https://github.com/ParBLiSS/FastANI/releases/download/v1.32/fastANI-Linux64-v1.32.zip -q && \
    unzip fastANI-Linux64-v1.32.zip -d /usr/bin && \
    rm fastANI-Linux64-v1.32.zip

RUN wget http://www.microbesonline.org/fasttree/FastTree.c && \
    gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm && mv FastTree /usr/bin && \
    gcc -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeMP FastTree.c -lm && mv FastTreeMP /usr/bin && \
    rm FastTree.c
