FROM ubuntu:bionic
WORKDIR /opt

# Ubuntu
# =============================================================
RUN apt-get update &&\
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata &&\
    rm -rf /var/lib/apt/lists/* &&\
    apt-get update &&\
    apt-get install -y \
                       ca-certificates \
                       gnupg2 \
                       apt-utils \
                       build-essential \
                       software-properties-common &&\

# Curl
# =============================================================
apt-get update &&\
apt-get install -y \
libssl-dev \
libxml2-dev \
libcurl4-openssl-dev \
libgmp-dev \
pandoc \
imagemagick \
curl &&\


# R-base
# =============================================================
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 &&\
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/"  &&\
apt-get update &&\
apt-get install -y r-base

# R-packages
# =============================================================
RUN R -e "install.packages('data.table',dependencies=TRUE)"
RUN R -e "install.packages('ggplot2',dependencies=TRUE)"
RUN R -e "install.packages('plotly',dependencies=TRUE)"
RUN R -e "install.packages('rjson',dependencies=TRUE)"
RUN R -e "install.packages('stringr',dependencies=TRUE)"
RUN R -e "install.packages('RColorBrewer',dependencies=TRUE)"
RUN R -e "install.packages('pander',dependencies=TRUE)"
RUN R -e "install.packages('rmarkdown',dependencies=TRUE)"
RUN R --no-save -e 'install.packages("BiocManager"); BiocManager::install();'
RUN R --no-save -e 'BiocManager::install(c("GenomicRanges"))'
#RUN R -e "source('https://bioconductor.org/biocLite.R'); biocLite(c(GenomicRanges));"
#RUN R -e "install.packages('vegan',dependencies=TRUE)"
#RUN R -e "install.packages('ape',dependencies=TRUE)"

# Install python deps
RUN apt-get install python3.6
RUN apt -y install python3-pip
RUN pip3 install pandas awscli


RUN R -e "install.packages('vegan',dependencies=TRUE)"
RUN R -e "install.packages('ape',dependencies=TRUE)"
