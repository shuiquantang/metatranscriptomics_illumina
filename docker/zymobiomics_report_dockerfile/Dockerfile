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
curl


# R-base
# =============================================================
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt-get update
RUN apt-get install -y r-base

# R-packages
# =============================================================
RUN Rscript -e 'install.packages("data.table",dependencies=TRUE, repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("ggplot2",dependencies=TRUE, repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("plotly",dependencies=TRUE, repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("rjson",dependencies=TRUE, repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("stringr",dependencies=TRUE, repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("RColorBrewer",dependencies=TRUE, repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("pander",dependencies=TRUE, repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("rmarkdown",dependencies=TRUE, repo="http://cran.us.r-project.org")'
RUN Rscript --no-save -e 'install.packages("BiocManager",, repo="http://cran.us.r-project.org"); BiocManager::install();'
RUN Rscript --no-save -e 'BiocManager::install(c("GenomicRanges"))'
#RUN R -e "source("https://bioconductor.org/biocLite.R"); biocLite(c(GenomicRanges));"


# Install python deps
RUN apt-get install python3.6
RUN apt -y install python3-pip
RUN pip3 install pandas awscli
RUN Rscript -e 'install.packages("vegan",dependencies=TRUE, repo="http://cran.us.r-project.org")'
