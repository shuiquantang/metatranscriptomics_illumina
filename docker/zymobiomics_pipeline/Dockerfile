# stang/zymobiomics_pipeline:v16


FROM stang/qiime1.9:v1

# install Java for trimmomatic adapter trimming
RUN apt-get --yes install default-jre

# install parallel perl
RUN perl -MCPAN -e "install Parallel::Loops"

# install krona plot
RUN apt-get install wget
RUN wget "https://github.com/marbl/Krona/releases/download/v2.6.1/KronaTools-2.6.1.tar"
RUN tar -xvf KronaTools-2.6.1.tar
RUN rm KronaTools-2.6.1.tar
RUN cd /KronaTools-2.6.1/ && perl install.pl

## required for Deep's heatmap script
RUN pip install seaborn==0.7.1

###############dependencis for lefse and TaxaHeatmap###########
## Install R

RUN Rscript -e 'install.packages("survival", repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("mvtnorm", repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("modeltools", repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("coin", repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("MASS", repo="http://cran.us.r-project.org")'


# install R packages
RUN Rscript -e 'install.packages("getopt", repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("gplots", repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("vegan", repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("RColorBrewer", repo="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("Heatplus", repo="http://cran.us.r-project.org")'

#install perl statistics module v15
RUN perl -MCPAN -e "install Statistics::Regression"

RUN apt-get -y install cpanminus
RUN perl -MCPAN -e "install Text::CSV_PP"

