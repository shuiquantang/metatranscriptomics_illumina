FROM ubuntu:18.04
ENV PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh &&\
    mkdir /root/.conda &&\
    bash Miniconda3-latest-Linux-x86_64.sh -b &&\
    rm -f Miniconda3-latest-Linux-x86_64.sh

RUN conda install python=3.7 &&\
    conda config --add channels defaults &&\
    conda config --add channels bioconda &&\
    conda config --add channels conda-forge 

RUN conda config --add channels biobakery
RUN conda install humann -c biobakery

RUN conda install perl-app-cpanminus

RUN conda install -c bioconda perl-parallel-loops

RUN apt-get update
# install aws-cli
RUN apt-get --yes install curl unzip git
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install
RUN rm -r awscliv2.zip aws/

RUN conda install -c cyclus java-jre

RUN humann_config --update database_folders nucleotide /home/ubuntu/shotgun_db/chocophlan
RUN humann_config --update database_folders protein /home/ubuntu/shotgun_db/uniref90
RUN humann_config --update database_folders utility_mapping /home/ubuntu/shotgun_db/util_mapping

#RUN mkdir -p /home/ubuntu/shotgun_db/chocophlan /home/ubuntu/shotgun_db/uniref90  /home/ubuntu/shotgun_db/util_mapping /home/ubuntu/shotgun_db/mf4

#RUN humann_databases --download chocophlan full /home/ubuntu/shotgun_db/chocophlan --update-config yes

#RUN humann_databases --download uniref uniref90_diamond /home/ubuntu/shotgun_db/uniref90 --update-config yes

#RUN humann_databases --download utility_mapping full /home/ubuntu/shotgun_db/util_mapping --update-config yes

#RUN metaphlan --install --bowtie2db /home/ubuntu/shotgun_db/mf4 --nproc 16
