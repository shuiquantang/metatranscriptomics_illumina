#FROM brunneis/python:3.7.0-ubuntu-18.04
FROM brunneis/python:3.8.3-ubuntu-20.04
# stang/metaillu:v1

RUN apt-get update

RUN apt-get install -y wget

# install diamond
RUN mkdir /diamond_db
RUN wget -O /diamond_db/NCBIATB.dmnd https://zymo-files.s3.amazonaws.com/WGS_Pipeline/shotgun_database/20200809_ref/diamond/NCBIATB.dmnd
RUN wget -O /diamond_db/VFDBS.dmnd https://zymo-files.s3.amazonaws.com/WGS_Pipeline/shotgun_database/20200809_ref/diamond/VFDBS.dmnd
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.15/diamond-linux64.tar.gz \
&& tar xzf diamond-linux64.tar.gz \
&& mv diamond /usr/bin
RUN chmod +x /usr/bin/diamond

#install others
RUN apt-get update
RUN apt-get install -y cpanminus build-essential
RUN cpanm Parallel::Loops


##install humann3
RUN apt-get install -y python3-pip
RUN pip3 install humann --no-binary :all:


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
RUN humann_config --update database_folders nucleotide /home/ubuntu/metaillu_database/humann3/chocophlan
RUN humann_config --update database_folders protein /home/ubuntu/metaillu_database/humann3/uniref90
RUN humann_config --update database_folders utility_mapping /home/ubuntu/metaillu_database/humann3/util_mapping

RUN apt-get install -y default-jre
RUN pip install khmer

RUN cpanm Text::CSV_PP

RUN apt-get install -y libbz2-dev

RUN pip install metaphlan=3.0


ENTRYPOINT [""]

