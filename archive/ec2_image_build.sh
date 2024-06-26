# Name: ZymoBIOMICS_shotgun_Pipeline_v1
# EC2 ubuntu instance base image: Ubuntu Server 18.04 LTS (HVM), SSD Volume Type - ami-0ac80df6eff0e70b5 (64-bit x86)
# 192Gb hardware storage is required

sudo apt-get update
sudo apt-get --yes upgrade

sudo apt-get --yes install zip
sudo apt-get --yes install python

# install AWS-cli
curl "https://s3.amazonaws.com/aws-cli/awscli-bundle.zip" -o "awscli-bundle.zip"
unzip awscli-bundle.zip
sudo ./awscli-bundle/install -i /usr/local/aws -b /usr/local/bin/aws
rm -r awscli-bundle/ awscli-bundle.zip

#install docker
sudo apt-get --yes install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common

curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

sudo apt-key fingerprint 0EBFCD88

sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
sudo apt-get update
sudo apt-get --yes install docker-ce docker-ce-cli containerd.io


#sudo apt-get install pigz
#sudo apt-get install tar
#looks like ubuntu 18.04 has pre-installed these two packages


#download shotgun database 1
database="20200327"
mkdir database/
aws s3 sync s3://zymo-files/WGS_Pipeline/shotgun_database/$database database/

#download shotgun database 2
mkdir /home/ubuntu/database/zymo_centrifuge/genomes
aws s3 sync s3://zymo-files/WGS_Pipeline/shotgun_database/20200327_genomes/ /home/ubuntu/database/zymo_centrifuge/genomes/

#download docker images
qiime_docker_image="stang/qiime1.9:v1"
centrifuge_docker_image="stang/shotgun_metagenomics:v1"
humann2_docker_image="stang/humann2:v2"
report_docker_image="aishaniprem/zymobiomics-shotgun-image:v1"
sudo docker pull $qiime_docker_image
sudo docker pull $centrifuge_docker_image
sudo docker pull $humann2_docker_image
sudo docker pull $report_docker_image

