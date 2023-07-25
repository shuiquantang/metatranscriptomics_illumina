#sudo su
#----------------------------------------
#set up ssh keys to git clone the pipeline scripts
aws s3 cp s3://zymo-files/Shuiquan/github_keys/id_rsa /root/.ssh/
chmod 400 /root/.ssh/id_rsa
ssh-keyscan -H github.com >> /root/.ssh/known_hosts
rm -r scripts/
git clone -b v1.1 git@github.com:shuiquantang/shotgun_metagenomics_with_Illumina.git
mv shotgun_metagenomics_with_Illumina scripts/
#----------------------------------------
# download reference database
rm -r /home/ubuntu/metaillu_database/ /home/ubuntu/database/
aws s3 cp s3://zymo-files/WGS_Pipeline/shotgun_database/metaillu_database/metaillu_database.tar /home/ubuntu/
tar xvf metaillu_database.tar
rm metaillu_database.tar
#----------------------------------------
# run_ID is for download QC data, and if you want to skip, put run_id="NA"
run_ID="nextseq1.230421"
projectID="in3743"
random_str="MRXGEYHIHSDNETYN"
sample_info_table="in3743.analysis.230421.csv"
seq_analysis_tag="5k"
report_tag="5k"
#----------------------------------------
#download data
/bin/bash /home/ubuntu/scripts/1.download.rawdata.sh -w=/home/ubuntu -p=$projectID -s=$random_str -r=$run_ID -m=$sample_info_table

# run the pipeline
# -R Trimmomatic controller, -M sourmash controller, -S StrainScan controller, -H humann3 controller, -D diamond controller, -V qiime and visualization controller, -k kmer size 51 by default
# 1, to run; 0, to skip and download the results from S3; 2, skip and do nothing
# if a customer wants new group analysis or you need to update a sample label, you only need to redo qiime and visualization,
# use: -R=0 -M=0 -S=2 -H=0 -D=0 -V=1

# if you want to only run sourmash for taxonomy profile and suppress humann and diamond,
# -R=1 -M=1 -S=2 -H=0 -D=0 -V=1

# do not try to skip read processing with -R=0.

# to set the ker size
# -k=21, non-human microbiome samples
# -k=51, human microbiome samples 

/bin/bash /home/ubuntu/scripts/2.run.analysis.sh -w=/home/ubuntu -p=$projectID -s=$random_str -i=$sample_info_table -x=$seq_analysis_tag -y=$report_tag -r=$run_ID -R=1 -M=1 -S=2 -H=1 -D=1 -V=1 -k=51 2>&1 > /home/ubuntu/$projectID/$projectID.run.log

