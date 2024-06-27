#sudo su
work_dir=/storage/metatrans
mkdir -p $work_dir
cd $work_dir
#----------------------------------------
#set up ssh keys to git clone the pipeline scripts
aws s3 cp s3://zymo-files/Shuiquan/github_keys/id_rsa /root/.ssh/
chmod 400 /root/.ssh/id_rsa
ssh-keyscan -H github.com >> /root/.ssh/known_hosts
rm -r scripts/
git clone -b v1.3 git@github.com:shuiquantang/metatranscriptomics_illumina.git
mv metatranscriptomics_illumina scripts/
#----------------------------------------
# download reference database
rm -r $work_dir/database/
mkdir $work_dir/metaillu_database/
aws s3 sync s3://zymo-files/WGS_Pipeline/shotgun_database/metaillu_database_20240524/ $work_dir/metaillu_database/
#----------------------------------------
# run_ID is for download QC data, and if you want to skip, put run_id="NA"
work_dir="/storage/metatrans"
run_ID="20230322" 
projectID="in1756" 
random_str="QLACGTDAED4NNFSCKWB3NEY3MP2BLYYE" 
sample_info_table="in1756.pe.analysis.csv" 
seq_analysis_tag="RNA_debug" 
report_tag="RNA_debug"
#----------------------------------------
#download data
/bin/bash $work_dir/scripts/1.download.rawdata.sh -w=$work_dir -p=$projectID -s=$random_str -r=$run_ID

# run the pipeline
# -R Trimmomatic controller, -M sourmash controller, -S StrainScan controller, -H humann3 controller, -D diamond controller, -V qiime and visualization controller, -k kmer size 51 by default
# 1, to run; 0, to skip and download the results from S3; 2, skip and do nothing
# if a customer wants new group analysis or you need to update a sample label, you only need to redo qiime and visualization,
# use: -R=0 -M=0 -S=2 -H=0 -D=0 -V=1

# if you want to only run sourmash for taxonomy profile and suppress humann and diamond,
# -R=1 -M=1 -S=2 -H=0 -D=0 -V=1

# do not try to skip read processing with -R=0.

 /bin/bash -x $work_dir/scripts/2.run.analysis.sh -w=$work_dir -p=$projectID -s=$random_str -i=$sample_info_table -x=$seq_analysis_tag -y=$report_tag -r=$run_ID -R=1 -M=1 -P=1 -S=2 -H=2 -D=1 -V=1 2>&1 > $work_dir/$projectID/$projectID.run.log

