#!/bin/bash
# ./run_qiime.sh /home/ubuntu/trial in517 WYMNGZUNFWJ56WUCJGWETYEHS3DQHNDD in517.rawdata.zip in517.master.table.csv trial

while getopts w:p:s:r: option; do
  case "${option}" in
    w) work_dir=${OPTARG:1};;
    p) project_ID=${OPTARG:1};;
    s) random_str=${OPTARG:1};;
    r) run_ID=${OPTARG:1};;
  esac
done

#create folders

cd $work_dir

mkdir -p $work_dir/$project_ID/rawdata

#download rawdata of controls
if [ -n "$run_ID" ]
  then
    #MiSeq run date e.g. 170728
    s3_path2="s3://zymo-files/Shotgun_QC/$run_ID"
    aws s3 sync $s3_path2 $project_ID/rawdata
fi

s3_path="s3://zymo-microbiomics-service/epiquest/epiquest_$project_ID/$random_str"

aws s3 sync $s3_path/metadata/ $project_ID

aws s3 sync $s3_path/rawdata/$run_ID $project_ID/rawdata
