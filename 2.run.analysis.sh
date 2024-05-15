#!/bin/bash
# ./run_qiime.sh /home/ubuntu/trial in517 WYMNGZUNFWJ56WUCJGWETYEHS3DQHNDD in517.rawdata.zip in517.master.table.csv trial

while getopts w:p:s:i:r:x:y:R:M:S:H:D:V:k: option; do
    case "${option}" in
       w) work_dir=${OPTARG:1};;
       p) projectID=${OPTARG:1};;
       s) random_str=${OPTARG:1};;
       i) sample_info=${OPTARG:1};;
       r) rawdata_tag=${OPTARG:1};;
       x) seq_analysis_tag=${OPTARG:1};;
       y) report_tag=${OPTARG:1};;
       R) read_processing_controller=${OPTARG:1};;
       M) sourmash_controller=${OPTARG:1};;
       S) strainscan_controller=${OPTARG:1};;
       H) humann_controller=${OPTARG:1};;
       D) diamond_controller=${OPTARG:1};;
       V) visualization_controller=${OPTARG:1};;
       k) kmer_size=${OPTARG:1};;
    esac
done

kmer_size=${kmer_size:-51}

#run qiime
s3_path="s3://zymo-microbiomics-service/epiquest/epiquest_$projectID/$random_str"
home_dir="/home/ubuntu"
script="/home/ubuntu/scripts"
workdir="/home/ubuntu/$projectID"
ref_database="/home/ubuntu/metaillu_database"
start=`date +%s`
bin_path="$script/bin"

# docker_images
read_processing_docker_image='stang/trimreads:v2'
sourmash_docker_image='stang/metaillu:v2'
strainscan_docker_image='stang/strainscan:v1'
humann_docker_image='stang/metaphlan4:v2'
diamond_docker_image='stang/strainscan:v1'
visualization_docker_image='stang/zymobiomics_pipeline:v16'
report_docker_image='stang/shotgun-report:v1'

# step perl scripts
read_processing_script="A.read_processing.pl"
sourmash_script="B.sourmash.pl"
strainscan_script="C.strainscan.pl"
humann_script="D.humann3.pl"
diamond_script="E.diamond.pl"
visualization_script="F.visualization.pl"
report_script="$script/zymobiomics-shotgun-report/exec.py"

# result zip files
read_processing_zipfile="$projectID.$seq_analysis_tag.read_processing.tar.gz"
sourmash_zipfile="$projectID.$seq_analysis_tag.sourmash.tar.gz"
strainscan_zipfile="$projectID.$seq_analysis_tag.strainscan.tar.gz"
humann_zipfile="$projectID.$seq_analysis_tag.humann3.tar.gz"
diamond_zipfile="$projectID.$seq_analysis_tag.diamond.tar.gz"
visualization_zipfile="$projectID.$report_tag.visualization.tar.gz"
report_zipfile="$projectID.$report_tag.report.zip"

# run time recording file
time_file="runtime.$report_tag.txt"
rm $time_file
touch $time_file


#### step function ####
function run_step() {
  local step_name=$1
  local step_controller_var_name="${step_name}_controller"
  local step_controller=${!step_controller_var_name}
  local step_docker_image_var_name="${step_name}_docker_image"
  local step_docker_image=${!step_docker_image_var_name}
  local step_script_var_name="${step_name}_script"
  local step_script=${!step_script_var_name}
  local step_zipfile_var_name="${step_name}_zipfile"
  local step_zipfile=${!step_zipfile_var_name}
  #sudo docker pull $step_docker_image
  if [[ $step_controller -eq 1 ]]
  then
    docker run -v $home_dir:$home_dir -w $workdir -i $step_docker_image perl -I $bin_path $script/$step_script -o $step_name -t $sample_info -r $ref_database -k $kmer_size -R $random_str -T $rawdata_tag -P $projectID 2>&1 > $step_name.log.txt
    cd $workdir
    tar -I pigz -cf $step_zipfile $step_name/ #--exclude='*.fastq.gz'
    aws s3 cp $step_zipfile $s3_path/intermediate_results/$step_zipfile --acl public-read-write 2>&1 > /dev/null
    rm $step_zipfile
  elif [[ $step_controller -eq 0 ]]
  then
    if aws s3 ls $s3_path/intermediate_results/$step_zipfile >/dev/null 2>&1; then
        aws s3 cp --quiet $s3_path/intermediate_results/$step_zipfile $step_zipfile 2>&1 > /dev/null
        tar -I pigz -xf $step_zipfile
        rm $step_zipfile
        echo "Downloaded the intermediate file of step \"$step_name\" from s3".
    else
        echo "No intermediate file of step \"$step_name\" exists; $step_name was skipped"
    fi
  else
    echo "Skip the $step_name step completely"
  fi
  echo -e "\nafter $step_name:" >> $time_file
  expr `date +%s` - $start >> $time_file
}


### run through processing steps ###
cd $workdir

steps=("read_processing" "sourmash" "strainscan" "humann" "diamond" "visualization")

for step in "${steps[@]}"; do
 rm -r -f $step
done

for step in "${steps[@]}"; do
 run_step $step
done


#### generate the html report ####
cd $workdir
rm -r -f $projectID.$report_tag

docker run -v $home_dir:$home_dir -w $workdir -i $report_docker_image python3 $report_script -p $projectID.$report_tag -m $sample_info > report.generation.log.txt

docker run -v $home_dir:$home_dir -w $workdir -i $visualization_docker_image perl -I $bin_path $script/G.cleanup.pl -t $sample_info -o $projectID.$report_tag 2>&1 > report.generation.log.txt

cd $workdir
zip -r -q $report_zipfile $projectID.$report_tag/
aws s3 cp $report_zipfile $s3_path/report/$report_zipfile --acl public-read-write 2>&1 > /dev/null
rm $report_zipfile
echo "analysis done!"

echo -e "\nafter HTML report:" >> $time_file
expr `date +%s` - $start >> $time_file

#### shut down instance ####
end=`date +%s`
runtime=$(expr $end - $start)
aws s3 cp $time_file $s3_path/logs/run_time.$report_tag.txt 2>&1 > /dev/null
echo "$runtime seconds" >> report.generation.log.txt
if [ $runtime -gt 300 ] 
then
 instance_id="$(ec2metadata --instance-id)"
 aws ec2 stop-instances --region us-east-1 --instance-ids $instance_id
fi
