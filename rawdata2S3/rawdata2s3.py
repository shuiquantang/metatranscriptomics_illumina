#!usr/bin/python3.6
__author__= 'Aishani Prem'
__email__='aprem@zymoresearch.com'
import argparse
import pandas as pd
import os
import re

#---Command Line Arguments----
parser=argparse.ArgumentParser(
      description="This script uses python3. Please make sure you have the correct version installed. This script is used to upload the fasta files to s3")

parser.add_argument('-m','--metadata', help='Input the metadata folder',required=True)
parser.add_argument('-r', '--rawdata', help='Input the rawdata folder' ,required=True)
parser.add_argument('-d', '--date', help='Enter the run date of sequencing' ,required=True)
parser.add_argument('-p', '--project', help='Enter the project info file, which contains the list of projects that were included in the run' ,required=True)
parser.add_argument('-s', '--sequencer', help='Enter the name of the sequencer used Ex MiniSeq1 MiSeq' ,required=True)


args = parser.parse_args()
metadata = args.metadata
rawdata = args.rawdata
date = args.date
pro = args.project
seq = args.sequencer

pro = pd.read_csv(pro, sep=',', header=0)
#print(pro)

fasta = os.listdir(rawdata)
#print(fasta)

os.system('mkdir new_data')


##Rename files
for name in fasta:
	if name.endswith('.fastq.gz'):
		newname = name.split('_',1)[0]
		newname= newname.replace('-','_')
		read = '_R' + name.split('_R')[1].replace('_001','')
		newname = newname.lower() + read
		os.system ('cp rawdata/%s new_data/%s' %(name, newname))
	else:
		print("%s is not formated properly" %name)


##Upload the rawdata to s3
new_fasta = os.listdir('new_data')
#os.system('mkdir new_data/QC')
for index, row in pro.iterrows():
	project = row["#projectIDs"]
	#print(row["#projectIDs"])
	os.system('mkdir new_data/%s' %project)
	os.system("mv new_data/%s*fastq.gz new_data/%s/" %(project, project))
	os.chdir("/home/ubuntu/new_data/%s/" %(project) )
	zip_file = "%s.rawdata.%s.zip" %(project,date)
	os.system("zip -r %s.rawdata.%s.zip *" %(project,date))
	if project.startswith('in') or project.startswith('zr'):
		s3_path= "s3://epiquest/epiquest_%s/%s/rawdata/%s.rawdata.%s.zip" %(project, row["random_strings"], project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print(s3_path)
		print("Rawdata from %s copied to %s" %(project, s3_path))
	elif project.startswith('m'):
		s3_path= "s3://midog/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
		#print(project)
	elif project.startswith('pbe'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Ear/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
	elif project.startswith('pbna'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Nail/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
	elif project.startswith('pbno'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Nose/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
	elif project.startswith('pbs'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Sputum/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
	elif project.startswith('pbt'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Throat/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
	elif project.startswith('pbu'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Urine/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
	elif project.startswith('pbv'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Vaginal/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
	elif project.startswith('pbw'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Wound/Projects/%s/rawdata/%s.rawdata.%s.zip"%(project, project, date)
		os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
		print("Rawdata from %s copied to %s" %(project, s3_path))
	else:
		pass
	os.chdir("/home/ubuntu/")



##Upload QC files:
os.system('mkdir new_data/QC')
os.system("mv new_data/*fastq.gz new_data/QC/")
zip_file = "QC.zip"
os.chdir("/home/ubuntu/new_data/QC/" )
os.system("zip -r %s *" %(zip_file))
if seq.startswith('MiniSeq'):
	ctvalues="/home/ubuntu/metadata/ctvalues_input.csv"
	s3_path= "s3://zymo-files/Shuiquan/MiSeq_QC/%s.%s/%s" %(seq,date,zip_file)
	os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
	s3_path= "s3://zymo-files/Shuiquan/MiSeq_QC/%s.%s/ctvalues_input.csv"%(seq,date)
	os.system("aws s3 cp %s %s --acl public-read-write" %(ctvalues, s3_path))	
else:
	s3_path= "s3://zymo-files/Shuiquan/MiSeq_QC/%s/%s" %(date,zip_file)
	os.system("aws s3 cp %s %s --acl public-read-write" %(zip_file, s3_path))
os.chdir("/home/ubuntu/")


#Save the analysis files to the right locations

#metadata = os.listdir (metadata)
#print(metadata)

for index, row in pro.iterrows():
	newname = "%s.master.table.%s.csv" %(row["#projectIDs"], date)
	#os.system("cp metadata/%s metadata/%s")
	project = row["#projectIDs"]
	if project.startswith('in') or project.startswith('zr'):
		s3_path =  "s3://epiquest/epiquest_%s/%s/metadata/%s" %(project, row["random_strings"], newname)
		os.system("aws s3 cp metadata/%s %s" %(row["metadata_table"], s3_path )) 
	elif project.startswith('m'):
		s3_path= "s3://midog/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
		s3_path= "s3://midog/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))
	elif project.startswith('pbe'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Ear/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
		s3_path= "s3://precisionbiome/PrecisionBIOME_Ear/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))
	elif project.startswith('pbu'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Urine/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
		s3_path= "s3://precisionbiome/PrecisionBIOME_Urine/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))
	elif project.startswith('pbw'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Wound/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
		s3_path= "s3://precisionbiome/PrecisionBIOME_Wound/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))		
	elif project.startswith('pbv'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Vaginal/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
		s3_path= "s3://precisionbiome/PrecisionBIOME_Vaginal/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))
	elif project.startswith('pbna'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Nail/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
		s3_path= "s3://precisionbiome/PrecisionBIOME_Nail/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))
	elif project.startswith('pbno'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Nose/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
		s3_path= "s3://precisionbiome/PrecisionBIOME_Nose/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))
	elif project.startswith('pbs'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Sputum/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
		s3_path= "s3://precisionbiome/PrecisionBIOME_Sputum/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))
	elif project.startswith('pbt'):
		s3_path= "s3://precisionbiome/PrecisionBIOME_Throat/Projects/%s/metadata/%s"%(row["#projectIDs"],newname)
		os.system("aws s3 cp metadata/%s %s"  %(row["metadata_table"], s3_path ))
		s3_path= "s3://precisionbiome/PrecisionBIOME_Throat/Projects/%s/metadata/%s"%(row["#projectIDs"],sample_info)
		os.system("aws s3 cp metadata/%s %s" %(row["sample_info_table"], s3_path ))

## Create manuals
#read the manual template

template = open('/home/ubuntu/scripts/s3.folders/template.manual', 'r')
template = template.readlines()
print(template)

os.system("mkdir manuals")

for index, row in pro.iterrows():
	proj_manual = open('manuals/%s.%s.manual' %(row['#projectIDs'], date), 'w')
	template_pro= template.copy()
	#print(template_pro)
	for i in range (0, len(template_pro)):
		if template_pro[i].startswith ('projectID='):
			line = template_pro[i].strip('\n')
			newline = "%s\"%s\" \n" %(line, row['#projectIDs'])
			newline = newline.replace('nan','')
			print(newline)
			print(template_pro[i])
			template_pro[i]= newline
		if template_pro[i].startswith ('runID='):	
			line = template_pro[i].strip('\n')
			if seq.startswith('MiniSeq'):
				newline = "%s\"%s.%s\" \n" %(line,seq, date)
				newline = newline.replace('nan','')
			else:
				newline = "%s\"%s\" \n" %(line, date)
				newline = newline.replace('nan','')
			print(newline)
			template_pro[i]= newline
		if template_pro[i].startswith ('random_str='):
			line = template_pro[i].strip('\n')
			newline = "%s\"%s\" \n" %(line, row['random_strings'])
			newline = newline.replace('nan','')
			print(newline)
			template_pro[i]= newline
		if template_pro[i].startswith ('rawdata='):
			line = template_pro[i].strip('\n')
			newline = "%s\"%s.rawdata.%s.zip\" \n" %(line, row['#projectIDs'], date)
			newline = newline.replace('nan','')
			print(newline)
			template_pro[i]= newline
		if template_pro[i].startswith ('sampleinfo='):
			sample_info= "%s.sample.info.%s.csv" %(row["#projectIDs"], date)
			line = template_pro[i].strip('\n')
			newline = "%s\"%s\" \n" %(line, sample_info)
			newline = newline.replace('nan','')
			print(newline)
			template_pro[i]= newline
		if template_pro[i].startswith ('metadata='):
			metadata_name= "%s.master.table.%s.csv" %(row["#projectIDs"], date)
			line = template_pro[i].strip('\n')
			newline = "%s\"%s\" \n" %(line, metadata_name)
			newline = newline.replace('nan','')
			print(newline)
			template_pro[i]= newline
		if template_pro[i].startswith ('tag='):
			line = template_pro[i].strip('\n')
			newline = "%s\"%s\" \n" %(line, date)
			newline = newline.replace('nan','')
			print(newline)
			template_pro[i]= newline
	#print(template_pro)
	template_pro = ''.join(template_pro)
	print(template_pro)
	proj_manual.write(template_pro)
	proj_manual.close()



##Upload the manuals to s3

for index, row in pro.iterrows():
	manual_path = "manuals/%s.%s.manual" %(row["#projectIDs"],date)
	project = row["#projectIDs"]
	if re.search('in|zr', project):
		s3_path = "s3://epiquest/epiquest_%s/%s/manual/%s.%s.manual" %(row["#projectIDs"],row["random_strings"],row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	if re.search('m', project):
		s3_path = "s3://midog/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	elif project.startswith('pbe'):
		s3_path = "s3://precisionbiome/PrecisionBIOME_Ear/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	elif project.startswith('pbu'):
		s3_path = "s3://precisionbiome/PrecisionBIOME_Urine/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	elif project.startswith('pbv'):
		s3_path = "s3://precisionbiome/PrecisionBIOME_Vaginal/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	elif project.startswith('pbs'):
		s3_path = "s3://precisionbiome/PrecisionBIOME_Sputum/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	elif project.startswith('pbt'):
		s3_path = "s3://precisionbiome/PrecisionBIOME_Throat/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	elif project.startswith('pbw'):
		s3_path = "s3://precisionbiome/PrecisionBIOME_Wound/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	elif project.startswith('pbna'):
		s3_path = "s3://precisionbiome/PrecisionBIOME_Nail/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
	elif project.startswith('pbno'):
		s3_path = "s3://precisionbiome/PrecisionBIOME_Nose/Projects/%s/manuals/%s.%s.manual" %(row["#projectIDs"], row["#projectIDs"], date)
		os.system("aws s3 cp %s %s" %(manual_path, s3_path))
