#!usr/bin/python3.6
__author__= 'Aishani Prem'
__email__='aprem@zymoresearch.com'
import pandas as pd
import re
import numpy as np
import argparse
import os


#---Command Line Arguments----
parser=argparse.ArgumentParser(
      description="This script uses python3. Please make sure you have the correct version installed. This script is used to add the sample metadata table to the analysis folder in the shotgun pipeline")

parser.add_argument('-i','--input', help='Input the path to the analysis file',required=True)
parser.add_argument('-f', '--folder', help='Input the name of the project folder',required=True)
parser.add_argument('-d', '--rundate', help='Input the name of the project folder',required=True)
parser.add_argument('-r', '--randomstr', help='Input the name of the project folder',required=True)
parser.add_argument('-t', '--tag', help='Input the name of the project folder',required=True)

args = parser.parse_args()

projectid = args.folder.split('.')[0]
input = pd.read_csv(args.input, sep=',', header=0, encoding="iso-8859-1")
input['sample_id'] = input['projectID'].astype(str) + '_' + input['#num'].astype(str)
input['customer_label'] = input['UniqueLabel']
input = input[~input.GroupID.str.contains("QC")]

input['analysis'] =  input['GroupID'].astype(str) + "." + input['SeqType']
input = input.drop (labels = ['#num', 'projectID','RunID','SeqType', 'UniqueLabel', 'GroupID'], axis=1)
input = input[['sample_id', 'customer_label']]
input = input.drop_duplicates(subset = ['sample_id', 'customer_label'])


if args.rundate != None and args.randomstr != None:
	for i,r in input.iterrows():
		input.at[i,'Read1 Download'] = "https://zymo-microbiomics-service/epiquest/epiquest_%s/%s/rawdata/%s/%s_R1.fastq.gz" %(projectid, args.randomstr, args.rundate, input.at[i,'sample_id'])
		input.at[i, 'Read2 Download'] = "https://zymo-microbiomics-service/epiquest/epiquest_%s/%s/rawdata/%s/%s_R2.fastq.gz" %(projectid, args.randomstr, args.rundate, input.at[i,'sample_id'])		
input.to_csv("RawdataLinks_%s.csv" %args.tag ,index=False)
os.system("aws s3 cp RawdataLinks_%s.csv s3://zymo-microbiomics-service/epiquest/epiquest_%s/%s/rawdata/%s/RawdataLinks_%s.csv --acl public-read-write" %(projectid, args.randomstr, args.rundate, args.tag))

RawdataFile = open("%s/rawdata_link.txt" % (args.folder),"w")
RawdataFile.write ("https://zymo-microbiomics-service/epiquest/epiquest_%s/%s/rawdata/%s/RawdataLinks_%s.csv" %(projectid, args.randomstr, args.rundate, args.tag))
RawdataFile.close()
#print(input)
