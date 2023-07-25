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

args = parser.parse_args()

projectid = args.folder.split('.')[0]
input = pd.read_csv(args.input, sep=',', header=0, encoding="iso-8859-1")
input['sample_id'] = input['projectID'].astype(str) + '_' + input['#num'].astype(str) 
input['customer_label'] = input['UniqueLabel']
input['analysis'] =  input['GroupID'].astype(str) + "." + input['SeqType']
input = input.drop (labels = ['#num', 'projectID','RunID','SeqType', 'UniqueLabel', 'GroupID'], axis=1)
ColumnsList =input.columns.tolist()
ColumnsList = ColumnsList[-3:] + ColumnsList[:-3]
input = input[ColumnsList]

"""
if args.rundate != None and args.randomstr != None:
	for i,r in input.iterrows():
		input.at[i,'Read1']="<a href=\"https://epiquest.s3.amazonaws.com/epiquest_%s/%s/rawdata/%s/%s_R1.fastq.gz\"target='_blank'>RawRead1</a>" %(projectid, args.randomstr, args.rundate, input.at[i,'sample_id'])
		input.at[i, 'Read2'] = "<a href=\"https://epiquest.s3.amazonaws.com/epiquest_%s/%s/rawdata/%s/%s_R2.fastq.gz\"target='_blank'>RawRead2</a>" %(projectid, args.randomstr, args.rundate, input.at[i,'sample_id'])
else:
	pass
"""

#print(ColumnsList)
#print(input[ColumnsList])

for value in input.analysis.unique():
	group = input[input["analysis"]== value]
	#print(value)
	group.to_csv("%s/%s/SampleMetadata.csv" %(args.folder, value),index=False)
	#print(folder)
