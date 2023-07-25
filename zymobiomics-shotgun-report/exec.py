#!usr/bin/python3.6
__author__= 'Aishani Prem'
__email__='aprem@zymoresearch.com'
import argparse
import os
import subprocess
import re
import glob

#---Command Line Arguments----
parser=argparse.ArgumentParser(
description="This script uses python3. Please make sure you have the correct version installed. This script is used to execute the report generation pipeline for ZymoBIOMICS Services")
parser.add_argument('-p','--project', help='Input the path to the  project folder',required=True)
parser.add_argument('-m','--metadata', help='Input the path to the metadata folder',required=True)

args = parser.parse_args()

project =  args.project
metadata = args.metadata
os.system("rm -r %s" %(project))
os.system("cp -r visualization %s" %(project))

groupfolders = os.listdir(args.project)

exec_path = os.path.dirname(os.path.realpath(__file__))
pwd = os.getcwd()

os.system("python3 %s/report_generator/SplitAnalysisFile.py -i %s -f %s" %(exec_path, metadata, project))
os.system("python3 %s/report_generator/RenameFiles.py -i %s " %(exec_path, project))

#Run the report page

os.chdir(project)
os.system("cp %s/report_generator/Report.Rmd Report.Rmd" %(exec_path))
os.system("cp %s/report_generator/*png ./" %(exec_path))
os.system("Rscript -e \"rmarkdown::render(\'Report.Rmd\',\'html_document\')\"")
os.system("rm Report.Rmd")
os.system("rm *png")
os.system("rm rawdata_link.txt")
os.chdir(pwd)


#Generate the Group Overview page and results page for the report

TaxonomyFolders = ['All', 'Prokaryote', 'Eukaryote', "Virus"]
for group in groupfolders:
	if group.endswith("illumina.pe") or group.endswith("illumina.se"):
		print ("GENERATING RESULTS PAGE FOR %s" %(group))
		os.system("cp %s/report_generator/GroupOverview.Rmd %s/%s/%s/"  %(exec_path, pwd, project,group))
		os.chdir("%s/%s/%s/" %(pwd, project,group))
		os.system("Rscript -e \"rmarkdown::render(\'GroupOverview.Rmd\',\'html_document\')\"")
		os.system("rm GroupOverview.Rmd")
		os.system("rm -r knitr_tmp")
		subfolders = os.listdir("%s/%s/%s/" %(pwd, project,group))
		for analysis in TaxonomyFolders:
			files = os.listdir("%s/%s/%s/%s/" %(pwd, project, group, analysis))
			if "AbundanceTables" in files:
				os.system("cp %s/report_generator/TaxonomyAnalysis.Rmd %s/%s/%s/%s/"  %(exec_path, pwd, project,group, analysis))
				os.chdir("%s/%s/%s/%s/" %(pwd,project,group,analysis))
				os.system("Rscript -e \"rmarkdown::render(\'TaxonomyAnalysis.Rmd\',\'html_document\')\"")
				os.system("rm TaxonomyAnalysis.Rmd")
				os.system("rm -r knitr_tmp")
			else:
				os.system("cp %s/report_generator/NoneDetected.Rmd %s/%s/%s/%s/"  %(exec_path, pwd, project,group, analysis))
				os.chdir("%s/%s/%s/%s/" %(pwd,project,group,analysis))
				os.system("Rscript -e \"rmarkdown::render(\'NoneDetected.Rmd\',\'html_document\')\"")
				os.system("rm NoneDetected.Rmd")
				os.system("mv NoneDetected.html TaxonomyAnalysis.html")
				os.system("rm -r knitr_tmp")
		#Generate the reports for functional analysis
		os.system("cp %s/report_generator/FunctionalResults.Rmd %s/%s/%s/FunctionalPathway/"  %(exec_path, pwd, project, group))
		os.chdir("%s/%s/%s/FunctionalPathway/" %(pwd,project,group))
		os.system("Rscript -e \"rmarkdown::render(\'FunctionalResults.Rmd\',\'html_document\')\"")
		os.system("rm FunctionalResults.Rmd")
		os.system("rm -r knitr_tmp")
