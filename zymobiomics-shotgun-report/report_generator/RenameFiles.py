#!usr/bin/python3.6
__author__= 'Aishani Prem'
__email__='aprem@zymoresearch.com'
import re
import argparse
import os
import pandas as pd

#---Command Line Arguments----
parser=argparse.ArgumentParser(
      description="This script uses python3. Please make sure you have the correct version installed. This script is used to convert the services output folder to a customer version.")

parser.add_argument('-i','--input', help='Input the name of the output folder',required=True)
args = parser.parse_args()


folders= os.listdir(args.input)


for folder in folders:
	if folder.endswith("illumina.pe") or folder.endswith("illumina.se"):
		os.chdir("%s/%s/"%(args.input,folder))
		#Delete unwanted files
		os.system("rm -r *.mapping.file.txt")
		os.system("rm -r sample.list.txt")
		#Renaming folders
		os.system("mv virulence_diamond VirulenceFactor")
		os.system("mv AMR_diamond AntibioticResistance")
		os.system("mv virus Virus")
		os.system("mv humann FunctionalPathway")		
		os.system("mv eukaryote Eukaryote")
		os.system("mv prokaryote Prokaryote")
		os.system("mv all All")
		#Abundance Tables
		os.system("mkdir AbundanceTables")
		#file = pd.read_csv("abun_table.tsv", sep='\t', header=0)
		#file.to_csv("AbundanceTables/AbundanceTable.csv" ,index=False)
		os.system("rm abun_table.tsv")
		file = pd.read_csv("host_dna_abun.csv", sep='\t', header=0)
		file.to_csv("AbundanceTables/ReadDistributionTable.csv" ,index=False)
		os.system("rm host_dna_abun.csv")
		#Sample information folder
		os.system("mkdir SampleInformation")
		os.system("mv SampleMetadata.csv SampleInformation/")
		os.system("cp -r SampleInformation/ FunctionalPathway/")
		os.system("cp -r SampleInformation/ All/")
		os.system("cp -r SampleInformation/ Eukaryote")
		os.system("cp -r SampleInformation/ Prokaryote")
		os.system("cp -r SampleInformation/ Virus")
		ReadProcessing = pd.read_csv("Trimmomatic/summary.tsv", sep='\t', header=0)
		ReadProcessing['both_surviving(%)'] = ReadProcessing['both_surviving(%)'].astype(str).str.strip("%")
		ReadProcessing['both_surviving(%)'] = ReadProcessing['both_surviving(%)'].astype(float)
		ReadProcessing['forward_only(%)'] = ReadProcessing['forward_only(%)'].astype(str).str.strip("%")
		ReadProcessing['forward_only(%)'] = ReadProcessing['forward_only(%)'].astype(float)
		ReadProcessing['reverse_only(%)'] = ReadProcessing['reverse_only(%)'].astype(str).str.strip("%")
		ReadProcessing['reverse_only(%)'] = ReadProcessing['reverse_only(%)'].astype(float)
		ReadProcessing["TM Surviving(%)"] = ReadProcessing[['both_surviving(%)', 'forward_only(%)', 'reverse_only(%)']].sum(axis=1).astype(str)
		ReadProcessing["TM Dropped(%)"] =  ReadProcessing['dropped(%)']
		ReadProcessing['RNAs Removed(%)'] = ReadProcessing['rRNAs_removed(%)']
		ReadProcessing['Final Reads(%)'] = ReadProcessing['final_reads(%)']
		ReadProcessing = ReadProcessing.drop(labels = ['both_surviving(%)','forward_only(%)', 'reverse_only(%)', 'dropped(%)', 'rRNAs_removed(%)', 'final_reads(%)'], axis= 1) 
		ReadProcessing.to_csv("SampleInformation/ReadProcessingSummary.csv"  ,index=False)
		os.system("rm -r Trimmomatic/")
		os.system("mv FastQC SampleInformation/")
		#Files in FunctionalPathway folder
		os.system("mkdir FunctionalPathway/RawData")
		os.system("mv FunctionalPathway/*tsv FunctionalPathway/RawData/")
		dirs = os.listdir("FunctionalPathway")
		for dir in dirs:
			if dir.startswith("heatmap_species_pathway_abun"):
				newname = re.sub("heatmap_species_pathway_abun", "Heatmap_SpeciesPathwayAbundance", dir)
				os.system("mv FunctionalPathway/%s FunctionalPathway/%s" %(dir,newname))
				os.system("mv FunctionalPathway/%s/new_abun_table.tsv FunctionalPathway/%s/Raw_Data.tsv" %(newname,newname))
				os.system("mv FunctionalPathway/%s/heatmap_with_sample_clustering.pdf FunctionalPathway/%s/Heatmap_with_SampleClustering.pdf" %(newname,newname))
				os.system("mv FunctionalPathway/%s/heatmap_without_sample_clustering.pdf FunctionalPathway/%s/Heatmap_without_SampleClustering.pdf" %(newname,newname))
			elif dir.startswith("heatmap_pathway_abun"):
				newname = re.sub("heatmap_pathway_abun", "Heatmap_PathwayAbundance", dir)
				os.system("mv FunctionalPathway/%s FunctionalPathway/%s" %(dir,newname))
				os.system("mv FunctionalPathway/%s/new_abun_table.tsv FunctionalPathway/%s/Raw_Data.tsv" %(newname,newname))
				os.system("mv FunctionalPathway/%s/heatmap_with_sample_clustering.pdf FunctionalPathway/%s/Heatmap_with_SampleClustering.pdf" %(newname,newname))
				os.system("mv FunctionalPathway/%s/heatmap_without_sample_clustering.pdf FunctionalPathway/%s/Heatmap_without_SampleClustering.pdf" %(newname,newname))
			elif dir.startswith("heatmap_gene_fam_cpm"):
				newname = re.sub("heatmap_gene_fam_cpm", "Heatmap_GeneFamily_CPM", dir)
				os.system("mv FunctionalPathway/%s FunctionalPathway/%s" %(dir,newname))
				os.system("mv FunctionalPathway/%s/new_abun_table.tsv FunctionalPathway/%s/Raw_Data.tsv" %(newname,newname))
				os.system("mv FunctionalPathway/%s/heatmap_with_sample_clustering.pdf FunctionalPathway/%s/Heatmap_with_SampleClustering.pdf" %(newname,newname))
				os.system("mv FunctionalPathway/%s/heatmap_without_sample_clustering.pdf FunctionalPathway/%s/Heatmap_without_SampleClustering.pdf" %(newname,newname))
			elif dir.startswith("lefse_gene_fam_cpm"):
				newname = re.sub("lefse_gene_fam_cpm", "LEfSe_GeneFamily_CPM", dir)
				os.system("mv FunctionalPathway/%s FunctionalPathway/%s" %(dir,newname))
				os.system("mv FunctionalPathway/%s/biomarkers FunctionalPathway/%s/Figures" %(newname,newname))
				os.system("rm FunctionalPathway/%s/map.txt" %(newname))
				os.system("rm FunctionalPathway/%s/lefse.in" %(newname))
				os.system("rm FunctionalPathway/%s/abun.txt" %(newname))
				os.system("mv FunctionalPathway/%s/lefse.input.txt  FunctionalPathway/%s/LEfSe_Input.txt" %(newname,newname))
				os.system("mv FunctionalPathway/%s/biomarkers.pdf FunctionalPathway/%s/Biomarkers.pdf" %(newname,newname))
				os.system("mv FunctionalPathway/%s/cladogram.pdf FunctionalPathway/%s/Cladogram.pdf" %(newname,newname))
				os.system("mv FunctionalPathway/%s/lefse.res.xls FunctionalPathway/%s/LEfSe_Results.tsv" %(newname,newname))
			elif dir.startswith("lefse_pathway_abun"):
				newname = re.sub("lefse_pathway_abun", "LEfSe_SpeciesPathwayAbundance", dir)
				os.system("mv FunctionalPathway/%s FunctionalPathway/%s" %(dir,newname))
				os.system("mv FunctionalPathway/%s/biomarkers FunctionalPathway/%s/Figures" %(newname,newname))
				os.system("rm FunctionalPathway/%s/map.txt" %(newname))
				os.system("rm FunctionalPathway/%s/lefse.in" %(newname))
				os.system("rm FunctionalPathway/%s/abun.txt" %(newname))
				os.system("mv FunctionalPathway/%s/lefse.input.txt  FunctionalPathway/%s/LEfSe_Input.txt" %(newname,newname))
				os.system("mv FunctionalPathway/%s/biomarkers.pdf FunctionalPathway/%s/Biomarkers.pdf" %(newname,newname))
				os.system("mv FunctionalPathway/%s/cladogram.pdf FunctionalPathway/%s/Cladogram.pdf" %(newname,newname))
				os.system("mv FunctionalPathway/%s/lefse.res.xls FunctionalPathway/%s/LEfSe_Results.tsv" %(newname,newname))
			else:
				pass
		#Rename folders for AMR and Virulence factor
		os.system("mv AntibioticResistance/summary AntibioticResistance/Summary")
		if len(os.listdir("AntibioticResistance/Summary")) > 0:
			ResultFiles =  os.listdir("AntibioticResistance/Summary")
			for file in ResultFiles:
				newfile = file.replace('.txt', '.csv')
				file = pd.read_csv("AntibioticResistance/Summary/%s" %file, sep='\t', header=0)
				#file = file.drop(labels = ['Percent_Unique_Protein','Percent_Unique_Bacteria', 'Hits', 'Average_Gene_Length'], axis= 1)
				file.to_csv("AntibioticResistance/Summary/%s" %newfile ,index=False)
			os.system("rm AntibioticResistance/Summary/*txt")
		else:
			os.system("rm -r AntibioticResistance/Summary")
		os.system("mv VirulenceFactor/summary VirulenceFactor/Summary")
		if len(os.listdir("VirulenceFactor/Summary")) > 0:
			ResultFiles =  os.listdir("VirulenceFactor/Summary")
			for file in ResultFiles:
				newfile = file.replace('.txt', '.csv')
				file = pd.read_csv("VirulenceFactor/Summary/%s" %file, sep='\t', header=0)
				#file = file.drop(labels = ['Percent_Unique_Protein','Percent_Unique_Bacteria', 'Hits', 'Average_Gene_Length'], axis= 1)
				file.to_csv("VirulenceFactor/Summary/%s" %newfile ,index=False)
			os.system("rm VirulenceFactor/Summary/*txt")
		else:
			os.system("rm -r VirulenceFactor/Summary")
		#Reanem folders for Analysis files
		analysis = ['All', 'Eukaryote', 'Prokaryote', 'Virus']
		for files in analysis:
			os.system("mv %s/barplots/ %s/CompositionBarplots/" %(files, files))
			os.system("mv %s/abun_tables/ %s/AbundanceTables/" %(files, files))
			#os.system("mv %s/abun_table.tsv %s/AbundanceTables/ReadAbundance.tsv" %(files, files))
			os.system("rm -r %s/alpha_div*" %(files))
			dirs = os.listdir(files)
			for dir in dirs:
				if dir == "AbundanceTables":
					os.system("mv %s/%s/superkingdom %s/%s/1.Superkingdom" %(files, dir, files, dir))
					os.system("mv %s/%s/phylum %s/%s/2.Phylum" %(files, dir, files, dir))
					os.system("mv %s/%s/order %s/%s/3.Order" %(files, dir, files, dir))
					os.system("mv %s/%s/family %s/%s/4.Family" %(files, dir, files, dir))
					os.system("mv %s/%s/genus %s/%s/5.Genus" %(files, dir, files, dir))
					os.system("mv %s/%s/species %s/%s/6.Species" %(files, dir, files, dir))
					os.system("mv %s/%s/strain %s/%s/7.Strain" %(files, dir, files, dir))
				elif dir.startswith('beta_div'):
					newname = re.sub("beta_div", "BetaDiversity", dir)
					os.system("mv %s/%s/genus %s/%s/1.Genus" %(files, dir, files, dir))
					os.system("mv %s/%s/1.Genus/cordinates %s/%s/1.Genus/Coordinates" %(files, dir, files, dir))
					os.system("mv %s/%s/1.Genus/dist %s/%s/1.Genus/DistanceMatrix" %(files, dir, files, dir))
					os.system("mv %s/%s/1.Genus/biplot %s/%s/1.Genus/Biplot" %(files, dir, files, dir))
					os.system("mv %s/%s/species %s/%s/2.Species" %(files, dir, files, dir))
					os.system("mv %s/%s/2.Species/cordinates %s/%s/2.Species/Coordinates" %(files, dir, files, dir))
					os.system("mv %s/%s/2.Species/dist %s/%s/2.Species/DistanceMatrix" %(files, dir, files, dir))
					os.system("mv %s/%s/2.Species/biplot %s/%s/2.Species/Biplot" %(files, dir, files, dir))
					os.system("mv %s/%s/strain %s/%s/3.Strain" %(files, dir, files, dir))
					os.system("mv %s/%s/3.Strain/cordinates %s/%s/3.Strain/Coordinates" %(files, dir, files, dir))
					os.system("mv %s/%s/3.Strain/dist %s/%s/3.Strain/DistanceMatrix" %(files, dir, files, dir))
					os.system("mv %s/%s/3.Strain/biplot %s/%s/3.Strain/Biplot" %(files, dir, files, dir))
					os.system("mv %s/%s %s/%s" %(files, dir, files, newname))
				elif dir.startswith('heatmaps'):
					newname = re.sub("heatmaps", "Heatmaps", dir)
					os.system("mv %s/%s/phylum %s/%s/1.Phylum" %(files, dir, files, dir))
					os.system("mv %s/%s/order %s/%s/2.Order" %(files, dir, files, dir))
					os.system("mv %s/%s/family %s/%s/3.Family" %(files, dir, files, dir))
					os.system("mv %s/%s/genus %s/%s/4.Genus" %(files, dir, files, dir))
					os.system("mv %s/%s/species %s/%s/5.Species" %(files, dir, files, dir))
					os.system("mv %s/%s/strain %s/%s/6.Strain" %(files, dir, files, dir))
					os.system("mv %s/%s %s/%s" %(files, dir, files, newname))
				elif dir.startswith('lefse_'):
					newname = dir.replace("lefse_", "LEfSe_")
					os.system("mv %s/%s %s/%s" %(files, dir, files, newname))
					os.system("mv %s/%s/biomarkers %s/%s/Figures" %(files, newname, files, newname))
					os.system("mv %s/%s/lefse.input.txt  %s/%s/LEfSe_Input.txt" %(files, newname, files, newname))
					os.system("mv %s/%s/biomarkers.pdf %s/%s/Biomarkers.pdf" %(files, newname, files, newname))
					os.system("mv %s/%s/cladogram.pdf %s/%s/Cladogram.pdf" %(files, newname, files, newname))
					try:
						file = pd.read_csv("%s/%s/lefse.res.xls" %(files,newname), sep='\t', header=0)
						file.to_csv("%s/%s/LEfSe_Results.csv" %(files, newname) ,index=False)
					except:
						pass
					os.system("rm %s/%s/lefse.res.xls" %(files, newname))
					os.system("rm %s/%s/map.txt" %(files, newname))
					os.system("rm %s/%s/lefse.in" %(files, newname))
					os.system("rm %s/%s/abun.txt" %(files, newname))
				else:
					pass
			os.system("rm %s/abun_table.tsv" %(files))
			os.system("rm %s/qiime.mapping.file.txt" %(files))
			os.system("rm %s/sample.list.txt" %(files))
		os.chdir("./../../")
		
		
