#!/usr/bin/python
__author__= 'Pedro J. Torres'
import argparse
import os
import pandas,numpy
import sys
#import shutil

#-----------Command Line Arguments-----------------
parser=argparse.ArgumentParser(description="Script allows you to convert your diamond output into a normalized relative abundance with percent uniqueness to protein and bacterial hits.")
parser.add_argument('-i','--input', help='Input txt file output from diamond',required=True)
parser.add_argument('-o','--out', help='Name of output file: jsut the mae e.g., "sample1_final"', required=True)
parser.add_argument('-db','--db', help='Database type, "e.g., -db atb"', required=True)

args = parser.parse_args()
o_file=str(args.out)
dbtype=str(args.db)
inputfile=str(args.input)

###########################################################################################################
#--This is part 1 of the script where we reformat and get unique bacteria and genes based on query id-----
###########################################################################################################
#this checks that there were results from the DIAMOND output,if not a blank file will be produced with no error
with open(inputfile) as friendsfile:
        first = friendsfile.read(1)
        if not first:
                print(inputfile+' is empty')
                df = pandas.read_table(inputfile, sep='\t',
                  names = ["Gene","Species","Percent_Unique_Protein","Percent_Unique_Bacteria","Hits","Average_Gene_Length","Relative_abundance"])
                df.to_csv(o_file+'.txt', sep='\t', index=False)
                sys.exit()
        else:pass


#take in file and read args will be above
df = pandas.read_table(inputfile, sep='\t',
                  names = ["Id_VF_Bac", "qseqid", "length"])

# set the index to be this and don't drop
df.set_index(keys=['qseqid'], drop=False,inplace=True)

#Set up dictionary and list for later use
gene_i={}
info=[]

#---- Reformatting the DIAMOND output to be able to manipulate it easier -----------

# will split the datafrane Id_VF_Bac into the qseqid and gene
df_id_vf_bac=df['Id_VF_Bac'].str.split(' ',1, expand=True).rename(columns={0:'ID', 1:'VF_Bacteria' })

# will split VF_Bacteria column into gene and bacteria column now
df_vf_bac=df_id_vf_bac['VF_Bacteria'].str.split('[',1, expand=True).rename(columns={0:'VF', 1:'Bacteria'})

#below will split the Bacteria column into multiple columns depending on how much info there is on the organism (e.g., G,,S, strain, ect..)
df_id_vf_bac=df_vf_bac['Bacteria'].str.split(' ',0, expand=True).rename(columns={0:'Genus', 1:'Species' })

# antibiotic resistant database is formatted differently than virulence factor. The following will remove some unwanted characters while
#  adding 'Spp.' if there is no species name present
if dbtype=='atb':
    df_id_vf_bac["Species"].fillna(value='Spp.', inplace=True)
    df_id_vf_bac['Genus'] = df_id_vf_bac['Genus'].str.replace(']',' ')
    df_id_vf_bac['Species'] = df_id_vf_bac['Species'].str.replace(']',' ')
else:pass

#concatinate columns/variables which are of use.
df_id_vf_bac['Bacteria']=df_id_vf_bac["Genus"].map(str) + " " + df_id_vf_bac["Species"] # this drops stuff
df_id_vf_bac['Bacteria']=df_id_vf_bac['Bacteria'].astype(str)

#concatinate ourVirulence factor, and new bacteria columns
df_qseqid_vf_bacteria=pandas.concat([df_vf_bac['VF'], df_id_vf_bac['Bacteria'], df['length']], axis=1)
df_qseqid_vf_bacteria=df_qseqid_vf_bacteria.reset_index(level=['qseqid']) # so we can pick out the bacteria column must reset its index

df_qseqid_vf_bacteria=df_qseqid_vf_bacteria.groupby('qseqid') #group by qseq id
df=[df_qseqid_vf_bacteria.get_group(x) for x in df_qseqid_vf_bacteria.groups]#split each individual qseqid into a df

print ('Reformatting complete')

#start gettig percent unique hits for gene and bacteria for each unique qseqid or raw fasta line
for qseqid in df:
    qid=qseqid[qseqid.columns[0]].iloc[0].strip()
    length=qseqid[qseqid.columns[3]].iloc[0]
    legnthofdataframe = len(qseqid.index)
    if len(qseqid.index) == 1:
        top_bacteria=qseqid[qseqid.columns[2]].iloc[0]
        qid=qseqid[qseqid.columns[0]].iloc[0]
        gene=qseqid[qseqid.columns[1]].iloc[0]
        length=qseqid[qseqid.columns[3]].iloc[0]
        merged=[gene,top_bacteria,str(100.0),str(100.0),str(length)]
        gene_i[qid]=merged
    else:
        qid=qseqid[qseqid.columns[0]].iloc[0]
        col=qseqid.groupby(["Bacteria"]).size().reset_index(name="Count")
        col=col.sort_values(by=['Count'], ascending=False) #after we sort based on most abundant bacteria then lets capture that in a list
        col['Percent_Unique_Hit_Bacteria']=(col['Count']/col['Count'].sum()*100)
        if col.empty==True:
            continue
        top_bacteria= col[col.columns[0]].iloc[0].strip()
        bacteria_uniq=col[col.columns[2]].iloc[0]
        col=qseqid.groupby(["VF"]).size().reset_index(name="Count")
        col=col.sort_values(by=['Count'], ascending=False) #after we sort based on most abundant bacteria then lets capture that in a list
        col['Percent_Unique_Hit_Gene']=(col['Count']/col['Count'].sum()*100)
        gene= col[col.columns[0]].iloc[0].strip()
        percent_uniq_gene=col[col.columns[2]].iloc[0]
        #add above info to list and then add that to dictionary where your query seq id is your key and results to that hit are your values in a list
        merged=[gene,top_bacteria,str(percent_uniq_gene),str(bacteria_uniq),str(length)]
        gene_i[qid]=merged

#write out file this could very well end up being a tmp file
#newpath=os.getcwd()
o=open("percent_abundanceatbtmp1.txt","w+")# will change this later
o.write("qseqid"+"\t"+"Gene"+"\t"+"Organism"+"\t"+"Percent_Unique_Hits_Gene" +"\t"+"Percent_Unique_Hits_Bacteria"+"\t"+"Gene_length"+ "\n")
for i in gene_i:
    o.write(i+"\t"+"\t".join([(xx) for xx in gene_i[i]])+"\n")
o.close()
print('Finished up part 1/3. Will now calculate total percent unique hits.')
########################################################################################################################
#---This is part 2 of the script where we start grouping based on genes and species abudnance in whole file--------------
########################################################################################################################
print('Starting part 2/3')
df = pandas.read_table("percent_abundanceatbtmp1.txt", sep='\t')

#dont care about qseqid anymore so we will drop that
df.drop(df.columns[[0]], axis=1, inplace=True)

#merge organism and gene column into one
df['Organism_Gene']=df['Organism'].map(str) + ";"+ df['Gene']

# set the index to be this and don't drop. Everything will revovle around the virulence gene as this is what we
# want to get abundace to as before
df.set_index(keys=['Organism_Gene'], drop=False,inplace=True)

# get a list of unique genes present in the directory
names=df['Organism_Gene'].unique().tolist()
print ('There are ' + str(len(names))+ ' unique organism/gene combinations')

# list and dictionary to add stuff to later
taxa={}
information=[]

df=df.groupby('Organism_Gene') #group by organism gene
df=[df.get_group(x) for x in df.groups]#split each individual qseqid into a df
#start getting ubique abundance of gene in the file
for org_gene in df:
    if len(org_gene.index) == 1:
        gene=org_gene[org_gene.columns[0]].iloc[0].strip()
        bacteria_test = org_gene[org_gene.columns[1]].iloc[0]
        if bacteria_test != bacteria_test:
                continue
        bacteria=org_gene[org_gene.columns[1]].iloc[0].strip()
        percnt_uniq_gene=org_gene[org_gene.columns[2]].iloc[0]
        percnt_uniq_bacteria=org_gene[org_gene.columns[3]].iloc[0]
        raw_hits=1.0
        average_length=org_gene[org_gene.columns[4]].iloc[0]
        normalized_abundance=(1.0/average_length)*100.0
        information=[bacteria,str(percnt_uniq_gene+0.0),str(percnt_uniq_bacteria+0.0),str(raw_hits),str(normalized_abundance),str(average_length+0.0)]
        information=[bacteria,str(percnt_uniq_gene),str(percnt_uniq_bacteria),str(raw_hits),str(normalized_abundance),str(average_length)]
        taxa[gene]=information
    else:
        pandas.options.mode.chained_assignment = None  # default='warn'
        gene=org_gene[org_gene.columns[0]].iloc[0].strip()
        bacteria=org_gene[org_gene.columns[1]].iloc[0]
        #Percent unique htis bacteria
        df_u_gene_per=org_gene[['Organism_Gene','Percent_Unique_Hits_Bacteria','Gene_length']]
        df_u_gene_per['Organism_Gene_Raw_Hits']=1
        col=df_u_gene_per.groupby(['Organism_Gene'])[["Percent_Unique_Hits_Bacteria","Organism_Gene_Raw_Hits"]].sum()
        col['Percent_Unique_Hits_Organism']=(col['Percent_Unique_Hits_Bacteria']/col['Organism_Gene_Raw_Hits'])
        col['Avg_len']=(df_u_gene_per['Gene_length'].sum()/col['Organism_Gene_Raw_Hits'])
        col['Normalized_count']=col['Organism_Gene_Raw_Hits']/col['Avg_len']*100
        # start organizing
        raw_hits= col[col.columns[1]].iloc[0]
        average_length=col[col.columns[3]].iloc[0] # this is the average length of he particular gene that was identified via blast
        normalized_abundance=col[col.columns[4]].iloc[0]
        col=col.reset_index(level=['Organism_Gene']) # so we can pick out the bacteria column must reset its index
        percnt_uniq_bacteria = col[col.columns[3]].iloc[0]
        # percent unique protein
        df_u_gene_per=org_gene[['Organism_Gene','Percent_Unique_Hits_Gene','Percent_Unique_Hits_Bacteria']]
        df_u_gene_per['Organism_Gene_Raw_Hits']=1
        col=df_u_gene_per.groupby(['Organism_Gene'])[["Percent_Unique_Hits_Gene","Organism_Gene_Raw_Hits"]].sum()
        col['Percent_Unique_Hits_GeneFinal']=(col['Percent_Unique_Hits_Gene']/col['Organism_Gene_Raw_Hits'])
        percnt_uniq_gene = col[col.columns[2]].iloc[0]
        #add above info to list and then add that to dictionary where your query seq id is your key and results to that hit are your values in a list
        information=[bacteria,str(percnt_uniq_gene),str(percnt_uniq_bacteria),str(raw_hits),str(normalized_abundance),str(average_length)]
        taxa[gene]=information

#write out file this could very well end up being a tmp file
#workdir=os.getcwd()
o=open("percent_abundanceatbtmp2.txt","w+")
o.write("Gene"+"\t"+"Species"+"\t"+"Percent_Unique_Protein"+"\t"+"Percent_Unique_Bacteria"+"\t"+"Hits" +"\t"+"Normalized_Abundance"+"\t"+"Average_Gene_Length"+ "\n")
for i in taxa:
    try:
        o.write(i+"\t"+"\t".join([(xx) for xx in taxa[i]])+"\n")
    except TypeError:
        continue
o.close()
print ('Finished part 2/3. Will now calculate gene relative abundance.')
########################################################################################################################
#-- This is the last part of the script pt 3 where we will get the relative abudnance of the gene in our sample ----
########################################################################################################################
df = pandas.read_table("percent_abundanceatbtmp2.txt", sep='\t')

# set the index to be this and don't drop. Everything will revovle around the virulence gene as this is what we
# want to get abundace to as before
df.set_index(keys=['Gene'], drop=False,inplace=True)
df['Relative_abundance']=(df['Normalized_Abundance'].div(df['Normalized_Abundance'].sum(axis=0)).multiply(100))

#lets drop unimportant columns in this case the normalized relative abudance
df.drop(df.columns[[5]], axis=1, inplace=True)
df.drop(df.columns[[0]], axis=1, inplace=True)
df=df.sort_values(by=['Species'], ascending=True)
df.to_csv(o_file+'.txt', sep='\t')

#shutil.rmtree(newpath)
os.remove('percent_abundanceatbtmp1.txt')
os.remove('percent_abundanceatbtmp2.txt')
print ('Finished part 3/3.')
print ('Done :)')
