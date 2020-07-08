'''
Annotation pipeline for Tsor2 data (~27k assembled contigs)

What does this program do?
	Parses the results from Tsor2 contigs blasted against the uniprot swissprot (protein)
	Select the best hit for each Tsor2 contig based on the percent of identity

Input:
	Tsor2_trembl_blast.result edited with a headline:
		contig_id
		Unirptot_Refseq_len
		Uniprot_Refseq_start
		Uniprot_Refseq_stop
		Uniprot_Refseq_query_seq
		Uniprot_Refseq_idt
		Uniprot_Refseq_percentage_idt
		Uniprot_Refseq_e-value
		Uniprot_Refseq_positive
		Uniprot_Refseq_percentage_pos
		Uniprot_Refseq_gaps
		Uniprot_Refseq_query_frame
		Uniprot_Refseq_acession

Output:
	Tsor2_uniprot_blast_best_hits, a .tsv file (tab separated file) with one contig per line and the best hit on the columns

@author: Waldeyr Mendes Cordeiro da Silva, June 2020
'''

# needed python libraries
# https://pandas.pydata.org
import pandas as pd
# https://github.com/biopython/DIST
from Bio import Entrez
# https://docs.python.org/3/library/re.html
import requests

# reading the file with refseq blast results
Tsor2_trembl_blast_df = pd.read_csv('Tsor2_trembl_blast.result', sep="\t", error_bad_lines=False, low_memory=False)
# selecting the desired columns from the Blast results
print("Working on selection of the desired columns from the Blast results...\n")
df_temp = Tsor2_trembl_blast_df.apply(lambda row: pd.Series(
	[
	 row['contig_id'],
	 row['sseqid'],
	 row['qseq'].replace("-", ""),  # translated query sequence without blast alignment inserted gaps
	 row['qlen'],
	 row['qframe'],
	 row['qstart'],
	 row['qend'],
	 row['nident'],
	 "{:.2f}".format(row['pident']),
	 row['evalue']
	 ],
	index=[
		'contig_id',
		'sseqid',
		'qseq',
		'qlen',
		'qframe',
		'qstart',
		'qend',
		'nident',
		'pident',
		'evalue'
	]), axis=1)
print("Done\n")
# selecting the best hit based on the percent of identity
print("Working on selection of best Blast results...\n")
df_sorted = df_temp.sort_values('pident', ascending=False).drop_duplicates(['contig_id'])
# write .tsv file with best results
df = df_sorted.to_csv("Tsor2_trembl_blast_best_hits.tsv", sep='\t', encoding='utf-8', index=False)
print("Done parsing.\n")

df = pd.read_csv('Tsor2_trembl_blast_best_hits.tsv', sep="\t", error_bad_lines=False, low_memory=False)
print(df.shape)
print(df.columns)
print("Working on annotation for the best Blast results...\n")
# list of protein names for each hit
proteinNamesList = []
# list of taxonomy species for each hit
taxonomyList = []
# providing an e-mail is mandatory to use the NCBI services
Entrez.email = "mendes@iscb.org"
for i in range(len(df)):
    row = str(df.iloc[i, 1])
    if (row.find('tr') != -1):
        query = str(row.split("|")[1])
        url = "https://www.uniprot.org/uniprot/?query=accession:" + query.strip() + "&format=tab"
        response = requests.get(url)
        if response.status_code == 200:
            result = response.text
            proteinNamesList.append(str(result.split("\n")[1].split("\t")[3]))
            taxonomyList.append(str(result.split("\n")[1].split("\t")[5]))
            print(str(result.split("\n")[1].split("\t")[3]) + "|" + str(result.split("\n")[1].split("\t")[5]))
        elif response.status_code == 404:
            proteinNamesList.append("Invalid Uniprot ID")
            taxonomyList.append("None")
            print("#", end='')
        else:
            proteinNamesList.append("None")
            taxonomyList.append("None")
            print("@", end='')
            pass
    else:
        proteinNamesList.append("None")
        taxonomyList.append("None")

print("Done\n")

print("Working on creating a file with compiled results...\n")
# join annotation to the original data
df.insert(len(df.columns), "trembl_annotation", proteinNamesList, True)
# join taxonomy to he original data
df.insert(len(df.columns), "trembl_taxonomy", taxonomyList, True)
# save annotations in a .tsv file
result = df.to_csv("Tsor2_trembl.tsv", sep='\t', encoding='utf-8', index=False)
print("Done\n")
