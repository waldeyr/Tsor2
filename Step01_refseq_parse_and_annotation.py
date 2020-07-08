'''
Annotation pipeline for Tsor2 data (~27k assembled contigs)

What does this program do?
	Inserts a headline to the Blast result file Tsor2_refseq_blast.result
	Parses the results from Tsor2 contigs blasted against the complete refseq (protein) release 200
	Select the best hit for each Tsor2 contig based on the percent of identity

Input:
	Tsor2_refseq_blast.result edited with a headline: contig_id	qgi	qacc	qlen	qstart	qend	qseq	nident	pident	evalue	positive	ppos	gaps	qframe	sseqid

Output:
	Tsor2_refseq_blast_best_hits, a .tsv file (tab separated file) with one contig per line and the best hit on the columns

@author: Waldeyr Mendes Cordeiro da Silva, June 2020
'''

# needed python libraries
# https://pandas.pydata.org
import pandas as pd
# https://github.com/biopython/DIST
from Bio import Entrez
# https://docs.python.org/3/library/re.html
import re

# reading the file with refseq blast results
Tsor2_refseq_blast_df = pd.read_csv('Tsor2_refseq_blast.result', sep="\t", error_bad_lines=False, low_memory=False)
# selecting the desired columns from the Blast results
print("Working on selection of the desired columns from the Blast results...\n")
df_temp = Tsor2_refseq_blast_df.apply(lambda row: pd.Series(
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
df = df_sorted.to_csv("Tsor2_refseq_blast_best_hits.tsv", sep='\t', encoding='utf-8', index=False)
print("Done parsing.\n")

print("Working on annotation for the best Blast results...\n")
df = pd.read_csv('Tsor2_refseq_blast_best_hits.tsv', sep="\t", error_bad_lines=False, low_memory=False)
# print(df.shape)
# print(df.columns)
# list of protein names for each hit
proteinNamesList = []
# list of taxonomy species for each hit
taxonomyList = []
# providing an e-mail is mandatory to use the NCBI services
Entrez.email = "mendes@iscb.org"
# For each best hit on results, get the annotation and the taxonomy from the NCBI
for i in range(len(df)):
	id_ncbi_for_query = str(df.iloc[i, 1])
	try:
		handle = Entrez.esummary(db="protein", id=str(id_ncbi_for_query))
		if (handle):
			query_info = Entrez.read(handle)
			print(str(df.iloc[i, 0]) + "|"+query_info[0]['Title'])
			if (query_info):
				match = re.search(r"(\[.*\])", str(query_info[0]['Title']), re.IGNORECASE)
				if (match):
					taxonomy = str(match.group()).rstrip("]").lstrip("[")
					annotation = str(query_info[0]['Title']).replace(str(match.group()), "").rstrip().lstrip()
				else:
					taxonomy = "None"
					annotation = str(query_info[0]['Title']).rstrip().lstrip()
			else:
				taxonomy = "None"
				annotation = "None"
		else:
			annotation = "None"
			taxonomy = "None"
	except Exception as e:
		print(e)
		annotation = "None"
		taxonomy = "None"
	proteinNamesList.append(annotation)
	taxonomyList.append(taxonomy)

print("Done\n")

print("Working on creating a file with compiled results...\n")
# join annotation to the original data
df.insert(len(df.columns), "refseq_annotation", proteinNamesList, True)
# join taxonomy to he original data
df.insert(len(df.columns), "taxonomy", taxonomyList, True)
# save annotations in a .tsv file
result = df.to_csv("Tsor2_refseq.tsv", sep='\t', encoding='utf-8', index=False)
print("Done\n")