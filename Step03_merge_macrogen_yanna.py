'''
This step join refseq with previous Macrogen annotation

What does this program do?
       It merges Macrogen and refseq annotations

Input:
       Tsor2_macrogen.tsv
       Tsor2_refseq.tsv
       Tsor2_uniprot.tsv

Output:
       Tsor_macrogen_plus_refseq_plus_uniprot.tsv

@author: Waldeyr Mendes Cordeiro da Silva, June 2020
'''

import pandas as pd

#Read files
macrogen =  pd.read_csv('Tsor2_macrogen.tsv', sep="\t")
refseq = pd.read_csv('Tsor2_refseq.tsv', sep="\t")
uniprot = pd.read_csv('Tsor2_uniprot.tsv', sep="\t")

#Join annotations
macrogen_plus_refseq = merged_inner = pd.merge(left=macrogen, right=refseq, left_on='contig_id', right_on='contig_id', how='outer')
print(macrogen_plus_refseq.shape)
print(macrogen_plus_refseq.columns)
macrogen_plus_refseq.columns = ['contig_id', 'Tsor2_Read_Count', 'Tsor2_FPKM', 'GO', 'UniProt', 'NR',
       'Pfam', 'EggNOG', 'NT', 'KO_EUK', 'del1', 'sseqid', 'qseq',
       'qlen', 'qframe', 'qstart', 'qend', 'nident', 'pident', 'evalue',
       'refseq_annotation', 'taxonomy']
del macrogen_plus_refseq['del1']
print(macrogen_plus_refseq.shape)
print(macrogen_plus_refseq.columns)

macrogen_plus_refseq_plus_uniprot = merged_inner = pd.merge(left=macrogen_plus_refseq, right=uniprot, left_on='contig_id', right_on='contig_id', how='outer')

print(macrogen_plus_refseq_plus_uniprot.shape)
print(macrogen_plus_refseq_plus_uniprot.columns)
result = macrogen_plus_refseq_plus_uniprot.to_csv("Tsor2_macrogen_plus_refseq_plus_uniprot.tsv", sep='\t', encoding='utf-8', index=False)
print("Done.\n")