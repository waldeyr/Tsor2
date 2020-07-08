'''
What does this program do?
	It merges the resulting files temp01.tsv (annotations plus categories) and temp02.tsv (signalp results)

Input:
    temp01.tsv
    temp02.tsv

Output:
	Tsor2_final_coding_contigs_annotation.tsv

@author: Waldeyr Mendes Cordeiro da Silva, June 2020
'''
import pandas as pd

#Read files
temp01 =  pd.read_csv('temp01.tsv', sep="\t")
temp02 = pd.read_csv('temp02.tsv', sep="\t")

#Join annotations
annotated_contigs = merged_inner = pd.merge(left=temp01, right=temp02, left_on='contig_id', right_on='contig_id', how='outer')
print(annotated_contigs.shape)
print(annotated_contigs.columns)

result = annotated_contigs.to_csv("Tsor2_final_coding_contigs_annotation.tsv", sep='\t', encoding='utf-8', index=False)
print("Done.\n")