import pandas as pd

#Read files
curated =  pd.read_csv('Tsor2_macrogen_plus_refseq_plus_uniprot.tsv', sep="\t")
trembl = pd.read_csv('Tsor2_trembl.tsv', sep="\t")

#Join annotations
curated_plus_trembl = merged_inner = pd.merge(left=curated, right=trembl, left_on='contig_id', right_on='contig_id', how='outer')
# print(curated_plus_trembl.shape)
# print(curated_plus_trembl.columns)

curated_plus_trembl.columns = [
       'contig_id', 
       'read_count', 
       'fpkm', 
       'macrogen_go', 
       'macrogen_uniprot', 
       'macrogen_nr',
       'macrogen_pfam', 
       'macrogen_eggnog', 
       'macrogen_nt', 
       'macrogen_ko_euk', 
       'refseq_sseqid', 
       'refseq_qseq', 
       'refseq_qlen',
       'refseq_qframe', 
       'refseq_qstart', 
       'refseq_qend', 
       'refseq_nident', 
       'refseq_pident', 
       'refseq_evalue',
       'refseq_annotation', 
       'refseq_taxonomy', 
       'uniprot_sseqid', 
       'uniprot_qseq', 
       'uniprot_qlen',
       'uniprot_qframe', 
       'uniprot_qstart', 
       'uniprot_qend', 
       'uniprot_nident', 
       'uniprot_pident', 
       'uniprot_evalue',
       'uniprot_annotation', 
       'uniprot_taxonomy', 
       'trembl_sseqid', 
       'trembl_qseq', 
       'trembl_qlen',
       'trembl_qframe', 
       'trembl_qstart', 
       'trembl_qend', 
       'trembl_nident', 
       'trembl_pident', 
       'trembl_evalue',
       'trembl_annotation', 
       'trembl_taxonomy']
#del curated_plus_trembl['del1']
# print(curated_plus_trembl.shape)
# print(curated_plus_trembl.columns)
result = curated_plus_trembl.to_csv("Tsor2_macrogen_plus_refseq_plus_uniprot_trembl.tsv", sep='\t', encoding='utf-8', index=False)
print("Done.\n")