import pandas as pd
from textwrap import wrap
from Bio import SeqIO

#Read files
annotation =  pd.read_csv('Tsor2_macrogen_plus_refseq_plus_uniprot_trembl.tsv', sep="\t")

# print(annotation.shape)
# print(annotation.columns)

del annotation['refseq_sseqid']
del annotation['refseq_qseq']
del annotation['refseq_qlen']
del annotation['refseq_qframe']
del annotation['refseq_qstart']
del annotation['refseq_qend']
del annotation['refseq_nident']
del annotation['refseq_pident']
del annotation['refseq_evalue']
del annotation['uniprot_sseqid']
del annotation['uniprot_qseq']
del annotation['uniprot_qlen']
del annotation['uniprot_qframe']
del annotation['uniprot_qstart']
del annotation['uniprot_qend']
del annotation['uniprot_nident']
del annotation['uniprot_pident']
del annotation['uniprot_evalue']
del annotation['trembl_sseqid']
del annotation['trembl_qseq']
del annotation['trembl_qlen']
del annotation['trembl_qframe']
del annotation['trembl_qstart']
del annotation['trembl_qend']
del annotation['trembl_nident']
del annotation['trembl_pident']
del annotation['trembl_evalue']

annotation['count_annotation'] = annotation.apply(lambda x: 0 if x.isna().sum() == 13 else 1, axis=1)
annotation = annotation[annotation.count_annotation != 0]
del annotation['count_annotation']

annotation.to_csv("Tsor2_annotated_contigs.tsv", sep='\t', encoding='utf-8', index=False)

#TODO Choose the best annotation for the contigs.

# write fasta files
def formatSequence(sequence, length):
    formatedSequence = ''
    listTemp = wrap(sequence, length)
    for temp in listTemp:
        formatedSequence = formatedSequence + str(temp) + '\n'
    return formatedSequence.rstrip('\n')

f= open("Tsor2_annotated_contigs.fna","w+")
for i, row in annotation.iterrows():
    for seq in SeqIO.parse("Tsor2_Trinity_Unigene.fasta", "fasta"):
        if str(row['contig_id']) == str(seq.id): #not annotated by refseq
            f.write(">" + str(seq.id)+"\n")
            f.write(formatSequence(str(seq.seq), 80) + "\n")
        else:
            pass
f.close()

print("Done.\n")