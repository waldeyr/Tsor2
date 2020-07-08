from Bio import SeqIO
from textwrap import wrap
import pandas as pd

def formatSequence(sequence, length):
    formatedSequence = ''
    listTemp = wrap(sequence, length)
    for temp in listTemp:
        formatedSequence = formatedSequence + str(temp) + '\n'
    return formatedSequence.rstrip('\n')

def getListfromTSV(filename: str):
    result = []
    with open(filename, 'r') as f:
        line = f.readline()
        while line:
            result.append(str(line).strip())
            line = f.readline()
    return result

notAnnotatedContigs = getListfromTSV('Tsor2_not_annotated_contigs_list.txt')
listOfNotAnnotated = []
for value in notAnnotatedContigs:
    for seq in SeqIO.parse("Tsor2_Trinity_Unigene.fasta", "fasta"):
        if str(value) == str(seq.id): #not annotated by refseq
            print(">" + str(seq.id))
            print(formatSequence(str(seq.seq), 80)) 
        else:
            pass
                 

# python3 select_contigs_for_TrEMBL_running.py > Tsor2_Trinity_Unigene_to_submit_to_TrEMBL.fasta