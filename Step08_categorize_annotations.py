'''
What does this program do?
	It categorizes the annotations acording to a mannually built repository of keywords, which are representatives of these categories

Input:
    Tsor2_annotated_contigs.tsv
    keywords_secreted_protein.txt
    keywords_housekeeping.txt
    keywords_immunity.txt
    keywords_enzyme.txt
    keywords_bacteria.txt
    keywords_virus.txt
    keywords_other.txt

Output:
	temp01.tsv

@author: Waldeyr Mendes Cordeiro da Silva, June 2020
'''

import re
import pandas as pd


def getListfromTSV(filename: str):
    result = []
    with open(filename, 'r') as f:
        line = f.readline()
        while line:
            result.append(str(line).strip())
            line = f.readline()
    return result

keywords_secreted_protein = getListfromTSV('keywords/keywords_secreted_protein.txt')
keywords_housekeeping = getListfromTSV('keywords/keywords_housekeeping.txt')
keywords_immunity = getListfromTSV('keywords/keywords_immunity.txt')
keywords_enzyme = getListfromTSV('keywords/keywords_enzyme.txt')
keywords_bacteria = getListfromTSV('keywords/keywords_bacteria.txt')
keywords_virus = getListfromTSV('keywords/keywords_virus.txt')
keywords_other = getListfromTSV('keywords/keywords_other.txt')

df = pd.read_csv('Tsor2_annotated_contigs.tsv', sep="\t")

# all frames need to have the same size, they were completed with @@@@@
keywords = pd.DataFrame({
    "secreted_protein": keywords_secreted_protein,
    "housekeeping": keywords_housekeeping,
    "enzyme": keywords_enzyme,
    "immunity": keywords_immunity,
    "bacteria": keywords_bacteria,
    "virus": keywords_virus,
    "other": keywords_other
})
result = pd.DataFrame(
    columns=['contig_id', 'secreted_protein', 'housekeeping', 'enzyme', 'immunity', 'bacteria', 'virus', 'other'])
# for each contig
for i, row in df.iterrows():
    rowResultDict = {}
    # for each keyword in each category
    secreted_protein = housekeeping = enzyme = immunity = bacteria = virus = other = ""
    annotations_to_search = \
        str(row['macrogen_go']) + " " + \
        str(row['macrogen_uniprot']) + " " + \
        str(row['macrogen_nr']) + " " + \
        str(row['macrogen_pfam']) + " " + \
        str(row['macrogen_eggnog']) + " " + \
        str(row['macrogen_nt']) + " " + \
        str(row['macrogen_ko_euk']) + " " + \
        str(row['refseq_annotation']) + " " + \
        str(row['uniprot_annotation']) + " " + \
        str(row['trembl_annotation'])
    tax_to_search = str(row['uniprot_taxonomy']) + " " + str(row['refseq_taxonomy'])
    for index, value in keywords.iterrows():
        if (re.search(str(value['secreted_protein']), annotations_to_search, re.IGNORECASE)):
            secreted_protein = str(value["secreted_protein"])
        if (re.search(str(value['housekeeping']), annotations_to_search, re.IGNORECASE)):
            housekeeping = str(value["housekeeping"])
        if (re.search(str(value['enzyme']), annotations_to_search, re.IGNORECASE)):
            enzyme = str(value["enzyme"])
        if (re.search(str(value['immunity']), annotations_to_search, re.IGNORECASE)):
            immunity = str(value["immunity"])
        if (re.search(str(value['bacteria']), tax_to_search, re.IGNORECASE)):
            bacteria = str(value["bacteria"])
        if (re.search(str(value['virus']), tax_to_search, re.IGNORECASE)):
            virus = str(value["virus"])
        if (re.search(str(value['other']), annotations_to_search, re.IGNORECASE)):
            other = str(value["other"])
        rowResultDict = {
            'contig_id': str(row['contig_id']),
            'secreted_protein': secreted_protein,
            'housekeeping': housekeeping,
            'enzyme': enzyme,
            'immunity': immunity,
            'bacteria': bacteria,
            'virus': virus,
            'other': other
        }
    result = result.append(rowResultDict, ignore_index=True)

print(df.shape)
Tsor2_annotations_plus_keywords = merged_inner = pd.merge(left=df, right=result,
                                                          left_on='contig_id',
                                                          right_on='contig_id',
                                                          how='outer')
print(Tsor2_annotations_plus_keywords.shape)

Tsor2_annotations_plus_keywords.to_csv("temp01.tsv", sep='\t', encoding='utf-8',
                                       index=False)
