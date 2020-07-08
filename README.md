# Tsor2 annotation workflow


## Dependencies

* NCBI-Blast 2.8+

* Signalp 5

* Python3.x libraries: 

** [pandas](https://pandas.pydata.org)

** [Biopython](https://github.com/biopython/DIST)

** [re](https://docs.python.org/3/library/re.html)

** [os](https://docs.python.org/3/library/os.html)



![Graphical overview of the Tsor2 workflow](https://github.com/waldeyr/Tsor2/blob/master/graphical_pipeline/Yanna.png)


## RefSeq Alignment

Database building:

`makeblastdb -in refseq_protein.fasta -out refseq_blast/RefSeq_Release_200 -dbtype prot -max_file_sz 2GB -parse_seqids -title RefSeq_Release_200 -logfile log.txt`

Blast running (Caution: it will take about 2 weeks to run):

`nohup blastx -db refseq_blast/RefSeq_Release_200 -num_threads 48 -evalue 10E-5 -query Tsor2_Trinity_Unigene.fasta -outfmt "6 qseqid qgi qacc qlen qstart qend qseq nident pident evalue positive ppos gaps qframe sseqid" -out Tsor2_RefSeq_blast.result &`

Insert headline on the result file (caution: execute it only once)

`echo -e "contig_id\tqgi\tqacc\tqlen\tqstart\tqend\tqseq\tnident\tpident\tevalue\tpositive\tppos\tgaps\tqframe\tsseqid" | cat - Tsor2_refseq_blast.result > /tmp/out && mv /tmp/out Tsor2_refseq_blast.result`

## RefSeq parse and annotation

Before run it, check if the Tsor2_refseq_blast.result has a headline with the command `head Tsor2_refseq_blast.result`

Caution: recommended 32 GB RAM

`python3 Step01_refseq_parse_and_annotation.py`


## Uniprot Alignment

Database building:

`makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot/uniprot_sprot -dbtype prot -max_file_sz 2GB -parse_seqids -title uniprot_sprot -logfile log.txt`

Blast running (Caution: it will take about 2 days to run):

`nohup blastx -db uniprot_sprot/uniprot_sprot -num_threads 48 -evalue 10E-5 -query Tsor2_Trinity_Unigene.fasta -outfmt "6 qseqid qgi qacc qlen qstart qend qseq nident pident evalue positive ppos gaps qframe sseqid" -out Tsor2_uniprot_blast.result &`

Insert headline on the result file (caution: execute it only once)

`echo -e "contig_id\tqgi\tqacc\tqlen\tqstart\tqend\tqseq\tnident\tpident\tevalue\tpositive\tppos\tgaps\tqframe\tsseqid" | cat - Tsor2_uniprot_blast.result > /tmp/out && mv /tmp/out Tsor2_uniprot_blast.result`


## TrEMBL Alignment (only for those contigs that were not annotated by any other means; see the script "Step04_build_input_for_trembl_running.py")

Database building:

`makeblastdb -in uniprot_trembl.fasta -out uniprot_trembl/uniprot_trembl -dbtype prot -max_file_sz 2GB -parse_seqids -title uniprot_trembl -logfile log.txt`

Blast running (Caution: it will take about 2 days to run):

`nohup blastx -db uniprot_trembl/uniprot_trembl -num_threads 36 -evalue 1E-5 -query Tsor2_input_for_TrEMBL.fasta -outfmt "6 qseqid qgi qacc qlen qstart qend qseq nident pident evalue positive ppos gaps qframe sseqid" -out Tsor2_trembl_blast.result &`

Insert headline on the result file (caution: execute it only once)

`echo -e "contig_id\tqgi\tqacc\tqlen\tqstart\tqend\tqseq\tnident\tpident\tevalue\tpositive\tppos\tgaps\tqframe\tsseqid" | cat - Tsor2_trembl_blast.result > /tmp/out && mv /tmp/out Tsor2_trembl_blast.result`


## SIGNALP
signalp -fasta Tsor2_Translated_Proteins.fasta -org euk -format short -prefix Tsor2_signalp -verbose -gff3