#!/bin/bash

# only applicable for Aedes samples 
# creating exogenous and endogenous Orthomyxovirus databases
/home/nbrait/TOOLS/diamond makedb --in genomic_EVEs_aa.fasta -d EVE_REFS
/home/nbrait/TOOLS/diamond makedb --in orthomyxo_references.fasta -d EXO_REFS

/home/nbrait/TOOLS/diamond blastx -d ./EVEs_REF -q .viral_contigs_RAW.fasta --sensitive --quiet -F 15-f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore btop\
       -k 1 --un unalignedEVE_matches.fa -o ./EVE_matches.tsv
/home/nbrait/TOOLS/diamond blastx -d ./EXO_REFS -q ./viral_contigs_RAW --sensitive --quiet -F 15 -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore btop\
       -k 1 --un unalignedEXO_matches.fa -o ./EXO_matches.tsv
