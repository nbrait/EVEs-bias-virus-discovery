#!/bin/bash

# making blast database
/home/nbrait/TOOLS/diamond makedb --in references_EXO_WD.fasta -d EXO_REFS

# translating corrected sequences 
home/nbrait/TOOLS/diamond blastx -d ./EXO_REFS -q nucleotide_corrected.fasta --sensitive --quiet -F 15 -f 6 qseqid sseqid pident length full_qseq qseq qseq_translated slen btop -k 1 --un unmatched.fa --masking 0 -o aa_sequences.tsv

# mafft alignments were created in geneious per segment

# creating trees
iqtree -nt AUTO -s NP_alignment.fasta -m TEST -bb 1000 -bnni -pre NP_tree
iqtree -nt AUTO -s PB1_alignment.fasta -m TEST -bb 1000 -bnni -pre PB1_tree
iqtree -nt AUTO -s PB2_alignment.fasta -m TEST -bb 1000 -bnni -pre PB2_tree
iqtree -nt AUTO -s PA_alignment.fasta -m TEST -bb 1000 -bnni -pre PA_tree
iqtree -nt AUTO -s GP_alignment.fasta -m TEST -bb 1000 -bnni -pre GP_tree
iqtree -nt AUTO -s HP1_alignment.fasta -m TEST -bb 1000 -bnni -pre HP1_tree
iqtree -nt AUTO -s HP2_alignment.fasta -m TEST -bb 1000 -bnni -pre HP2_tree
iqtree -nt AUTO -s HP3_alignment.fasta -m TEST -bb 1000 -bnni -pre HP3_tree
