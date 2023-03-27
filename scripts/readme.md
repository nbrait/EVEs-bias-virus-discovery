# Scripts

In order to run the project, please install the packages that are required and check their version, it may occur that a different version of a package does not work as intended.

## Databases

Note that the databases are not present within this repository, these can be downloaded individually.

NCBI nr database: downloaded 9.11.2021 [current release](https://ftp.ncbi.nlm.nih.gov/blast/db/)

most tools were installed with `sudo apt-get install` in bash environment or from CRAN-like repositories through R-Studio - otherwise installation guides are given below. 
## Prerequisites
* SRA-toolkit (v2.10.9): (NCBI 2012), [(installation guide)](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) - download SRA datasets
* efetch (v14.6), fastq-dump (v2.10.9): NCBI 2012 - download and processing of NCBI data
* Trimmomatic (v0.39): Bolger, Lohse, and Usadel 2014 [(installation guide)](http://www.usadellab.org/cms/?page=trimmomatic) - quality and length trimming
* Bowtie2 (v2.3.5.1): Langmead and Salzberg 2012; - read mapping
* Samtools (v1.10), htslib (v1.10.2-3): (H. Li et al. 2009) - processing, indexing alignment files
* bedtools (v2.27.1): Quinlan and Hall 2010 - coverage analysis 
* Ray-tools (v2.3.1): [(installation guide)](https://github.com/sebhtml/Ray-Releases) - read assembly into contigs
* Diamond blast (v2.0.13.151):Buchfink, Xie, and Huson 2015 [(installation guide)](https://github.com/bbuchfink/diamond) - blastx search
* NCBI BLAST+ (v.2.13.0): Altschul et al. 1997; Johnson et al. 2008 - tblastn, tblastx and blastx search for EVE screening
* Lofreq* (2.1.6): Wilm et al 2012 [(installation guide)](https://github.com/CSB5/lofreq/tree/master/dist) - variant caller
* MAFFT (v7.490): Katoh and Standley 2013 - multiple sequence alignment
* IQTREE, Modelfinder, UFBoot2 (v1.6.12): (Nguyen et al. 2015), (Kalyaanamoorthy et al. 2017) - tree construction
* R (v4.1.2): R. C. Team 2021
* R-Studio (v2022.07.2): Rs. Team 2021
* ggtree (v3.2.1): Yu et al. 2018 - tree visualization
* ggplot2 (v3.3.5): Wickham 201 - data visualization

## Manual steps 

After the pre-processing script, contigs were mapped to a non-redundant referencing sequence list corresponding to their best blastx hits in Geneious Prime (v.2022.1.1) [https://www.geneious.com/](https://www.geneious.com/ "see here"). Contigs spanning only parts of the Orthomyxovirus segments were artificially extended to the full segment length by merging them with their reference sequence determined by their best blast hit. 

Options: `Mapper: Geneious` `Sensitivity: Highest Sensitivity/Medium` `Interation: Up to 5 times` `Map multiple best matches: Randomly`

Consensus sequences of the extended virus sequences from Geneious are used as input for the post-processing script. 
Fasta headers have to look as followed: `[SRA-accesion]_[Segment].fasta` and were concatenated into a single file called `references.fasta`

