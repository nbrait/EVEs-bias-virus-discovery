# A Tale of Caution: Endogenous viral elements might bias virus discovery in transcriptomic datasets

Repository containing source code and data supporting our manuscript on EVEs from Orthomyxoviruses in mosquitoes

Nadja Brait, Sebastian Lequime

### Abstract
With the advent of large-scale metagenomic studies, viral genetic diversity and abundance have been characterised at an unparalleled speed and depth. However, the viral genomes detected are often incomplete:  while this can be partly attributed to suboptimal viral sequence enrichment, especially in publicly available datasets that were not generated for viral discovery, we argue that some of these sequences may not be exogenous viruses but likely represent endogenous viral elements (EVEs). Non-retroviral EVEs are a highly underestimated bias in metagenomics/transcriptomics, and seldom removed in viral discovery bioinformatics pipelines. However, EVEs can be transcribed, which could potentially be falsely assigned to an exogenous virus in traditional virus discovery pipelines. Here, we used publicly available genomic and transcriptomic datasets of Culicinae mosquitoes to characterise EVEs and exogenous viral sequences associated with Orthomyxoviridae, a family of negative-sense segmented RNA viruses. 
After screening 13 publicly available mosquito genomes, we describe the currently most extended EVE collection in Orthomyxoviridae, with a total of 224 detected EVE sequences. All EVEs corresponded solely to Aedes datasets and 5 out of 8 known Orthomyxovirus segments, with the majority (81%) being nucleoprotein sequences, were detected. Further on, we screened over 500 publicly available transcriptomic datasets and detected 95 with positive virus hits, including 23 samples with complete viral genomes. Still, we also observed a vast number of datasets with only some of the segments available, alone or in addition to the fully segmented sequences. A quantity of these detected sequences contain frameshift mutations and display low amino acid similarities to known exogenous viral reference sequences. Yet, some of these transcribed EVEs have full segment lengths without frameshift and nonsense mutations as well as high enough mean read depths to be considered as positive virus hits for classical virus discovery pipelines. 
Our study highlights that our knowledge of the genetic diversity of viruses gained through large-scale metavirome studies can be altered by the underestimated presence of EVEs in transcriptomic datasets. This has two major implications: (i) potential false positives results in computational virus discovery pipelines based on full or partial sequences, which might derive either from exogenous viruses or EVEs, and (ii) in case of co-occurence, chimeric assemblies of EVEs and exogenous viruses, leading to false or biased discoveries. In both cases, these could alter downstream analyses, such as molecular clock dating or the characterization of evolutionary pressures.

### Citation

### Code
All code used to generate the figures and analyses in this manuscript are available in this repository, in the `scripts` repo.

### Work in Progress
