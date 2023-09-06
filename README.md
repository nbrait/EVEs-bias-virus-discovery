# A Tale of Caution: Endogenous viral elements might bias virus discovery in transcriptomic datasets

**Nadja Brait, Thomas Hackl, Côme Morel, Antoni Exbrayat, Serafín Gutierrez, Sebastian Lequime**

Groningen Institue for Evolutionary Life Sciences, University of Groningen, Netherlands

### Purpose 
Repository containing source code and data supporting our manuscript on EVEs from Orthomyxoviruses in mosquitoes

### Abstract
Large-scale metagenomic and -transcriptomic studies have revolutionized our understanding of viral diversity and abundance. In contrast, endogenous viral elements (EVEs), remnants of viral sequences integrated into host genomes, have received limited attention in the context of virus discovery, especially in RNA-seq data. As EVEs can retain similarity to the exogenous virus from which they originated, this could lead to misidentification and challenges in distinguishing between active viral infections and integrated viral remnants. This poses a notable challenge, affecting the accurate classification of viruses and introducing biases during the bioinformatic process. By assessing the influence of EVEs at multiple stages of a virus discovery pipeline, this study aimed to explore their impact on data integrity and classification accuracy. 
We conducted a case study utilizing 13 public genomic and 538 transcriptomic datasets of Culicinae mosquitoes to examine EVEs and exogenous viral sequences linked to Orthomyxoviridae, a diverse family of negative-sense segmented RNA viruses. Our analysis revealed a substantial number of viral sequences in transcriptomic datasets. However, upon comparing these hits to beforehand detected genomic EVE sequences, a significant portion of them appeared to be transcribed EVEs rather than exogenous viruses. Classifying contigs into EVE transcripts or exogenous associated sequences turned out especially difficult in samples with low viral abundance. Furthermore, mapping reads on a host genome containing EVEs before viral sequence assembly led to a drastic reduction of initial detected viral hits. This led to reduced read depth and coverage, particularly in regions displaying a minimum similarity of 74% to EVEs. However, upon removal of EVEs from the reference genome, we discovered three transcribed EVEs with full-length segments, devoid of frameshift and nonsense mutations, and exhibiting high enough mean read depths, that would qualify them as exogenous virus hits. 
Our study highlights that our knowledge of the genetic diversity of viruses can be altered by the underestimated presence of EVEs in transcriptomic datasets, leading to false positives, altered or missing sequence information. Recognizing and addressing the influence of EVEs in virus discovery pipelines enhances our ability to capture the full spectrum of viral diversity.


### Citation

### Project Workflow
![Alt text](/final_workflow.jpg?raw=true "Project workflow")

### Code
All code used to generate the figures and analyses in this manuscript are available in this repository, in the `scripts` repo.

### Work in Progress
