# A tale of caution: How endogenous viral elements affect virus discovery in transcriptomic data

**Nadja Brait, Thomas Hackl, Côme Morel, Antoni Exbrayat, Serafín Gutierrez, Sebastian Lequime**

### Purpose 
Repository containing source code and data supporting our manuscript on EVEs from Orthomyxoviruses in mosquitoes

### Abstract
Large-scale metagenomic and -transcriptomic studies have revolutionized our understanding of viral diversity and abundance. In contrast, endogenous viral elements (EVEs), remnants of viral sequences integrated into host genomes, have received limited attention in the context of virus discovery, especially in RNA-Seq data. EVEs resemble their original viruses, a challenge that makes distinguishing between active infections and integrated remnants difficult, affecting virus classification and biases downstream analyses. Here, we systematically assessed the effects of EVEs on a prototypical virus discovery pipeline, evaluated their impact on data integrity and classification accuracy, and provide some recommendations for better practices.
We examined EVEs and exogenous viral sequences linked to Orthomyxoviridae, a diverse family of negative-sense segmented RNA viruses, in 13 genomic and 538 transcriptomic datasets of Culicinae mosquitoes. Our analysis revealed a substantial number of viral sequences in transcriptomic datasets. However, a significant portion appeared not to be exogenous viruses but transcripts derived from EVEs. Distinguishing between transcribed EVEs or exogenous virus sequences was especially difficult in samples with low viral abundance. For example, three transcribed EVEs showed full-length segments, devoid of frameshift and nonsense mutations, and exhibiting sufficient mean read depths that would qualify them as exogenous virus hits. Mapping reads on a host genome containing EVEs before assembly somewhat alleviated the EVE burden but it led to drastic reduction of viral hits and reduced quality of assemblies especially in regions of the viral genome relatively similar to EVEs.
Our study highlights that our knowledge of the genetic diversity of viruses can be altered by the underestimated presence of EVEs in transcriptomic datasets, leading to false positives and altered or missing sequence information. Thus, recognizing and addressing the influence of EVEs in virus discovery pipelines will be key to enhancing our ability to capture the full spectrum of viral diversity.

### Citation
not yet peer-reviewed

### Project Workflow
![Alt text](/Workflow_figure_github.png?raw=true "Project workflow")

### Code
All code used to generate the figures and analyses in this manuscript are available in this repository, in the `scripts` repo.
