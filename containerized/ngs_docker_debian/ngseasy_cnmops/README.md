# cn-MOPS

**** 

- http://www.bioinf.jku.at/software/cnmops/  


## cn.MOPS: Mixture Of PoissonS for discovering Copy Number variations in next generation sequencing data

cn.MOPS is an algorithm that accurately detects copy number variations in next generation sequencing data in a study of multiple samples. 
Quantitative analyses of next generation sequencing (NGS) data, such as the detection of copy number variations (CNVs), remain challenging. 
Current methods detect CNVs as changes in the depth of coverage along chromosomes. 
Technological or genomic variations in the depth of coverage thus lead to a high false discovery rate (FDR), even upon correction for GC content. 
In the context of association studies between CNVs and disease, a high FDR means many false CNVs, thereby decreasing the discovery power of the 
study after correction for multiple testing. 
We propose “Copy Number estimation by a Mixture Of PoissonS” (cn.MOPS), a data processing pipeline for CNV detection in NGS data. 
In contrast to previous approaches, cn.MOPS incorporates modeling of depths of coverage across samples at each genomic position. 
Therefore, cn.MOPS is not affected by read count variations along chromosomes. 
Using a Bayesian approach, cn.MOPS decomposes variations in the depth of coverage across samples into integer copy numbers and noise by means of its mixture components and Poisson distributions, respectively. 
The noise estimate allows for reducing the FDR by filtering out detections having high noise which are likely to be false detections. 
We compared cn.MOPS with the five most popular methods for CNV detection methods in 
NGS data using four benchmark data sets: (1) simulated data, (2) NGS data from a male HapMap individual with implanted CNVs from the X chromosome, (3) data from HapMap individuals with known CNVs, 
(4) high coverage data from the 1000 Genomes Project. cn.MOPS outperformed its five competitors in terms of precision (1–FDR) and recall for both gains and losses in all benchmark data sets.  

Please cite:  

Günter Klambauer, Karin Schwarzbauer, Andreas Mayr, Djork-Arné Clevert, Andreas Mitterecker, Ulrich Bodenhofer, Sepp Hochreiter. "cn.MOPS: mixture of Poissons for discovering copy number variations in next generation sequencing data with a low false discovery rate." Nucleic Acids Research 2012 40(9); doi:10.1093/nar/gks003.

