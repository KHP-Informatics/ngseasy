# ExomeDepth

********

- http://cran.r-project.org/web/packages/ExomeDepth/index.html  
- http://plagnol-lab.blogspot.co.uk/2013/11/faq-and-clarifications-for-exomedepth-r.html  

## A robust model for read count data in exome sequencing experiments and implications for copy number variant calling

- http://bioinformatics.oxfordjournals.org/content/28/21/2747.long  
- http://cran.r-project.org/web/packages/ExomeDepth/vignettes/ExomeDepth-vignette.pdf

**Motivation:** Exome sequencing has proven to be an effective tool to discover the genetic basis of Mendelian disorders. It is well established that copy number variants (CNVs) contribute to the etiology of these disorders. However, calling CNVs from exome sequence data is challenging. A typical read depth strategy consists of using another sample (or a combination of samples) as a reference to control for the variability at the capture and sequencing steps. However, technical variability between samples complicates the analysis and can create spurious CNV calls.

**Results:** Here, we introduce ExomeDepth, a new CNV calling algorithm designed to control for this technical variability. ExomeDepth uses a robust model for the read count data and uses this model to build an optimized reference set in order to maximize the power to detect CNVs. As a result, ExomeDepth is effective across a wider range of exome datasets than the previously existing tools, even for small (e.g. one to two exons) and heterozygous deletions. We used this new approach to analyse exome data from 24 patients with primary immunodeficiencies. Depending on data quality and the exact target region, we find between 170 and 250 exonic CNV calls per sample. Our analysis identified two novel causative deletions in the genes GATA2 and DOCK8.

**Availability:** The code used in this analysis has been implemented into an R package called ExomeDepth and is available at the Comprehensive R Archive Network (CRAN).

Contact: v.plagnol@ucl.ac.uk