# Trimmomatic
- http://www.usadellab.org/cms/?page=trimmomatic  

## Trimmomatic: A flexible read trimming tool for Illumina NGS data

### Citations
```
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.
Description
```

Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection of trimming steps and their associated parameters are supplied on the command line.
The current trimming steps are:

- ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.  
- SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.  
- LEADING: Cut bases off the start of a read, if below a threshold quality  
- TRAILING: Cut bases off the end of a read, if below a threshold quality  
- CROP: Cut the read to a specified length  
- HEADCROP: Cut the specified number of bases from the start of the read  
- MINLEN: Drop the read if it is below a specified length  
- TOPHRED33: Convert quality scores to Phred-33  
- TOPHRED64: Convert quality scores to Phred-64  

It works with FASTQ (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used), either uncompressed or gzipp'ed FASTQ. Use of gzip format is determined based on the .gz extension.
For single-ended data, one input and one output file are specified, plus the processing steps. For paired-end data, two input files are specified, and 4 output files, 2 for the 'paired' output where both reads survived the processing, and 2 for corresponding 'unpaired' output where a read survived, but the partner read did not.