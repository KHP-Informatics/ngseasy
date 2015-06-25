# Stampy #

Investigators: 	Gerton Lunter and Martin Goodson  
 	 
Description: 	Stampy is a package for the mapping of short reads from illumina sequencing machines onto a reference genome. It's recommended for most workflows, including those for genomic resequencing, RNA-Seq and Chip-seq. Stampy excels in the mapping of reads containing that contain sequence variation relative to the reference, in particular for those containing insertions or deletions. It can map reads from a highly divergent species to a reference genome for instance. Stampy achieves high sensitivity and speed by using a fast hashing algorithm and a detailed statistical model. Stampy has the following features:
Maps single, paired-end and mate pair Illumina reads to a reference genome  

Fast: about 20 Gbase per hour in hybrid mode (using BWA)  
Low memory footprint: 2.7 Gb shared memory for a 3Gbase genome  
High sensitivity for indels and divergent reads, up to 10-15%  
Low mapping bias for reads with SNPs  
Well calibrated mapping quality scores  
Input: Fastq and Fasta; gzipped or plain  
Output: SAM, Maq's map file  
Optionally calculates per-base alignment posteriors  
Optionally processes part of the input  
Handles reads of up to 4500 bases  

Reference: 	Lunter and Goodson. Stampy: a statistical algorithm for sensitive and fast mapping of Illumina sequence reads. Genome Res. 2011. 21:936-939. [Link to paper](http://genome.cshlp.org/content/21/6/936)
Resources: 	Stampy documentation http://www.well.ox.ac.uk/~gerton/README.txt  
Stampy registration and download http://www.well.ox.ac.uk/software-download-registration  
Contact: 	Gerton Lunter gerton.lunter@well.ox.ac.uk  