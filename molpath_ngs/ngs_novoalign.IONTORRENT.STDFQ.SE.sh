#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

####################
## Call Novoalign ##
####################

###################################################################################################################################
##  IONTORRENT :- http://www.novocraft.com/wiki/tiki-index.php?page=Benchmarking+Ion+Torrent+PGM+aligners&highlight=Ion%20Torrent 
## HaloPlex cancer panel


## note -g -x -H changes 
## needs exploring


N_CPU=${novo_cpu}
fastq_prefix=${1}
sample_name=${2}
sample_dir=${3}

cd ${sample_dir}


## start novo 
echo "----------------------------------------------------------------------------------------" 
echo " Fastq prefix" ${1}
echo " Sample Name" ${2}
echo " Output dir" ${3}
echo "starting novoalign " 
echo " number of cpus " $N_CPU
echo "----------------------------------------------------------------------------------------" 

cd ${sample_dir}

fastqc ${fastq_dir}/${fastq_prefix}_1.fastq 

${ngs_novo}/novoalign \
-d ${reference_genome_novoindex} \
-f ${fastq_dir}/${fastq_prefix}_1.fastq \
-F STDFQ \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-a GCTGAGGATCACCGACTGCCCATAGAGAGGCTGAGAC \
-a GTGCTCTTCCGATCT \
--Q2Off --3Prime  -g 15 -x 4 -r All -H -c ${N_CPU} \
-k -K ${sample_dir}/${sample_name}.novoalign.K.stats \
-o SAM > ${sample_dir}/${sample_name}.aln.sam;

echo "end novoalign"

echo "----------------------------------------------------------------------------------------" 

ls -l 

echo "----------------------------------------------------------------------------------------" 

samtools view -hS ${sample_dir}/${sample_name}.aln.sam | head -500
 
chmod 765 ${sample_dir}/${sample_name}.aln.sam;






















