#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1
#$ -V

################################
## covert sam to bam and sort ##
################################
echo "convert novoaligne sam to bam and index"

sample_name=${1}
sample_dir=${2}

cd ${sample_dir}

samtools view -bS ${sample_dir}/${sample_name}.aln.sam > ${sample_dir}/${sample_name}.aln.bam;
samtools index    ${sample_dir}/${sample_name}.aln.bam;


