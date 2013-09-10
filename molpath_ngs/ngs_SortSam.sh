#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1
#$ -V
#############
## SortSam ##
#############

echo " Sort BAM file and index"

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

cd ${sample_dir}


java -XX:ParallelGCThreads=1 -Xmx6g -jar ${ngs_picard}/SortSam.jar \
TMP_DIR=${sample_temp} \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
CREATE_INDEX=true \
SORT_ORDER=coordinate \
INPUT=${sample_dir}/${sample_name}.aln.bam \
OUTPUT=${sample_dir}/${sample_name}.alnSrt.bam;


# ngs_SortSam.sh
