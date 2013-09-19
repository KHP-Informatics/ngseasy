#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1
#$ -V

####################
## MarkDuplicates ##
####################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

cd ${sample_dir}


${java_v1_7}/java  -XX:ParallelGCThreads=1 -Xmx${java_mem}g -jar ${ngs_picard}/MarkDuplicates.jar \
TMP_DIR=${sample_temp} \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
CREATE_INDEX=true \
REMOVE_DUPLICATES=false \
ASSUME_SORTED=true \
INPUT=${sample_dir}/${sample_name}.alnSrtRG.bam \
METRICS_FILE=${sample_dir}/${sample_name}.novoraw.MarkDuplicates \
OUTPUT=${sample_dir}/${sample_name}.novoraw.bam;


# ngs_MarkDuplicates.sh
