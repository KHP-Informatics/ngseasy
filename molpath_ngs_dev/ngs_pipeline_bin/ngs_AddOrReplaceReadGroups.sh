#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

##############################
## AddOrReplaceReadGroups   ##
##############################

echo "Start AddOrReplaceReadGroups"

sample_name=${1}
sample_dir=${2}
sample_temp=${3}
mRGID=${4}
mRGLB=${5}
mRGPL=${6}
mRGPU=${7}
mRGSM=${8}
mRGCN=${9}
mRGDS=${10}
mRGDT=${11}

cd ${sample_dir}

${java_1.7}/java  -XX:ParallelGCThreads=1 -Xmx${java_mem}g -jar ${ngs_picard}/AddOrReplaceReadGroups.jar \
TMP_DIR=${sample_temp} \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
CREATE_INDEX=true \
SORT_ORDER=coordinate \
RGID=${mRGID} \
RGLB=${mRGLB} \
RGPL=${mRGPL} \
RGPU=${mRGPU} \
RGSM=${mRGSM} \
RGCN=${mRGCN} \
RGDS=${mRGDS} \
RGDT=${mRGDT} \
INPUT=${sample_dir}/${sample_name}.alnSrt.bam \
OUTPUT=${sample_dir}/${sample_name}.alnSrtRG.bam;
# ngs_AddOrReplaceReadGroups.sh
