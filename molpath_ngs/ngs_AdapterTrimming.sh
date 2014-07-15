#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

#########################################
## adapter trimming with trimmomatic   ##
#########################################

echo "Start Trimming and Quality filtering"

sample_name=${1}
sample_dir=${2}
sample_temp=${3}
# additional parameters
sample_PE=${4}
sample_adapterfile=${5}


# set PE/SE variable
sample_librarytype='SE'
if [[ sample_PE ==  ]]; then
    sample_librarytype='PE'
fi

echo "Reads are ${sample_librarytype}"



cd ${sample_dir}

${java_v1_7}/java  -XX:ParallelGCThreads=1 -Xmx${java_mem}g -jar ${ngs_trimmomatic} \
PE \
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

java -jar $TRIM/trimmomatic PE seqfile_R1_001.fastq seqfile_R2_001.fastq s1_pe s1_se s2_pe s2_se ILLUMINACLIP:/mnt/research/common-data/Bio/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36