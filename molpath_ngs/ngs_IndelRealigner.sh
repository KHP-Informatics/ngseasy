#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

#############################
## IndelRealigner ###
#############################

echo "get env var IndelRealigner"
echo ".................................................."

sample_name=${1};
sample_dir=${2};
sample_temp=${3};

cd ${sample_dir}


## run IndelRealigner ##

${java_v1_7}/java  -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference_genome_seq}  \
-I ${sample_dir}/${sample_name}.novoraw.bam -o ${sample_dir}/${sample_name}.novorealn.bam -targetIntervals ${sample_dir}/${sample_name}.novoraw.output.intervals \
-known ${b37_1000G_indels} \
-known ${b37_Mills_Devine_2hit_indels} \
-model USE_READS -LOD 3 -compress 0;


## index
samtools index ${sample_dir}/${sample_name}.novorealn.bam 



















