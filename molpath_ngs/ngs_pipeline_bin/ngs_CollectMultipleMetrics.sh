#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

####################
## Call CollectMultipleMetrics ##
####################
sample_name=${1}
sample_dir=${2}
sample_temp=${3}


java -XX:ParallelGCThreads=1 -Xmx${java_mem}g -jar ${ngs_picard}/CollectMultipleMetrics.jar \
TMP_DIR=${sample_temp} \
INPUT=${sample_dir}/${sample_name}.novorecal.bam \
ASSUME_SORTED=true \
REFERENCE_SEQUENCE=${reference_genome_seq} \
OUTPUT=${sample_dir}/${sample_name}.novorecal.qc \
PROGRAM=CollectAlignmentSummaryMetrics \
PROGRAM=QualityScoreDistribution \
PROGRAM=CollectInsertSizeMetrics \
PROGRAM=MeanQualityByCycle;

java -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T FlagStat -R ${reference_genome_seq} \
-I ${sample_dir}/${sample_name}.novorecal.bam \
-o ${{sample_dir}/${sample_name}.novorecal.FlagStat;

java -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T FindCoveredIntervals -R ${reference_genome_seq} \
-I ${sample_dir}/${sample_name}.novorecal.bam \
-o ${{sample_dir}/${sample_name}.novorecal.overedIntervals.list \
--coverage_threshold 10 \



