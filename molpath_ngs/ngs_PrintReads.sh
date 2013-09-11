#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

###############################################
## PrintReads 
################################################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

cd ${sample_dir}

java -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T PrintReads -R ${reference_genome_seq} -S SILENT -l ERROR \
-I ${sample_dir}/${sample_name}.novorealn.bam \
-o ${sample_dir}/${sample_name}.novorecal.bam \
-BQSR ${sample_dir}/${sample_name}.novorealn.recal_data.table \
-log ${sample_dir}/${sample_name}.novorecal.PrintReads.log;

## index
samtools index ${sample_dir}/${sample_name}.novorecal.bam

