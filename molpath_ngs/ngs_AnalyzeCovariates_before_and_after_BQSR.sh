#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

###############################################
## AnalyzeCovariates with default covars
################################################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

cd ${sample_dir}

${java_v1_7}/java -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ${reference_genome_seq}  \
-before ${sample_dir}/${sample_name}.novorealn.recal_data.table \
-after  ${sample_dir}/${sample_name}.novorecal.recal_data.table \
-csv    ${sample_dir}/${sample_name}.novoRealnRecal.AnalyzeCovariates.BQSR.csv \
-plots  ${sample_dir}/${sample_name}.novoRealnRecal.AnalyzeCovariates.BQSR.pdf;
#-log    ${sample_dir}/${sample_name}.novoRealnRecal.AnalyzeCovariates.log;


##########mv -v AnalyzeCovariates_before_and_after_BQSR.${sample_name}.* ${sample_dir}/sge_out/
