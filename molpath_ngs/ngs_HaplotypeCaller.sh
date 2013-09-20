#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

###############################################
## HaplotypeCaller ############################
###############################################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

## HaplotypeCaller option 
stand_call_conf=${4} ## 40.0 ; default 30
stand_emit_conf=${5} ## 10.0 ; default 30

cd ${sample_dir}

## run HaplotypeCaller

${java_v1_7}/java  -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${reference_genome_seq} \
-I ${sample_dir}/${sample_name}.novorecal.bam \
-o ${sample_dir}/${sample_name}.novorecal.HaplotypeCaller.raw.snps.indels.vcf \
--dbsnp ${b37_dbsnp} \
-stand_call_conf ${stand_call_conf} \
-stand_emit_conf ${stand_emit_conf} \
--bamOutput ${sample_dir}/${sample_name}.novorecal.HaplotypeCaller.bam \
--bamWriterType CALLED_HAPLOTYPES;

##########mv -v HaplotypeCaller.${sample_name}.* ${sample_dir}/sge_out/
