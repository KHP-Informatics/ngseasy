#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

###############################################
## UnifiedGenotyper ############################
###############################################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

## UnifiedGenotyper option 
stand_call_conf=${4} ## 40.0 ; default 30
stand_emit_conf=${5} ## 10.0 ; default 30

cd ${sample_dir}

${java_v1_7}/java  -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${reference_genome_seq} \
-I ${sample_dir}/${sample_name}.novorecal.bam \
-o ${sample_dir}/${sample_name}.novorecal.UnifiedGenotyper.raw.snps.indels.vcf \
--dbsnp ${b37_dbsnp} \
-stand_call_conf ${stand_call_conf} \
-stand_emit_conf ${stand_emit_conf} \
--genotype_likelihoods_model BOTH \
-baq CALCULATE_AS_NECESSARY \
-baqGOP 30;

# check
if [ -e ${sample_dir}/${sample_name}.novorecal.UnifiedGenotyper.raw.snps.indels.vcf ]
then
echo " Checking VCF"
head ${sample_dir}/${sample_name}.novorecal.UnifiedGenotyper.raw.snps.indels.vcf;
wc -l ${sample_dir}/${sample_name}.novorecal.UnifiedGenotyper.raw.snps.indels.vcf;
fi





##########mv -v UnifiedGenotyper.${sample_name}.* ${sample_dir}/sge_out/
