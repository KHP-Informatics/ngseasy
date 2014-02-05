#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999999999999999999
#$ -V

###############################################
## VCFtoolsSiteFilter ###
###############################################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

cd ${sample_dir}

INPUT=${sample_dir}/${sample_name}.novorecal.UnifiedGenotyper.raw.snps.indels \
OUTPUT=${sample_dir}/${sample_name}.novorecal.UnifiedGenotyper.raw.snps.indels \

#######################################
## Run VCFSiteFilter ##
########################################

vcftools --vcf ${INPUT}.vcf \
--recode-INFO-all \
--min-meanDP 60 \
--max-meanDP 99999 \
--minQ 30 \
--out ${OUTPUT} \
--recode;






