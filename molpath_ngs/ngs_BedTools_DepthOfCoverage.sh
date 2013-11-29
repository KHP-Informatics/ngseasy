#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

####################
## Call BedTools_DepthOfCoverage ##
####################
sample_name=${1}
sample_dir=${2}
sample_temp=${3}
bed_list=${ngs_pipeline}/${4}  ## cancer_genes; exome_genes, genome_genes, etc 
bed_type=${5}  ## region eg target bed, exons,  or whole_genome

cd ${sample_dir}

## human_g1k_v37.bed
## head ${bed_list}.bed
## head ${bed_list}_exons.bed
## head ${bed_list}_exons_50.bed ## padded with 50 base paires


if [ $bed_type == "region" ]; then
for i in \
${bed_list};do
##${bed_list}_genes_exons \
##${bed_list}_genes_exons_50;do
coverageBed -a ${sample_dir}/${sample_name}.novorecal.bed -b ${i}.bed > ${sample_dir}/${sample_name}.novorecal.coverage;
coverageBed -hist -a ${sample_dir}/${sample_name}.novorecal.bed -b ${i}.bed  > ${sample_dir}/${sample_name}.novorecal.coverage_hist;
coverageBed -hist -d -a ${sample_dir}/${sample_name}.novorecal.bed -b ${i}.bed > ${sample_dir}/${sample_name}.novorecal.coverage_hist.d;
done

fi

## default is to to genome 
## human_g1k_v37.bed

if [ $bed_type == "whole_genome" ]; then

genome_bed="/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/human_g1k_v37"

coverageBed -a ${sample_dir}/${sample_name}.novorecal.bed -b ${genome_bed}.bed > ${sample_dir}/${sample_name}.novorecal.coverage_${genome_bed};

coverageBed -hist -a ${sample_dir}/${sample_name}.novorecal.bed -b ${genome_bed}.bed  > ${sample_dir}/${sample_name}.novorecal.coverage_hist.${genome_bed};

## coverageBed -hist -d -a ${sample_dir}/${sample_name}.novorecal.bed -b ${genome_bed}.bed > ${sample_dir}/${sample_name}.novorecal.coverage_hist.d.${genome_bed};

fi

##########mv -v DepthOfCoverage.${sample_name}.* ${sample_dir}/sge_out/ 

