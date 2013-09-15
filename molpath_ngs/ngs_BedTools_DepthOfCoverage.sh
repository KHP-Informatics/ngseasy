#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

####################
## Call BedTools_DepthOfCoverag ##
####################
sample_name=${1}
sample_dir=${2}
sample_temp=${3}
bed_list=${ngs_pipeline}/${4}  ## cancer; exome, genome, etc 

for i in ${bed_list} ${bed_list}_exons ${bed_list}_exons_50;do

coverageBed -a ${sample_dir}/${sample_name}.novorecal.bed -b ${i}.bed > ${sample_dir}/${sample_name}.novorecal.coverage_${i};

coverageBed -hist -a ${sample_dir}/${sample_name}.novorecal.bed -b ${i}.bed  > ${sample_dir}/${sample_name}.novorecal.coverage_hist.${i};

coverageBed -hist -d -a ${sample_dir}/${sample_name}.novorecal.bed -b ${i}.bed > ${sample_dir}/${sample_name}.novorecal.coverage_hist.d.${i};

coverageBed -d -a ${sample_dir}/${sample_name}.novorecal.bed -b ${i}.bed  > ${sample_dir}/${sample_name}.novorecal.coverage_cancer_d.${i};


done
