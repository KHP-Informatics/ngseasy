#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

################################
## covert sam to bam and sort ##
################################
echo "convert novoaligne sam to bam and index"

sample_name=${1}
sample_dir=${2}

cd ${sample_dir}

ngs="/scratch/project/pipelines/ngs_pipeline_dev/ngs_dev_sjn_tmp/ngs_bin"

${ngs}/samtools view -bS ${sample_dir}/${sample_name}.aln.sam > ${sample_dir}/${sample_name}.aln.bam;
${ngs}/samtools index    ${sample_dir}/${sample_name}.aln.bam;


