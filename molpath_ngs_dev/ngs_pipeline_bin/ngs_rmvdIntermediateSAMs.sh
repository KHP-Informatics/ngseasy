#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1
#$ -V
##################################
## clean up intermediate file ####
##################################

sample_name=$1
sample_dir=$2
sample_temp=$3

cd ${sample_dir}


mkdir ${sample_dir}/sge_out

echo "move sge out to '${sample_dir}'/sge_out"

mv -v ${sample_dir}/rmvIntermediateSAMs.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/MarkDuplicates.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/SortSam.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/sam2bam.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/novoalign.${sample_name}.* ${sample_dir}/sge_out/

echo " clean up files"

rm -v ${sample_dir}/${sample_name}.alnSrtRG.*
rm -v ${sample_dir}/${sample_name}.alnSrt.*
rm -v ${sample_dir}/${sample_name}.aln.sam
rm -v ${sample_dir}/${sample_name}.aln.sam.bai

echo "clean temp dir"

rm -f -v ${sample_temp}/*.*
