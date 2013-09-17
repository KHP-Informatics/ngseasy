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

echo "move sge output to '${sample_dir}'/sge_out"


mv -v *.o[0-9][0-9]* ${sample_dir}/sge_out/
mv -v *.e[0-9][0-9]* ${sample_dir}/sge_out/
mv -v *.po[0-9][0-9]* ${sample_dir}/sge_out/
mv -v *.pe[0-9][0-9]* ${sample_dir}/sge_out/

echo "clean up temp dir "

rm -r -f -v ${sample_temp}/*.*



                       
