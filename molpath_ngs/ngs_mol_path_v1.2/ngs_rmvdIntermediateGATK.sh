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

mv -v ${sample_dir}/novoalign.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/sam2bam.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/SortSam.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/AddOrReplaceReadGroups.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/MarkDuplicates.${sample_name}.* ${sample_dir}/sge_out/

mv -v ${sample_dir}/RealignerTargetCreator.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/IndelRealigner.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/BaseRecalibrator_before.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/PrintReads_BQSR.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/BaseRecalibrator_after.${sample_name}.* ${sample_dir}/sge_out/
mv -v ${sample_dir}/AnalyzeCovariates_before_and_after_BQSR.${sample_name}.* ${sample_dir}/sge_out/

echo " clean up files"

#rm -v ${sample_dir}/${sample_name}.aln.sam
rm -v ${sample_dir}/${sample_name}.alnSrt.ba*
rm -v ${sample_dir}/${sample_name}.alnSrtRG.ba*
#rm -v ${sample_dir}/${sample_name}.novoraw.ba*
rm -v ${sample_dir}/${sample_name}.novorealn.ba*


echo "clean temp dir"
rm -r -f -v ${sample_temp}/*.*

mv -v ${sample_dir}/rmvIntermediateGATK.${sample_name}.* ${sample_dir}/sge_out/


                       
