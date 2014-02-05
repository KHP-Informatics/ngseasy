#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

################################################
## PrintReads ##################################
################################################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

cd ${sample_dir}

echo "PrintReads"

${java_v1_7}/java  -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T PrintReads -R ${reference_genome_seq} \
-I ${sample_dir}/${sample_name}.novorealn.bam \
-o ${sample_dir}/${sample_name}.novorecal.bam \
-baq RECALCULATE \
-baqGOP 40 \
-BQSR ${sample_dir}/${sample_name}.novorealn.recal_data.table;

###########
## index ##
###########

samtools index ${sample_dir}/${sample_name}.novorecal.bam;

######################################
## make bed file [needs header etc] #
#####################################
echo "bam to bed file" 

samtools view -b -h -q 20 -F 1796 ${sample_dir}/${sample_name}.novorecal.bam | bamToBed -i stdin > ${sample_dir}/${sample_name}.novorecal.bed;

##########################
## FindCoveredIntervals ##
##########################
echo "FindCoveredIntervals"

${java_v1_7}/java  -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T FindCoveredIntervals -R ${reference_genome_seq} \
-I ${sample_dir}/${sample_name}.novorecal.bam \
-o ${sample_dir}/${sample_name}.novorecal.CoveredIntervals.list \
--coverage_threshold 10;


##echo "................................................................................................"
##echo "clean up intermediate SAM/BAMs"

##rm -v ${sample_dir}/${sample_name}.aln.sam;
##rm -v ${sample_dir}/${sample_name}.aln.bam;
##rm -v ${sample_dir}/${sample_name}.alnSrt.ba*;
##rm -v ${sample_dir}/${sample_name}.alnSrtRG.ba*;
##rm -v ${sample_dir}/${sample_name}.novoraw.ba*;
##rm -v ${sample_dir}/${sample_name}.novorealn.ba*;

## 
##echo "................................................................................................"
##echo "clean up sge output "

##########mv -v novoalign.${sample_name}.* ${sample_dir}/sge_out/
##########mv -v sam2bam.${sample_name}.* ${sample_dir}/sge_out/
##########mv -v SortSam.${sample_name}.* ${sample_dir}/sge_out/
##########mv -v AddOrReplaceReadGroups.${sample_name}.* ${sample_dir}/sge_out/
##########mv -v MarkDuplicates.${sample_name}.* ${sample_dir}/sge_out/

##########mv -v RealignerTargetCreator.${sample_name}.* ${sample_dir}/sge_out/ 
##########mv -v IndelRealigner.${sample_name}.* ${sample_dir}/sge_out/
##########mv -v BaseRecalibrator_before.${sample_name}.* ${sample_dir}/sge_out/
##########mv -v PrintReads_BQSR.${sample_name}.* ${sample_dir}/sge_out/




















