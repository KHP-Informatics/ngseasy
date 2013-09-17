#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=3G
#$ -q short.q
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 8


####################
## Call Novoalign ##
####################

N_CPU=8

echo " number of cpus " $N_CPU

echo "starting novoalign " 

${ngs_novo}/novoalign \
-d ${reference_genome_novoindex} \
-f ${fastq_dir}/${fastq_prefix}_1.fastq  ${fastq_dir}/${fastq_prefix}_2.fastq  \
-F STDFQ \
--Q2Off \
--3Prime  \
-g 65 \
-x 7 \
-r All \
-i PE 300,150 \
-c ${N_CPU} \
-k -K ${sample_dir}/${sample_name}.novoalign.K.stats \
-o SAM > ${sample_dir}/${sample_name}.aln.sam;
 
echo "end novoalign"

###########################################################################################################


#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -q short.q
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1

################################
## covert sam to bam and sort ##
################################
echo "convert sam to bam and sort"

samtools view -bS ${sample_dir}/${sample_name}.aln.sam > ${sample_dir}/${sample_name}.aln.bam;
samtools index    ${sample_dir}/${sample_name}.aln.bam;

##########################################################################################################

#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -q short.q
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1
#############
## SortSam ##
#############

${java_1.7}/java  -XX:ParallelGCThreads=1 -Xmx6g -jar ${ngs_picard}/SortSam.jar \
TMP_DIR=${sample_temp} \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
CREATE_INDEX=true \
SORT_ORDER=coordinate \
INPUT=${sample_dir}/${sample_name}.aln.bam \
OUTPUT=${sample_dir}/${sample_name}.alnSrt.bam;


##############################################################################################################

#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -q short.q
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1
##############################
## AddOrReplaceReadGroups   ##
##############################

${java_1.7}/java  -XX:ParallelGCThreads=1 -Xmx6g -jar ${ngs_picard}/AddOrReplaceReadGroups.jar \
TMP_DIR=${sample_temp} \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
CREATE_INDEX=true \
SORT_ORDER=coordinate \
RGID=${mRGID} \
RGLB=${mRGLB} \
RGPL=${mRGPL} \
RGPU=${mRGPU} \
RGSM=${mRGSM} \
RGCN=${mRGCN} \
RGDS=${mRGDS} \
RGDT=${mRGDT} \
INPUT=${sample_dir}/${sample_name}.alnSrt.bam \
OUTPUT=${sample_dir}/${sample_name}.alnSrtRG.bam;

###################################################################################################

#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -q short.q
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1
####################
## MarkDuplicates ##
####################

${java_1.7}/java  -XX:ParallelGCThreads=1 -Xmx6g -jar ${ngs_picard}/MarkDuplicates.jar \
TMP_DIR=${sample_temp} \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
CREATE_INDEX=true \
REMOVE_DUPLICATES=false \
ASSUME_SORTED=true \
INPUT=${sample_dir}/${sample_name}.alnSrtRG.bam \
METRICS_FILE=${sample_dir}/${sample_name}.novo \
OUTPUT=${sample_dir}/${sample_name}.novo.bam;


##################################################################


#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -l h_vmem=8G
#$ -q short.q
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -pe multi_thread 1
##################################
## clean up intermediate file ####
##################################

mkdir ${sample_dir}/sge_out

mv -v rmvIntermediateSAMs.${sample_name} ${sample_dir}/sge_out/
mv -v MarkDuplicates.${sample_name} ${sample_dir}/sge_out/
mv -v SortSam.${sample_name} ${sample_dir}/sge_out/
mv -v sam2bam.${sample_name} ${sample_dir}/sge_out/
mv -v novoalign.${sample_name} ${sample_dir}/sge_out/

rm -v ${sample_dir}/${sample_name}.alnSrtRG.*
rm -v ${sample_dir}/${sample_name}.alnSrt.*
rm -v ${sample_dir}/${sample_name}.aln.sam
rm -v ${sample_dir}/${sample_name}.aln.sam.bai






















