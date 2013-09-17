#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

#############################
## IndelRealigner ###
#############################

echo "get env var IndelRealigner"
echo ".................................................."

sample_name=${1};
sample_dir=${2};
sample_temp=${3};

echco ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"


echo ".................................................."
echo "fastq_dir : " ${fastq_dir}
echo "project dir : " ${aln_dir}
echo ".................................................."
echo "Sample working dir : " ${sample_dir}
echo "Sample temp dir : " ${sample_temp}
###echo "Sample information : "  ${fastq_prefix} ${sample_name} ${mRGID} ${mRGLB} ${mRGPL} ${mRGPU} ${mRGSM} ${mRGCN} ${mRGDS} ${mRGDT}
echo ".................................................."
echo "sge ques : " ${queue_name}
echo "novo_mem : " ${novo_mem}G
echo "sge_h_vmem : " ${sge_h_vmem}G
echo "gatk_java_mem : " ${gatk_java_mem}
echo "novo_cpu : " ${novo_cpu}
echo "email_contact: " ${email_contact}
echo "ngs_pipeline : " ${ngs_pipeline}
echo ".................................................."
ls -l ${ngs_pipeline}
echo ".................................................."
echo "ngs_picard : " ${ngs_picard} 
echo "ngs_novo : " ${ngs_novo} 
echo "ngs_gatk : " ${ngs_gatk} 
echo "ngs_samtools: " ${ngs_samtools} 
echo ".................................................."
ls -l ${ngs_picard}
echo ".................................................."
ls -l ${ngs_novo}
echo ".................................................."
ls -l ${ngs_gatk}
echo ".................................................."
ls -l ${ngs_samtools}
echo ".................................................."
echo "fastq seq : " ${reference_genome_seq}
head -5  ${reference_genome_seq}
echo ".................................................."
head ${b37_1000G_biallelic_indels}
echo ".................................................."
head ${b37_Mills_Devine_2hit_indels_sites}
echo ".................................................."
head ${b37_Mills_Devine_2hit_indels}
echo ".................................................."
head ${b37_1000G_omni2_5}
echo ".................................................."
head ${b37_hapmap_3_3}
echo ".................................................."
head ${b37_dbsnp_132_excluding_sites_after_129}
echo ".................................................."
head ${b37_dbsnp_132}
echo ".................................................."





echco ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

echo "move to sample dir"
echo ".................................................."
cd ${sample_dir}
ls -ls 
echo ".................................................."
echo "sample Id :" ${sample_name}
echo ".................................................."
echo "sample dir :" ${sample_dir}
echo ".................................................."
echo "sample temp :" ${sample_temp}
echo ".................................................."
echo "path to gatk :" ${ngs_gatk}
echo ".................................................."
echo "vcf1 :" ${b37_1000G_biallelic_indels}
echo ".................................................."
echo "vcf2 :" ${b37_Mills_Devine_2hit_indels_sites}
echo ".................................................."
echo "intervals file:" ${sample_dir}/${sample_name}.novoraw.output.intervals
echo ".................................................."
head ${b37_1000G_biallelic_indels}
echo ".................................................."
head ${b37_Mills_Devine_2hit_indels_sites}
echo ".................................................."
head  ${sample_dir}/${sample_name}.novoraw.output.intervals
echo ".................................................."
samtools view ${sample_dir}/${sample_name}.novoraw.bam | head -50 ;
echo ".................................................."

## run IndelRealigner ##

${java_1.7}/java  -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T IndelRealigner -R ${reference_genome_seq}  \
-I ${sample_dir}/${sample_name}.novoraw.bam -o ${sample_dir}/${sample_name}.novorealn.bam -targetIntervals ${sample_dir}/${sample_name}.novoraw.output.intervals \
-known ${b37_1000G_biallelic_indels} -known ${b37_Mills_Devine_2hit_indels_sites} \
-model USE_READS -LOD 3 -compress 0;


## index
samtools index ${sample_dir}/${sample_name}.novorealn.bam 

ls -ls
echo ".................................................."
echo "sample Id x:" ${sample_name}
echo ".................................................."
echo "sample dir x:" ${sample_dir}
echo ".................................................."
echo "sample temp x:" ${sample_temp}
echo ".................................................."
echo "path to gatk x:" ${ngs_gatk}
echo ".................................................."
echo "vcf1 x:" ${b37_1000G_biallelic_indels}
echo ".................................................."
echo "vcf2 x:" ${b37_Mills_Devine_2hit_indels_sites}
echo ".................................................."
echo "intervals file x:" ${sample_dir}/${sample_name}.novoraw.output.intervals
echo ".................................................."
head ${b37_1000G_biallelic_indels}
echo ".................................................."
head ${b37_Mills_Devine_2hit_indels_sites}
echo ".................................................."
head  ${sample_dir}/${sample_name}.novoraw.output.intervals
echo ".................................................."
samtools view ${sample_dir}/${sample_name}..novorealn.bam | head -50 ;
echo ".................................................."

echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

echo ".................................................."
echo "fastq_dir : " ${fastq_dir}
echo "project dir : " ${aln_dir}
echo ".................................................."
echo "Sample working dir : " ${sample_dir}
echo "Sample temp dir : " ${sample_temp}
#####echo "Sample information : "  ${fastq_prefix} ${sample_name} ${mRGID} ${mRGLB} ${mRGPL} ${mRGPU} ${mRGSM} ${mRGCN} ${mRGDS} ${mRGDT}
echo ".................................................."
echo "sge ques : " ${queue_name}
echo "novo_mem : " ${novo_mem}G
echo "sge_h_vmem : " ${sge_h_vmem}G
echo "gatk_java_mem : " ${gatk_java_mem}
echo "novo_cpu : " ${novo_cpu}
echo "email_contact: " ${email_contact}
echo "ngs_pipeline : " ${ngs_pipeline}
echo ".................................................."
ls -l ${ngs_pipeline}
echo ".................................................."
echo "ngs_picard : " ${ngs_picard} 
echo "ngs_novo : " ${ngs_novo} 
echo "ngs_gatk : " ${ngs_gatk} 
echo "ngs_samtools: " ${ngs_samtools} 
echo ".................................................."
ls -l ${ngs_picard}
echo ".................................................."
ls -l ${ngs_novo}
echo ".................................................."
ls -l ${ngs_gatk}
echo ".................................................."
ls -l ${ngs_samtools}
echo ".................................................."
echo "fastq seq : " ${reference_genome_seq}
head -5  ${reference_genome_seq}
echo ".................................................."
head ${b37_1000G_biallelic_indels}
echo ".................................................."
head ${b37_Mills_Devine_2hit_indels_sites}
echo ".................................................."
head ${b37_Mills_Devine_2hit_indels}
echo ".................................................."
head ${b37_1000G_omni2_5}
echo ".................................................."
head ${b37_hapmap_3_3}
echo ".................................................."
head ${b37_dbsnp_132_excluding_sites_after_129}
echo ".................................................."
head ${b37_dbsnp_132}
echo ".................................................."





















