#!/bin/bash

# Stephen Newhouse <stephen.j.newhouse@gmail.com>
# UNDER HEAVY DEVELOPMENT
echo ""
echo ""
echo  "#------------------------------------------------------------#"
echo  " Starting NGSeasy : NGS made easy!!!!"
echo  " Version 1.0"
echo  " Authors: Amos Folarin <amosfolarin@gmail.com>"
echo  " Authors: Stephen Newhouse <stephen.j.newhouse@gmail.com>"
echo  " Run Date : `date +"%d-%m-%y"`"
echo  "#------------------------------------------------------------#"
echo ""

# refgenomes and contaminant lists
  REFGenomes='/home/pipeman/reference_genomes_b37'
  adapter_fa='/home/pipeman/reference_genomes_b37/contaminant_list.fa'
  gatk_resources='/home/pipeman/gatk_resources'
  FASTQDIR='/home/pipeman/fastq_raw'
  KNOWN_SNPS_b138=/usr/local/pipeline/gatk_resources/dbsnp_138.b37.vcf
  echo " NGSeasy: Note:- Using ${KNOWN_SNPS_b138} for SNP annotations in GATK.\n Edit bash script if you want to change this " `date`
  
#---------------------------------------------CONFIG------------------------------------------------#
# get the following from the config file

# set varibales
  DATE=`date +"%d%m%y"`
  POJECT_ID=${1}
  SAMPLE_ID=${2}
  FASTQ1=${3}
  FASTQ2=${4}
  PROJECT_DIR=${5}
  DNA_PREP_LIBRARY_ID=${6}
  NGS_PLATFORM=${7}
  NGS_TYPE=${8}
  BED_ANNO=${9}
  PIPELINE=${10}
  ALIGNER=${11}
  VARCALLER=${12}
  GTMODEGATK=${13}
  CLEANUP=${14}
  NCPU=${15}

  
echo ""
echo "................................................"
echo " NGSeasy: START NGSeasy Pipleine [${PIPELINE}]" `date`
echo "................................................"
echo ""

##############################
## BEGING FULL NGS PIPELINE ##
##############################

## BAM PREFIX 
  BAM_PREFIX=${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}
  echo " NGSeasy: Setting BAM_PREFIX directory [$BAM_PREFIX]"
  echo ""
 
## OUTPUT SAMPLE DIR
 SOUT=${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
  
if [ ! -e ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID} ]
then
  echo " NGSeasy: Cant Find Project directory. This is then end. Please Stop and check everything is ok " `date`
  exit 1

else 
  echo " NGSeasy: Setting OUTPUT directory [${SOUT}]"
fi
  
##---------------------- POST ALIGNMENT QC ----------------------##

echo "................................................"
echo " NGSeasy: START Post Alignmnet QC " `date`
echo "................................................"
echo ""

if [ ! -s ${SOUT}/reports/${BAM_PREFIX}.alignment_summary_metrics ]
then 

  echo " NGSeasy: Running Post Alignmnet QC " `date`
# CollectMultipleMetrics
  echo " NGSeasy: CollectMultipleMetrics " `date`
 
 java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectMultipleMetrics.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
  OUTPUT=${SOUT}/reports/${BAM_PREFIX} \
  REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
  PROGRAM=CollectAlignmentSummaryMetrics \
  PROGRAM=CollectInsertSizeMetrics \
  PROGRAM=QualityScoreDistribution \
  PROGRAM=MeanQualityByCycle;

# CollectAlignmentSummaryMetrics
  echo " NGSeasy: CollectAlignmentSummaryMetrics " `date`
  
  java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectAlignmentSummaryMetrics.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
  OUTPUT=${SOUT}/reports/${BAM_PREFIX}.alignment_summary_metrics_alt \
  REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
  ASSUME_SORTED=true \
  METRIC_ACCUMULATION_LEVEL=SAMPLE;

# CollectWgsMetrics
  echo " NGSeasy: CollectWgsMetrics " `date`
  
  java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectWgsMetrics.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
  OUTPUT=${SOUT}/reports/${BAM_PREFIX}.wgs_coverage \
  REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
  MINIMUM_MAPPING_QUALITY=20 \
  MINIMUM_BASE_QUALITY=20 \
  COVERAGE_CAP=1000;

  awk 'NR>9' ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage > ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage_hist
  sed '8q' ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage > ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage_stats

# FlagStat
  echo " NGSeasy: FlagStat " `date`
  
    java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FlagStat -R ${REFGenomes}/human_g1k_v37.fasta \
    -I ${SOUT}/alignments/${BAM_PREFIX}.bam  \
    -o ${SOUT}/reports/${BAM_PREFIX}.FlagStat;
    
else
    echo " NGSeasy: Post Alignmnet QC Already Run " `date` 
fi

echo ""
echo "................................................"
echo " NGSeasy: END Post Alignmnet QC " `date`
echo "................................................"
echo ""
