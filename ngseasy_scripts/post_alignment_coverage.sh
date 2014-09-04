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

echo "................................................"
echo " NGSeasy: START Post Alignmnet Coverage Calculations " `date`
echo "................................................"
echo ""

if [ ! -s ${SOUT}/reports/${BAM_PREFIX}.genomecov.bed ]
then
echo " NGSeasy: Running bedtools genomecov [-bga] " `date`

  /usr/local/pipeline/bedtools2/bin/bedtools genomecov -ibam ${SOUT}/alignments/${BAM_PREFIX}.bam -bga > ${SOUT}/reports/${BAM_PREFIX}.genomecov.bed;

else
  echo " NGSeasy: bedtools genomecov " `date`
fi

# Picardt tools CollectTargetedPcrMetric
echo " NGSeasy: Run CollectTargetedPcrMetric " `date`
echo " NGSeasy: This requires a Custom or Exome BED File provided by the user or manufacturer of your NGS Exome Library " 

if [ "${NGS_TYPE}" == "TGS" ] && [ -s ${BED_ANNO} ] && [ ! -s ${SOUT}/reports/${BAM_PREFIX}.target_coverage ]
then
  echo " NGSeasy: Finding TGS coverage over target regions in ${BED_ANNO} " `date`
  
    java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectTargetedPcrMetrics.jar \
    TMP_DIR=${SOUT}/tmp \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=100000 \
    INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
    OUTPUT=${SOUT}/reports/${BAM_PREFIX}.target_coverage \
    REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
    AMPLICON_INTERVALS=${BED_ANNO} \
    TARGET_INTERVALS=${BED_ANNO} \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    PER_TARGET_COVERAGE=${SOUT}/reports/${BAM_PREFIX}.per_target_coverage;
    
elif [ "${NGS_TYPE}" == "WEX" ]  && [ -s ${BED_ANNO} ] && [ ! -s ${SOUT}/reports/${BAM_PREFIX}.target_coverage ]
then
  echo " NGSeasy: Finding WEX Coverage over annotated exons/genes from ensembl " `date`
    java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectTargetedPcrMetrics.jar \
    TMP_DIR=${SOUT}/tmp \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=100000 \
    INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
    OUTPUT=${SOUT}/reports/${BAM_PREFIX}.target_coverage \
    REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
    AMPLICON_INTERVALS=${REFGenomes}/${BED_ANNO} \
    TARGET_INTERVALS=${REFGenomes}/${BED_ANNO} \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    PER_TARGET_COVERAGE=${SOUT}/reports/${BAM_PREFIX}.per_target_coverage;
    
elif [ "${NGS_TYPE}" == "WGS" ] && [ ! -s ${SOUT}/reports/${BAM_PREFIX}.target_coverage ]
then
  echo " NGSeasy: Finding WGS Coverage over 500bp windows " `date`
    java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectTargetedPcrMetrics.jar \
    TMP_DIR=${SOUT}/tmp \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=100000 \
    INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
    OUTPUT=${SOUT}/reports/${BAM_PREFIX}.target_coverage \
    REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
    AMPLICON_INTERVALS=${REFGenomes}/human_g1k_v37_0.5Kwindows.bed \
    TARGET_INTERVALS=${REFGenomes}/human_g1k_v37_0.5Kwindows.bed \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    PER_TARGET_COVERAGE=${SOUT}/reports/${BAM_PREFIX}.per_target_coverage;

else

  echo " NGSeasy: Whoops! Something Went wrong! NGS_TYPE and Annotation Files not found! Check your config file and data...or not..."
  echo " NGSeasy: Skipping Collect Targeted Pcr Metrics...."
  echo " NGSeasy: You may want to run this manually later....it is quite nice"
  echo ""
fi

echo "................................................"
echo " NGSeasy: END Post Alignmnet Coverage Calculations " `date`
echo "................................................"
echo ""
