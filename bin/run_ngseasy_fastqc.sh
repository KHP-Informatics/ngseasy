#!/bin/bash

# run_ngseasy_fastqc_pre_trimm.sh 
pre_or_post=${1}
config_tsv=${2}

# usage function 
usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker fastqc container:
  See NGSEasy containerized instructions.

  EXAMPLE USAGE:

  bash run_ngseasy_fastqc_pre_trimm.sh <1> <config.file.tsv>

EOF 
}

#-------------------------------------------------------------------#

#check exists.
 if [[ ! -e ${config_tsv} ]] 
 then
	  echo " ${config_tsv} does not exist "
	  usage;
	  exit 1;
 fi

#-------------------------------------------------------------------#

# read config file 

# begin reading config file line by line
  
while read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15
do

# set varibales  
  PROJECT_ID=$f1
  SAMPLE_ID=$f2
  FASTQ1=$f3
  FASTQ2=$f4
  PROJECT_DIR=$f5 
  DNA_PREP_LIBRARY_ID=$f6
  NGS_PLATFORM=$f7
  NGS_TYPE=$f8
  BED_ANNO=$f9
  PIPELINE=$f10
  ALIGNER=$f11
  VARCALLER=$f12
  GTMODEGATK=$f13
  CLEANUP=$f14
  NCPU=$f15

#-------------------------------------------------------------------#
#BAM PREFIX 
  BAM_PREFIX=${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}
  echo " NGSeasy: Setting BAM_PREFIX directory [$BAM_PREFIX]"
 
#OUTPUT SAMPLE DIR /home/pipeman/
 SOUT=${PROJECT_DIR}/${PROJECT_ID}/${SAMPLE_ID}
#-------------------------------------------------------------------#

#check dir exists
if [ ! -d ${PROJECT_DIR}/${PROJECT_ID}/${SAMPLE_ID} ]
then
  echo " NGSeasy: Cant Find Project directory. This is then end. Please Stop and check everything is ok " `date`
  exit 1

else 
  echo " NGSeasy: Setting OUTPUT directory [${SOUT}]"
fi

#check for fastq files exist
if [ ! -s ${SOUT}/fastq/${FASTQ1} ] && [ ! -s ${SOUT}/fastq/${FASTQ2} ]
then
  echo " NGSeasy: Can't Find fastq files [${SOUT}/fastq/${FASTQ1}] and [${SOUT}/fastq/${FASTQ2}] in ${SOUT}/fastq/ " `date`
  exit 1
fi

#set new names for copied fastq files
  rawFASTQ1=`basename ${SOUT}/fastq/${FASTQ1} _1.fq.gz`
  rawFASTQ2=`basename ${SOUT}/fastq/${FASTQ2} _2.fq.gz`
    
echo " NGSeasy: Fastq Basename : [$rawFASTQ1] "

#-------------------------------------------------------------------#
#FASTQC on raw files
echo " NGSeasy: START FASTQC " `date`

#check if qc'd data alread exists 
if [ -s ${SOUT}/fastq/${rawFASTQ1}_1.fq_fastqc.zip ] && [ -s ${SOUT}/fastq/${rawFASTQ2}_2.fq_fastqc.zip ]
then
  echo " FastQC Data already exists...skipping"
else
fi


# if before trimmomatic
if [ "${pre_or_post}" == "1" ]
then
  echo " NGSeasy: Run Pre-Alignment QC on raw Fastq files " `date`

  /usr/local/pipeline/FastQC/fastqc \
    --threads ${NCPU} \
    --extract \
    --quiet \
    --dir ${SOUT}/tmp \
    --outdir ${SOUT}/fastq \
    ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz \
    ${SOUT}/fastq/${rawFASTQ2}_2.fq.gz
fi

#if after trimmomatic
if [ "${pre_or_post}" == "0" ]
then
# Trimmomatic paired output
  qcdPeFASTQ1=${SOUT}/fastq/${rawFASTQ1}_1.filtered.fq.gz;
  qcdPeFASTQ2=${SOUT}/fastq/${rawFASTQ2}_2.filtered.fq.gz;
# Trimmomatic unpaired ouput
  qcdSeFASTQ1=${SOUT}/fastq/${rawFASTQ1}_1.unpaired.fq.gz;
  qcdSeFASTQ2=${SOUT}/fastq/${rawFASTQ2}_2.unpaired.fq.gz;
   
#check if fastq files exist
if [ ! -s ${SOUT}/fastq/${rawFASTQ1}_1.filtered.fq.gz ] && [ ! -s ${SOUT}/fastq/${rawFASTQ2}_2.filtered.fq.gz ]
then
  echo " NGSeasy: Can'd Find Trimmed Fastq files in ${SOUT}/fastq/ ie [${rawFASTQ1}_1.filtered.fq.gz] or [${rawFASTQ1}_1.unpaired.fq.gz] " `date`
  exit 1
fi

echo " NGSeasy: Run FastQC on Trimmed Fastq files " `date`

  /usr/local/pipeline/FastQC/fastqc \
    --threads ${NCPU} \
    --extract \
    --quiet \
    --dir ${SOUT}/tmp \
    --outdir ${SOUT}/fastq \
    ${SOUT}/fastq/${rawFASTQ1}_1.filtered.fq.gz \
    ${SOUT}/fastq/${rawFASTQ2}_2.filtered.fq.gz \
    ${SOUT}/fastq/${rawFASTQ1}_1.unpaired.fq.gz \
    ${SOUT}/fastq/${rawFASTQ2}_2.unpaired.fq.gz
fi

echo " NGSeasy: END FASTQC  " `date`

done < ${config_tsv}

