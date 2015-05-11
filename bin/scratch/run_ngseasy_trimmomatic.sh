#!/bin/bash

# run_ngseasy_trimmomatic.sh

config_tsv=${1}

# usage function 
usage()
  {
  cat << EOF
  This script sets up the NGSeasy docker fastqc container:
  See NGSEasy containerized instructions.

  ARGUMENTS:
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  EXAMPLE USAGE:
  bash run_ngseasy_trimmomatic.sh <config.file.tsv>
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

#-------------------------------------------------------------------#
#OUTPUT SAMPLE DIR
 SOUT=${PROJECT_DIR}/${PROJECT_ID}/${SAMPLE_ID}
 
if [ ! -e ${PROJECT_DIR}/${PROJECT_ID}/${SAMPLE_ID} ]
then
  echo " NGSeasy: Cant Find Project directory. This is then end. Please Stop and check everything is ok " `date`
  exit 1

else 
  echo " NGSeasy: Setting OUTPUT directory [${SOUT}]"
fi

#-------------------------------------------------------------------#
#check fot fastq files to sample directory
if [ ! -s ${SOUT}/fastq/${FASTQ1} ] && [ ! -s ${SOUT}/fastq/${FASTQ2} ]
then
  echo " NGSeasy: Cant find fastq files from ${FASTQDIR}/ in ${SOUT}/fastq/ " `date`
  exit 1
else
  echo " NGSeasy: Fastq Files exist in  ${SOUT}/fastq/ " `date`
fi

#-------------------------------------------------------------------#
#set new names for copied fastq files
  rawFASTQ1=`basename ${SOUT}/fastq/${FASTQ1} _1.fq.gz`
  rawFASTQ2=`basename ${SOUT}/fastq/${FASTQ2} _2.fq.gz`
    
echo " NGSeasy: Fastq Basename : [$rawFASTQ1] "

#-------------------------------------------------------------------#

# Trimmomatic paired output
  qcdPeFASTQ1=${SOUT}/fastq/${rawFASTQ1}_1.filtered.fq.gz;
  qcdPeFASTQ2=${SOUT}/fastq/${rawFASTQ2}_2.filtered.fq.gz;
# Trimmomatic unpaired ouput
  qcdSeFASTQ1=${SOUT}/fastq/${rawFASTQ1}_1.unpaired.fq.gz;
  qcdSeFASTQ2=${SOUT}/fastq/${rawFASTQ2}_2.unpaired.fq.gz;
 
# if qc'd files exits then skip this step
if [ ! -s ${qcdPeFASTQ1} ] && [ ! -s ${qcdPeFASTQ2} ]
then

  # skip next step if using novoalign or ya really cant be bothered - saves about an hour
  echo " NGSeasy: Running Trimmomatic: A flexible read trimming tool for NGS data " `date`

    # run Trimmomatic
    java -jar /usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar PE \
    -threads ${NCPU} \
    -phred33 \
    -trimlog ${SOUT}/fastq/trimm_qc.log \
    ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz ${SOUT}/fastq/${rawFASTQ2}_2.fq.gz \
    ${qcdPeFASTQ1} ${qcdSeFASTQ1} \
    ${qcdPeFASTQ2} ${qcdSeFASTQ2} \
    ILLUMINACLIP:${adapter_fa}:2:30:10:5:true \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50

  echo " NGSeasy: Run Pre-Alignment QC on Filtered/Trimmed Fastq files " `date`
    # FASTQC on paired trimmed files

    /usr/local/pipeline/FastQC/fastqc \
    --threads ${NCPU} \
    --extract \
    --quiet \
    --dir ${SOUT}/tmp \
    --outdir ${SOUT}/fastq \
    ${qcdPeFASTQ1} ${qcdPeFASTQ2};    
else
  echo " NGSeasy: QC'd Fastq files Already Exist " `date`
  echo "................................................"
  zcat ${qcdPeFASTQ1} | head -4
  echo "................................................"
fi
    
    
done < ${config_file}

