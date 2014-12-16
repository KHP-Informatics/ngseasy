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
  
##---------------------- FASTQ-QC ----------------------## 

  if [ ! -s ${SOUT}/fastq/${FASTQ1} ] && [ ! -s ${SOUT}/fastq/${FASTQ2} ]
then
  echo " NGSeasy: Copying fastq files from ${FASTQDIR}/ to ${SOUT}/fastq/ " `date`
    cp ${FASTQDIR}/${FASTQ1} ${SOUT}/fastq/${FASTQ1}
    cp ${FASTQDIR}/${FASTQ2} ${SOUT}/fastq/${FASTQ2}

else
  echo " NGSeasy: Fastq Files exist in  ${SOUT}/fastq/ " `date`
  ls ${SOUT}/fastq/
fi

## set new names for copied fastq files
  rawFASTQ1=`basename ${SOUT}/fastq/${FASTQ1} _1.fq.gz`
  rawFASTQ2=`basename ${SOUT}/fastq/${FASTQ2} _2.fq.gz`
    
echo " NGSeasy: Fastq Basename : [$rawFASTQ1] "

# FASTQC on raw files

echo ""
echo "................................................"
echo " NGSeasy: START Pre-Alignment QC " `date`
echo "................................................"
echo ""

# check if qc'd data alread exists 
if [ ! -s ${SOUT}/fastq/${rawFASTQ1}_1.fq_fastqc.zip ] && [ ! -s ${SOUT}/fastq/${rawFASTQ2}_2.fq_fastqc.zip ]
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
    
else
  echo " NGSeasy: Pre-Alignment QC on raw Fastq files Already run " `date`
fi

echo ""
echo "................................................"
echo " NGSeasy: START Filtering/Trimming " `date`
echo "................................................"
echo ""

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

echo ""
echo "................................................"
echo " NGSeasy: END Filtering/Trimming " `date`
echo "................................................"
echo ""