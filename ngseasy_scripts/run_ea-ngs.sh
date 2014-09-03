#!/bin/bash

#------------------------------------------------------------#
# NGSeasy
# Version 1.0
#------------------------------------------------------------#

#------------------------------------------------------------#
# AUTHORS
#------------------------------------------------------------#
#
# Stephen Newhouse
# stephen.j.newhouse@gmail.com
# Amos Folarin
# amosfolarin@gmail.com
#
#------------------------------------------------------------#



#------------------------------------------------------------#
# USEAGE
#------------------------------------------------------------#
#
# sh run_ea-ngs.sh <config.file>
#
#------------------------------------------------------------#


echo  "#------------------------------------------------------------#"
echo  " Starting NGSeasy : NGS made easy!!!!"
echo  " Version 1.0"
echo  " Authors: Amos Folarin <amosfolarin@gmail.com>"
echo  " Authors: Stephen Newhouse <stephen.j.newhouse@gmail.com>"
echo  " USAGE:  sh run_ea-ngs.sh <config.file>"
echo  " Run Date : `date +"%d-%m-%y"`"
echo  "#------------------------------------------------------------#"
echo ""
  
#------------------------------------------------------------#
# Config file
#------------------------------------------------------------#

  config_file=${1}
 
echo  " NGS config_file : [${config_file}]"

#------------------------------------------------------------#
# Read ngs config file 
#------------------------------------------------------------#
  
echo  " Checking [${config_file}] format"

  numfeilds=`awk '{print NF}' ${config_file} | sort -gr | sed '1q'`
  numsamples=`wc -l ${config_file} | awk '{print $1}'`
echo  " Number of samples : [$numsamples]"
    
  if [ $numfeilds -ne 15 ]
    then
    echo  " WARNING! Number of fields not equal to 15 for one of your samples You appear to have missing data in [${config_file}]"
    echo  " Check  [${config_file}] format"
    echo  " The [config.file] should contain the following 15 columns for each sample:-"
    echo  " 
      POJECT_ID
      SAMPLE_ID
      FASTQ1
      FASTQ2
      PROJECT_DIR
      DNA_PREP_LIBRARY_ID
      NGS_PLATFORM
      NGS_TYPE
      BED_ANNO
      PIPELINE
      ALIGNER
      VARCALLER
      GTMODEGATK
      CLEANUP
      NCPU"
    echo  " Data are [TAB] delimited and each line should end in a [NEW LINE]"
      exit 1
    else
   echo  " All Good...Proceding to read : [${config_file}]"
  fi

# begin reading config file line by line
  
while read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15
do

# set varibales  
  POJECT_ID=$f1
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
  
# refgenomes and contaminant lists
  REFGenomes='/home/pipeman/reference_genomes_b37'
  adapter_fa='/home/pipeman/reference_genomes_b37/contaminant_list.fa'
  gatk_resources='/home/pipeman/gatk_resources'
  FASTQDIR='/home/pipeman/fastq_raw'

  
#------------------------------------------------------------#
# Checking Project Directory 
#------------------------------------------------------------#
  
echo  "#------------------------------------------------------------#"
echo  " Checking Project/Sample Directory Structure"
echo  "#------------------------------------------------------------#"
echo ""
  
  if [ ! -d $PROJECT_DIR/${POJECT_ID} ]
    then
    echo  " Project Directory Does not exist. Setting up correct Project Directory"
    echo ""
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/]"
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/config_files]"
    echo  " Making Directory : [$PROJECT_DIR/${POJECT_ID}/cohort_vcfs]"
        mkdir $PROJECT_DIR/${POJECT_ID}/
        mkdir ${PROJECT_DIR}/${POJECT_ID}/config_files
        mkdir ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs
  else
    echo  " Project Directory Exists"
  fi

  if [ ! -d $PROJECT_DIR/${POJECT_ID}/$SAMPLE_ID ]  
    then
    echo  " Sample Directory Does not exist. Setting up correct Sample Directory"
    echo ""
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/fastq
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/tmp
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/alignments
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/vcf
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/reports
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/fastq]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/tmp]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/alignments]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/vcf]"
    echo  " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/reports]"
      else
    echo  " Sample Directory Exists"
  fi
  
#------------------------------------------------------------#
# Pipeline Settings print to screen
#------------------------------------------------------------#
echo ""
echo  "#------------------------------------------------------------#"
echo  " Pipeline Settings "
echo  "#------------------------------------------------------------#"
echo ""
  DATE=`date +"%d%m%y"`
echo  " Run Date : [$DATE]"
echo  " Project_Id : [$POJECT_ID]" 
echo  " Sample_Id  : [$SAMPLE_ID]"
echo  " Fastq File 1 : [$FASTQ1]"
echo  " Fastq File 2 : [$FASTQ2]"
echo  " Project Directory : [$PROJECT_DIR]"
echo  " DNA Library : [$DNA_PREP_LIBRARY_ID]"
echo  " Platform : [$NGS_PLATFORM]"
echo  " NGS Type : [$NGS_TYPE]"
echo  " BED file annotation : [$BED_ANNO]"
echo  " NGS Pipeline : [$PIPELINE]"
echo  " NGS Aligner  : [$ALIGNER]"
echo  " Variant Caller : [$VARCALLER]"
echo  " Genotyping Mode (GATK Specific) : [$GTMODEGATK]"
echo  " Clean up TRUE/FALE : [$CLEANUP]"
echo  " Number of cpu : [$NCPU]"
echo  " Output Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}]"
echo  " Sample Prefix : [${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}]"
echo  " Fastq Directory : [$FASTQDIR]"
echo ""
echo  " NOTE: Paths are relative to the Structure set up in the Docker Image."
echo ""
echo " These will be [/home/pipeman] and [/usr/local/pipeline] on the image and mounted to your local directories"
echo ""
echo " See Doc's for details."
echo ""
 
#------------------------------------------------------------#
# run pipeline
#------------------------------------------------------------#

echo  "#------------------------------------------------------------#"
echo  " Run NGS Pipeline "
echo  "#------------------------------------------------------------#"
echo ""
  
echo  " Saving Sample configuration settings to : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE}.config.file]"
echo  "$DATE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 ${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE} ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}" \
  >  ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE}.config.file
echo ""
  
echo  " Saving Project configuration settings to : [${PROJECT_DIR}/${POJECT_ID}/config_files/master.config.file]"
echo  "$DATE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 ${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE} ${PROJECT_DIR}/${POJECT_ID}" >> ${PROJECT_DIR}/${POJECT_ID}/config_files/project.config.file;
echo ""
    
echo  " Processing Sample : [$SAMPLE_ID]"
echo  " Sample Data will be saved into : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/]"
echo  " Sample Prefix : [${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}]"
echo ""
echo  " Selected Pipeline : [$PIPELINE]"
echo ""
echo  " SYSTEM COMMAND: bash /usr/local/pipeline/ngseasy_scripts/${PIPELINE} $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 " 
date
echo ""
echo ""
#------------------------------------------------------------#
# run command 
#------------------------------------------------------------#  

bash /home/pipeman/ngseasy_scripts/${PIPELINE} $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15
  
done < ${config_file}


















