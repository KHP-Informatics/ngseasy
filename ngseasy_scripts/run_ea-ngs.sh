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


  echo "\n\n#------------------------------------------------------------#\n"
  echo " Starting NGSeasy : NGS made easy!!!!"
  echo " Version 1.0\n"
  echo " Authors: Amos Folarin <amosfolarin@gmail.com>"
  echo " Authors: Stephen Newhouse <stephen.j.newhouse@gmail.com>\n"
  echo " USAGE:  sh run_ea-ngs.sh <config.file>\n"
  echo " Run Date : `date +"%d-%m-%y"`"
  echo "#------------------------------------------------------------#\n\n"
  
#------------------------------------------------------------#
# Config file
#------------------------------------------------------------#

  config_file=${1}
 
  echo " NGS config_file : [${config_file}]"

#------------------------------------------------------------#
# Read ngs config file 
#------------------------------------------------------------#
  
  echo " Checking [${config_file}] format"

  numfeilds=`awk '{print NF}' ${config_file} | sort -gr | sed '1q'`
  numsamples=`wc -l ${config_file} | awk '{print $1}'`
  echo " Number of samples : [$numsamples]"
    
  if [ $numfeilds -ne 15 ]
    then
      echo " WARNING!\n Number of fields not equal to 15 for one of your samples\n You appear to have missing data in [${config_file}]\n"
      echo " Check  [${config_file}] format\n"
      echo " The [config.file] should contain the following 15 columns for each sample:-\n\n"
      echo " 
      POJECT_ID\n
      SAMPLE_ID\n
      FASTQ1\n
      FASTQ2\n
      PROJECT_DIR\n
      DNA_PREP_LIBRARY_ID\n
      NGS_PLATFORM\n
      NGS_TYPE\n
      BED_ANNO\n
      PIPELINE\n
      ALIGNER\n
      VARCALLER\n
      GTMODEGATK\n
      CLEANUP\n
      NCPU\n\n"
      echo " Data are [TAB] delimited and each line should end in a [NEW LINE]\n"
      exit 1
    else
     echo " All Good\nProceding to read : [${config_file}]\n"
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

#------------------------------------------------------------#
# Checking Project Directory 
#------------------------------------------------------------#
  
  echo "#------------------------------------------------------------#"
  echo " Checking Project/Sample Directory Structure"
  echo "#------------------------------------------------------------#\n"
  
  if [ ! -d $PROJECT_DIR/${POJECT_ID} ]
    then
      echo " Project Directory Does not exist\n Setting up correct Project Directory\n"
      echo " Making Directory : [$PROJECT_DIR/${POJECT_ID}/]"
      echo " Making Directory : [$PROJECT_DIR/${POJECT_ID}/config_files]"
      echo " Making Directory : [$PROJECT_DIR/${POJECT_ID}/cohort_vcfs]\n"
        mkdir $PROJECT_DIR/${POJECT_ID}/
        mkdir ${PROJECT_DIR}/${POJECT_ID}/config_files
        mkdir ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs
  else
      echo " Project Directory Exists\n"
  fi

  if [ ! -d $PROJECT_DIR/${POJECT_ID}/$SAMPLE_ID ]  
    then
      echo " Sample Directory Does not exist\n Setting up correct Sample Directory\n"
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/fastq
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/tmp
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/alignments
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/vcf
        mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/reports
      echo " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}]"
      echo " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/fastq]"
      echo " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/tmp]"
      echo " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/alignments]"
      echo " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/vcf]"
      echo " Making Sample Directory : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/reports]"
      else
      echo " Sample Directory Exists\n"
  fi
  
#------------------------------------------------------------#
# Pipeline Settings print to screen
#------------------------------------------------------------#

  echo "\n#------------------------------------------------------------#"
  echo " Pipeline Settings "
  echo "#------------------------------------------------------------#\n"
  DATE=`date +"%d%m%y"`
  echo " Run Date : [$DATE]"
  echo " Project_Id : [$POJECT_ID]" 
  echo " Sample_Id  : [$SAMPLE_ID]"
  echo " Fastq File 1 : [$FASTQ1]"
  echo " Fastq File 2 : [$FASTQ2]"
  echo " Project Directory : [$PROJECT_DIR]"
  echo " DNA Library : [$DNA_PREP_LIBRARY_ID]"
  echo " Platform : [$NGS_PLATFORM]"
  echo " NGS Type : [$NGS_TYPE]"
  echo " BED file annotation : [$BED_ANNO]"
  echo " NGS Pipeline : [$PIPELINE]"
  echo " NGS Aligner  : [$ALIGNER]"
  echo " Variant Caller : [$VARCALLER]"
  echo " Genotyping Mode (GATK Specific) : [$GTMODEGATK]"
  echo " Clean up TRUE/FALE : [$CLEANUP]"
  echo " Number of cpu : [$NCPU]"
  echo " Sample Prefix : [${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE}]\n"
 
#------------------------------------------------------------#
# run pipeline
#------------------------------------------------------------#

  echo "\n#------------------------------------------------------------#"
  echo " Run NGS Pipeline "
  echo "#------------------------------------------------------------#\n"
  
  echo " Saving Sample configuration settings to : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE}.config.file]\n"
  echo "$DATE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 ${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE} ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}" \
  >  ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE}.config.file
  
  echo " Saving Project configuration settings to : [${PROJECT_DIR}/${POJECT_ID}/config_files/master.config.file]\n"
  echo "$DATE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 ${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE} ${PROJECT_DIR}/${POJECT_ID}" >> ${PROJECT_DIR}/${POJECT_ID}/config_files/project.config.file;
  
  echo " Processing Sample : [$SAMPLE_ID]\n"
  echo " Sample Data will be saved into : [${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/]\n"
  echo " Sample Prefix : [${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE}]\n"
  
  echo " Run Pipeline : [$PIPELINE]\n"
  echo " SYSTEM COMMAND: sh /usr/local/pipeline/ngs_scripts/${PIPELINE} $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 \n"

#------------------------------------------------------------#
# run command 
#------------------------------------------------------------#  

	sh /usr/local/pipeline/ngs_scripts/${PIPELINE} $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15
  
done < ${config_file}


















