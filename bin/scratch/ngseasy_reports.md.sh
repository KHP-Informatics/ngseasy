#!/bin/bash

################################################################
# Program: ngseasy_sample_report
# Version 1.0 
# Author: Stephen Newhouse (stephen.j.newhouse@gmail.com) http://stackoverflow.com/questions/2953081/how-can-i-write-a-here-doc-to-a-file-in-bash-script
#################################################################

echo -e "\n################################################################"
echo -e "# Program: ngseasy_sample_report"
echo -e "# Version 1.0"
echo -e "# Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)"
echo -e "#################################################################\n"


########################################################################################################
## Set version and run date
#
NGSEASYVERSION="1.0"
RUNDATE=`date +"%d%m%y"`
ngseasy_step="ngseasy_sample_report"

########################################################################################################
## global logging fuction
#
function logger_ngseasy() {
 message=${1}
 mylogfile=${2}
 echo -e [`date`]":[NGSEASY:${NGSEASYVERSION}]:"${message}":[${USER}]:[`uname -a`]" >> ${mylogfile}.log;
 echo -e [`date`]":[NGSEASY:${NGSEASYVERSION}]:"${message}":[${USER}]:[`uname -a`]"
}

########################################################################################################
## global usage
#
function ngseasy_alignment_usage() {
    echo "
Program: ngseasy_sample_report
Version 1.0
Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)

usage:   ngseasy_sample_report -c <config_file> -d <project_directory> -f <First Name> -l <Last Name>

options:  -c  STRING	configuration file
          -d  STRING	project directory
          -f  STRING	First Name
          -l  STRING	Last Name
          -h  NULL	    show this message

ngseasy_sample_report generates basic sample report (.md .html .pdf) using pandoc
"
}
########################################################################################################
## Check options passed in.
#
if test -z "$2"
then
  logger_ngseasy "[${ngseasy_step}]:ERROR:No options found"
  ngseasy_alignment_usage
  exit 1
fi

########################################################################################################
## get options for command line args
#
  while  getopts "hc:d:f:l:" opt
  do

      case ${opt} in
	  h)
	  ngseasy_alignment_usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  echo -e "CONFIG FILE [-c] = ${config_tsv}"
	  ;;

	  d)
	  project_directory=${OPTARG}
	  echo -e "PROJECT DIR [-d] = ${project_directory}"
	  ;;
	  
	  f)
	  FirstName=${OPTARG}
	  echo -e "First Name [-f] = ${FirstName}"
	  ;; 
	  
	  l)
	  LastName=${OPTARG}
	  echo -e "First Name [-f] = ${LastName}"
	  ;;
      esac
  done


########################################################################################################
## Read config file 

## check ${config_tsv}. is this a batch file or the orginal config file 
#
logger_ngseasy "[${ngseasy_step}]:Checking [${config_tsv}] format" ${HOME}/ngseasy_logs/ngseasy.${config_tsv}.${RUNDATE}

hasheader=`sed 1q ${config_tsv} | grep PROJECT_ID | wc -l | awk '{print $1}'`

if [[ "${config_tsv}" == *.batch.* ]]
then
  logger_ngseasy "[${ngseasy_step}]:[${config_tsv}] is a BACTH file ie a subset of the original config file" ${HOME}/ngseasy_logs/ngseasy.${config_tsv}.${RUNDATE}
  RUNFILE="${config_tsv}"
  logger_ngseasy "[${ngseasy_step}]:Setting RUNFILE to [${RUNFILE}]" ${HOME}/ngseasy_logs/ngseasy.${config_tsv}.${RUNDATE}
elif [[ "$hasheader" -eq 1 ]]
then
  logger_ngseasy "[${ngseasy_step}]:[${config_tsv}] header present. Removing this" ${HOME}/ngseasy_logs/ngseasy.${config_tsv}.${RUNDATE}
  logger_ngseasy "[${ngseasy_step}]:[cmd]:sed 1d \${config_tsv} > ${config_tsv}.tmp" ${HOME}/ngseasy_logs/ngseasy.${config_tsv}.${RUNDATE}
  sed 1d ${config_tsv} > ${config_tsv}.tmp
  RUNFILE="${config_tsv}.tmp"
  logger_ngseasy "[${ngseasy_step}]:Setting RUNFILE to [${RUNFILE}]" ${HOME}/ngseasy_logs/ngseasy.${config_tsv}.${RUNDATE}
else
  RUNFILE="${config_tsv}"
  logger_ngseasy "[${ngseasy_step}]:[${RUNFILE}] is seemingly perfect" ${HOME}/ngseasy_logs/ngseasy.${config_tsv}.${RUNDATE}
  logger_ngseasy "[${ngseasy_step}]:Setting RUNFILE to [${RUNFILE}]" ${HOME}/ngseasy_logs/ngseasy.${config_tsv}.${RUNDATE}
fi

########################################################################################################
## Read config and loop through all lines calling fastqc docker
#
while read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19 f20 f21 f22 f23
do
PROJECT_ID=$f1;
SAMPLE_ID=$f2;
FASTQ1=$f3;
FASTQ2=$f4;
PROJECT_DIR=$f5;
DNA_PREP_LIBRARY_ID=$f6;
NGS_PLATFORM=$f7;
NGS_TYPE=$f8;
BAIT=$f9;
CAPTURE=$f10;
GENOMEBUILD=$f11;
FASTQC=$f12;
TRIM=$f13;
BSQR=$f14;
REALN=$f15;
ALIGNER=$f16;
VARCALLER=$f17;
CNV=$f18;
ANNOTATOR=$f19;
CLEANUP=$f20;
NCPU=$f21;
NGSEASYVERSION=$f22;
NGSUSER=$f23;
DATE=`date +"%d%m%y"`

########################################################################################################
## sample report name
#
SAMPLE_REPORT="${PROJECT_DIR}/${PROJECT_ID}/${SAMPLE_REPORT}/reports/${SAMPLE_ID}.NGSeasy.Report.md"
SAMPLE_REPORT_HTML="${PROJECT_DIR}/${PROJECT_ID}/${SAMPLE_REPORT}/reports/${SAMPLE_ID}.NGSeasy.Report.html"
SAMPLE_REPORT_PDF="${PROJECT_DIR}/${PROJECT_ID}/${SAMPLE_REPORT}/reports/${SAMPLE_ID}.NGSeasy.Report.pdf"
########################################################################################################
## NGSeasy list of files and stuff (keep dupemk.bam & final.bam)
#
BAM_PREFIX="${SAMPLE_ID}.${NGS_TYPE}.${DNA_PREP_LIBRARY_ID}.${NGS_PLATFORM}.${TRIM}.${ALIGNER}.${GENOMEBUILD}"
rawBAM="${BAM_PREFIX}.dupemk.bam"
finalBAM="${BAM_PREFIX}.dupemk.bam-bsqr.bam-realn.bam" 
# gatk-bsqr.gatk-realn.bam or no-bsqr.no-realn.bam based on config options or whats in dir then array of filenames
VCF="${BAM_PREFIX}.dupemk.bam-bsqr.bam-realn.bam.raw.platypus.vcf" 
# array of all files if multi callers used
FQ1=
FQ2=
peQCFQ1=
peQCFQ2=
seQCFQ1=
seQCFQ2=
# list qc reports ; check if exits and fail if they dont

########################################################################################################
## start report
#
tee ${SAMPLE_REPORT} <<EOF

NGSeasy Sample Report
======================

Author: ${FirstName} ${LastName}  
Email: ${NGSUSER}  
Date: ${RUNDATE}  
NGSeasy Version: ${NGSEASYVERSION}  

********

## Project

Project: ${PROJECT_ID}
Sample: ${SAMPLE_ID}
Location: ${PROJECT_DIR}
Config File: ${config_tsv}

## Config Options

Table options inc genome build

********

## Raw Fastq Files

Location

PE Illumina Platform Unit Quality Type

${FQ1}
${FQ2}

URL to FASTQC Report

### QC Report
Trimmomatic version
fastx version

summary of commands 

R1 ${NUMREADS1} ${NQCR1}
R2 ${NUMREADS2} ${NQCR2}

URL to FASTQC Report

**********

## Alignment Report

Aligner: bwa version 
Options[cmd]: bwa -M blah 
Location:
Bam:

Table QC report - picard tools & depth et
Bed Files
Plots
url link ucsc uplaod???

**********

## Variant Report
Variant Caller: platypus
Options[cmd]: 
Location:
VCF:
Misc Output:

**********

## Annotation Report
Annotator: snpeff version
Options[cmd]: 
Location:
VCF:
Misc Output:

**********

## Log Files

## END REPORT
#
EOF

########################################################################################################
## pandoc use rocker versions with new pandoc intalled
#
pandoc -s -o ${SAMPLE_REPORT_HTML} ${SAMPLE_REPORT}
pandoc ${SAMPLE_REPORT} -o ${SAMPLE_REPORT_PDF}

########################################################################################################
## END
done < ${RUNFILE}
chmod -R 765 ${PROJECT_DIR}/*