


logger_ngseasy
logger_ngseasy_new
ngseasy

ngseasy_alignment_v0.9




ngseasy_functions



ngseasy_qcfiler_bam
ngseasy_filter_recalbam


ngs_full_gatk
ngs_full_no_gatk



## Step 1

## get contianers 
get_containers.sh

## step up directory structure 

ngseasy_initiate_project 

## copy fastq file from storage folder to project and sample folders
ngseasy_initiate_fastq

## start volumes container
ngseasy_volumes_container

### TO DO 
add checks for containers/images 

****

## NGS pipeline

```
#!/bin/bash -e -x

################################################################
# Program: ngseasy
# Version 1.0 
# Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)
#################################################################

## NGSeasy version
NGSEASYVERSION="1.0"

## check and make ~/ngseasy_logs if needed
if [[ ! -e  ${HOME}/ngseasy_logs ]]
then
  mkdir ${HOME}/ngseasy_logs
  global_run_logs="${HOME}/ngseasy_logs"
fi

## global logging fuction
function logger_ngseasy() {
 message=${1}
 mylogfile=${2}
 echo -e [`date`]":[NGSEASY:${NGSEASYVERSION}]:"${message}":[`pwd`]:[${USER}]:[`uname -a`]" >> ${mylogfile}.log >> ${global_run_logs}/ngseasy-run.log;
 echo -e [`date`]":[NGSEASY:${NGSEASYVERSION}]:"${message}":[`pwd`]:[${USER}]:[`uname -a`]"
}

## global usage
function usage_ngseasy() {
    echo "
Program: ngseasy
Version 1.0
Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)

usage:   ngseasy_full -c <config_file> -d <project_directory>

options:  -c  configuration file
          -d  project directory
          -h  show this message 
"
}

## example config: https://docs.google.com/spreadsheets/d/1VWqmMffkVDnvOtRJGlPqOYzXWnIN_IONXQHDAawaN5Q/edit#gid=0

## Check options passed in.
if test -z "$2"
  then
  usage_ngseasy
  exit 1
fi


## get options for command line args
while  getopts "hc:d:" opt
do

  case ${opt} in

   h)
   usage_ngseasy #print help
   exit 0
   ;;

   c)
   config_tsv=${OPTARG}
   echo "-c = ${config_tsv}"
   ;;

   d)
   project_directory=${OPTARG}
   echo "-d = ${project_directory}"
   ;;


   esac
done

## check file and directory exist.
if [[ ! -e "${config_tsv}" ]] 
  then
	  usage_ngseasy;
	  echo -e "ERROR : ${config_tsv} does not exist\n"
	  exit 1;
fi

## check exists.
if [[ ! -d "${project_directory}" ]] 
  then
	  usage_ngseasy;
	  echo -e "ERROR :  ${project_directory} does not exist\n"
	  exit 1;
fi


# --- Start NGS Pipeline -------------------------------------------------------- #



## params to get
## trim gatk bsqr realign 
##

while read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17
do

## set varibales  
  DATE=`date +"%d%m%y"`
  
POJECT_ID=f1
SAMPLE_ID f2
FASTQ1  f3
FASTQ2  f4
PROJECT_DIR f5
DNA_PREP_LIBRARY_ID f6
NGS_PLATFORM  f7
NGS_TYPE  f8
BAIT  f9
CAPTURE f10
TRIM  f11
BSQR  f12
REALN f13
ALIGNER f14
VARCALLER f15
CNV f16
ANNOTATOR f17
CLEANUP f18
NCPU  f19
VERSION f20
NGSUSER f21

  POJECT_ID=$f1
  SAMPLE_ID=$f2
  FASTQ1=$f3
  FASTQ2=$f4
  PROJECT_DIR=$f5 
  DNA_PREP_LIBRARY_ID=$f6
  NGS_PLATFORM=$f7
  NGS_TYPE=$f8
  BAIT=
  CAPTURE=

  TRIM=$f11
  BSQR=$f13
  REALN=$f14
  ALIGNER=$f15
  VARCALLER=$f16
  CNV=
  CLEANUP=$f18
  NCPU=$f19
  VERSION=$f20
  NGSUSER=$f21


# Read config file 
logger_ngseasy "[ngseasy]:Reading [${config_tsv}] " ${HOME}/ngseasy_logs/ngseasy.${POJECT_ID}.${USER}.$(date +"%d%m%y"




## fastqc
ngseasy_fastqc -c ${config_tsv} -d ${project_directory}

##-------------------------------------------------------------------------##
## adapter and read/base quality trimming
if [[ "${TRIM}" -eq 1 ]]
then
  ngseasy_trimmomatic -c ${config_tsv} -d ${project_directory}
fi

##-------------------------------------------------------------------------##
## alignment : includes addition of read groups at alignment stage 
## and then duplicate marking (samblaster), indexing and sorting with sambamba
ngseasy_alignment -c ${config_tsv} -d ${project_directory}


##-------------------------------------------------------------------------##
## NGS Processing : Indel realignment and base quality score reclibration using GATK or BamUtil/ogap
if [[ "${GATK}" -eq 1 ]] && [[ "${REALN}" -eq 1 ]] && [[ "${BSQR}" -eq 1 ]]
then

  ngseasy_indel_realn -c ${config_tsv} -d ${project_directory}
  ngseasy_base_recal -c ${config_tsv} -d ${project_directory}
    
elif [[ "${GATK}" -eq 1 ]] && [[ "${REALN}" -eq 0 ]] && [[ "${BSQR}" -eq 1 ]]
then

  ngseasy_base_recal -c ${config_tsv} -d ${project_directory}

elif [[ "${GATK}" -eq 0 ]] && [[ "${REALN}" -eq 1 ]] && [[ "${BSQR}" -eq 1 ]]
then

  ngseasy_ogap_realn -c ${config_tsv} -d ${project_directory}
  ngseasy_bamutil_base_recal -c ${config_tsv} -d ${project_directory}

elif [[ "${GATK}" -eq 0 ]] && [[ "${REALN}" -eq 0 ]] && [[ "${BSQR}" -eq 1 ]]
then

  ngseasy_bamutil_base_recal -c ${config_tsv} -d ${project_directory}
  
fi

##-------------------------------------------------------------------------##
## Alignment statistics
ngseasy_alignment_qc -c ${config_tsv} -d ${project_directory}

##-------------------------------------------------------------------------##
## SNP/INDEL calling
ngseasy_variant_calling -c ${config_tsv} -d ${project_directory}
 
ngseasy_variant_calling_fast_ensemble -c ${config_tsv} -d ${project_directory}

##-------------------------------------------------------------------------##
## CNV Calling

##-------------------------------------------------------------------------##
## Annotation

##-------------------------------------------------------------------------##
## NGS Report

```

## pipelines
gatk_realn_recab
gatk_recab




## options
POJECT_ID
SAMPLE_ID
FASTQ1
FASTQ2
PROJECT_DIR
DNA_PREP_LIBRARY_ID
NGS_PLATFORM
NGS_TYPE
BAIT
CAPTURE
TRIM
GATK
BSQR
REALN
ALIGNER
VARCALLER
CNV # if WGS cn.MOPs mHMM DELLY and LUMPY if WEX ExomeDepth DELLY and LUMPY if TGS SLOPE DELLY and LUMPY 
GTMODEGATK
CLEANUP
NCPU
VERSION
NGSUSER


****
ngseasy_full
ngseasy_fastqc
ngseasy_trimmomatic
ngseasy_alignment
ngseasy_indel_realn
ngseasy_base_recal
ngseasy_ogap_realn
ngseasy_bamutil_base_recal
ngseasy_variant_calling

## Dumped

```
ngseasy_addreadgroup
ngseasy_markduplicates
```


*****


```
#usage printing func
usage()
{
cat << EOF
  This script calls the NGSeasy pipeline : <full_gatk/full_no_gatk/fastqc/fastq_trimm/alignment/var_call/cnv_call/var_annotate/alignment_qc>
  
  See NGSEasy containerized instructions.

  ARGUMENTS:
  
  -h      Flag: Show this help message
  -c      NGSeasy project and run configureation file
  -d      Base directory for (fastq_raw, reference_genomes_b37, gatk_resources, ngs_projects, ngseasy_scripts)

  EXAMPLE USAGE:
    
    ngseasy -c config.file.tsv -d /media/ngs_projects 

EOF
}

#get options for command line args
  while  getopts "hc:d:" opt
  do

      case ${opt} in
	  h)
	  usage #print help
	  exit 0
	  ;;
	  
	  c)
	  config_tsv=${OPTARG}
	  echo "-c = ${config_tsv}"
	  ;;
	 
	  d)
	  project_directory=${OPTARG}
	  echo "-d = ${project_directory}"
	  ;; 
      esac
  done

#check exists.
if [ ! -e "${config_tsv}" ] 
  then
	  echo "ERROR : ${config_tsv} does not exist "
	  usage;
	  exit 1;
fi

#check exists.
if [ ! -d "${project_directory}" ] 
  then
	  echo "ERROR :  ${project_directory} does not exist "
	  usage;
	  exit 1;
fi
```






