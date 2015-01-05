


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
# Program: ngseasy_full
# Version 1.0 
# Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)
#################################################################

## global usage
function usage_ngseasy_full() {
    echo "
Program: ngseasy_full
Version 1.0
Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)

usage:   ngseasy_full -c <config_file> -d <project_directory> [options]

options:  -t  1=quality trim [default], 0=no quality trimming 
          -g  1=GATK, 0=No GATK [default] 
          -r  1=Indel Realignment, 0=No Indel Realignment [default]
          -h  show this message 
"
}


## set defaults
TRIM=1
GATK=0
REALN=1
BSQR=1
project_directory=""
config_tsv=""

## Check options passed in.
    if test -z "$2"
    then
	usage_ngseasy_full
	exit 1
    f


## get options for command line args
while  getopts "hc:d:tgr" opt
do

  case ${opt} in

   h)
   usage_ngseasy_full #print help
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

   t)
   TRIM=${OPTARG}
   echo "-t = ${TRIM}"
   ;; 
 
   g)
   GATK=${OPTARG}
   echo "-d = ${GATK}"
   ;; 
 
   r)
   REALN=${OPTARG}
   echo "-d = ${REALN}"
   ;; 
   
  esac

done

## check file and directory exist.
if [[ ! -e "${config_tsv}" ]] 
  then
	  usage_ngseasy_full;
	  echo -e "ERROR : ${config_tsv} does not exist\n"
	  exit 1;
fi

## check exists.
if [[ ! -d "${project_directory}" ]] 
  then
	  usage_ngseasy_full;
	  echo -e "ERROR :  ${project_directory} does not exist\n"
	  exit 1;
fi


# --- Start NGS Pipeline --- #

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
CNV
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






