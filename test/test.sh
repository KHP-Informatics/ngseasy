#!/bin/bash 

############################################
## RUN NGSEASY                             ##
#############################################
#
# Default install dir is user's /home/${USER}

#############################################
## 0. Move to config file dir

NGSEASY_DIR=${1}

cd ${NGSEASY_DIR}/ngs_projects/config_files/

#############################################
## 1. Set up project and sample directories

ngseasy_initiate_project -c ngseasy_test.config.tsv -d ${NGSEASY_DIR}//ngs_projects 

#############################################
## 2. Move project/sample fastq from raw_fastq
## to project and sample directories
 
ngseasy_initiate_fastq -c ngseasy_test.config.tsv -d ${NGSEASY_DIR}/ngs_projects 

#############################################
## 3. Run basic test 

ngseasy -c ngseasy_test.config.tsv -d ${NGSEASY_DIR}/ngs_projects 

