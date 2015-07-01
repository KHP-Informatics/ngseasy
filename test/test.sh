#!/bin/bash 

############################################
## RUN NGSEASY                             ##
#############################################
#
# Default install dir is user's /home/${USER}

NGSEASY_DIR=${1}

#############################################
## 0. Move to config file dir

cd ${NGSEASY_DIR}/ngs_projects/config_files/

#############################################
## 1. Run basic test 

ngseasy -c ngseasy_test.config.tsv -d ${NGSEASY_DIR}/ngs_projects 

