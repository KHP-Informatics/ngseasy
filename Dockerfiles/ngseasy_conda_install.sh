#!/usr/bin/env bash
set -o errexit
set -o pipefail
#set -o nounset
set -o xtrace

## version
VERSION="ngseasy-v0.1-conda"

## home and user
MYOS=`uname`
MYHOME=${HOME}
GETUSER=${USER}
ARCH=`uname -m`
echo "User: ${GETUSER}"
echo "Home: ${MYHOME}"
echo "OS: ${MYOS} ${ARCH}"

## Set Install directory. Default is HOME
if [[ -z "${1}"  ]]; then
  INSTALL_DIR="${HOME}"
  echo "No arguments supplied"
  echo "Install directory set to default /home/user [${INSTALL_DIR}]"
else
  INSTALL_DIR="${1}"
  echo "Install directory set to [${INSTALL_DIR}]"
fi


## Install anaconda2
cd ${INSTALL_DIR}

if [[ "${MYOS}" == "Linux" ]]; then
  CONDA=""
  CONDA="Anaconda2-4.0.0-Linux-x86_64.sh"
  echo "Anaconda2-4.0.0 is being installed to [${INSTALL_DIR}/anaconda2]"
  wget https://repo.continuum.io/archive/${CONDA} && \
  /bin/bash ./${CONDA} -b -p ${INSTALL_DIR}/anaconda2 && \
  rm -v ./${CONDA}
  unset CONDA
  # add conda bin to path
  echo "running: export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin"
  export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin

elif [[  "${MYOS}" == "Darwin"  ]]; then
  CONDA=""
  CONDA="Anaconda2-4.0.0-MacOSX-x86_64.sh"
  echo "Anaconda2-4.0.0 is being installed to [${INSTALL_DIR}/anaconda2]"
  wget https://repo.continuum.io/archive/${CONDA} && \
  /bin/bash ./${CONDA} -b -p ${INSTALL_DIR}/anaconda2 && \
  rm -v ./${CONDA}
  unset CONDA=
  # add conda bin to path
  echo "running: export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin"
  export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin

else
   echo "Error: No OS detected"
   echo "Exiting"
   sleep 2s
   exit 1
fi

# source dotfile
echo "Souce dotfile"
if [[ "${MYOS}" == "Linux" ]]; then
  touch ${INSTALL_DIR}/.conda_bin
  echo 'PATH=$PATH:${INSTALL_DIR}/anaconda2/bin' >> ${INSTALL_DIR}/.conda_bin ## linux
  /bin/bash -c "source ${INSTALL_DIR}/.conda_bin"

elif [[  "${MYOS}" == "Darwin"  ]]; then
   touch ${INSTALL_DIR}/.conda_bin
   echo 'PATH=$PATH:${INSTALL_DIR}/anaconda2/bin' >> ${INSTALL_DIR}/.conda_bin ## mac
  /bin/bash -c "source ${INSTALL_DIR}/.conda_bin"

else
   echo "Error: No OS detected"
   echo "Exiting"
   sleep 2s
   exit 1
fi


# setup conda
echo "Start conda set up"
${INSTALL_DIR}/anaconda2/bin/conda update -y conda
${INSTALL_DIR}/anaconda2/bin/conda update -y conda-build
${INSTALL_DIR}/anaconda2/bin/conda update -y --all
mkdir -p ${INSTALL_DIR}/anaconda2/conda-bld/linux-64 ${INSTALL_DIR}/anaconda2/conda-bld/osx-64
${INSTALL_DIR}/anaconda2/bin/conda index ${INSTALL_DIR}/anaconda2/conda-bld/linux-64 ${INSTALL_DIR}/anaconda2/conda-bld/osx-64

## add channels
echo "add channels"
${INSTALL_DIR}/anaconda2/bin/conda config --add channels bioconda
${INSTALL_DIR}/anaconda2/bin/conda config --add channels r
${INSTALL_DIR}/anaconda2/bin/conda config --add channels sjnewhouse


## ngs tools
echo "get ngs tool list"
wget https://raw.githubusercontent.com/KHP-Informatics/ngseasy/f1000_dev/Dockerfiles/ngs_conda_tool_list.txt

## create ngseasy environment python >=2.7 for ngs tools
echo "create ngseasy environment"
${INSTALL_DIR}/anaconda2/bin/conda create --yes --name ngseasy --file ngs_conda_tool_list.txt
rm -v ./ngs_conda_tool_list.txt

## activate
echo "activate ngseasy environment"
 /bin/bash -c "source ${INSTALL_DIR}/anaconda2/bin/activate ngseasy"

## update nextflow
nextflow self-update

# list
TIME_STAMP=`date +"%d-%m-%y"`
${INSTALL_DIR}/anaconda2/bin/conda list -e > ${INSTALL_DIR}/anaconda2/ngseasy-spec-file-${TIME_STAMP}.txt
unset INSTALL_DIR
unset MYOS
unset MYHOME
unset GETUSER
unset ARCH
## The end
echo "Done installing ngseasy tools [Version: ${VERSION}]"
unset VERSION
