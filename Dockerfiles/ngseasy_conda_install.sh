#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset
#set -o xtrace

##----------------------------------------------------------------------------##
## version
VERSION="ngseasy-v0.1-conda"

##----------------------------------------------------------------------------##
## home and user
MYOS=$(uname)
MYHOME=${HOME}
GETUSER=$(whoami)
ARCH=$(uname -m)
CONDA_URL="https://repo.continuum.io/archive"
CONDA_VERSION="4.2.0"
PYTHON_VERSION="2"
BIN_DIR="anaconda${PY_VERSION}"

##----------------------------------------------------------------------------##
## opening message
echo ""
echo "-------------------------------------------------------------------------"
echo "NGSeasy (conda) Version : ${VERSION}"
echo "-------------------------------------------------------------------------"
echo "Installs Anaconda2-4 to local system along with a suite of NGS tools"
echo "WARNING: Currently only supports x86_64 Linux"
echo "NOTE: Some important conda builds are missing for mac osx-64"
echo "Contact: stephen.j.newhouse@gmail.com"
echo "-------------------------------------------------------------------------"
echo ""
echo "INFO: User: ${GETUSER}"
echo "INFO: Home: ${MYHOME}"

##----------------------------------------------------------------------------##
## check OS
if [[ "${MYOS}" == "Linux" ]]; then
  echo "INFO: ${MYOS}" == "Linux"
elif [[  "${MYOS}" == "Darwin"  ]]; then
  echo "INFO: ${MYOS}" == "Darwin"
  echo "INFO: Currently only supports x86_64 Linux"
  echo "INFO: Exiting"
  exit 1
else
  echo "ERROR: No OS detected"
  echo "INFO: Exiting"
  exit 1
fi

##----------------------------------------------------------------------------##
## check arch
if [[ "${ARCH}" == "x86_64" ]]; then
  echo "INFO: OS: ${MYOS} ${ARCH}"
else
  echo "ERROR: Currently only supports x86_64 Linux"
  echo "Exiting"
  exit 1
fi


##----------------------------------------------------------------------------##
## Set Install directory. Default is /opt
if [[ $# -eq 0  ]]; then
  echo "WARNING: no arguments set"
  INSTALL_DIR="/opt"
  if [ "${INSTALL_DIR}" == "/home/root" ] || [ "${INSTALL_DIR}" == "/root" ] ; then
    echo "ERROR: trying to install to /home/root not permitted"
    echo "WARNING: are you really root?"
    echo "INFO: Exiting"
    exit 1
  fi
  echo "INFO: RUNNING: ngseasy_conda_install.sh ${INSTALL_DIR}"
  echo "WARNING: No INSTALL_DIR specified"
  echo "WARNING: Install directory  will be set to default /opt"
  echo "USEAGE: bash ngseasy_conda_install.sh /PATH/TO/INSTALL/DIR"
  sleep 1s
else
  INSTALL_DIR="${1}"
  echo "INFO: Install directory set to [${INSTALL_DIR}]"
fi

##----------------------------------------------------------------------------##
## Install anaconda

## move to install dir
cd ${INSTALL_DIR}

## set version name
CONDA="Anaconda${PY_VERSION}-${CONDA_VERSION}-Linux-x86.sh"
echo "Anaconda${PY_VERSION} is being installed to [${INSTALL_DIR}/anaconda${PY_VERSION}]"

## get install script from https://repo.continuum.io/archive/
  wget -N https://repo.continuum.io/archive/${CONDA} && \
  /bin/bash ./${CONDA} -b -p ${INSTALL_DIR}/anaconda${PY_VERSION} && \
  rm -v ./${CONDA}
  unset CONDA

## add conda bin to path
echo "INFO: export PATH=$PATH:${INSTALL_DIR}/anaconda${PY_VERSION}/bin"
export PATH=$PATH:${INSTALL_DIR}/anaconda${PY_VERSION}/bin


##----------------------------------------------------------------------------##
# setup conda
if [[ -x  "${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda" ]]; then
which conda
echo "INFO: Start conda set up and updates"

# run cmd
${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda update -y conda
${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda update -y conda-build
${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda update -y --all

## add channels
echo "INFO: add channels bioconda r sjnewhouse"

${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda config --add channels bioconda
${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda config --add channels r
${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda config --add channels sjnewhouse

##----------------------------------------------------------------------------##
## ngs tools
echo "INFO: get ngs tool list"
curl -L https://raw.githubusercontent.com/KHP-Informatics/ngseasy/f1000_dev/Dockerfiles/ngs_conda_tool_list.txt \
-o ${INSTALL_DIR}//ngs_conda_tool_list.txt

echo "INFO: Install ngseasy tools"
${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda install --yes --file ${INSTALL_DIR}/ngs_conda_tool_list.txt
wait
${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda clean -tipsy
wait
rm -v ./ngs_conda_tool_list.txt

else
  echo "ERROR: can not find ${INSTALL_DIR}/anaconda${PY_VERSION}/bin/conda"
  echo "WARNING: did you install Anaconda"
  echo "INFO: Exiting"
  exit 1
fi

##----------------------------------------------------------------------------##
# list tools and versions
TIME_STAMP=`date +"%d-%m-%y"`
conda info > ${INSTALL_DIR}/anaconda2/ngseasy-spec-file-${TIME_STAMP}.txt
${INSTALL_DIR}/anaconda2/bin/conda list -e >> ${INSTALL_DIR}/anaconda2/ngseasy-spec-file-${TIME_STAMP}.txt
unset INSTALL_DIR
unset MYOS
unset MYHOME
unset GETUSER
unset ARCH

##----------------------------------------------------------------------------##
## The end
echo "INFO: Done installing ngseasy tools [Version: ${VERSION}]"
unset VERSION
