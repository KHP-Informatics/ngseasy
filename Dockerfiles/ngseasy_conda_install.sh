#!/usr/bin/env bash

set -o errexit
set -o pipefail
#set -o nounset
#set -o xtrace

##----------------------------------------------------------------------------##
## version
VERSION="ngseasy-v0.1-conda"

##----------------------------------------------------------------------------##
## home and user
MYOS=`uname`
MYHOME=${HOME}
GETUSER=`whoami`
ARCH=`uname -m`
echo ""
echo "-------------------------------------------------------------------------"
echo "NGSeasy (conda) Version : ${VERSION}"
echo "-------------------------------------------------------------------------"
echo "Installs Anaconda2-4 to local system along with a suite of NGS tools"
echo "WARNING: Currently only supports x86_64 Linux and Darwin (MAC OSX)"
echo "WARNING: Some important tools are missing for mac osx-64"
echo "Contact: stephen.j.newhouse@gmail.com"
echo "-------------------------------------------------------------------------"
echo ""
echo "User: ${GETUSER}"
echo "Home: ${MYHOME}"

##----------------------------------------------------------------------------##
## check arch
if [[ "${ARCH}" == "x86_64" ]]; then
  echo "OS: ${MYOS} ${ARCH}"
  echo ""
else
  echo "ERROR: Currently only supports x86_64 Linux and Darwin (MAC OSX)"
  echo "Exiting"
  sleep 2s
  exit 1
fi

##----------------------------------------------------------------------------##
## Set Install directory. Default is HOME
if [[ -z "${1}"  ]]; then
  INSTALL_DIR="${HOME}"
  if [ "${INSTALL_DIR}" == "/home/root" ]; then
    echo "ERROR: trying to install to /home/root not permitted"
    echo "are you really root?"
    echo "Exiting"
    sleep 3s
    exit 1
  fi
  echo "RUNNING: ngseasy_conda_install.sh ${1}"
  echo "WARNING: No INSTALL_DIR specified"
  echo "WARNING:Install directory  will be set to default /home/user [${INSTALL_DIR}]"
  echo ""
  echo "Usage: bash ngseasy_conda_install.sh /PATH/TO/INSTALL/DIR"
  sleep 1s
else
  INSTALL_DIR="${1}"
  echo "Install directory set to [${INSTALL_DIR}]"
fi

##----------------------------------------------------------------------------##
## Install anaconda2
cd ${INSTALL_DIR}

if [[ "${MYOS}" == "Linux" ]]; then
  CONDA=""
  CONDA="Anaconda2-4.0.0-Linux-x86_64.sh"
  echo "Anaconda2-4.0.0 is being installed to [${INSTALL_DIR}/anaconda2]"
  wget --quiet https://repo.continuum.io/archive/${CONDA} && \
  /bin/bash ./${CONDA} -b -p ${INSTALL_DIR}/anaconda2 && \
  rm -v ./${CONDA}
  unset CONDA
  # add conda bin to path
  export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin

elif [[  "${MYOS}" == "Darwin"  ]]; then
  CONDA=""
  CONDA="Anaconda2-4.0.0-MacOSX-x86_64.sh"
  echo "Anaconda2-4.0.0 is being installed to [${INSTALL_DIR}/anaconda2]"
  wget --quiet  https://repo.continuum.io/archive/${CONDA} && \
  /bin/bash ./${CONDA} -b -p ${INSTALL_DIR}/anaconda2 && \
  rm -v ./${CONDA}
  unset CONDA
  # add conda bin to path
  export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin

else
  echo "ERROR: No OS detected"
  echo "Exiting"
  sleep 2s
  exit 1
fi

if [[ -x  "${INSTALL_DIR}/anaconda2/bin/conda" ]]; then
# setup conda
which conda
echo "Start conda set up and updates"
${INSTALL_DIR}/anaconda2/bin/conda update -y conda
${INSTALL_DIR}/anaconda2/bin/conda update -y conda-build
${INSTALL_DIR}/anaconda2/bin/conda update -y --all
mkdir -p ${INSTALL_DIR}/anaconda2/conda-bld/linux-64 ${INSTALL_DIR}/anaconda2/conda-bld/osx-64
${INSTALL_DIR}/anaconda2/bin/conda index ${INSTALL_DIR}/anaconda2/conda-bld/linux-64 ${INSTALL_DIR}/anaconda2/conda-bld/osx-64

## add channels
echo "add channels bioconda r sjnewhouse"
${INSTALL_DIR}/anaconda2/bin/conda config --add channels bioconda
${INSTALL_DIR}/anaconda2/bin/conda config --add channels r
${INSTALL_DIR}/anaconda2/bin/conda config --add channels sjnewhouse

##----------------------------------------------------------------------------##
## ngs tools
echo "get ngs tool list"
wget https://raw.githubusercontent.com/KHP-Informatics/ngseasy/f1000_dev/Dockerfiles/ngs_conda_tool_list.txt
wget https://raw.githubusercontent.com/KHP-Informatics/ngseasy/f1000_dev/Dockerfiles/ngs_conda_tool_list_osx_64_may2016.txt
# some tools are not built for osx yet. quikc fix is to remove them from list
# I will see if we can build thiese tools for osx-64 conda...coming soon
# some pretty important tools
#  - parallel
#  - biobambam
#  - samblaster
#  - vt
#  - vawk
#  - pysamstats
#  - trim-galore
#  - snap-aligner
#  - khmer

echo "Install ngseasy tools"
${INSTALL_DIR}/anaconda2/bin/conda install --yes --file ngs_conda_tool_list.txt
rm -v ./ngs_conda_tool_list.txt

else
  echo "ERROR: can not find ${INSTALL_DIR}/anaconda2/bin/conda"
  which conda
  echo "WARNING: did you install Anaconda2-4"
  echo "Exiting"
  sleep 3s
  exit 1
fi

# list
TIME_STAMP=`date +"%d-%m-%y"`
${INSTALL_DIR}/anaconda2/bin/conda list -e > ${INSTALL_DIR}/anaconda2/ngseasy-spec-file-${TIME_STAMP}.txt
conda info
unset INSTALL_DIR
unset MYOS
unset MYHOME
unset GETUSER
unset ARCH

##----------------------------------------------------------------------------##
## The end
echo "Done installing ngseasy tools [Version: ${VERSION}]"
unset VERSION
