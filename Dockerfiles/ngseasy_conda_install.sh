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
echo "User: ${GETUSER}"
echo "Home: ${MYHOME}"

##----------------------------------------------------------------------------##
## check OS
if [[ "${MYOS}" == "Linux" ]]; then
  echo "${MYOS}" == "Linux"
elif [[  "${MYOS}" == "Darwin"  ]]; then
  echo "${MYOS}" == "Darwin"
  echo "Currently only supports x86_64 Linux"
  echo "Exiting"
  exit 1
else
  echo "ERROR: No OS detected"
  echo "Exiting"
  exit 1
fi

##----------------------------------------------------------------------------##
## check arch
if [[ "${ARCH}" == "x86_64" ]]; then
  echo "OS: ${MYOS} ${ARCH}"
else
  echo "ERROR: Currently only supports x86_64 Linux"
  echo "Exiting"
  exit 1
fi


##----------------------------------------------------------------------------##
## Set Install directory. Default is /opt
if [[ $# -eq 0  ]]; then
  echo "WARNING: no arguments set. Using opt as default install directory"
  INSTALL_DIR="/opt"
  if [ "${INSTALL_DIR}" == "/home/root" ] || [ "${INSTALL_DIR}" == "/root" ] ; then
    echo "ERROR: trying to install to /home/root not permitted"
    echo "are you really root?"
    echo "Exiting"
    exit 1
  fi
  echo "RUNNING: ngseasy_conda_install.sh ${1}"
  echo "WARNING: No INSTALL_DIR specified"
  echo "WARNING:Install directory  will be set to default /opt"
  echo "USEAGE: bash ngseasy_conda_install.sh /PATH/TO/INSTALL/DIR"
  sleep 1s
else
  INSTALL_DIR="${1}"
  echo "Install directory set to [${INSTALL_DIR}]"
fi

##----------------------------------------------------------------------------##
## Install anaconda2
cd ${INSTALL_DIR}

  CONDA=""
  CONDA="Anaconda2-4.1.1-Linux-x86_64.sh"
  echo "Anaconda2-4.0.0 is being installed to [${INSTALL_DIR}/anaconda2]"
  echo "wget -N https://repo.continuum.io/archive/${CONDA} && \
    /bin/bash ./${CONDA} -b -p ${INSTALL_DIR}/anaconda2 && \
    rm -v ./${CONDA}"

  wget -N --quiet https://repo.continuum.io/archive/${CONDA} && \
  /bin/bash ./${CONDA} -b -p ${INSTALL_DIR}/anaconda2 && \
  rm -v ./${CONDA}
  unset CONDA
  # add conda bin to path
  echo "export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin"
  export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin


##----------------------------------------------------------------------------##
# setup conda
if [[ -x  "${INSTALL_DIR}/anaconda2/bin/conda" ]]; then
which conda
echo "Start conda set up and updates"

# run cmd
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
curl -L https://raw.githubusercontent.com/KHP-Informatics/ngseasy/f1000_dev/Dockerfiles/ngs_conda_tool_list.txt \
-o ${INSTALL_DIR}//ngs_conda_tool_list.txt

echo "Install ngseasy tools"
${INSTALL_DIR}/anaconda2/bin/conda install --yes --file ${INSTALL_DIR}/ngs_conda_tool_list.txt
wait
${INSTALL_DIR}/anaconda2/bin/conda clean -tipsy
wait
rm -v ./ngs_conda_tool_list.txt

else
  echo "ERROR: can not find ${INSTALL_DIR}/anaconda2/bin/conda"
  echo "WARNING: did you install Anaconda2-4"
  echo "Exiting"
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
echo "Done installing ngseasy tools [Version: ${VERSION}]"
unset VERSION
