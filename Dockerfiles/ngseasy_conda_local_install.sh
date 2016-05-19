#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset
# set -o xtrace

## version
VERSION="0.1"

INSTALL_DIR="/home/ngseasy"

echo "Installing to [${INSTALL_DIR}]"

cd ${INSTALL_DIR} && pwd

wget https://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh && \
/bin/bash ./Anaconda2-4.0.0-Linux-x86_64.sh -b -p ${INSTALL_DIR}/anaconda2 && \
rm -v ./Anaconda2-4.0.0-Linux-x86_64.sh
echo "Anaconda is being installed at [${INSTALL_DIR}/anaconda2]"

# add conda bin to path
export PATH=$PATH:${INSTALL_DIR}/anaconda2/bin
echo 'PATH=$PATH:${INSTALL_DIR}/anaconda2/bin' >> /home/ngseasy/.bashrc

# setup conda
${INSTALL_DIR}/anaconda2/bin/conda update -y conda
${INSTALL_DIR}/anaconda2/bin/conda update -y conda-build
${INSTALL_DIR}/anaconda2/bin/conda update -y --all
mkdir -p ${INSTALL_DIR}/anaconda2/conda-bld/linux-64 ${INSTALL_DIR}/anaconda2/conda-bld/osx-64
${INSTALL_DIR}/anaconda2/bin/conda index ${INSTALL_DIR}/anaconda2/conda-bld/linux-64 ${INSTALL_DIR}/anaconda2/conda-bld/osx-64

## add channels
${INSTALL_DIR}/anaconda2/bin/conda config --add channels bioconda
${INSTALL_DIR}/anaconda2/bin/conda config --add channels r
${INSTALL_DIR}/anaconda2/bin/conda config --add channels sjnewhouse


## ngs tools
wget https://raw.githubusercontent.com/KHP-Informatics/ngseasy/f1000_dev/Dockerfiles/ngs_conda_tool_list.txt

## create ngseasy environment python >=2.7 for ngs tools
${INSTALL_DIR}/anaconda2/bin/conda create --yes --name ngseasy --file ngs_conda_tool_list.txt

## activate
 /bin/bash -c "source ${INSTALL_DIR}/anaconda2/bin/activate ngseasy"

rm -v ./ngs_conda_tool_list.txt

## update nextflow
nextflow self-update

# list
TIME_STAMP=`date +"%d-%m-%y"`
${INSTALL_DIR}/anaconda2/bin/conda list -e > ${INSTALL_DIR}/anaconda2/ngseasy-spec-file-${TIME_STAMP}.txt
