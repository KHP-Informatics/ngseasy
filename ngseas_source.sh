#!/usr/bin/env bash
set -o errexit
set -o pipefail
#set -o nounset
#set -o xtrace

function ngseasy_fastqc(){

local EXEC=`which fastqc || true`

cd $PROJECT/${SAMPLE}

mkdir -p $PROJECT/TEMP-FASTQC
mkdir -p fastq

local TEMP_DIR=$PROJECT/TEMP-FASTQC

local FQ1=${FASTQ1}
local FQ2=${FASTQ2}
local OUTDIR="./fastq"

local OPTIONS=`--threads 2 --extract --dir ${TEMP_DIR} --outdir ${OUTDIR}`

local CMD=`${${EXEC} ${OPTIONS} ${FASTQ1} ${FASTQ2}`

echo "RUN: ${CMD}"
${CMD}

rm -rf $PROJECT/TEMP-FASTQC
