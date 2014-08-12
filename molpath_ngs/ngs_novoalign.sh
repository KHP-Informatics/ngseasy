#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

#
# ONE SCRIPT FOR ALL DATA TYPES
#
# PE/SE (autodetect)
# illumina only (for now)
# no trimming (done by trimmomatic)
#

####################
## Call Novoalign ##
####################

# get filtered fastq from sample directory
sample_name=${1}
sample_dir=${2}

cd ${sample_dir}

# get input fastQ (filtered and paired if appropriate)
FILES=(`ls *.filtered.fastq`)
if [[ ${#FILES} -eq 2 ]];then
    sample_librarytype='PE'
    IFS=$' ' SORTED=($(sort <<< "${FILES[0]} ${FILES[1]}"))
    FQ1=${SORTED[0]}
    FQ2=${SORTED[1]}
elif [[ ${#FILES} -eq 1 ]];then
    sample_librarytype='SE'
    FQ1=${FILES[0]}
else
    # is something else
    echo "Didn't find the right FastQ files"
    exit 1
fi

# outputfile
OUTFILE=${sample_name}.aln.sam

echo "----------------------------------------------------------------------------------------"
echo " sampleDir        " $sample_dir
echo " CWD              " `pwd`
echo " in FASTQ         " $FQ1 "<>" $FQ2
echo " Reads are        " ${sample_librarytype}
echo " alignment CPUs   " ${novo_cpu}
echo " Reference Genome " ${reference_genome_novoindex}
echo " output SAM       " ${OUTFILE}
echo "----------------------------------------------------------------------------------------"

if [[ $sample_librarytype -eq 'PE' ]];then
    ${ngs_novo}/novoalign \
    -c ${novo_cpu} \
    -d ${reference_genome_novoindex} \
    -f ${FQ1} ${FQ2} \
    -F STDFQ \
    --Q2Off \
    --3Prime \
    -g 40 \
    -x 6 \
    -r All \
    -i PE 200,50 \
    -k -K ${sample_name}.novoalign_calibration.stats \
    -o SAM > ${OUTFILE};
elif [[ $sample_librarytype -eq 'SE']];then
    ${ngs_novo}/novoalign \
    -c ${novo_cpu} \
    -d ${reference_genome_novoindex} \
    -f ${FQ1} \
    -F STDFQ \
    --Q2Off \
    --3Prime  \
    -g 40 \
    -x 6 \
    -r All \
    -k -K s${sample_name}.novoalign_calibration.stats \
    -o SAM > ${OUTFILE};
else
    echo "unkown PE/SE?"
    exit 1
fi


















