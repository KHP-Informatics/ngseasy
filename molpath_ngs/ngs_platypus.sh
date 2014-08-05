#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V


###TODO
# check insert size
# is the proper paired flag used by GATK/platypus





###################################
### Platypus genotype calling
###################################
#
# Takes direct output of stampy mapping
#
#python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fa --output=out.vcf --filterDuplicates=0

# get filtered fastq from sample directory
sample_name=${1}
sample_dir=${2}

cd ${sample_dir}

# outputfile
INFILE=${sample_name}.stampy.sam
OUTFILE=${sample_name}.platypus.vcf

echo "----------------------------------------------------------------------------------------"
echo " sampleDir        " $sample_dir
echo " CWD              " `pwd`
echo " Reference Genome " ${reference_genome_novoindex}
echo " BWA output BAM   " ${BWAFILE}
echo " FINAL output SAM " ${OUTFILE}
echo "----------------------------------------------------------------------------------------"

# for targeted data
if [[ ${analysisType} == *panel* ]];then
    python ${ngs_platypus}/Platypus.py callVariants \
    --bamFiles=${INFILE} \
    --refFile=${reference_genome_seq} \
    --output=${OUTFILE} \
    --filterDuplicates=0;
# for exome/whole genome data (no duplicate filtering)
else
    python ${ngs_platypus}/Platypus.py callVariants \
    --bamFiles=${INFILE} \
    --refFile=${reference_genome_seq} \
    --output=${OUTFILE};

















