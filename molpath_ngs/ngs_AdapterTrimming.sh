#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

#########################################
## adapter trimming with trimmomatic   ##
#########################################

# TARGET and SOURCE
sample_name=${1}
sample_dir=${2}
# fastq file(s) with path

if [[ $# -eq 4 ]];then
    sample_librarytype='PE'
    IFS=$' ' SORTED=($(sort <<< "${3} ${4}"))
    # inputs
    FQ1=${SORTED[0]}
    FQ2=${SORTED[1]}
    # outputs
    OUT1pe=${sample_dir}/`basename $FQ1 .fastq`.filtered.fastq
    OUT1se=${sample_dir}/`basename $FQ1 .fastq`.unpaired.fastq
    OUT2pe=${sample_dir}/`basename $FQ2 .fastq`.filtered.fastq
    OUT2se=${sample_dir}/`basename $FQ2 .fastq`.unpaired.fastq
elif [[ $# -eq 3 ]];then
    sample_librarytype='SE'
    # input
    FQ1=${fastq_dir}/${3}
    # output
    OUT1se=${sample_dir}/`basename $FQ1 .fastq`.filtered.fastq
else
    # is something else
    echo "argument list error"
    exit 1
fi

echo "----------------------------------------------------------------------------------------"
echo " in FASTQ    " $FQ1 "<>" $FQ2
echo " Sample Name " $sample_name # target directory name (sample)
echo " Sample Dir  " $sample_dir  # target directory
echo " Reads are   " ${sample_librarytype}
echo " out FASTQ   " ${OUT1pe} ${OUT1se} ${OUT2pe} ${OUT2se}
echo "----------------------------------------------------------------------------------------"

echo "PreFiltering QC"

fastqc --noextract --outdir=${sample_dir} ${FQ1};
if [[${sample_librarytype} -eq 'PE']];then
    fastqc --noextract --outdir=${sample_dir} ${FQ2};
fi

echo "Start Filtering/Trimming"

cd ${sample_dir}

${java_v1_7}/java  -XX:ParallelGCThreads=4 -Xmx${java_mem}g -jar ${ngs_trimmomatic} \
${sample_librarytype} \
-phred64 \
-trimlog ${sample_name}.trimming.log \
-threads 4 \
${FQ1} ${FQ2} \
${OUT1pe} ${OUT1se} \
${OUT2pe} ${OUT2se} \
ILLUMINACLIP:${ngs_trimmomatic}/adapters/TruSeq3-${sample_librarytype}.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36

### CHANGE THIS TO READ FROM PIPELINE CONFIGURATION
#${ngs_filtering}

echo "PostFiltering QC"

fastqc --noextract --outdir=${sample_dir} ${OUT1pe};
if [[${sample_librarytype} -eq 'PE']];then
    fastqc --noextract --outdir=${sample_dir} ${OUT1se};
    fastqc --noextract --outdir=${sample_dir} ${OUT2pe};
    fastqc --noextract --outdir=${sample_dir} ${OUT2se};
fi



