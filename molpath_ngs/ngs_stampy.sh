#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V


###TODO
# check insert size
# is the proper paired flag used by GATK/platypus





###################################
### Stampy mapping
###################################
#
# ONE SCRIPT FOR ALL DATA TYPES
#
# PE/SE (autodetect)
# illumina only (for now)
# no trimming (done by trimmomatic)
#

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
OUTFILE=${sample_name}.stampy.sam
BWAFILE=${sample_name}.bwa.bam

echo "----------------------------------------------------------------------------------------"
echo " sampleDir        " $sample_dir
echo " CWD              " `pwd`
echo " in FASTQ         " $FQ1 "<>" $FQ2
echo " Reads are        " ${sample_librarytype}
echo " alignment CPUs   " ${novo_cpu}
echo " Reference Genome " ${reference_genome_novoindex}
echo " BWA output BAM   " ${BWAFILE}
echo " FINAL output SAM " ${OUTFILE}
echo "----------------------------------------------------------------------------------------"

if [[ -f ${ngs_bwa}/bwa && ${reference_genome_bwaindex}.bwt && ${reference_genome_bwaindex}.rbwt ]]; then
    echo "BWA mapping"
    ${ngs_bwa}/bwa aln -q10 -t8 ${reference_genome_bwaindex} ${FQ1} > ${FQ1}.sai && ${ngs_bwa}/bwa aln -q10 -t8 hg18 ${FQ2} > ${FQ2}.sai;
    echo "BWA merging mappings"
    ${ngs_bwa}/bwa sampe hg18 1.sai 2.sai ${FQ1} ${FQ2} | \
    ${ngs_samtools}/samtools view -Sb - > ${BWAFILE} && rm -f ${FQ1}.sai ${FQ2}.sai
    echo "STAMPY mapping"
    ${ngs_stampy}/stampy.py \
    -g ${reference_genome_stampyindex} \
    -h ${reference_genome_stampyindex} \
    -M ${BWAFILE} \
    #-t ${novo_cpu} \
    --bamkeepgoodreads \
    --sanger \
    --gapopen=40 \
    --gapextend=3 \
    --noautosense \
    --insertsize=200 \
    --insertsd=50 \
    > ${OUTFILE};
else
    echo "BWA not available - STAMPY only mapping (slower)"
    ${ngs_stampy}/stampy.py \
    -g ${reference_genome_stampyindex} \
    -h ${reference_genome_stampyindex} \
    -M ${FQ1} ${FQ2} \
    #-t ${novo_cpu} \
    --sanger \
    --gapopen=40 \
    --gapextend=3 \
    --noautosense \
    --insertsize=200 \
    --insertsd=50 \
    > ${OUTFILE};
fi

















