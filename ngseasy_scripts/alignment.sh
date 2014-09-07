#!/bin/bash

# Stephen Newhouse <stephen.j.newhouse@gmail.com>
# UNDER HEAVY DEVELOPMENT
echo ""
echo ""
echo  "#------------------------------------------------------------#"
echo  " Starting NGSeasy : NGS made easy!!!!"
echo  " Version 1.0"
echo  " Authors: Amos Folarin <amosfolarin@gmail.com>"
echo  " Authors: Stephen Newhouse <stephen.j.newhouse@gmail.com>"
echo  " Run Date : `date +"%d-%m-%y"`"
echo  "#------------------------------------------------------------#"
echo ""

# refgenomes and contaminant lists
  REFGenomes='/home/pipeman/reference_genomes_b37'
  adapter_fa='/home/pipeman/reference_genomes_b37/contaminant_list.fa'
  gatk_resources='/home/pipeman/gatk_resources'
  FASTQDIR='/home/pipeman/fastq_raw'
  KNOWN_SNPS_b138=/usr/local/pipeline/gatk_resources/dbsnp_138.b37.vcf
  echo " NGSeasy: Note:- Using ${KNOWN_SNPS_b138} for SNP annotations in GATK.\n Edit bash script if you want to change this " `date`
  
#---------------------------------------------CONFIG------------------------------------------------#
# get the following from the config file

# set varibales
  DATE=`date +"%d%m%y"`
  POJECT_ID=${1}
  SAMPLE_ID=${2}
  FASTQ1=${3}
  FASTQ2=${4}
  PROJECT_DIR=${5}
  DNA_PREP_LIBRARY_ID=${6}
  NGS_PLATFORM=${7}
  NGS_TYPE=${8}
  BED_ANNO=${9}
  PIPELINE=${10}
  ALIGNER=${11}
  VARCALLER=${12}
  GTMODEGATK=${13}
  CLEANUP=${14}
  NCPU=${15}

  
echo ""
echo "................................................"
echo " NGSeasy: START NGSeasy Pipleine [${PIPELINE}]" `date`
echo "................................................"
echo ""

##############################
## BEGING FULL NGS PIPELINE ##
##############################

## BAM PREFIX 
  BAM_PREFIX=${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}
  echo " NGSeasy: Setting BAM_PREFIX directory [$BAM_PREFIX]"
  echo ""
 
## OUTPUT SAMPLE DIR
 SOUT=${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
  
if [ ! -e ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID} ]
then
  echo " NGSeasy: Cant Find Project directory. This is then end. Please Stop and check everything is ok " `date`
  exit 1

else 
  echo " NGSeasy: Setting OUTPUT directory [${SOUT}]"
fi

##---------------------- ALIGNMENT ----------------------------------------------------------------------------------------------------------------------------##


echo "................................................"
echo " NGSeasy: Start Alignment [${ALIGNER}] " `date`
echo "................................................"
echo ""

if [ "${ALIGNER}" == "bwa" ] && [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.sort.bam ]
then
  # BWA alignment
  echo " NGSeasy: Running bwa " `date`
    /usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
  echo " NGSeasy: SAM to BAM and INDEX  " `date`  
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   
elif [ "${ALIGNER}" == "bowtie2" ] && [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.sort.bam ]
then
  echo " NGSeasy: Running bowtie2 " `date`
    # Bowtie2 alignment
   /usr/local/pipeline/bowtie2-2.2.3/bowtie2 --local --threads ${NCPU} -x ${REFGenomes}/human_g1k_v37 -1 ${qcdPeFASTQ1} -2 ${qcdPeFASTQ2} -S ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
  echo " NGSeasy: SAM to BAM and INDEX " `date`
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam > ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

elif [ "${ALIGNER}" == "novoalign" ] && [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.sort.bam ]
then
  echo " NGSeasy: Running novolaign " `date`
    # Novoalign alignment 
    # TO DO: Use raw${SOUT}/fastq/${FASTQ1} ${SOUT}/fastq/${FASTQ2} as novoalign does this all internally if provideing adapter lists
    # ADD -a @ADAPTERS
   /usr/local/pipeline/novocraft/novoalign -d ${REFGenomes}/human_g1k_v37.fasta -f ${qcdPeFASTQ1} ${qcdPeFASTQ2} -F STDFQ --Q2Off --3Prime -g 40 -x 6 -r All -i PE 300,150 -c ${NCPU} -k -K ${SOUT}/alignments/${BAM_PREFIX}.K.stats -o SAM > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
  echo " NGSeasy: SAM to BAM and INDEX  " `date`
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

elif [ "${ALIGNER}" == "stampy" ] && [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.sort.bam ]
then
  echo " NGSeasy: Running stampy " `date`

  echo " NGSeasy: Running stampy bwa "
    /usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} \
    > ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam;

  echo " NGSeasy: Running sam to bam stampy bwa "
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam;
    rm ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam

  echo " NGSeasy: Running sort bam stampy bwa "
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam ${SOUT}/alignments/${BAM_PREFIX}.tmpsort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.tmpsort.bam;
   rm ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam

  echo " NGSeasy: Running Stampy aligner on bwa bam "
  
python  /usr/local/pipeline/stampy-1.0.23/stampy.py \
    -g ${REFGenomes}/human_g1k_v37 \
    -h ${REFGenomes}/human_g1k_v37 \
    -t ${NCPU} \
    --bamsortprefix ${SOUT}/tmp \
    --bamkeepgoodreads \
    --sanger \
    --gapopen=40 \
    --gapextend=3 \
    --noautosense \
    --insertsize=300 \
    --insertsd=50 \
    -M ${SOUT}/alignments/${BAM_PREFIX}.tmpsort.bam \
    -o ${SOUT}/alignments/${BAM_PREFIX}.raw.sam \
    -f sam

  echo " NGSeasy: Running Stampy sam to bam "
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

  echo " NGSeasy: Remove Stampy intermediate bam files " 
    rm ${SOUT}/alignments/${BAM_PREFIX}.raw.sam
    rm ${SOUT}/alignments/${BAM_PREFIX}.raw.bam
    rm ${SOUT}/alignments/${BAM_PREFIX}.tmp

fi

echo " NGSeasy: Basic $ALIGNER Complete " `date`

##---------------------- RAW ALIGNMENT PROCESSING ---------------------------------------------------------------------------------------##

# AddOrReplaceReadGroups

echo " NGSeasy: START AddOrReplaceReadGroups " `date`

if [ ! -e ${SOUT}/alignments/${BAM_PREFIX}.addrg.bam ]
then
echo " NGSeasy: Adding Read Group Information " `date`
echo " NGSeasy: Getting Platform Unit Information "  `date`
  platform_unit=`zcat ${qcdPeFASTQ1} | head -1 | perl -p -i -e 's/:/\t/' | cut -f 1 | perl -p -i -e 's/@//g'`
echo " NGSeasy: Platform Unit: [$platform_unit]"  `date`

java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/AddOrReplaceReadGroups.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  CREATE_INDEX=true \
  SORT_ORDER=coordinate \
  RGID=${BAM_PREFIX} \
  RGLB=${DNA_PREP_LIBRARY_ID} \
  RGPL=${NGS_PLATFORM} \
  RGPU=${platform_unit} \
  RGSM=${SAMPLE_ID} \
  RGDT=${DATE} \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.sort.bam \
  OUTPUT=${SOUT}/alignments/${BAM_PREFIX}.addrg.bam;
fi

# MarkDuplicates
echo " NGSeasy: START MarkDuplicates " `date`  

if [ ! -e ${SOUT}/alignments/${BAM_PREFIX}.dupemk_metrics ]
then
echo " NGSeasy: Marking Duplicate Reads " `date`

  java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/MarkDuplicates.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  CREATE_INDEX=true \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.addrg.bam \
  OUTPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
  METRICS_FILE=${SOUT}/alignments/${BAM_PREFIX}.dupemk_metrics;
fi

  # FindCoveredIntervals: these are used in GATK Calling to help speed things up

echo " NGSeasy: START FindCoveredIntervals " `date`  
if [ ! -s ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list ]
then
echo " NGSeasy: Finding Covered Intervals : minimum coverage of 4 " `date`
  java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FindCoveredIntervals -R ${REFGenomes}/human_g1k_v37.fasta \
  -I ${SOUT}/alignments/${BAM_PREFIX}.bam  \
  -o ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list \
  --coverage_threshold 4;
fi
  
  
# BAM to BED
echo " NGSeasy: START bamtobed " `date`  

if [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.merged.bed ]
then
echo " NGSeasy: Converting Aligned BAM To BED File " `date`
# pulls out paired only reads with min qual 10. Set low as donwstream tools filter these regions
 /usr/local/pipeline/samtools/samtools view -b -h -q 10 -F 1796 ${SOUT}/alignments/${BAM_PREFIX}.bam | /usr/local/pipeline/bedtools2/bin/bedtools bamtobed -i stdin > ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.bed;

echo " NGSeasy: Converting Aligned BED To MERGED BED File " `date` 
 /usr/local/pipeline/bedtools2/bin/bedtools merge -i ${SOUT}/alignments/${BAM_PREFIX}.bed > ${SOUT}/alignments/${BAM_PREFIX}.merged.bed;
fi

echo ""
echo "................................................"
echo " NGSeasy: END Alignment [${ALIGNER}] " `date`
echo "................................................"
echo ""
