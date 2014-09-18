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
  KNOWN_SNPS_b138=${gatk_resources}/dbsnp_138.b37.vcf
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
  
##---------------------- FASTQ-QC ----------------------## 

  if [ ! -s ${SOUT}/fastq/${FASTQ1} ] && [ ! -s ${SOUT}/fastq/${FASTQ2} ]
then
  echo " NGSeasy: Copying fastq files from ${FASTQDIR}/ to ${SOUT}/fastq/ " `date`
    cp ${FASTQDIR}/${FASTQ1} ${SOUT}/fastq/${FASTQ1}
    cp ${FASTQDIR}/${FASTQ2} ${SOUT}/fastq/${FASTQ2}

else
  echo " NGSeasy: Fastq Files exist in  ${SOUT}/fastq/ " `date`
  ls ${SOUT}/fastq/
fi

## set new names for copied fastq files
  rawFASTQ1=`basename ${SOUT}/fastq/${FASTQ1} _1.fq.gz`
  rawFASTQ2=`basename ${SOUT}/fastq/${FASTQ2} _2.fq.gz`
    
echo " NGSeasy: Fastq Basename : [$rawFASTQ1] "

# FASTQC on raw files

echo ""
echo "................................................"
echo " NGSeasy: START Pre-Alignment QC " `date`
echo "................................................"
echo ""

# check if qc'd data alread exists 
if [ ! -s ${SOUT}/fastq/${rawFASTQ1}_1.fq_fastqc.zip ] && [ ! -s ${SOUT}/fastq/${rawFASTQ2}_2.fq_fastqc.zip ]
then
  echo " NGSeasy: Run Pre-Alignment QC on raw Fastq files " `date`

  /usr/local/pipeline/FastQC/fastqc \
    --threads ${NCPU} \
    --extract \
    --quiet \
    --dir ${SOUT}/tmp \
    --outdir ${SOUT}/fastq \
    ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz \
    ${SOUT}/fastq/${rawFASTQ2}_2.fq.gz
    
else
  echo " NGSeasy: Pre-Alignment QC on raw Fastq files Already run " `date`
fi

echo ""
echo "................................................"
echo " NGSeasy: START Filtering/Trimming " `date`
echo "................................................"
echo ""

  # Trimmomatic paired output
  qcdPeFASTQ1=${SOUT}/fastq/${rawFASTQ1}_1.filtered.fq.gz;
  qcdPeFASTQ2=${SOUT}/fastq/${rawFASTQ2}_2.filtered.fq.gz;
  # Trimmomatic unpaired ouput
  qcdSeFASTQ1=${SOUT}/fastq/${rawFASTQ1}_1.unpaired.fq.gz;
  qcdSeFASTQ2=${SOUT}/fastq/${rawFASTQ2}_2.unpaired.fq.gz;
 
# if qc'd files exits then skip this step
if [ ! -s ${qcdPeFASTQ1} ] && [ ! -s ${qcdPeFASTQ2} ]
then

  # skip next step if using novoalign or ya really cant be bothered - saves about an hour
  echo " NGSeasy: Running Trimmomatic: A flexible read trimming tool for NGS data " `date`

    # run Trimmomatic
    java -jar /usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar PE \
    -threads ${NCPU} \
    -phred33 \
    -trimlog ${SOUT}/fastq/trimm_qc.log \
    ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz ${SOUT}/fastq/${rawFASTQ2}_2.fq.gz \
    ${qcdPeFASTQ1} ${qcdSeFASTQ1} \
    ${qcdPeFASTQ2} ${qcdSeFASTQ2} \
    ILLUMINACLIP:${adapter_fa}:2:30:10:5:true \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50

  echo " NGSeasy: Run Pre-Alignment QC on Filtered/Trimmed Fastq files " `date`
    # FASTQC on paired trimmed files

    /usr/local/pipeline/FastQC/fastqc \
    --threads ${NCPU} \
    --extract \
    --quiet \
    --dir ${SOUT}/tmp \
    --outdir ${SOUT}/fastq \
    ${qcdPeFASTQ1} ${qcdPeFASTQ2};

else
  echo " NGSeasy: QC'd Fastq files Already Exist " `date`
  echo "................................................"
  zcat ${qcdPeFASTQ1} | head -4
  echo "................................................"
fi

echo ""
echo "................................................"
echo " NGSeasy: END Filtering/Trimming " `date`
echo "................................................"
echo ""

##---------------------- ALIGNMENT ----------------------------------------------------------------------------------------------------------------------------##


echo "................................................"
echo " NGSeasy: Start Alignment [${ALIGNER}] " `date`
echo "................................................"
echo ""

if [ "${ALIGNER}" == "bwa" ] && [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.bam ]
then
  # BWA alignment
  echo " NGSeasy: Running bwa " `date`
    /usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;

 echo " NGSeasy: SAM to BAM and INDEX  " `date` 
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam > ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   
elif [ "${ALIGNER}" == "bowtie2" ] && [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.bam ]
then
  echo " NGSeasy: Running bowtie2 " `date`
    # Bowtie2 alignment
   /usr/local/pipeline/bowtie2-2.2.3/bowtie2 --local --threads ${NCPU} -x ${REFGenomes}/human_g1k_v37 -1 ${qcdPeFASTQ1} -2 ${qcdPeFASTQ2} -S ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;

  echo " NGSeasy: SAM to BAM and INDEX " `date`
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam > ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

elif [ "${ALIGNER}" == "novoalign" ] && [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.bam ]
then
  echo " NGSeasy: Running novolaign " `date`
    # Novoalign alignment 
   /usr/local/pipeline/novocraft/novoalign \
   -d ${REFGenomes}/human_g1k_v37.novoIndex \
   -f ${qcdPeFASTQ1} ${qcdPeFASTQ2} \
   -F STDFQ \
   --Q2Off \
   --3Prime \
   -g 40 -x 6 \
   -r All \
   -i PE 300,150 -c ${NCPU} -k -K ${SOUT}/alignments/${BAM_PREFIX}.K.stats -o SAM > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;

  echo " NGSeasy: SAM to BAM and INDEX  " `date`
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam > ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

elif [ "${ALIGNER}" == "stampy" ] && [ ! -s ${SOUT}/alignments/${BAM_PREFIX}.bam ]
then
  echo " NGSeasy: Running stampy " `date`

  echo " NGSeasy: Running stampy bwa "
    /usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} \
    > ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam;

  echo " NGSeasy: Running sam to bam stampy bwa "
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam > ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam;
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
   /usr/local/pipeline/samtools/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam > ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
   /usr/local/pipeline/samtools/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
   /usr/local/pipeline/samtools/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

  echo " NGSeasy: Remove Stampy intermediate bam files " 
    rm ${SOUT}/alignments/${BAM_PREFIX}.raw.sam
    rm ${SOUT}/alignments/${BAM_PREFIX}.raw.bam
    rm ${SOUT}/alignments/${BAM_PREFIX}.tmp
    rm ${SOUT}/alignments/${BAM_PREFIX}.tmpsort.bam


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

if [ ! -e ${SOUT}/reports/${BAM_PREFIX}.dupemk_metrics ]
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
  METRICS_FILE=${SOUT}/reports/${BAM_PREFIX}.dupemk_metrics;
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

if [ ! -s ${SOUT}/reports/${BAM_PREFIX}.merged.bed ]
then
echo " NGSeasy: Converting Aligned BAM To BED File " `date`
# pulls out paired only reads with min qual 10. Set low as donwstream tools filter these regions
 /usr/local/pipeline/samtools/samtools view -b -h -q 10 -F 1796 ${SOUT}/alignments/${BAM_PREFIX}.bam | /usr/local/pipeline/bedtools2/bin/bedtools bamtobed -i stdin > ${SOUT}/reports/${BAM_PREFIX}.bed;

echo " NGSeasy: Converting Aligned BED To MERGED BED File " `date` 
 /usr/local/pipeline/bedtools2/bin/bedtools merge -i ${SOUT}/reports/${BAM_PREFIX}.bed > ${SOUT}/reports/${BAM_PREFIX}.merged.bed;
fi

echo ""
echo "................................................"
echo " NGSeasy: END Alignment [${ALIGNER}] " `date`
echo "................................................"
echo ""

##---------------------- POST ALIGNMENT QC ----------------------##

echo "................................................"
echo " NGSeasy: START Post Alignmnet QC " `date`
echo "................................................"
echo ""

if [ ! -s ${SOUT}/reports/${BAM_PREFIX}.alignment_summary_metrics ]
then 

  echo " NGSeasy: Running Post Alignmnet QC " `date`
# CollectMultipleMetrics
  echo " NGSeasy: CollectMultipleMetrics " `date`
 
 java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectMultipleMetrics.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
  OUTPUT=${SOUT}/reports/${BAM_PREFIX} \
  REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
  PROGRAM=CollectAlignmentSummaryMetrics \
  PROGRAM=CollectInsertSizeMetrics \
  PROGRAM=QualityScoreDistribution \
  PROGRAM=MeanQualityByCycle;

# CollectAlignmentSummaryMetrics
  echo " NGSeasy: CollectAlignmentSummaryMetrics " `date`
  
  java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectAlignmentSummaryMetrics.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
  OUTPUT=${SOUT}/reports/${BAM_PREFIX}.alignment_summary_metrics_alt \
  REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
  ASSUME_SORTED=true \
  METRIC_ACCUMULATION_LEVEL=SAMPLE;

# CollectWgsMetrics
  echo " NGSeasy: CollectWgsMetrics " `date`
  
  java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectWgsMetrics.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
  OUTPUT=${SOUT}/reports/${BAM_PREFIX}.wgs_coverage \
  REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
  MINIMUM_MAPPING_QUALITY=20 \
  MINIMUM_BASE_QUALITY=20 \
  COVERAGE_CAP=1000;

  awk 'NR>9' ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage > ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage_hist
  sed '8q' ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage > ${SOUT}/reports/${BAM_PREFIX}.wgs_coverage_stats

# FlagStat
  echo " NGSeasy: FlagStat " `date`
  
    java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FlagStat -R ${REFGenomes}/human_g1k_v37.fasta \
    -I ${SOUT}/alignments/${BAM_PREFIX}.bam  \
    -o ${SOUT}/reports/${BAM_PREFIX}.FlagStat;
    
else
    echo " NGSeasy: Post Alignmnet QC Already Run " `date` 
fi

echo ""
echo "................................................"
echo " NGSeasy: END Post Alignmnet QC " `date`
echo "................................................"
echo ""

## Depending in NGS type Run PicardTools CollectTargetedPcrMetrics
## NB This is rough and ready quicl fix to get coverage ver genome bins, annotated exomes and targetd sequene files

echo "................................................"
echo " NGSeasy: START Post Alignmnet Coverage Calculations " `date`
echo "................................................"
echo ""

if [ ! -s ${SOUT}/reports/${BAM_PREFIX}.genomecov.bed ]
then
echo " NGSeasy: Running bedtools genomecov [-bga] " `date`

  /usr/local/pipeline/bedtools2/bin/bedtools genomecov -ibam ${SOUT}/alignments/${BAM_PREFIX}.bam -bga > ${SOUT}/reports/${BAM_PREFIX}.genomecov.bed;

else
  echo " NGSeasy: bedtools genomecov " `date`
fi

# Picardt tools CollectTargetedPcrMetric
echo " NGSeasy: Run CollectTargetedPcrMetric " `date`
echo " NGSeasy: This requires a Custom or Exome BED File provided by the user or manufacturer of your NGS Exome Library " 

if [ "${NGS_TYPE}" == "TGS" ] && [ -s ${BED_ANNO} ] && [ ! -s ${SOUT}/reports/${BAM_PREFIX}.target_coverage ]
then
  echo " NGSeasy: Finding TGS coverage over target regions in ${BED_ANNO} " `date`
  
    java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectTargetedPcrMetrics.jar \
    TMP_DIR=${SOUT}/tmp \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=100000 \
    INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
    OUTPUT=${SOUT}/reports/${BAM_PREFIX}.target_coverage \
    REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
    AMPLICON_INTERVALS=${BED_ANNO} \
    TARGET_INTERVALS=${BED_ANNO} \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    PER_TARGET_COVERAGE=${SOUT}/reports/${BAM_PREFIX}.per_target_coverage;
    
elif [ "${NGS_TYPE}" == "WEX" ]  && [ -s ${BED_ANNO} ] && [ ! -s ${SOUT}/reports/${BAM_PREFIX}.target_coverage ]
then
  echo " NGSeasy: Finding WEX Coverage over annotated exons/genes from ensembl " `date`
    java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectTargetedPcrMetrics.jar \
    TMP_DIR=${SOUT}/tmp \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=100000 \
    INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
    OUTPUT=${SOUT}/reports/${BAM_PREFIX}.target_coverage \
    REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
    AMPLICON_INTERVALS=${REFGenomes}/${BED_ANNO} \
    TARGET_INTERVALS=${REFGenomes}/${BED_ANNO} \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    PER_TARGET_COVERAGE=${SOUT}/reports/${BAM_PREFIX}.per_target_coverage;
    
elif [ "${NGS_TYPE}" == "WGS" ] && [ ! -s ${SOUT}/reports/${BAM_PREFIX}.target_coverage ]
then
  echo " NGSeasy: Finding WGS Coverage over 500bp windows " `date`
    java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectTargetedPcrMetrics.jar \
    TMP_DIR=${SOUT}/tmp \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=100000 \
    INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
    OUTPUT=${SOUT}/reports/${BAM_PREFIX}.target_coverage \
    REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
    AMPLICON_INTERVALS=${REFGenomes}/human_g1k_v37_0.5Kwindows.bed \
    TARGET_INTERVALS=${REFGenomes}/human_g1k_v37_0.5Kwindows.bed \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    PER_TARGET_COVERAGE=${SOUT}/reports/${BAM_PREFIX}.per_target_coverage;

else

  echo " NGSeasy: Whoops! Something Went wrong! NGS_TYPE and Annotation Files not found! Check your config file and data...or not..."
  echo " NGSeasy: Skipping Collect Targeted Pcr Metrics...."
  echo " NGSeasy: You may want to run this manually later....it is quite nice"
  echo ""
fi

echo "................................................"
echo " NGSeasy: END Post Alignmnet Coverage Calculations " `date`
echo "................................................"
echo ""


##---------------------- VARIANT CALLING SINGLE SAMPLE ----------------------------------------------------------------##
echo "................................................"
echo " NGSeasy: START SNP and Small INDEL Calling " `date`
echo "................................................"
echo ""

if [ "${VARCALLER}" == "freebayes" ]
then

echo " NGSeasy: Starting Variant Calling using Freebayes " `date`
  /usr/local/pipeline/freebayes/bin/freebayes \
    -f ${REFGenomes}/human_g1k_v37.fasta \
    -b ${SOUT}/alignments/${BAM_PREFIX}.bam \
    --min-coverage 4 \
    --min-mapping-quality 30 \
    --min-base-quality 20 \
    --genotype-qualities > ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ;
    
  # copy vcf to cohort vcf directory
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs/;
  
elif [ "${VARCALLER}" == "platypus" ]
then

  echo " NGSeasy: Starting Variant Calling using Platypus " `date`
    if [ ${NGS_TYPE} == "TGS" ]
    then
    echo " NGSeasy: NGS_TYPE is Targeted so no duplicate filtering  " `date`
    # for exome/whole genome data no duplicate filtering
      python /usr/local/pipeline/Platypus_0.7.9.1/Platypus.py callVariants \
      --nCPU ${NCPU} \
      --bamFiles=${SOUT}/alignments/${BAM_PREFIX}.bam \
      --refFile=${REFGenomes}/human_g1k_v37.fasta \
      --output=${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
      --filterDuplicates=0;
      
  # copy vcf to cohort vcf directory
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs/;
      
     else
	python /usr/local/pipeline/Platypus_0.7.9.1/Platypus.py callVariants \
	  --nCPU ${NCPU} \
	  --bamFiles=${SOUT}/alignments/${BAM_PREFIX}.bam \
	  --refFile=${REFGenomes}/human_g1k_v37.fasta \
	  --output=${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf;
	  
	    # copy vcf to cohort vcf directory
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs/;
  
    fi
	  
elif [ ${VARCALLER} == "gatk_ug" ]
then
  # UnifiedGenotyper EMIT_ALL_CONFIDENT_SITES
  echo " NGSeasy: Running GATK UnifiedGenotyper " `date`
  java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
  -I ${SOUT}/alignments/${BAM_PREFIX}.bam \
  -o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
  -stand_call_conf 30 \
  -stand_emit_conf 10 \
  --dbsnp ${KNOWN_SNPS_b138} \
  -dcov 250 -minPruning 4 \
  --unsafe ALL \
  --genotype_likelihoods_model BOTH \
  --genotyping_mode DISCOVERY \
  --output_mode ${GTMODEGATK} \
  --annotation AlleleBalance \
  --annotation BaseCounts \
  --annotation BaseQualityRankSumTest \
  --annotation ChromosomeCounts \
  --annotation ClippingRankSumTest \
  --annotation Coverage \
  --annotation FisherStrand \
  --annotation GCContent \
  --annotation HaplotypeScore \
  --annotation HomopolymerRun \
  --annotation InbreedingCoeff \
  --annotation LikelihoodRankSumTest \
  --annotation LowMQ \
  --annotation MappingQualityRankSumTest \
  --annotation MappingQualityZero \
  --annotation QualByDepth \
  --annotation RMSMappingQuality \
  --annotation ReadPosRankSumTest \
  --annotation SpanningDeletions \
  --annotation TandemRepeatAnnotator \
  --annotation VariantType;
  
  # copy vcf to cohort vcf directory
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs/;

elif [ "${VARCALLER}" == "gatk_hc" ]
then 
  echo " NGSeasy: Running GATK HaplotypeCaller THIS TAKES A LOOOONG TIME " `date`
  ## HaplotypeCaller Standard EMIT_ALL_CONFIDENT_SITES EMIT_VARIANTS_ONLY
  java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
  -I ${SOUT}/alignments/${BAM_PREFIX}.bam \
  -o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
  -stand_call_conf 30 \
  -stand_emit_conf 10 \
  --dbsnp ${KNOWN_SNPS_b138} \
  -dcov 250 -minPruning 4 \
  --unsafe ALL \
  -pairHMM VECTOR_LOGLESS_CACHING \
  --genotyping_mode DISCOVERY \
  --output_mode ${GTMODEGATK} \
  --annotation AlleleBalance \
  --annotation BaseCounts \
  --annotation BaseQualityRankSumTest \
  --annotation ChromosomeCounts \
  --annotation ClippingRankSumTest \
  --annotation Coverage \
  --annotation FisherStrand \
  --annotation GCContent \
  --annotation HaplotypeScore \
  --annotation HomopolymerRun \
  --annotation InbreedingCoeff \
  --annotation LikelihoodRankSumTest \
  --annotation LowMQ \
  --annotation MappingQualityRankSumTest \
  --annotation MappingQualityZero \
  --annotation QualByDepth \
  --annotation RMSMappingQuality \
  --annotation ReadPosRankSumTest \
  --annotation SpanningDeletions \
  --annotation TandemRepeatAnnotator \
  --annotation VariantType;
  
  # copy vcf to cohort vcf directory
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs/;
  
elif [ ${VARCALLER} == "gatk_hc_gvcf" ]
then
  echo " NGSeasy: Running GATK HaplotypeCaller GVCF THIS TAKES A VERY LOOOONG TIME" `date` 
  ## HaplotypeCaller GVCF
  java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
  -I ${SOUT}/alignments/${BAM_PREFIX}.bam \
  -o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.g.vcf \
  -stand_call_conf 30 \
  -stand_emit_conf 10 \
  --dbsnp ${KNOWN_SNPS_b138} \
  -dcov 250 -minPruning 4 \
  --unsafe ALL \
  -pairHMM VECTOR_LOGLESS_CACHING \
  --emitRefConfidence GVCF \
  --variant_index_type LINEAR \
  --variant_index_parameter 128000 \
  --annotation AlleleBalance \
  --annotation BaseCounts \
  --annotation BaseQualityRankSumTest \
  --annotation ChromosomeCounts \
  --annotation ClippingRankSumTest \
  --annotation Coverage \
  --annotation FisherStrand \
  --annotation GCContent \
  --annotation HaplotypeScore \
  --annotation HomopolymerRun \
  --annotation InbreedingCoeff \
  --annotation LikelihoodRankSumTest \
  --annotation LowMQ \
  --annotation MappingQualityRankSumTest \
  --annotation MappingQualityZero \
  --annotation QualByDepth \
  --annotation RMSMappingQuality \
  --annotation ReadPosRankSumTest \
  --annotation SpanningDeletions \
  --annotation TandemRepeatAnnotator \
  --annotation VariantType;
  ## -minPruning 10 -dcov 250

  # copy vcf to cohort vcf directory
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs/;

else 
  echo " NGSeasy: Something Went wrong! Check your config file "
  exit 1 
fi
## TO ADD:- LCR filters and HIGH COV LOW MAP Q FILTERS AND ANNOTATIONS

##---------------------- CLEANUP ----------------------------------------------------------------##

if [ "${CLEANUP}" == "TRUE" ]
then
echo " NGSeasy: Removing Intermediate SAM/BAM Files " `date`

  rm -v ${SOUT}/alignments/${BAM_PREFIX}.realn.bam
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.dupemk.bam
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.addrg.bam
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.sort.bam
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.raw.sam
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.raw.bam
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.realn.bam.bai
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.dupemk.bam.bai
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.addrg.bam.bai
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.sort.bam.bai
  
else
  echo " NGSeasy: Keepping all Intermediate SAM/BAM Files " `date`
fi
 

echo ""
echo "................................................"
echo " NGSeasy: END SNP and Small INDEL Calling " `date`
echo "................................................"
echo ""
echo ""
#-----------------------------------------THE END------------------------------------------------------------#
echo ""
echo ""
echo "................................................"
echo " NGSeasy: END NGSeasy Pipleine...Have a Nice Day " `date`
echo "................................................"
echo ""














