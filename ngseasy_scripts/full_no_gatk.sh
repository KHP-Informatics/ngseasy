#!/bin/bash

## FOR DEBUGGING
# set -x
##################################################################################

echo ""
echo "................................................"
echo " NGSeasy: START NGSeasy Pipleine " `date`
echo "................................................"
echo ""

# refgenomes and contaminant lists
  REFGenomes='/home/pipeman/reference_genomes_b37'
  adapter_fa='/home/pipeman/reference_genomes_b37/contaminant_list.fa'
  gatk_resources='/home/pipeman/gatk_resources'
  FASTQDIR='/home/pipeman/fastq_raw'

#---------------------------------------------CONFIG------------------------------------------------#
# get the following from the config file

# set varibales
  DATE=`date +"%d%m%y"`
  POJECT_ID=$1
  SAMPLE_ID=$2
  FASTQ1=$3
  FASTQ2=$4
  PROJECT_DIR=$5
  DNA_PREP_LIBRARY_ID=$6
  NGS_PLATFORM=$7
  NGS_TYPE=$8
  BED_ANNO=$9
  PIPELINE=$10
  ALIGNER=$11
  VARCALLER=$12
  GTMODEGATK=$13
  CLEANUP=$14
  NCPU=$15


##############################
## BEGING FULL NGS PIPELINE ##
##############################

## OUTPUT SAMPLE DIR
  SOUT=${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
  echo " NGSeasy: Setting OUTPUT directory [${SOUT}]"

## BAM PREFIX 
  BAM_PREFIX=${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}
  echo " NGSeasy: Setting BAM_PREFIX directory [$BAM_PREFIX]"
  echo ""
##---------------------- FASTQ-QC ----------------------## 

## symlink orginal fastq file for sample x to sample x fastq directory
echo " NGSeasy: Copying fastq files from ${FASTQDIR}/ to ${SOUT}/fastq/ " `date`

  cp -v ${FASTQDIR}/${FASTQ1} ${SOUT}/fastq/${FASTQ1}
  cp -v ${FASTQDIR}/${FASTQ2} ${SOUT}/fastq/${FASTQ2}

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

echo " NGSeasy: Run Pre-Alignment QC on raw Fastq files " `date`
  /usr/local/pipeline/FastQC/fastqc \
  --threads ${NCPU} \
  --extract \
  --quiet \
  --dir ${SOUT}/tmp \
  --outdir ${SOUT}/fastq \
  ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz \
  ${SOUT}/fastq/${rawFASTQ2}_2.fq.gz

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
 

# skip next step if using novoalign or ya really cant be bothered - saves about an hour
if [ "$ALIGNER" != "novoalign" ]
then

echo " NGSeasy: Run Trimmomatic " `date`
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
echo " NGSeasy: Skipping Filtering/Trimming as Aligner is [novoalign] " `date`
  cp ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz ${qcdPeFASTQ1}
  cp ${SOUT}/fastq/${rawFASTQ1}_2.fq.gz ${qcdPeFASTQ2}
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

if [ "${ALIGNER}" == "bwa" ]
then
  # BWA alignment
  echo " NGSeasy: Running bwa " `date`
    /usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} \
    > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
    /usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
    
elif [ "${ALIGNER}" == "bowtie2" ]
then
  echo " NGSeasy: Running bowtie2 " `date`
    # Bowtie2 alignment
    /usr/local/pipeline/bowtie2-2.2.3/bowtie2 --local -x ${REFGenomes}/human_g1k_v37 -1 ${qcdPeFASTQ1} -2 ${qcdPeFASTQ2} -S ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
    /usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam > ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

elif [ "${ALIGNER}" == "novoalign" ]
then
  echo " NGSeasy: Running novolaign " `date`
    # Novoalign alignment 
    # TO DO: Use raw${SOUT}/fastq/${FASTQ1} ${SOUT}/fastq/${FASTQ2} as novoalign does this all internally if provideing adapter lists
    # ADD -a @ADAPTERS
    /usr/local/pipeline/novocraft -d ${REFGenomes}/human_g1k_v37.fasta -f ${qcdPeFASTQ1} ${qcdPeFASTQ2} -F STDFQ --Q2Off --3Prime -g 40 -x 6 -r All -i PE 300,150 -c ${NCPU} -k -K ${SOUT}/alignments/${BAM_PREFIX}.K.stats -o SAM > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
    /usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

elif [ "${ALIGNER}" == "stampy" ]
then
  echo " NGSeasy: Running stampy " `date`
  # stampy alignment
  echo " NGSeasy: Running stampy bwa "
    /usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} \
    > ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam;

  echo " NGSeasy: Running sam to bam stampy bwa "
    /usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam;
    rm ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam

  echo " NGSeasy: Running sort bam stampy bwa "
    /usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam ${SOUT}/alignments/${BAM_PREFIX}.tmpsort.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.tmpsort.bam;
    rm ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam

  echo " NGSeasy: Running Stampy aligner on bwa bam"
    /usr/local/pipeline/stampy-1.0.23/stampy.py \
    -g human_g1k_v37 \
    -h human_g1k_v37 \
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
    /usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
    /usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

  echo " NGSeasy: Remove Stampy intermediate bam files " 
    rm ${SOUT}/alignments/${BAM_PREFIX}.raw.sam
    rm ${SOUT}/alignments/${BAM_PREFIX}.raw.bam
    rm ${SOUT}/alignments/${BAM_PREFIX}.tmp
    
  echo " NGSeasy: $ALIGNER Complete " ``date``

else
  echo " NGSeasy: Something Went wrong! Aligner not found! Check your config file "
  exit 1
fi

##---------------------- RAW ALIGNMENT PROCESSING ---------------------------------------------------------------------------------------##

echo " NGSeasy: Getting Platform Unit Information " 
  
  platform_unit=`zcat ${qcdPeFASTQ1} | head -1 | perl -p -i -e 's/:/\t/' | cut -f 1 | perl -p -i -e 's/@//g'`

echo " NGSeasy: Getting Platform Unit Information = [$platform_unit]" 
  
echo " NGSeasy: Adding Read Group Information " `date`

  # AddOrReplaceReadGroups
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

echo " NGSeasy: Marking Duplicate Reads " `date`
  # MarkDuplicates
  java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/MarkDuplicates.jar \
  TMP_DIR=${SOUT}/tmp \
  VALIDATION_STRINGENCY=SILENT \
  MAX_RECORDS_IN_RAM=100000 \
  CREATE_INDEX=true \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  INPUT=${SOUT}/alignments/${BAM_PREFIX}.addrg.bam \
  OUTPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
  METRICS_FILE=${SOUT}/alignments/${BAM_PREFIX}.dupemk;

echo " NGSeasy: Removing Intermediate SAM/BAM Files " `date`

  rm -v ${SOUT}/alignments/${BAM_PREFIX}.addrg.ba*
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.sort.ba*
  rm -v ${SOUT}/alignments/${BAM_PREFIX}.raw.sam

# FindCoveredIntervals: these are used in GATK Calling to help speed things up
echo " NGSeasy: Finding Covered Intervals : minimum coverage of 4 " `date`
  java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FindCoveredIntervals -R ${REFGenomes}/human_g1k_v37.fasta \
  -I ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.bam  \
  -o ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list \
  --coverage_threshold 4;

# BAM to BED
echo " NGSeasy: Converting BAM To BED " `date`

# pulls out paired only reads with min qual 10. Set low as donwstream tools filter these regions
  /usr/local/pipeline/samtools-0.1.19/samtools \
  view \
  -b -h -q 10 -F 1796 \
  ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.bam  | /usr/local/pipeline/bedtools2/bin/bedtools bamtobed  -i stdin \
  > ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.bed;
  
  /usr/local/pipeline/bedtools2/bin/bedtools merge -i ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.bed > ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.merged.bed
  

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

## Depending in NGS type Run PicardTools CollectTargetedPcrMetrics
## NB This is rough and ready quicl fix to get coverage ver genome bins, annotated exomes and targetd sequene files

echo " NGSeasy: Run CollectTargetedPcrMetric " `date`

if [ "${NGS_TYPE}" == "TGS" ]
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
    
elif [ "${NGS_TYPE}" == "WEX" ]
then
  echo " NGSeasy: Finding WEX Coverage over annotated exons/genes from ensembl " `date`
    java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/CollectTargetedPcrMetrics.jar \
    TMP_DIR=${SOUT}/tmp \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=100000 \
    INPUT=${SOUT}/alignments/${BAM_PREFIX}.bam \
    OUTPUT=${SOUT}/reports/${BAM_PREFIX}.target_coverage \
    REFERENCE_SEQUENCE=${REFGenomes}/human_g1k_v37.fasta \
    AMPLICON_INTERVALS=${REFGenomes}/ensgenes_protein_coding_genes.bed \
    TARGET_INTERVALS=${REFGenomes}/ensgenes_protein_coding_genes.bed \
    METRIC_ACCUMULATION_LEVEL=ALL_READS \
    PER_TARGET_COVERAGE=${SOUT}/reports/${BAM_PREFIX}.per_target_coverage;
    
elif [ "${NGS_TYPE}" == "WGS" ]
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

  echo " NGSeasy: Something Went wrong! NGS_TYPE not found! Check your config file "
  exit 1 

fi

echo ""
echo "................................................"
echo " NGSeasy: END Post Alignmnet QC " `date`
echo "................................................"
echo ""
 

##---------------------- VARIANT CALLING SINGLE SAMPLE ----------------------------------------------------------------##
echo "................................................"
echo " NGSeasy: START SNP and Small INDEL Calling " `date`
echo "................................................"
echo ""
echo " NGSeasy: NOTE: All tools are set to call variants over targeted regions created from BAM file to help speed things up (BED file format or GATK interval lists) "
echo ""

  KNOWN_SNPS_b138=/usr/local/pipeline/gatk_resources/dbsnp_138.b37.vcf
  
echo " NGSeasy: Using ${KNOWN_SNPS_b138} for SNP annotations in GATK " `date`

if [ "${VARCALLER}" == "freebayes" ]
then

echo " NGSeasy: Starting Variant Calling using Freebayes " `date`
  /usr/local/pipeline/freebayes/bin/freebayes \
    -f ${REFGenomes}/human_g1k_v37.fasta \
    -b ${SOUT}/alignments/${BAM_PREFIX}.bam \
    --min-coverage 4 \
    --min-mapping-quality 30 \
    --min-base-quality 20 \
    --targets ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.merged.bed \
    --genotype-qualities > ${SOUT}/alignments/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ;
    
    cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/
  
elif [ "${VARCALLER}" == "platypus" ]
then

  echo " NGSeasy: Starting Variant Calling using Platypus " `date`
    if [ ${NGS_TYPE} == "TGS" ]
    then
    echo " NGSeasy: NGS_TYPE is Targeted so no duplicate filtering  " `date`
    # for exome/whole genome data (no duplicate filtering)
      python /usr/local/pipeline/Platypus_0.7.4/Platypus.py callVariants \
      --nCPU ${NCPU} \
      --bamFiles=${SOUT}/alignments/${BAM_PREFIX}.bam \
      --refFile=${REFGenomes}/human_g1k_v37.fasta \
      --output=${SOUT}/alignments/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
      --filterDuplicates=0 \
      --regions=${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list;
      
      cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/
      
     else
	python /usr/local/pipeline/Platypus_0.7.4/Platypus.py callVariants \
	  --nCPU ${NCPU} \
	  --bamFiles=${SOUT}/alignments/${BAM_PREFIX}.bam \
	  --refFile=${REFGenomes}/human_g1k_v37.fasta \
	  --output=${SOUT}/alignments/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
	  --regions=${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list;
	  
	  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/;
    fi
	  
elif [ ${VARCALLER} == "gatk_ug" ]
then
  # UnifiedGenotyper EMIT_ALL_CONFIDENT_SITES
  echo " NGSeasy: Running GATK UnifiedGenotyper " `date`
  java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
  -I ${SOUT}/alignments/${BAM_PREFIX}.bam \
  -L ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list \
  -o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
  -stand_call_conf 30 \
  -stand_emit_conf 10 \
  --dbsnp ${KNOWN_SNPS_b138} \
  -dcov 250 \
  --genotype_likelihoods_model BOTH \
  --genotyping_mode DISCOVERY \
  --output_mode ${GTMODE} \
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
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/

elif [ "${VARCALLER}" == "gatk_hc" ]
then 
  echo " NGSeasy: Running GATK HaplotypeCaller (THIS TAKES A LOOOONG TIME) " `date`
  ## HaplotypeCaller Standard EMIT_ALL_CONFIDENT_SITES EMIT_VARIANTS_ONLY
  java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
  -I ${SOUT}/alignments/${BAM_PREFIX}.bam \
  -L ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list \
  -o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
  -stand_call_conf 30 \
  -stand_emit_conf 10 \
  --dbsnp ${KNOWN_SNPS_b138} \
  -dcov 250 \
  --genotyping_mode DISCOVERY \
  --output_mode ${GTMODE} \
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
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/
  
elif [ ${VARCALLER} == "gatk_hc_gvcf" ]
then
  echo " NGSeasy: Running GATK HaplotypeCaller GVCF (THIS TAKES A VERY LOOOONG TIME)" `date` 
  ## HaplotypeCaller GVCF
  java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
  -I ${SOUT}/alignments/${BAM_PREFIX}.bam \
  -L ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list \
  -o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.g.vcf \
  -stand_call_conf 30 \
  -stand_emit_conf 10 \
  --dbsnp ${KNOWN_SNPS_b138} \
  -dcov 250 \
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
  cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.g.vcf ${PROJECT_DIR}/cohort_vcfs/

else 
  echo " NGSeasy: Something Went wrong! Check your config file "
  exit 1 
fi

## TO ADD:- LCR filters and HIGH COV LOW MAP Q FILTERS AND ANNOTATIONS
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














