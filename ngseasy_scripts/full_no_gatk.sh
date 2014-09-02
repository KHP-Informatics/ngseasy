#!/bin/bash
echo  "................................................"
echo  " START NGSeasy Pipleine " `date`
echo  "................................................"

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
echo " Setting OUTPUT directory ${SOUT}"
echo ""

## BAM PREFIX 
BAM_PREFIX=${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE}

##---------------------- FASTQ-QC ----------------------## 

## symlink orginal fastq file for sample x to sample x fastq directory
echo  " Copying fastq files from ${FASTQDIR}/ to ${SOUT}/fastq/ " `date`

cp -v ${FASTQDIR}/${FASTQ1} ${SOUT}/fastq/${FASTQ1}
cp -v ${FASTQDIR}/${FASTQ2} ${SOUT}/fastq/${FASTQ2}


## set new names for copied fastq files
rawFASTQ1=`basename ${SOUT}/fastq/${FASTQ1} _1.fq.gz`
rawFASTQ2=`basename ${SOUT}/fastq/${FASTQ2} _2.fq.gz`

echo  " Fastq Basename : [$rawFASTQ1] "

# FASTQC on raw files

echo  " Run Pre-Alignment QC on raw Fastq files " `date`
/usr/local/pipeline/FastQC/fastqc \
--threads ${NCPU} \
--extract \
--quiet \
--dir ${SOUT}/tmp \
--outdir ${SOUT}/fastq \
${SOUT}/fastq/${rawFASTQ1}_1.fq.gz \
${SOUT}/fastq/${rawFASTQ2}_2.fq.gz


echo  " Start Filtering/Trimming" `date`
# Trimmomatic paired output
qcdPeFASTQ1=${SOUT}/fastq/${rawFASTQ1}_1.filtered.fq.gz
qcdPeFASTQ2=${SOUT}/fastq/${rawFASTQ2}_2.filtered.fq.gz
# Trimmomatic unpaired ouput
qcdSeFASTQ1=${SOUT}/fastq/${rawFASTQ1}_1.unpaired.fq.gz
qcdSeFASTQ2=${SOUT}/fastq/${rawFASTQ2}_2.unpaired.fq.gz

if [ ${ALIGNER} != "novoalign" ]; then
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
MINLEN:50;
else
cp ${SOUT}/fastq/${rawFASTQ1}_1.fq.gz ${qcdPeFASTQ1}
cp ${SOUT}/fastq/${rawFASTQ1}_2.fq.gz ${qcdPeFASTQ2}
fi

echo  " Run Pre-Alignment QC on Filtered/Trimmed Fastq files " `date`
# FASTQC on paired trimmed files
/usr/local/pipeline/FastQC/fastqc \
--threads ${NCPU} \
--extract \
--quiet \
--dir ${SOUT}/tmp \
--outdir ${SOUT}/fastq \
${qcdPeFASTQ1} ${qcdPeFASTQ2};

##---------------------- ALIGNMENT ----------------------##

echo  " Start Alignment" `date`

if [ ${ALIGNER} == "bwa" ]; then
# BWA alignment
echo  " Running bwa " `date`
/usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} \
> ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
elif [ ${ALIGNER} == "bowtie2" ]; then
echo  " Running bowtie2 " `date`
# Bowtie2 alignment
/usr/local/pipeline/bowtie2-2.2.3/bowtie2 -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x ${REFGenomes}/human_g1k_v37 -1 ${qcdPeFASTQ1} -2 ${qcdPeFASTQ2} -S ${SOUT}/alignments/${BAM_PREFIX}.raw.sam; # 5GB RAM
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam > ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
elif [ ${ALIGNER} == "novoalign" ]; then
echo  " Running novolaign " `date`
# Novoalign alignment 
# TO DO: Use raw${SOUT}/fastq/${FASTQ1} ${SOUT}/fastq/${FASTQ2} as novoalign does this all internally if provideing adapter lists
# ADD -a @ADAPTERS
/usr/local/pipeline/novocraft -d ${REFGenomes}/human_g1k_v37.fasta -f ${qcdPeFASTQ1} ${qcdPeFASTQ2} -F STDFQ --Q2Off --3Prime -g 40 -x 6 -r All -i PE 300,150 -c ${NCPU} -k -K ${SOUT}/alignments/${BAM_PREFIX}.K.stats -o SAM > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
elif [ ${ALIGNER} == "stampy" ]; then
echo  " Running stampy " `date`
# stampy alignment
echo  " Running stampy bwa "
/usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} \
> ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam;

echo  " Running sam to bam stampy bwa "
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam;
rm ${SOUT}/alignments/${BAM_PREFIX}.tmp.sam

echo  " Running sort bam stampy bwa "
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam ${SOUT}/alignments/${BAM_PREFIX}.tmpsort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.tmpsort.bam;
rm ${SOUT}/alignments/${BAM_PREFIX}.tmp.bam

echo  " Running Stampy aligner on bwa bam"
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

echo  " Running Stampy sam to bam "
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;

echo  " Remove Stampy intermediate bam files " 
rm ${SOUT}/alignments/${BAM_PREFIX}.raw.sam
rm ${SOUT}/alignments/${BAM_PREFIX}.raw.bam
rm ${SOUT}/alignments/${BAM_PREFIX}.tmp

fi

echo  " $ALIGNER Complete " ``date``
 
##---------------------- RAW ALIGNMENT PROCESSING ----------------------##

echo  " Adding Read Group Information " `date`

platform_unit=`cat ${qcdPeFASTQ1} | head -1 | perl -p -i -e 's/:/\t/' | cut -f 1 | perl -p -i -e 's/@//g'`

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

echo  " Marking Duplicate Reads " `date`
# MarkDuplicates
java -XX:ParallelGCThreads=${NCPU} -Xmx6g -jar /usr/local/pipeline/picardtools/picard-tools-1.115/MarkDuplicates.jar \
TMP_DIR=${SOUT}/tmp \
VALIDATION_STRINGENCY=SILENT \
MAX_RECORDS_IN_RAM=100000 \
CREATE_INDEX=true \
REMOVE_DUPLICATES=false \
ASSUME_SORTED=true \
INPUT=${SOUT}/alignments/${BAM_PREFIX}.addrg.bam \
OUTPUT=${SOUT}/alignments/${BAM_PREFIX}.dupemk.bam \
METRICS_FILE=${SOUT}/alignments/${BAM_PREFIX}.dupemk;

echo  " Removing Intermediate SAM/BAM Files " `date`

rm -v ${SOUT}/alignments/${BAM_PREFIX}.addrg.ba*
rm -v ${SOUT}/alignments/${BAM_PREFIX}.sort.ba*
rm -v ${SOUT}/alignments/${BAM_PREFIX}.raw.sam

# FindCoveredIntervals: these are used in GATK Calling to help speed things up
echo  " Finding Covered Intervals : minimum coverage of 4 " `date`
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FindCoveredIntervals -R ${REFGenomes}/human_g1k_v37.fasta \
-I ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.dupemk.bam  \
-o ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list \
--coverage_threshold 4;

echo  " Alignment Complete " `date`


##---------------------- POST ALIGNMENT QC ----------------------##

echo  " Start Post Alignmnet QC " `date`

# CollectMultipleMetrics
echo  " CollectMultipleMetrics " `date`
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
echo  " CollectAlignmentSummaryMetrics " `date`
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
echo  " CollectWgsMetrics " `date`
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
echo  " FlagStat " `date`
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FlagStat -R ${REFGenomes}/human_g1k_v37.fasta \
-I ${SOUT}/alignments/${BAM_PREFIX}.bam  \
-o ${SOUT}/reports/${BAM_PREFIX}.FlagStat;

# BAM to BED
echo  " BAM To BED " `date`
/usr/local/pipeline/samtools-0.1.19/samtools view -b -h -q 20 -F 1796  ${SOUT}/alignments/${BAM_PREFIX}.bam  | /usr/local/pipeline/bedtools2/bin/bedtools bamtobed  -i stdin > ${SOUT}/reports/${BAM_PREFIX}.bed;

## Depending in NGS type Run PicardTools CollectTargetedPcrMetrics
## NB This is rough and ready quicl fix to get coverage ver genome bins, annotated exomes and targetd sequene files

echo  " Run CollectTargetedPcrMetric " `date`

if [ ${NGS_TYPE} == "TGS" ];then
echo  " Finding TGS coverage over target regions in ${BED_ANNO} " `date`
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
elif [ ${NGS_TYPE} == "WEX" ]; then
echo  " Finding WEX Coverage over annotated exons/genes from ensembl " `date`
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
elif [ ${NGS_TYPE} == "WGS" ]; then
echo  " Finding WGS Coverage over 500bp windows " `date`
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
fi

echo  " Basic Post Alignmnet QC Complete " `date`
 

##---------------------- GATK VARIANT CALLING SINGLE SAMPLE ----------------------##
echo  " Start SNP and Small INDEL Calling " `date`

KNOWN_SNPS_b138=/usr/local/pipeline/gatk_resources/dbsnp_138.b37.vcf

echo  " Starting Variant Calling using Freebayes " `date`
if [ ${VARCALLER} == "freebayes" ]; then
/usr/local/pipeline/freebayes/bin/freebayes \
-f ${REFGenomes}/human_g1k_v37.fasta \
-b ${SOUT}/alignments/${BAM_PREFIX}.bam \
--min-coverage 4 \
--min-mapping-quality 30 \
--min-base-quality 20 \
--genotype-qualities > ${SOUT}/alignments/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ;
cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/
elif [ ${VARCALLER} == "platypus" ]; then
echo  " Starting Variant Calling using Platypus " `date`
if [ ${NGS_TYPE} == "TGS" ];then
python /usr/local/pipeline/Platypus_0.7.4/Platypus.py callVariants \
--bamFiles=${SOUT}/alignments/${BAM_PREFIX}.bam \
--refFile=${REFGenomes}/human_g1k_v37.fasta \
--output=${SOUT}/alignments/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
--filterDuplicates=0;
cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/
# for exome/whole genome data (no duplicate filtering)
else
python /usr/local/pipeline/Platypus_0.7.4/Platypus.py callVariants \
--bamFiles=${SOUT}/alignments/${BAM_PREFIX}.bam \
--refFile=${REFGenomes}/human_g1k_v37.fasta \
--output=${SOUT}/alignments/${BAM_PREFIX}..raw.snps.indels.${VARCALLER}.vcf;
cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/
elif [ ${VARCALLER} == "gatk_ug" ]; then
# UnifiedGenotyper EMIT_ALL_CONFIDENT_SITES
echo  " Running GATK UnifiedGenotyper " `date`
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
elif [ ${VARCALLER} == "gatk_hc" ]; then 
echo  " Running GATK HaplotypeCaller " `date`
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
elif [ ${VARCALLER} == "gatk_hc_gvcf" ]; then
echo  " Running GATK HaplotypeCaller GVCF " `date` 
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
fi

echo  " Variant Calling Complete " `date` 

echo  "......................"
echo  " END NGSeasy Pipeline " `date`
exit


