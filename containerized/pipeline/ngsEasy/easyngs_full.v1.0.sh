# refgenomes and contaminant lists
REFGenomes='/home/pipeman/reference_genomes'
adapter_fa='/home/pipeman/reference_genomes/contaminant_list.fa'
gatk_resources='/home/pipeman/gatk_resources'
FASTQDIR='/home/pipeman/FASTQ_STAGING'
##----------------------------------------------------------------------------------------------##


#---------------------------------------------CONFIG------------------------------------------------#
# get the following from a config file

# SAMPLE INFORMATION
POJECT_ID='example'
SAMPLE_ID='test'
FASTQ1='WEXs_R1.fastq'
FASTQ2='WEXs_R2.fastq'
PROJECT_DIR='/home/pipeman/ngs_projects'
DNA_PREP_LIBRARY_ID="EXOME"

## NGS TYPE
NGS_TYPE='WEX'
NGS_PLATFORM='ILLUMINA'
BED_ANNO='NULL'

## Pipeline to call
PIPELINE='full'

## ALIGNER
ALIGNER='Bowtie2'

# GENOTYPING MODE
GTMODE='EMIT_VARIANTS_ONLY'

## CLEANUP TRUE/FALSE
CLEANUP='TRUE'

## NCPU
NCPU=4

## GATK GENOTYPE Settings
stand_call_conf=30
stand_emit_conf=10

## DATE
DATE=`date +"%m%d%y"`
#---------------------------------------------CONFIG------------------------------------------------#



##----------------------------------- SET UP PROJECT DIRECTORY ---------------------------##
## FORU UNIX THIS CAN BE CALLED FROM SCRIPT
## NEED TO DO THIS MANUALLY IN CMD/DOS PROMPT
## MAKE PROJECT DIRs
mkdir ${PROJECT_DIR}
mkdir ${PROJECT_DIR}/${POJECT_ID}
mkdir ${PROJECT_DIR}/${POJECT_ID}/config_files
mkdir ${PROJECT_DIR}/${POJECT_ID}/cohort_vcfs
mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}
mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/fastq
mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/tmp
mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/alignments
mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/vcf
mkdir ${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}/reports
echo " made project dir"
tree ${PROJECT_DIR}/${POJECT_ID}
##----------------------------------- SET UP PROJECT DIRECTORY ---------------------------##



##############################
## BEGING FULL NGS PIPELINE ##
##############################

## OUTPUT SAMPLE DIR
SOUT=${PROJECT_DIR}/${POJECT_ID}/${SAMPLE_ID}

## BAM PREFIX
BAM_PREFIX=${SAMPLE_ID}.${NGS_TYPE}.${NGS_PLATFORM}.${ALIGNER}.${DATE}

##---------------------- FASTQ-QC ----------------------##

## copy orginal fastq file for sample x to sample x fastq directory
cp -v ${FASTQDIR}/${FASTQ1} ${SOUT}/fastq/
cp -v ${FASTQDIR}/${FASTQ2} ${SOUT}/fastq/

## set new names for copied fastq files
rawFASTQ1=${SOUT}/fastq/${FASTQ1}
rawFASTQ2=${SOUT}/fastq/${FASTQ2}

# FASTQC on raw files
/usr/local/pipeline/FastQC/fastqc --threads ${NCPU} --extract --quiet --dir ${SOUT}/tmp --outdir ${SOUT}/fastq ${rawFASTQ1} ${rawFASTQ2};

# Trimmomatic paired output
qcdPeFASTQ1=${rawFASTQ1}.qcd.pe.fq
qcdPeFASTQ2=${rawFASTQ2}.qcd.pe.fq
# Trimmomatic unpaired ouput
qcdSeFASTQ1=${rawFASTQ1}.qcd.se.fq
qcdSeFASTQ2=${rawFASTQ2}.qcd.se.fq
# run Trimmomatic
# set min lengh to 35bp
java -jar /usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads ${NCPU} -phred64 -trimlog ${SOUT}/fastq/trimm_qc.log \
${rawFASTQ1} ${rawFASTQ2} \
${qcdPeFASTQ1} ${qcdSeFASTQ1} \
${qcdPeFASTQ2} ${qcdSeFASTQ2} \
ILLUMINACLIP:${adapter_fa}:2:30:10:5:true LEADING:3 TRAILING:3 MINLEN:50;

# FASTQC on paired trimmed files
/usr/local/pipeline/FastQC/fastqc --threads ${NCPU} --extract --quiet --dir ${SOUT}/tmp --outdir ${SOUT}/fastq ${qcdPeFASTQ1} ${qcdPeFASTQ2};

##---------------------- ALIGNMENT ----------------------##

if [ "${ALIGNER}" == "bwa" ]; then
# BWA alignment
/usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
fi

if [ "${ALIGNER}" == "Bowtie2" ]; then
# Bowtie2 alignment
/usr/local/pipeline/bowtie2-2.2.3/bowtie2 -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 -x ${REFGenomes}/human_g1k_v37 -1 ${qcdPeFASTQ1} -2 ${qcdPeFASTQ2} -S ${SOUT}/alignments/${BAM_PREFIX}.raw.sam; # 5GB RAM
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam > ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
fi

if [ "${ALIGNER}" == "novoalign" ]; then
# Novoalign alignment
ngsNovoalign -d ${REFGenomes}/human_g1k_v37.fasta -f ${qcdPeFASTQ1} ${qcdPeFASTQ2} -F STDFQ --Q2Off --3Prime -g 40 -x 6 -r All -i PE 300,150 -c ${NCPU} -k -K ${SOUT}/alignments/${BAM_PREFIX}.K.stats -o SAM > ${SOUT}/alignments/${BAM_PREFIX}.raw.sam;
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
fi

## TO DO : NEED TO TEST STAMPY
if [ "${ALIGNER}" == "stampy" ]; then
# stampy alignment
/usr/local/pipeline/bwa-0.7.10/bwa mem -M -t ${NCPU} ${REFGenomes}/human_g1k_v37.fasta ${qcdPeFASTQ1} ${qcdPeFASTQ2} > ${SOUT}/alignments/${BAM_PREFIX}.raw0.sam;
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw0.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort0.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort0.bam;
/usr/local/pipeline/samtools-0.1.19/samtools view -H ${SOUT}/alignments/${BAM_PREFIX}.sort0.bam > ${SOUT}/alignments/${BAM_PREFIX}.sort0.bam.header;
/usr/local/pipeline/stampy-1.0.23/stampy.py -g human_g1k_v37 -h human_g1k_v37 -t${NCPU} --bamsortprefix ${SOUT}/tmp --bamkeepgoodreads -M ${SOUT}/alignments/${BAM_PREFIX}.sort0.bam -o ${SOUT}/alignments/${BAM_PREFIX}.sort0_unmapped.bam -f sam
/usr/local/pipeline/stampy-1.0.23/stampy.py -g human_g1k_v37 -h human_g1k_v37 -t${NCPU} --bamsortprefix ${SOUT}/tmp -M ${SOUT}/alignments/${BAM_PREFIX}.sort0.bam ${SOUT}/alignments/${BAM_PREFIX}.sort0_unmapped.bam -o ${SOUT}/alignments/${BAM_PREFIX}.raw.sam -f sam
/usr/local/pipeline/samtools-0.1.19/samtools view -bhS ${SOUT}/alignments/${BAM_PREFIX}.raw.sam ${SOUT}/alignments/${BAM_PREFIX}.raw.bam;
/usr/local/pipeline/samtools-0.1.19/samtools sort -f   ${SOUT}/alignments/${BAM_PREFIX}.raw.bam ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
/usr/local/pipeline/samtools-0.1.19/samtools index     ${SOUT}/alignments/${BAM_PREFIX}.sort.bam;
fi

##---------------------- RAW ALIGNMENT PROCESSING ----------------------##

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
RGPU=`cat ${qcdPeFASTQ1} | head -1 | perl -p -i -e 's/:/\t/' | cut -f 1 | perl -p -i -e 's/@//g'` \
RGSM=${SAMPLE_ID} \
RGDT=${DATE} \
INPUT=${SOUT}/alignments/${BAM_PREFIX}.sort.bam \
OUTPUT=${SOUT}/alignments/${BAM_PREFIX}.addrg.bam;

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


##---------------------- GATK CLEANING ----------------------##

# RealignerTargetCreator

## gunzip /usr/local/pipeline/gatk_resources/Mills_and_1000G_gold_standard.indels.b37***
KNOWN_INDELS=/usr/local/pipeline/gatk_resources/Mills_and_1000G_gold_standard.indels.b37.vcf # CHECK IT CAN READ gz'd files

# RealignerTargetCreator
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${REFGenomes}/human_g1k_v37.fasta -nt ${NCPU} \
--known ${KNOWN_INDELS} \
-I ${SOUT}/alignments/${BAM_PREFIX}.dupemk.bam \
-o ${SOUT}/alignments/${BAM_PREFIX}.dupemk.bam.IndelRealigner.intervals;

# IndelRealigner
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R ${REFGenomes}/human_g1k_v37.fasta \
-known ${KNOWN_INDELS} \
-I ${SOUT}/alignments/${BAM_PREFIX}.dupemk.bam \
-targetIntervals ${SOUT}/alignments/${BAM_PREFIX}.dupemk.bam.IndelRealigner.intervals \
-o ${SOUT}/alignments/${BAM_PREFIX}.realn.bam \
--consensusDeterminationModel USE_READS \
-LOD 0.4;

# index bam file
/usr/local/pipeline/samtools-0.1.19/samtools index ${SOUT}/alignments/${BAM_PREFIX}.realn.bam;

##---------------------- GATK BASE RECALIBRATION ----------------------##

KNOWN_INDELS=/usr/local/pipeline/gatk_resources/Mills_and_1000G_gold_standard.indels.b37.vcf # CHECK IT CAN READ gz'd files
KNOWN_SNPS_1000G=/usr/local/pipeline/gatk_resources/1000G_phase1.snps.high_confidence.b37.vcf
KNOWN_SNPS_b138=/usr/local/pipeline/gatk_resources/dbsnp_138.b37.vcf

# RealignerTargetCreator pre recal
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU}  \
-I ${SOUT}/alignments/${BAM_PREFIX}.realn.bam \
-o ${SOUT}/alignments/${BAM_PREFIX}.realn.bam.BaseRecalibrator.table \
--knownSites ${KNOWN_INDELS} \
--knownSites ${KNOWN_SNPS_1000G} \
--knownSites ${KNOWN_SNPS_b138};

# PrintReads BQSR recalibration
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T PrintReads -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
--baq RECALCULATE \
--baqGapOpenPenalty 40 \
--BQSR ${SOUT}/alignments/${BAM_PREFIX}.realn.bam.BaseRecalibrator.table \
-I ${SOUT}/alignments/${BAM_PREFIX}.realn.bam \
-o ${SOUT}/alignments/${BAM_PREFIX}.recal.bam;
# index bam file
/usr/local/pipeline/samtools-0.1.19/samtools index ${SOUT}/alignments/${BAM_PREFIX}.recal.bam;

# RealignerTargetCreator post recal
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
-I ${SOUT}/alignments/${BAM_PREFIX}.recal.bam \
-o ${SOUT}/alignments/${BAM_PREFIX}.recal.bam.BaseRecalibrator.table \
--knownSites ${KNOWN_INDELS} \
--knownSites ${KNOWN_SNPS_1000G} \
--knownSites ${KNOWN_SNPS_b138};

cp -v ${SOUT}/alignments/${BAM_PREFIX}.recal.bam ${SOUT}/alignments/${BAM_PREFIX}.bam
# index bam file
/usr/local/pipeline/samtools-0.1.19/samtools index ${SOUT}/alignments/${BAM_PREFIX}.bam;



##---------------------- GATK VARIANT CALLING SINGLE SAMPLE ----------------------##

KNOWN_SNPS_b138=/usr/local/pipeline/gatk_resources/dbsnp_138.b37.vcf

# UnifiedGenotyper EMIT_ALL_CONFIDENT_SITES
if [ "${VARCALLER}" == "UG" ]; then
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
-I ${SOUT}/alignments/${BAM_PREFIX}.bam \
-o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.ug.vcf \
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
cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.ug.vcf ${PROJECT_DIR}/cohort_vcfs/
fi


## HaplotypeCaller Standard EMIT_ALL_CONFIDENT_SITES EMIT_VARIANTS_ONLY
if [ "${VARCALLER}" == "HC" ]; then
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
-I ${SOUT}/alignments/${BAM_PREFIX}.bam \
-o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.hc.vcf \
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
cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.hc.vcf ${PROJECT_DIR}/cohort_vcfs/
fi


## HaplotypeCaller GVCF
if [ "${VARCALLER}" == "HCgVCF" ]; then
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T HaplotypeCaller -R ${REFGenomes}/human_g1k_v37.fasta -nct ${NCPU} \
-I ${SOUT}/alignments/${BAM_PREFIX}.bam \
-o ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.hc.g.vcf \
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
# copy vcf to cohort vcf directory
cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.hc.g.vcf ${PROJECT_DIR}/cohort_vcfs/
fi
## -minPruning 10 -dcov 250

##---------------------- ALIGNMENT QC ----------------------##

# CollectMultipleMetrics
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

# FindCoveredIntervals
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FindCoveredIntervals -R ${REFGenomes}/human_g1k_v37.fasta \
-I ${SOUT}/alignments/${BAM_PREFIX}.bam  \
-o ${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals.list \
--coverage_threshold 4;

# FlagStat
java -Xmx6g -Djava.io.tmpdir=${SOUT}/tmp -jar /usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar -T FlagStat -R ${REFGenomes}/human_g1k_v37.fasta \
-I ${SOUT}/alignments/${BAM_PREFIX}.bam  \
-o ${SOUT}/reports/${BAM_PREFIX}.FlagStat;

# BAM to BED
/usr/local/pipeline/samtools-0.1.19/samtools view -b -h -q 20 -F 1796  ${SOUT}/alignments/${BAM_PREFIX}.bam  | /usr/local/pipeline/bedtools2/bin/bedtools bamtobed  -i stdin > ${SOUT}/reports/${BAM_PREFIX}.bed;

## TO ADD ####################################################################################################
CollectTargetedPcrMetrics
fastq to BAM
plots of coverage
vcf reports
var annotation
multi sample calling
bed files annotations 
calling given snps only
