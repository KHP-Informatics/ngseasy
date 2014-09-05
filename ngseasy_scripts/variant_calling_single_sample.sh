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
  
##---------------------- VARIANT CALLING SINGLE SAMPLE ----------------------------------------------------------------##
echo "................................................"
echo " NGSeasy: START SNP and Small INDEL Calling " `date`
echo "................................................"
echo ""
echo " NGSeasy: NOTE: All tools are set to call variants over targeted regions created from BAM file to help speed things up  "
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
    --targets ${SOUT}/alignments/${SOUT}/alignments/${BAM_PREFIX}.merged.bed \
    --genotype-qualities > ${SOUT}/alignments/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ;
    
    cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/
  
elif [ "${VARCALLER}" == "platypus" ]
then

  echo " NGSeasy: Starting Variant Calling using Platypus " `date`
    if [ ${NGS_TYPE} == "TGS" ]
    then
    echo " NGSeasy: NGS_TYPE is Targeted so no duplicate filtering  " `date`
    # for exome/whole genome data no duplicate filtering
      python /usr/local/pipeline/Platypus_0.7.8/Platypus.py callVariants \
      --nCPU ${NCPU} \
      --bamFiles=${SOUT}/alignments/${BAM_PREFIX}.bam \
      --refFile=${REFGenomes}/human_g1k_v37.fasta \
      --output=${SOUT}/alignments/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
      --filterDuplicates=0 \
      --regions=${SOUT}/reports/${BAM_PREFIX}.CoveredIntervals_x4.list;
      
      cp -v ${SOUT}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf ${PROJECT_DIR}/cohort_vcfs/
      
     else
	python /usr/local/pipeline/Platypus_0.7.8/Platypus.py callVariants \
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
  echo " NGSeasy: Running GATK HaplotypeCaller THIS TAKES A LOOOONG TIME " `date`
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
  echo " NGSeasy: Running GATK HaplotypeCaller GVCF THIS TAKES A VERY LOOOONG TIME" `date` 
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