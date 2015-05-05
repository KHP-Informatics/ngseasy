

########################################################################################################
## Set version and run date
#
NGSEASYVERSION="1.0"
RUNDATE=`date +"%d%m%y"`
NGSEASY_STEP="ngseasy_variant_calling"

########################################################################################################
## docker run command
#
DOCKER_RUN="docker run -P -w /home/pipeman -e HOME=/home/pipeman -e USER=pipeman --user pipeman"


########################################################################################################
## PLATYPUS_OPTIONS
#
PLATYPUS_OPTIONS=" --assemble=1 --assembleAll=1 \
--assemblyRegionSize=1500 \
--minReads=4 \
--maxGOF=30 \
--bufferSize=50000 \
--maxReads=10000000 \
--minPosterior=5 \
--minMapQual=20 \
--minBaseQual=10 \
--maxSize 10000 \
--maxVariants=8 \
--maxHaplotypes=50 \
--filterReadsWithDistantMates=1 \
--hapScoreThreshold 10 \
--scThreshold 0.99 \
--filteredReadsFrac 0.9 \
--rmsmqThreshold 20 \
--qdThreshold 0 \
--abThreshold 0.0001 \
--minVarFreq 0.0 "

########################################################################################################
## PLATYPUS RUN

NCPU=30

    ${DOCKER_RUN} \
    -v ${PROJECT_DIR}:/home/pipeman/ngs_projects \
    --name platypus_${BAM_PREFIX} \
    -t compbio/ngseasy-platypus:${NGSEASYVERSION} /bin/bash -c \
    "time python /usr/local/pipeline/Platypus/bin/Platypus.py callVariants \
      --nCPU ${NCPU} \
      --bamFiles=${BAMFILE} \
      --refFile=${REFFASTA} \
      --output=${SOUTDocker}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf \
      --logFileName=${SOUTDocker}/vcf/${BAM_PREFIX}.raw.snps.indels.${VARCALLER}.vcf.log \
      --filterDuplicates=0 \
      ${PLATYPUS_OPTIONS}



