#!/bin/bash

#################################################################
## Tool Configurations

## Aligners

local BWA_OPTIONS=" mem -M "

local STAMPY_OPTIONS=" --bamkeepgoodreads --sanger --bwamark --baq --gapopen=40 --gapextend=6 --noautosense --insertsize=500 --insertsd=100 "

local SNAP_OPTIONS=" -b -M -s 50 1000 -H 300000 -h 300 -d 15 -mcp 6000000 -map -pre -I "

local NOVOALIGN_OPTIONS=" -F STDFQ --3Prime -g 40 -x 6 -i PE 500,100 "

local BOWTIE2_OPTIONS=" -I 50 -X 10000 "

#################################################################
## Tool Configurations


## Variant Callers

local FREEBAYES_OPTIONS=" --min-coverage 4 --min-mapping-quality 20 --min-base-quality 10 --min-repeat-entropy 1 --genotype-qualities "

local PLATYPUS_OPTIONS=" --assemble=1 --assembleAll=1 \
--assemblyRegionSize=1500 \
--minReads=4 \
--maxGOF=30 \
--bufferSize=50000 \
--maxReads=10000000 \
--minPosterior=5 \
--minMapQual=20 --minBaseQual=10 \
--maxSize 10000 \
--maxVariants=8 \
--maxHaplotypes=50 \
--filterReadsWithDistantMates=1 \
--hapScoreThreshold 10 \
--scThreshold 0.99 \
--filteredReadsFrac 0.7 \
--rmsmqThreshold 20 \
--qdThreshold 0 \
--abThreshold 0.0001 \
--minVarFreq 0.0 "

local HAPLOTYPECALLER_OPTIONS=" -stand_call_conf 30 -stand_emit_conf 10 \
--output_mode EMIT_VARIANTS_ONLY \
-dcov 250 \
-minPruning 10 \
-pairHMM VECTOR_LOGLESS_CACHING \
--genotyping_mode DISCOVERY \
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
--annotation VariantType"
