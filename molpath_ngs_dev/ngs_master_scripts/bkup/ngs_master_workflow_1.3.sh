#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -pe multi_thread 1
#$ -l h_vmem=1G
#$ -p -0.99999999999999999999999999999999999999999999999999999999999999999
#$ -j y
#------------------------------------------------------------------------#


#############################################################################################
# -- Authors: Stepgen Newhouse, Amos Folarin, Aditi Gulati                                  #
# -- Organisation: KCL/SLaM/NHS                                                             #
# -- Email: stephen.j.newhouse@gmail.com, amosfolarin@gmail.com,aditi.gulati@nhs.ne         #
# -- Verion: 1.3                                                                            #
# -- Date: 11/09/2013                                                                       #
# -- DESC: NGS pipeline to perform SE/PE Alignments & GATK cleaning                         #
#############################################################################################

# called using "Rscript call_ngs_master_workflow.R <config_file>"

# this calls the following with options read in from R and fired off using R's system() command

# qsub ngs_master_workflow.sh <fastq_prefix> <sample_name> <qual_type> <RGID> <RGLB> <RGPL> <RGPU> <RGSM> <RGCN> <RGDS> <RGDT> <PE> <bed_list> [to add: options for pipelines steps, bed targets and geno all sites and knowns]

# ALL OPTIONS REQUIRED
# -- fastq_prefix=string	fastq prefix for Sample ID ie <fastq_prefix>_1.fastq <fastq_prefix>_2.fastq
# -- sample_name=string		Sample ID
# -- qual_type=string		Base quality coding for novoalign ie STFQ, ILMFQ, ILM1.8
# -- RGID=String		Read Group ID Required. <PROJECT_NAME>
# -- RGLB=String		Read Group Library Required. <PATIENT_ID>.<RGPL>.<>
# -- RGPL=String		Read Group platform (e.g. illumina, solid) Required.
# -- RGPU=String		Read Group platform unit (eg. run barcode) Required.
# -- RGSM=String		Read Group sample name Required. PATIENT_ID
# -- RGCN=String		Read Group sequencing center name Required.
# -- RGDS=String		Read Group description Required.
# -- RGDT=Iso8601Date		Read Group run date Required.
# -- PE=1 or 0			Indicates PE or SE
# -- bed_list=string		name/prefix of target bedfile

#------------------------------------------------------------------------#
# Set environmental variables required for child processes (all )
#------------------------------------------------------------------------#

########################################
## target bed file for coverage calcs ##
########################################

## export bed_list="cancer"  ## TO Add:- genes, exons and other custom targets from ENSEMBL

#######
# QUE #
#######
export queue_name="short.q,long.q"
export email_contact="stephen.newhouse@kcl.ac.uk"

#####################
# mem and cpu vars ##
#####################
## Novoalign
export novo_cpu=8
export novo_mem=3

## Java & Picardtools
export sge_h_vmem=8
export java_mem=6

## Java & GATK
export gatk_h_vmem=8
export gatk_java_mem=6

###############
## ngs tools ##
###############
export ngs_picard="/share/apps/picard-tools_1.91/jar"
export ngs_gatk="/share/apps/gatk_2.5-2"  ## gatk_2.7-2 needs java 1.7
export ngs_novo="/share/apps/novocraft_current/bin/" ## Novoalign V3.01.01
export ngs_samtools="/share/apps/samtools_0.1.18/bin"

##################
## pipeline dir ##
##################
export ngs_pipeline="/home/snewhousebrc/scratch/pipelines/ngs/molpath_ngs"

######################
## reference genomes #
######################
export reference_genome_novoindex="/scratch/data/reference_genomes/gatk_resources/b37/human_g1k_v37.fasta.novoindex"
export reference_genome_seq="/scratch/data/reference_genomes/gatk_resources/b37/human_g1k_v37.fasta"

############################
## ref vcf files for gatk ##
############################
## needs updating!!!!!!
# indels
export b37_1000G_biallelic_indels="/scratch/data/reference_genomes/gatk_resources/b37/1000G_biallelic.indels.b37.vcf"
export b37_Mills_Devine_2hit_indels_sites="/scratch/data/reference_genomes/gatk_resources/b37/Mills_Devine_2hit.indels.b37.sites.vcf"
export b37_Mills_Devine_2hit_indels="/scratch/data/reference_genomes/gatk_resources/b37/Mills_Devine_2hit.indels.b37.vcf"
# snps
export b37_1000G_omni2_5="/scratch/data/reference_genomes/gatk_resources/b37/1000G_omni2.5.b37.vcf"
export b37_hapmap_3_3="/scratch/data/reference_genomes/gatk_resources/b37/hapmap_3.3.b37.sites.vcf"
export b37_dbsnp_132_excluding_sites_after_129="/scratch/data/reference_genomes/gatk_resources/b37/dbsnp_132.b37.excluding_sites_after_129.vcf"
export b37_dbsnp="/scratch/data/reference_genomes/gatk_resources/b37/dbsnp_132.b37.vcf"

#######################
## Path to fastq dir ##
#######################
export fastq_dir="/scratch/project/pipelines/ngs_pipeline_dev/aditi_fastq"

###########################
## path to final aln dir ##
###########################
## export aln_dir="/scratch/project/pipelines/ngs_pipeline_dev/ngs_molpath_sjn"
export aln_dir="/home/snewhousebrc/scratch/pipelines/ngs"

###########################
## path to tmp aln dir   ##
###########################
export ngstmp="/home/snewhousebrc/scratch/ngs_temp"

#############################
## get and set all options ##
#############################
fastq_prefix=${1}
sample_name=${2}.${5}.${6}
qual_type=${3}  ## Base quality coding for novoalign ie STFQ, ILMFQ, ILM1.8
mRGID=${4}	#Read Group ID Required.
mRGLB=${5}	#Read Group Library Required.
mRGPL=${6}	#Read Group platform (e.g. illumina, solid,IONTORRENT) Required.
mRGPU=${7}	#Read Group platform unit (eg. run barcode) Required.
mRGSM=${8}	#Read Group sample name Required.
mRGCN=${9}	#Read Group sequencing center name Required.
mRGDS=${10}	#Read Group description Required.
mRGDT=${11}	#Read Group run date Required.
mPE=${12} 	#Indicates PE or SE
bed_list=${13}	#target bed for coverage

#------------------------------------------------------------------------------#
# END setting environmental variables required for child processes (all )
#------------------------------------------------------------------------------#


##############################
## make sample dir for data ##
##############################

## temp dir
mkdir ${ngstmp}/${sample_name}_temp

## final data dir
mkdir ${aln_dir}/${sample_name}
mkdir ${aln_dir}/${sample_name}/sge_out


sample_temp=${ngstmp}/${sample_name}_temp

sample_dir=${aln_dir}/${sample_name}

#------------------------------------------------------------------------------#

######################################
# moving to sample directory #########
######################################

echo "moving to sample directory "  ${sample_dir}

cd ${sample_dir}

#------------------------------------------------------------------------------#
echo ".................................................."
echo "fastq_dir : " ${fastq_dir}
echo "project dir : " ${aln_dir}
echo ".................................................."
echo "Sample working dir : " ${sample_dir}
echo "Sample temp dir : " ${sample_temp}
echo "Sample information : "  ${fastq_prefix} ${sample_name} ${mRGID} ${mRGLB} ${mRGPL} ${mRGPU} ${mRGSM} ${mRGCN} ${mRGDS} ${mRGDT}
echo ".................................................."
echo "sge ques : " ${queue_name}
echo "novo_mem : " ${novo_mem}G
echo "sge_h_vmem : " ${sge_h_vmem}G
echo "gatk_java_mem : " ${gatk_java_mem}
echo "novo_cpu : " ${novo_cpu}
echo "email_contact: " ${email_contact}
echo "ngs_pipeline : " ${ngs_pipeline}
echo ".................................................."
ls -l ${ngs_pipeline}
echo ".................................................."
echo "ngs_picard : " ${ngs_picard} 
echo "ngs_novo : " ${ngs_novo} 
echo "ngs_gatk : " ${ngs_gatk} 
echo "ngs_samtools: " ${ngs_samtools} 
echo ".................................................."
ls -l ${ngs_picard}
echo ".................................................."
ls -l ${ngs_novo}
echo ".................................................."
ls -l ${ngs_gatk}
echo ".................................................."
ls -l ${ngs_samtools}
echo ".................................................."
echo "fastq seq : " ${reference_genome_seq}
head -5  ${reference_genome_seq}
echo ".................................................."
head ${b37_1000G_biallelic_indels}
echo ".................................................."
head ${b37_Mills_Devine_2hit_indels_sites}
echo ".................................................."
head ${b37_Mills_Devine_2hit_indels}
echo ".................................................."
head ${b37_1000G_omni2_5}
echo ".................................................."
head ${b37_hapmap_3_3}
echo ".................................................."
head ${b37_dbsnp_132_excluding_sites_after_129}
echo ".................................................."
head ${b37_dbsnp_132}
echo ".................................................."
#------------------------------------------------------------------------------#

##############################
## START ALIGNMENT PIPELINE ##
##############################


#----------------------------------------------------------------------#
# 1. Align PE or SE data 
#----------------------------------------------------------------------#
##echo ">>>>>" `date` " :-> " "Aligning PE data "

if [ ${mPE} -eq 1 ] && [ ${mRGPL} -eq "illumina" ]; then

echo "Reads are PE ILLUMINA"

qsub -q ${queue_name} -N novoalign.${sample_name} -l h_vmem=${novo_mem}G -pe multi_thread ${novo_cpu} -M ${email_contact} -m beas ${ngs_pipeline}/ngs_novoalign.illumina.${qual_type}.PE.sh ${fastq_prefix} ${sample_name} ${sample_dir};

elif [ ${mPE} -eq 0 ] && [ ${mRGPL} -eq "illumina" ]; then

echo "Reads are SE ILLUMINA"

qsub -q ${queue_name} -N novoalign.${sample_name} -l h_vmem=${novo_mem}G -pe multi_thread ${novo_cpu} -M ${email_contact} -m beas ${ngs_pipeline}/ngs_novoalign.illumina.${qual_type}.SE.sh ${fastq_prefix} ${sample_name} ${sample_dir};

elif [ ${mPE} -eq 0 ] && [ ${mRGPL} -eq "IONTORRENT" ]; then

echo " Reads are SE IONTORRENT"

qsub -q ${queue_name} -N novoalign.${sample_name} -l h_vmem=${novo_mem}G -pe multi_thread ${novo_cpu} -M ${email_contact} -m beas ${ngs_pipeline}/ngs_novoalign.IONTORRENT.${qual_type}.SE.sh ${fastq_prefix} ${sample_name} ${sample_dir};

fi


#----------------------------------------------------------------------#
# 2. sam2bam
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Converting Novoalign SAM to BAM and indexing"

qsub -q ${queue_name} -N sam2bam.${sample_name} -hold_jid novoalign.${sample_name} -l h_vmem=${sge_h_vmem}G -pe multi_thread 1 -M ${email_contact} -m beas ${ngs_pipeline}/ngs_sam2bam.sh ${sample_name} ${sample_dir};


#----------------------------------------------------------------------#
# 3. SortSam
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N SortSam.${sample_name} -hold_jid sam2bam.${sample_name}  -l h_vmem=${sge_h_vmem}G -pe multi_thread 1 -M ${email_contact} -m beas ${ngs_pipeline}/ngs_SortSam.sh ${sample_name} ${sample_dir} ${sample_temp};

#----------------------------------------------------------------------#
# 4. AddOrReplaceReadGroups 
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running AddOrReplaceReadGroups" 

qsub -q ${queue_name} -N AddOrReplaceReadGroups.${sample_name} -hold_jid SortSam.${sample_name} -l h_vmem=${sge_h_vmem}G -pe multi_thread 1 -M ${email_contact} -m beas ${ngs_pipeline}/ngs_AddOrReplaceReadGroups.sh ${sample_name} ${sample_dir} \
${sample_temp} ${mRGID} ${mRGLB} ${mRGPL} ${mRGPU} ${mRGSM} ${mRGCN} ${mRGDS} ${mRGDT};


#----------------------------------------------------------------------#
# 5. MarkDuplicates
#----------------------------------------------------------------------#

###echo ">>>>>" `date` " :-> " "Running MarkDuplicates"

qsub -q ${queue_name} -N MarkDuplicates.${sample_name} -hold_jid AddOrReplaceReadGroups.${sample_name} -l h_vmem=${sge_h_vmem}G  -pe multi_thread 1 -M ${email_contact} -m beas ${ngs_pipeline}/ngs_MarkDuplicates.sh ${sample_name} ${sample_dir} ${sample_temp};

##############################
## END ALIGNMENT PIPELINE   ## 
##############################

#########################
## BEGIN GATK CLEANING ##
#########################

#----------------------------------------------------------------------#
# 7. RealignerTargetCreator
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running RealignerTargetCreator"

qsub -q ${queue_name} -N RealignerTargetCreator.${sample_name} -hold_jid MarkDuplicates.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_RealignerTargetCreator.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 8. IndelRealigner
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running IndelRealigner"

qsub -q ${queue_name} -N IndelRealigner.${sample_name} -hold_jid RealignerTargetCreator.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_IndelRealigner.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 9. BaseRecalibrator before recal
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running BaseRecalibrator before QUAL SCORE RECALIBRATION"

qsub -q ${queue_name} -N BaseRecalibrator_before.${sample_name} -hold_jid IndelRealigner.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_BaseRecalibrator_before.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 10. PrintReads = QUAL SCORE RECALIBRATION
#----------------------------------------------------------------------#

###echo ">>>>>" `date` " :-> " "Running PrintReads > QUAL SCORE RECALIBRATION "

qsub -q ${queue_name} -N PrintReads_BQSR.${sample_name} -hold_jid BaseRecalibrator_before.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_PrintReads_BQSR.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 11. BaseRecalibrator after recal
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running BaseRecalibrator after QUAL SCORE RECALIBRATION"

qsub -q ${queue_name} -N BaseRecalibrator_after.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_BaseRecalibrator_after.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 12. AnalyzeCovariates before & after recal
#----------------------------------------------------------------------#

##echo ">>>>>" `date` " :-> " "Running AnalyzeCovariates" 

qsub -q ${queue_name} -N AnalyzeCovariates_before_and_after_BQSR.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_AnalyzeCovariates_before_and_after_BQSR.sh \
${sample_name} ${sample_dir} ${sample_temp};

##############################################
##  CALL VARIANTS SINGLE SAMPLE ##############
##############################################

#----------------------------------------------------------------------#
# 13. HaplotypeCaller 
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N HaplotypeCaller.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_HaplotypeCaller.sh \
${sample_name} ${sample_dir} ${sample_temp} 30 10;

#----------------------------------------------------------------------#
# 14. UnifiedGenotyper 
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N UnifiedGenotyper.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_UnifiedGenotyper.sh \
${sample_name} ${sample_dir} ${sample_temp} 30 10;

## disovery 


#############################################################
## summary metrics ##########################################
#############################################################

#----------------------------------------------------------------------#
# 15. BedTools_DepthOfCoverage 
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N DepthOfCoverage.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_BedTools_DepthOfCoverage.sh \
${sample_name} ${sample_dir} ${sample_temp} ${bed_list};


#----------------------------------------------------------------------#
# 16. CollectMultipleMetrics
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N CollectMultipleMetrics.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_CollectMultipleMetrics.sh \
${sample_name} ${sample_dir} ${sample_temp}


#########################
## END GATK CLEANING   ##
#########################






















