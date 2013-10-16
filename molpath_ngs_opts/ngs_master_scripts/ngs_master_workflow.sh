#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
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

# called by the main dispatch: 
#   Rscript call_ngs_master_workflow.R <patient.config> <pipeline_env.config> 

# For instructions on running the pipeline see the README file in this this directory:




#------------------------------------------------------------------------#
# Set environmental variables required for child processes (all )
#------------------------------------------------------------------------#

##################
## pipeline dir ##
##################
export ngs_pipeline=${18} #"/home/snewhousebrc/scratch/pipelines/ngs/molpath_ngs_dev" 

#######
# QUE #
#######
export queue_name=${19} #"short.q,long.q"

###############
## ngs tools ##
###############
export ngs_picard=${20} #"/share/apps/picard-tools_1.91/jar"
export ngs_gatk=${21} #"/share/apps/gatk_2.7-2"  ## gatk_2.7-2 needs java 1.7
export ngs_novo=${22} #"/share/apps/novocraft_current/bin/" ## Novoalign V3.01.01
export ngs_samtools=${23} #"/share/apps/samtools_0.1.18/bin"

######################
## reference genomes #
## movef refs to isilon and all b37 now!
######################
export reference_genome_novoindex=${24} #"/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/human_g1k_v37.fasta.novoindex"
export reference_genome_seq=${25} #"/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/human_g1k_v37.fasta"

############################
## ref vcf files for gatk ##
############################
# indels #
export b37_1000G_indels=${26} #"/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/1000G_phase1.indels.b37.vcf"
export b37_Mills_Devine_2hit_indels=${27} #"/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
# snps #
export b37_1000G_omni2_5=${28} #"/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/1000G_omni2.5.b37.vcf"
export b37_dbsnp=${29} #"/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/dbsnp_137.b37.vcf"
export b37_hapmap_3_3=${30} #"/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/hapmap_3.3.b37.vcf"
export b37_1000G_snps=${31} #"/isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/1000G_phase1.snps.high_confidence.b37.vcf"

#############
## annovar ##
#############
export annovar=${32} #"/isilon/irods_a/datasets_res/Vault/annovar_2013Aug23"
export annovar_humandb=${33} #"/isilon/irods_a/datasets_res/Vault/annovar_2013Aug23/humandb"

#######################
## JAVA 1.7 for GATK ##
#######################
export java_v1_7=${34} #"/share/java/jdk1.7.0/bin"

#####################
# mem and cpu vars ##
#####################

## Novoalign ##
export novo_cpu=${35} #8
export novo_mem=${36} #3

## Java & Picardtools ##
export sge_h_vmem=${37} #8
export java_mem=${38} #6

## Java & GATK ##
export gatk_h_vmem=${39} #8
export gatk_java_mem=${40} #6

###########################
## path to tmp aln dir   ##
###########################
export ngstmp=${41}  #"/home/snewhousebrc/scratch/ngs_temp"



#---------------------------------------------------------------------------------------------------------------------------------------------#
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
bed_type=${14}  # region or whole_genome : needed for coverage


############################
## email contact for qsub ##
############################
export email_contact=${15}

#######################
## Path to fastq dir ##
#######################
export fastq_dir=${16}

###########################
## path to final aln dir ##
###########################
export aln_dir=${17}


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

#########################################
## set sample temp and output dirs ######
#########################################

sample_temp=${ngstmp}/${sample_name}_temp

sample_dir=${aln_dir}/${sample_name}

#------------------------------------------------------------------------------#

######################################
# moving to sample directory #########
######################################

echo "moving to sample directory "  ${sample_dir}

cd ${sample_dir}

##############################
## START ALIGNMENT PIPELINE ##
##############################


#----------------------------------------------------------------------#
# 1. Align PE or SE data 
#----------------------------------------------------------------------#
##echo ">>>>>" `date` " :-> " "Aligning PE data "

if [ ${mPE} -eq 1 ] && [ ${mRGPL} == "illumina" ]; then

echo "Reads are PE ILLUMINA"

qsub -q ${queue_name} -N novoalign.${sample_name} -l h_vmem=${novo_mem}G -pe multi_thread ${novo_cpu} -M ${email_contact} -m beas ${ngs_pipeline}/ngs_novoalign.illumina.${qual_type}.PE.sh ${fastq_prefix} ${sample_name} ${sample_dir};

elif [ ${mPE} -eq 0 ] && [ ${mRGPL} == "illumina" ]; then

echo "Reads are SE ILLUMINA"

qsub -q ${queue_name} -N novoalign.${sample_name} -l h_vmem=${novo_mem}G -pe multi_thread ${novo_cpu} -M ${email_contact} -m beas ${ngs_pipeline}/ngs_novoalign.illumina.${qual_type}.SE.sh ${fastq_prefix} ${sample_name} ${sample_dir};

elif [ ${mPE} -eq 0 ] && [ ${mRGPL} == "IONTORRENT" ]; then

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

############ qsub -q ${queue_name} -N BaseRecalibrator_after.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_BaseRecalibrator_after.sh ${sample_name} ${sample_dir} ${sample_temp};


#----------------------------------------------------------------------#
# 12. AnalyzeCovariates before & after recal
#----------------------------------------------------------------------#
##echo ">>>>>" `date` " :-> " "Running AnalyzeCovariates" 
#### qsub -q ${queue_name} -N AnalyzeCovariates_before_and_after_BQSR.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_AnalyzeCovariates_before_and_after_BQSR.sh \
#### ${sample_name} ${sample_dir} ${sample_temp};

##############################################
##  CALL VARIANTS SINGLE SAMPLE ##############
##############################################

#----------------------------------------------------------------------#
# 13. HaplotypeCaller 
#----------------------------------------------------------------------#
i
qsub -q ${queue_name} -N HaplotypeCaller.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_HaplotypeCaller.sh \
${sample_name} ${sample_dir} ${sample_temp} 30 10;

#----------------------------------------------------------------------#
# 14. UnifiedGenotyper 
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N UnifiedGenotyper.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_UnifiedGenotyper.sh \
${sample_name} ${sample_dir} ${sample_temp} 30 10;

## to do : discovery > list > merge novels with knowns and re-genotype. Include a genotype all bases option


#############################################################
## summary metrics ##########################################
#############################################################

#----------------------------------------------------------------------#
# 15. BedTools_DepthOfCoverage 
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N DepthOfCoverage.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_BedTools_DepthOfCoverage.sh \
${sample_name} ${sample_dir} ${sample_temp} ${bed_list} ${bed_type};


#----------------------------------------------------------------------#
# 16. CollectMultipleMetrics
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N CollectMultipleMetrics.${sample_name} -hold_jid PrintReads_BQSR.${sample_name} -l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_CollectMultipleMetrics.sh \
${sample_name} ${sample_dir} ${sample_temp}



#######################################################################
## Variant annotations ################################################
#######################################################################

#----------------------------------------------------------------------#
# 17. Table Annovar
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N annovar_UnifiedGenotyper.${sample_name} -hold_jid UnifiedGenotyper.${sample_name} \
-l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_table_annovar_UnifiedGenotyper_hg19.sh ${sample_name} ${sample_dir} ${sample_temp} "UnifiedGenotyper";

qsub -q ${queue_name} -N annovar_HaplotypeCaller.${sample_name} -hold_jid HaplotypeCaller.${sample_name} \
-l h_vmem=${gatk_h_vmem}G -M ${email_contact} -m beas ${ngs_pipeline}/ngs_table_annovar_HaplotypeCaller_hg19.sh ${sample_name} ${sample_dir} ${sample_temp} "HaplotypeCaller";


#########################
## END ##################
#########################






















