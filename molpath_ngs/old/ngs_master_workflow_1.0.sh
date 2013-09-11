#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N ngs_master_workflow_v1.0
#$ -M stephen.newhouse@kcl.ac.uk
#$ -m beas
#$ -pe multi_thread 1
#$ -l h_vmem=2G
#$ -p -0.9999
#------------------------------------------------------------------------#


#########################################################################
# -- Authors: Stepgen Newhouse, Amos Folarin, Aditi Gulati              #
# -- Organisation: KCL/SLaM/NHS                                         #
# -- Email: stephen.j.newhouse@gmail.com, amosfolarin@gmail.com         #
# -- Verion: 1.00                                                       #
# -- Date: 20/07/2013                                                   #
# -- DESC: NGS pipeline to perform PE Alignments & GATK cleaning        #
#########################################################################

#------------------------------------------------------------------------#
# Set environmental variables required for child processes (all )
#------------------------------------------------------------------------#

#######
# QUE #
#######

queue_name="short.q,long.q"
email_contact="stephen.newhouse@kcl.ac.uk"

export novo_cpu=8
export novo_mem=3
export sge_h_vmem=8
export java_mem=6


###############
## ngs tools ##
###############

export ngs_picard="/share/apps/picard-tools_1.91/jar"
export ngs_gatk="/share/apps/gatk_2.5-2"
export ngs_novo="/share/apps/novocraft_20130415/bin" ## Novoalign V3.00.03
export ngs_samtools="/share/apps/samtools_0.1.18/bin"
#export PATH=$PATH:${ngs_picard}:${ngs_gatk}:${ngs_novo}:${ngs_samtools}
#ngs="/scratch/project/pipelines/ngs_pipeline_dev/ngs_dev_sjn_tmp/ngs_bin"
#export PATH=$PATH:${ngs}

##################
## pipeline dir ##
##################

export ngs_pipeline="/scratch/project/pipelines/ngs_pipeline_dev/ngs_dev_sjn_tmp"

######################
## reference genomes #
######################

export reference_genome_novoindex="/scratch/data/reference_genomes/gatk_resources/b37/human_g1k_v37.fasta.novoindex"

export reference_genome_seq="/scratch/data/reference_genomes/gatk_resources/b37/human_g1k_v37.fasta"

#######################
## Path to fastq dir ##
#######################

export fastq_dir="/scratch/project/pipelines/ngs_pipeline_dev/aditi_fastq"

###########################
## path to final aln dir ##
###########################

export aln_dir="/scratch/project/pipelines/ngs_pipeline_dev/ngs_molpath_sjn"

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
mPE=${12}

#------------------------------------------------------------------------------#
# END setting environmental variables required for child processes (all )
#------------------------------------------------------------------------------#

# called using "Rscript call_ngs_master_workflow.R <config_file>"
# this calls the following with options read in from R and fired off using R's system() command
# qsub ngs_master_workflow.sh <fastq_prefix> <sample_name> <qual_type> <RGID> <RGLB> <RGPL> <RGPU> <RGSM> <RGCN> <RGDS> <RGDT>
# ALL 8 OPTIONS REQUIRED
# -- RGID=String	Read Group ID Required. <PROJECT_NAME>
# -- RGLB=String	Read Group Library Required. <PATIENT_ID>.<RGPL>.<>
# -- RGPL=String	Read Group platform (e.g. illumina, solid) Required.
# -- RGPU=String	Read Group platform unit (eg. run barcode) Required.
# -- RGSM=String	Read Group sample name Required. PATIENT_ID
# -- RGCN=String	Read Group sequencing center name Required.
# -- RGDS=String	Read Group description Required.
# -- RGDT=Iso8601Date	Read Group run date Required.


##############################
## make sample dir for data ##
##############################

## temp dir
mkdir ${ngstmp}/${sample_name}_temp

## final data dir
mkdir ${aln_dir}/${sample_name}

sample_temp=${ngstmp}/${sample_name}_temp

sample_dir=${aln_dir}/${sample_name}

echo "moving to sample directory "  ${sample_dir}

cd ${sample_dir}

#----------------------------------------------------------------------#
# 0. Check for proper number of command line args and get settings [TO DO]
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# 1. Align PE or SE data 
#----------------------------------------------------------------------#
echo ">>>>>" `date` " :-> " "Aligning PE data "

if [ ${mPE} -eq 1 ]; then

	qsub -q ${queue_name} -N novoalign.${sample_name} -l h_vmem=${novo_mem}G -pe multi_thread ${novo_cpu} -M ${email_contact} -m beas \
	${ngs_pipeline}/ngs_novoalign.${qual_type}.PE.sh \
	${fastq_prefix} \
	${sample_name} \
	${sample_dir};
	
else

	qsub -q ${queue_name} -N novoalign.${sample_name} -l h_vmem=${novo_mem}G -pe multi_thread ${novo_cpu} -M ${email_contact} -m beas \
	${ngs_pipeline}/ngs_novoalign.${qual_type}.SE.sh \
	${fastq_prefix} \
	${sample_name} \
	${sample_dir};
	
fi

#----------------------------------------------------------------------#
# 2. sam2bam
#----------------------------------------------------------------------#

echo ">>>>>" `date` " :-> " "Converting Novoalign SAM to BAM and indexing"

qsub -q ${queue_name} -N sam2bam.${sample_name} -hold_jid novoalign.${sample_name} -l h_vmem=${sge_h_vmem}G -pe multi_thread 1 -M ${email_contact} -m beas \
${ngs_pipeline}/ngs_sam2bam.sh ${sample_name} ${sample_dir};


#----------------------------------------------------------------------#
# 3. SortSam
#----------------------------------------------------------------------#

qsub -q ${queue_name} -N SortSam.${sample_name} -hold_jid sam2bam.${sample_name}  -l h_vmem=${sge_h_vmem}G -pe multi_thread 1 -M ${email_contact} -m beas \
${ngs_pipeline}/ngs_SortSam.sh \
${sample_name} \
${sample_dir} \
${sample_temp};

#----------------------------------------------------------------------#
# 4. AddOrReplaceReadGroups 
#----------------------------------------------------------------------#

echo ">>>>>" `date` " :-> " "Running AddOrReplaceReadGroups" 

qsub -q ${queue_name} -N AddOrReplaceReadGroups.${sample_name} -hold_jid SortSam.${sample_name} -l h_vmem=${sge_h_vmem}G -pe multi_thread 1 -M ${email_contact} -m beas \
${ngs_pipeline}/ngs_AddOrReplaceReadGroups.sh \
${sample_name} \
${sample_dir} \
${sample_temp} \
${mRGID} ${mRGLB} ${mRGPL} ${mRGPU} ${mRGSM} ${mRGCN} ${mRGDS} ${mRGDT};


#----------------------------------------------------------------------#
# 5. MarkDuplicates
#----------------------------------------------------------------------#

echo ">>>>>" `date` " :-> " "Running MarkDuplicates"

qsub -q ${queue_name} -N MarkDuplicates.${sample_name} -hold_jid AddOrReplaceReadGroups.${sample_name} -l h_vmem=${sge_h_vmem}G  -pe multi_thread 1 -M ${email_contact} -m beas \
${ngs_pipeline}/ngs_MarkDuplicates.sh \
${sample_name} \
${sample_dir} \
${sample_temp};

#----------------------------------------------------------------------#
# 6. Clean up sample dir
#----------------------------------------------------------------------#

echo ">>>>>" `date` " :-> " "Running Clean Up 01"

qsub -q ${queue_name} -N rmvIntermediateSAMs.${sample_name} -hold_jid MarkDuplicates.${sample_name} -l h_vmem=1G -M ${email_contact} -m beas \
${ngs_pipeline}/ngs_rmvdIntermediateSAMs.sh \
${sample_name} \
${sample_dir} \
${sample_temp};

#########################
## BEGIN GATK CLEANING ##
#########################

#----------------------------------------------------------------------#
# 6. Clean up
#----------------------------------------------------------------------#


















