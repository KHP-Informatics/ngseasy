#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M aditi.gulati@nhs.net
#$ -m beas
#$ -pe multi_thread 3
#$ -l h_vmem=8G
#$ -p -0.99999999999999999999999999999999999999999999999999999999999999999


###############################################
## UnifiedGenotyper ############################
###############################################

##sample_name=${1}
##sample_dir=${2}
##sample_temp=${3}

## UnifiedGenotyper option 
##stand_call_conf=${4} ## 40.0 ; default 30
##stand_emit_conf=${5} ## 10.0 ; default 30

##cd ${sample_dir}

/share/java/jdk1.7.0/bin/java -Xmx8g -jar /share/apps/gatk_2.7-2/GenomeAnalysisTK.jar -T VariantsToTable -R /isilon/irods_a/datasets_res/Vault/ngs_ref_resources_b37/human_g1k_v37.fasta -V Sample_barneylib3samp8.Bait_Capture_LiverGenes_Panel.illumina/Sample_barneylib3samp8.Bait_Capture_LiverGenes_Panel.illumina.novorecal.UnifiedGenotyper.raw.snps.indels.recode.vcf --allowMissingData -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -GF GT -GF AD -GF DP -GF GQ -o Sample_barneylib3samp8.Bait_Capture_LiverGenes_Panel.illumina/Sample_barneylib3samp8.Bait_Capture_LiverGenes_Panel.illumina.novorecal.UnifiedGenotyper.raw.snps.indels.recode.table;

###ngs_VariantsToTable.sh###


##########mv -v UnifiedGenotyper.${sample_name}.* ${sample_dir}/sge_out/
