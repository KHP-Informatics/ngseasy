#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -p -0.99999999999999999999999999999999999999999999999999
#$ -V

###############################################
## BaseRecalibrator with default covars
################################################

sample_name=${1}
sample_dir=${2}
sample_temp=${3}

cd ${sample_dir}

${java_1_7}/java  -Xmx${gatk_java_mem}g -Djava.io.tmpdir=${sample_temp} -jar ${ngs_gatk}/GenomeAnalysisTK.jar -T BaseRecalibrator -R ${reference_genome_seq}  \
-I ${sample_dir}/${sample_name}.novorealn.bam \
-o ${sample_dir}/${sample_name}.novorealn.recal_data.table \
-knownSites ${b37_1000G_indels} \
-knownSites ${b37_Mills_Devine_2hit_indels} \
-knownSites ${b37_1000G_omni2_5} \
-knownSites ${b37_1000G_snps} \
-knownSites ${b37_hapmap_3_3} \
-knownSites ${b37_dbsnp};
##-nct 1 \
##-log ${sample_dir}/${sample_name}.novorealn.BaseRecalibrator.log;
# indels #
##export b37_1000G_indels="/isilon/irods_a/data_resources/ngs_ref_resources_b37/1000G_phase1.indels.b37.vcf"
## export b37_Mills_Devine_2hit_indels="/isilon/irods_a/data_resources/ngs_ref_resources_b37/Mills_and_1000G_gold_standard.indels.b37.vcf"

# snps #
## export b37_1000G_omni2_5=="/isilon/irods_a/data_resources/ngs_ref_resources_b37/1000G_omni2.5.b37.vcf"
## export b37_dbsnp="/isilon/irods_a/data_resources/ngs_ref_resources_b37/dbsnp_137.b37.vcf"
## export b37_hapmap_3_3="/isilon/irods_a/data_resources/ngs_ref_resources_b37/hapmap_3.3.b37.vcf"
## export b37_1000G_snps="/isilon/irods_a/data_resources/ngs_ref_resources_b37/1000G_phase1.snps.high_confidence.b37.vcf"
###--covariate / -cov ( String[] )
###One or more covariates to be used in the recalibration. 
###Can be specified multiple times. 
##Note that the ReadGroup and QualityScore covariates are required and do not need to be specified. 
##Also, unless --no_standard_covs is #specified, the Cycle and Context covariates are standard and are included by default.
##Use the --list argument to see the available covariates.
