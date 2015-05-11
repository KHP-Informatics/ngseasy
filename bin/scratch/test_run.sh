


#Good Run on Google Cloud

ngseasy_initiate_project -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects

ngseasy_volumes_container -d /media/container-vol/ngs_projects

ngseasy_initiate_fastq -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects

# get annovar databases
sudo docker run \
	--rm=true \
	--name annovar_db_get \
	--volumes-from volumes_container \
	-i \
	-t \
	compbio/ngseasy-annovar:v0.9.2 \
	/bin/bash \
	/usr/local/pipeline/annovar/get_annovar_gene_databases.sh

sudo docker run \
	--rm=true \
	--volumes-from volumes_container \
	-i \
	-t \
	compbio/ngseasy-annovar:v0.9.2 \
	/bin/bash \
	/usr/local/pipeline/annovar/get_annovar_databases.sh

# PIPELINE SET TO ngseasy_fastq in config file
# this creates new sample confog file from project confif file and feeds as input to pipeline ngseasy_fastq
ngseasy -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects


#--------------Step through------------------------------------#

# 1. copy fastq to project/sample/fastq
# sep call to ngseasy_fastqc
ngseasy_fastqc -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 

# 2 Trimm
# sep call to ngseasy_trimmomatic
ngseasy_trimmomatic -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 

# 3 Align
# run aligment 
ngseasy_alignment -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 

# 4 Add read groups
# Add Read Groups
ngseasy_addreadgroup -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 

# 5 Mark Duplicates
#Mark Dupes
ngseasy_markduplicates -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 

# 6 Realign Around Indels
#GATK Indel Realn
ngseasy_indel_realn -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 

# 7 Base recalibration
#GATK Base Recal
ngseasy_base_recal -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 

# 8 Filter BAM to good paried reads only
ngseasy_filter_recalbam -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 


# 9 Call SNPs and Small Indels (ensembl)

ngseasy_variant_calling_fast_ensemble -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 

ngseasy_variant_calling -c /media/container-vol/ngs_projects/config_files/example.config.tsv -d /media/container-vol/ngs_projects 


#Alignment qc
ngseasy_alignment_qc


#
ngseasy_sv_call




*****


http://basecallbio.wordpress.com/2013/04/23/base-quality-score-rebinning/
