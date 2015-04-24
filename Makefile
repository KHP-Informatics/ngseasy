## ngseasy Makefile
## Version 1.0 
## Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)

## MEM Min 32GB

## Edit this to reflect version if ya need
VERSION=1.0
WHOAMI=$(shell whoami)

## This is where we will make ngs_projects and download metadata to etc etc
## Edit this if you want to install all somewhere else
# eg:
# make INSTALLDIR="/medida/scratch" all
#
INSTALLDIR=$(HOME)

## Current working dir
DIR=$(shell pwd)

## Install bin path - edit at will
TARGET_BIN=/usr/local/bin

## relative path to ngseasy scripts
SRC=./bin

## Basic install - no annotation data bases or manual build tools
all: scripts ngsprojectdir dockerimages b37 hg19 testdata

## install scripts to target bin
scripts:
	chmod 775 $(SRC)/*
	cp -v $(SRC)/* $(TARGET_BIN)/

rmscripts:
	rm -fv $(TARGET_BIN)/ngseasy* && rm -fv $(TARGET_BIN)/ngseasy 

updatescripts:
	git pull && rm -fv $(TARGET_BIN)/ngseasy* && rm -fv $(TARGET_BIN)/ngseasy && chmod 777 $(SRC)/ && cp -v $(SRC)/* $(TARGET_BIN)/

## Make Top level project directories
ngsprojectdir: 
	mkdir -v -p $(INSTALLDIR)/ngs_projects && \
	mkdir -v -p $(INSTALLDIR)/ngs_projects/raw_fastq && \
	mkdir -v -p $(INSTALLDIR)/ngs_projects/config_files && \
	mkdir -v -p $(INSTALLDIR)/ngs_projects/run_logs && \
	mkdir -v -p $(HOME)/ngseasy_logs && \
	mkdir -v -p $(HOME)/ngseasy_tmp

purgengsprojectsdir: 
	rm -rfv $(INSTALLDIR)/ngs_projects

## Get all docker images 
dockerimages:	
	docker pull compbio/ngseasy-base:$(VERSION) && \
	docker pull compbio/ngseasy-fastqc:$(VERSION) && \
	docker pull compbio/ngseasy-trimmomatic:$(VERSION) && \
	docker pull compbio/ngseasy-snap:$(VERSION) && \
	docker pull compbio/ngseasy-bwa:$(VERSION) && \
	docker pull compbio/ngseasy-bowtie2:$(VERSION) && \
	docker pull compbio/ngseasy-stampy:$(VERSION) && \
	docker pull compbio/ngseasy-picardtools:$(VERSION) && \
	docker pull compbio/ngseasy-freebayes:$(VERSION) && \
	docker pull compbio/ngseasy-platypus:$(VERSION) && \
	docker pull compbio/ngseasy-delly:$(VERSION) && \
	docker pull compbio/ngseasy-lumpy:$(VERSION) && \
	docker pull compbio/ngseasy-bcbiovar:$(VERSION) && \
	docker pull compbio/ngseasy-cnmops:$(VERSION) && \
	docker pull compbio/ngseasy-mhmm:$(VERSION) && \
	docker pull compbio/ngseasy-exomedepth:$(VERSION) && \
	docker pull compbio/ngseasy-slope:$(VERSION)

indexgenomes: dockerimages
	bwaindex bowtie2index stampyindex snapindex

baseimage:
	docker pull compbio/ngseasy-base:$(VERSION)

fastqc: baseimage
	docker pull compbio/ngseasy-fastqc:$(VERSION)

trimmomatic: baseimage
	docker pull compbio/ngseasy-trimmomatic:$(VERSION)

bwa: baseimage
	docker pull compbio/ngseasy-bwa:$(VERSION)

bowtie2: baseimage
	docker pull compbio/ngseasy-bowtie2:$(VERSION)

snap: baseimage
	docker pull compbio/ngseasy-snap:$(VERSION)
	
stampy: bwa baseimage
	docker pull compbio/ngseasy-stampy:$(VERSION)
	
picardtools: baseimage
	docker pull compbio/ngseasy-picardtools:$(VERSION)

freebayes: baseimage
	docker pull compbio/ngseasy-freebayes:$(VERSION)

platypus: baseimage
	docker pull compbio/ngseasy-platypus:$(VERSION)

delly: baseimage
	docker pull compbio/ngseasy-delly:$(VERSION)

lumpy: baseimage
	docker pull compbio/ngseasy-lumpy:$(VERSION)

bcbiovar:
	docker pull compbio/ngseasy-bcbiovar:$(VERSION)

cnmops: baseimage
	docker pull compbio/ngseasy-cnmops:$(VERSION)

mhmm: baseimage
	docker pull compbio/ngseasy-mhmm:$(VERSION)

exomedepth: baseimage
	docker pull compbio/ngseasy-exomedepth:$(VERSION)

slope: baseimage
	docker pull compbio/ngseasy-slope:$(VERSION)

# b37 Genomes idexed and resources	
b37: 
	cd $(INSTALLDIR)/ngs_projects && \
	mkdir -p reference_genomes_b37 && \
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b37 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/1000G_omni2.5.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/1000G_omni2.5.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/1000G_phase1.indels.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/1000G_phase1.indels.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/1000G_phase1.snps.high_confidence.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/1000G_phase1.snps.high_confidence.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/CEUTrio.HiSeq.WGS.b37.NA12878.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/CEUTrio.HiSeq.WGS.b37.NA12878.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/CEUTrio.HiSeq.WGS.b37.bestPractices.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/CEUTrio.HiSeq.WGS.b37.bestPractices.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/Genome && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/GenomeIndex && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/GenomeIndexHash && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/LCR_hg19_rmsk.bed && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/Mills_and_1000G_gold_standard.indels.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/Mills_and_1000G_gold_standard.indels.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.sites.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.sites.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/NA12878.knowledgebase.snapshot.20131119.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/NA12878.knowledgebase.snapshot.20131119.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/OverflowTable && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/S03723314_Regions.bed && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/S06588914_Regions_trimmed.bed && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/SeqCap_EZ_Exome_v3_capture.bed && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/SeqCap_EZ_Exome_v3_primary.bed && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/b37.genome && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/contaminant_list.fa && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/dbsnp_138.b37.excluding_sites_after_129.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/dbsnp_138.b37.excluding_sites_after_129.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/dbsnp_138.b37.recab && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/dbsnp_138.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/dbsnp_138.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/hapmap_3.3.b37.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/hapmap_3.3.b37.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.1.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.2.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.3.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.4.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.chrom_lengths && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta.amb && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta.ann && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta.bwt && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta.fai && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta.pac && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.fasta.sa && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.novoindex && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.rev.1.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.rev.2.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.sthash && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37.stidx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37_0.5Kwindows.bed && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/human_g1k_v37_1Kwindows.bed && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/nexterarapidcapture_exome_targetedregions_v1.2.bed && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37/nexterarapidcapture_expandedexome_targetedregions.bed && \
	chmod -R 775 $(INSTALLDIR)/ngs_projects/reference_genomes_b37/

# hg19 Genomes idexed and resources	
hg19: 
	cd $(INSTALLDIR)/ngs_projects && \
	mkdir -p reference_genomes_hg19 && \
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_hg19 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/1000G_omni2.5.hg19.sites.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/1000G_omni2.5.hg19.sites.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/1000G_phase1.indels.hg19.sites.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/1000G_phase1.indels.hg19.sites.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/Genome && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/GenomeIndex && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/GenomeIndexHash && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz.tbi && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.sites.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.sites.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/NA12878.knowledgebase.snapshot.20131119.hg19.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/NA12878.knowledgebase.snapshot.20131119.hg19.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/OverflowTable && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/contaminant_list.fa && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/dbsnp_138.hg19.recab && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/dbsnp_138.hg19.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/dbsnp_138.hg19.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/get_hg19.sh && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/get_hg19_others.sh && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/hapmap_3.3.hg19.sites.vcf && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/hapmap_3.3.hg19.sites.vcf.idx && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/hg19.genome && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/index_bowtie.sh && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/index_bwa.sh && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/index_novo.sh && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/index_snap.sh && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/index_stampy.sh && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19-bs.umfa && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.1.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.2.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.3.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.4.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.dict && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fai && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.amb && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.ann && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.bwt && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.fai && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.fai.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.novoindex && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.pac && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.fasta.sa && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.rev.1.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.rev.2.bt2 && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.sthash && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_hg19/ucsc.hg19.stidx && \
	chmod -R 775 $(INSTALLDIR)/ngs_projects/reference_genomes_hg19/

##  Test data and stick it in raw_fastq
testdata: ngsprojectdir
	cd $(INSTALLDIR)/ngs_projects/raw_fastq && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/ && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/NA12878D_HiSeqX_R1.fastq.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/NA12878D_HiSeqX_R2.fastq.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/illumina.100bp.pe.wex.150x_1.fastq.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/illumina.100bp.pe.wex.150x_2.fastq.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/illumina.100bp.pe.wex.30x_1.fastq.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/illumina.100bp.pe.wex.30x_2.fastq.gz && \
	chmod -R 775 $(INSTALLDIR)/ngs_projects/raw_fastq/ && \
	cd $(INSTALLDIR)/ngs_projects/config_files && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/ngseasy_test.config.tsv && \
	chmod -R 775 $(INSTALLDIR)/ngs_projects/config_files/ 

## Manual Builds
gatk:
	cd $(DIR)/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_gatk && \
	docker build --rm=true compbio/ngseasy-gatk:$(VERSION) .

novoalign:
	cd $(DIR)/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_novoalin && \
	docker build --rm=true compbio/ngseasy-novoalign:$(VERSION) .

vep:
	cd $(DIR)/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_vep && \
	docker build --rm=true compbio/ngseasy-vep:$(VERSION) .

snpeff:
	cd $(DIR)/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_snpeff && \
	docker build --rm=true compbio/ngseasy-snpeff:$(VERSION) .

annovar:
	cd $(DIR)/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_annovar && \
	docker build --rm=true compbio/ngseasy-annovar:$(VERSION) .

annovardb:	
	docker run \
	--volume $(INSTALLDIR)/ngs_projects/annovardb/:/home/annovardb \
	--name get_annovardb \
	--rm=true \
	-i -t compbio/ngseasy-annovar:$(VERSION) \
	 /bin/bash -c \
	"/bin/bash /usr/local/pipeline/annovar/get_annovar_gene_databases.sh && /bin/bash /usr/local/pipeline/annovar/get_annovar_databases.sh"

## Cleanups
clean:
	rm -f -v $(TARGET_BIN)/ngseas* && \
	rm -f -v $(TARGET_BIN)/ensembl****yaml

purgeall:
	rm -f -v $(TARGET_BIN)/ngseas* && \
	rm -f -v $(TARGET_BIN)/ensembl****yaml && \
	docker kill $(shell docker ps -a | awk '(print $$1)') && \
	docker rm -f $(shell docker ps -a | awk '(print $$1)') && \
	docker rmi -f $(shell docker images -a |  grep ngseasy | awk '(print $$3)') && \

purgegenomes: 
	rm -rfv $(INSTALLDIR)/ngs_projects/reference_genomes_b37 && \
	rm -rfv $(INSTALLDIR)/ngs_projects/reference_genomes_hg19

############################################################################	
## Make sep chroms
#chrb37:
#	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b37 && \
#	mkdir chroms && \
#	cd chroms && \
#	awk 'BEGIN { CHROM="" } { if ($$1~"^>") CHROM=substr($$1,2); print $$0 > CHROM".fasta" }' ${INSTALLDIR}/ngs_projects/reference_genomes_b37/human_g1k_v37.fasta

#chrhg19:
#	cd $(INSTALLDIR)/ngs_projects/reference_genomes_hg19 && \
#	mkdir chroms && \
#	cd chroms && \
#	awk 'BEGIN { CHROM="" } { if ($$1~"^>") CHROM=substr($$1,2); print $$0 > CHROM".fasta" }' ${INSTALLDIR}/ngs_projects/reference_genomes_hg19/ucsc.hg19.fasta










