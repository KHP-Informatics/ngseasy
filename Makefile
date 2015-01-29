## ngseasy Makefile
## Version 1.0 
## Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)

## MEM Min 32GB

## Edit this to reflect version if ya need
VERSION=1.0
GENOMEBUILD=37

## This is where we will make ngs_projects and download metadata to etc etc
## Edit this if you want to install all somewhere else
INSTALLDIR=$(HOME)

## Current working dir
DIR=$(shell pwd)

## Install bin path - edit at will
TARGET_BIN=/bin

## relative path to ngseasy scripts
SRC=./bin

all:
	scripts ngsprojectdir dockerimages genomes annotations bwaindex bowtie2index stampyindex snapindex resources vep snpeff testdata

scripts:
	cp -v $(SRC)/* $(TARGET_BIN)/

ngsprojectdir: 
	mkdir -v $(INSTALLDIR)/ngs_projects && \
	mkdir -v $(INSTALLDIR)/ngs_projects/raw_fastq && \
	mkdir -v $(INSTALLDIR)/ngs_projects/config_files && \
	mkdir -v $(INSTALLDIR)/ngs_projects/annovardb && \
	mkdir -v $(HOME)/ngseasy_logs

purgengsprojectsdir: 
	rm -rfv $(INSTALLDIR)/ngs_projects

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

baseimage:
	docker pull compbio/ngseasy-base:$(VERSION)

fastqc:
	docker pull compbio/ngseasy-fastqc:$(VERSION)

trimmomatic:
	docker pull compbio/ngseasy-trimmomatic:$(VERSION)

bwa:
	docker pull compbio/ngseasy-bwa:$(VERSION)

bowtie2:
	docker pull compbio/ngseasy-bowtie2:$(VERSION)

snap:
	docker pull compbio/ngseasy-snap:$(VERSION)
	
stampy:
	docker pull compbio/ngseasy-stampy:$(VERSION)
	
picardtools:
	docker pull compbio/ngseasy-picardtools:$(VERSION)

freebayes:
	docker pull compbio/ngseasy-freebayes:$(VERSION)

platypus:
	docker pull compbio/ngseasy-platypus:$(VERSION)

delly:
	docker pull compbio/ngseasy-delly:$(VERSION)

lumpy:
	docker pull compbio/ngseasy-lumpy:$(VERSION)

bcbiovar:
	docker pull compbio/ngseasy-bcbiovar:$(VERSION)

cnmops:
	docker pull compbio/ngseasy-cnmops:$(VERSION)

mhmm:
	docker pull compbio/ngseasy-mhmm:$(VERSION)

exomedepth:
	docker pull compbio/ngseasy-exomedepth:$(VERSION)

slope:
	docker pull compbio/ngseasy-slope:$(VERSION)
	
genomes: 
	cd $(INSTALLDIR)/ngs_projects && \
	mkdir reference_genomes_b$(GENOMEBUILD) && \
	cd reference_genomes_b$(GENOMEBUILD) && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta.fai.tar.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta.tar.gz && \
	chmod -R 777 $(INSTALLDIR)/ngs_projects/ && \
	tar -xvf $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta.tar.gz && \
	tar -xvf $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta.fai.tar.gz && \
	rm -v $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta.tar.gz && \
	rm -v $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta.fai.tar.gz

annotations:
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD) && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/LCR_hg19_rmsk.bed.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/SeqCap_EZ_Exome_v3_capture.bed.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/SeqCap_EZ_Exome_v3_primary.bed.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/dbsnp_138.b37.recab.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/human_g1k_v37.chrom_lengths.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/human_g1k_v37_0.5Kwindows.bed.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/human_g1k_v37_1Kwindows.bed.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/nexterarapidcapture_exome_targetedregions_v1.2.bed.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/reference_genomes_b37_annotations/nexterarapidcapture_expandedexome_targetedregions.bed.gz

chromosomes:
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD) && \
	mkdir chroms && \
	cd chroms && \
	awk 'BEGIN { CHROM="" } { if ($$1~"^>") CHROM=substr($$1,2); print $$0 > CHROM".fasta" }' ${INSTALLDIR}/ngs_projects/reference_genomes_b${GENOMEBUILD}/human_g1k_v${GENOMEBUILD}.fasta

purgegenomes: 
	rm -rfv $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)

bwaindex:
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD) && \
	docker run \
	--volume $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/:/home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD) \
	--name bwa_indexing \
	--rm=true \
	-e USER=pipeman \
	--user=pipeman \
	-i -t compbio/ngseasy-bwa:$(VERSION) \
	 /bin/bash -c \
        "/usr/local/pipeline/bwa-0.7.12/bwa index \
        -a bwtsw \
        /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta" && \
	chmod -R 777 $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/*

stampyindex:
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD) && \
	docker run \
	--volume $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/:/home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD) \
	--name stampy_indexing \
	--rm=true \
	-e USER=pipeman \
	--user=pipeman \
	-i -t compbio/ngseasy-stampy:$(VERSION) \
	 /bin/bash -c \
        "python  \
        /usr/local/pipeline/stampy-1.0.23/stampy.py \
        -G /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD) \
        /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta && \
        chmod -R 777 /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/* && \
        python  \
        /usr/local/pipeline/stampy-1.0.23/stampy.py \
        -g /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD) \
        -H /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD) && \
        chmod -R 777 /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/*" && \
	chmod -R 777 $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/*

bowtie2index:
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD) && \
	docker run \
	--volume $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/:/home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD) \
	--name bowtie2_indexing \
	--rm=true \
	-e USER=pipeman \
	--user=pipeman \
	-i -t compbio/ngseasy-bowtie2:$(VERSION) \
	 /bin/bash -c \
         "/usr/local/pipeline/bowtie2-2.2.4/bowtie2-build \
        /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta \
        /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD)" && \
	chmod -R 777 $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/*

snapindex:
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD) && \
	docker run \
	--volume $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/:/home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD) \
	--name snap_indexing \
	--rm=true \
	-e USER=pipeman \
	--user=pipeman \
	-i -t compbio/ngseasy-snap:$(VERSION) \
	 /bin/bash -c \
        "/usr/local/bin/snap index \
        /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta \
        /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/ -hg19" && \
	chmod -R 777 $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/*
	
resources: ngsprojectdir
	cd $(INSTALLDIR)/ngs_projects && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/gatk_resources.tar.gz && \
	tar -xvf gatk_resources.tar.gz && \
	chmod -R 766 gatk_resources && \
	rm gatk_resources.tar.gz

testdata: ngsprojectdir
	cd $(INSTALLDIR)/ngs_projects && \
	mkdir $(INSTALLDIR)/ngs_projects/fastq_test_data && \
	cd $(INSTALLDIR)/ngs_projects/fastq_test_data && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/NA12878s.WEX_1.fq.gz && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/fastq_test_data/NA12878s.WEX_2.fq.gz


gatk:
	cd $(DIR)/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_gatk && \
	docker build --rm=true compbio/ngseasy-gatk:$(VERSION) .

novoalign:
	cd $(DIR)/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_novoalin && \
	docker build --rm=true compbio/ngseasy-novoalign:$(VERSION) .


novoalignindex:	
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD) && \
	docker run \
	--volume $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/:/home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD) \
	--name novoalign_indexing \
	--rm=true \
	-e USER=pipeman \
	--user=pipeman \
	-i -t compbio/ngseasy-novoalign:$(VERSION) \
	 /bin/bash -c \
    "/usr/local/pipeline/novocraft/novoindex \
    -k 13 \
    -s 4 \
    /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).novoindex \
    /home/pipeman/ngs_projects/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta" && \
	chmod -R 777 $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/*

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

clean:
	rm -f -v $(TARGET_BIN)/ngseas* && \
	rm -f -v $(TARGET_BIN)/ensembl****yaml

purgeall:
	rm -f -v $(TARGET_BIN)/ngseas* && \
	rm -f -v $(TARGET_BIN)/ensembl****yaml && \
	docker kill $(shell docker ps -a | awk '(print $1)') && \
	docker rm -f $(shell docker ps -a | awk '(print $1)') && \
	docker rmi -f $(shell docker images -a |  grep ngseasy | awk '(print $3)')

## to do: add options to download and build reference genome builds 
## indexing of genome 












