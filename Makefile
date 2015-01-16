## ngseasy Makefile
## Version 1.0 
## Author: Stephen Newhouse (stephen.j.newhouse@gmail.com)


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
	ngseasybin ngsprojects dockerimages genomes bwaindex bowtie2index stampyindex resources vep snpeff testdata

ngseasybin:
	cp -v $(SRC)/* $(TARGET_BIN)/

ngsprojects: 
	mkdir -v $(INSTALLDIR)/ngs_projects && \
	mkdir -v $(INSTALLDIR)/ngs_projects/raw_fastq && \
	mkdir -v $(INSTALLDIR)/ngs_projects/config_files && \
	mkdir -v $(INSTALLDIR)/ngs_projects/annovardb && \
	mkdir -v $(HOME)/ngseasy_logs

dockerimages:	
	docker pull compbio/ngseasy-base:$(VERSION) && \
	docker pull compbio/ngseasy-fastqc:$(VERSION) && \
	docker pull compbio/ngseasy-trimmomatic:$(VERSION) && \
	docker pull compbio/ngseasy-bwa:$(VERSION) && \
	docker pull compbio/ngseasy-bowtie2:$(VERSION) && \
	docker pull compbio/ngseasy-stampy:$(VERSION) && \
	docker pull compbio/ngseasy-picardtools:$(VERSION) && \
	docker pull compbio/ngseasy-freebayes:$(VERSION) && \
	docker pull compbio/ngseasy-platypus:$(VERSION) && \
	docker pull compbio/ngseasy-delly:$(VERSION) && \
	docker pull compbio/ngseasy-lumpy:$(VERSION) && \
	docker pull compbio/ngseasy-cnmops:$(VERSION) && \
	docker pull compbio/ngseasy-mhmmm:$(VERSION) && \
	docker pull compbio/ngseasy-exomedepth:$(VERSION) && \
	docker pull compbio/ngseasy-bcbiovar:$(VERSION)

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

purgegenomes: 
	rm -rfv $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)

bwaindex:
	cd $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD) && \
	docker run \
	--volume $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/:/home/reference_genomes_b$(GENOMEBUILD) \
	--name bwa_indexing \
	--rm=true \
	-i -t compbio/ngseasy-bwa:$(VERSION) \
	sudo /bin/bash -c "/usr/local/pipeline/bwa-0.7.12/bwa index -a bwtsw /home/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta" && \
	chmod -R 777 $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/*

stampyindex:

bowtie2index:

resources:
	cd $(INSTALLDIR)/ngs_projects && \
	wget https://s3-eu-west-1.amazonaws.com/ngseasy.data/gatk_resources.tar.gz && \
	tar -xvf gatk_resources.tar.gz && \
	chmod -R 766 gatk_resources && \
	rm gatk_resources.tar.gz

testdata:
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
	--volume $(INSTALLDIR)/ngs_projects/reference_genomes_b$(GENOMEBUILD)/:/home/reference_genomes_b$(GENOMEBUILD) \
	--name novoalign_indexing \
	--rm=true \
	-i -t compbio/ngseasy-novoalign:$(VERSION) \
	sudo /bin/bash -c \
    "/usr/local/pipeline/novocraft/novoindex \
    -k 14 \
    -s 3 \
    /home/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).novoindex \
    /home/reference_genomes_b$(GENOMEBUILD)/human_g1k_v$(GENOMEBUILD).fasta" && \
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
	sudo /bin/bash -c \
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




