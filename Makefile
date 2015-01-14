## Edit this to reflect version if ya need
VERSION=1.0

## Dowload dir
DIR=`pwd`

## Install bin path - edit at will
TARGET_BIN=/bin

## relative path to ngseasy scripts
SRC=./bin

all:
	install dockerimages genomes resources testdata vep snpeff

install:
	cp -v ${SRC}/* ${TARGET_BIN}

dockerimages:	
	docker pull compbio/ngseasy-base:${VERSION} && \
	docker pull compbio/ngseasy-fastqc:${VERSION} && \
	docker pull compbio/ngseasy-trimmomatic:${VERSION} && \
	docker pull compbio/ngseasy-bwa:${VERSION} && \
	docker pull compbio/ngseasy-bowtie2:${VERSION} && \
	docker pull compbio/ngseasy-stampy:${VERSION} && \
	docker pull compbio/ngseasy-picardtools:${VERSION} && \
	docker pull compbio/ngseasy-freebayes:${VERSION} && \
	docker pull compbio/ngseasy-platypus:${VERSION} && \
	docker pull compbio/ngseasy-delly:${VERSION} && \
	docker pull compbio/ngseasy-lumpy:${VERSION} && \
	docker pull compbio/ngseasy-cnmops:${VERSION} && \
	docker pull compbio/ngseasy-mhmmm:${VERSION} && \
	docker pull compbio/ngseasy-exomedepth:${VERSION} && \
	docker pull compbio/ngseasy-bcbiovar:${VERSION}

vep:
	cd ${DIR}/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_vep && \
	docker build --rm=true compbio/ngseasy-vep:${VERSION} 

snpeff:
	cd ${DIR}/containerized/ngs_docker_debian/ngs_variant_annotators/ngseasy_snpeff && \
	docker build --rm=true compbio/ngseasy-snpeff:${VERSION} 

genomes:
	wget http://s3.amazonaws.com/ && \
	tar -xvf 

resources:
	wget http://s3.amazonaws.com/ && \
	tar -xvf

testdata:
	wget http://s3.amazonaws.com/ && \
	tar -xvf


clean:
	rm -fv ${TARGET_BIN}/ngseas* && \
	rm -fv ${TARGET_BIN}/ensembl****yaml

purge:
	rm -f -v ${TARGET_BIN}/ngseas* && \
	rm -f -v ${TARGET_BIN}/ensembl****yaml && \
	docker kill $(docker ps -a | awk '{print $1}') && \
	docker rm -f $(docker ps -a | awk '{print $1}') && \
	docker rmi -f $(docker images -a | awk '{print $3}')




