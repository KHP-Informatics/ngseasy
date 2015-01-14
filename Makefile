## Edit this to reflect version if ya need
VERSION=1.0

## Install bin path - edit at will
TARGET_BIN=/bin

## relative path to ngseasy scripts
SRC=./bin

all:
	install dockerimages genomes resources testdata

install:
	cp -v ${SRC}/* ${TARGET_BIN}

dockerimages:	
	sudo docker pull compbio/ngseasy-base:${VERSION} && \
	sudo docker pull compbio/ngseasy-fastqc:${VERSION} && \
	sudo docker pull compbio/ngseasy-trimmomatic:${VERSION} && \
	sudo docker pull compbio/ngseasy-bwa:${VERSION} && \
	sudo docker pull compbio/ngseasy-bowtie2:${VERSION} && \
	sudo docker pull compbio/ngseasy-stampy:${VERSION} && \
	sudo docker pull compbio/ngseasy-picardtools:${VERSION} && \
	sudo docker pull compbio/ngseasy-freebayes:${VERSION} && \
	sudo docker pull compbio/ngseasy-platypus:${VERSION} && \
	sudo docker pull compbio/ngseasy-delly:${VERSION} && \
	sudo docker pull compbio/ngseasy-lumpy:${VERSION} && \
	sudo docker pull compbio/ngseasy-cnmops:${VERSION} && \
	sudo docker pull compbio/ngseasy-mhmmm:${VERSION} && \
	sudo docker pull compbio/ngseasy-exomedepth:${VERSION} && \
	sudo docker pull compbio/ngseasy-bcbiovar:${VERSION} && \
	sudo docker pull compbio/ngseasy-snpeff:${VERSION}

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



