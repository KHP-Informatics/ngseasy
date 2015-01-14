
VERSION=1.0

TARGET_BIN=/bin

SRC=./bin

install:
	cp -v ${SRC}/* ${TARGET_BIN} && \
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

clean:
	rm /bin/ngseas* && rm /bin/ensembl****yaml


