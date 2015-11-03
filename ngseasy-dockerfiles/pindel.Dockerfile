# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Update
RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install libffi-dev libxml2-dev libxslt1-dev python-dev

## PINDEL #####################################################
## http://gmt.genome.wustl.edu/packages/pindel/quick-start.html
# samtools, htslib, bcftools
RUN cd /usr/local/pipeline && \
    git clone --branch=develop git://github.com/samtools/htslib.git && \
    git clone --branch=develop git://github.com/samtools/bcftools.git && \
    git clone --branch=develop git://github.com/samtools/samtools.git && \
    cd /usr/local/pipeline/htslib && \
    autoconf && \
    ./configure  && \
    make && \
    make install && \
    cd /usr/local/pipeline/bcftools && \
    make && \
    make install 
    ##&& \
    ##cd /usr/local/pipeline/samtools && \
    ##make && \
    ##make install

ADD samtools-0.1.19.tar.bz2 /usr/local/pipeline/

RUN cd /usr/local/pipeline/ && \
	cd samtools-0.1.19 && \
	chmod -R 777 ./* && \
	make

RUN cd /usr/local/pipeline/ && \
	git clone git://github.com/genome/pindel.git && \
	chmod -R 777 ./pindel && \
	cd pindel && \
	/bin/sh ./INSTALL /usr/local/pipeline/samtools-0.1.19 && \
    chmod -R 777 /usr/local/pipeline/pindel && \
    cp -puv Adaptor.pm bam2pindel.pl pindel pindel2vcf sam2pindel /usr/local/bin/

#-------------------------------PERMISSIONS-------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R ngseasy:ngseasy /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/USER ngseasy
WORKDIR /home/ngseasy
