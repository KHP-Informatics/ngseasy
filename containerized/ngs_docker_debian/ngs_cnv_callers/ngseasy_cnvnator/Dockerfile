# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# general purpose tools
RUN yum update
RUN yum -y install cmake gcc-c++ gcc git make python27 python-devel python-yaml ncurses-devel zlib-devel numpy python-pip

# PERL modules for VEP
RUN yum -y install "perl(Archive::Extract)" "perl(CGI)" "perl(DBI)" "perl(Time::HiRes)" "perl(Archive::Tar)" "perl(Archive::Zip)"

# install dependencies
RUN yum -y install libX11-devel libXpm-devel libXft-devel libXext-devel

# get ROOT
RUN curl -OL ftp://root.cern.ch/root/root_v5.34.20.source.tar.gz && \
    tar -xvf root_v5.34.20.source.tar.gz && \
    cd root && \
    ./configure && \
    make && \
    cd .. && \
    sudo mv root /usr/local

RUN source /usr/local/root/bin/thisroot.sh && \
    . ~/.bashrc

# speedseq    
RUN git clone --recursive https://github.com/cc2qe/speedseq && \
    cd speedseq && \
    python speedseq_setup.py && \
    sudo cp -r bin/* /usr/local/bin/ && \ 
    make cnvnator-multi && \
    sudo cp bin/* /usr/local/bin
    
RUN mkdir -p /root/genomes/GRCh37/chroms && \
    cd /root/genomes/GRCh37/chroms && \
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz && \
    wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai && \
    gunzip ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz && \
    cat human_g1k_v37.fasta | awk 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > CHROM".fa" }'

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/    