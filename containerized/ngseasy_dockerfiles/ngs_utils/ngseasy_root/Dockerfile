################################################################################
# Base image
################################################################################
FROM compbio/ngseasy-base:1.0-r001

################################################################################
# Maintainer
################################################################################
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

################################################################################
# Set correct environment variables.
################################################################################
##ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

################################################################################
# general purpose tools
################################################################################
RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y python python-dev python-yaml ncurses-dev python-numpy python-pip
RUN pip install scipy pysam

################################################################################
# install dependencies
################################################################################
RUN apt-get install -y libX11-dev libXpm-dev libXft-dev libXext-dev

################################################################################
# get ROOT
################################################################################
RUN cd /tmp && \
    curl -OL ftp://root.cern.ch/root/root_v5.34.20.source.tar.gz && \
    tar -xvf root_v5.34.20.source.tar.gz && \
    cd root && \
    ./configure && \
    make && \
    cd .. && \
    mv root /usr/local

RUN  echo "source /usr/local/root/bin/thisroot.sh" >>  /root/.bashrc  && \
    /bin/bash -c ". /root/.bashrc"

################################################################################
# speedseq
################################################################################
ENV ROOTSYS /usr/local/root

RUN cd /usr/local/pipeline && \
  git clone --recursive https://github.com/cc2qe/speedseq && \
  chmod -R 777 /usr/local/pipeline/speedseq && \
  cd /usr/local/pipeline/speedseq && \
  make ROOTFLAGS=" -pthread -m64"

################################################################################
# PERMISSIONS
################################################################################
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

################################################################################
#Cleanup the temp dir
################################################################################
RUN rm -rf /tmp/*

################################################################################
#open ports private only
################################################################################
EXPOSE 8080

################################################################################
# Use baseimage-docker's bash.
################################################################################
CMD ["/bin/bash"]

################################################################################
#Clean up APT when done.
################################################################################
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
