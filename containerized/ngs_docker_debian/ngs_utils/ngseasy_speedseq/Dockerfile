# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
#ENV HOME /home/root
ENV DEBIAN_FRONTEND noninteractive

# general purpose tools
RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y python python-dev python-yaml ncurses-dev python-numpy python-pip

# install dependencies
RUN apt-get install -y libX11-dev libXpm-dev libXft-dev libXext-dev

# make .bashrc 
RUN touch ${HOME}/.bashrc 

# get ROOT
RUN cd /tmp && \
    curl -OL ftp://root.cern.ch/root/root_v5.34.20.source.tar.gz && \
    tar -xvf root_v5.34.20.source.tar.gz && \
    chmod -R 777 ./root && \
    cd ./root && \
    ./configure && \
    make -j30 && \
    cd .. && \
    mv root /usr/local

#RUN echo "source /usr/local/root/bin/thisroot.sh" >> ${HOME}/.bashrc  && \
#    source ${HOME}/.bashrc

   
# speedseq-------------------------------------------------------------
ENV ROOTSYS /usr/local/root

RUN cd /usr/local/pipeline && \
  git clone --recursive https://github.com/cc2qe/speedseq && \
  chmod -R 777 ./speedseq && \
  cd /usr/local/pipeline/speedseq/src/ && \
  sed -i 's/ROOTFLAGS = -m64 -O3/#ROOTFLAGS = -m64 -O3/g' Makefile && \
  sed -i 's/#ROOTFLAGS = -pthread -m64/ROOTFLAGS = -pthread -m64/g' Makefile && \ 
  cd /usr/local/pipeline/speedseq && \
  make

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