# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.-------------------------------------
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update ---------------------------------------------------------------
RUN apt-get update -y && apt-get upgrade -y

# get ROOT & speedseq/cnvnator-multi-------------------------------------------------------------
RUN cd /usr/local/pipeline && \
  curl -OL ftp://root.cern.ch/root/root_v5.34.20.source.tar.gz && \
  tar -xvf root_v5.34.20.source.tar.gz && \
  cd root && \
  ./configure && \
  make && \
  cd .. && \
  sudo mv root /usr/local  && \
  cat "source /usr/local/root/bin/thisroot.sh" >> .~/bashrc && \
  . ~/.bashrc && \
  cd /usr/local/pipeline && \
  git clone --recursive https://github.com/cc2qe/speedseq && \
  cd /usr/local/pipeline/speedseq/src/ && \
  sed -i s/ROOTFLAGS = -m64 -O3/#ROOTFLAGS = -m64 -O3/g Makefile && \
  sed -i s/#ROOTFLAGS = -pthread -m64/ROOTFLAGS = -pthread -m64/g Makefile && \ 
  make cnvnator-multi && \
  sudo cp bin/* /usr/local/bin

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 80

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/

