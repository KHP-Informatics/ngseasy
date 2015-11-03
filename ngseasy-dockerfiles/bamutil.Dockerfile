# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y && apt-get install -y g++ libssl-dev zlib1g-dev

#--------------------------------libStatGen/BamUtil-----------------------------
# Install libStatGen and BamUtil

RUN cd /usr/local/pipeline && \
  git clone https://github.com/statgen/libStatGen.git && \
  cd libStatGen && \
  make all && \
  cd /usr/local/pipeline && \
  git clone https://github.com/statgen/bamUtil.git && \
  cd bamUtil && \
  make cloneLib && \
  make all

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R ngseasy:ngseasy /usr/local/pipeline

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





