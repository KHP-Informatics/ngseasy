FROM compbio/ngseasy-base:1.0

# Maintainer from Tobias Rausch rausch@embl.de
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# update package repository
RUN apt-get update

# install g++, git, zlib, cmake, boost, ...
RUN apt-get install -y build-essential g++ git cmake zlib1g-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev

# set environment
ENV BOOST_ROOT /usr

# install delly
RUN cd /usr/local/pipeline && \
    wget https://github.com/tobiasrausch/delly/releases/download/v0.6.1/delly_v0.6.1_CentOS5.4_x86_64bit && \
    wget https://github.com/tobiasrausch/delly/releases/download/v0.6.1/delly_v0.6.1_parallel_linux_x86_64bit && \
    mv -v delly_v0.6.1_CentOS5.4_x86_64bit delly && \
    mv -v delly_v0.6.1_parallel_linux_x86_64bit delly-parallel && \
    chmod -R 777 /usr/local/pipeline && \
    cp -v delly /usr/local/bin && \
    cp -v delly-parallel /usr/local/bin

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





