# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# + bwa 
RUN wget -O /tmp/bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2 \
    && tar xjvf /tmp/bwa-0.7.10.tar.bz2 -C /usr/local/pipeline/ \
    && cd /usr/local/pipeline/bwa-0.7.10 && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bwa-0.7.10 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bwa-0.7.10' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bwa-0.7.10' ~/.bashrc \
    && echo "alias ngsBWA='/usr/local/pipeline/bwa-0.7.10'" >>  /home/pipeman/.bashrc
    

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
