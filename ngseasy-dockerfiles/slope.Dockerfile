# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# general purpose tools
RUN apt-get update -y && apt-get upgrade -y

RUN cd /usr/local/pipeline && \
    wget http://www-genepi.med.utah.edu/suppl/SLOPE/slope_align && \
    wget http://www-genepi.med.utah.edu/suppl/SLOPE/slope_align.cpp && \
    wget http://www-genepi.med.utah.edu/suppl/SLOPE/slope_cluster && \
    wget http://www-genepi.med.utah.edu/suppl/SLOPE/slope_cluster.cpp && \
    chmod 777 ./* && \
    cp -v /usr/local/pipeline/slope* /usr/local/bin

#-------------------------------PERMISSIONS--------------------------
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
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/    
USER ngseasy
WORKDIR /home/ngseasy
