# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root

# Update
RUN apt-get update -y && apt-get install -y python-dev
 
#-----------------------------PLATYPUS-------------------------------
RUN wget -O /tmp/Platypus-latest.tgz http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz \
  && tar xvf /tmp/Platypus-latest.tgz -C /usr/local/pipeline/ \
  && cd /usr/local/pipeline/Platypus_0.7.9.1 \
  && /bin/bash /usr/local/pipeline/Platypus_0.7.9.1/buildPlatypus.sh \
  && chmod -R 755 /usr/local/pipeline/Platypus_0.7.9.1 \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Platypus_0.7.9.1' /home/pipeman/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Platypus_0.7.9.1' ~/.bashrc

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 8080

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
