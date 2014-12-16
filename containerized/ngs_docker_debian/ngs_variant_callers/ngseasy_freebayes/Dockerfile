# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

#---------------------- FREEBAYES  -------------------------------- 
RUN cd /usr/local/pipeline \
  && git clone --recursive git://github.com/ekg/freebayes.git \
  && cd /usr/local/pipeline/freebayes \
  && make \
  && chmod -R 755 /usr/local/pipeline/freebayes \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/freebayes \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/freebayes/bin' /home/pipeman/.bashrc \
  && echo "alias ngsFreebayes='/usr/local/pipeline/freebayes/bin'" >>  /home/pipeman/.bashrc


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





