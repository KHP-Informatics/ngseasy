# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

#---------------------- FREEBAYES  -------------------------------- 
RUN cd /usr/local/pipeline \
  && git clone --recursive git://github.com/ekg/freebayes.git \
  && cd /usr/local/pipeline/freebayes \
  && make \
  && chmod -R 777 /usr/local/pipeline/freebayes \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/freebayes/bin' /home/ngseasy/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/freebayes/bin' ~/.bashrc \
  && cp -v /usr/local/pipeline/freebayes/bin/* /usr/local/bin

ADD fix_ambiguous /usr/local/bin/

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





