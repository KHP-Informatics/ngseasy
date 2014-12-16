# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# + novoalign (registration required,  get compressed file and put in context dir for the build)
# + get novoalign.lic from novoalign ~ $1000 and put in context dir for the build
# + get updated version of novosort novosortV1.03.01.Linux3.0.tar.gz and put in context dir for the build

ADD novocraftV3.02.07.Linux3.0.tar.gz /usr/local/pipeline/
ADD novosortV1.03.01.Linux3.0.tar.gz /usr/local/pipeline/

RUN sed  -i '$aPATH=${PATH}:/usr/local/pipeline/novocraft' /home/pipeman/.bashrc \
      && chown -R pipeman:ngsgroup /usr/local/pipeline/novocraft \
      && echo "alias ngsNovoalign='/usr/local/pipeline/novocraft/novoalign'" >>  /home/pipeman/.bashrc

# + get novoalign.lic 
COPY novoalign.lic /usr/local/pipeline/novocraft/

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
