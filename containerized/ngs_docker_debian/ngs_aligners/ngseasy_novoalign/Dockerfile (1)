# Base image
FROM compbio/ngseasy-base:1.0

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

ADD novocraftV3.02.11.Linux3.0.tar.gz /usr/local/pipeline/
ADD novosortV1.03.01.Linux3.0.tar.gz /tmp/

RUN cp -v /tmp/novocraft/novosort /usr/local/pipeline/novocraft/ \
      && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/novocraft' /home/pipeman/.bashrc \
      && chmod -R 777 /usr/local/pipeline/novocraft \
      && cp -rfv /usr/local/pipeline/novocraft/* /usr/local/bin

# + get novoalign.lic 
COPY novoalign.lic /usr/local/pipeline/novocraft/

RUN cp -v /usr/local/pipeline/novocraft/novoalign.lic /usr/local/bin 

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
