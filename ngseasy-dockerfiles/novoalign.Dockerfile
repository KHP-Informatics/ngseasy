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

ADD novocraftV3.02.11.Linux3.0.tar.gz /usr/local/ngs/bin/
ADD novosortV1.03.01.Linux3.0.tar.gz /tmp/

RUN cp -v /tmp/novocraft/novosort /usr/local/ngs/bin/novocraft/ \
      && sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/novocraft' /home/ngseasy/.bashrc \
      && chmod -R 777 /usr/local/ngs/bin/novocraft \
      && cp -rfv /usr/local/ngs/bin/novocraft/* /usr/local/bin

# + get novoalign.lic 
COPY novoalign.lic /usr/local/ngs/bin/novocraft/

RUN cp -v /usr/local/ngs/bin/novocraft/novoalign.lic /usr/local/bin 

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/ngs/bin/***
RUN chown -R ngseasy:ngseasy /usr/local/ngs/bin

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 8080

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
USER ngseasy
WORKDIR /home/ngseasy
