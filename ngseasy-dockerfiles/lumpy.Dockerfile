# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

#-------------------------------lumpy--------------------------   
# + lumpy split read PE mapping

RUN apt-get install -y python-numpy

RUN cd /usr/local/ngs/bin \
    && git clone git://github.com/arq5x/lumpy-sv.git \
    && cd /usr/local/ngs/bin/lumpy-sv \
    && make \
    && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/lumpy-sv/bin' /home/ngseasy/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/lumpy-sv/bin' ~/.bashrc

ADD fix_ambiguous /usr/local/bin/
#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 755 /usr/local/ngs/bin
RUN chown -R ngseasy:ngseasy /usr/local/ngs/bin

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
