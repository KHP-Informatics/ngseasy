# NGSeasy Base Image

FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Remain current
RUN apt-get update &&  apt-get upgrade -y && apt-get dist-upgrade -y

# glia
RUN cd /usr/local/pipeline && \
    git clone --recursive https://github.com/ekg/glia.git && \
    chmod -R 777 ./glia && \
    cd ./glia && \
    make && \
    cp -v ./glia /usr/local/bin/ && \
    cd /usr/local/pipeline && \
    git clone --branch v0.9.20 --recursive git://github.com/ekg/freebayes.git  && \
    chmod -R 777 freebayes && \
    cd freebayes && \
    make && \
    cp -v ./bin/* /usr/local/bin/
  
#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline 
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rvf /tmp/*

#open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/ && \
    rm -rf /usr/local/pipeline/*



