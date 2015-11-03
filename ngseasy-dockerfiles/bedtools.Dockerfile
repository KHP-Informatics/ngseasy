# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# + BEDtools
RUN cd /usr/local/pipeline \
    && git clone https://github.com/arq5x/bedtools2.git \
    && cd bedtools2 && make clean && make all \
    && chown -R ngseasy:ngseasy /usr/local/pipeline/bedtools2 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bedtools2/bin' /home/ngseasy/.bashrc \
    && echo "alias ngsBedtools='/usr/local/pipeline/bedtools2/bin'" >> /home/ngseasy/.bashrc


#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R ngseasy:ngseasy /usr/local/pipeline

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

