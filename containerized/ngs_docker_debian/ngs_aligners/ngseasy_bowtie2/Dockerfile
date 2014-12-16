# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y


# + bowtie
RUN wget -O /tmp/bowtie2-2.2.4-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip \
    && unzip /tmp/bowtie2-2.2.4-linux-x86_64.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bowtie2-2.2.4 \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.4:/usr/local/pipeline/bowtie2-2.2.4/scripts' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.4:/usr/local/pipeline/bowtie2-2.2.4/scripts' ~/.bashrc \
    && echo "alias ngsBowtie2='/usr/local/pipeline/bowtie2-2.2.4'" >>  /home/pipeman/.bashrc

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


