# Base image
FROM compbio/ngseasy-base:1.0-r001

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# + bowtie
RUN wget -O /tmp/bowtie2-2.2.5-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip \
    && unzip /tmp/bowtie2-2.2.5-linux-x86_64.zip -d /usr/local/pipeline/ \
    && chmod -R 777 /usr/local/pipeline/ \  
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.5:/usr/local/pipeline/bowtie2-2.2.5/scripts' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.5:/usr/local/pipeline/bowtie2-2.2.5/scripts' ~/.bashrc \
    && cp -v /usr/local/pipeline/bowtie2-2.2.5/bowtie* /usr/local/bin/

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

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


