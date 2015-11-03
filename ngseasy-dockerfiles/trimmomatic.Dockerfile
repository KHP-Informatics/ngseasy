#------------------------------------------------#
# FASTQC TOOL: Dockerfile
#------------------------------------------------#

# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Remain current
RUN apt-get update && apt-get upgrade -y

#-------------------------NGS-TOOL---------------------------------------

# + Trimmomatic
    RUN wget -O /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip \
        && unzip /tmp/Trimmomatic-0.32.zip -d /usr/local/pipeline/ \
        && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar' /home/ngseasy/.bashrc \
        && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Trimmomatic-0.32' /home/ngseasy/.bashrc \
        && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Trimmomatic-0.32' /root/.bashrc \
        && cp -v /usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar /usr/local/bin
        
#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R ngseasy:ngseasy /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/

