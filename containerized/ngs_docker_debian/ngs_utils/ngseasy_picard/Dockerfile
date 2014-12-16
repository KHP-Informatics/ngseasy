# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# + Picard
RUN wget -O /tmp/picard-tools-1.119.zip http://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip \
    && mkdir /usr/local/pipeline/picardtools \
    && unzip /tmp/picard-tools-1.119.zip -d /usr/local/pipeline/picardtools/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/picardtools \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/picardtools/picard-tools-1.119/snappy-java-1.0.3-rc3.jar' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/picardtools/picard-tools-1.119' /home/pipeman/.bashrc \
    && echo "alias ngsPicard='/usr/local/pipeline/picardtools/picard-tools-1.119'" >>  /home/pipeman/.bashrc


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

