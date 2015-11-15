# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# + Picard
RUN wget -O /tmp/picard-tools-1.129.zip https://github.com/broadinstitute/picard/releases/download/1.129/picard-tools-1.129.zip \
    && mkdir /usr/local/ngs/bin/picardtools \
    && unzip /tmp/picard-tools-1.129.zip -d /usr/local/ngs/bin/picardtools/ \
    && chown -R ngseasy:ngseasy /usr/local/ngs/bin/picardtools \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/picardtools/picard-tools-1.129/snappy-java-1.0.3-rc3.jar' /home/ngseasy/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/picardtools/picard-tools-1.129' /home/ngseasy/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/picardtools/picard-tools-1.129' ~/.bashrc

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/ngs/bin
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
