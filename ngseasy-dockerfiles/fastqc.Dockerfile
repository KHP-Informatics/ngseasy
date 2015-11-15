# Base image
FROM compbio/ngseasy-ubuntu-14.04.3:1.0-r002
# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com
# + FastQC
    RUN wget -O /tmp/fastqc_v0.11.2.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip \
        && unzip /tmp/fastqc_v0.11.2.zip -d /usr/local/ngs/bin/ \
        && chmod -R 766 /usr/local/ngs/bin/ \
        && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/FastQC/jbzip2-0.9.jar:/usr/local/ngs/bin/FastQC/sam-1.103.jar' /home/pipeman/.bashrc \
        && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/FastQC' /home/pipeman/.bashrc \
        && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/FastQC' /root/.bashrc \
        && ln -s /usr/local/ngs/bin/FastQC/fastqc /usr/local/bin/fastqc

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/ngs/bin
RUN chown -R pipeman:ngsgroup /usr/local/ngs/bin

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
