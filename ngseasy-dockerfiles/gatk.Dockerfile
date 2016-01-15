# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

#-------------------------NGS-TOOL---------------------------------------
ADD GenomeAnalysisTK-3.3-0.tar.bz2 /usr/local/ngs/bin/GenomeAnalysisTK-3.3-0/

RUN cd /usr/local/ngs/bin/GenomeAnalysisTK-3.3-0 && \
	sed -i '$aCLASSPATH=:${CLASSPATH}:/usr/local/ngs/bin/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar' /home/ngseasy/.bashrc && \
	sed -i '$aCLASSPATH=:${CLASSPATH}:/usr/local/ngs/bin/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar' /root/.bashrc && \
	sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/GenomeAnalysisTK-3.3-0' /home/ngseasy/.bashrc && \
	sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/GenomeAnalysisTK-3.3-0' /root/.bashrc && \
	chmod 755 /usr/local/ngs/bin/GenomeAnalysisTK-3.3-0/* && \	
	cp -rfv /usr/local/ngs/bin/GenomeAnalysisTK-3.3-0/* /usr/local/bin/

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
