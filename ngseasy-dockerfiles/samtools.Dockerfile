# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y && \
	apt-get install -y debhelper dpkg-dev g++-multilib  libgpm-dev pkg-config ncurses-base libncursesw5 libncurses5-dev
 
# + samtools, htslib and bcftools
RUN cd /usr/local/pipeline \
	&& git clone --branch=develop git://github.com/samtools/htslib.git \
	&& git clone --branch=develop git://github.com/samtools/bcftools.git \
	&& git clone --branch=develop git://github.com/samtools/samtools.git \
	&& cd /usr/local/pipeline/bcftools \
	&& make \
	&& cd /usr/local/pipeline/samtools \
	&& make \
	&& cd /usr/local/pipeline/htslib \
	&& make \
	&& chown -R ngseasy:ngseasy /usr/local/pipeline/samtools \
	&& chown -R ngseasy:ngseasy /usr/local/pipeline/bcftools \
	&& chown -R ngseasy:ngseasy /usr/local/pipeline/htslib \
	&& sed  -i '$aPATH=${PATH}:/usr/local/pipeline/samtools' /home/ngseasy/.bashrc \
	&& sed  -i '$aPATH=${PATH}:/usr/local/pipeline/samtools' ~/.bashrc \
	&& echo "alias ngsSAMtools='/usr/local/pipeline/samtools'" >>  /home/ngseasy/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bcftools' /home/ngseasy/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bcftools' ~/.bashrc \        
        && echo "alias ngsBCFtools='/usr/local/pipeline/bcftools'" >>  /home/ngseasy/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/htslib' /home/ngseasy/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/htslib' ~/.bashrc \        
        && echo "alias ngsHTSlib='/usr/local/pipeline/htslib'" >>  /home/ngseasy/.bashrc 


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
