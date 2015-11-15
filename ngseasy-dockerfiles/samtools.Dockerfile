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
RUN cd /usr/local/ngs/bin \
	&& git clone --branch=develop git://github.com/samtools/htslib.git \
	&& git clone --branch=develop git://github.com/samtools/bcftools.git \
	&& git clone --branch=develop git://github.com/samtools/samtools.git \
	&& cd /usr/local/ngs/bin/bcftools \
	&& make \
	&& cd /usr/local/ngs/bin/samtools \
	&& make \
	&& cd /usr/local/ngs/bin/htslib \
	&& make \
	&& chown -R ngseasy:ngseasy /usr/local/ngs/bin/samtools \
	&& chown -R ngseasy:ngseasy /usr/local/ngs/bin/bcftools \
	&& chown -R ngseasy:ngseasy /usr/local/ngs/bin/htslib \
	&& sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/samtools' /home/ngseasy/.bashrc \
	&& sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/samtools' ~/.bashrc \
	&& echo "alias ngsSAMtools='/usr/local/ngs/bin/samtools'" >>  /home/ngseasy/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/bcftools' /home/ngseasy/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/bcftools' ~/.bashrc \        
        && echo "alias ngsBCFtools='/usr/local/ngs/bin/bcftools'" >>  /home/ngseasy/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/htslib' /home/ngseasy/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/htslib' ~/.bashrc \        
        && echo "alias ngsHTSlib='/usr/local/ngs/bin/htslib'" >>  /home/ngseasy/.bashrc 


#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/ngs/bin
RUN chown -R ngseasy:ngseasy /usr/local/ngs/bin

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
USER ngseasy
WORKDIR /home/ngseasy
