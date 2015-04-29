# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

## CHANJO #####################################################

RUN cd /usr/local/pipeline/ && \
	wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
	bash Miniconda-latest-Linux-x86_64.sh -b && \
	sed -i '$aPATH=${PATH}:${HOME}/miniconda/bin' /home/pipeman/.bashrc && \
	sed -i '$aPATH=${PATH}:${HOME}/miniconda/bin' ~/.bashrc && \
	export PATH=${PATH}:${HOME}/miniconda/bin

RUN	/root/miniconda/bin/conda create -n chanjo2 python=2 cython numpy sqlalchemy pip && \
	source activate chanjo2 && \
	/root/miniconda/bin/conda install -c https://conda.binstar.org/robinandeer pysam && \
	pip install chanjo && \
	pip install pymysql && \
	pip install chanjo-report

#-------------------------------PERMISSIONS-------------------------
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