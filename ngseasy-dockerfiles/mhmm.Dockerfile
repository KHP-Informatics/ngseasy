# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

RUN apt-get install -y \
	libatlas3-base \
	libblas3 \
	liblzma5 \
	libpango1.0-0 \
	libpaper-utils \
	libpcre3 \
	libpng12-0 \
	libquadmath0  \
	libreadline6 \
	libsm6 \
	libx11-6 \
	libxext6 \
	libxss1 \
	libxt6 \
	tcl8.5 \
	tk8.5 \
	ucf \
	unzip \
	xdg-utils \
	zip \
	zlib1g \
	ed \
	less \
	littler \
	locales \
	r-base-core \
	r-base-dev \
	r-recommended \
	r-cran-vgam \
	r-cran-rsqlite

#  upgrade
RUN apt-get update -y && apt-get upgrade -y

#---------------------------------m-HMM -------------------------------------#
  
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite()' && \
	cd /usr/local/pipeline/ && \
	curl -k -L -O https://www.stt.msu.edu/users/hengwang/mHMM%20Data/mHMM_1.0.tar.gz && \
	R CMD INSTALL mHMM_1.0.tar.gz && \
	rm mHMM_1.0.tar.gz;

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R ngseasy:ngseasy /usr/local/pipeline

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

