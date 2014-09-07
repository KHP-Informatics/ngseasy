# Base image will be Trusty
FROM afolarin/base-config


# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com

# Set correct environment variables.
ENV HOME /root

# install basic dependencies
RUN apt-get update
RUN apt-get install -y gcc 
RUN apt-get install -y g++
RUN apt-get install -y python-dev
RUN apt-get install -y make 
RUN apt-get install -y build-essential 
RUN apt-get install -y cmake
RUN apt-get install -y tabix
RUN apt-get install -y dpkg-dev 
RUN apt-get install -y binutils
RUN apt-get install -y gfortran
RUN apt-get install -y git-core 
RUN apt-get install -y wget 
RUN apt-get install -y unzip bzip2
RUN apt-get install -y zlib1g-dev 
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y sysvbanner
RUN apt-get install -y openjdk-7-jre
RUN apt-get install -y openjdk-7-jdk
RUN apt-get install -y curl 
RUN apt-get install -y libcurl3
RUN apt-get install -y php5-curl 
RUN apt-get install -y libxml2 
RUN apt-get install -y libxslt-dev 
RUN apt-get install -y libxml2-dev 
RUN apt-get install -y fastx-toolkit
RUN apt-get install -y htmldoc 
RUN apt-get install -y html2text
RUN apt-get install -y perl
RUN apt-get install -y bioperl
RUN apt-get install -y r-base-dev



#-------------------------------------v-callers R----------------------------------#

## bioconductor and R libs ##
## needed fo R CNV tools  
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite()' \
	&& /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("GenomicRanges",dependencies=TRUE)' \
	&& /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("IRanges",dependencies=TRUE)' \
	&& /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("Biostrings",dependencies=TRUE)' \
	&& /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("genomeIntervals",dependencies=TRUE)' \
	&& /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("rtracklayer",dependencies=TRUE)' \
	&& /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome",dependencies=TRUE)' \
	&& /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("Rsamtools",dependencies=TRUE)' \
	&& /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("ShortRead",dependencies=TRUE)'  \
## R ExomeDepth - Exome Seq \
	&& /usr/bin/Rscript --no-save --no-restore -e 'install.packages("ExomeDepth",dependencies=TRUE ,repos="http://cran.r-project.org")'

## R m-HMM - WGS 
# https://www.stt.msu.edu/users/hengwang/mHMM.html
RUN cd /usr/local/pipeline/ \
&& wget https://www.stt.msu.edu/users/hengwang/mHMM%20Data/mHMM_1.0.tar.gz \
&& sudo R CMD INSTALL mHMM_1.0.tar.gz \
&& rm mHMM_1.0.tar.gz;

#---------------------------------sv-callers-------------------------------------#
# + DELLY - split read PE mapping
# install bamtools
RUN cd /usr/local/pipeline \
&& git clone https://github.com/pezmaster31/bamtools.git \
&& cd bamtools/\
&& mkdir build \
&& cd build/ \
&& cmake ..\
&& make;

# set environment
ENV BOOST_ROOT /usr
ENV BAMTOOLS_ROOT /usr/local/pipeline/bamtools
ENV SEQTK_ROOT /usr/local/pipeline/seqtk
ENV LD_LIBRARY_PATH /usr/local/pipeline/bamtools/lib

RUN echo "export BOOST_ROOT=/usr" >>  /home/pipeman/.bashrc \
&& echo "export BAMTOOLS_ROOT=/usr/local/pipeline/bamtools" >>  /home/pipeman/.bashrc \
&& echo "export SEQTK_ROOT=/usr/local/pipeline/seqtk" >>  /home/pipeman/.bashrc \
&& echo "export LD_LIBRARY_PATH=/usr/local/pipeline/bamtools/lib" >>  /home/pipeman/.bashrc \
&& source /home/pipeman/.bashrc

# install delly
RUN cd /usr/local/pipeline \
&& git clone https://github.com/tobiasrausch/delly.git \
&& cd delly/ ; make -B src/delly \
&& sed -i '$aPATH=${PATH}:/usr/local/pipeline/delly/src' /home/pipeman/.bashrc \
&& echo "alias ngsDelly='/usr/local/pipeline/delly/src'" >> /home/pipeman/.bashrc \
&& source /home/pipeman/.bashrc

# + lumpy  - split read PE mapping
RUN  apt-get install python-numpy \
&& cd /usr/local/pipeline \
&& git clone git://github.com/arq5x/lumpy-sv.git \
&& cd cd lumpy-sv \
&& make \
&& sed -i '$aPATH=${PATH}:/usr/local/pipeline/lumpy-sv/bin' /home/pipeman/.bashrc \
&& echo "alias ngsLumpy='/usr/local/pipeline/lumpy-sv/bin'" >> /home/pipeman/.bashrc \
&& source /home/pipeman/.bashrc

## CNVnator - readdepth mapping

# install root
# ./configure --help
RUN cd /usr/local/pipeline/ \
&& git clone http://root.cern.ch/git/root.git \
&& cd root \
&& git tag -l \
&& git checkout -b v5-34-00-patch \
&& ./configure linuxx8664gcc  \
&& make \
&& source bin/thisroot.sh \
&& source /home/pipeman/.bashrc \
&& echo "export ROOTSYS=/usr/local/pipeline/root" >>  /home/pipeman/.bashrc \
&& sed -i '$aLD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib' /home/pipeman/.bashrc \
&& source /home/pipeman/.bashrc \


# install CNVnator
RUN cd /usr/local/pipeline/ \
&& sudo apt-get install root-system-common \
&& sudo apt-get install libroot-dev \
&& wget http://sv.gersteinlab.org/cnvnator/CNVnator_v0.3.zip \
&& unzip CNVnator_v0.3.zip \
&& export ROOTSYS=/usr/local/pipeline/root \
&& export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib \
&& source /home/pipeman/.bashrc \
&& cd /usr/local/pipeline/CNVnator_v0.3/src/samtools \
&& make \
&& cd /usr/local/pipeline/CNVnator_v0.3/src \
&& cp -v Makefile Makefile.orginal \
&& sed -i 's/g++ /g++ -std=c++0x -pthread /' Makefile \
&& make \
&& sed -i '$aPATH=${PATH}:/usr/local/pipeline/CNVnator_v0.3/src' /home/pipeman/.bashrc \

