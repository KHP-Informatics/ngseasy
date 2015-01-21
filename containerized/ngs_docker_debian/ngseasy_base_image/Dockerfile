# NGSeasy Base Image

# base image
FROM debian:jessie

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Remain current
RUN apt-get update && apt-get dist-upgrade -y

# Basic dependencies
RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  apt-utils \
  automake \
  ant \
  bash \
  binutils \
  perl \
  bioperl \
  build-essential \
  bzip2 \
  c++11 \
  cdbs \
  cmake \
  cron \
  curl \
  dkms \
  dpkg-dev \
  g++ \
  gpp \
  gcc \
  gfortran \
  git \
  git-core \
  libblas-dev \
  libatlas-dev \
  libbz2-dev \
  liblzma-dev \
  libpcre3-dev \
  libreadline-dev \
  make \
  mercurial \
  php5-curl \
  python python-dev python-yaml ncurses-dev zlib1g-dev python-numpy python-pip \
  sudo \
  subversion \
  tabix \
  tree \
  unzip \
  vim \
  wget \
  python-software-properties \
  libc-bin \
  llvm \
  libconfig-dev \
  ncurses-dev \
  zlib1g-dev \
  yum \
  libX11-dev libXpm-dev libXft-dev libXext-dev \
  asciidoc

#---------------------------------JAVA-------------------------------------------------------------------------------------#  
# upgrade java
RUN apt-get install -y openjdk-7-jdk openjdk-7-doc openjdk-7-jre-lib

#set java
ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java
RUN sed -i 'aPATH=$PATH:/usr/lib/jvm/java-7-openjdk-amd64/jre/bin' /root/.bashrc

#-------------------------------Add user----------------------------------------------------------------------------------#
# Create a pipeline user:pipeman and group:ngsgroup

RUN useradd -m -s /bin/bash pipeman && \
  cd /home/pipeman && \
  echo "#bash config file for user pipeman" >> /home/pipeman/.bashrc

RUN groupadd ngsgroup && \
  usermod -aG ngsgroup pipeman && \
  usermod -aG sudo pipeman

#-----------------------------NGS TOOLS DIRECTORY------------------------------------------------------------------------#  
#make pipeline install dirs
RUN mkdir /usr/local/pipeline && \
  chown pipeman:ngsgroup /usr/local/pipeline && \
  chmod 775 /usr/local/pipeline

#--------------STANDARD NGS TOOLS----------------------------------------------------------------------------------------------#
# Tools used for processing SAM/BAM/BED/VCF files
# samtools,htslib,bcftools,parallel,bamUtil,sambamba,samblaster,vcftools,vcflib,seqtk,ogap,bamleftalign,bedtools2,libStatGen

# ngs tools     

# samtools, htslib, bcftools
RUN cd /usr/local/pipeline && \
    git clone --branch=develop git://github.com/samtools/htslib.git && \
    git clone --branch=develop git://github.com/samtools/bcftools.git && \
    git clone --branch=develop git://github.com/samtools/samtools.git && \
    cd /usr/local/pipeline/htslib && \
    make && \
    make install && \
    cd /usr/local/pipeline/bcftools && \
    make && \
    make install && \
    cd /usr/local/pipeline/samtools && \
    make && \
    make install    

# parallel    
RUN cd /usr/local/pipeline && \
    wget http://ftpmirror.gnu.org/parallel/parallel-20140222.tar.bz2 && \
    bzip2 -dc parallel-20140222.tar.bz2 | tar xvf - && \
    cd parallel-20140222 && \
    ./configure && \
    make && \
    make install

# libStatGen    
RUN cd /usr/local/pipeline && \
    git clone https://github.com/statgen/libStatGen.git && \
    cd libStatGen && \
    make all

# bamUtil    
RUN cd /usr/local/pipeline && \
    git clone https://github.com/statgen/bamUtil.git && \
    cd bamUtil && \
    make cloneLib && \
    make all && \
    make install

# samblaster    
RUN cd /usr/local/pipeline && \ 
    git clone git://github.com/GregoryFaust/samblaster.git && \ 
    cd samblaster && \ 
    make && \ 
    cp samblaster /usr/local/bin/
    
# sambamba
RUN cd /usr/local/pipeline && \ 
    curl -OL https://github.com/lomereiter/sambamba/releases/download/v0.5.0-alpha/sambamba_v5.0.0-alpha_linux.tar.bz2 && \ 
    tar -xvf sambamba_v5.0.0-alpha_linux.tar.bz2 && \ 
    cp sambamba_v5.0.0-alpha sambamba && \
    chmod +rwx sambamba && \
    cp sambamba /usr/local/bin/

# vcftools 
RUN cd /usr/local/pipeline/ && \
    wget -O /tmp/vcftools_0.1.12b.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12b.tar.gz && \
    tar xzvf /tmp/vcftools_0.1.12b.tar.gz -C /usr/local/pipeline/  && \
    export PERL5LIB=/usr/local/pipeline/vcftools_0.1.12b/perl/  && \
    cd /usr/local/pipeline/vcftools_0.1.12b/ && \
    make && \
    sed  -i '$aPATH=${PATH}:/usr/local/pipeline/vcftools_0.1.12b/bin' /root/.bashrc

# seqtk and trimadap
RUN cd /usr/local/pipeline/ && \
    git clone https://github.com/lh3/seqtk.git && \
    chown -R pipeman:ngsgroup /usr/local/pipeline/seqtk && \
    cd seqtk/ && \
    chmod -R 777 ./* && \
    make && \
    cp -v seqtk /usr/local/bin/ && \
    cp -v trimadap /usr/local/bin/ && \
    sed  -i '$aPATH=${PATH}:/usr/local/pipeline/seqtk' /root/.bashrc

# ogap    
RUN cd /usr/local/pipeline/ && \
    git clone --recursive https://github.com/ekg/ogap.git && \
    cd ogap && \
    make all && \
    chmod -R 777 ./* && \
    cp -v ogap /usr/local/bin/

# bamleftalign (from Freebayes)
RUN cd /usr/local/pipeline/ && \
    git clone --recursive git://github.com/ekg/freebayes.git && \
    cd freebayes && \
    make all && \
    chmod -R 777 ./* && \
    cp bin/bamleftalign /usr/local/bin/ && \
    rm -frv /usr/local/pipeline/freebayes

# vcflib    
RUN cd /usr/local/pipeline/ && \
    rm -rfv /usr/local/pipeline/vcflib && \
    git clone --recursive git://github.com/ekg/vcflib.git && \
    cd vcflib && \
    make && \
    cp ./bin/* /usr/local/bin/
 
# bedtools2
RUN cd /usr/local/pipeline && \
    git clone https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && \
    make clean && \
    make all && \
    make install
  
#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline 
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rvf /tmp/*

#open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get remove -y asciidoc && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/



