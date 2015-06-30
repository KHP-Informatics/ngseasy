#------------------------------------------------#
# NGSeasy Dev Tool Box
#------------------------------------------------#

# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# FastQC & Trimmomatic ---------------------------------------------------------------------

RUN wget -O /tmp/fastqc_v0.11.2.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip \
  && unzip /tmp/fastqc_v0.11.2.zip -d /usr/local/pipeline/ \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/FastQC \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/FastQC/jbzip2-0.9.jar:/usr/local/pipeline/FastQC/sam-1.103.jar' /home/pipeman/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/FastQC' /home/pipeman/.bashrc \
  && echo "alias ngsFastqc='/usr/local/pipeline/FastQC'" >>  /home/pipeman/.bashrc \
  && cp -vr /usr/local/pipeline/FastQC/fastqc /usr/local/bin \
  && wget -O /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip \
  && unzip /tmp/Trimmomatic-0.32.zip -d /usr/local/pipeline/ \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/Trimmomatic-0.32 \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar' /home/pipeman/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Trimmomatic-0.32' /home/pipeman/.bashrc \
  && echo "alias ngsTrimfq='/usr/local/pipeline/Trimmomatic-0.32'" >>  /home/pipeman/.bashrc \
  && cp -vr /usr/local/pipeline/Trimmomatic-0.32/* /usr/local/bin \
   
# picard ----------------------------------------------------------------------
RUN wget -O /tmp/picard-tools-1.119.zip http://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip \
    && mkdir /usr/local/pipeline/picardtools \
    && unzip /tmp/picard-tools-1.119.zip -d /usr/local/pipeline/picardtools/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/picardtools \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/picardtools/picard-tools-1.119/snappy-java-1.0.3-rc3.jar' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/picardtools/picard-tools-1.119' /home/pipeman/.bashrc \
    && echo "alias ngsPicard='/usr/local/pipeline/picardtools/picard-tools-1.119'" >>  /home/pipeman/.bashrc \
    && cp -vr /usr/local/pipeline/picardtools/picard-tools-1.119/* /usr/local/bin
  
# Aligners ---------------------------------------------------------------------
# bwa, stampy, bowtie2
# manually build novoalign and commit local image

RUN wget -O /tmp/bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2 \
    && tar xjvf /tmp/bwa-0.7.10.tar.bz2 -C /usr/local/pipeline/ \
    && cd /usr/local/pipeline/bwa-0.7.10 && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bwa-0.7.10 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bwa-0.7.10' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bwa-0.7.10' ~/.bashrc \
    && echo "alias ngsBWA='/usr/local/pipeline/bwa-0.7.10'" >>  /home/pipeman/.bashrc \
    && cp -vr /usr/local/pipeline/bwa-0.7.10/* /usr/local/bin \
    && wget -O /tmp/stampy-latest.tgz  http://www.well.ox.ac.uk/~gerton/software/Stampy/stampy-latest.tgz \
    && tar xvf /tmp/stampy-latest.tgz -C /usr/local/pipeline/ \
    && sed -i 's/-Wl//' /usr/local/pipeline/stampy-1.0.23/makefile \
    && chmod -R 755 /usr/local/pipeline/stampy-1.0.23 \
    && cd /usr/local/pipeline/stampy-1.0.23 && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/stampy-1.0.23 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/stampy-1.0.23' /home/pipeman/.bashrc \
    && echo "alias ngsStampy='/usr/local/pipeline/stampy-1.0.23'" >>  /home/pipeman/.bashrc \
    && cp -vr /usr/local/pipeline/stampy-1.0.23/* /usr/local/bin \
    && wget -O /tmp/bowtie2-2.2.4-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip \
    && unzip /tmp/bowtie2-2.2.4-linux-x86_64.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bowtie2-2.2.4 \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.4:/usr/local/pipeline/bowtie2-2.2.4/scripts' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.4:/usr/local/pipeline/bowtie2-2.2.4/scripts' ~/.bashrc \
    && echo "alias ngsBowtie2='/usr/local/pipeline/bowtie2-2.2.4'" >>  /home/pipeman/.bashrc \
    && cp -vr /usr/local/pipeline/bowtie2-2.2.4/* /usr/local/bin
    
# Variant callers -------------------------------------------------------------------------------------

RUN cd /usr/local/pipeline \
  && git clone --recursive git://github.com/ekg/freebayes.git \
  && cd /usr/local/pipeline/freebayes \
  && make \
  && chmod -R 755 /usr/local/pipeline/freebayes \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/freebayes \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/freebayes/bin' /home/pipeman/.bashrc \
  && echo "alias ngsFreebayes='/usr/local/pipeline/freebayes/bin'" >>  /home/pipeman/.bashrc \
  && cp -vr /usr/local/pipeline/freebayes/bin/* /usr/local/bin \
  && wget -O /tmp/Platypus-latest.tgz http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz \
  && tar xvf /tmp/Platypus-latest.tgz -C /usr/local/pipeline/ \
  && cd /usr/local/pipeline/Platypus_0.7.9.1 \
  && sh ./buildPlatypus.sh \
  && chmod -R 755 /usr/local/pipeline/Platypus_0.7.9.1 \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/Platypus_0.7.9.1 \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Platypus_0.7.9.1' /home/pipeman/.bashrc \
  && echo "alias ngsPlatypus='/usr/local/pipeline/Platypus_0.7.9.1'" >>  /home/pipeman/.bashrc \
  && cp -vr /usr/local/pipeline/Platypus_0.7.9.1/* /usr/local/bin \
  &&  mkdir /usr/local/pipeline/bcbio \
  && cd /usr/local/pipeline/bcbio \
  && wget https://github.com/chapmanb/bcbio.variation/releases/download/v0.1.9/bcbio.variation-0.1.9-standalone.jar \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/bcbio/bcbio.variation-0.1.9-standalone.jar' /home/pipeman/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bcbio/' /home/pipeman/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bcbio/' ~/.bashrc \
  && echo "alias ngsBCBIOVAR'/usr/local/pipeline/bcbio/'" >>  /home/pipeman/.bashrc \
  && cp -vr /usr/local/pipeline/bcbio/* /usr/local/bin 


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


















    
    
    
    
    
    
    