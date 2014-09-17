# Base image
FROM snewhouse/ngseasy-base:v1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Regenerate SSH host keys. baseimage-docker does not contain any, so you
# have to do that yourself. You may also comment out this instruction; the
# init system will auto-generate one during boot.
RUN /etc/my_init.d/00_regen_ssh_host_keys.sh

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update

# + Trimmomatic
    RUN wget -O /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip \
        && unzip /tmp/Trimmomatic-0.32.zip -d /usr/local/pipeline/ \
        && chown -R pipeman:ngsgroup /usr/local/pipeline/Trimmomatic-0.32 \
        && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar' /home/pipeman/.bashrc \
        && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Trimmomatic-0.32' /home/pipeman/.bashrc \
        && echo "alias ngsTrimfq='/usr/local/pipeline/Trimmomatic-0.32'" >>  /home/pipeman/.bashrc
        
# + FastQC
    RUN wget -O /tmp/fastqc_v0.11.2.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip \
        && unzip /tmp/fastqc_v0.11.2.zip -d /usr/local/pipeline/ \
        && chown -R pipeman:ngsgroup /usr/local/pipeline/FastQC \
        && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/FastQC/jbzip2-0.9.jar:/usr/local/pipeline/FastQC/sam-1.103.jar' /home/pipeman/.bashrc \
        && sed -i '$aPATH=${PATH}:/usr/local/pipeline/FastQC' /home/pipeman/.bashrc \
        && echo "alias ngsFastqc='/usr/local/pipeline/FastQC'" >>  /home/pipeman/.bashrc    

# + seqtk
    RUN cd /usr/local/pipeline/ \
    && git clone https://github.com/lh3/seqtk.git \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/seqtk \
    && cd seqtk/ \
    && make \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/seqtk' /home/pipeman/.bashrc \
    && echo "alias ngsSeqtk='/usr/local/pipeline/seqtk'" >>  /home/pipeman/.bashrc \
    && cd

# + novoalign (registration required,  get compressed file and put in context dir for the build)
# + get novoalign.lic from novoalign ~ $1000 and put in context dir for the build
# + get updated version of novosort novosortV1.03.01.Linux3.0.tar.gz and put in context dir for the build
    RUN wget -O /tmp/novocraftV3.02.07.Linux3.0.tar.gz https://www.dropbox.com/s/kbgds9e6ucog6oj/novocraftV3.02.07.Linux3.0.tar.gz \
      && tar xvf /tmp/novocraftV3.02.07.Linux3.0.tar.gz  -C /usr/local/pipeline/ \
      && cp -v /usr/local/pipeline/novocraft/novosort /usr/local/pipeline/novocraft/novosort_original \
      && rm /usr/local/pipeline/novocraft/novosort \
      && wget -O /tmp/novosortV1.03.01.Linux3.0.tar.gz https://www.dropbox.com/s/1j4w9rpyoug1pbo/novosortV1.03.01.Linux3.0.tar.gz \
      && tar xvf /tmp/novosortV1.03.01.Linux3.0.tar.gz -C /usr/local/pipeline/novocraft \
      && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/novocraft' /home/pipeman/.bashrc \
      && chown -R pipeman:ngsgroup /usr/local/pipeline/novocraft \
      && echo "alias ngsNovoalign='/usr/local/pipeline/novocraft/novoalign'" >>  /home/pipeman/.bashrc

# + get novoalign.lic 
ADD novoalign.lic /usr/local/pipeline/novocraft/
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/novocraft
  
# + bowtie
RUN wget -O /tmp/bowtie2-2.2.3-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip \
    && unzip /tmp/bowtie2-2.2.3-linux-x86_64.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bowtie2-2.2.3 \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.3:/usr/local/pipeline/bowtie2-2.2.3/scripts' /home/pipeman/.bashrc \
    && echo "alias ngsBowtie2='/usr/local/pipeline/bowtie2-2.2.3'" >>  /home/pipeman/.bashrc

# + bwa 
RUN wget -O /tmp/bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2 \
    && tar xjvf /tmp/bwa-0.7.10.tar.bz2 -C /usr/local/pipeline/ \
    && cd /usr/local/pipeline/bwa-0.7.10 && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bwa-0.7.10 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bwa-0.7.10' /home/pipeman/.bashrc \
    && echo "alias ngsBWA='/usr/local/pipeline/bwa-0.7.10'" >>  /home/pipeman/.bashrc
    
# + stampy (registration required, get compressed file and put in context dir for the build)
RUN wget -O /tmp/stampy-latest.tgz  http://www.well.ox.ac.uk/~gerton/software/Stampy/stampy-latest.tgz \
    && tar xvf /tmp/stampy-latest.tgz -C /usr/local/pipeline/ \
    && sed -i 's/-Wl//' /usr/local/pipeline/stampy-1.0.23/makefile \
    && chmod -R 755 /usr/local/pipeline/stampy-1.0.23 \
    && cd /usr/local/pipeline/stampy-1.0.23 && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/stampy-1.0.23 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/stampy-1.0.23' /home/pipeman/.bashrc \
    && echo "alias ngsStampy='/usr/local/pipeline/stampy-1.0.23'" >>  /home/pipeman/.bashrc

# + GATK (see licence, and registration) 
# use ADD to copy the downloaded GATK binary from <buildcontext_path>/GenomeAnalysisTK-3.2-2.tar.bz2 to the /tmp dir in the container
# NOTE: ADD automatically unpacks compressed files
RUN wget -O /tmp/GenomeAnalysisTK-3.2-2.tar.bz2 https://www.dropbox.com/s/wey1edv9tqdc7so/GenomeAnalysisTK-3.2-2.tar.bz2 \
  && mkdir /usr/local/pipeline/GenomeAnalysisTK-3.2-2 \ 
  && tar xjvf /tmp/GenomeAnalysisTK-3.2-2.tar.bz2 -C /usr/local/pipeline/GenomeAnalysisTK-3.2-2 \
  && chmod -R 755 /usr/local/pipeline/GenomeAnalysisTK-3.2-2 \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/GenomeAnalysisTK-3.2-2 \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar' /home/pipeman/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/GenomeAnalysisTK-3.2-2' /home/pipeman/.bashrc \
  && echo "alias ngsGATK='/usr/local/pipeline/GenomeAnalysisTK-3.2-2'" >>  /home/pipeman/.bashrc
 
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
	&& chown -R pipeman:ngsgroup /usr/local/pipeline/samtools \
	&& chown -R pipeman:ngsgroup /usr/local/pipeline/bcftools \
	&& chown -R pipeman:ngsgroup /usr/local/pipeline/htslib \
	&& sed  -i '$aPATH=${PATH}:/usr/local/pipeline/samtools' /home/pipeman/.bashrc \
	&& echo "alias ngsSAMtools='/usr/local/pipeline/samtools'" >>  /home/pipeman/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bcftools' /home/pipeman/.bashrc \
        && echo "alias ngsBCFtools='/usr/local/pipeline/bcftools'" >>  /home/pipeman/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/htslib' /home/pipeman/.bashrc \
        && echo "alias ngsHTSlib='/usr/local/pipeline/htslib'" >>  /home/pipeman/.bashrc 

# + Picard
RUN wget -O /tmp/picard-tools-1.115.zip http://sourceforge.net/projects/picard/files/picard-tools/1.115/picard-tools-1.115.zip \
    && mkdir /usr/local/pipeline/picardtools \
    && unzip /tmp/picard-tools-1.115.zip -d /usr/local/pipeline/picardtools/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/picardtools \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/picardtools/picard-tools-1.115/snappy-java-1.0.3-rc3.jar' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/picardtools/picard-tools-1.115' /home/pipeman/.bashrc \
    && echo "alias ngsPicard='/usr/local/pipeline/picardtools/picard-tools-1.115'" >>  /home/pipeman/.bashrc

# + BEDtools
RUN cd /usr/local/pipeline \
    && git clone https://github.com/arq5x/bedtools2.git \
    && cd bedtools2 && make clean && make all \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bedtools2 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bedtools2/bin' /home/pipeman/.bashrc \
    && echo "alias ngsBedtools='/usr/local/pipeline/bedtools2/bin'" >> /home/pipeman/.bashrc

# + VCFtools: http://vcftools.sourceforge.net/index.html
RUN wget -O /tmp/vcftools_0.1.12a.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12a.tar.gz \
    && tar xzvf /tmp/vcftools_0.1.12a.tar.gz -C /usr/local/pipeline/  \
    && cd /usr/local/pipeline/vcftools_0.1.12a/ && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/vcftools_0.1.12a \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/vcftools_0.1.12a/bin' /home/pipeman/.bashrc \
    && echo "alias ngsVCFtools='/usr/local/pipeline/vcftools_0.1.12a/bin/'" >>  /home/pipeman/.bashrc    
    
#---------------------- FREEBAYES  -------------------------------- 
RUN cd /usr/local/pipeline \
  && git clone --recursive git://github.com/ekg/freebayes.git \
  && cd /usr/local/pipeline/freebayes \
  && make \
  && chmod -R 755 /usr/local/pipeline/freebayes \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/freebayes \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/freebayes/bin' /home/pipeman/.bashrc \
  && echo "alias ngsFreebayes='/usr/local/pipeline/freebayes/bin'" >>  /home/pipeman/.bashrc

#-----------------------------PLATYPUS-------------------------------
RUN wget -O /tmp/Platypus-latest.tgz http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz \
  && tar xvf /tmp/Platypus-latest.tgz -C /usr/local/pipeline/ \
  && cd /usr/local/pipeline/Platypus_0.7.8 \
  && sh ./buildPlatypus.sh \
  && chmod -R 755 /usr/local/pipeline/Platypus_0.7.8 \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/Platypus_0.7.8 \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/Platypus_0.7.8' /home/pipeman/.bashrc \
  && echo "alias ngsPlatypus='/usr/local/pipeline/Platypus_0.7.8'" >>  /home/pipeman/.bashrc

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

