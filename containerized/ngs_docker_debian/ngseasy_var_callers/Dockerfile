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

#---------------------------------GATK-------------------------------------#
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
 
#---------------------------------sam/vcf processing-------------------------------------#
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

# + VCFtools: http://vcftools.sourceforge.net/index.html
RUN wget -O /tmp/vcftools_0.1.12a.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12a.tar.gz \
    && tar xzvf /tmp/vcftools_0.1.12a.tar.gz -C /usr/local/pipeline/  \
    && cd /usr/local/pipeline/vcftools_0.1.12a/ && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/vcftools_0.1.12a \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/vcftools_0.1.12a/bin' /home/pipeman/.bashrc \
    && echo "alias ngsVCFtools='/usr/local/pipeline/vcftools_0.1.12a/bin/'" >>  /home/pipeman/.bashrc

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

