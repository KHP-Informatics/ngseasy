# Base image will be Trusty
FROM snewhouse/base-config

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com

# Set correct environment variables.
ENV HOME /root

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

#------------------ Variant Annotation & Reports -------------------------

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
    
# + SnpEff
RUN wget -O /tmp/snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip \
    && unzip /tmp/snpEff_latest_core.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/snpEff \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/snpEff/snpEff.jar' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/snpEff' /home/pipeman/.bashrc \
    && echo "alias ngsSNPeff='/usr/local/pipeline/snpEff'" >> /home/pipeman/.bashrc

# + ANNOVAR (see licence, and registration)
# Available on reg:  http://www.openbioinformatics.org/annovar/download/Ht8qRwQSTi/annovar.latest.tar.gz
RUN wget -O /tmp/annovar.latest.tar.gz http://www.openbioinformatics.org/annovar/download/mP628pfL21/annovar.latest.tar.gz \
  && tar xzvf /tmp/annovar.latest.tar.gz -C /usr/local/pipeline/ \
  && cd /usr/local/pipeline/annovar \
  && cp -v /usr/local/pipeline/annovar/table_annovar.pl table_annovar.pl.original \
  && wget http://www.openbioinformatics.org/annovar/download/table_annovar.pl \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/annovar \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/annovar' /home/pipeman/.bashrc \
   && echo "alias ngsSNPeff='/usr/local/pipeline/annovar'" >> /home/pipeman/.bashrc

   # + ENSEMBL API
RUN cd \
  && mkdir src \
  && cd src \
  && wget ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz \
  && wget http://bioperl.org/DIST/old_releases/bioperl-1.2.3.tar.gz \
  && tar zxvf ensembl-api.tar.gz \
  && tar zxvf bioperl-1.2.3.tar.gz \
  && cd \
  && chmod -R 766 /src

ENV PERL5LIB ${PERL5LIB}:${HOME}/src/bioperl-1.2.3
ENV PERL5LIB ${PERL5LIB}:${HOME}/src/ensembl/modules
ENV PERL5LIB ${PERL5LIB}:${HOME}/src/ensembl-compara/modules
ENV PERL5LIB ${PERL5LIB}:${HOME}/src/ensembl-variation/modules
ENV PERL5LIB ${PERL5LIB}:${HOME}/src/ensembl-funcgen/modules

# + VEP
# Download latest version (75)
# http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html
RUN wget -O /tmp/75.zip https://github.com/Ensembl/ensembl-tools/archive/release/75.zip \
  && unzip /tmp/75.zip -d /usr/local/pipeline/ 
#  && cd ensembl-tools-release-75/scripts/variant_effect_predictor/
#  && perl INSTALL.pl -a acf -t -s homo_sapiens


#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***

# CleanRUN chmod -R 766 /usr/local/pipeline/***up the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80
