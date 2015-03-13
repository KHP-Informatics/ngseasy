# NGSeasy Base Image

# base image
# FROM debian:jessie
FROM compbio/debian-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Remain current
RUN apt-get update &&  apt-get upgrade -y && apt-get dist-upgrade -y && apt-get install -y ldc

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
    autoconf && \
    ./configure  && \
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

# libStatGen and bamUtil
RUN cd /usr/local/pipeline && \
    git clone https://github.com/statgen/libStatGen.git && \
    cd libStatGen && \
    make all && \
    cd /usr/local/pipeline && \
    git clone https://github.com/statgen/bamUtil.git && \
    cd bamUtil && \
    make cloneLib && \
    make all && \
    make install

# samblaster and sambamba
RUN cd /usr/local/pipeline && \ 
    git clone git://github.com/GregoryFaust/samblaster.git && \ 
    cd samblaster && \ 
    make && \ 
    cp samblaster /usr/local/bin/ && \
    cd /usr/local/pipeline && \ 
    curl -OL https://github.com/lomereiter/sambamba/releases/download/v0.5.1/sambamba_v0.5.1_linux.tar.bz2 && \ 
    tar -xjvf sambamba_v0.5.1_linux.tar.bz2 && \ 
    mv sambamba_v0.5.1 sambamba && \
    chmod +rwx sambamba && \
    cp sambamba /usr/local/bin/

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

# ogap  and bamleftalign  
RUN cd /usr/local/pipeline/ && \
    git clone --recursive https://github.com/ekg/ogap.git && \
    cd ogap && \
    make all && \
    chmod -R 777 ./* && \
    cp -v ogap /usr/local/bin/ && \
    cd /usr/local/pipeline/ && \
    git clone --recursive git://github.com/ekg/freebayes.git && \
    cd freebayes && \
    make all && \
    chmod -R 777 ./* && \
    cp bin/bamleftalign /usr/local/bin/ && \
    rm -frv /usr/local/pipeline/freebayes

# vcftools and vcflib and bedtools2 and vt
RUN cd /usr/local/pipeline/ && \
    wget -O /tmp/vcftools_0.1.12b.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12b.tar.gz && \
    tar xzvf /tmp/vcftools_0.1.12b.tar.gz -C /usr/local/pipeline/  && \
    export PERL5LIB=/usr/local/pipeline/vcftools_0.1.12b/perl/  && \
    cd /usr/local/pipeline/vcftools_0.1.12b/ && \
    make && \
    cp -vrf /usr/local/pipeline/vcftools_0.1.12b/bin/*  /usr/local/bin/ && \
    cd /usr/local/pipeline/ && \
    rm -rfv /usr/local/pipeline/vcflib && \
    git clone --recursive git://github.com/ekg/vcflib.git && \
    cd vcflib && \
    make && \
    cp ./bin/* /usr/local/bin/ && \
    cd /usr/local/pipeline && \
    git clone https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && \
    make clean && \
    make all && \
    make install && \
    cd /usr/local/pipeline && \
    git clone https://github.com/atks/vt.git && \
    chmod -R 777 vt/ && \
    cd vt && \
    make && \
    cp -v vt /usr/local/bin 

# vawk and bioawk
RUN cd /usr/local/pipeline && \
    git clone https://github.com/cc2qe/vawk.git && \
    chmod -R 777 vawk/ && \
    cp -v vawk/vawk /usr/local/bin && \
    apt-get install -y bison flex byacc && \
    && \cd /usr/local/pipeline && \
    git clome https://github.com/lh3/bioawk.git && \
    chmod -R 777 bioawk/ && \
    cd bioawk && \
    make && \
    cp -v vawk/vawk /usr/local/bin

  
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
    rm -rf /var/lib/{apt,dpkg,cache,log}/ && \
    rm -rf /usr/local/pipeline/*



