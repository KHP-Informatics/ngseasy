# base image
FROM compbio/debian:r1.0-002
# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com
LABEL Description="This is the base image for all ngseasy tools images; Contains SAM/BAM/VCF/BED Parsers" Version="r1.0-002"
# Remain current
RUN apt-get update && \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install -y \
# needed for htslib build
  zlib1g-dev \
  libncurses5-dev

# Create a user:ngseasy and group:ngseasy
RUN useradd -m -U -s /bin/bash ngseasy && \
  cd /home/ngseasy && \
  usermod -aG sudo ngseasy && \

# make /usr/local/ngs/bin install dirs and sort permissions out
  mkdir /usr/local/ngs && \
  mkdir /usr/local/ngs/bin && \
  chown ngseasy:ngseasy /usr/local/ngs/bin  && \
  chmod -R 777 /usr/local/ngs/bin  && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin

# STANDARD NGS TOOLS
# Tools used for processing SAM/BAM/BED/VCF files
# samtools,htslib,bcftools,parallel,bamUtil,sambamba,samblaster,vcftools,vcflib,seqtk,bedtools2,libStatGen

# samtools, htslib, bcftools
RUN cd /usr/local/ngs/bin && \
  git clone --branch=develop git://github.com/samtools/htslib.git && \
  git clone --branch=develop git://github.com/samtools/bcftools.git && \
  git clone --branch=develop git://github.com/samtools/samtools.git && \
  cd /usr/local/ngs/bin/htslib && \
  autoconf && \
  ./configure  && \
  make && \
  make install && \
  cd /usr/local/ngs/bin/bcftools && \
  make && \
  make install && \
  cd /usr/local/ngs/bin/samtools && \
  make && \
  make install && \

# parallel
  cd /usr/local/ngs/bin && \
  wget http://ftpmirror.gnu.org/parallel/parallel-20140222.tar.bz2 && \
  bzip2 -dc parallel-20140222.tar.bz2 | tar xvf - && \
  cd parallel-20140222 && \
  ./configure && \
  make && \
  make install && \

# libStatGen and bamUtil
  cd /usr/local/ngs/bin && \
  git clone https://github.com/statgen/libStatGen.git && \
  cd libStatGen && \
  make all && \
  cd /usr/local/ngs/bin && \
  git clone https://github.com/statgen/bamUtil.git && \
  cd bamUtil && \
  make cloneLib && \
  make all && \
  make install && \

# samblaster and sambamba
  cd /usr/local/ngs/bin && \
  git clone git://github.com/GregoryFaust/samblaster.git && \
  cd samblaster && \
  make && \
  cp samblaster /usr/local/bin/ && \
  cd /usr/local/ngs/bin && \
  curl -OL https://github.com/lomereiter/sambamba/releases/download/v0.5.1/sambamba_v0.5.1_linux.tar.bz2 && \
  tar -xjvf sambamba_v0.5.1_linux.tar.bz2 && \
  mv sambamba_v0.5.1 sambamba && \
  chmod +rwx sambamba && \
  cp sambamba /usr/local/bin/ && \

# seqtk and trimadap
  cd /usr/local/ngs/bin/ && \
  git clone https://github.com/lh3/seqtk.git && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin/seqtk && \
  cd seqtk/ && \
  chmod -R 777 ./* && \
  make && \
  cp -v seqtk /usr/local/bin/ && \
  cp -v trimadap /usr/local/bin/ && \

# vcftools
  cd /usr/local/ngs/bin/ && \
  wget -O /tmp/vcftools_0.1.12b.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12b.tar.gz && \
  tar xzvf /tmp/vcftools_0.1.12b.tar.gz -C /usr/local/ngs/bin/  && \
  export PERL5LIB=/usr/local/ngs/bin/vcftools_0.1.12b/perl/  && \
  cd /usr/local/ngs/bin/vcftools_0.1.12b/ && \
  make && \
  cp -vrf /usr/local/ngs/bin/vcftools_0.1.12b/bin/*  /usr/local/bin/ && \

# vcflib
  cd /usr/local/ngs/bin/ && \
  rm -rfv /usr/local/ngs/bin/vcflib && \
  git clone --recursive git://github.com/ekg/vcflib.git && \
  cd vcflib && \
  make && \
  cp ./bin/* /usr/local/bin/ && \

# bedtools
  cd /usr/local/ngs/bin && \
  git clone https://github.com/arq5x/bedtools2.git && \
  cd bedtools2 && \
  make clean && \
  make all && \
  make install && \

# vt
  cd /usr/local/ngs/bin && \
  git clone https://github.com/atks/vt.git && \
  chmod -R 777 vt/ && \
  cd vt && \
  make && \
  cp -v vt /usr/local/bin && \

# vawk
  cd /usr/local/ngs/bin && \
  git clone https://github.com/cc2qe/vawk.git && \
  chmod -R 777 vawk/ && \
  cp -v vawk/vawk /usr/local/bin && \
  apt-get install -y bison flex byacc && \

# bioawk
  cd /usr/local/ngs/bin && \
  git clone https://github.com/lh3/bioawk.git && \
  chmod -R 777 bioawk/ && \
  cd bioawk && \
  make && \
  cp -v bioawk /usr/local/bin && \
  cp -v maketab /usr/local/bin && \

# Clean up APT when done.
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/ && \
  rm -rf /usr/local/ngs/bin/*

# PERMISSIONS
RUN chmod -R 777 /usr/local/ngs/bin && chown -R ngseasy:ngseasy /usr/local/ngs/bin

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

USER ngseasy
WORKDIR /home/ngseasy
