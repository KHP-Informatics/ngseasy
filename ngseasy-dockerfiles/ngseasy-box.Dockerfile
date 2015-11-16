# base image
FROM compbio/debian:r1.0-002

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

LABEL Description="This is the NGSeasy Big-Kahuna Tool Box" NickName="Big-Kahuna" Version="0.5"

# Remain current
RUN apt-get update && \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install -y --no-install-recommends \
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
  chown -R ngseasy:ngseasy /usr/local/ngs/bin && \

# NGSeasy Tools
# samtools, htslib, bcftools
  cd /usr/local/ngs/bin && \
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

# parallel (this is now in little-fatty-deb)
#  cd /usr/local/ngs/bin && \
#  wget http://ftpmirror.gnu.org/parallel/parallel-20140222.tar.bz2 && \
#  bzip2 -dc parallel-20140222.tar.bz2 | tar xvf - && \
#  cd parallel-20140222 && \
#  ./configure && \
#  make && \
#  make install && \

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

# samblaster
  cd /usr/local/ngs/bin && \
  git clone git://github.com/GregoryFaust/samblaster.git && \
  cd samblaster && \
  make && \
  cp samblaster /usr/local/bin/ && \

# sambamba
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

# gVCFtools
  cd /usr/local/ngs/bin && \
  git clone --recursive https://github.com/sequencing/gvcftools.git && \
  cd gvcftools && \
  make && \
  sed -i '$aPATH=$PATH:/usr/local/ngs/bin/gvcftools' /home/ngseasy/.bashrc && \
  ln -s /usr/local/ngs/bin/gvcftools/* /usr/local/bin/ && \

# fastc
  cd /usr/local/ngs/bin \
  && wget -O /tmp/fastqc_v0.11.2.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip \
  && unzip /tmp/fastqc_v0.11.2.zip -d /usr/local/ngs/bin/ \
  && chmod -R 766 /usr/local/ngs/bin/ \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/FastQC/jbzip2-0.9.jar:/usr/local/ngs/bin/FastQC/sam-1.103.jar' /home/ngseasy/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/FastQC' /home/ngseasy/.bashrc \
  && ln -s /usr/local/ngs/bin/FastQC/fastqc /usr/local/bin/fastqc && \

# Trimmomatic
  cd /usr/local/ngs/bin \
  && wget -O /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip \
  && unzip /tmp/Trimmomatic-0.32.zip -d /usr/local/ngs/bin/ \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/Trimmomatic-0.32/trimmomatic-0.32.jar' /home/ngseasy/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/Trimmomatic-0.32' /home/ngseasy/.bashrc \
  && cp -v /usr/local/ngs/bin/Trimmomatic-0.32/trimmomatic-0.32.jar /usr/local/bin && \

# Picard
  cd /usr/local/ngs/bin \
  && wget -O /tmp/picard-tools-1.129.zip https://github.com/broadinstitute/picard/releases/download/1.129/picard-tools-1.129.zip \
  && mkdir /usr/local/ngs/bin/picardtools \
  && unzip /tmp/picard-tools-1.129.zip -d /usr/local/ngs/bin/picardtools/ \
  && chown -R ngseasy:ngseasy /usr/local/ngs/bin/picardtools \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/picardtools/picard-tools-1.129/snappy-java-1.0.3-rc3.jar' /home/ngseasy/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/picardtools/picard-tools-1.129' /home/ngseasy/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/picardtools/picard-tools-1.129' ~/.bashrc && \

# bwa
  cd /usr/local/ngs/bin \
  && wget -O /tmp/bwa-0.7.12.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2 \
  && tar xjvf /tmp/bwa-0.7.12.tar.bz2 -C /usr/local/ngs/bin/ \
  && chmod -R 777 /usr/local/ngs/bin \
  && cd /usr/local/ngs/bin/bwa-0.7.12 && make \
  && cp -v /usr/local/ngs/bin/bwa-0.7.12/bwa /usr/local/bin && \

# snap
  cd /usr/local/ngs/bin && \
  git clone https://github.com/amplab/snap.git && \
  chmod -R 777 snap && \
  cd snap && \
  make && \
  chmod -R 777 /usr/local/ngs/bin && \
  cp -v snap* /usr/local/bin && \

# freebayes
  cd /usr/local/ngs/bin \
  && git clone --recursive git://github.com/ekg/freebayes.git \
  && cd /usr/local/ngs/bin/freebayes \
  && make \
  && chmod -R 777 /usr/local/ngs/bin/freebayes \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/freebayes/bin' /home/ngseasy/.bashrc \
  && cp -v /usr/local/ngs/bin/freebayes/bin/* /usr/local/bin && \

# get python (now in little-fatty-deb)
#  apt-get install -y python-dev && \
# cython for platypus
#  cd /tmp && \
#  wget https://github.com/cython/cython/archive/0.22.1.tar.gz && \
#  chmod 777 0.22.1.tar.gz && \
#  tar xvf 0.22.1.tar.gz && \
#  cd cython-0.22.1 && \
#  python setup.py install && \

# Platypus
  cd /usr/local/ngs/bin && \
  git clone --recursive https://github.com/andyrimmer/Platypus.git && \
  chmod -R 777 Platypus && \
  cd Platypus && \
  make && \
  chmod -R 777 ./* && \
  cp -vrf ./bin/* /usr/local/bin && \

# VarDict
  cd /usr/local/ngs/bin/ && \
  git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git && \
  cd /usr/local/ngs/bin/VarDictJava && \
  ./gradlew clean installApp && \
  chmod -R 777 /usr/local/ngs/bin/VarDictJava && \
  sed  -i '$aPATH=$PATH:/usr/local/ngs/bin/VarDictJava/VarDict' /home/ngseasy/.bashrc && \

# ABRA - Assembly Based ReAligner https://github.com/mozack/abra
  cd /usr/local/ngs/bin && \
  wget https://github.com/mozack/abra/releases/download/v0.94/abra-0.94-SNAPSHOT-jar-with-dependencies.jar && \
  chmod 775 abra-0.94-SNAPSHOT-jar-with-dependencies.jar && \
  mv -v abra-0.94-SNAPSHOT-jar-with-dependencies.jar /usr/local/bin/abra-0.94 && \

# glia
  cd /usr/local/ngs/bin && \
  git clone --recursive https://github.com/ekg/glia.git && \
  chmod -R 777 ./glia && \
  cd ./glia && \
  make && \
  cp -v ./glia /usr/local/bin/ && \

# ogap
  cd /usr/local/ngs/bin/ && \
  git clone --recursive https://github.com/ekg/ogap.git && \
  cd ogap && \
  make all && \
  chmod -R 777 ./* && \
  cp -v ogap /usr/local/bin/ && \

# lein for chapmanb/bcbio.variation
  cd  /usr/local/ngs/bin/ && \
  wget -O lein https://raw.githubusercontent.com/technomancy/leiningen/stable/bin/lein && \
  chmod 777 lein && \
  cp -v lein /usr/local/bin && \
  lein && \

# chapmanb/bcbio.variation
  cd  /usr/local/ngs/bin/ && \
  wget https://github.com/chapmanb/bcbio.variation/releases/download/v0.2.6/bcbio.variation-0.2.6-standalone.jar && \
  chmod -R 777 bcbio.variation-0.2.6-standalone.jar && \
  ln -s bcbio.variation-0.2.6-standalone.jar /usr/local/bin/bcbio.variation.jar && \

# CNVkit
# https://github.com/chapmanb/cnvkit and https://github.com/etal/cnvkit
  cd  /usr/local/ngs/bin/ && \
  pip install pandas biopython pysam pyvcf --upgrade && \
  Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite()' && \
  Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("PSCBS", "cghFLasso")' && \
  pip install cnvkit==0.7.3 && \

# lumpy-sv
  cd  /usr/local/ngs/bin/ && \
  git clone --recursive https://github.com/arq5x/lumpy-sv.git && \
  cd lumpy-sv && \
  make && \
  chmod -R 777 ./bin/ && \
  cp -v bin/* /usr/local/bin/ && \

# bamkit
  cd  /usr/local/ngs/bin/ && \
  git clone --recursive https://github.com/hall-lab/bamkit.git && \
  chmod -R 777 ./bamkit/ && \
  cp -v ./bamkit/bam* /usr/local/bin/ && \
  cp -v ./bamkit/sectosupp /usr/local/bin/ && \

# svtyper
  cd  /usr/local/ngs/bin/ && \
  git clone --recursive  https://github.com/hall-lab/svtyper.git && \
  chmod -R 777 ./svtyper && \
  cp -v ./svtyper/svtyper /usr/loca/bin/ && \

# Clean up APT when done.
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/ && \
  rm -rf /usr/local/ngs/bin/*

ADD fix_ambiguous /usr/local/bin/

# PERMISSIONS
RUN chmod -R 777 /usr/local/ngs/bin && chown -R ngseasy:ngseasy /usr/local/ngs/bin

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

# run as user ngseasy
HOME /home/ngseasy
USER ngseasy
WORKDIR /home/ngseasy
