# base image
FROM compbio/debian:v1.0-r002

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

LABEL Description="This is the NGSeasy Big-Kahuna Tool Box" NickName="Big-Kahuna" Version="0.5"

# Remain current and get random dependencies
RUN apt-get update && \
  DEBIAN_FRONTEND=noninteractive \
  apt-get install -y --no-install-recommends \
  man \
  bison \
  flex \
  byacc \
  zlib1g-dev \
  libncurses5-dev && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \

# Create a user:ngseasy and group:ngseasy
  useradd -m -U -s /bin/bash ngseasy && \
  cd /home/ngseasy && \
  usermod -aG sudo ngseasy && \

# make dirs: /usr/local/ngs/bin and sort permissions out
  mkdir /usr/local/ngs && \
  mkdir /usr/local/ngs/bin && \
  chown ngseasy:ngseasy /usr/local/ngs/bin  && \
  chmod -R 777 /usr/local/ngs/bin  && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin && \
  sed -i '$aPATH=$PATH:/usr/local/ngs/bin' /home/ngseasy/.bashrc && \
  bash -c "source /home/ngseasy/.bashrc" && \

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
  rm -r bcftools htslib samtools && \

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
  make install && \
  cd /usr/local/ngs/bin && \
  git clone https://github.com/statgen/bamUtil.git && \
  cd bamUtil && \
  make cloneLib && \
  make all && \
  make install && \
  rm -r /usr/local/ngs/bin/libStatGen && \
  rm -r /usr/local/ngs/bin/bamUtil && \

# samblaster
  cd /usr/local/ngs/bin && \
  git clone git://github.com/GregoryFaust/samblaster.git && \
  cd samblaster && \
  make && \
  cp samblaster /usr/local/bin/ && \
  cd /usr/local/bin/ && \
  rm -r /usr/local/ngs/bin/samblaster && \

# sambamba
  cd /usr/local/ngs/bin && \
  curl -OL https://github.com/lomereiter/sambamba/releases/download/v0.5.1/sambamba_v0.5.1_linux.tar.bz2 && \
  tar -xjvf sambamba_v0.5.1_linux.tar.bz2 && \
  mv sambamba_v0.5.1 sambamba && \
  chmod +rwx sambamba && \
  cp -v sambamba /usr/local/bin/ && \
  cd /usr/local/ngs/bin && \
  rm sambamba_v0.5.1_linux.tar.bz2 && \
  rm /usr/local/ngs/bin/sambamba && \

# seqtk and trimadap
  cd /usr/local/ngs/bin/ && \
  git clone https://github.com/lh3/seqtk.git && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin/seqtk && \
  cd seqtk/ && \
  chmod -R 777 ./* && \
  make && \
  chmod 777 seqtk && \
  chmod 777 trimadap && \
  cp -v seqtk /usr/local/bin/ && \
  cp -v trimadap /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r /usr/local/ngs/bin/seqtk && \

# vcftools https://github.com/vcftools/vcftools
  cd /usr/local/ngs/bin/ && \
  git clone https://github.com/vcftools/vcftools.git && \
  cd vcftools && \
  ./autogen.sh && \
  ./configure && \
  make && \
  make install && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vcftools && \

# vcflib
  cd /usr/local/ngs/bin/ && \
  rm -rfv /usr/local/ngs/bin/vcflib && \
  git clone --recursive git://github.com/ekg/vcflib.git && \
  cd vcflib && \
  make && \
  chmod -R 777 ./bin/ && \
  cp -v ./bin/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vcflib && \

# bedtools
  cd /usr/local/ngs/bin && \
  git clone https://github.com/arq5x/bedtools2.git && \
  cd bedtools2 && \
  make all && \
  chmod -R 777 ./* && \
  make install && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./bedtools2 && \

# vt
  cd /usr/local/ngs/bin && \
  git clone https://github.com/atks/vt.git && \
  chmod -R 777 vt/ && \
  cd vt && \
  make && \
  chmod -R 777 vt && \
  cp -v vt /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vt && \

# vawk
  cd /usr/local/ngs/bin && \
  git clone https://github.com/cc2qe/vawk.git && \
  chmod -R 777 vawk/ && \
  cp -v vawk/vawk /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vawk && \

# bioawk
  cd /usr/local/ngs/bin && \
  git clone https://github.com/lh3/bioawk.git && \
  chmod -R 777 bioawk/ && \
  cd bioawk && \
  make && \
  chmod -R 777 ./ && \
  cp -v bioawk /usr/local/bin && \
  cp -v maketab /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./bioawk && \

# gVCFtools: https://sites.google.com/site/gvcftools/home/strelka-workflow-installation/installation-prerequisites
  GVCFTOOLS_VERSION="0.16" && \
  cd /usr/local/ngs/bin && \
  wget https://sites.google.com/site/gvcftools/home/download/gvcftools-${GVCFTOOLS_VERSION}.tar.gz && \
  tar -xvf gvcftools-${GVCFTOOLS_VERSION}.tar.gz && \
  cd gvcftools-${GVCFTOOLS_VERSION} && \
  make && \
  chmod -R 777 ./ && \
  cp -v ./bin/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./gvcftools-${GVCFTOOLS_VERSION}* && \

# fastc
  FASTQC_VERSION="0.11.4" && \
  cd /usr/local/ngs/bin \
  && wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip \
  && unzip fastqc_v${FASTQC_VERSION}.zip \
  && chmod -R 777 ./FastQC/ \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/FastQC/jbzip2-0.9.jar:/usr/local/ngs/bin/FastQC/sam-1.103.jar' /home/ngseasy/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/FastQC' /home/ngseasy/.bashrc \
  && ln -s /usr/local/ngs/bin/FastQC/fastqc /usr/local/bin/fastqc && \
  cd /usr/local/ngs/bin/ && \
  rm -r fastqc_v${FASTQC_VERSION}.zip && \

# Trimmomatic
  TRIMMOMATIC_VERSION="0.32" && \
  cd /usr/local/ngs/bin \
  && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
  && unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip \
  && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/Trimmomatic-0.32/trimmomatic-0.32.jar' /home/ngseasy/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/Trimmomatic-0.32' /home/ngseasy/.bashrc \
  && cp -v /usr/local/ngs/bin/Trimmomatic-0.32/trimmomatic-0.32.jar /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \

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

# Platypus
  cd /usr/local/ngs/bin && \
  git clone --recursive https://github.com/andyrimmer/Platypus.git && \
  chmod -R 777 Platypus && \
  cd Platypus && \
  make && \
  chmod -R 777 ./* && \
  cp -vrf ./bin/* /usr/local/bin && \

# scalpel
  SCALPEL_VERSION="0.5.2" && \
  cd /usr/local/ngs/bin && \
  http://sourceforge.net/projects/scalpel/files/scalpel-${SCALPEL_VERSION}.tar.gz && \
  tar -xvf scalpel-${SCALPEL_VERSION}.tar.gz && \
  cd scalpel-${SCALPEL_VERSION} && \
  make && \
  chmod -R 777 ./* && \
  sed -i '$aPATH=$PATH:/usr/local/ngs/bin/scalpel-0.5.2/' /home/ngseasy/.bashrc && \
  rm scalpel-${SCALPEL_VERSION}.tar.gz && \

# varscan
  VARSCAN_VERSION="v2.3.9" && \
  cd /usr/local/ngs/bin && \
  wget http://sourceforge.net/projects/varscan/files/VarScan.${VARSCAN_VERSION}.jar && \
  chmod -R 777 VarScan.${VARSCAN_VERSION}.jar && \
  ln -s VarScan.${VARSCAN_VERSION}.jar VarScan.jar && \

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
  BCBIO_VAR_VERSION="0.2.6" && \
  wget https://github.com/chapmanb/bcbio.variation/releases/download/v0.2.6/bcbio.variation-${BCBIO_VAR_VERSION}-standalone.jar && \
  chmod -R 777 bcbio.variation-${BCBIO_VAR_VERSION}--standalone.jar && \
  ln -s bcbio.variation-${BCBIO_VAR_VERSION}--standalone.jar /usr/local/bin/bcbio.variation.jar && \

# CNVkit
# https://github.com/chapmanb/cnvkit and https://github.com/etal/cnvkit
  cd  /usr/local/ngs/bin/ && \
  CNVKIT_VERSION="0.7.3" && \
  pip install pandas biopython pysam pyvcf --upgrade && \
  Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite()' && \
  Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("PSCBS", "cghFLasso")' && \
  pip install cnvkit==${CNVKIT_VERSION} && \

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
  cp -v ./svtyper/svtyper /usr/local/bin/ && \

# source .bashrc
  bash -c "source /home/ngseasy/.bashrc" && \

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
