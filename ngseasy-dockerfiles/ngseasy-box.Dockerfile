# base image
FROM compbio/ngseasy-debian:v1.0

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

LABEL Description="This is the NGSeasy Tool Box" NickName="Big-Kahuna"

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
  apt-get purge

## NGS Tools -------------------------------------------------------------------
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
  cd /usr/local/ngs/bin && \
  rm -r bcftools && \
  rm -r htslib && \
  rm -r samtools && \

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
  SAMBAMBA_VERSION="v0.5.9" && \
  curl -OL https://github.com/lomereiter/sambamba/releases/download/${SAMBAMBA_VERSION}/sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  tar -xjvf sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  mv sambamba_${SAMBAMBA_VERSION} sambamba && \
  chmod +rwx sambamba && \
  cp -v sambamba /usr/local/bin/ && \
  cd /usr/local/ngs/bin && \
  rm sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  rm /usr/local/ngs/bin/sambamba && \

# seqtk and trimadap
  cd /usr/local/ngs/bin/ && \
  git clone https://github.com/lh3/seqtk.git && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin/seqtk && \
  cd seqtk/ && \
  chmod -R 755 ./* && \
  make && \
  chmod 755 seqtk && \
  chmod 755 trimadap && \
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
  chmod -R 755 ./bin/ && \
  cp -v ./bin/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vcflib && \

# bedtools
  cd /usr/local/ngs/bin && \
  git clone https://github.com/arq5x/bedtools2.git && \
  cd bedtools2 && \
  make all && \
  chmod -R 755 ./* && \
  make install && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./bedtools2 && \

# vt
  cd /usr/local/ngs/bin && \
  git clone https://github.com/atks/vt.git && \
  chmod -R 755 vt/ && \
  cd vt && \
  make && \
  chmod -R 755 vt && \
  cp -v vt /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vt && \

# vawk
  cd /usr/local/ngs/bin && \
  git clone https://github.com/cc2qe/vawk.git && \
  chmod -R 755 vawk/ && \
  cp -v vawk/vawk /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vawk && \

# bioawk
  cd /usr/local/ngs/bin && \
  git clone https://github.com/lh3/bioawk.git && \
  chmod -R 755 bioawk/ && \
  cd bioawk && \
  make && \
  chmod -R 755 ./ && \
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
  chmod -R 755 ./ && \
  cp -v ./bin/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./gvcftools-${GVCFTOOLS_VERSION}* && \

# fastc
  FASTQC_VERSION="0.11.4" && \
  cd /usr/local/ngs/bin && \
  wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip && \
  unzip fastqc_v${FASTQC_VERSION}.zip && \
  chmod -R 755 ./FastQC/ && \
  sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/FastQC/jbzip2-0.9.jar:/usr/local/ngs/bin/FastQC/sam-1.103.jar' /home/ngseasy/.bashrc && \
  sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/FastQC' /home/ngseasy/.bashrc && \
  ln -s /usr/local/ngs/bin/FastQC/fastqc /usr/local/bin/fastqc && \
  cd /usr/local/ngs/bin/ && \
  rm -r fastqc_v${FASTQC_VERSION}.zip && \

# Trimmomatic
  TRIMMOMATIC_VERSION="0.32" && \
  cd /usr/local/ngs/bin && \
  wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \
  unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \
  chmod -R 755 Trimmomatic-${TRIMMOMATIC_VERSION}/* && \
  sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/Trimmomatic-0.32/trimmomatic-0.32.jar' /home/ngseasy/.bashrc && \
  sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/Trimmomatic-0.32' /home/ngseasy/.bashrc && \
  cp -v /usr/local/ngs/bin/Trimmomatic-0.32/trimmomatic-0.32.jar /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \

# skewer https://github.com/relipmoc/skewer
  cd /usr/local/ngs/bin && \
  git clone https://github.com/relipmoc/skewer.git && \
  cd skewer && \
  make && \
  chmod -R 755 ./* && \
  make install && \
  cd /usr/local/ngs/bin/ && \

# Picard
  PICARD_VERSION="1.141" && \
  cd /usr/local/ngs/bin && \
  wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip && \
  unzip picard-tools-${PICARD_VERSION}.zip && \
  chmod -R 755 /usr/local/ngs/bin/picard-tools-${PICARD_VERSION} && \
  mv -v /usr/local/ngs/bin/picard-tools-${PICARD_VERSION}/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r picard-tools-${PICARD_VERSION} && \
  rm -r picard-tools-${PICARD_VERSION}.zip && \

# bamkit
  cd  /usr/local/ngs/bin/ && \
  git clone --recursive https://github.com/hall-lab/bamkit.git && \
  chmod -R 755 ./bamkit/ && \
  cp -v ./bamkit/bam* /usr/local/bin/ && \
  cp -v ./bamkit/sectosupp /usr/local/bin/ && \
  rm -r ./bamkit/

## Aligners --------------------------------------------------------------------
# bwakit
# http://sourceforge.net/projects/bio-bwa/files/bwakit/
RUN cd /usr/local/ngs/bin && \
  BWA_VERSION="0.7.12" && \
  wget \
  http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-${BWA_VERSION}_x64-linux.tar.bz2 && \
  tar xjvf bwakit-${BWA_VERSION}_x64-linux.tar.bz2 && \
  chmod -R 755 /usr/local/ngs/bin && \
  ln -s /usr/local/ngs/bin/bwa.kit/bwa /usr/local/bin/bwa && \
  cd /usr/local/ngs/bin && \
  rm bwakit-${BWA_VERSION}_x64-linux.tar.bz2 && \

# snap v1.0beta.18 snap-aligner
  cd /usr/local/ngs/bin && \
  git clone https://github.com/amplab/snap.git && \
  chmod -R 755 snap && \
  cd snap && \
  make && \
  chmod -R 755 /usr/local/ngs/bin/snap && \
  cp -v snap-aligner /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r snap && \

# novoalign: need to add novoalign.lic
  NOVOALIGN_VERSION="V3.03.02" && \
  cd /usr/local/ngs/bin && \
  wget https://s3-eu-west-1.amazonaws.com/novoalign/novocraft${NOVOALIGN_VERSION}.Linux3.0.tar.gz && \
  tar xvf novocraft${NOVOALIGN_VERSION}.Linux3.0.tar.gz && \
  chmod -R 755 novocraft && \
  cp -vr novocraft/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r novocraft && \
  rm novocraft${NOVOALIGN_VERSION}.Linux3.0.tar.gz && \

# bowtie2
  BOWTIE2_VERSION="2.2.6" && \
  cd /usr/local/ngs/bin && \
  wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
  unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
  cd bowtie2-${BOWTIE2_VERSION} && \
  chmod -R 755 ./ && \
  cp -v bowtie2* /usr/local/bin/ && \
  cp -r scripts /usr/local/ngs/bin/bowtie2_scripts && \
  cd /usr/local/ngs/bin/ && \
  rm -r bowtie2-${BOWTIE2_VERSION} && \
  rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \

# stampy
  STAMPY_VERSION="1.0.28" && \
  cd /usr/local/ngs/bin && \
  wget http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz && \
  tar -xvf Stampy-latest.tgz && \
  chmod -R 755 stampy-${STAMPY_VERSION} && \
  cd stampy-${STAMPY_VERSION} && \
  mv makefile Makefile && \
  make && \
  export PATH=$PATH:/usr/local/ngs/bin/stampy-${STAMPY_VERSION} && \
  cd /usr/local/ngs/bin/ && \
  rm -r Stampy-latest.tgz && \

# mrsfast
# https://github.com/sfu-compbio/mrsfast/
  cd /usr/local/ngs/bin && \
  git clone https://github.com/sfu-compbio/mrsfast && \
  chmod -R 755 mrsfast && \
  cd mrsfast && \
  make && \
  chmod 755 snp_indexer && \
  chmod 755 mrsfast && \
  cp -v snp_indexer /usr/local/bin && \
  cp -v mrsfast /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r mrsfast

## Variant Calling -------------------------------------------------------------
# freebayes
RUN cd /usr/local/ngs/bin && \
  git clone --recursive git://github.com/ekg/freebayes.git && \
  cd /usr/local/ngs/bin/freebayes && \
  make && \
  chmod -R 755 /usr/local/ngs/bin/freebayes && \
  sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/freebayes/bin' /home/ngseasy/.bashrc && \
  cp -v /usr/local/ngs/bin/freebayes/bin/* /usr/local/bin && \
  cp -v /usr/local/ngs/bin/freebayes/scripts/* /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r freebayes && \

# Platypus VERSION := 0.8.1
  cd /usr/local/ngs/bin && \
  git clone --recursive https://github.com/andyrimmer/Platypus.git && \
  chmod -R 755 Platypus && \
  cd Platypus && \
  make && \
  chmod -R 755 ./* && \
  cp -vrf ./bin/* /usr/local/bin && \
  cp -vrf ./scripts/* /usr/local/bin && \
  cp -v ./*.py /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r Platypus && \

# scalpel
  SCALPEL_VERSION="0.5.2" && \
  cd /usr/local/ngs/bin && \
  wget http://sourceforge.net/projects/scalpel/files/scalpel-${SCALPEL_VERSION}.tar.gz && \
  tar -xvf scalpel-${SCALPEL_VERSION}.tar.gz && \
  cd scalpel-${SCALPEL_VERSION} && \
  make && \
  chmod -R 755 ./* && \
  sed -i '$aPATH=$PATH:/usr/local/ngs/bin/scalpel-0.5.2/' /home/ngseasy/.bashrc && \
  cd /usr/local/ngs/bin/ && \
  rm scalpel-${SCALPEL_VERSION}.tar.gz && \

# varscan
  VARSCAN_VERSION="v2.3.9" && \
  cd /usr/local/ngs/bin && \
  wget http://sourceforge.net/projects/varscan/files/VarScan.${VARSCAN_VERSION}.jar && \
  chmod -R 755 VarScan.${VARSCAN_VERSION}.jar && \
  ln -s VarScan.${VARSCAN_VERSION}.jar varscan.jar && \

# VarDict cmd: VarDict
  cd /usr/local/ngs/bin/ && \
  export JAVA_HOME="/usr/lib/jvm/java-1.7.0-openjdk-amd64" && \
  git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git && \
  cd /usr/local/ngs/bin/VarDictJava && \
  ./gradlew clean installApp && \
  chmod -R 755 /usr/local/ngs/bin/VarDictJava && \
  sed -i '$aPATH=$PATH:/usr/local/ngs/bin/VarDictJava/VarDict' /home/ngseasy/.bashrc && \
  sed -i '$aPATH=$PATH:/usr/local/ngs/bin/VarDictJava/build/install/VarDict/bin' /home/ngseasy/.bashrc && \
  sed -i '$aPATH=$PATH:/usr/local/ngs/bin/VarDictJava/build/install/VarDict/lib' /home/ngseasy/.bashrc && \

## chapmanb/bcbio.variation ----------------------------------------------------
# lein for chapmanb/bcbio.variation
  cd  /usr/local/ngs/bin/ && \
  wget -O lein https://raw.githubusercontent.com/technomancy/leiningen/stable/bin/lein && \
  chmod 755 lein && \
  cp -v lein /usr/local/bin && \
  export LEIN_ROOT=yes && \
  lein && \

# chapmanb/bcbio.variation
  cd /usr/local/ngs/bin/ && \
  BCBIO_VAR_VERSION="0.2.6" && \
  wget https://github.com/chapmanb/bcbio.variation/releases/download/v0.2.6/bcbio.variation-${BCBIO_VAR_VERSION}-standalone.jar && \
  chmod -R 755 bcbio.variation-${BCBIO_VAR_VERSION}-standalone.jar && \
  ln -s bcbio.variation-${BCBIO_VAR_VERSION}-standalone.jar bcbio.variation.jar && \

##  bcbio.variation.recall (https://github.com/chapmanb/bcbio.variation.recall)
# GATK_FRAMEWORK, io_lib, scramble, staden for bcbio.variation.recall.
# http://sourceforge.net/p/staden/discussion/347718/thread/88e0e73b/
  cd /usr/local/ngs/bin/ && \
  GATK_FRAMEWORK_VERSION="3.4-46" && \
  wget https://github.com/chapmanb/gatk/releases/download/v${GATK_FRAMEWORK_VERSION}-framework/gatk-framework-${GATK_FRAMEWORK_VERSION}.tar.gz && \
  tar -xvf gatk-framework-${GATK_FRAMEWORK_VERSION}.tar.gz && \
  chmod -R 755 gatk-framework-${GATK_FRAMEWORK_VERSION}/* && \
  cp -v gatk-framework-${GATK_FRAMEWORK_VERSION}/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r gatk-framework-${GATK_FRAMEWORK_VERSION}.tar.gz && \
  rm -r gatk-framework-${GATK_FRAMEWORK_VERSION} && \
  apt-get install -y --no-install-recommends libstaden-read-dev staden-io-lib-utils staden staden-common && \
  IO_LIB_VERSION="1.14.6" && \
  wget http://sourceforge.net/projects/staden/files/io_lib/${IO_LIB_VERSION}/io_lib-${IO_LIB_VERSION}.tar.gz && \
  tar -xvf io_lib-${IO_LIB_VERSION}.tar.gz && \
  cd io_lib-${IO_LIB_VERSION} && \
  chmod -R 755 ./* && \
  ./configure && \
  make && \
  make install && \
  rm -r /usr/local/ngs/bin/io_lib-${IO_LIB_VERSION}.tar.gz && \
  cd /usr/local/ngs/bin/ && \
  git clone --recursive https://github.com/chapmanb/bcbio.variation.recall.git && \
  cd bcbio.variation.recall && \
  chmod -R 755 ./* && \
  export LEIN_ROOT=yes && \
  make && \
  chmod -R 755 ./* && \
  cp -v /usr/local/ngs/bin/bcbio.variation.recall/bin/* /usr/local/bin

## SV CALLERS ------------------------------------------------------------------
# CNVkit
# https://github.com/chapmanb/cnvkit and https://github.com/etal/cnvkit
RUN cd  /usr/local/ngs/bin/ && \
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
  chmod -R 755 ./* && \
  cp -v bin/* /usr/local/bin/ && \
  cp -vr scripts/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r lumpy-sv && \

# svtyper
  cd  /usr/local/ngs/bin/ && \
  git clone --recursive  https://github.com/hall-lab/svtyper.git && \
  chmod -R 755 ./svtyper && \
  cp -v ./svtyper/svtyper /usr/local/bin/ && \
  cp -vr ./svtyper/scripts/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r /usr/local/ngs/bin/svtyper && \

# wham: https://github.com/zeeev/wham
  cd  /usr/local/ngs/bin/ && \
  pip install -U numpy scipy scikit-learn argparse && \
  git clone --recursive https://github.com/jewmanchue/wham.git && \
  chmod -R 755 ./wham && \
  cd wham && \
  make && \
  cp -v ./wham/WHAM-BAM /usr/local/bin/ && \
  cp -v ./wham/mergeIndvs /usr/local/bin/ && \
  cp -v ./wham/WHAM-GRAPHENING/usr/local/bin/ && \

# manta
  MANTA_VERSON="0.29.4" && \
  cd  /usr/local/ngs/bin/ && \
  wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSON}/manta-${MANTA_VERSON}.centos5_x86_64.tar.bz2 && \
  tar -xjvf manta-${MANTA_VERSON}.centos5_x86_64.tar.bz2 && \


## ReAligners ------------------------------------------------------------------
# ABRA - Assembly Based ReAligner https://github.com/mozack/abra
RUN cd /usr/local/ngs/bin && \
  ABRA_VERSION="0.94" && \
  wget https://github.com/mozack/abra/releases/download/v${ABRA_VERSION}/abra-${ABRA_VERSION}-SNAPSHOT-jar-with-dependencies.jar && \
  chmod 755 abra-${ABRA_VERSION}-SNAPSHOT-jar-with-dependencies.jar && \
  ln -s abra-0.94-SNAPSHOT-jar-with-dependencies.jar abra.jar && \

# glia
  cd /usr/local/ngs/bin && \
  git clone --recursive https://github.com/ekg/glia.git && \
  chmod -R 755 ./glia && \
  cd ./glia && \
  make && \
  chmod -R 755 ./* && \
  cp -v ./glia /usr/local/bin/ && \
  cp -v fastahack/fastahack /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r glia && \

# ogap
  cd /usr/local/ngs/bin/ && \
  git clone --recursive https://github.com/ekg/ogap.git && \
  cd ogap && \
  make all && \
  chmod -R 755 ./* && \
  cp -v ogap /usr/local/bin/ && \
  cp -v ./smithwaterman/smithwaterman /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r ogap

## Nextflow --------------------------------------------------------------------
RUN cd /usr/local/bin/ && curl -fsSL get.nextflow.io | bash && cd

# source .bashrc
RUN  bash -c "source /home/ngseasy/.bashrc" && \

# Clean up APT when done.
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

ADD fix_ambiguous /usr/local/bin/

# PERMISSIONS
RUN chmod -R 755 /usr/local/ngs/bin && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin

# run as non root user ngseasy
ENV HOME /home/ngseasy
ENV PATH /usr/local/ngs/bin:/usr/local/bin:$PATH
USER ngseasy
WORKDIR /home/ngseasy
VOLUME /home/ngseasy/ngs_projects
VOLUME /home/ngseasy/scripts
CMD ["/bin/bash"]
