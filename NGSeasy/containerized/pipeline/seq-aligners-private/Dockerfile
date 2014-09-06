#########################################################################
# -- Author: Amos Folarin                                               #
# -- Author: Stephen J Newhouse                                         #
# -- Organisation: KCL/SLaM                                             #
# -- Email: amosfolarin@gmail.com                                       #
#########################################################################


#------------------------------------------------------------------------
#
#                            #######
#     #    #   ####    ####  #          ##     ####    #   #
#     ##   #  #    #  #      #         #  #   #         # #
#     # #  #  #        ####  #####    #    #   ####      #
#     #  # #  #  ###       # #        ######       #     #
#     #   ##  #    #  #    # #        #    #  #    #     #
#     #    #   ####    ####  #######  #    #   ####      #
#
# This dockerfile should build the enviroment required for the ngs 
# pipeline. Alternatively you can just get the docker container pre-built
# from our repository
#------------------------------------------------------------------------

#As a multi-component system, NGS pipeline setup is traditionally heavy on
#configuration. Our idea is to provide this in a simple encapsulated container.
#Users also typically wish to configure their own environments and run the
#pipeline on a wide range of hardware (workstations to clusters to cloud), being
#able to stand-up a pipeline with minimal fuss is made straightforward with this
#container.
#
# The containers for this NGS pipeline:
# 
# ngseasy pipeline 
#   https://github.com/KHP-Informatics/ngs/blob/master/containerized/pipeline/Dockerfile
# ngs_storage container 
#   https://github.com/KHP-Informatics/ngs/blob/master/containerized/storage/Dockerfile
#
#- Pipeline components are stored in /usr/local/pipeline/
#- The user for running the pipeline is: pipeman
#- Intermediary pipeine files are stored in the container under this dir
#   structure: TODO
#
# The system looks something like this:
#
#>     [ ngseasy  container ]... x n
#>     |
#>     |____[ storage container ]  #reference genomes
#>     |
#>     |____{ mounted volume pwd }  #pass config file from host
#>     |
#>     |____{ volume //data } #container output
#

#------------------------------------------------------------------------
# BUILDING THE DOCKER IMAGE FROM THIS Dockerfile:
# download this Dockerfile into a clean directory e.g. build_dir
# Due to licencing, you should get your own versions of:
# 
# 1) Download the compressed files into the build_dir, then edit the lines
# of the Dockerfile to the correct versions of ANNOVAR and GATK
#   * novoalign http://www.novocraft.com/
#   * Stampy http://www.well.ox.ac.uk/project-stampy
#   * GATK https://www.broadinstitute.org/gatk/
#   * ANNOVAR http://www.openbioinformatics.org/annovar/
#
# 2) $ cd build_dir
# 3) $ sudo docker build --tag <repo-name:tag> .
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# VOLUMES and VOLUME CONTAINERS
# Mount a volume external to the container to store stuff that needs to persist
#------------------------------------------------------------------------
#------ Pass in the config file via another volume, mount host pwd
# You should mount the pwd when you launch docker, and then make sure the config file 
# resides in pwd.
# Then the docker process can read the file from /tmp/config/<configfile>
#   $ docker run -v .:/tmp/config
#
#------ RUN mkdir /pipeln
# FASTQ Staging Area
#RUN bash mkdir /pipeln/staging/ ## !!! this should be in an external volume !!!!
# REF Data
#RUN bash mkdir /pipeln/refdata ## !!! this should be in an external volume!!!!!
# Results Folder by projetc and sample

#------ Reference sequence data TODO:
# Option1. Mount a volume in a dedicated Volume Container *Probably Preferred*
#   $ docker run volumes-from <storage_container> 
# the storage_container should have a directory with the data
# /media/ngsEasy/data/reference
#
# Or...
#
# Option2.  mount shared directory on host
#   $ docker run -v /*host*/data/reference:/media/ngsEasy/data/reference

# Base image will be Trusty
FROM afolarin/base-config

# Maintainer Amos Folarin
    MAINTAINER Amos Folarin amosfolarin@gmail.com

# Set correct environment variables.
    ENV HOME /root

# install basic dependencies
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


#---------------------------------------------------------------------------------------
#***************** INSTALL NGS ALIGNMENTY PIPELINE COMPONENTS **************************
#---------------------------------------------------------------------------------------

#---------------------- Quality Control ----------------------------------

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

#---------------------- Sequence Aligners --------------------------------
# Aligner options: 1)*novoalign*  2)bowtie, 3)stampy, 4)bwa

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

#---------------------- GATK --------------------------------

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
  
#---------------------- b37 gatk resource -----------------------------
#RUN wget -O /tmp/gatk_resources.tar.gz https://www.dropbox.com/s/ff1b5l2xbhisidx/gatk_resources.tar.gz \
#  && tar xvf /tmp/gatk_resources.tar.gz -C /usr/local/pipeline/ \
#  && gunzip /usr/local/pipeline/gatk_resources/*.gz \ 
#  && chown -R pipeman:ngsgroup /usr/local/pipeline/gatk_resources \
#  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/gatk_resources' /home/pipeman/.bashrc \
#  && echo "alias ngsGATKresources='/usr/local/pipeline/gatk_resources'" >>  /home/pipeman/.bashrc


#---------------------- SAM/BED/VCF Processing tools  --------------------------------    
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


#RUN wget -O /tmp/samtools-0.1.19.tar.bz2 http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2 \
#    && tar xjvf /tmp/samtools-0.1.19.tar.bz2 -C /usr/local/pipeline/ \
#    && cd /usr/local/pipeline/samtools-0.1.19 && make \
#    && chown -R pipeman:ngsgroup /usr/local/pipeline/samtools-0.1.19 \
#    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/samtools-0.1.19/bin' /home/pipeman/.bashrc \
#    && echo "alias ngsSamtools='/usr/local/pipeline/samtools-0.1.19/bin'" >>  /home/pipeman/.bashrc

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

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/*** \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# change to bash for the default entrypoint. You can override ENTRYPOINT and CMD now in docker run.
# ENTRYPOINT ["/bin/bash"]

# Launch the pipeline through a command this will be a default executaion of the pipeline 
# and can be overridden with # docker run
# CMD ["/usr/local/pipeline/molpath_ngs/ngs_master_scripts/xxxxx"]



