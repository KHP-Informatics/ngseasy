#########################################################################
# -- Author: Amos Folarin                                               #
# -- Author: Stephen J Newhouse                                         #
# -- Organisation: KCL/SLaM                                             #
# -- Email: amosfolarin@gmail.com                                       #
#########################################################################


#------------------------------------------------------------------------
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

# Base image will be Trusty
FROM ubuntu:trusty

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com

# Required basic stuff, make, gcc, wget etc
RUN apt-get update
RUN apt-get install -y gcc g++
RUN apt-get install -y python-dev
RUN apt-get install -y make tabix
RUN apt-get install -y git-core wget unzip
RUN apt-get install -y zlib1g-dev libncurses5-dev
RUN apt-get install -y sysvbanner

# delly requirements 
RUN apt-get install -y build-essential cmake zlib1g-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev

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



#------------------------------------------------------------------------
# USER SETUP
#------------------------------------------------------------------------
# Create a pipeline user:pipeman and group:ngsgroup
RUN useradd -m -s /bin/bash pipeman && cd /home/pipeman && echo "#bash config file for user pipeman" >> /home/pipeman/.bashrc
RUN groupadd ngsgroup
RUN usermod -G ngsgroup pipeman


#------------------------------------------------------------------------
#***************** INSTALL PIPELINE COMPONENTS **************************
#------------------------------------------------------------------------

# make FASTQ_STAGGING dir
RUN mkdir /usr/local/pipeline/FASTQ_STAGGING
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/FASTQ_STAGGING
# make reference_data dir
RUN mkdir /usr/local/pipeline/reference_data
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/reference_data

#---------------------- Quality Control ----------------------------------

# + Trimmomatic
RUN wget -O /tmp/Trimmomatic-0.32.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip \
    && unzip /tmp/Trimmomatic-0.32.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/Trimmomatic-0.32 \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/Trimmomatic-0.32/trimmomatic-0.32.jar' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/Trimmomatic-0.32' /home/pipeman/.bashrc

# + FastQC
RUN wget -O /tmp/fastqc_v0.11.2.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip \
    && unzip /tmp/fastqc_v0.11.2.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/FastQC \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/FastQC/jbzip2-0.9.jar:/usr/local/pipeline/FastQC/sam-1.103.jar' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/FastQC' /home/pipeman/.bashrc
    

#---------------------- Sequence Aligners --------------------------------
# Aligner options: 1)*novoalign*  2)bowtie, 3)stampy, 4)bwa

# + novoalign (registration required,  get compressed file and put in context dir for the build)
ADD novocraftV3.02.06.Linux3.0.tar.gz /usr/local/pipeline/novocraftV3.02.06.Linux3.0
RUN sed  -i '$aPATH=${PATH}:/usr/local/pipeline/novocraftV3.02.06.Linux3.0/novocraft' /home/pipeman/.bashrc

# + bowtie
RUN wget -O /tmp/bowtie2-2.2.3-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/download \
    && unzip /tmp/bowtie2-2.2.3-linux-x86_64.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bowtie2-2.2.3 \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.3:/usr/local/pipeline/bowtie2-2.2.3/scripts' /home/pipeman/.bashrc 

# + stampy (registration required, get compressed file and put in context dir for the build)
ADD Stampy-latest.tgz /usr/local/pipeline/
# (below) require to replace the -Wl flags for building in ubuntu, see stampy README file
RUN sed -i 's/-Wl//' /usr/local/pipeline/stampy-1.0.23/makefile \
    && chmod -R 755 /usr/local/pipeline/stampy-1.0.23 \
    && cd /usr/local/pipeline/stampy-1.0.23 && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/stampy-1.0.23 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/stampy-1.0.23' /home/pipeman/.bashrc

# stampy install alternative if registration already completed -- seems to have a static download page.
#RUN wget Stampy-latest.tgz -O http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
#    && tar xzvf Stampy-latest.tgz -C /usr/local/pipeline/
#    && sed -i "$aPATH=${PATH}:"$(readlink -f stampy-*) /home/pipeman/.bashrc

#---------------------- SAM/BAM Processing -------------------------------
# + VCFtools: http://vcftools.sourceforge.net/index.html
RUN wget -O /tmp/vcftools_0.1.12a.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12a.tar.gz/download \
    && tar xzvf /tmp/vcftools_0.1.12a.tar.gz -C /usr/local/pipeline/  \
    && cd /usr/local/pipeline/vcftools_0.1.12a/ && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/vcftools_0.1.12a \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/vcftools_0.1.12a/bin' /home/pipeman/.bashrc



#---------------------- SNP/INDEL Calling --------------------------------

# + GATK (see licence, and registration) 
# use ADD to copy the downloaded GATK binary from <buildcontext_path>/GenomeAnalysisTK-3.1-1.tar.bz2 to the /tmp dir in the container
# NOTE: ADD automatically unpacks compressed files
ADD GenomeAnalysisTK-3.1-1.tar.bz2 /usr/local/pipeline/GenomeAnalysisTK-3.1-1 
RUN chmod -R 755 /usr/local/pipeline/GenomeAnalysisTK-3.1-1 \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/GenomeAnalysisTK-3.1-1 \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/GenomeAnalysisTK-3.1-1' /home/pipeman/.bashrc

# + BEDtools
RUN cd /usr/local/pipeline \
    && git clone https://github.com/arq5x/bedtools2.git \
    && cd bedtools2 && make clean && make all \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bedtools2 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bedtools2/bin' /home/pipeman/.bashrc

# + Picard
RUN wget -O /tmp/picard-tools-1.115.zip http://sourceforge.net/projects/picard/files/picard-tools/1.115/picard-tools-1.115.zip/download \
    && mkdir /usr/local/pipeline/picardtools \
    && unzip /tmp/picard-tools-1.115.zip -d /usr/local/pipeline/picardtools/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/picardtools \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/picardtools/picard-tools-1.115/snappy-java-1.0.3-rc3.jar' /home/pipeman/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/picardtools/picard-tools-1.115' /home/pipeman/.bashrc

# + SAMtools
RUN wget -O /tmp/samtools-0.1.19.tar.bz2 http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download \
    && tar xjvf /tmp/samtools-0.1.19.tar.bz2 -C /usr/local/pipeline/ \
    && cd /usr/local/pipeline/samtools-0.1.19 && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/samtools-0.1.19 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/samtools-0.1.19/bin' /home/pipeman/.bashrc


#---------------------------------- DELLY ---------------------------------------

# + install bamtools
RUN cd /usr/local/pipeline
RUN git clone https://github.com/pezmaster31/bamtools.git
RUN cd bamtools/ ; mkdir build ; cd build/ ; cmake .. ; make
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/bamtools

# + install seqtk
RUN cd /usr/local/pipeline
RUN git clone https://github.com/lh3/seqtk.git
RUN cd seqtk/ ; make
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/seqtk

# + set environment
ENV BOOST_ROOT /usr
ENV BAMTOOLS_ROOT /bamtools
ENV SEQTK_ROOT /seqtk
ENV LD_LIBRARY_PATH /bamtools/lib

# + install delly
RUN cd /usr/local/pipeline
RUN git clone https://github.com/tobiasrausch/delly.git
RUN cd delly/ ; make -B src/delly
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/delly



#---------------------- Visualisations -----------------------------------

# + R (replace with installing the dev, if there is need to compile packages)
#RUN sudo apt-get install r-base-dev
RUN sudo apt-get install -y r-base
#RUN R CMD INSTALL caret # doesn't do dependancies
#RUN /usr/bin/Rscript --no-save --no-restore <( cat <<< 'install.packages("caret", repos="http://cran.r-project.org")' )
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("caret",dependencies=TRUE ,repos="http://cran.r-project.org")'
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("gdata",dependencies=TRUE ,repos="http://cran.r-project.org")'
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("gplots",dependencies=TRUE ,repos="http://cran.r-project.org")'
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("ggplot2",dependencies=TRUE ,repos="http://cran.r-project.org")'
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("RColorBrewer",dependencies=TRUE ,repos="http://cran.r-project.org")'
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("doMC",dependencies=TRUE ,repos="http://cran.r-project.org")'
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("caTools",dependencies=TRUE ,repos="http://cran.r-project.org")'
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("utils",dependencies=TRUE ,repos="http://cran.r-project.org")'
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("lattice",dependencies=TRUE ,repos="http://cran.r-project.org")'


#------------------------------------bioconductor---------------------------------
# install bioconductor -- mod the biocLite.R script so it fires default bioc installation
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite()'
# GenomicRanges
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("GenomicRanges",dependencies=TRUE)'
# IRanges
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("IRanges",dependencies=TRUE)'
# Biostrings
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("Biostrings",dependencies=TRUE)'
# genomeIntervals
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("genomeIntervals",dependencies=TRUE)'
# rtracklayer
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("rtracklayer",dependencies=TRUE)'
# workflows.R: variants
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://bioconductor.org/workflows.R");workflowInstall("variants")'
# BSgenome Rsamtools ShortRead
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome",dependencies=TRUE)'
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("Rsamtools",dependencies=TRUE)'
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("ShortRead",dependencies=TRUE)'


#------------------ Variant Annotation & Reports -------------------------
# + ANNOVAR (see licence, and registration)
# Available on reg:  http://www.openbioinformatics.org/annovar/download/Ht8qRwQSTi/annovar.latest.tar.gz
# hmmm... no version on tar.
ADD annovar.latest.tar.gz /usr/local/pipeline/
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/annovar \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/annovar' /home/pipeman/.bashrc

# + SnpEff
RUN wget -O /tmp/snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip \
    && unzip /tmp/snpEff_latest_core.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/snpEff \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/snpEff/snpEff.jar' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/snpEff' /home/pipeman/.bashrc


#----------------------------- vcf phaseing -----------------------------------

# used to phase vcf files for family and population based studies 
# + BEAGLE
RUN wget -O /usr/local/pipeline/beagle/beagle.jar http://faculty.washington.edu/browning/beagle/beagle.jar \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/beagle/ \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/beagle/beagle.jar' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/beagle' /home/pipeman/.bashrc



#------------------------- NGS PIPELINE ----------------------------------
RUN cd /usr/local/pipeline \
    && git clone https://github.com/KHP-Informatics/ngs.git \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/ngs/molpath_ngs:/usr/local/pipeline/ngs/molpath_ngs/ngs_master_scripts' /home/pipeman/.bashrc \

#----------------------  --------------------------------



#----------------------  --------------------------------




# Cleanup the temp dir
RUN rm -rf /tmp/*


# change to bash for the default entrypoint. You can override ENTRYPOINT and CMD now in docker run.
# ENTRYPOINT /bin/bash

# Launch the pipeline through a command this will be a default executaion of the pipeline 
# and can be overridden with # docker run
#CMD ["/usr/local/pipeline/molpath_ngs/ngs_master_scripts/xxxxx"]



