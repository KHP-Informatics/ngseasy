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


#------------------------------------------------------------------------
# VOLUMES and VOLUME CONTAINERS
# Mount a volume external to the container to store stuff that needs to persist
#------------------------------------------------------------------------
# VOLUME ["/media/ngsEasy/data/reference"]
# VOLUME ["/media/ngsEasy/data/reference"]
# Shared data should sit in a Volume Container, 

# RUN mkdir /pipeln
# FASTQ Staging Area
#RUN bash mkdir /pipeln/staging/ ## !!! this should be in an external volume !!!!
# REF Data
#RUN bash mkdir /pipeln/refdata ## !!! this should be in an external volume!!!!!
# Results Folder by projetc and sample



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

#---------------------- Sequence Aligners --------------------------------
# Aligner options: 1)*novoalign*  2)bowtie, 3)stampy, 4)bwa

# + novoalign
ADD novocraftV3.02.06.Linux3.0.tar.gz /usr/local/pipeline/novocraftV3.02.06.Linux3.0
RUN sed  -i '$aPATH=${PATH}:/usr/local/pipeline/novocraftV3.02.06.Linux3.0/novocraft' /home/pipeman/.bashrc

# + bowtie
RUN wget -O /tmp/bowtie2-2.2.3-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/download \
    && unzip /tmp/bowtie2-2.2.3-linux-x86_64.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bowtie2-2.2.3 \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.3' /home/pipeman/.bashrc 

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
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar' /home/pipeman/.bashrc

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
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/picardtools/picard-tools-1.115' /home/pipeman/.bashrc

# + SAMtools
RUN wget -O /tmp/samtools-0.1.19.tar.bz2 http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download \
    && tar xjvf /tmp/samtools-0.1.19.tar.bz2 -C /usr/local/pipeline/ \
    && cd /usr/local/pipeline/samtools-0.1.19 && make \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/samtools-0.1.19 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/samtools-0.1.19/bin' /home/pipeman/.bashrc



#---------------------- Visualisations -----------------------------------

# + R (replace with installing the dev, if there is need to compile packages)
#RUN sudo apt-get install r-base-dev
RUN sudo apt-get install -y r-base
#RUN R CMD INSTALL caret # doesn't do dependancies
#RUN /usr/bin/Rscript --no-save --no-restore <( cat <<< 'install.packages("caret", repos="http://cran.r-project.org")' )
RUN /usr/bin/Rscript --no-save --no-restore -e 'install.packages("caret", repos="http://cran.r-project.org")'
# install bioconductor -- mod the biocLite.R script so it fires default bioc installation
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite()'


#------------------ Variant Annotation & Reports -------------------------
# + ANNOVAR (see licence, and registration)
#ADD
#RUN



#---------------------- NGS PIPELINE  --------------------------------
RUN cd /usr/local/pipeline \
    && git clone https://github.com/KHP-Informatics/ngs.git \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/ngs/molpath_ngs' /home/pipeman/.bashrc
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/ngs/molpath_ngs/ngs_master_scripts' /home/pipeman/.bashrc

#----------------------  --------------------------------



#----------------------  --------------------------------







# Cleanup the temp dir
RUN rm -rf /tmp/*



# Add an entrypoint or CMD that runs the pipeline, taking in the config file to parameterize
#ENTRYPOINT ["/app/init"]
#CMD ["app:start"]
