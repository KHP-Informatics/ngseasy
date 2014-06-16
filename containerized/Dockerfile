#########################################################################
# -- Author: Amos Folarin                                               #
# -- Organisation: KCL/SLaM                                             #
# -- Email: amosfolarin@gmail.com                                       #
#########################################################################


#------------------------------------------------------------------------
# This dockerfile should build the enviroment required for the ngs 
# pipeline. Of course you can just get the docker container pre-built
# from our repository
#------------------------------------------------------------------------


# Base image will be Trusty
FROM ubuntu:trusty

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com

# Required basic stuff, make, gcc, wget etc
RUN apt-get update
RUN apt-get install -y gcc g++
RUN apt-get install -y python
RUN apt-get install -y make tabix
RUN apt-get install -y git-core wget unzip
RUN apt-get install -y zlib1g-dev libncurses5-dev

# Mount a volume external to the container to store stuff that needs to persist
VOLUME ["/media/ngsEasy/data/reference"] # reference indexed genomes

# Set up file structure
# RUN mkdir /pipeln
# FASTQ Staging Area
#RUN bash mkdir /pipeln/staging/ ## !!! this should be in an external volume !!!!
# REF Data
#RUN bash mkdir /pipeln/refdata ## !!! this should be in an external volume!!!!!
# Results Folder by projetc and sample


# Create a pipeline user
RUN useradd -m -s /bin/bash pipeman && cd /home/pipeman && echo "#bash config file for user pipeman" >> /home/pipeman/.bashrc
RUN groupadd ngsgroup

#------------------ INSTALL PIPELINE COMPONENTS --------------------------

#---------------------- Sequence Aligners --------------------------------
# 1)*novoalign*  2)bowtie, 3)bwa, 4) stampy

# + bowtie
RUN wget -O /tmp/bowtie2-2.2.3-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/download \
    && unzip /tmp/bowtie2-2.2.3-linux-x86_64.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bowtie2-2.2.3 \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bowtie2-2.2.3' /home/pipeman/.bashrc 


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
ADD GenomeAnalysisTK-3.1-1.tar.bz2 /usr/local/pipeline/GenomeAnalysisTK-3.1-1.tar.bz2 
RUN chmod 755 /usr/local/pipeline/GenomeAnalysisTK-3.1-1.tar.bz2 \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/GenomeAnalysisTK-3.1-1.tar.bz2/GenomeAnalysisTK.jar' /home/pipeman/.bashrc

# + BEDtools
RUN cd /usr/local/pipeline \
    && git clone https://github.com/arq5x/bedtools2.git \
    && cd bedtools2 && make clean && make all \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/bedtools2 \
    && sed -i '$aPATH=${PATH}:/usr/local/pipeline/bedtools2/bin' /home/pipeman/.bashrc

# + Picard
RUN wget -O /tmp/picard-tools-1.115.zip http://sourceforge.net/projects/picard/files/picard-tools/1.115/picard-tools-1.115.zip/download \
    && mkdir /usr/local/pipeline/picardtools
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



#------------------ Variant Annotation & Reports -------------------------
# + ANNOVAR (see licence, and registration)





#----------------------  --------------------------------



#----------------------  --------------------------------



#----------------------  --------------------------------











# Cleanup the temp dir
RUN rm -rf /tmp/*



# Add an entrypoint or CMD that runs the pipeline, taking in the config file to parameterize
#ENTRYPOINT ["/app/init"]
#CMD ["app:start"]
