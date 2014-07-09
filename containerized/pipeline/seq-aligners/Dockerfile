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
#******************** IMAGE: afolarin/seq-aligners *******************************
#
# This dockerfile should build the sequence alignment conmponents and 
# quality control tools
#------------------------------------------------------------------------

# Starting image
FROM afolarin/base-config

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com

#------------------------------------------------------------------------
#******** INSTALL PIPELINE SEQUENCE ALIGNMENT COMPONENTS ****************
#------------------------------------------------------------------------

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
# + get novoalign.lic from novoalign ~ $1000 and put in context dir for the build
# + get updated version of novosort novosortV1.03.01.Linux3.0.tar.gz and put in context dir for the build
ADD novocraftV3.02.06.Linux3.0.tar.gz /usr/local/pipeline/novocraftV3.02.06.Linux3.0
ADD novoalign.lic /usr/local/pipeline/novocraftV3.02.06.Linux3.0/novocraft/
ADD novosortV1.03.01.Linux3.0.tar.gz  /usr/local/pipeline/novocraftV3.02.06.Linux3.0/novocraft/
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



#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 755 /usr/local/pipeline/*


# Cleanup the temp dir
RUN rm -rf /tmp/*


# change to bash for the default entrypoint. You can override ENTRYPOINT and CMD now in docker run.
# ENTRYPOINT /bin/bash

# Launch the pipeline through a command this will be a default executaion of the pipeline 
# and can be overridden with # docker run
#CMD ["/usr/local/pipeline/molpath_ngs/ngs_master_scripts/xxxxx"]



