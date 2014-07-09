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
#******************* IMAGE: afolarin/var-callers ************************
#
# This dockerfile should build the container image for varient caller
# components
#------------------------------------------------------------------------

# Starting image
FROM afolarin/seq-aligners

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com



#------------------------------------------------------------------------
#*********** INSTALL PIPELINE VARIENT CALLER COMPONENTS *****************
#------------------------------------------------------------------------




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
##RUN cd /usr/local/pipeline \
RUN git clone https://github.com/pezmaster31/bamtools.git \
&& cd bamtools/ \
&& mkdir build \
&& cd build/ \
&& cmake .. \
&& make
#&& chown -R pipeman:ngsgroup /usr/local/pipeline/bamtools

# + install seqtk
#RUN cd /usr/local/pipeline \
RUN git clone https://github.com/lh3/seqtk.git \
&& cd seqtk/ \
&& make
##&& chown -R pipeman:ngsgroup /usr/local/pipeline/seqtk

# + set environment
ENV BOOST_ROOT /usr
ENV BAMTOOLS_ROOT /bamtools
ENV SEQTK_ROOT /seqtk
ENV LD_LIBRARY_PATH /bamtools/lib

# + install delly
#RUN cd /usr/local/pipeline \
RUN git clone https://github.com/tobiasrausch/delly.git \
&& cd delly/ \
&& make -B src/delly
#&& chown -R pipeman:ngsgroup /usr/local/pipeline/delly


#----------------------------PINDEL-------------------------------------#
# http://gmt.genome.wustl.edu/pindel/current/install.html
RUN cd /usr/local/pipeline \
&& git clone git://github.com/genome/pindel.git \
&& cd pindel \
&& ./INSTALL /usr/local/pipeline/samtools-0.1.19;

# TEST FOR UBUNTU BUILD ISSUES!

#----------------------------breakdancer-------------------------------------#
#RUN cd /usr/local/pipeline \
#&& git clone --recursive https://github.com/genome/breakdancer.git \
#&& cd breakdancer \
#&& mkdir build \
#&& cd build \
#&& cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local/pipeline \
#&& make \
#&& sudo make install
#Dependence over other Perl Modules
#------------------------------------
#use Statistics::Descriptive;
#use Math::CDF;#
#These are available at CPAN.
#http://search.cpan.org/~colink/Statistics-Descriptive-2.6/Descriptive.pm
#http://search.cpan.org/~callahan/Math-CDF-0.1/CDF.pm
#use Poisson;
#This is provided. Please make sure the "use lib" at the beginning of the perl scripts contains the correct path.

# http://sv.gersteinlab.org/breakseq/
# http://root.cern.ch/drupal/content/production-version-534


#---------------------------CNVnator--------------------------------#
#RUN wget -O /tmp/CNVnator_v0.3.zip http://sv.gersteinlab.org/cnvnator/CNVnator_v0.3.zip \
#&& unzip /tmp/CNVnator_v0.3.zip -d /usr/local/pipeline/ \
#&& chown -R pipeman:ngsgroup /usr/local/pipeline/CNVnator_v0.3

#------------------------m-HMM-----------------------------------------#
# https://www.stt.msu.edu/users/hengwang/mHMM.html
RUN cd /usr/local/pipeline/ \
&& wget https://www.stt.msu.edu/users/hengwang/mHMM%20Data/mHMM_1.0.tar.gz \
&& R CMD INSTALL mHMM_1.0.tar.gz \
&& rm mHMM_1.0.tar.gz;

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 755 /usr/local/pipeline/*


# Cleanup the temp dir
RUN rm -rf /tmp/*


# change to bash for the default entrypoint. You can override ENTRYPOINT and CMD now in docker run.
# ENTRYPOINT /bin/bash

# Launch the pipeline through a command this will be a default executaion of the pipeline 
# and can be overridden with # docker run
#CMD ["/usr/local/pipeline/molpath_ngs/ngs_master_scripts/xxxxx"]



