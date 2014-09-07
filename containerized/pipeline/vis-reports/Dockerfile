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
#******************* IMAGE: afolarin/vis-reports ************************
#
# This dockerfile should build the container image for visualization and 
# report generation
#------------------------------------------------------------------------
# Starting image
FROM afolarin/var-anno

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com


#------------------------------------------------------------------------
#*** INSTALL PIPELINE VISUALIZATION AND REPORT GENERATION COMPONENTS ****
#------------------------------------------------------------------------

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
###RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://bioconductor.org/workflows.R");workflowInstall("variants")'
# BSgenome Rsamtools ShortRead
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("BSgenome",dependencies=TRUE)'
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("Rsamtools",dependencies=TRUE)'
RUN /usr/bin/Rscript --no-save --no-restore -e 'source("http://www.bioconductor.org/biocLite.R"); biocLite("ShortRead",dependencies=TRUE)'


#------------------ Variant Annotation & Reports -------------------------
# + ANNOVAR (see licence, and registration)
# Available on reg:  http://www.openbioinformatics.org/annovar/download/Ht8qRwQSTi/annovar.latest.tar.gz
# hmmm... no version on tar.
ADD annovar.latest.tar.gz /usr/local/pipeline/
##RUN tar xvf /usr/local/pipeline/annovar.latest.tar.gz
RUN cd /usr/local/pipeline/annovar \
&& cp -v /usr/local/pipeline/annovar/table_annovar.pl table_annovar.pl.original \
&& wget http://www.openbioinformatics.org/annovar/download/table_annovar.pl \
&& chown -R pipeman:ngsgroup /usr/local/pipeline/annovar \
&& sed -i '$aPATH=${PATH}:/usr/local/pipeline/annovar' /home/pipeman/.bashrc

# TO Do : ADD scripts to then dowload ANNOVAR DATABASES
# Alternative is mount annovar database from local machine -- prefer this option


# + SnpEff
RUN wget -O /tmp/snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip \
    && unzip /tmp/snpEff_latest_core.zip -d /usr/local/pipeline/ \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/snpEff \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/snpEff/snpEff.jar' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/snpEff' /home/pipeman/.bashrc


#----------------------------- vcf phaseing -----------------------------------

# used to phase vcf files for family and population based studies 
# + BEAGLE
RUN mkdir /usr/local/pipeline/beagle
RUN wget -O /usr/local/pipeline/beagle/beagle.jar http://faculty.washington.edu/browning/beagle/beagle.jar \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/beagle/ \
    && sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/pipeline/beagle/beagle.jar' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/beagle' /home/pipeman/.bashrc


#---------------------- b37 gatk resource -----------------------------
# HASH this out if you want to include gatk resource data files
ADD gatk_resources.tar.gz /usr/local/pipeline/gatk_resources \ 
&& chown -R pipeman:ngsgroup /usr/local/pipeline/gatk_resources


#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 755 /usr/local/pipeline/*



# Cleanup the temp dir
RUN rm -rf /tmp/*


# change to bash for the default entrypoint. You can override ENTRYPOINT and CMD now in docker run.
# ENTRYPOINT /bin/bash

# Launch the pipeline through a command this will be a default executaion of the pipeline 
# and can be overridden with # docker run
#CMD ["/usr/local/pipeline/molpath_ngs/ngs_master_scripts/xxxxx"]



