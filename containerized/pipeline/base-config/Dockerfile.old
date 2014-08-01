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
#******************* IMAGE: afolarin/base-config ************************
# 
# This dockerfile should build the base image with configuration setup
# required for the pipeline container
#------------------------------------------------------------------------

# Starting image 
FROM ubuntu:trusty

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com

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
# FILE SYS SETUP
#------------------------------------------------------------------------
# make FASTQ_STAGGING dir
RUN mkdir -p /usr/local/pipeline/FASTQ_STAGGING
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/FASTQ_STAGGING
# make reference_data dir
RUN mkdir -p /usr/local/pipeline/reference_data
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/reference_data


#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 755 /usr/local/pipeline/*



# Cleanup the temp dir
RUN rm -rf /tmp/*



