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

# Base image will be latest
FROM ubuntu

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com

# Set correct environment variables.
ENV HOME /root

RUN apt-get update


#------------------------------------------------------------------------
# USER SETUP
#------------------------------------------------------------------------
# Create a pipeline user:pipeman and group:ngsgroup
RUN useradd -m -s /bin/bash pipeman && cd /home/pipeman && echo "#bash config file for user pipeman" >> /home/pipeman/.bashrc
RUN groupadd ngsgroup
RUN usermod -G ngsgroup pipeman
# USER pipeman

# give root access
# echo "pipeman	ALL=(ALL)	ALL" >>  /etc/sudoers && echo "%ngsgroup	ALL=(ALL)	ALL" >>  /etc/sudoers




