#########################################################################
# -- Author: Amos Folarin                                               #
# -- Author: Stephen J Newhouse                                         #
# -- Organisation: KCL/SLaM                                             #
# -- Email: amosfolarin@gmail.com                                       #
#########################################################################

#------------------------------------------------------------------------
# CONTAINER INVENTORY: Ruffus                                           |
#   A pipeline manager for NGSEasy                                      |
#                                                                       |
#------------------------------------------------------------------------
#      ####   #    #  ######  ######  #    #   ####
#     #    #  #    #  #       #       #    #  #
#     #    #  #    #  #####   #####   #    #   ####
#     #####   #    #  #       #       #    #       #
#     #   #   #    #  #       #       #    #  #    #
#     #    #   ####   #       #        ####    ####

#------------------------------------------------------------------------
# BUILDING THE DOCKER IMAGE FROM THIS Dockerfile:
# download this Dockerfile into a clean directory e.g. build_dir
# 1) Download the Dockerfile
# 2) $ cd build_dir
# 3) $ sudo docker build --tag <repo-name:tag> .
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# VOLUMES and VOLUME CONTAINERS
# Mount a volume external to the container to store stuff that needs to persist
#------------------------------------------------------------------------

# Base image will be Trusty
FROM compbio/ngseasy-base:wheezy

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com


# install basic dependencies
RUN apt-get update && \
    apt-get install -y gcc \
    g++ \
    python-dev

        
#-------------------Install Ruffus Container-------------------------
RUN cd /usr/local/pipeline/ && \
    git clone https://bunbun68@code.google.com/p/ruffus/  && \
    cd /usr/local/pipeline/ruffus && \
    python ./setup.py install
    

# Cleanup the temp dir
RUN rm -rf /tmp/*

# set container user
#USER doesn't seem to set HOME dir? bug? manually set it...
ENV HOME /home/pipeman  
USER pipeman
WORKDIR /home/pipeman
