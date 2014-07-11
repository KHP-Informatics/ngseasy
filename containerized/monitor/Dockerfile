########################################################################
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
# This dockerfile should construct a monitor for the ngsEasy pipeline
# this includes logging and cAdvisor among other things.
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# BUILDING THE DOCKER IMAGE FROM THIS Dockerfile:
# 2) $ cd build_dir
# 3) $ sudo docker build --tag <repo-name:tag> .
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# Example usage
# sudo docker run \
#  --volume=/var/run:/var/run:rw \
#  --volume=/sys/fs/cgroup/:/sys/fs/cgroup:ro \
#  --volume=/var/lib/docker/:/var/lib/docker:ro \
#  --publish=8080:8080 \
#  --detach=true \
#------------------------------------------------------------------------

# Base image will be Trusty
FROM ubuntu:trusty

# Maintainer Amos Folarin
MAINTAINER Amos Folarin amosfolarin@gmail.com

# Required basic stuff, make, gcc, wget etc
RUN apt-get update


#------------------------------------------------------------------------------
# Install Supervisord and add config file to manage all monitors etc.
#------------------------------------------------------------------------------
RUN apt-get install -y openssh-server apache2 supervisor \
    && mkdir -p /var/run/sshd \
    && mkdir -p /var/log/supervisor
#supervisord.conf req in build context
ADD supervisord.conf /etc/supervisor/conf.d/supervisord.conf

#------------------------------------------------------------------------------
# PORTs exposed on the monitor container
#    -> 8080 for cAdvisor
#    -> 22, 80 for supervisord
#------------------------------------------------------------------------------
EXPOSE 8080 22 80

#------------------------------------------------------------------------------
# Run supervisord, you can use this to manage any microservices required to 
# handle things like logging and monitoring etc..
#------------------------------------------------------------------------------

#default command for entrypoint
CMD ["/usr/bin/supervisord"]

#default entrypoint report value of CMD to terminal, can override with --entrypoint
#ENTRYPOINT



