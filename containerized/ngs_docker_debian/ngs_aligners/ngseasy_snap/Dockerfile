# Base image
FROM compbio/ngseasy-base:1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# + snap
RUN wget -O /tmp/snap-1.0beta.16-linux.tar.gz http://snap.cs.berkeley.edu/downloads/snap-1.0beta.16-linux.tar.gz \
    && tar xvf /tmp/snap-1.0beta.16-linux.tar.gz -C /usr/local/pipeline/ \
    && cd /usr/local/pipeline/snap-1.0beta.16-linux \
    && chmod -R 777 /usr/local/pipeline/snap-1.0beta.16-linux \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/snap-1.0beta.16-linux' /home/pipeman/.bashrc \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/snap-1.0beta.16-linux' ~/.bashrc \
    && cp -v /usr/local/pipeline/snap-1.0beta.16-linux/* /usr/local/bin/

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/


