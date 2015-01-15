# Base image
FROM compbio/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update

# NOTE: ADD automatically unpacks compressed files

#--------------------------GEM------------------------------------

RUN wget -O /tmp/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2 http://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2 \
  && tar xjvf /tmp/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2 -C /usr/local/pipeline/ \
  && chmod -R 755 /usr/local/pipeline/GEM-binaries-Linux-x86_64-core_i3-20130406-045632 \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/GEM-binaries-Linux-x86_64-core_i3-20130406-045632 \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin' /home/pipeman/.bashrc \
  && echo "alias ngsGEM='/usr/local/pipeline/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin'" >>  /home/pipeman/.bashrc

RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline
 
# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
