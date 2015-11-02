# Base image
FROM compbio/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y

# + seqtk
    RUN cd /usr/local/pipeline/ \
    && git clone https://github.com/lh3/seqtk.git \
    && chown -R pipeman:ngsgroup /usr/local/pipeline/seqtk \
    && cd seqtk/ \
    && make \
    && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/seqtk' /home/pipeman/.bashrc \
    && echo "alias ngsSeqtk='/usr/local/pipeline/seqtk'" >>  /home/pipeman/.bashrc \


#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
