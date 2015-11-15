# Base image
FROM compbio/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y

# + seqtk
    RUN cd /usr/local/ngs/bin/ \
    && git clone https://github.com/lh3/seqtk.git \
    && chown -R ngseasy:ngseasy /usr/local/ngs/bin/seqtk \
    && cd seqtk/ \
    && make \
    && sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/seqtk' /home/ngseasy/.bashrc \
    && echo "alias ngsSeqtk='/usr/local/ngs/bin/seqtk'" >>  /home/ngseasy/.bashrc \


#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/ngs/bin/***
RUN chown -R ngseasy:ngseasy /usr/local/ngs/bin

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
USER ngseasy
WORKDIR /home/ngseasy
