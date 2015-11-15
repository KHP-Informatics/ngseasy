# Base image
FROM compbio/ngseasy-base:1.0-r001
# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

RUN apt-get install -y gnuplot libnuma-dev libz-dev markdown zlib1g-dev doxygen

RUN cd /usr/local/ngs/bin && \
    git clone --recursive https://github.com/sequencing/isaac_aligner.git && \
    mkdir ISSAC && \
    /usr/local/ngs/bin/isaac_aligner/src/configure --prefix=/usr/local/ngs/bin/ISSAC && \
    make && \
    make install

RUN cd /usr/local/ngs/bin && \
    wget https://github.com/sequencing/isaac_variant_caller/archive/v1.0.7.tar.gz && \
    tar -xvf v1.0.7.tar.gz && \

#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
USER ngseasy
WORKDIR /home/ngseasy
