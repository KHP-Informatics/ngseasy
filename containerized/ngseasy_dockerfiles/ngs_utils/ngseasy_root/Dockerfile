################################################################################
# Base image
################################################################################
FROM compbio/ngseasy-base:1.0-r001

################################################################################
# Maintainer
################################################################################
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

################################################################################
# Set correct environment variables.
################################################################################
##ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

################################################################################
# general purpose tools
################################################################################
RUN apt-get update -y && apt-get upgrade -y

# python and stuff
RUN apt-get install -y \
                    make automake cmake gcc gcc c++ ncurses-dev \
                    gfortran libghc-pipes-zlib-dev libghc-zlib-dev libghc-bzlib-dev \
                    libblas-dev liblapack3 liblapack-dev libatlas-base-dev \
                    libblas3 libblas-dev \
                    python python-dev python-yaml python-pip \
                    python-numpy python-scipy python-pysam

# cython
RUN cd /tmp && \
    wget http://cython.org/release/Cython-0.22.tar.gz && \
    chmod 777 Cython-0.22.tar.gz && \
    tar xvf Cython-0.22.tar.gz && \
    chmod -R 777 Cython-0.22 && \
    cd Cython-0.22 && \
    python setup.py install

################################################################################
# install dependencies
################################################################################
RUN apt-get install -y libX11-dev libXpm-dev libXft-dev libXext-dev

################################################################################
# get ROOT
################################################################################
RUN cd /tmp && \
    curl -OL ftp://root.cern.ch/root/root_v5.34.20.source.tar.gz && \
    tar -xvf root_v5.34.20.source.tar.gz && \
    cd root && \
    ./configure && \
    make -j `ncpu` && \
    cd .. && \
    mv root /usr/local

RUN  echo "source /usr/local/root/bin/thisroot.sh" >>  /root/.bashrc  && \
    /bin/bash -c ". /root/.bashrc && source /root/.bashrc"

################################################################################
# speedseq
################################################################################
ENV ROOTSYS /usr/local/root

RUN cd /usr/local/pipeline && \
    git clone --recursive https://github.com/cc2qe/speedseq && \
    chmod -R 777 /usr/local/pipeline/speedseq && \
    cd /usr/local/pipeline/speedseq && \
    make -j `ncpu` ROOTFLAGS=" -pthread -m64" && \
    chmod -R 777 /usr/local/pipeline/speedseq

RUN echo -e "export PATH=/usr/local/pipeline/speedseq/bin:\$PATH" >> ~/.bashrc && \
    echo -e "export PATH=/usr/local/pipeline/speedseq/bin:\$PATH" >> /home/pipman/.bashrc && \
    /bin/bash -c "source ~/.bashrc"

################################################################################
# PERMISSIONS
################################################################################
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

################################################################################
#Cleanup the temp dir
################################################################################
RUN rm -rf /tmp/*

################################################################################
#open ports private only
################################################################################
EXPOSE 8080

################################################################################
# Use baseimage-docker's bash.
################################################################################
CMD ["/bin/bash"]

################################################################################
#Clean up APT when done.
################################################################################
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
