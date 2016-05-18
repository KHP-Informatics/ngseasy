FROM debian:7.4

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# update and get some things
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get install -y \
    wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion \
    make cmake gcc build-essential zlib1g-dev \
    curl grep sed dpkg \
    unzip tree && \
    apt-get autoremove -y && \
    apt-get autoclean && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# http://bugs.python.org/issue19846
# > At the moment, setting "LANG=C" on a Linux system *fundamentally breaks Python 3*, and that's not OK.
ENV LANG C.UTF-8

# User set up
# Create a user:ngseasy and group:ngseasy
RUN useradd -m -U -s /bin/bash ngseasy && \
  cd /home/ngseasy && \
  chown ngseasy:ngseasy /home/ngseasy && \
  usermod -aG sudo ngseasy && \
  echo 'export PATH=/home/ngseasy/bin:$PATH' > /etc/profile.d/conda.sh

# switch to ngseasy user
USER ngseasy

# Anaconda2 install
RUN cd /home/ngseasy && \
    wget --quiet https://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh && \
    /bin/bash ./Anaconda2-4.0.0-Linux-x86_64.sh -b -p /home/ngseasy/conda && \
    rm ./Anaconda2-4.0.0-Linux-x86_64.sh

# add conda bin to path
ENV PATH /home/ngseasy/conda/bin:$PATH

# volumes
VOLUME /home/ngseasy
VOLUME /home/ngseasy/reference_genomes
VOLUME /home/ngseasy/ngs_projects
VOLUME /home/ngseasy/scratch

# set Home
ENV HOME /home/ngseasy
WORKDIR /home/ngseasy

# Expose for future
EXPOSE 80
EXPOSE 8080

# starting point
CMD [ "/bin/bash" ]
