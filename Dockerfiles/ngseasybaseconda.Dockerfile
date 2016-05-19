FROM debian:stable

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# update and get some things
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    git \
    mercurial \
    subversion \
    make \
    cmake \
    gcc \
    gpp \
    build-essential \
    zlib1g-dev \
    patch \
    perl \
    libextutils-pkgconfig-perl \
    byacc \
    wget \
    curl \
    grep \
    sed \
    dpkg \
    debconf \
    unzip \
    tree \
    libc6 \
    llvm \
    locales \
    time && \
    apt-get autoremove -y && \
    apt-get autoclean && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# http://bugs.python.org/issue19846
# > At the moment, setting "LANG=C" on a Linux system *fundamentally breaks Python 3*, and that's not OK.
# configure locales
RUN dpkg-reconfigure locales && \
  locale-gen C.UTF-8 && \
  /usr/sbin/update-locale LANG=C.UTF-8 && \
  echo 'en_US.UTF-8 UTF-8' >> /etc/locale.gen && \
  locale-gen

# Set default locale for the environment
ENV LC_ALL C.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# User set up
# Create a user:ngseasy and group:ngseasy
RUN useradd -m -U -s /bin/bash ngseasy && \
  cd /home/ngseasy && \
  chown ngseasy:ngseasy /home/ngseasy && \
  usermod -aG sudo ngseasy && \
  echo 'export PATH=/home/ngseasy/bin:$PATH' > /etc/profile.d/conda.sh

RUN chown -R ngseasy:ngseasy /home/ngseasy/

# switch to ngseasy user
USER ngseasy

# Anaconda2 install
RUN cd /home/ngseasy && \
 wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh  && \
 /bin/bash /home/ngseasy/Miniconda3-latest-Linux-x86_64.sh -b -p /home/ngseasy/conda && \
 rm -v /home/ngseasy/Miniconda3-latest-Linux-x86_64.sh

RUN mkdir -p /home/ngseasy/conda/conda-bld/linux-64 /home/ngseasy/conda/conda-bld/osx-64

# add conda bin to path
ENV PATH /home/ngseasy/conda/bin:$PATH

# setup conda
ADD requirements.txt requirements.txt
RUN conda update -y conda
RUN conda install -y --file requirements.txt
RUN conda update -y conda-build
RUN conda index /home/ngseasy/conda/conda-bld/linux-64 /home/ngseasy/conda/conda-bld/osx-64

## add channels
RUN conda config --add channels bioconda
RUN conda config --add channels r
RUN conda config --add channels sjnewhouse

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
