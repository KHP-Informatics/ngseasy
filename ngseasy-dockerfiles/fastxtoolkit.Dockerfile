# NGSeasy Base Image

# FROM 
FROM debian:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# basic dependencies
RUN apt-get update \
  && apt-get install -y \
  automake \
  ant \
  bash \
  binutils \
  bioperl \
  build-essential \
  bzip2 \
  c++11 \
  cdbs \
  cmake \
  cron \
  curl \
  dkms \
  dpkg-dev \
  g++ \
  gcc \
  gfortran \
  git \
  git-core \
  libblas-dev libatlas-dev libbz2-dev liblzma-dev \
  libpcre3-dev \
  libreadline-dev \
  zlib1g-dev \  
  make \
  mercurial \
  openjdk-7-jdk \
  openjdk-7-jre \
  perl \
  php5-curl \
  python \
  subversion \
  tabix \
  tree \
  unzip \
  vim \
  wget \
  fastx-toolkit


#------------------------------------------------------------------------
# USER SETUP
#------------------------------------------------------------------------
# Create a pipeline user:ngseasy and group:ngseasy
RUN useradd -m -s /bin/bash ngseasy && cd /home/ngseasy && echo "#bash config file for user ngseasy" >> /home/ngseasy/.bashrc
RUN groupadd ngseasy
RUN usermod -G ngseasy ngseasy
RUN mkdir /usr/local/pipeline && chown ngseasy:ngseasy /usr/local/pipeline && chmod 775 /usr/local/pipeline

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



USER ngseasy
WORKDIR /home/ngseasy
