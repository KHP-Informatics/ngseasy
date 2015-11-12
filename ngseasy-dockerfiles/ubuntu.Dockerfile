# base image
FROM ubuntu:14.04.3
# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com
LABEL Description="This is the base ubuntu:14.04.3 image for compbio" Version="1.0"
# Remain current
RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get upgrade -y
# Basic dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y \
  apt-utils \
  automake \
  bash \
  binutils \
  build-essential \
  bzip2 \
  c++11 \
  cmake \
  cron \
  curl \
  dkms \
  dpkg-dev \
  g++ \
  gpp \
  gcc \
  git \
  git-core \
  libreadline-dev \
  make \
  php5-curl \
  sudo \
  tabix \
  tree \
  unzip \
  vim \
  wget \
  libc-bin \
  llvm \
  libconfig-dev \
  ncurses-dev \
  zlib1g-dev \
  libX11-dev libXpm-dev libXft-dev libXext-dev \
  openjdk-7-jdk \
  openjdk-7-doc \
  openjdk-7-jre-lib

# set JAVA_HOME
ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java

# Clean up APT when done.
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# Use baseimage-docker's bash.
  CMD ["/bin/bash"]
