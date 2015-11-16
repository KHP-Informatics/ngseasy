FROM debian:jessie
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com
LABEL Description="This is the base image for compbio. Based on debian:jessie. Fat Image" NickName="little-fatty-deb" URL="https://hub.docker.com/r/library/debian/" Version="1.0"
# Remain current and upgrade apt-get and add additional repos
RUN sed -i '$adeb http://cran.ma.imperial.ac.uk/bin/linux/debian jessie-cran3/' /etc/apt/sources.list && \
  apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 && \
  DEBIAN_FRONTEND=noninteractive \
  apt-get update && apt-get dist-upgrade -y && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y \
  ant \
  asciidoc \
  automake \
  bash \
  binutils \
  build-essential \
  bzip2 \
  cmake \
  curl \
  dkms \
  dpkg-dev \
  gcc \
  g++ \
  gpp \
  gcc \
  git \
  git-core \
  gfortran \
  gnuplot \
  gradle \
  graphviz \
  htop \
  libatlas-dev \
  libblas-dev \
  libbz2-dev \
  libc-bin \
  libconfig-dev \
  libcurl4-openssl-dev \
  libfreetype6-dev \
  liblapack-dev \
  liblzma-dev \
  libpcre3-dev \
  libpng-dev \
  libreadline-dev \
  libssl-dev \
  libxml2-dev \
  libxml2-utils \
  llvm \
  locales \
  make \
  ncurses-dev \
  openssl \
  openssl-blacklist \
  parallel \
  pkg-config \
  python \
  python-dev \
  python2.7-dev \
  python-pip \
  python-yaml \
  ssl-cert \
  sudo \
  tabix \
  time \
  tree \
  unzip \
  vim \
  wget \
  zlib1g \
# java
  openjdk-7-jdk \
  openjdk-7-doc \
  openjdk-7-jre-lib \
# R
  libatlas3-base \
  libopenblas-base \
  r-base \
  r-base-dev \
  littler && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
# python modules
  pip install \
  Cython \
  pytabix \
  python-dateutil \
  numpy \
  matplotlib \
  regex \
  deepdiff \
  intervaltree \
  pygr \
  XlsxWriter && \

# Ensure permissions are set for update in place by arbitrary users
# From: https://github.com/chapmanb/bcbio-nextgen/blob/master/Dockerfile#L68
  find /usr/local -perm /u+x -execdir chmod a+x {} \; && \
  find /usr/local -perm /u+w -execdir chmod a+w {} \; && \

# clean up
  apt-get clean && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /tmp/* && \
  rm -rf /var/tmp/* && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

# set JAVA_HOME
  ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java

# Use baseimage-docker's bash.
CMD ["/bin/bash"]
