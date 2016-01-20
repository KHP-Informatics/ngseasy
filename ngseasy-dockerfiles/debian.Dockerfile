FROM debian:jessie

MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

LABEL Description="This is the debian:jessie base image for compbio ngseasy builds. This is Fat Dev Box Image to help get around the dependency mess of all the components" NickName="little-fatty-deb"

ENV DEBIAN_FRONTEND noninteractive

# set R version
# https://github.com/rocker-org/rocker/blob/master/r-base/Dockerfile#L38
ENV R_BASE_VERSION 3.2.2

# Remain current, upgrade apt-get, add build tools, R, JAVA and python
RUN sed -i '$adeb http://cran.ma.imperial.ac.uk/bin/linux/debian jessie-cran3/' /etc/apt/sources.list && \
  apt-key adv --keyserver keys.gnupg.net --recv-key 381BA480 && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get dist-upgrade -y && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
  apt-get install -y --no-install-recommends \
  apt-utils \
  ant \
  asciidoc \
  automake \
  bash \
  binutils \
  build-essential \
  bzip2 \
  ca-certificates \
  cmake \
  curl \
  dkms \
  dpkg-dev \
  debconf \
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
  less \
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
  zlib1g && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
# java
  apt-get install -y --no-install-recommends \
  openjdk-7-jdk \
  openjdk-7-doc \
  openjdk-7-jre-lib && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
# Install R
# https://github.com/rocker-org/rocker/blob/master/r-base/Dockerfile#L45
  apt-get install -y --no-install-recommends \
  libatlas3-base \
  libopenblas-base \
  r-base=${R_BASE_VERSION}* \
  r-base-dev=${R_BASE_VERSION}* \
  r-recommended=${R_BASE_VERSION}* \
  littler && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
# python modules scipy stack
  apt-get install -y  --no-install-recommends \
  python-biopython \
  python-numpy \
  python-scipy \
  python-matplotlib \
  python-reportlab \
  python-pandas \
  python-sympy \
  python-tk \
  python-nose && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
# install cython
  CYTHON_VERSION="0.23" && \
  cd /tmp && \
  wget http://cython.org/release/Cython-${CYTHON_VERSION}.tar.gz && \
  tar xvf Cython-${CYTHON_VERSION}.tar.gz && \
  chmod -R 755 Cython-${CYTHON_VERSION} && \
  cd Cython-${CYTHON_VERSION} && \
  python setup.py install && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \

# Ensure permissions are set for update in place by arbitrary users
# From: https://github.com/chapmanb/bcbio-nextgen/blob/master/Dockerfile#L68
  find /usr/local -perm /u+x -execdir chmod a+x {} \; && \
  find /usr/local -perm /u+w -execdir chmod a+w {} \; && \

# clean up
  apt-get update && \
  apt-get update -y && \
  apt-get clean && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /tmp/* && \
  rm -rf /var/tmp/* && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/ && \

# Create a user:ngseasy and group:ngseasy
  useradd -m -U -s /bin/bash ngseasy && \
  cd /home/ngseasy && \
  usermod -aG sudo ngseasy

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

# set JAVA_HOME
ENV JAVA_HOME /usr/lib/jvm/java-1.7.0-openjdk-amd64

# Use baseimage-docker's bash.
CMD ["/bin/bash"]
