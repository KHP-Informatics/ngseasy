FROM debian:jessie

MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

LABEL Description="This is the debian:jessie base image for compbio ngseasy builds. This is Fat Dev Box Image to help get around the dependency mess of all the components" NickName="little-fatty-deb"

ENV DEBIAN_FRONTEND noninteractive

# set R version
# https://github.com/rocker-org/rocker/blob/master/r-base/Dockerfile#L38
ENV R_BASE_VERSION 3.2.3

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
  apt-utils && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
# ant \
# asciidoc \
  apt-get install -y --no-install-recommends \
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
  parallel \
  pkg-config \
  python \
  python-dev \
  python2.7-dev \
  python-pip \
  python-yaml \
  tabix \
  time \
  tree \
  unzip \
  vim \
  wget \
  zlib1g && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge

# java
# http://www.webupd8.org/2014/03/how-to-install-oracle-java-8-in-debian.html
RUN JAVA_INSTALL_VERSION="8" && \
  echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" | tee /etc/apt/sources.list.d/webupd8team-java.list && \
  echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" | tee -a /etc/apt/sources.list.d/webupd8team-java.list && \
  apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886 && \
  apt-get update && \
  echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections && \
  apt-get install -y --no-install-recommends \
  oracle-java8-installer && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
  apt-get install -y --no-install-recommends oracle-java8-set-default && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge

# Install R
# https://github.com/rocker-org/rocker/blob/master/r-base/Dockerfile#L45
RUN apt-get install -y --no-install-recommends \
  libatlas3-base \
  libopenblas-base && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
  apt-get install -y --no-install-recommends \
  r-base=${R_BASE_VERSION}* \
  r-base-dev=${R_BASE_VERSION}* \
  r-recommended=${R_BASE_VERSION}* && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge

# python modules scipy stack
RUN apt-get install -y  --no-install-recommends \
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
  usermod -aG sudo ngseasy && \

# make dirs: /usr/local/ngs/bin and sort permissions out
  mkdir /usr/local/ngs && \
  mkdir /usr/local/ngs/bin && \
  chown ngseasy:ngseasy /usr/local/ngs/bin  && \
  chmod -R 755 /usr/local/ngs/bin  && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin

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
ENV JAVA_HOME /usr/lib/jvm/java-8-oracle
# Use baseimage-docker's bash.
CMD ["/bin/bash"]
