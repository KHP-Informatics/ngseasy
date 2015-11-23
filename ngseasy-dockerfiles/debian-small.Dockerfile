FROM debian:jessie

MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

LABEL Description="This is the base image for compbio. Based on debian:jessie. This is a Slim Image" NickName="little-deb" URL="https://hub.docker.com/r/library/debian/" Version="1.0"

ENV DEBIAN_FRONTEND noninteractive

# Remain current, upgrade apt-get, add build tools and java
RUN apt-get update && \
  apt-get dist-upgrade -y && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
  apt-get install -y --no-install-recommends \
  apt-utils \
  automake \
  bash \
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
  locales \
  make \
  parallel \
  pkg-config \
  sudo \
  tabix \
  time \
  tree \
  unzip \
  vim \
  wget && \
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

# Ensure permissions are set for update in place by arbitrary users
# From: https://github.com/chapmanb/bcbio-nextgen/blob/master/Dockerfile#L68
  find /usr/local -perm /u+x -execdir chmod a+x {} \; && \
  find /usr/local -perm /u+w -execdir chmod a+w {} \; && \

# configure locales
dpkg-reconfigure locales && \
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

# clean up
RUN apt-get clean && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /tmp/* && \
  rm -rf /var/tmp/* && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

# Use baseimage-docker's bash.
CMD ["/bin/bash"]
