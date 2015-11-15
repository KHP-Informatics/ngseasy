FROM debian:jessie
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com
LABEL Description="This is the base debian:jessie image for compbio images" URL="https://hub.docker.com/r/library/debian/" Version="1.0"
# Remain current
RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y \
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
  libconfig-dev \
  zlib1g \
  openjdk-7-jdk \
  openjdk-7-doc \
  openjdk-7-jre-lib && \

# clean up
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /tmp/* && \
  rm -rf /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

# set JAVA_HOME
  ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java

# Use baseimage-docker's bash.
CMD ["/bin/bash"]
