FROM debian:jessie
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com
LABEL Description="This is the base debian:jessie image for compbio images" URL="https://hub.docker.com/r/library/debian/" Version="1.0"
# Set correct environment variables.
ENV DEBIAN_FRONTEND noninteractive
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
  libc-bin \
  llvm \
  libconfig-dev \
  ncurses-dev \
  zlib1g-dev \
  libX11-dev libXpm-dev libXft-dev libXext-dev \
  openjdk-7-jdk \
  openjdk-7-doc \
  openjdk-7-jre-lib && \
  apt-get clean && \

  ## clean up
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /tmp/* && \
  rm -rf /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

# set JAVA_HOME
  ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java

#open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]
