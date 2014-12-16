# NGSeasy Base Image

# FROM 
FROM debian:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Remain current
RUN apt-get update && \
  apt-get dist-upgrade -y

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
  libblas-dev \
  libatlas-dev \
  libbz2-dev \
  liblzma-dev \
  libpcre3-dev \
  libreadline-dev \
  zlib1g-dev \  
  make \
  mercurial \
  perl \
  php5-curl \
  python \
  sudo \
  subversion \
  tabix \
  tree \
  unzip \
  vim \
  wget \
  python-software-properties

# Remain current and upgrade java
RUN apt-get update -y && \
  apt-get upgrade -y && \
  apt-get purge -y openjdk-6-jdk openjdk-6-doc openjdk-6-jre-lib && \
  apt-get install -y openjdk-7-jdk openjdk-7-doc openjdk-7-jre-lib && \
  apt-get update -y && \
  apt-get upgrade -y

#set java
ENV JAVA_HOME /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java
RUN sed -i 'aPATH=$PATH:/usr/lib/jvm/java-7-openjdk-amd64/jre/bin' /root/.bashrc

# Create a pipeline user:pipeman and group:ngsgroup
RUN useradd -m -s /bin/bash pipeman && cd /home/pipeman && echo "#bash config file for user pipeman" >> /home/pipeman/.bashrc

RUN groupadd ngsgroup && \
  usermod -aG ngsgroup pipeman && \
  usermod -aG sudo pipeman

#make pipeline install dirs
RUN mkdir /usr/local/pipeline && chown pipeman:ngsgroup /usr/local/pipeline && chmod 775 /usr/local/pipeline

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



