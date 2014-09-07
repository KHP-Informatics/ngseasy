# FROM 
FROM phusion/baseimage:0.9.13

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Regenerate SSH host keys. baseimage-docker does not contain any, so you
# have to do that yourself. You may also comment out this instruction; the
# init system will auto-generate one during boot.
RUN /etc/my_init.d/00_regen_ssh_host_keys.sh

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update

# basic dependencies
RUN apt-get install -y build-essential  gcc g++ \
python-dev make dkms cmake tabix dpkg-dev binutils gfortran git-core wget unzip bzip2 zlib1g-dev \
libncurses5-dev sysvbanner openjdk-7-jre openjdk-7-jdk \
curl libcurl3 php5-curl libxml2 libxslt-dev libxml2-dev fastx-toolkit htmldoc html2text perl bioperl cron  libzip-dev libbz2-dev  \
libcurl4-openssl-dev  libcurl3-gnutls libjpeg62 libjpeg62-dev libpng12-0 libpng12-dev libmcrypt-dev \
libmcrypt4 libxslt1-dev mercurial subversion vim tree libssl1.0.0 libssl-dev libcrypto++9 libcrypto++-dev libbam-dev \
libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev \
git dpkg-dev make g++ gcc binutils libx11-dev libxpm-dev libxft-dev libxext-dev libmpc3 libmpc-dev \
gfortran libssl-dev libpcre3-dev xlibmesa-glu-dev libglew1.5-dev libftgl-dev libmysqlclient-dev libfftw3-dev cfitsio-dev \
graphviz-dev libavahi-compat-libdnssd-dev libldap2-dev python-dev libxml2-dev libkrb5-dev libgsl0-dev libqt4-dev c++11 libroot-dev root-system-common automake;

RUN apt-get install -y r-base-dev 

RUN apt-get update -y
RUN apt-get upgrade -y

RUN apt-get install -y  linux-libc-dev

# Create a pipeline user:pipeman and group:ngsgroup
RUN mkdir /usr/local/pipeline \
  && chmod 766 /usr/local/pipeline
  
RUN useradd -m -s /bin/bash pipeman && cd /home/pipeman && echo "#bash config file for user pipeman" >> /home/pipeman/.bashrc
RUN groupadd ngsgroup
RUN usermod -G ngsgroup pipeman
RUN adduser pipeman sudo
RUN usermod -a -G sudo pipeman


# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*



