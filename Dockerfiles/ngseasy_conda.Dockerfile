FROM debian:stable

MAINTAINER Stephen Newhouse <stephen.j.newhouse@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

# update and get some things
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    git \
    mercurial \
    subversion \
    make \
    cmake \
    gcc \
    gpp \
    build-essential \
    zlib1g-dev \
    patch \
    perl \
    libextutils-pkgconfig-perl \
    byacc \
    wget \
    curl \
    grep \
    sed \
    dpkg \
    debconf \
    unzip \
    tree \
    libc6 \
    llvm \
    locales \
    time \
    sudo && \
    apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get autoremove -y && \
    apt-get autoclean && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

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

# User set up
RUN useradd -m -U -s /bin/bash ngseasy && \
  cd /home/ngseasy && \
  chown ngseasy:ngseasy /home/ngseasy && \
  usermod -aG sudo ngseasy && \
# echo 'export PATH=/home/ngseasy/anaconda2/bin:$PATH' > /etc/profile.d/conda.sh && \
  groupadd docker && \
  usermod -aG docker ngseasy

# get ngseasy_conda_install script
COPY ngseasy_conda_install.sh /home/ngseasy/

# install Anaconda and NGS tools
RUN /bin/bash /home/ngseasy/ngseasy_conda_install.sh /home/ngseasy && \
# Ensure permissions are set for update in place by arbitrary users
# From: https://github.com/chapmanb/bcbio-nextgen/blob/master/Dockerfile#L68
  find /usr/local -perm /u+x -execdir chmod a+x {} \; && \
  find /usr/local -perm /u+w -execdir chmod a+w {} \; && \
  chown -R ngseasy:ngseasy /home/ngseasy/ && \
  chmod -R 755 /home/ngseasy/anaconda2
  
# switch to ngseasy user
USER ngseasy

# volumes
RUN mkdir -p /home/ngseasy/resources && \
    mkdir -p /home/ngseasy/reference_genomes && \
    mkdir -p /home/ngseasy/ngs_projects && \
    mkdir -p /home/ngseasy/ngs_projects/fastq && \
    mkdir -p /home/ngseasy/scratch && \
    export "export PATH=$PATH:/home/ngseasy/anaconda2/bin" >> /home/ngseasy/.bashrc && \
    /bin/bash -c "source /home/ngseasy/.bashrc"

# set Home
ENV HOME /home/ngseasy
WORKDIR /home/ngseasy

# Expose for future
EXPOSE 80
EXPOSE 8080

# starting point
CMD [ "/bin/bash" ]
