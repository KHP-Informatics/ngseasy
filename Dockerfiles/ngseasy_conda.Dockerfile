FROM debian:8.5

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
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    # configure locales
    && dpkg-reconfigure locales && \
    locale-gen C.UTF-8 && \
    /usr/sbin/update-locale LANG=C.UTF-8 && \
    echo 'en_US.UTF-8 UTF-8' >> /etc/locale.gen && \
    locale-gen

# Set default locale for the environment
ENV LC_ALL C.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# get ngseasy_conda_install script
COPY ngseasy_conda_install.sh /opt

# install Anaconda and NGS tools
RUN /bin/bash /opt/ngseasy_conda_install.sh /opt && \
# Ensure permissions are set for update in place by arbitrary users
# From: https://github.com/chapmanb/bcbio-nextgen/blob/master/Dockerfile#L68
  find /usr/local -perm /u+x -execdir chmod a+x {} \; && \
  find /usr/local -perm /u+w -execdir chmod a+w {} \; && \
  chown -R ngseasy:ngseasy /home/ngseasy/


  RUN apt-get install -y curl grep sed dpkg && \
      TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
      curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
      dpkg -i tini.deb && \
      rm tini.deb && \
      apt-get clean

  ENV PATH /opt/conda/bin:$PATH

  ENTRYPOINT [ "/usr/bin/tini", "--" ]
  CMD [ "/bin/bash" ]


# Expose for future
EXPOSE 80
EXPOSE 8080

# starting point
CMD [ "/bin/bash" ]
