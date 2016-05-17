
#FROM debian:jessie
FROM java:8-jdk


# MAINTAINER Kamil Kwiek <kamil.kwiek@continuum.io>
# sjn additions for local dev

MAINTAINER Stephen J Newhouse (sjn) <stephen.j.newhouse@gmail.com>
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update --fix-missing && \
    apt-get install -y && \
    apt-utils && \
    apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    apt-get install -y && \
    wget \
    bzip2 \
    ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion \
    ## sjn additions
    make \
    cmake \
    curl \
    build-essential \
    unzip \
    xz-utils && \
    rm -rf /var/lib/apt/lists/*

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/archive/Anaconda2-4.0.0-Linux-x86_64.sh && \
    /bin/bash /Anaconda2-4.0.0-Linux-x86_64.sh -b -p /opt/conda && \
    rm /Anaconda2-4.0.0-Linux-x86_64.sh && \
    mkdir -p /home/sjn

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENV PATH /opt/conda/bin:$PATH

# sjn : clean up
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/* && \
    rm -rf /var/tmp/* && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# sjn : volumes added for use in kitematic etc
VOLUME /
VOLUME /home
VOLUME /opt

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
