FROM  continuumio/anaconda3:4.0

MAINTAINER Stephen J Newhouse <stephen.j.newhouse@gmail.com>

RUN apt-get update --fix-missing && \
    apt-get install -y \
    curl \
    wget \
    make \


    apt-get autoremove -y && \
    apt-get autoclean && \
    apt-get clean && \
    apt-get purge && \



# clean up
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get clean && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /tmp/* && \
    rm -rf /var/tmp/* && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/


## Taken from baseimage ontinuumio/anaconda3:4.0
## ENTRYPOINT [ "/usr/bin/tini", "--" ]
## CMD [ "/bin/bash" ]
