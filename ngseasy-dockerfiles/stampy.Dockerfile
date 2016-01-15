# Base image
FROM compbio/ngseasy-base:1.0-r001

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Remain current
RUN apt-get update && apt-get upgrade -y

# + stampy (registration required, get compressed file and put in context dir for the build)

RUN wget -O /tmp/stampy-latest.tgz http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz \
    && tar xvf /tmp/stampy-latest.tgz -C /usr/local/ngs/bin/ \
    && sed -i 's/-Wl//' /usr/local/ngs/bin/stampy-1.0.27/makefile \
    && chmod -R 755 /usr/local/ngs/bin/stampy-1.0.27 \
    && cd /usr/local/ngs/bin/stampy-1.0.27 \
    && make \ 
    && chmod -R 755 /usr/local/ngs/bin/ \
    && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/stampy-1.0.27' /home/ngseasy/.bashrc \
    && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/stampy-1.0.27' ~/.bashrc

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/ngs/bin/***
RUN chown -R ngseasy:ngseasy /usr/local/ngs/bin

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
USER ngseasy
WORKDIR /home/ngseasy
