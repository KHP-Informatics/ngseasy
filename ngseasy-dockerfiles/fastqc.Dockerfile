# base image
FROM compbio/ngseasy-debian:small

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# fastc
RUN  FASTQC_VERSION="0.11.4" && \
  cd /usr/local/ngs/bin && \
  wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip && \
  unzip fastqc_v${FASTQC_VERSION}.zip && \
  chmod -R 755 ./FastQC/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r fastqc_v${FASTQC_VERSION}.zip && \
  bash -c "source /home/ngseasy/.bashrc"

# Clean up APT when done.
RUN apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

# user and env vars
ENV PATH /usr/local/ngs/bin/FastQC:$PATH
ENV HOME /home/ngseasy
USER ngseasy
WORKDIR /home/ngseasy
VOLUME /home/ngseasy/ngs_projects

# base command
CMD ["/bin/bash"]
