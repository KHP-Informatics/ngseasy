# Base image
FROM compbio/debian:small-0.1-0bef85c

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Picard
# htsjdk-1.141.jar  libIntelDeflater.so  picard.jar  picard-lib.jar
RUN PICARD_VERSION="1.141" && \
  cd /usr/local/ngs/bin && \
  wget https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip && \
  unzip picard-tools-${PICARD_VERSION}.zip && \
  chmod -R 777 /usr/local/ngs/bin/picard-tools-${PICARD_VERSION} && \
  mv -v /usr/local/ngs/bin/picard-tools-${PICARD_VERSION}/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r picard-tools-${PICARD_VERSION} && \
  rm -r picard-tools-${PICARD_VERSION}.zip && \
  bash -c "source /home/ngseasy/.bashrc" && \

# Clean up APT when done.
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

USER ngseasy
WORKDIR /home/ngseasy

CMD ["java","-jar",""/usr/local/bin/picard.jar"]
