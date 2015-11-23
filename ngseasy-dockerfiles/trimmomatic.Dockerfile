# Base image
FROM compbio/debian:small-0.1-0bef85c

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Trimmomatic
RUN TRIMMOMATIC_VERSION="0.32" && \
  cd /usr/local/ngs/bin && \
  wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \
  unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \
  chmod -R 777 Trimmomatic-${TRIMMOMATIC_VERSION}/* && \
  sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/Trimmomatic-${TRIMMOMATIC_VERSION}/trimmomatic-${TRIMMOMATIC_VERSION}.jar' /home/ngseasy/.bashrc && \
  sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/Trimmomatic-${TRIMMOMATIC_VERSION}' /home/ngseasy/.bashrc && \
  cp -v /usr/local/ngs/bin/Trimmomatic-${TRIMMOMATIC_VERSION}/trimmomatic-${TRIMMOMATIC_VERSION}.jar /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip && \
  bash -c "source /home/ngseasy/.bashrc" && \

# Clean up APT when done.
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

USER ngseasy
WORKDIR /home/ngseasy

CMD ["java","-jar","trimmomatic-${TRIMMOMATIC_VERSION}.jar","--help"]
