# Base image
FROM compbio/debian:small-0.1-0bef85c

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# FastQC
RUN wget -O /tmp/fastqc_v0.11.2.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip && \
  unzip /tmp/fastqc_v0.11.2.zip -d /usr/local/ngs/bin/ && \
  chmod -R 755 /usr/local/ngs/bin/ && \
  sed -i '$aCLASSPATH=.:${CLASSPATH}:/usr/local/ngs/bin/FastQC/jbzip2-0.9.jar:/usr/local/ngs/bin/FastQC/sam-1.103.jar' /home/ngseasy/.bashrc && \
  sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/FastQC' /home/ngseasy/.bashrc && \
  ln -s /usr/local/ngs/bin/FastQC/fastqc /usr/local/bin/fastqc && \
  chmod -R 755 /usr/local/ngs/bin && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin && \
  bash -c "source /home/ngseasy/.bashrc"

#Clean up APT when done.
RUN apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
  apt-get autoclean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

USER ngseasy
WORKDIR /home/ngseasy

# base command
# CMD ["/bin/bash"]
CMD ["fastqc", "--help"]
