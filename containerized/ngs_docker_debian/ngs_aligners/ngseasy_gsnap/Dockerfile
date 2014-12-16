# Base image
FROM compbio/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y

## http://research-pub.gene.com/gmap/
# GMAP: A Genomic Mapping and Alignment Program for mRNA and EST Sequences, and
# GSNAP: Genomic Short-read Nucleotide Alignment Program
RUN cd  /usr/local/pipeline && \
  wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-10-09.tar.gz && \
  tar -xvf gmap-gsnap-2014-10-09.tar.gz && \
  rm gmap-gsnap-2014-10-09.tar.gz && \
  cd gmap-2014-10-09 && \
  ./configure && \
  make && \
  chmod -R 755 /usr/local/pipeline/gmap-2014-10-09 && \
  chown -R pipeman:ngsgroup /usr/local/pipeline/gmap-2014-10-09 && \
  sed -i '$aPATH=${PATH}:/usr/local/pipeline/gmap-2014-10-09/src' /home/pipeman/.bashrc && \
  sed -i '$aPATH=${PATH}:/usr/local/pipeline/gmap-2014-10-09/util' /home/pipeman/.bashrc && \
  echo "alias ngsGMAP='/usr/local/pipeline/gmap-2014-10-09'" >>  /home/pipeman/.bashrc

RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline
 
# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

