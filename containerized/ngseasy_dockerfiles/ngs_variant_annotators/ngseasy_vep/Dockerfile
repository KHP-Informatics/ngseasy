# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Remain current
RUN apt-get update && \
  apt-get upgrade -y

# VEP GRCh37
RUN cd /usr/local/pipeline && \ 
  wget https://github.com/Ensembl/ensembl-tools/archive/release/78.zip && \
  unzip 78.zip && \
  chmod -R 777 /usr/local/pipeline/* && \
  cd /usr/local/pipeline/ensembl-tools-release-78/scripts/variant_effect_predictor/ && \
  perl INSTALL.pl -a ac -s homo_sapiens -y GRCh37 --CACHEDIR ~/.vep

# Ports
EXPOSE 80

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


############################################################################################
# Notes: 
#
# If using human GRCh37 add "--port 3337" to use the GRCh37 database, or --offline to avoid database connection entirely
#
# -h | --help        Display this message and quit
#-d | --DESTDIR     Set destination directory for API install (default = './')
#-v | --VERSION     Set API version to install (default = 78)
#-c | --CACHEDIR    Set destination directory for cache files (default = '/.vep/')
#-n | --UPDATE      EXPERIMENTAL! Check for and download new VEP versions
#-a | --AUTO        Run installer without user prompts. Use a (API), c (cache),
#                  f (FASTA) to specify parts to install e.g. -a ac for API and
#                   cache
#-s | --SPECIES     Comma-separated list of species to install when using --AUTO
#-y | --ASSEMBLY    Assembly name to use if more than one during --AUTO
#-q | --QUIET       Don't write any status output when using --AUTO
#-p | --PREFER_BIN  Use this if the installer fails with out of memory errors
#-t | --CONVERT     Convert downloaded caches to use tabix for retrieving
#                   co-located variants (requires tabix)
#-u | --CACHEURL    Override default cache URL; this may be a local directory or
#                   a remote (e.g. FTP) address.
#-f | --FASTAURL    Override default FASTA URL; this may be a local directory or
#                   a remote (e.g. FTP) address. The FASTA URL/directory must have
#                   gzipped FASTA files under the following structure:
#                   [species]/[dna]/
