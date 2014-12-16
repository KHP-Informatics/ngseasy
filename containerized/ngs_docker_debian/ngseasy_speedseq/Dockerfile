# Base image
FROM compbio/ngseasy-base:latest

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

#----------------------------------------------------------------------------------#
# SPEEDSEQ https://github.com/cc2qe/speedseq
#----------------------------------------------------------------------------------#

# get ROOT
RUN curl -OL ftp://root.cern.ch/root/root_v5.34.20.source.tar.gz && \
  tar -xvf root_v5.34.20.source.tar.gz && \
  cd root && \
  ./configure && \
  make && \
  cd ../ && \
  sudo mv -v root /usr/local/ && \
  echo "source /usr/local/root/bin/thisroot.sh" >> ~/.bashrc && \
  . ~/.bashrc

# speedseq 
RUN cd /usr/local/pipeline && \
  git clone --recursive https://github.com/cc2qe/speedseq && \
  cd /usr/local/pipeline/speedseq && \
  python speedseq_setup.py && \
  sudo cp -r bin/* /usr/local/bin/ && \
  cd /usr/local/pipeline/speedseq  && \
  make cnvnator-multi  && \
  sudo cp bin/* /usr/local/bin
  
# chroms data  
RUN cd /usr/local/pipeline/speedseq && \
  mkdir chroms && \
  cd chroms && \
  wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz && \
  wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai && \
  gunzip human_g1k_v37.fasta.gz && \
  cat human_g1k_v37.fasta | awk 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > CHROM".fa" }'  
  
# VEP  
RUN sudo yum -y install "perl(Archive::Extract)" "perl(CGI)" "perl(DBI)" "perl(Time::HiRes)" "perl(Archive::Tar)" "perl(Archive::Zip)"  

RUN cd /usr/local/pipeline && \
  curl -OL https://github.com/Ensembl/ensembl-tools/archive/release/76.zip  && \
  unzip 76.zip  && \
  perl ensembl-tools-release-76/scripts/variant_effect_predictor/INSTALL.pl -a ac -s homo_sapiens -y GRCh37 --CACHEDIR ~/.vep && \
  sudo cp ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl /usr/local/bin && \
  sudo cp -r ensembl-tools-release-76/scripts/variant_effect_predictor/Bio /usr/local/bin
  
# GEMINI  
RUN cd /usr/local/pipeline && \
  wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py && \
  sudo python2.7 gemini_install.py /usr/local /usr/local/share/gemini && \
  export PATH=$PATH:/usr/local/gemini/bin && \
  echo "export PATH=$PATH:/usr/local/gemini/bin" >> ~/.bashrc && \
  . ~/.bashrc && \
  gemini update
  
#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 80

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
  
  
  
  
  
  
  
