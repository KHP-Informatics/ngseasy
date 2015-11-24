# Base image
FROM compbio/debian:small-0.1-0bef85c

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# dependencies 
RUN apt-get update && apt-get install -y cmake llvm zlib1g zlib1g-dev && \
  apt-get autoremove -y && \
  apt-get autoclean && \
  apt-get clean && \
  apt-get purge && \
# sambamba
  cd /usr/local/ngs/bin && \
  SAMBAMBA_VERSION="v0.5.9" && \
  curl -OL https://github.com/lomereiter/sambamba/releases/download/${SAMBAMBA_VERSION}/sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  tar -xjvf sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  mv sambamba_${SAMBAMBA_VERSION} sambamba && \
  chmod +rwx sambamba && \
  cp -v sambamba /usr/local/bin/ && \
  cd /usr/local/ngs/bin && \
  rm sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  rm /usr/local/ngs/bin/sambamba && \

# vcftools https://github.com/vcftools/vcftools
  cd /usr/local/ngs/bin/ && \
  git clone https://github.com/vcftools/vcftools.git && \
  cd vcftools && \
  ./autogen.sh && \
  ./configure && \
  make && \
  make install && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vcftools && \

# vcflib
  cd /usr/local/ngs/bin/ && \
  rm -rfv /usr/local/ngs/bin/vcflib && \
  git clone --recursive git://github.com/ekg/vcflib.git && \
  cd vcflib && \
  make && \
  chmod -R 777 ./bin/ && \
  cp -v ./bin/* /usr/local/bin/ && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vcflib && \

# bedtools
  cd /usr/local/ngs/bin && \
  git clone https://github.com/arq5x/bedtools2.git && \
  cd bedtools2 && \
  make all && \
  chmod -R 777 ./* && \
  make install && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./bedtools2 && \

# vt
  cd /usr/local/ngs/bin && \
  git clone https://github.com/atks/vt.git && \
  chmod -R 777 vt/ && \
  cd vt && \
  make && \
  chmod -R 777 vt && \
  cp -v vt /usr/local/bin && \
  cd /usr/local/ngs/bin/ && \
  rm -r ./vt && \

## Variant caller
# freebayes
  cd /usr/local/ngs/bin \
  && git clone --recursive git://github.com/ekg/freebayes.git \
  && cd /usr/local/ngs/bin/freebayes \
  && make \
  && chmod -R 777 /usr/local/ngs/bin/freebayes \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/freebayes/bin' /home/ngseasy/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/ngs/bin/freebayes/bin' ~/.bashrc \
  && cp -v /usr/local/ngs/bin/freebayes/bin/* /usr/local/bin && \

# source
  bash -c "source /home/ngseasy/.bashrc"

ADD fix_ambiguous /usr/local/bin/

# user and wd
USER ngseasy
WORKDIR /home/ngseasy

# command
CMD ["freebayes"]
