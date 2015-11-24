# Base image
FROM compbio/debian:small-0.1-0bef85c

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# bwakit
# http://sourceforge.net/projects/bio-bwa/files/bwakit/
# http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.12_x64-linux.tar.bz2
RUN cd /usr/local/ngs/bin && \
  BWA_VERSION="0.7.12" && \
  wget http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-${BWA_VERSION}_x64-linux.tar.bz2 && \
  tar xjvf bwakit-${BWA_VERSION}_x64-linux.tar.bz2 && \
  chmod -R 777 /usr/local/ngs/bin && \
  ln -s /usr/local/ngs/bin/bwa.kit/bwa /usr/local/bin/bwa && \
  cd /usr/local/ngs/bin && \
  rm bwakit-${BWA_VERSION}_x64-linux.tar.bz2 && \
  chown -R ngseasy:ngseasy /usr/local/ngs/bin && \

# samblaster
  cd /usr/local/ngs/bin && \
  git clone git://github.com/GregoryFaust/samblaster.git && \
  cd samblaster && \
  make && \
  chmod -R 777 /usr/local/ngs/bin && \
  cp -v samblaster /usr/local/bin/ && \
  cd /usr/local/ngs/bin && \
  rm -r /usr/local/ngs/bin/samblaster && \

# sambamba
  cd /usr/local/ngs/bin && \
  SAMBAMBA_VERSION="v0.5.1" && \
  curl -OL https://github.com/lomereiter/sambamba/releases/download/${SAMBAMBA_VERSION}/sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  tar -xjvf sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  mv sambamba_${SAMBAMBA_VERSION} sambamba && \
  chmod +rwx sambamba && \
  cp -v sambamba /usr/local/bin/ && \
  cd /usr/local/ngs/bin && \
  rm sambamba_${SAMBAMBA_VERSION}_linux.tar.bz2 && \
  rm /usr/local/ngs/bin/sambamba && \

# source
  bash -c "source /home/ngseasy/.bashrc"

# user and wd
USER ngseasy
WORKDIR /home/ngseasy

# command
CMD ["bwa","mem"]
