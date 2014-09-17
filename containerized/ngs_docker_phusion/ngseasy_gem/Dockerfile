# Base image
FROM snewhouse/ngseasy-base:v1.0

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Regenerate SSH host keys. baseimage-docker does not contain any, so you
# have to do that yourself. You may also comment out this instruction; the
# init system will auto-generate one during boot.
RUN /etc/my_init.d/00_regen_ssh_host_keys.sh

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update

# NOTE: ADD automatically unpacks compressed files

RUN cd /usr/local/pipeline \
	&& git clone --branch=develop git://github.com/samtools/htslib.git \
	&& git clone --branch=develop git://github.com/samtools/bcftools.git \
	&& git clone --branch=develop git://github.com/samtools/samtools.git \
	&& cd /usr/local/pipeline/bcftools \
	&& make \
	&& cd /usr/local/pipeline/samtools \
	&& make \
	&& cd /usr/local/pipeline/htslib \
	&& make \
	&& chown -R pipeman:ngsgroup /usr/local/pipeline/samtools \
	&& chown -R pipeman:ngsgroup /usr/local/pipeline/bcftools \
	&& chown -R pipeman:ngsgroup /usr/local/pipeline/htslib \
	&& sed  -i '$aPATH=${PATH}:/usr/local/pipeline/samtools' /home/pipeman/.bashrc \
	&& echo "alias ngsSAMtools='/usr/local/pipeline/samtools'" >>  /home/pipeman/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/bcftools' /home/pipeman/.bashrc \
        && echo "alias ngsBCFtools='/usr/local/pipeline/bcftools'" >>  /home/pipeman/.bashrc \
        && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/htslib' /home/pipeman/.bashrc \
        && echo "alias ngsHTSlib='/usr/local/pipeline/htslib'" >>  /home/pipeman/.bashrc 


RUN wget -O /tmp/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2 http://sourceforge.net/projects/gemlibrary/files/gem-library/Binary%20pre-release%203/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2 \
  && tar xjvf /tmp/GEM-binaries-Linux-x86_64-core_i3-20130406-045632.tbz2 -C /usr/local/pipeline/ \
  && chmod -R 755 /usr/local/pipeline/GEM-binaries-Linux-x86_64-core_i3-20130406-045632 \
  && chown -R pipeman:ngsgroup /usr/local/pipeline/GEM-binaries-Linux-x86_64-core_i3-20130406-045632 \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin' /home/pipeman/.bashrc \
  && echo "alias ngsGEM='/usr/local/pipeline/GEM-binaries-Linux-x86_64-core_i3-20130406-045632/bin'" >>  /home/pipeman/.bashrc




RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline
 
# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

