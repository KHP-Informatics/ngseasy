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

 
#---------------------------------sam/vcf processing-------------------------------------#
# + samtools, htslib and bcftools
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

#-----------------------------VarScan2-------------------------------
RUN cd /usr/local/pipeline && \
	mkdir varscan && \
	cd varscan && \
	wget http://sourceforge.net/projects/varscan/files/VarScan.v2.3.7.jar


#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

