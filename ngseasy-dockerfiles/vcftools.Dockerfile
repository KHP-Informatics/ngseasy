# Base image
FROM compbio/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Update
RUN apt-get update -y

# + VCFtools: http://vcftools.sourceforge.net/index.html
RUN wget -O /tmp/vcftools_0.1.12b.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12b.tar.gz \
    && tar xzvf /tmp/vcftools_0.1.12b.tar.gz -C /usr/local/ngs/bin/  \
    && cd /usr/local/ngs/bin/vcftools_0.1.12b/ && make \
    && chown -R ngseasy:ngseasy /usr/local/ngs/bin/vcftools_0.1.12b \
    && sed  -i '$aPATH=${PATH}:/usr/local/ngs/bin/vcftools_0.1.12b/bin' /home/ngseasy/.bashrc \
    && echo "alias ngsVCFtools='/usr/local/ngs/bin/vcftools_0.1.12b/bin/'" >>  /home/ngseasy/.bashrc    

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/ngs/bin/***
RUN chown -R ngseasy:ngseasy /usr/local/ngs/bin

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
USER ngseasy
WORKDIR /home/ngseasy
