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


# + novoalign (registration required,  get compressed file and put in context dir for the build)
# + get novoalign.lic from novoalign ~ $1000 and put in context dir for the build
# + get updated version of novosort novosortV1.03.01.Linux3.0.tar.gz and put in context dir for the build
    RUN wget -O /tmp/novocraftV3.02.07.Linux3.0.tar.gz https://www.dropbox.com/s/kbgds9e6ucog6oj/novocraftV3.02.07.Linux3.0.tar.gz \
      && tar xvf /tmp/novocraftV3.02.07.Linux3.0.tar.gz  -C /usr/local/pipeline/ \
      && cp -v /usr/local/pipeline/novocraft/novosort /usr/local/pipeline/novocraft/novosort_original \
      && rm /usr/local/pipeline/novocraft/novosort \
      && wget -O /tmp/novosortV1.03.01.Linux3.0.tar.gz https://www.dropbox.com/s/1j4w9rpyoug1pbo/novosortV1.03.01.Linux3.0.tar.gz \
      && tar xvf /tmp/novosortV1.03.01.Linux3.0.tar.gz -C /usr/local/pipeline/novocraft \
      && sed  -i '$aPATH=${PATH}:/usr/local/pipeline/novocraft' /home/pipeman/.bashrc \
      && chown -R pipeman:ngsgroup /usr/local/pipeline/novocraft \
      && echo "alias ngsNovoalign='/usr/local/pipeline/novocraft/novoalign'" >>  /home/pipeman/.bashrc

# + get novoalign.lic 
ADD novoalign.lic /usr/local/pipeline/novocraft/
RUN chown -R pipeman:ngsgroup /usr/local/pipeline/novocraft
 

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

