# Base image
FROM compbio/ngseasy-base:wheezy

# Maintainer 
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
ENV HOME /root

# Use baseimage-docker's init system.
CMD ["bash"]

# Update
RUN apt-get update -y

#-------------------------------TOOL--------------------------

# + SAVANT
<<<<<<< HEAD
wget -O /tmp/Savant-full-latest.tgz http://www.well.ox.ac.uk/bioinformatics/Software/Savant-full-latest.tgz \
  && tar xvf /tmp/Savant-full-latest.tgz -C /usr/local/pipeline \
  && cd /usr/local/pipeline/savant-v1.1.0-complete \
  && chmod -R 755 /usr/loca/pipeline/savant-v1.1.0-complete \
  && sed -i '$aPATH=${PATH}:$INSTALL_TARGET/savant-v1.1.0-complete' $PATH_FILE \
  && echo "alias ngsSavant='$INSTALL_TARGET/savant-v1.1.0-complete'" >>  $PATH_FILE

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/
=======
RUN wget -O /tmp/Savant-full-latest.tgz http://www.well.ox.ac.uk/bioinformatics/Software/Savant-full-latest.tgz \
  && tar xvf /tmp/Savant-full-latest.tgz -C /usr/local/pipeline/ \
  && cd /usr/local/pipeline/savant-v1.1.0-complete \
  && chmod -R 755 /usr/local/pipeline/savant-v1.1.0-complete \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/savant-v1.1.0-complete' /home/pipeman/.bashrc \
  && echo "alias ngsSavant='/usr/local/pipeline/savant-v1.1.0-complete'" >>  /home/pipeman/.bashrc \
  && sed -i '$aPATH=${PATH}:/usr/local/pipeline/savant-v1.1.0-complete' /root/.bashrc

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 766 /usr/local/pipeline/***
>>>>>>> 7ea2ff7bd9089cf957f5097e745f7b6adf4f873b
RUN chown -R pipeman:ngsgroup /usr/local/pipeline

# Cleanup the temp dir
RUN rm -rf /tmp/*

# open ports private only
EXPOSE 80

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
