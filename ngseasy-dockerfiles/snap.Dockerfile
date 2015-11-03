# Base image
FROM compbio/ngseasy-base:1.0-r001

# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

# Set correct environment variables.
# ENV HOME /root
# ENV DEBIAN_FRONTEND noninteractive

# Update
RUN apt-get update -y && apt-get upgrade -y

# + snap
RUN cd /usr/local/pipeline && \
	git clone https://github.com/amplab/snap.git && \
	chmod -R 777 snap && \
	cd snap && \
	make && \
	chmod -R 777 /usr/local/pipeline && \
	cp -v snap* /usr/local/bin

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 777 /usr/local/pipeline
RUN chown -R ngseasy:ngseasy /usr/local/pipeline

#---------------------------------------------------------------------
#Cleanup the temp dir
RUN rm -rf /tmp/*

#open ports private only
EXPOSE 8080

# Use baseimage-docker's bash.
CMD ["/bin/bash"]

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get autoclean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/
