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
RUN cd /usr/local/ngs/bin && \
	git clone https://github.com/amplab/snap.git && \
	chmod -R 755 snap && \
	cd snap && \
	make && \
	chmod -R 755 /usr/local/ngs/bin && \
	cp -v snap* /usr/local/bin

#-------------------------------PERMISSIONS--------------------------
RUN chmod -R 755 /usr/local/ngs/bin
RUN chown -R ngseasy:ngseasy /usr/local/ngs/bin

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
USER ngseasy
WORKDIR /home/ngseasy
