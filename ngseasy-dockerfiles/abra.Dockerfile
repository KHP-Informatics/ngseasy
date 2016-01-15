FROM compbio/ngseasy-base:1.0
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com

## ABRA - Assembly Based ReAligner https://github.com/mozack/abra
RUN cd /usr/local/ngs/bin && \
    wget https://github.com/mozack/abra/releases/download/v0.94/abra-0.94-SNAPSHOT-jar-with-dependencies.jar && \
    chmod 775 abra-0.94-SNAPSHOT-jar-with-dependencies.jar && \
    mv -v abra-0.94-SNAPSHOT-jar-with-dependencies.jar /usr/local/bin/abra-0.94

## PERMISSIONS
RUN chmod -R 755 /usr/local/ngs/bin
RUN chown -R ngseasy:ngseasy /usr/local/ngs/bin

## Cleanup the temp dir
RUN rm -rvf /tmp/*

## open ports private only
EXPOSE 8080

## Use baseimage-docker's bash.
CMD ["/bin/bash"]

## Clean up APT when done.
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/ && \
    rm -rf /usr/local/ngs/bin/*
USER ngseasy
WORKDIR /home/ngseasy
