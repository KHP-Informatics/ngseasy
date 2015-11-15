# base image
FROM compbio/debian:r1.0-002
# Maintainer
MAINTAINER Stephen Newhouse stephen.j.newhouse@gmail.com
LABEL Description="This is the base image for http://www.nextflow.io/" Version="r1.0-002"
RUN curl -fsSL get.nextflow.io | bash
