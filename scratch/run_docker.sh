#!/bin/bash

# playing and blah

docker run --rm -it \
-v ${HOME}/ngs_projects:/home/pipeman/ngs_projects: \
-w /home/pipeman \
-e HOME=/home/pipeman \
-e USER=pipeman \
--name ngseasy_base \
--user=pipeman compbio/ngseasy-base:1.0