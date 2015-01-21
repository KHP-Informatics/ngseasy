#!/bin/bash

# playing and blah

docker run --rm -it \
-v ${HOME}/ngs_projects:/home/pipeman/ngs_projects: \
-w /home/pipeman \
-e HOME=/root \
--name ngseasy_test compbio/ngseasy-stampy:1.0