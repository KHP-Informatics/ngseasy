#!/bin/bash

# playing and blah

docker run --rm -it \
-v $(HOME)/ngs_projects:/home/pipeman \
-w /home/pipeman \
-e HOME=/home/pipeman \
-e USER=pipeman \
--user=pipeman compbio/ngseasy-base:1.0