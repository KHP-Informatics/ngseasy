#!/usr/bin/env bash

TAG=$1

docker build \
--rm=true \
--tag=compbio/ngseasy-box:${TAG} \
--file=ngseasy-box.Dockerfile .


