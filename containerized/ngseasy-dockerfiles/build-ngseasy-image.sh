#!/bin/bash
TOOL=${1}
VERSION=${2}
docker build \
--no-cache=true \
--rm=true \
--tag=compbio/ngseasy-${TOOL}:${VERSION} \
--file=${TOOL}.Dockerfile .
