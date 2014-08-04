#!/bin/sh
docker build -v /home/admin/temp:/tmp --force-rm=true --no-cache --rm -t ngseasy/var-anno:1.0 .
