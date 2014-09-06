#!/bin/sh

docker ps -a  | grep 'Exit' | awk '{print $1}' | xargs -r docker rm

docker images -a  | grep none | awk '{print $3}' | xargs -r docker rmi
