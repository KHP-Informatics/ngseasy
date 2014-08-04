#!/bin/sh

docker ps -a | grep 'Exit' | awk '{print }' | xargs -r docker rm

docker images | grep none | awk '{print }' | xargs -r docker rmi
