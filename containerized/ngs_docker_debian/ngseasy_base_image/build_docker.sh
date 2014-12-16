#!/bin/bash
repo=${1}
sudo docker build --rm=true -t ${repo}/ngseasy-base:wheezy .
