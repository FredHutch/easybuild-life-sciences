#!/bin/bash

# script to download sources too large to fit into a git repo
curl -o /app/sources/jdk-8u121-linux-x64.tar.gz https://s3-us-west-2.amazonaws.com/ls2-sources/jdk-8u121-linux-x64.tar.gz
curl -o /app/sources/tensorflow-1.5.0rc0-cp36-cp36m-linux_x86_64.whl https://s3-us-west-2.amazonaws.com/ls2-sources/tensorflow-1.5.0rc0-cp36-cp36m-linux_x86_64.whl
