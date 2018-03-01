#!/bin/bash

set -x
set -e

# variables used: EB_NAME

# try to preserve group write here
umask 002

# load modules
source /app/lmod/lmod/init/bash
module use /app/easybuid/modules/all

# load Easybuild
module load EasyBuild

# build the easyconfig file
eb -l ${EB_NAME}.eb --robot
