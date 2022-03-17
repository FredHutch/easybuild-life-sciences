#!/bin/bash

files='fh_easyconfigs/f/FlashPCA2/FlashPCA2-2.0-GCCcore-11.2.0.eb
fh_easyconfigs/j/JAGS/JAGS-4.3.0-foss-2021b.eb
fh_easyconfigs/l/librsvg/librsvg-2.53.2-GCCcore-11.2.0.eb
fh_easyconfigs/m/MUMPS/
fh_easyconfigs/p/Python/fhPython-3.9.6-foss-2021b.eb
fh_easyconfigs/p/pybedtools/Dockerfile
fh_easyconfigs/q/QUAST/
fh_easyconfigs/s/SCOTCH/'

for f in $files; do
    base=`basename $f`
    name="${base%.*}"
    echo commit $name
    git add $f
    git commit -m "update $name"
    git push
done

