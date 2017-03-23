#! /bin/bash

if [[ -d easybuild/easyconfigs ]]; then
  rsync -av /app/easybuild/ebfiles_repo/ easybuild/easyconfigs/
  cp  /app/easybuild/fh_easyconfigs/*.patch easybuild/easyconfigs/
  echo "deleting unwanted easyconfigs ..."
  oldstuff=".*\(2014a\|2014b\|2015a\|2015b\|goolf\|ictce\|iimpi\|ifort\|icc-\|CrayGNU\|iomkl\|gimkl\).*.eb"
  find easybuild/easyconfigs -regex "${oldstuff}"  -exec rm {} \;
  ./sandbox/osdep.py
  mv os-dependencies.* ./sandbox
  git add easybuild
else
  echo "folder easybuild/easyconfigs does not exist. Please switch to the root of the easybuild-life-sciences repos."
fi
