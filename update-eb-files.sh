#! /bin/bash

if [[ -d easybuild/easyconfigs ]]; then
  rsync -av /app/easybuild/ebfiles_repo/ easybuild/easyconfigs/
  mkdir -p easybuild/easyconfigs/-PATCHES
  cp  /app/easybuild/fh_easyconfigs/*.patch easybuild/easyconfigs/-PATCHES
  echo "deleting unwanted easyconfigs ..."
  oldstuff=".*\(2014a\|2014b\|2015a\|2015b\|goolf\|ictce\|iimpi\|ifort\|icc-\|CrayGNU\|iomkl\|gimkl\).*.eb"
  find easybuild/easyconfigs -regex "${oldstuff}"  -exec rm {} \;
else
  echo "folder easybuild/easyconfigs does not exist."
fi
