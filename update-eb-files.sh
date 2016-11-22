#! /bin/bash

if [[ -d easybuild/easyconfigs ]]; then 
  rsync -av /app/easybuild/ebfiles_repo/ easybuild/easyconfigs/
  find easybuild/easyconfigs -name *2014a*.eb -o -name *2014b*.eb -o -name *2015a*.eb -o -name *2015b*.eb -exec rm {} \;
  find easybuild/easyconfigs -name *goolf*.eb -exec rm {} \;
else
  echo "folder easybuild/easyconfigs does not exist."
fi
