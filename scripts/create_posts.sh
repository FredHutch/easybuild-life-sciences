#!/bin/bash

# create_posts.sh
#
# create Jekyll post from Easyconfig file 
# usage: eb_file
#
# John Dey Sept 2018

if [ "$#" -ne 1 ]; then
    echo "Must supply one easybuild-easyconfig as an argument"
    exit 1
fi
pkg=$1
doc_date=`date +'%Y-%m-%d'`

# where is the script being run from? The script can only be run from inside
# the "easybuild-life-sciences" repo, but it can be run from any 
# location within the repo. 
repo='easybuild-life-sciences'
if [[ $PWD = *${repo}* ]]; then
   base_dir=${PWD%${repo}*}${repo}
else
   echo 'can only be run from inside the ${repo} repo'
   exit 1
fi
docs_dir=${base_dir}/docs
posts_dir=${base_dir}/docs/_posts/
scripts_dir=${base_dir}/scripts
eb_dir=${base_dir}/fh_easyconfigs/

  short=`echo $pkg | sed 's/-.*//'`
  pp=${pkg%.eb}
  pkg_name=${pp/-//}
  cc=`grep moduleclass ${eb_dir}/${pkg}`
  ccc=${cc%\'*}
  module_class=${ccc#*\'}
  # locate home URL from easyconfig
  gg=`grep homepage ${eb_dir}${pkg}`
  ggg=${gg%\'*}
  url=${ggg#*\'}
  # get description from easyconfig
  dd=`awk -v RS= -v FPAT="'''.*'''|"'""".*"""' '{print $1}' ${eb_dir}${pkg}`
  description=`echo $dd | sed 's/["'\'']["'\'']["'\'']//g'`
cat << EOF  >${posts_dir}${doc_date}-${short}.md
---
layout: post
title: $short
catagory: $module_class 
homepage: $url
---
$description
\`\`\`
module load $pkg_name
\`\`\`
EOF

exit
