#!/bin/bash

# build basic template for Markdown Post for new easybuild easyconfigs
# John Dey Sept 2018

date='2018-09-21'
eb_dir='../fh_easyconfigs/'
doc_dir='../docs/_post/'

list='BCFtools-1.9-foss-2018b.eb
BEDTools-2.27.1-foss-2018b.eb
BamTools-2.5.1-foss-2018b.eb
HTSlib-1.9-foss-2018b.eb
SAMtools-1.9.0-foss-2018b.eb'

for pkg in $list; do
  short=`echo $pkg | sed 's/-.*//'`
  pkg_name=${pkg%.eb}
  # locate home URL from easyconfig
  gg=`grep homepage ${eb_dir}${pkg}`
  ggg=${gg%\'*}
  url=${ggg#*\'}
  # get description from easyconfig
  description=`awk -v RS= -v FPAT="'''.*'''|"'""".*"""' '{print $1}' ${eb_dir}${pkg}`
cat >${doc_dir}${date}-${short}.md <<EOF
---
layout: post
title: $pkg_name
catagory: bio
homepage: $url
___
$description
EOF
done
