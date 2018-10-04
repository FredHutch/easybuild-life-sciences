#!/bin/bash

# build basic template for Markdown Post for new easybuild easyconfigs
# John Dey Sept 2018

doc_date='2018-10-03'
eb_dir='../fh_easyconfigs/'
doc_dir='../docs/_posts/'

list='VEP-94-foss-2016b-Perl-5.24.1.eb
Bio-DB-HTS-2.11-foss-2016b-Perl-5.24.0.eb
MariaDB-10.1.17-foss-2016b.eb
DBD-mysql-4.033-foss-2016b-Perl-5.24.0.eb
BioPerl-2.1.9-foss-2016b-Perl-5.24.0.eb'

for pkg in $list; do
  short=`echo $pkg | sed 's/-.*//'`
  pp=${pkg%.eb}
  pkg_name=${pp/-//}
  cc=`grep moduleclass ${eb_dir}${pkg}`
  ccc=${cc%\'*}
  module_class=${ccc#*\'}
  # locate home URL from easyconfig
  gg=`grep homepage ${eb_dir}${pkg}`
  ggg=${gg%\'*}
  url=${ggg#*\'}
  # get description from easyconfig
  dd=`awk -v RS= -v FPAT="'''.*'''|"'""".*"""' '{print $1}' ${eb_dir}${pkg}`
  description=`echo $dd | sed 's/["'\'']["'\'']["'\'']//g'`
cat << EOF  >${doc_dir}${doc_date}-${short}.md
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

done
exit
