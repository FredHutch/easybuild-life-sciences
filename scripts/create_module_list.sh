#!/bin/bash

# create_module_list.sh
#
# Create a list of LMOD modules to publish a Software Inventory
#   modules all or Bio 
#

function usage {
   echo usage: create_module_list.sh [all, bio] lable
   echo Lable should be [chorus, gizmo, ermine, etc]
   echo Lable for Nobel on Gizmo is 'skylake'
   exit
}

repo='easybuild-life-sciences'

if [[ ! -z "${PWD##*${repo}*}" ]]; then
    echo "Can not find github pages docs directory."
    echo "Run script from github repo: ${repo}"
    exit 1
fi

label=""
if [[ $# -eq 2 ]]; then
   inventory_type="${1}"
   label="${2}"
   echo Lable: $label
else
  usage
fi

case $inventory_type in 
  "bio") moduleclass="ai bio chem data math"
  ;;
  "all") moduleclass="all"
  ;;
  *)
  usage
  ;;
esac

# get VERSION_ID from /etc/os-release 
. /etc/os-release

# remove the $repo from left everything else; add the $repo and the "/docs" dir back
base_dir=${PWD%${repo}*}${repo}
docs_dir=${base_dir}/docs/sw_inventory
scripts_dir=${base_dir}/scripts

# Location of LMOD modules
spider=/app/lmod/lmod/libexec/spider
module_dir=/app/modules

echo Collecting Inventory
cd $base_dir
module_search_path=''
for class in $moduleclass; do
   if [ -z $module_search_path ]; then
       module_search_path=${module_dir}/${class}
   else
       module_search_path=${module_search_path}:${module_dir}/${class}
   fi
done
echo $module_search_path
json_in=${docs_dir}/${label}-${inventory_type}-modules-${VERSION_ID}.json
$spider -o spider-json $module_search_path | \
    python3 -mjson.tool | \
    sed 's/4.release-/4./' |\
    sed 's/004.\*release/004/' > ${json_in}
    # 04.*release.

echo Generating Markdown
md_file=${label}-${inventory_type}-modules-${VERSION_ID}
md_out=${docs_dir}/${md_file}.md
csv_out=${docs_dir}/${md_file}.csv

echo '---' > ${md_out}
echo "title: $label Bio Modules" $VERSION_ID >> ${md_out}
echo 'layout: single' >> ${md_out}
echo "permalink: /sw_inventory/${md_file}/" >> ${md_out}
echo 'created: '`date +"%Y-%m-%d"` >> ${md_out}
echo 'toc: true' >> ${md_out}
echo 'toc_label: "On This Page"' >> ${md_out}
echo 'sidebar:' >> ${md_out}
echo '  nav: "docs"' >> ${md_out}
echo '---' >> ${md_out}
echo '' >> ${md_out}

# convert json modules spider to Markdown
cat ${json_in} | ${scripts_dir}/spider2post.py >> ${md_out}
cat ${json_in} | ${scripts_dir}/spider2csv.py > ${csv_out}

echo Wrote inventory to ${md_out}

