#!/bin/bash

# create_module_list.sh
#
# Create a list of LMOD modules. Filter modules based on module class 'bio', 'chem', 'phys' and
# 'math'
# Used to create software Inventory

repo='easybuild-life-sciences'
marker='<!---DO NOT EDIT BELOW HERE--->'


if [[ ! -z "${PWD##*${repo}*}" ]]; then
    echo "Can not find github pages docs directory."
    echo "Run script from github repo: ${repo}"
    exit 1
fi

label=""
if [[ $# -eq 1 ]]; then
   label="${1}-"
   echo Lable: $label
fi

# get VERSION_ID from /etc/os-release 
. /etc/os-release

# remove the $repo from left everything else; add the $repo and the "/docs" dir back
base_dir=${PWD%${repo}*}${repo}
docs_dir=${base_dir}/docs
scripts_dir=${base_dir}/scripts

# Location of LMOD modules
spider=/app/lmod/lmod/libexec/spider
module_dir=/app/modules

echo Collecting Inventory
cd $base_dir
moduleclass='ai bio chem math'
module_search_path=''
for class in $moduleclass; do
   if [ -z $module_search_path ]; then
       module_search_path=${module_dir}/${class}
   else
       module_search_path=${module_search_path}:${module_dir}/${class}
   fi
done
echo $module_search_path
json_in=${docs_dir}/${label}bio-modules-${VERSION_ID}.json
$spider -o spider-json $module_search_path | python3 -mjson.tool > ${json_in} 

echo Generating Markdown
md_file=${label}bio-modules-${VERSION_ID}
md_out=${docs_dir}/${md_file}.md

echo '---' > ${md_out}
echo 'title: Bio Modules' $VERSION_ID >> ${md_out}
echo 'layout: single' >> ${md_out}
echo "permalink: /${md_file}/" >> ${md_out}
echo 'created: '`date +"%Y-%m-%d"` >> ${md_out}
echo 'toc: true' >> ${md_out}
echo 'toc_label: "On This Page"' >> ${md_out}
echo 'sidebar:' >> ${md_out}
echo '  nav: "docs"' >> ${md_out}
echo '---' >> ${md_out}
echo '' >> ${md_out}

# convert json modules spider to Markdown
cat ${json_in} | ${scripts_dir}/spider2post.py >> ${md_out}

echo Wrote inventory to ${md_out}

