#!/bin/bash

# create_module_list.sh
#
# Create a list of modules. Filter modules based on module class 'bio', 'chem', 'phys' and
# 'math'
# Used to create software Inventory

repo='easybuild-life-sciences'

if [[ ! -z "${PWD##*${repo}*}" ]]; then
    echo "Can not find github pages docs directory."
    echo "Run script from github repo: ${repo}"
    exit 1
fi

label=""
if [[ $# -eq 1 ]]; then
    label="${1}-"
fi

# extract major.minor from OS relase: VERSION="16.04.3 LTS (Xenial Xerus)"
os_ver=`grep '^VERSION=' /etc/os-release | sed 's/.*="\([0-9]*\.[0-9]*\).*$/\1/'`
inventory=all-modules-${os_ver}.md

# remove the $repo from left everything else; add the $repo and the "/docs" dir back
base_dir=${PWD%${repo}*}${repo}
docs_dir=${base_dir}/docs
scripts_dir=${base_dir}/scripts

if [[ $os_ver == '14.04' ]]; then
    spider=/app/Lmod/lmod/lmod/libexec/spider
    module_dir=/app/easybuild/modules
else
    spider=/app/lmod/lmod/libexec/spider
    module_dir=/app/modules
fi

echo Collecting Inventory
cd $base_dir
json_in=${docs_dir}/${label}all-modules-${os_ver}.json
$spider -o spider-json ${module_dir}/all | python3 -mjson.tool >${json_in}

echo Generating Markdown
md_file=${label}all-modules-${os_ver}
md_out=${docs_dir}/${md_file}.md

echo '---' > ${md_out}
echo 'title: All Modules' $os_ver >> ${md_out}
echo "permalink: /${md_file}/" >> ${md_out}
echo 'layout: single' >> ${md_out}
echo 'toc: true' >> ${md_out}
echo 'toc_label: "On This Page"' >> ${md_out}
echo 'sidebar:' >> ${md_out}
echo '  nav: "docs"' >> ${md_out}
echo '---' >> ${md_out}
echo '' >> ${md_out}

# convert json modules spider to Markdown
cat ${json_in} | ${scripts_dir}/spider2post.py >> ${md_out}

echo Wrote inventory to ${md_out}

