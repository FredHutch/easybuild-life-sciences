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

# extract major.minor from OS relase: VERSION="16.04.3 LTS (Xenial Xerus)"
os_ver=`grep '^VERSION=' /etc/os-release | sed 's/.*="\([0-9]*\.[0-9]*\).*$/\1/'`
inventory=bio-modules-${os_ver}.md

# remove the $repo from left everything else; add the $repo and the "/docs" dir back
base_dir=${PWD%${repo}*}${repo}
docs_dir=${base_dir}/docs
scripts_dir=${base_dir}/scripts

if [[ $os_ver == '16.04' ]]; then
    spider=/app/lmod/lmod/libexec/spider
    module_dir=/app/modules
elif [[ $os_ver == '14.04' ]]; then
    spider=/app/Lmod/lmod/lmod/libexec/spider
    module_dir=/app/easybuild/modules
fi

echo Collecting Inventory
cd $base_dir
$spider -o spider-json ${module_dir}/bio:${module_dir}/math | python -mjson.tool >${docs_dir}/modules-${os_ver}.json

echo Generating Markdown
json_in=${docs_dir}/modules-${os_ver}.json
md_out=${docs_dir}/bio-modules-${os_ver}.md

echo '---' > ${md_out}
echo 'layout: post' >> ${md_out}
echo 'title: Bio Modules ' $os_ver >> ${md_out}
echo 'date: '`date +'%Y-%m-%d'` >> ${md_out}
echo '---' >> ${md_out}
echo '' >> ${md_out}

# convert json modules spider to Markdown
cat ${json_in} | ${scripts_dir}/spider2post.py >> ${md_out}

echo Wrote inventory to ${md_out}

