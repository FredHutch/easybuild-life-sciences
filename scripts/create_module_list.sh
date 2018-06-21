#!/bin/bash

# create_module_list.sh
#
# Create a list of modules. Filter modules based on module class 'bio', 'chem', 'phys' and
# 'math'

repo='easybuild-life-sciences'

list='bio
chem
phys
math'

if [[ ! -z "${PWD##*${repo}*}" ]]; then
    echo "Can not find github pages docs directory."
    echo "Run script from github repo: ${repo}"
    exit 1
fi

cwd=$PWD
# remove the $repo from left everything else; add the $repo and the "/docs" dir back
base_dir=${cwd%${repo}*}${repo}
docs_dir=${base_dir}/docs
scripts_dir=${base_dir}/scripts

for mtype in $list; do
    /app/Lmod/lmod/lmod/libexec/spider -o /app/easybuild/modules/${mtype} >> $docs_dir/mod_list.txt
done
spider=/app/Lmod/lmod/lmod/libexec/spider
module_dir=/app/easybuild/modules
$spider -o spider-json ${module_dir}/bio | python -mjson.tool >bio.json


