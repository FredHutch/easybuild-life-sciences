#!/bin/bash

# create_module_list.sh
#
# Create a list of modules. Filter modules based on module class 'bio', 'chem', 'phys' and
# 'math'

repo='easybuild-life-sciences'

if [[ ! -z "${PWD##*${repo}*}" ]]; then
    echo "Can not find github pages docs directory."
    echo "Run script from github repo: ${repo}"
    exit 1
fi

# remove the $repo from left everything else; add the $repo and the "/docs" dir back
base_dir=${PWD%${repo}*}${repo}
docs_dir=${base_dir}/docs
scripts_dir=${base_dir}/scripts

spider=/app/Lmod/lmod/lmod/libexec/spider
module_dir=/app/easybuild/modules
cd $base_dir
$spider -o spider-json ${module_dir}/bio:${module_dir}/math | python -mjson.tool >${docs_dir}/modules.json

echo '---' > ${docs_dir}/bio-modules.md
echo 'layout: post' >> ${docs_dir}/bio-modules.md
echo 'title: Bio Modules' >> ${docs_dir}/bio-modules.md
echo 'date: '`date +'%Y-%m-%d'` >> ${docs_dir}/bio-modules.md
echo '---' >> ${docs_dir}/bio-modules.md
echo '' >> ${docs_dir}/bio-modules.md 

python - <<EOF
import json
jfile = open('docs/modules.json', 'r')
outfile = open('docs/bio-modules.md', 'a')
data = json.load(jfile)
packages = data.keys()
slist = sorted(packages)
for p in slist:
   pac = data[p]
   for ver in pac.keys():
       if '-2015' in pac[ver]['full']: 
            continue
       descrp = ''
       url = ''
       if 'Description' in pac[ver]:
           text = pac[ver]['Description'].split(' - ')[0]
           descrp = text.encode('utf8', 'replace')
       if 'whatis' in pac[ver]:
           entry = [x for x in pac[ver]['whatis'] if 'Homepage: ' in x]
           text = entry[0].split('Homepage: ')[1]
           url = text.encode('utf8', 'replace')
       outfile.write(' - [' + pac[ver]['full'] + ']')
       outfile.write('(' + url + ')  ')
       outfile.write(descrp + '\n')
EOF

