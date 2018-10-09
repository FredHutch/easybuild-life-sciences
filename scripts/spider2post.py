#!/usr/bin/env python

import json
import platform

(name, version, id) = platform.dist()

jfile = open('../docs/modules-%s.json' % version, 'r')
outfile = open('../docs/bio-modules-%s.md' % version, 'a')

data = json.load(jfile)
packages = data.keys()
slist = sorted(packages)
for p in slist:
   pac = data[p]
   for ver in pac.keys():
       if 'fullName' in pac[ver] and '-2015' in pac[ver]['fullName']:
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
       outfile.write(' - [' + pac[ver]['fullName'] + ']')
       outfile.write('(' + url + ')  ')
       outfile.write(descrp + '\n')
