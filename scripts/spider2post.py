#!/usr/bin/env python

import json
import sys

data = json.load(sys.stdin)
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
       print(' - [' + pac[ver]['fullName'] + ']')
       print('(' + url + ')  ')
       print(descrp + '\n')
