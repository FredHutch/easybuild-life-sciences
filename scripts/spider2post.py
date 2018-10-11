#!/usr/bin/env python

from pprint import pprint
import json
import sys

data = json.load(sys.stdin)
packages = data.keys()
slist = sorted(packages)
for p in slist:
   pac = data[p]
   paths = list(pac.keys())
   if len(paths) == 1:
      release = pac[paths[0]]
   else:
      vers = sorted(paths)
      release = pac[vers[-1]]
   if 'fullName' in release and '-2015' in release['fullName']:
        continue
   descrp = ''
   url = ''
   if 'Description' in release:
       text = release['Description'].split(' - ')[0]
       descrp = text.encode('utf8', 'replace')
   if 'whatis' in release:
       entry = [x for x in release['whatis'] if 'Homepage: ' in x]
       text = entry[0].split('Homepage: ')[1]
       url = text.encode('utf8', 'replace')
   print(' - [' + release['fullName'] + ']')
   print('(' + url + ')  ')
   print(descrp + '\n')
