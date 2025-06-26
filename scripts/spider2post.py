#!/usr/bin/env python3

""" convert output from module spider to a single Markdown page.
Not all modules can be located in the Fred Hutch Easybuild-life-scienes repo.
check for module type of bio, then search for easyconfig

Note: 'repo_path' needs to be relative to ../easyconfigs, to create
URL links to EasyCconfigs
"""

from pprint import pprint
import json
import sys
import os
from packaging import version

data = json.load(sys.stdin)
packages = data.keys()
slist = sorted(packages)
repo_path = os.path.dirname(os.path.abspath(__file__)) + '/../easyconfigs/'
github_repo = 'https://github.com/FredHutch/easybuild-life-sciences/'
github_repo += 'blob/master/easyconfigs/'
if not os.path.isdir(repo_path):
   print("don't know where I am")
   sys.exit(1)
for p in slist:
   pac = data[p]
   print(p, file=sys.stderr)
   paths = list(pac.keys())
   if len(paths) == 1:
      latest = pac[paths[0]]
   else:
      maxVal = version.parse("0.0.0")
      for release in pac.keys():
         verVal = version.parse(pac[release]['Version'].split('-')[0].replace('_', '.'))
         if verVal > maxVal:
             latest = pac[release]
             maxVal = verVal
             latestVersion = pac[release]['Version']
   if 'fullName' in latest and '-2015' in latest['fullName']:
        continue
   descrp = ''
   url = ''
   easyconfig_url = None
   fullName = latest['fullName']
   eb_filename = fullName.replace('/', '-') + '.eb'
   topdir = fullName[0].lower()
   projdir = p
   eb_path = repo_path +  topdir +'/'+ projdir +'/'+ eb_filename
   # search for eb_filename
   if os.path.isfile(eb_path):
       easyconfig_url = github_repo + topdir +'/'+ projdir +'/'+ eb_filename
   eb_filename = None
   if 'Description' in latest:
       text = latest['Description'].split(' - ')[0]
       if isinstance(text, bytes):
           descrp = text.decode()
       elif isinstance(text, str):
           descrp = text
   if 'whatis' in latest:
       entry = None
       entry = [x for x in latest['whatis'] if 'Homepage: ' in x]
       if entry and 'Homepage: ' in entry:
           text = entry[0].split('Homepage: ')[1]
       if isinstance(text, bytes):
           url = text.decode()
       elif isinstance(text, str):
           url = text
   print(' - [{}]({})'.format(latest['fullName'], url))
   if easyconfig_url:
       print('[easyconfig]({})'.format(easyconfig_url) )
   else:
       print
   print('{}'.format(descrp))
