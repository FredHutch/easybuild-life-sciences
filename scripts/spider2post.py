#!/usr/bin/env python

""" convert output from module spider to a single Markdown page.
Not all modules can be located in the Fred Hutch Easybuild-life-scienes repo.
check for module type of bio, then search for easyconfig

Note: 'repo_path' needs to be relative to ../fh_easyconfigs, to create
URL links to EasyCconfigs
"""

from pprint import pprint
import json
import sys
import os

data = json.load(sys.stdin)
packages = data.keys()
slist = sorted(packages)
repo_path = os.path.dirname(os.path.abspath(__file__)) + '/../fh_easyconfigs/'
github_repo = 'https://github.com/FredHutch/easybuild-life-sciences/'
github_repo += 'blob/master/fh_easyconfigs/'
if not os.path.isdir(repo_path):
   print("don't know where I am")
   sys.exit(1)
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
   easyconfig_url = None
   fullName = release['fullName']
   eb_filename = fullName.replace('/', '-') + '.eb'
   topdir = fullName[0].lower()
   projdir = p
   eb_path = repo_path +  topdir +'/'+ projdir +'/'+ eb_filename
   # search for eb_filename
   if os.path.isfile(eb_path):
       easyconfig_url = github_repo + topdir +'/'+ projdir +'/'+ eb_filename
   eb_filename = None
   if 'Description' in release:
       text = release['Description'].split(' - ')[0]
       descrp = text.encode('utf8', 'replace')
   if 'whatis' in release:
       entry = [x for x in release['whatis'] if 'Homepage: ' in x]
       text = entry[0].split('Homepage: ')[1]
       url = text.encode('utf8', 'replace')
   print(' - [' + release['fullName'] + '](' + url + ')  '),
   if easyconfig_url:
       print('[easyconfig](' + easyconfig_url + ')  ')
   else:
       print
   print(descrp)
