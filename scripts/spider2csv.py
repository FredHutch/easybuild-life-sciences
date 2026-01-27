#!/usr/bin/env python3

""" convert JSON output from module spider to a single [markdown/csv] file.
Not all modules can be located in the Fred Hutch Easybuild-life-scienes repo.
check for module type of bio, then search for easyconfig

Note: 'repo_path' needs to be relative to ../easyconfigs, to create
URL links to EasyCconfigs
"""

from pprint import pprint
import json
import sys
import os
import re
from packaging import version
from tc import Toolchain as TC

__author__ = "John Dey"
__date__ = "Jan 2026"
__version__ = "1.0.0"

data_path = os.path.dirname(sys.argv[1])
file_name = os.path.basename(sys.argv[1])
base_name = file_name.removesuffix('.json')
markdown_name =os.path.join(data_path , base_name + ".md") 
csv_name = os.path.join(data_path , base_name + ".csv")
csv_f = open(csv_name, "w")
md_f = open(markdown_name, 'a')

with open(sys.argv[1], 'r') as file:
    json_data = file.read()
    data = json.loads(json_data)

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
   latestVersion = ""
   if len(paths) == 1:
      latest = pac[paths[0]]
      latestVersion = pac[paths[0]]['Version']
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
   software_version = TC.tc_trim(latestVersion)
   cleaned_descrp = "none"
   url = ''
   fullName = latest['fullName']
   if 'Description' in latest:
       text = latest['Description'].split(' - ')[0]
       if isinstance(text, bytes):
           descrp = text.decode()
       elif isinstance(text, str):
           descrp = text
       cleaned_descrp = re.sub(r'[\x00-\x1F\x7F]', ' ', descrp)  # Remove control characters
       if '\r' in cleaned_descrp:
           print(f"{p} has an issue")
   url = latest.get('URL', "none")
   print(f'"{p}","{software_version}","{latestVersion}","{url}","{cleaned_descrp}"', file=csv_f)
   print(f" - [{latest['fullName']}]({url})", file=md_f)
   print(f"{cleaned_descrp}", file=md_f)

csv_f.close()
md_f.close()

if __file__ == '__main__':
    if len(sys.argv) != 2:
       print('usage: need one filename argument')
       sys.exit(1)
