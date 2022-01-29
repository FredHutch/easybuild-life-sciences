#!/usr/bin/env python3

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
from packaging import version

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
   print(p, file=sys.stderr)
