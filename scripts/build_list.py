#!/usr/bin/env python3

"""
Create a list of software that needs to be built.

Filter by toolchain
input is output from spider *.json file
"""

from toolchain import toolchain
import json
import sys

data = json.load(sys.stdin)
packages = data.keys()
slist = sorted(packages)
tc = toolchain('2022b')

for p in slist:
    pac = data[p]
    paths = list(pac.keys())
    if len(paths) == 1:
       latest = pac[paths[0]]
    else:
       maxVal = None
       for release in pac:
          verVal = pac[release]['pV']
          if maxVal and verVal > maxVal:
              latest = pac[release]
              maxVal = verVal
          if not maxVal:
              latest = pac[release]
              maxVal = verVal
              latestVersion = pac[release]['Version']
    if tc.tc_ge(latest['Version']):
        print('{}:{}'.format(p,latest['Version']))
    else:
        print('Needs updated: {}:{}.eb'.format(p,latest['Version']), file=sys.stderr)

