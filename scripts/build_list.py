#!/usr/bin/env python3

"""
file: build_list.py
Create a list of software that needs to be built.

Build list is based on a existing system

input is output from spider *.json file
"""

from toolchain import Toolchain
import json
import sys
import os

""" <cutoff> is the lowest matching tool chain to be used """
cutoff = '2022b'

if len(sys.argv) > 1:
    with open(sys.argv[1]) as f:
        data = json.load(f)
else:
    print(f'usage: build_list json_file')
    sys.exit(1)

tc_tools = Toolchain(cutoff)

for p in data.keys():
    for name in data[p]:
        print(f"name: {name}")
        if tc_tools.cutoff(data[p][name]['fullName']):
            if os.path.isfile(name):
                print(f"Installed: {name}")
            else:
                eb = data[p][name]['fullName'].replace('/', '-') + ".eb"
                print(f"eb {eb} --robot")
        else:
            print(f"== cutoff: {data[p][name]['fullName']}")
