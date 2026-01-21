#!/usr/bin/env python3

"""
file: build_list.py

  Compare two LMOD systems to create a list of
  software that needs to be built.

  - The comparison only works one way, from old to new
    creating a list of what is missing on the <new> LMOD
    system when comparing to the <old>.
  - Comparison is done in two steps, since you cant' run
    two LMOD systems on one host.
  - Step One - Run create_module_list.py on the <old> LMOD system.
    This will generate JSON file from Spider of all LMOD modules.
  - Step Two - from <new> system run this script with the JSON
    file as the only argument.

"""

from toolchain import Toolchain
import json
import sys
import os
from os import listdir
from os.path import isfile, join
import pathlib

""" <cutoff> is the lowest matching tool chain to be used """
cutoff = '2022b'

local_eb_path = None
eb_path = None

def find_easyconfig_paths():
    """find the paths to EasyConfigs, search within the path given by the easyconfig
    given to update.
    This script exists ../scripts and asumes the `easyconfigs` directory is one level up
    """
    global local_eb_path, eb_path

    script_path = pathlib.Path(__file__).parent.resolve()
    (head, tail) = os.path.split(script_path)
    local_eb_path = os.path.join(head, 'easyconfigs')
    if not os.path.exists(local_eb_path):
        print('You are not working in an EB repository, Quiting because I can not find easyconfig directory.')
        sys.exit(1)
    eb_root = os.getenv('EBROOTEASYBUILD')
    if eb_root is None:
        print('$EBROOTEASYBUILD environment variable must be defined.')
        sys.exit(1)
    else:
        eb_path = os.path.join(eb_root, 'easybuild/easyconfigs')
    print(f"eb path: {eb_path}")

def find_easyconfig(package_name, easyconfig):
    """ search base_paths for easyconfig filename """
    found = None
    first_letter = package_name[0].lower()
    for base_path in [eb_path, local_eb_path]:
        narrow = os.path.join(base_path, first_letter, package_name)
        if not os.path.isdir(narrow):
            if base_path == eb_path:
                print(f"# package is not upstream: {package_name}")
            continue
        for file in listdir(narrow):
            if file == easyconfig:
               if base_path == eb_path:
                   return
               else:
                  full_ec_path = os.path.join(narrow, file)
                  print(f"# staging {easyconfig}")
                  print(f"cp {full_ec_path} /loc/eb_build")
                  return

if len(sys.argv) > 1:
    with open(sys.argv[1]) as f:
        data = json.load(f)
else:
    print(f'usage: build_list json_file')
    sys.exit(1)

find_easyconfig_paths()
tc_tools = Toolchain(cutoff)

"""
  JSON output from LMOD spider has full pathname to lua module
"""
for package_name in data.keys():
    for name in data[package_name]:
        if tc_tools.cutoff(data[package_name][name]['fullName']):
            if not os.path.isfile(name):
                eb = data[package_name][name]['fullName'].replace('/', '-') + ".eb"
                find_easyconfig(package_name, eb)
                print(f"echo ' == install {package_name}'")
                print(f"eb {eb} --robot")
