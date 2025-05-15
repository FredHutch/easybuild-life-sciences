#!/usr/bin/env python3

"""
file: update_easyconfig.py 

   modules is a file with a single module name one per line.

   Processes multiple output for modules.
   - is two versions with and without CUDA
   - two versions per toolchain
   - Python 2 and 3 versions

   todo, Pick the highest version
   add a --CUDA flag, no flag no CUDA 
"""

from toolchain import Toolchain
import json
import sys
import os

def read_project(project_file):
    modules = []
    with open(project_file, "r") as f:
        for line in f:
            modules.append(line.rstrip('\n'))
    return modules

def set_easyconfig_path():
    """
    get the easyconfig file path
    """
    ebroot = os.environ.get('EBROOTEASYBUILD')
    if not ebroot:
        print(f"environment EBROOTEASYBUILD must be set")
        sys.exit(1)
    return ebroot

def find_modules(modules, tc):
    global ebroot
    for mod  in modules:
        letter_dir = mod[0].lower()
        mod_dir = os.path.join(ebroot, 'easybuild/easyconfigs', letter_dir, mod)
        found = []
        for eb_config in os.listdir(mod_dir):
            if tc.tc_filter(eb_config):
                found.append(eb_config)
                print(f"using: {eb_config}")
        if len(found) == 0:
            print(f"missing toolchain for {mod}")

def main():
    global ebroot

    if len(sys.argv) < 3:
        print("Usage: update_toolchain.py <project_file> <toolchain> <eb_path>")
        sys.exit(1)

    project_file = sys.argv[1]
    toolchain = sys.argv[2]
    ebroot = set_easyconfig_path()
    tc_tools = Toolchain(toolchain)
    modules = read_project(project_file)
    find_modules(modules, tc_tools)

if __name__ == '__main__':
    main()
