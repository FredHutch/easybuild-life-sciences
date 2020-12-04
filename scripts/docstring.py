#!/usr/bin/env python
import sys
import ast


"""some handy ast code to parse static Python and extract elements
   print to stdout the value of 'description' from an easyconfig
"""
__date__ = 'Dec 2020'
__author__ = 'John Dey jfdey@fredhutch.org'
__version__ = '0.1'

if len(sys.argv) != 2:
   sys.stderr.write('argument needs to be easyconfig file')
   sys.ext(1)

with open(sys.argv[1]) as f:
    code = ast.parse(f.read())

for node in ast.walk(code):
    if isinstance(node, (ast.Assign)):
        if len(node.targets) == 1:
            name = node.targets[0]
            if isinstance(name, ast.Name):
                #print('Name: {}'.format(name.id))
                v = node.value
                if name.id == 'description' and isinstance(v, ast.Str):
                    print(v.s)
                #elif isinstance(v, ast.List):
                #    for e in v.elts:
                #       if isinstance(e, ast.Str):
                #           print('  list val: {}'.format(e.s))
