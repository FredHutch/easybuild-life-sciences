#! /usr/bin/env python

import os, sys, ast

def _convert(node):
        if isinstance(node, ast.Str):
            return node.s
        elif isinstance(node, ast.Num):
            return node.n
        elif isinstance(node, ast.Tuple):
            return tuple(map(_convert, node.elts))
        elif isinstance(node, ast.List):
            return list(map(_convert, node.elts))
        elif isinstance(node, ast.Dict):
            return dict((_convert(k), _convert(v)) for k, v
                        in zip(node.keys, node.values))
        elif isinstance(node, ast.Name):
            return node.id
        else:
            return None

def _pkgadd (aptpkg, yumpkg, osdepts):
    for pkglist in osdepts:
        if isinstance(pkglist, basestring):
            pkglist = (pkglist,)        
        for pkg in pkglist:
            if pkg.endswith('-dev'):
                if not pkg in aptpkg:
                    aptpkg.append(pkg)
                    print ("added %s to os-dependencies.apt" % pkg)
            elif pkg.endswith('-devel'):
                if not pkg in yumpkg:
                    yumpkg.append(pkg)
                    print ("added %s to os-dependencies.yum" % pkg)
            else:
                if not pkg in aptpkg:
                    aptpkg.append(pkg)
                    print ("added %s to os-dependencies.apt" % pkg)
                if not pkg in yumpkg:
                    yumpkg.append(pkg)
                    print ("added %s to os-dependencies.yum" % pkg)
    return aptpkg, yumpkg

# ----------------------------------------------------------------------

ebfolder = 'easybuild/easyconfigs'
aptpkg = []
yumpkg = []

for root, dirs, files in os.walk(ebfolder):
    for f in files:
        if f.endswith(".eb"):
             #print(os.path.join(root, f))
            getnext=False
            p = os.path.join(root, f)
            t=ast.parse(open(p).read())
            for node in ast.walk(t):
                if getnext:
                    osdepts=(_convert(node))
                    aptpkg, yumpkg = _pkgadd(aptpkg, yumpkg, osdepts)
                    getnext=False
                if _convert(node) == 'osdependencies':
                    getnext=True

aptf = open('os-dependencies.apt','w')
yumf = open('os-dependencies.yum','w')
aptf.write(' '.join(aptpkg))
yumf.write(' '.join(yumpkg))

print ("os-dependencies.apt and os-dependencies.yum written.")
