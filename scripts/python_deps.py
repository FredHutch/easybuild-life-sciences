#!/usr/bin/env python

import re
import os
import sys
import imp
import json
import urllib2
import xmlrpclib

def parse_pypi_requires(pkg_name, requires):
    """ pip requirement specifier is defined in full in PEP 508
        The project name is the only required portion of a requirement string. 

        Test that <python_version> and <sys_platform> conform.
        If <extra> is present and required check that extra is contained in "exts_list".
        
        We only care about the package name so ignore all version information
        input: 'numpy (>=1.7.1)'  output: ['numpy', '']
        badness:
            psutil (<2.0.0,>=1.1.1)
            netaddr (!=0.7.16,>=0.7.13)
            requests>=2.0    # no delimiters at all
            apscheduler  != APScheduler
    """
    python_version = '2.7'
    sys_platform = 'Linux'
    extra = ''
    require_re = '^([A-Za-z0-9_\-\.]+)(?:.*)$'
    extra_re =   "and\sextra\s==\s'([A-Za-z0-9_\-\.]+)'"  # only if the 
    targets = ['python_version', 'sys_platform', 'extra'] 
    ans = re.search(require_re, requires)
    name = ans.group(1)
    test = False    # do we need to test extra requires field?
    state = True

    version = requires.split(';')
    if len(version) > 1:
        for target in targets:
            if target in version[1]:
                 test = True
                 if target == 'extra':
                     extra = re.search(extra_re, version[1])
                     extra = extra.group(1) if extra else None
        if test:
            state = eval(version[1])
    if state:
        print('Install!    depend: package: %s, eval: %r, Expression: %s' % (
                name, state, requires) )
        return name
    else:
        print('No install: depend: package: %s, eval: %r, Expression: %s' % (
                name, state, requires) )
        return None

def get_package_info(client, pkg_name):
    """Python pypi API for package version and dependancy list 
       pkg is a list; ['package name', 'version', 'other stuff']
       return the version number for the package and a list of dependancie
    """
    indent4 = '    '
    depends = []
    xml_vers = client.package_releases(pkg_name)
    if xml_vers:
        pkg_ver = xml_vers[0]
        xml_info = client.release_data(pkg_name, pkg_ver)
        print("%s %s \n" % (pkg_name, pkg_ver)),
        url_info = client.release_urls(pkg_name, pkg_ver)
        found = False
        for url in url_info:
            if (url['url'].endswith('.gz') or 
                url['url'].endswith('.zip') or
                url['url'].endswith('whl')):
                URL = url['url']
                print('wget %s' % URL)
                found = True
            if url['filename'].endswith('.gz') or (
                url['filename'].endswith('.zip')):
                print('filename: %s' % url['filename'])
        if not found:
            print('Error url not found!')
        print('Keys: %s' % xml_info.keys() )
        if 'requires_dist' in xml_info.keys():
            for requires in xml_info['requires_dist']:
                pkg_requires = parse_pypi_requires(pkg_name, requires)
                if pkg_requires:
                    print("  required: %s" % pkg_requires) 

        print("\n%s('%s', '%s', {" % (indent4, pkg_name, pkg_ver))
        print("%s%s'source_url': [%s]," % (indent4, indent4, URL))
        print("%s)}," % indent4 )
    else:
        print("Warning: %s Not in PyPi. No depdancy checking performed" % pkg_name)
        



if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("usage: %s[Python module" % sys.argv[0])
        sys.exit(0)

    client = xmlrpclib.ServerProxy('https://pypi.python.org/pypi')
    get_package_info(client, sys.argv[1])
