#!/usr/bin/env python

import re
import os
import sys
import argparse
import imp
import json
import requests
from pprint import pprint
from pprint import pformat

"""
update to use pypi.org JSON API
  Project API:  GET /pypi/<project_name>/json
  Release API: GET /pypi/<project_name>/<version>/json
"""

__version__ = '1.3.0'
__author__ = 'John Dey jfdey@fredhutch.org'
__date__ = 'Jul 11, 2018'

class ExtsList(object):
    """ Extension List Update is a utilty program for maintaining EasyBuild
    easyconfig files for R, Python and Bioconductor.  R, Python and
    Bioconductor langueges support package extensions.  Easyconfig supports
    the building of the language with a list of extions. This program
    automates the the updating of extension lists for R and Python by using
    API for resolving current version for each package.

    command line arguments -> arguments
    --add [file] ->add_packages. add_packages is a list of package names to
          be added to exts_list.
    --check [package name] -> package.  If package is not None just check
            this one package and exit.

    Note: At the FredHutch BioCondutor packages are merged with R.
    An earlier version of easy_update built BioCondutor in a seperate package.

    Issues
       There are many small inconsistancies with PyPI which make it difficult
       to fully automate building of easyconfig files.
       - dependancy checking - check for extras=='all'
       - pypi projects names do not always match module names and or file names
         project: liac-arff, module: arff,  file name: liac_arff.zip
    """
    def __init__(self, args):  # Old:  file_name, add_packages, package, verbose
        self.verbose = args.verbose
        self.debug = False
        self.meta = args.meta
        self.code = None
        self.biocver = None 
        self.ext_counter = 0
        self.pkg_update = 0
        self.pkg_new = 0
        self.pkg_duplicate = 0

        self.ext_list_len = 1 
        self.exts_processed = []  # single list of package names
        self.depend_exclude = []  # built in packages not be added to exts_list
        self.bioc_data = {}
        self.prolog = '## remove ##\n'
        self.ptr_head = 0
        self.indent_n = 4
        self.indent = ' ' * self.indent_n

        # update easyconfig exts_list or check single package
        if args.easyconfig:
            eb = self.parse_eb(args.easyconfig, primary=True)
            self.search_pkg = None 
            self.exts_orig = eb.exts_list
            self.toolchain = eb.toolchain
            self.version = eb.version
            self.name = eb.name
            self.eb_filename = eb.name + '-' + eb.version
            self.dependencies = eb.dependencies
            self.eb_filename += '-' + eb.toolchain['name']
            self.eb_filename += '-' + eb.toolchain['version']
            try:
                self.eb_filename += eb.versionsuffix
            except AttributeError, NameError:
                print('versionsuffix not defined')
            try: 
                self.biocver = eb.biocver
                if self.debug: print('biocver: %s' % self.biocver)
            except:
                pass
            self.pkg_version = eb.version
            self.check_eb_package_name(args.easyconfig)
            if args.add_pkg:
                self.get_package_list(args.add_pkg) 
            self.out = open(self.eb_filename + ".update", 'w')
        elif args.search_pkg:
            self.search_pkg = args.search_pkg
            if args.biocver: self.biocver = args.biocver
            if args.pyver:
                self.name = "Python"
                self.version = args.pyver
            elif args.rver:
                self.name = "R"
                self.version = args.rver
            else:
                print('Languange and version must be specified with ' +
                      '[--pyver x.x | --rver x.x | --biocver x.x]')
            sea_pkg = {'name': args.search_pkg, 'version': 'x', 'type': 'add'}
            self.search_ext = sea_pkg


    def parse_eb(self, file_name, primary):
        """ interpret easyconfig file with 'exec'.  Interperting fails if
        constants that are not defined within the easyconfig file.  Add
        undefined constants to <header>.
        """
        header = 'SOURCE_TGZ  = "%(name)s-%(version)s.tgz"\n'
        header += 'SOURCE_TAR_GZ = "%(name)s-%(version)s.tar.gz"\n'
        header += self.prolog
        code = header

        eb = imp.new_module("easyconfig")
        with open(file_name, "r") as f:
            code += f.read()
        try:
            exec (code, eb.__dict__)
        except Exception as err:
            print("interperting easyconfig error: %s" % err)
            sys.exit(1)
        if primary:     # save original text of source code
            self.code = code
            self.ptr_head = len(header)
        return eb

    def get_package_list(self, fname):
        """read package names from <fname>
        add packages to exts_list to be checked. 
        """
        with open(fname, "r") as pkg_file:
            for pkg in pkg_file:
                d = dict()
                d = {'type': 'add'}
                self.exts_orig.append((pkg, 'x', d))


    def check_eb_package_name(self, easyconfig):
        """" check that easybuild filename matches package name
        easyconfig is filename of easyconfig file
        """
        f_name = os.path.basename(easyconfig)[:-3]
        if f_name != self.eb_filename:
            sys.stderr.write("Warning: file name does not match easybuild " +
                             "module name\n"),
            sys.stderr.write(" file name: %s, module name: %s\n" % (
                             f_name, self.eb_filename))
            sys.stderr.write('Writing output to: %s' % self.eb_filename +
                             '.update\n')


    def processed(self, pkg):
        """ save Processed packages
        save a normalize version of packae name to <exts_search> for Python
        updated July 2018
        """
        pkg['processed'] = True 
        self.exts_processed.append(pkg)


    def is_processed(self, pkg):
        """ check if package has been previously processed 
            if package exists AND is in the original exts_lists
                Mark as 'duplicate' 
        updated July 2018
        """
        name = pkg['name']
        for p_pkg in self.exts_processed:
            if 'modulename' in p_pkg['spec']:
                modulename = p_pkg['spec']['modulename']
            else:
                modulename = ''
            if name == p_pkg['name'] or name == modulename: 
                if pkg['type'] == 'orig':
                    pkg['action'] = 'duplicate'
                    self.pkg_duplicate += 1
                    self.processed(pkg)
                    if self.verbose: self.print_status(pkg)
                return True
            else:
                return False


    def get_package_info(self, pkg):
        """base class 
        subclasses of <get_package_info> search CRAN, Pypi or Bioconductor
        based on easyconfig.
        input: pkg - dict['name', 'version', 'spec': dict]
        return: msg, dict 
        <depends> list of package names that are dependancies
        look hard for dependant packages, try to match different cases
        and interchange underscore and dash.  If input package_name
        does not match rewrite pkg[0] with name found from repository.
        """
        pass


    def print_status(self, pkg):
        """ print one line status for each package if --verbose
        updated July 2018
        """
        if pkg['action'] == 'update':
            version = '%s -> %s' % (pkg['orig_ver'], pkg['version'])
        else:
            version = pkg['version']
        action = '(%s)' % pkg['action']
        tmpl = "%20s : %-20s %12s [%2d, %d]"
        print(tmpl % (pkg['name'], version, action,
                      self.ext_list_len, self.ext_counter))

    def print_meta(self, info):
        """ print meta data from repository
        this is broken for R 
        """
        pass
       
    def check_package(self, pkg):
        """query package authority [Pypi, CRAN, Bio] to get the latest version
        information for a package. This is heart of the program.

        input: pkg - list [pkg_name, pkg_version, dict]
        check that all dependancies are meet for each package.
        check_package can be called recursivly.
        dict['type'] is used to track status.
          - 'orig' is also used to track recursion
          - 'dep' package that is added as result of dependancy
          - 'add' packages read from file
        dict['action'] What action will be take to exts_list.
          - 'add'; new package
          - 'keep'; no update required
          - 'update'; version change
          - 'duplicate' package appears twice
        """
        if self.debug: print('check_package: %s' % pkg['name'])
        if self.is_processed(pkg):
            return
        msg, info = self.get_package_info(pkg)
        #msg, result, depends = self.get_package_info(pkg)
        if msg in ["error", 'not found']:
            if pkg['type'] == 'orig':
                pkg['action'] = 'keep'
                self.processed(pkg)
                return
            else:
                msg = "Warning: %s is dependency, but can't be found!"
                print(msg % pkg['name'])
                return

        if pkg['version'] == info['version']:
            pkg['action'] = 'keep'
        else:
            pkg['orig_ver'] = pkg['version']
            pkg['version'] = info['version']
            pkg['checksums'] = info['digests']['sha256']
            self.pkg_update += 1
            if pkg['type'] == 'orig': 
                pkg['action'] = 'update'
            elif pkg['type'] in ['dep', 'add']:
                if self.debug: print('check_package; dep or add')
                pkg['action'] = 'add'
                if self.name == "Python":
                    pkg['source_urls'] = "['PYPI_SOURCE']" 
                    self.pkg_new += 1

        if 'requests_dist' in info:
            for depend in info['requests_dist']:
                if depend not in self.depend_exclude:
                    depPkg = dict()
                    depPkg = {'name': depnd, 'type': 'dep'}
                    self.check_package(depPkg)
        self.processed(pkg)
        self.ext_counter += 1
        if self.search_pkg:
            output = self.output_module(pkg)
            print(output)
        if self.verbose:
            self.print_status(pkg)
        if self.meta:
            self.print_meta(info)
          

    def update_exts(self):
        """Loop through exts_list and check which packages need to be updated.
        this is an external method for the class
        """
        if self.search_pkg:
            self.check_package(self.search_ext)
        else:
            self.ext_list_len = len(self.exts_orig)
            for ext in self.exts_orig:
                if isinstance(ext, tuple):
                    pkg = {'name': ext[0],
                           'version': ext[1],
                           'spec': ext[2],
                           'type': 'orig'}
                    self.check_package(pkg)
                else:
                    self.processed({'name': ext})

    def write_chunk(self, indx):
        self.out.write(self.code[self.ptr_head:indx])
        self.ptr_head = indx

    def rewrite_extension(self, pkg):
        name = pkg['name']
        name_indx = self.code[self.ptr_head:].find(name)
        name_indx += self.ptr_head + len(name) + 1
        indx = self.code[name_indx:].find("'") + name_indx + 1
        self.write_chunk(indx)
        self.out.write("%s'," % pkg['version'])  # write version Number
        self.ptr_head = self.code[self.ptr_head:].find(',') + (
                        self.ptr_head + 1)
        if 'checksums' in pkg['spec']:
            indx = self.code[self.ptr_head:].find('checksums') + self.ptr_head + 10 
            self.write_chunk(indx)
            indx = self.code[self.ptr_head:].find("'") + self.ptr_head + 1 
            self.write_chunk(indx)
            self.out.write(pkg['checksums'])
            self.ptr_head = self.code[self.ptr_head:].find("'") + (
                        self.ptr_head)
        indx = self.code[self.ptr_head:].find('),') + self.ptr_head + 3
        self.write_chunk(indx)


    def output_module(self, pkg):
        """
        """
        
    def print_update(self):
        """ this needs to be re-written in a Pythonesque manor
        """
        indx = self.code.find('exts_list')
        indx += self.code[indx:].find('[')
        indx += self.code[indx:].find('\n') + 1
        self.write_chunk(indx)

        for extension in self.exts_processed:
            name = extension['name']
            if 'orig_ver' not in extension:  # no change skip
                indx = self.code[self.ptr_head:].find(name)
                indx += self.ptr_head + len(name) + 2
                self.write_chunk(indx)
                continue
            if extension['action'] in ['keep', 'update']:
                self.rewrite_extension(extension)
                # sys.exit(0)
            elif extension['action'] == 'duplicate':
                name_indx = self.code[self.ptr_head:].find(name)
                name_indx += self.ptr_head + len(name)
                indx = self.code[name_indx:].find('),') + name_indx + 3
                self.ptr_head = indx
                continue
            elif extensio['action'] in ['add','dep']:
                output = self.output_module(extension)
                self.out.write("%s\n" % output) 
        self.out.write(self.code[self.ptr_head:])
        print("Updated Packages: %d" % self.pkg_update)
        print("New Packages: %d" % self.pkg_new)
        print("Dropped Packages: %d" % self.pkg_duplicate)


    def download_url(self, filename, url):
        print('downloading: %s' % filename)
        #handle = open(self.source_dir + filename, "w")
        response = requests.get(url, stream=True)
        if response.status_code < 200 or response.status_code > 299:
            sys.stderr.write('download failed: %s filename: %s' % (
                resp.status_code, filename))
            return
        with open(self.source_dir + filename, "wb") as handle:
            for data in tqdm(response.iter_content()):
                handle.write(data)
        handle.close()
        

    def unpack_package(self, filename):
        content = None
        print('unpack: %s' % self.source_dir + filename)
        tar = tarfile.open(self.source_dir + filename, "r:gz")
        for member in tar.getmembers():
            f = tar.extractfile(member)
            print('extract: %s' % f)
            if f is not None and f == 'setup.py':
                content = f.read()
                break
        if content and 'install_requires' in content:
            print(content)
        


class R(ExtsList):
    """extend ExtsList class to update package names from CRAN and BioCondutor
    """
    def __init__(self, args):
        ExtsList.__init__(self, args)
        self.depend_exclude = ['R', 'base', 'compiler', 'datasets', 'graphics',
                               'grDevices', 'grid', 'methods', 'parallel', 
                               'splines', 'stats', 'stats4', 'tcltk', 'tools',
                               'utils',] 
        if self.biocver:
            self.read_bioconductor_pacakges()
        else:
            print('BioCondutor verserion: biocver not set')

    def read_bioconductor_pacakges(self):
        """ read the Bioconductor package list into bio_data dict
        """
        base_url = 'https://bioconductor.org/packages/json/%s' % self.biocver
        bioc_urls = ['%s/bioc/packages.json' % base_url,
                     '%s/data/annotation/packages.json' % base_url,
                     '%s/data/experiment/packages.json' % base_url]
        for url in bioc_urls:
            resp = requests.get(url)
            if resp.status_code != 200: 
                print('Error: %s %s' % (resp.status_code, url))
                sys.exit(1)
            self.bioc_data.update(resp.json())
            if self.debug:
                print('reading Bioconductor Package inf: %s' % url)
                pkgcount = len(self.bioc_data.keys())
                print('size: %s' % pkgcount)
             

    def get_CRAN_info(self, pkg):
        cran_list = "http://crandb.r-pkg.org/"
        resp = requests.get(url=cran_list + pkg['name'])

        cran_info = json.loads(resp.text)
        if 'error' in cran_info and cran_info['error'] == 'not_found':
            return "not found", []
        try:
            pkg_ver = cran_info[u'Version']
        except KeyError:
            return "error", []
        depends = []
        if u'License' in cran_info and u'Part of R' in cran_info[u'License']:
            return 'base package', []
        if u"Depends" in cran_info:
            depends.extend(cran_info[u"Depends"].keys())
        if u"Imports" in cran_info:
            depends.extend(cran_info[u"Imports"].keys())
        return pkg_ver, depends


    def get_BioC_info(self, pkg):
        """Extract <Depends> and <Imports> from BioCondutor json metadata
        Example:
        bioc_data['pkg']['Depends']
                    [u'R (>= 2.10)', u'BiocGenerics (>= 0.3.2)', u'utils']
        interesting fields from BioCoductor:
        bioc_data['pkg']['Depends', 'Imports', 'Biobase', 'graphics', 'URL']
        """
        depends = []
        if pkg['name'] in self.bioc_data:
            pkg_ver = self.bioc_data[pkg['name']]['Version']
            source_file = self.bioc_data[pkg['name']]['source.ver']
            if 'Depends' in self.bioc_data[pkg['name']]:
                depends.extend([re.split('[ (><=,]', s)[0]
                           for s in self.bioc_data[pkg['name']]['Depends']])
            if 'Imports' in self.bioc_data[pkg['name']]:
                depends.extend([re.split('[ (><=,]', s)[0]
                           for s in self.bioc_data[pkg['name']]['Imports']])
        else:
            pkg_ver = "not found"
        if self.debug: print('%s: Depends on: %s' % (pkg['name'], depends))
        return pkg_ver, depends

    def print_depends(self, pkg, depends):
        """ used for debugging """
        for p in depends:
            if p not in self.depend_exclude:
                print("%20s : requires %s" % (pkg, p))

    def get_package_info(self, pkg):
        """R version, check CRAN and BioConductor for version information
        """
        depends = []
        pkg_name = pkg['name']
        print('get_package_info: %s' % pkg_name)
        pkg_ver, depends = self.get_BioC_info(pkg)
        if pkg_ver == 'not found':
            pkg_ver, depends = self.get_CRAN_info(pkg)
            if pkg_ver == 'not found':
                return 'not found', [], []
            else:
                pkg['R_source'] = 'ext_options'
        else:
            pkg['R_source'] = 'bioconductor_options'
        return 'ok', (pkg_name, pkg_ver), depends


    def output_module(self, pkg):
        """R version: format a pkg for output"""
        output = "%s('%s', '%s', " % (self.indent, pkg['name'], pkg['Version'])
        for source in pkg[2].keys():
           if source == 'R_source': 
               output += pkg['R_source'] 
           elif source == 'ext_options':
               output += 'ext_options),'
           elif source == 'bioconductor_options':
               output += 'bioconductor_options),'
        output += '),'
        return output


class PythonExts(ExtsList):
    """extend ExtsList class to update package names from PyPI
    """
    def __init__(self, args):
        ExtsList.__init__(self, args)
        self.pkg_dict = None
        (nums) = self.version.split('.')
        self.python_version = "%s.%s" % (nums[0], nums[1])
        # Python >3.3 has additional built in modules
        if nums[0] == 3 and nums[1] > 3:
            self.depend_exclude = ['argparse', 'asyncio', ]
        if self.debug and self.search_pkg:
            print('Python Search PyPi: %s' % self.search_pkg)

    def parse_pypi_requires(self, eb_filename, requires):
        """pip requirement specifier is defined in full in PEP 508
        The project name is the only required portion of a requirement string.

        Only install the latest version so ignore all version information
        input: 'numpy (>=1.7.1)'  output: 'numpy'

        Test that <python_version> and <sys_platform> conform.
        If <extra> is present and required check that extra is contained
        in "exts_list".
        wincertstore (==0.2); sys_platform=='win32' and extra == 'ssl'
        futures (>=3.0); (python_version=='2.7' or python_version=='2.6')
        requests-kerberos (>=0.6); extra == 'kerberos'
        trollius; python_version == "2.7" and extra == 'asyncio'
        asyncio; python_version == "3.3" and extra == 'asyncio'
        """
        sys_platform = 'Linux'
        python_version = self.python_version
        extra = ''
        platform_python_implementation = 'CPython'
        require_re = '^([A-Za-z0-9_\-\.]+)(?:.*)$'
        extra_re = "and\sextra\s==\s'([A-Za-z0-9_\-\.]+)'"  # only if the
        targets = ['python_version', 'sys_platform', 'extra']
        ans = re.search(require_re, requires)
        name = ans.group(1)
        state = True    # result of eval(requires)

        version = requires.split(';')
        if len(version) > 1:
            for target in targets:
                if target in version[1]:
                    # eval Example: python_version == "3.3" and sys_platform == 'Linux'
                    try:
                        state = eval(version[1])
                    except NameError:
                        print('Error: NameError on eval, ignoring')
                    except SyntaxError:
                        pass
                    if target == 'extra':
                        extra = re.search(extra_re, version[1])
                        extra = extra.group(1) if extra else None
                        if extra not in [i[0] for i in self.exts_processed]:
                            extra = None
        if state:
            if self.verbose > 1:
                if name not in [i[0] for i in self.exts_processed] and (
                   name not in self.depend_exclude):
                    print('Add dependent package: %s ' % name),
                    print('for: %s, Expression: %s' % (pkg_name, requires))
            return name
        else:
            if self.verbose > 1:
                print('Do not install: %s, ' % name +
                      'for package: %s, Expression: %s' % (pkg_name, requires))
            return None

    def print_meta(self, info):
        """ Display meta from pypi.org
        """
        tags = ['filename', 'packagetype', 'url', 'python_version',
                'requires_dist', 'requires_python']
        for tag in tags:
            if tag in info:
                print("%s'%s': '%s'" % (self.indent, tag, info[tag]))

    def get_package_depends(self, filenam, url):
        """ inspect setup.py for requires_dist and install_requires"""
        if not os.path.isfile(self.source_dir + filename):
            self.download_url(filename, url)
        self.unpack_package(filename)
        # do more stuff

    def get_package_info(self, pkg):
        """Python version
           Get package version and dependency list via PyPi API.
           pkg is a list; ['package name', 'version', dict]
           return the version number for the package and a list of dependencies

           TODO; whl file should check versions. Output for whl should also
           if version == '3.6':
               pyver = 'cp36'
               arch = 'linux' ['manylinux', 'anylinux', 'linux']
           if arch == url_info['python_version']:
               use whl
        """
        (msg, info) = self.get_pypi_info(pkg)
        # (pkg_ver, filename, URL, digest)  = self.get_pypi_info(pkg[0])
        if msg == 'not found':
            return 'not found', {} 
        else:
            return 'ok', info 


    def print_python_exts_list(self, pkg, ver, source_tmpl, url_info):
        """ output <pkg> information in easyconfig exts_list format
            Used when <self.checkpage> is True; for checking a single
            package; Used with flag --pyver
        """
        indent = self.indent
        url = "https://pypi.python.org/packages/source/"
        url += "%s/%s" % (pkg[0], pkg)
        print("%s('%s', '%s', {" % (indent, pkg, ver))
        print("%s%s'source_urls': ['%s']," % (indent, indent, url))
        if source_tmpl:
            print("%s%s'source_tmpl': '%s'," % (indent, indent, source_tmpl))
        print("%s})," % indent )


    def get_pypi_release(self, pkg, version):
        """ if source dist is not available from pypi search
        the release for for a whl file.
        """
        req = 'https://pypi.org/pypi/%s/%s/json' % (pkg['name'], version)
        resp = requests.get(req)
        if resp.status_code != 200:
            msg = "API error: %s GET release %s"
            sys.stderr.write(msg % (resp.status_code, pkg['name'])) 
            return 'not found', {}
        release = resp.json()
        cplist = ['cp35', 'cp36', 'cp37']
        for rel in release['releases'][version]:
            if any(cver in rel['python_version'] for cver in cplist):
                if 'manylinux' in rel['filename']:
                    return 'ok', rel 
        return 'not found', {}


    def get_pypi_info(self, pkg):
        """get version information from pypi.  If <pkg_name> is not processed
           seach pypi. pkg_name is now case sensitive and must match  
           info['digests']['sha256'], 'summary', 'url', 'filename', 'home_page'
        """
        status = 'not found'
        resp = requests.get('https://pypi.org/pypi/%s/json' % pkg['name'])
        if resp.status_code != 200:
                msg = "API error: %s GET %s" % (resp.status_code, pkg['name'])
                sys.stderr.write(msg)
                return 'not found', {}
        project = resp.json()
        info = project['info']
        version = info['version']
        for ver in project['releases'][version]:
            if 'packagetype' in ver and ver['packagetype'] == 'sdist':
                info.update(ver)
                status = 'ok'
                break
        if status != 'ok':
             status, rel = self.get_pypi_release(pkg, version)
             if status == 'ok':
                 info.update(rel)
        return status, info 


    def output_module(self, pkg):
        """Python version: format single pkg for output. 
        Used if --search argument is used.
        if self.search_pkg:
        """
        output = "%s('%s', '%s', {\n" % (self.indent, pkg['name'], pkg['version'])
        for item in pkg.keys():
           if item in ['name', 'version', 'action', 'type', 'orig_ver',
                       'processed']:
               continue 
           if item == 'checksums':
               output += "%s%s'%s': ['%s'],\n" % (self.indent, self.indent, item, pkg[item])
           else:
               output += "%s%s'%s': %s,\n" % (self.indent, self.indent, item, pkg[item])
        output += "%s})," % self.indent
        return output


def help():
    print("usage: easy_update  easyconfig.eb [flags]")
    print("easy_update Updates ext_list information of EasyBuild"),
    print(" easyconfig files")
    print("easy_update works with R, Python and R-bioconductor"),
    print(" easyconfig files")
    print("  --verbose  diplay status for each package")
    print("  --add [filename]  filename contains list of package"),
    print(" names to add")
    sys.exit() 


def main():
    """ main """
    parser = argparse.ArgumentParser(description='Update easyconfig extslist')
    
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + __version__)
    parser.add_argument('-v', '--verbose', dest='verbose', required=False,
                        action='store_true',
                        help='Verbose; print lots of extra stuff, (default: false)')
    parser.add_argument('--rver', dest='rver', required=False, action='store',
                        help='Set R version (major.minor) example 3.4')
    bioc_help = 'Set BioConductor version (major.minor) example 3.6. Use with --rver'
    parser.add_argument('--biocver', dest='biocver', required=False, action='store',
                        help=bioc_help)
    parser.add_argument('--pyver', dest='pyver', required=False, action='store',
                        help='Set Python version [2.7 or 3.6]')
    parser.add_argument('--add', dest='add_pkg', required=False, action='store',
                        help='File that contains additional packages to be added')
    sea_help = 'Search for single package. requires --rver or --pyber' 
    parser.add_argument('--search', dest='search_pkg', required=False, action='store',
                        help=sea_help)
    parser.add_argument('--meta', dest='meta', required=False, action='store_true',
                        help='output all meta data keys from Pypi, (default: false)')
    parser.add_argument('easyconfig', nargs='?') 
    args = parser.parse_args()
  
    eb_name = '' 
    if args.easyconfig:
        eb_name = os.path.basename(args.easyconfig)
    elif args.search_pkg:
        pass
    else:
        print('If no easyconfig is given, a module name must be ' +
              'specified with --search pkg_name') 
        sys.exit()


    if args.rver or eb_name[:2] == 'R-':
        module = R(args)
    elif args.pyver or eb_name[:7] == 'Python-':
        module = PythonExts(args)
    else:
        msg = "EasyConfig file name must begin with R- or Python- OR"
        msg += "  --search argument must be used"
        print(msg)
        sys.exit(1)
    module.update_exts()
    if args.easyconfig:
        module.print_update()

if __name__ == '__main__':
    main()
