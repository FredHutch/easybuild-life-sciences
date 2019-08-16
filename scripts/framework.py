#!/usr/bin/env python

import re
import os
import sys
import argparse
import imp
import requests

""" 1.0.1 Aug, 15, 2019
    fix parse_dependencies
    For the case of reading Python dependancies, conver the
    case of 'Biopython-1.74-foss-2016b-Python-3.7.4'
    Search dependcies for versionsuffix == '-Python-%(pyver)s'
    add dep_exts are exts_list from dependent packages

    - remove the variable dep_eb
    - All to resolve dependancie in the FrameWork, FrameWork only
      needs a single argument. It had three.

    1.0.0 July 8, 2019
    framework.py becomes seperate package. Share code
    between easy_update and easy_annotate

    Read exts_list for R and Python listed in dependencies.
"""

__version__ = '1.0.1'
__maintainer__ = 'John Dey jfdey@fredhutch.org'
__date__ = 'July 9, 2019'


class FrameWork:
    """provide access to EasyBuild Config file variables
    name, version, toolchain, eb.exts_list, dependencies, modulename, biocver,
    methods:
        print_update()
    """
    def __init__(self, args):
        self.debug = False
        self.code = None
        self.pyver = None
        self.search_pkg = None
        self.indent_n = 4
        self.indent = ' ' * self.indent_n
        self.ptr_head = 0
        self.modulename = None

        full_path = os.path.dirname(args.easyconfig)
        (head, tail) = os.path.split(full_path)
        while tail:
            if 'easyconfig' in tail:
                self.base_path = os.path.join(head, tail)
                break
            (head, tail) = os.path.split(head)

        # update EasyConfig exts_list or check single package
        if args.easyconfig:
            eb = self.parse_eb(args.easyconfig, primary=True)
            self.exts_list = eb.exts_list
            self.toolchain = eb.toolchain
            self.name = eb.name
            self.version = eb.version
            if eb.name == 'Python':
                self.pyver = eb.version
            else:
                self.pyver = None
            self.modulename = eb.name + '-' + eb.version
            self.modulename += '-' + eb.toolchain['name']
            self.modulename += '-' + eb.toolchain['version']
            self.interpolate = {'name': eb.name, 'namelower': eb.name.lower(),
                                'version': eb.version,
                                'pyver': None,
                                'rver': None}
            self.dep_exts = self.parse_dependencies(eb, self.name)
            # exts_defaultclass = 'PythonPackage' | 'RPackage' | 'PerlModule'
            try:
                self.versionsuffix = eb.versionsuffix
                self.modulename += eb.versionsuffix
            except (AttributeError, NameError):
                self.versionsuffix = None
            self.modulename = self.modulename % self.interpolate
            if self.debug:
                sys.stderr.write('debug - modulename: %s\n' % self.modulename)
                sys.stderr.write('debug -       file: %s\n' % filename[:-3])
            try:
                self.dependencies = eb.dependencies
            except (AttributeError, NameError):
                self.dependencies = None
            try:
                self.biocver = eb.biocver
            except (AttributeError, NameError):
                pass
            self.check_eb_package_name(args.easyconfig)
            self.out = open(args.easyconfig[:-3] + ".update", 'w')

    def parse_eb(self, file_name, primary):
        """ interpret EasyConfig file with 'exec'.  Interperting fails if
        constants that are not defined within the EasyConfig file.  Add
        undefined constants to <header>.
        """
        header = 'SOURCE_TGZ  = "%(name)s-%(version)s.tgz"\n'
        header += 'SOURCE_TAR_GZ = "%(name)s-%(version)s.tar.gz"\n'
        header += 'SOURCELOWER_TAR_GZ = "%(namelower)s-%(version)s.tar.gz"\n'
        header += ('PYPI_SOURCE = "https://pypi.python.org/packages/' +
                   'source/%(nameletter)s/%(name)s"\n')
        header += ('SOURCEFORGE_SOURCE = "https://download.sourceforge.net/' +
                   '%(namelower)s"\n')
        eb = imp.new_module("EasyConfig")
        try:
            with open(file_name, "r") as f:
                code = f.read()
        except IOError as err:
            print("opening %s: %s" % (file_name, err))
            sys.exit(1)
        try:
            exec (header + code, eb.__dict__)
        except Exception as err:
            print("interperting EasyConfig error: %s" % err)
            sys.exit(1)
        if primary:     # save original text of source code
            self.code = code
        return eb

    def build_dep_filename(self, eb, dep, pyver=None):
        """build a filename from a dependencie objecwt"""
        dep_filename = '{}-{}-{}-{}'.format(dep[0], dep[1],
                                            eb.toolchain['name'],
                                            eb.toolchain['version'])
        if pyver and len(dep) > 2:
            versionsuffix = dep[2] % {'pyver': pyver}
            dep_filename += '{}'.format(versionsuffix)
        dep_filename += '.eb'
        sys.stderr.write(" - dependency: {}\n".format(dep_filename))
        return dep_filename

    def find_easyconfig(self, easyconfig):
        """ search base_path for easyconfig filename """
        found = None
        for r,d,f in os.walk(self.base_path):
            for filename in f:
                 if filename == easyconfig:
                      found = os.path.join(r,filename)
                      break
        return found

    def parse_dependencies(self, eb, lang):
        """ inspect dependencies for R and Python easyconfigs,
        if found add the exts_list to the list of dependent
        exts  <dep_exts>
        """
        try:
            dependencies = eb.dependencies
        except NameError:
            return None
        dep_exts = []
        for dep in dependencies:
            dep_filename = None
            if dep[0] in ['R', 'Python']:
                if lang == dep[0]:
                   dep_filename = self.build_dep_filename(eb, dep)
            if lang == 'Python' and len(dep) > 2 and dep[2] == '-Python-%(pyver)s':
                """ this is a PythonBundle """
                # dep_exts.extend([(dep[0], dep[1])]) # Explictly add the module
                dep_filename = self.build_dep_filename(eb, dep, self.pyver)
            if dep_filename:
                easyconfig = self.find_easyconfig(dep_filename)
                if easyconfig:
                    eb = self.parse_eb(str(easyconfig), False)
                    try:
                        dep_exts.extend(eb.exts_list)
                    except (AttributeError, NameError):
                        dep_exts.extend([dep])
        return dep_exts

    def parse_python_deps(self, eb):
        """ Python EasyConfigs can have other Python packages in the
        dependancy field. check 3rd element for -Python-%(pyver)s"
        """
        try:
            dependencies = eb.dependencies
        except NameError:
            return None
        for dep in eb.dependencies:
            pass

    def check_eb_package_name(self, filename):
        """" check that easybuild filename matches package name
        easyconfig is original filename
        """
        eb_name = os.path.basename(filename)[:-3]
        if eb_name != self.modulename:
            sys.stderr.write("Warning: file name does not match easybuild " +
                             "module name\n"),
        if eb_name != self.modulename or self.debug:
            sys.stderr.write("   file name: %s\n module name: %s\n" % (
                eb_name, self.modulename))

    def write_chunk(self, indx):
        self.out.write(self.code[self.ptr_head:indx])
        self.ptr_head = indx

    def rewrite_extension(self, pkg):
        name = pkg['name']
        name_indx = self.code[self.ptr_head:].find(name)
        name_indx += self.ptr_head + len(name) + 1
        indx = self.code[name_indx:].find("'") + name_indx + 1
        self.write_chunk(indx)
        self.out.write("%s'" % pkg['version'])
        self.ptr_head = self.code[self.ptr_head:].find("'") + self.ptr_head + 1
        indx = self.code[self.ptr_head:].find('),') + self.ptr_head + 3
        self.write_chunk(indx)

    def output_module(self, lang, pkg):
        """write exts_list entry
        """
        output = None
        if lang == 'R':
            output = "%s('%s', '%s')," % (self.indent, pkg['name'], pkg['version'])
        elif lang == 'Python':
            pkg_fmt = self.indent + "('%s', '%s', {\n"
            item_fmt = self.indent + self.indent + "'%s': %s,\n"
            output = pkg_fmt % (pkg['name'], pkg['version'])
            for item in pkg.keys():
                if item in ['name', 'version', 'action', 'type', 'orig_ver',
                            'processed', 'meta', 'level', 'spec']:
                    continue
                output += item_fmt % (item, pkg[item])
            for item in pkg['spec'].keys():
                output += item_fmt % (item, pkg['spec'][item])
            output += self.indent + "}),"
        return output

    def print_update(self, lang, exts_list):
        """ this needs to be re-written in a Pythonesque manor
        if module name matches extension name then skip
        """
        indx = self.code.find('exts_list')
        indx += self.code[indx:].find('[')
        indx += self.code[indx:].find('\n') + 1
        self.write_chunk(indx)

        for extension in exts_list:
            name = extension['name']
            if 'action' not in extension:
                sys.stderr.write('No action: %s\n' % name)
                extension['action'] = 'keep'

            if self.name.lower() == name.lower():
                # special case for bundles, if "name" is used in exts_list
                indx = self.code[self.ptr_head:].find('),') + 2
                indx += self.ptr_head
                self.write_chunk(indx)
            elif extension['type'] == 'base':  # base library with no version
                indx = self.code[self.ptr_head:].find(name)
                indx += self.ptr_head + len(name) + 2
                self.write_chunk(indx)
            elif extension['action'] in ['keep', 'update']:
                self.rewrite_extension(extension)
            elif extension['action'] == 'duplicate':
                print('Duplicate: %s' % name)
                name_indx = self.code[self.ptr_head:].find(name)
                name_indx += self.ptr_head + len(name)
                indx = self.code[name_indx:].find('),') + name_indx + 3
                self.ptr_head = indx
                continue
            elif extension['action'] in ['add', 'dep']:
                output = self.output_module(lang, extension)
                self.out.write("%s\n" % output)
        self.out.write(self.code[self.ptr_head:])


class UpdateExts:
    """
    """
    def __init__(self, args, eb, dep_eb):
        """
        """
        self.verbose = args.verbose
        self.debug = True 
        self.tree = args.tree
        self.meta = args.meta
        self.search_pkg = args.search_pkg
        self.ext_counter = 0
        self.pkg_update = 0
        self.pkg_new = 0
        self.pkg_duplicate = 0
        self.indent_n = 4
        self.indent = ' ' * self.indent_n
        self.ext_list_len = 1
        self.exts_dep = list()
        self.checking = list()  # pytest -> attrs -> pytest
        self.depend_exclude = list()
        self.exts_processed = list()

        if dep_eb:
            for exten in dep_eb.exts_list:
                if isinstance(exten, tuple):
                    if len(exten) == 3 and 'modulename' in exten[2]:
                        self.exts_dep.append(exten[2]['modulename'])
                    else:
                        self.exts_dep.append(exten[0])
                else:
                    self.exts_dep.append(exten)
        if args.easyconfig:
            self.exts_orig = eb.exts_list
            self.interpolate = {'name': eb.name, 'namelower': eb.name.lower(),
                                'version': eb.version}
        if self.search_pkg:
            self.name = args.name
            self.search_pkg = args.search_pkg
            if args.biocver:
                self.biocver = args.biocver
            if args.pyver:
                self.version = args.pyver
            elif args.rver:
                self.version = args.rver
            else:
                print('Language and version must be specified with ' +
                      '[--pyver x.x | --rver x.x | --biocver x.x]')
            self.sea_pkg = {'name': args.search_pkg, 'version': '', 'type': 'orig',
                            'level': 0, 'meta': {}}
            # move exts from easyconfig into processed list
            if args.easyconfig:
                for ext in self.exts_orig:
                    if isinstance(ext, tuple):
                        name = ext[0] % self.interpolate
                        version = ext[1] % self.interpolate
                        pkg = {'name': name, 'version': version}
                        self.processed(pkg)
                    else:
                        self.processed({'name': ext, 'type': 'base'}) 

    def is_processed(self, pkg):
        """ check if package has been previously processed
            if package exists AND is in the original exts_lists
                Mark as 'duplicate'
        updated July 2018
        """
        name = pkg['name']
        found = False
        if name in self.exts_dep:
            found = True
        elif name in self.checking:
            found = True
        else:
            for p_pkg in self.exts_processed:
                if 'spec' in p_pkg and 'modulename' in p_pkg['spec']:
                    modulename = p_pkg['spec']['modulename']
                else:
                    modulename = ''
                if (str(name) == str(p_pkg['name'])) or (name == modulename):
                    found = True
                    break
        if found:
            if pkg['type'] == 'orig':
                pkg['action'] = 'duplicate'
                self.pkg_duplicate += 1
                self.ext_counter -= 1
                self.processed(pkg)
                if self.verbose:
                    self.print_status(pkg)
                    return True
        return found

    def processed(self, pkg):
        """ save Processed packages
        save a normalize version of package name to <exts_search> for Python
        updated July 2018
        """
        pkg['processed'] = True
        pkg2 = dict(pkg)
        self.exts_processed.append(pkg2)

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
        """ print meta data from CRAN or Pypi
        :param info: dict
        """
        pass

    def check_package(self, pkg):
        """query package authority [Pypi, CRAN, Bio] to get the latest version
        information for a package. This is the heart of the program.

        input: pkg{}
        check that all dependencies are meet for each package.
        check_package can be called recursively.
        pkg['type'] is used to track status.
          - 'orig' is also used to track recursion
          - 'dep' package that is added as result of dependancy
          - 'add' packages read from file
        pkg['action'] What action will be take to exts_list.
          - 'add'; new package
          - 'keep'; no update required
          - 'update'; version change
          - 'duplicate' package appears twice
        """
        if self.debug:
            sys.stderr.write('check_package: %s\n' % pkg['name'])
        if self.is_processed(pkg):
            return
        self.checking.append(pkg['name'])
        status = self.get_package_info(pkg)
        if status in ["error", 'not found']:
            if pkg['type'] == 'orig':
                pkg['action'] = 'keep'
                self.processed(pkg)
                return
            else:
                msg = " Warning: %s is dependency, but can't be found!"
                print(msg % pkg['name'])
                return
        if 'version' in pkg['meta']:
            version = pkg['meta']['version']
        else:
            print('version not in %s' % pkg['name'])
            version = pkg['version']
        if pkg['version'] == version:
            pkg['action'] = 'keep'
        else:
            pkg['orig_ver'] = pkg['version']
            pkg['version'] = pkg['meta']['version']
            if pkg['type'] == 'orig':
                pkg['action'] = 'update'
                self.pkg_update += 1
            elif pkg['type'] in ['dep', 'add']:
                if self.debug:
                    print('check_package - action %s: dep or add' % pkg['name'])
                pkg['action'] = 'add'
                self.pkg_new += 1

        if 'requires' in pkg['meta'] and pkg['meta']['requires'] is not None:
            # sys.stderr.write('%s: %s\n' % (pkg['name'], pkg['meta']['requires']))
            for depend in pkg['meta']['requires']:
                if depend not in self.depend_exclude:
                    dep_pkg = {'name': depend, 'version': 'x', 'type': 'dep',
                               'spec': {}, 'meta': {}, 'level': pkg['level']+1}
                    self.check_package(dep_pkg)
        self.processed(pkg)
        if self.search_pkg:
            output = self.output_module(pkg)
            print(output)
        if self.verbose:
            self.print_status(pkg)
        if self.meta:
            self.print_meta(pkg['meta'])
        if pkg['action'] == 'add':
            self.ext_counter += 1

    def updateexts(self):
        """Loop through exts_list and check which packages need to be updated.
        this is an external method for the class
        """
        if self.search_pkg:
            self.check_package(self.sea_pkg)
        else:
            self.ext_list_len = len(self.exts_orig)
            for ext in self.exts_orig:
                self.ext_counter += 1
                if isinstance(ext, tuple):
                    name = ext[0] % self.interpolate
                    version = ext[1] % self.interpolate
                    pkg = {'name': name, 'version': version, 'type': 'orig',
                           'level': 0, 'meta': {}}
                    if len(ext) > 2:
                        pkg['spec'] = ext[2]
                    pkg['meta'] = {}
                    self.check_package(pkg)
                else:
                    self.processed({'name': ext, 'type': 'base'})
            if self.verbose:
                self.stats()

    def stats(self):
        sys.stderr.write("Updated Packages: %d\n" % self.pkg_update)
        sys.stderr.write("New Packages: %d\n" % self.pkg_new)
        sys.stderr.write("Dropped Packages: %d\n" % self.pkg_duplicate)

    def get_package_info(self, pkg):
        pass


class UpdateR(UpdateExts):
    """extend UpdateExts class to update package names from CRAN and BioCondutor
    """
    def __init__(self, args, eb, deps_eb):
        UpdateExts.__init__(self, args, eb, deps_eb)
        self.debug = False
        self.bioc_data = {}
        self.depend_exclude = ['R', 'base', 'compiler', 'datasets', 'graphics',
                               'grDevices', 'grid', 'methods', 'parallel',
                               'splines', 'stats', 'stats4', 'tcltk', 'tools',
                               'utils', ]
        try:
            self.biocver = args.biocver
        except NameError:
            self.biocver = None
        if not self.biocver:
            try:
                self.biocver = eb.biocver
            except (AttributeError, NameError):
                self.biocver = None
                print('BioCondutor version: biocver not set')
        if self.biocver:
            self.read_bioconductor_packages()
        self.updateexts()
        if not self.search_pkg:
            eb.print_update('R', self.exts_processed)

    def read_bioconductor_packages(self):
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

    def get_cran_info(self, pkg):
        """ MD5sum, Description, Package, releases[]
        """
        cran_list = "http://crandb.r-pkg.org/"
        resp = requests.get(url=cran_list + pkg['name'])
        if resp.status_code != 200:
            return "not found"
        cran_info = resp.json()
        pkg['meta']['info'] = cran_info
        if self.meta:
            self.print_meta(cran_info)
        pkg['meta']['version'] = cran_info['Version']
        if u'License' in cran_info and u'Part of R' in cran_info[u'License']:
            return 'base package'
        pkg['meta']['requires'] = []
        if u"LinkingTo" in cran_info:
            pkg['meta']['requires'].extend(cran_info[u"LinkingTo"].keys())
        if u"Depends" in cran_info:
            pkg['meta']['requires'].extend(cran_info[u"Depends"].keys())
        if u"Imports" in cran_info:
            pkg['meta']['requires'].extend(cran_info[u"Imports"].keys())
        return 'ok'

    def get_bioc_info(self, pkg):
        """Extract <Depends> and <Imports> from BioCondutor json metadata
        Example:
        bioc_data['pkg']['Depends']
                    [u'R (>= 2.10)', u'BiocGenerics (>= 0.3.2)', u'utils']
        interesting fields from BioCoductor:
        bioc_data['pkg']['Depends', 'Imports', 'Biobase', 'graphics', 'URL']
        """
        status = 'ok'
        if pkg['name'] in self.bioc_data:
            pkg['meta']['version'] = self.bioc_data[pkg['name']]['Version']
            if 'LinkingTo' in self.bioc_data[pkg['name']]:
                pkg['meta']['requires'].extend(
                    [re.split('[ (><=,]', s)[0]
                     for s in self.bioc_data[pkg['name']]['LinkingTo']])
            if 'Depends' in self.bioc_data[pkg['name']]:
                pkg['meta']['requires'].extend(
                    [re.split('[ (><=,]', s)[0]
                     for s in self.bioc_data[pkg['name']]['Depends']])
            if 'Imports' in self.bioc_data[pkg['name']]:
                pkg['meta']['requires'].extend(
                    [re.split('[ (><=,]', s)[0]
                     for s in self.bioc_data[pkg['name']]['Imports']])
        else:
            status = "not found"
        return status

    def print_depends(self, pkg):
        """ used for debugging """
        for p in pkg['meta']['requires']:
            if p not in self.depend_exclude:
                print("%20s : requires %s" % (pkg['name'], p))

    def get_package_info(self, pkg):
        """R version, check CRAN and BioConductor for version information
        """
        if self.debug:
            print('get_package_info: %s' % pkg['name'])
        pkg['meta']['requires'] = []
        status = self.get_bioc_info(pkg)
        if status == 'not found':
            status = self.get_cran_info(pkg)
        if self.debug:
            self.print_depends(pkg)
        return status

    def print_meta(self, meta):
        """Display metadata from CRAN"""
        for tag in meta:
            if tag == 'info':
                for md in meta['info']:
                    print("%s: %s" % (md, meta['info'][md]))
            else:
                print("%s: %s" % (tag, meta[tag]))

    def output_module(self, pkg):
        output = "%s('%s', '%s')," % (self.indent, pkg['name'],
                                      pkg['version'])
        return output


class UpdatePython(UpdateExts):
    """extend ExtsList class to update package names from PyPI
    Python Issues
       There are many small inconsistancies with PyPi which make it difficult
       to fully automate building of EasyConfig files.
       - dependancy checking - check for extras=='all'
       - pypi projects names do not always match module names and or file names
         project: liac-arff, module: arff,  file name: liac_arff.zip
    """
    def __init__(self, args, eb, deps_eb):
        UpdateExts.__init__(self, args, eb, deps_eb)
        self.debug = True
        self.pkg_dict = None
        if eb:
            (nums) = eb.version.split('.')
        else:
            (nums) = args.pyver.split('.')
        self.python_version = "%s.%s" % (nums[0], nums[1])
        self.pymajornum = nums[0]
        self.pyminor = nums[1]
        # Python >3.3 has additional built in modules
        if nums[0] == 3 and nums[1] > 3:
            self.depend_exclude.extends(['argparse', 'asyncio'])
        if self.debug and self.search_pkg:
            print('Python Search PyPi: %s' % self.search_pkg)
        self.updateexts()
        if not self.search_pkg:
            eb.print_update('Python', self.exts_processed)

    def get_pypi_pkg_data(self, pkg, version=None):
        """
        return meta data from PyPi.org
        """
        if version:
            req = 'https://pypi.org/pypi/%s/%s/json' % (pkg['name'], version)
        else:
            req = 'https://pypi.org/pypi/%s/json' % pkg['name']
        resp = requests.get(req)
        if resp.status_code != 200:
            msg = "API error: %s GET release %s\n"
            sys.stderr.write(msg % (resp.status_code, pkg['name']))
            return 'not found'
        project = resp.json()
        # verify that package name is correct
        if pkg['name'] != project['info']['name']:
            sys.stderr.write('package name mismatch: %s -> %s\n' % (
                pkg['name'], project['info']['name']))
            pkg['name'] = project['info']['name']
        return project

    def check_package_name(self, pkg_name):
        """
        verify that package name from EasyConfig
        matches package name from PyPi
        """
        pkg = {'name': pkg_name}
        response = self.get_pypi_pkg_data(pkg)
        if response == 'not found':
            return response
        else:
            return response['info']['name']

    def parse_pypi_requires(self, requires):
        """pypi requires_dist PEP - 345
        https://dustingram.com/articles/2018/03/05/why-pypi-doesnt-know-dependencies/
        The project name must be as specified at pypi.org.
        requires_dist: <name> <version>[; Environment Markers]

        Only install the latest version so ignore all version information
        input: 'numpy (>=1.7.1)'  output: 'numpy'

        Test that <python_version> and <sys_platform> conform.
        example:  (sys_platform=='win32')
        If <extra> is present use ignore_list to validate
        """
        ignore_list = ['dev', 'tests', 'docs']
        if requires is None:
            return []
        dists = []
        #  sys_platform = 'Linux'
        #  python_version = self.python_version
        extra_re = re.compile("extra\s?==\s?'([\w\-\.\*]+)'")
        platform_re = re.compile('sys_platform\s?==\s?\'(\w+)')
        for req in requires:
            ignore = False
            pkg_name = req.split()[0]
            pkg_name = pkg_name.rstrip(';')
            pypi_name = self.check_package_name(pkg_name)
            if pypi_name == 'not found':
                continue
            # test for Environment Marker (stuff after ;)
            fields = req.split(';')
            if len(fields) > 1:
                platform = re.findall(platform_re, fields[1])
                extras = re.findall(extra_re, fields[1])
                for extra in ignore_list:
                    if extra in extras:
                        ignore = True
                for x in platform:
                    if x == 'win32':
                        ignore = True
            if self.is_processed(pkg={'name': pypi_name,
                                      'type': 'dep'}):
                continue
            if not ignore and pypi_name not in dists:
                dists.append(pypi_name)
        return dists

    def print_meta(self, meta):
        """ Display info from pypi.org
        """
        if self.verbose:
            tags = meta.keys()
        else:
            tags = ['filename', 'packagetype', 'url', 'python_version',
                    'requires_dist', 'summary', 'requires_python',
                    'classifiers', 'description', 'platform']
        for key in tags:
            if key == 'description':
                print("%s: %s" % (key, meta[key][1:60]))
            elif key == 'requires_dist':
                print('%s:' % 'requires_dist')
                for req in meta[key]:
                    print("   %s" % req)
            else:
                print("%s: %s" % (key, meta[key]))

    def get_package_info(self, pkg):
        """Python version
           Get package meta data via PyPi.org
           pkg is a dict; {'name', 'version', 'spec'}
           return metadata as dict
             pkg['meta']['version']
           if source package is not found look for whl
           if pyver == ['3.5, '3.6', '3.7']:
               arch = 'linux' ['manylinux', 'anylinux', 'linux']
        """
        status = self.get_pypi_info(pkg)
        return status

    def get_pypi_release(self, pkg, version):
        """if source dist is not available from pypi search
        the release for a wheel file.
        """
        release = self.get_pypi_pkg_data(pkg, version)
        if release == 'not found':
            return 'not found'
        cplist = ['cp35', 'cp36', 'cp37']
        for rel in release['releases'][version]:
            if any(cver in rel['python_version'] for cver in cplist):
                if 'manylinux' in rel['filename']:
                    pkg['meta'].update(rel)
                    return 'ok'
        return 'not found'

    def get_pypi_info(self, pkg):
        """get version information from pypi.  If <pkg_name> is not processed
        seach pypi. pkg_name is now case sensitive and must match
        info['digests']['sha256'], 'summary', 'url', 'filename', 'home_page'
        """
        project = self.get_pypi_pkg_data(pkg)
        if project == 'not found':
            return 'not found'
        status = 'not found'
        pkg['meta'].update(project['info'])
        new_version = pkg['meta']['version']
        requires = project['info']['requires_dist']
        pkg['meta']['requires'] = self.parse_pypi_requires(requires)
        for ver in project['releases'][new_version]:
            if 'packagetype' in ver and ver['packagetype'] == 'sdist':
                pkg['meta']['url'] = ver['url']
                pkg['meta']['filename'] = ver['filename']
                status = 'ok'
                break
        # one last try to find package release data
        if status != 'ok':
            status = self.get_pypi_release(pkg, new_version)
        # only set this if not set
        if ('spec' in pkg and 'source_urls' not in pkg['spec'] and
                new_version != pkg['version']):
            url = "['https://pypi.io/packages/source/%s/%s']"
            pkg['spec']['source_urls'] = url % (pkg['name'][0], pkg['name'])
        if self.meta:
            self.print_meta(project['info'])
            sys.exit(0)
        return status

    def output_module(self, pkg):
        pkg_fmt = self.indent + "('%s', '%s', {\n"
        item_fmt = self.indent + self.indent + "'%s': %s,\n"
        if self.tree:
            spaces = self.indent * pkg['level']
            output = '%s%s' % (spaces, pkg['name'])
        else:
            output = pkg_fmt % (pkg['name'], pkg['version'])
            for item in pkg.keys():
                if item in ['name', 'version', 'action', 'type', 'orig_ver',
                            'processed', 'meta', 'level', 'spec']:
                    continue
                output += item_fmt % (item, pkg[item])
            if 'spec' in pkg:
                for item in pkg['spec'].keys():
                    output += item_fmt % (item, pkg['spec'][item])
            output += self.indent + "}),"
        return output


if __name__ == '__main__':
    """ create unit test """
    print('none') 
