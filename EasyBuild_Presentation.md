# EasyBuild

---

# EasyBuild is...

## building software with ease
 
### Features
 
- multiple versions and builds of many software packages
- automated dependency system
- re-producible software package builds
- automated environment modules managment
- curated toolchains

---

# EasyBuild Terms

- **environment modules:**

   a system of managing environment variables
- **easyconfig:**

   a file that describes a software package version
- **toolchain:**

   a defined collection of compilers and support libraries

---

# Environment Modules

---

# Environment Modules
##Overview
 
- sets and unsets environment variables
- avoids conflicts between software packages
- provides administrative hooks
 
---
 
# Environment Modules
##Shell integration
 
Environment Variables integrates with the users shell

- shell function/alias module()
- calls `modulecmd` with parameters and shell name from the alias
- modulecmd outputs commands for specified shell to set/unset variables
- we recommend the use of Lmod as it is the most widely used by EasyBuilders

---

# Environment Modules
##Example of module use
 
Here is a quick example of a module in use:
 
   !bash
   $ which R
   /usr/bin/R
   $ module load R/3.3.0-intel-2016a
   $ which R
   /app/easybuild/software/R/3.3.0-intel-2016a/bin/R

---

# Example

   !bash
   $ module list

---

# Easyconfigs

---

# Easyconfig Simple Example

- Easyconfigs are text files
- Easyconfigs define a software package + version + toolchain

   !python
   easyblock = 'ConfigureMake'
   name = 'make'
   version = '4.1'
   homepage = 'http://www.gnu.org/software/make/make.html'
   description = "GNU version of make utility"
   toolchain = {'name': 'GCC', 'version': '4.9.2'}
   source_urls = [GNU_SOURCE]
   sources = [SOURCE_TAR_BZ2]
   moduleclass = 'devel'

---

# Easyconfig Extended Example

extended R example

---

# What are toolchains exactly?

## Toolchains

- are defined in easyconfigs
- are a collection of compilers and support libraries
- provide consistent build parameters/env
- optimize module stacks
- aid in *re-producibility*

---

# EasyBuild Benefits

## Metrics

- 1058 Software Packages (6600+ package versions)
- 54 toolchains
- 150+ contributors

---

# Toolchains and performance

Toolchains are curated and have been built by the EasyBuild community to provide optimized performance

##rbench

- Base R: secs
- EasyBuild foss-2016a toolchain R: secs (nn% faster)
- EasyBuild intel-2016a toolchain R: secs (nn% faster)

---

# Example of an EasyBuild

Demo a quick EasyBuild in shell

---

# Using EasyBuild

We talked about toolchains, how and why you might want to build with EasyBuild
 
---
 
# EasyBuild @FredHutch
 
Fred Hutch goals:

- centralized, shared packages
- non-root, multiple, individual builders in a group

Toolchains:

- **foss-n:** Free Open Source Software - GCC and friends
- **intel-n:** Intel C and Fortran
 
---

# EasyBuild @FredHutch Details

- NFS mounted /app on all systems (ro)
- EasyBuild PREFIX owned by POSIX group
- EASYBUILD ENV VARS to support group building
- EasyBuild is easy-built; everything in /app
- EasyBuilds done on a build host

---

# EasyBuild @FredHutch Metrics

- nnn/mmm software packages/versions built
- 4 builders
- nn% of our old software stack re-built

---

# FredHutch Next Step Goals

- build EB life-sciences community
- provide new versions quickly/publicly
- share tools help promote code

---

# Next Step Details

- publish easyconfigs
  - upstream
  - life-sciences github repo
- published detailed implementation example
- take ownership of R easyconfig

# Bootstrapping EasyBuild

    !bash
       # pick an installation prefix to install EasyBuild to (change this to your liking)
    $ EASYBUILD_PREFIX=$HOME/.local/easybuild
       # download script
    $ curl -O https://raw.githubusercontent.com/hpcugent/easybuild-framework/develop/easybuild/scripts/bootstrap_eb.py
       # bootstrap EasyBuild
    $ python bootstrap_eb.py $EASYBUILD_PREFIX
       # update $MODULEPATH, and load the EasyBuild module
    $ module use $EASYBUILD_PREFIX/modules/all
    $ module load EasyBuild

---

# EasyBuild Environment

EasyBuild environment and options - paths, shared paths, Lmod hooks, logging

---

CHANGES

X move toolchain difference earlier - toolchain a nd reproducibility
X better transition between env modules and toolchains

X bootstrap slide - quick steps... "getting easybuild going quickly" is python

X easybuild environment: easybuid at FH - link to easybuild-lifesciences
 multi user
 shared paths
 no root

X swap zlib for a binary with 'which' example

new slide: R dependency example

new slides: easybuild metrics/stats - totals, FH totals + packages for labs, life sciences

X new slide: easyconfig simple example

X new slide: easyconfig complex example

X new slide: collaboration and community

new slide: new easyconfig example

new slide: easybuild github and easyconfigs search

X new slide: contribution, mailing list, community, processes, testing, github-based

X new slide: outlook/next steps - lifesciences easybuild community for collaboration using github,

new slide: containers, dockers, contrast and compare

new slide: reproducibility/citation lxc/docker containers
