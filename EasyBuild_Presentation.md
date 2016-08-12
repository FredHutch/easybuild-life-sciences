# EasyBuild

.footer: *Presented at Seattle Information Technology Exchange 2016-08-17*

---

# EasyBuild is building software with ease
 
*Features*
 
- multiple versions and builds of many software packages
- automated dependency system
- re-producible software package builds
- automated environment modules managment
- curated toolchains


---

# EasyBuild Terms

- **environment modules:** a system of managing environment variables
- **easyconfig:** a file that describes a software package version
- **toolchain:** a defined collection of compilers and support libraries

---

# Environment Modules
 
Easybuild uses Environment Modules to automatically "wrap" software packages it builds.

This is how it makes software packages available, isolated from eachother, and able to be used together.

Environment Modules is a system that:
 
- sets and unsets environment variables
- avoids conflicts
- provides administrative hooks
 
---
 
# Environment Variables: shell integration
 
Environment Variables integrates with the users shell

- shell function module()
- calls modulecmd with options and shell name
- modulecmd outputs commands for specified shell to set/unset variables
- we recommend the use of Lmod as it is the most widely used by EasyBuilders

---
 
# Example of module use
 
Here is a quick example of a module in use:
 
> $ which R
> /usr/bin/R
> $ module load R/3.3.0-intel-2016a
> $ which R
> /app/easybuild/software/R/3.3.0-intel-2016a/bin/R

---

# Easyconfigs

Easyconfigs are text files that define a software package version + toolchain.

> easyblock = 'ConfigureMake'
>
> name = 'make'
> version = '4.1'
>
> homepage = 'http://www.gnu.org/software/make/make.html'
> description = "GNU version of make utility"
>
> toolchain = {'name': 'GCC', 'version': '4.9.2'}
>
> source_urls = [GNU_SOURCE]
> sources = [SOURCE_TAR_BZ2]
>
> moduleclass = 'devel'

---

# Markdown font size options

# header 1
## header 2
### header 3
#### header 4
##### header5
###### header 6

Normal text

# Why do I care about toolchains?

Toolchains are simply grouping of EasyBuilt software packages used to supply a consistent set of compilers, options, and libraries with which to build additional software packages.

---

# Toolchain consistency

Using the same toolchain enables:

- more efficient loading of modules for multiple software packages
- package stacking (bundles)
- **re-producibility**

---

# Toolchains and performance

Toolchains are curated and have been built by the EasyBuild community to provide optimized performance

*rbench:*

R | EB GCC R | EB Intel R
---|---|---
res | res | res

# Example of an EasyBuild

Demo a quick EasyBuild in shell

---

# Using EasyBuild

We talked about toolchains, how and why you might want to build with EasyBuild
 

 
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
- EasyBuild PREFIX owned by group, EASYBUILD ENV VARS to support group building
- EasyBuilds done on a single machine (Intel license, sanity)

# Bootstrapping EasyBuild

> # pick an installation prefix to install EasyBuild to (change this to your liking)
> EASYBUILD_PREFIX=$HOME/.local/easybuild
>
> # download script
> curl -O https://raw.githubusercontent.com/hpcugent/easybuild-framework/develop/easybuild/scripts/bootstrap_eb.py
>
> # bootstrap EasyBuild
> python bootstrap_eb.py $EASYBUILD_PREFIX
>
> # update $MODULEPATH, and load the EasyBuild module
> module use $EASYBUILD_PREFIX/modules/all
> module load EasyBuild

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

swap zlib for a binary with 'which' example

new slide: R dependency example

new slides: easybuild metrics/stats - totals, FH totals + packages for labs, life sciences

new slide: easyconfig simple example

new slide: easyconfig complex example

new slide: collaboration and community

new slide: new easyconfig example

new slide: easybuild github and easyconfigs search

new slide: contribution, mailing list, community, processes, testing, github-based

new slide: outlook/next steps - lifesciences easybuild community for collaboration using github,

new slide: containers, dockers, contrast and compare

new slide: reproducibility/citation lxc/docker containers
