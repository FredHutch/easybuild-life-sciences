# EasyBuild

---

# Fred Hutch Researchers

## What they want:
1. reliable systems and hardware
2. reliable software:
  - current versions built consistently and quickly
  - software and version accounting for all jobs
  - reliable citation sources -> re-producibility

---

# Reliable Systems & Hardware

## Out of scope (call us)

---

# Reliable Software

## EasyBuild

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

# EasyBuild Benefits

- 1058 Software Packages (6600+ package versions)
- 54 toolchains
- 150+ contributors
- 87% coverage of our existing software stack

---

# EasyBuild Understanding

- Terms & Components
- Building software
- Fred Hutch Implementation
- Community and Collaboration

---

# EasyBuild Terms

- **environment modules**
- **easyconfigs**
- **toolchains**

---

# Environment Modules

---

# Environment Module Overview
 
- sets and unsets environment variables
- avoids conflicts between software packages
- provides dependency resolution
- provides administrative hooks
 
---
 
# How It Works
 
Environment Variables integrates with the users shell

- shell function/alias module()
- calls `modulecmd`
- `modulecmd` echoes shell cmds to be eval'd

---

# Example of module use
 
     !bash
     $ which R
     /usr/bin/R
     $ module load R/3.3.0-intel-2016a
     $ which R
     /app/easybuild/software/R/3.3.0-intel-2016a/bin/R

---

# Dependency Example

![Module List](module_list.png)

---

# Shell Env Example

![Module List](path.png)

---

# Easyconfigs

---

# Easyconfigs

Easyconfigs...

- are text files (python)
- define a software package + version + toolchain
- can define extentions/modules/libraries
- can define options and environment variables

---

# Easyconfig Example

## Python

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

![Module List](r_easyconfig.png)

---

# Toolchains

---

# What are toolchains exactly?

- are defined in easyconfig files
- are a collection of compilers and support libraries
- provide consistent build parameters/env
- are a great place for optimization fan-out

---

# How far does it go?

Down:

- just above kernel/hardware
- is inconsistent

Up:

- build tools (compilers, static base libraries)

---

# Toolchain example

![Toolchain Example](foss_toolchain.png)

---

# Why do you care about Toolchains?

- optimize module stacks
- optimize software
- aid in *re-producibility*

---

# Toolchains and performance

rbench:

- *Ubuntu R*: secs
- EasyBuild *foss-2016a R*: secs (nn% faster)
- EasyBuild *intel-2016a R*: secs (nn% faster)

---

# Re-producibility

---

# Common citation

Software is commonly cited using some combination of:

   - Author Name
   - Project Name
   - URL

---

# EasyBuild possible citation

- software + version + toolchain + release
- compiler + libraries + parameters + options

The easyconfig will allow one to re-build precisely.

---

# EasyBuilding Software

---

# Prerequisites

- Python
- Environment Modules

Use *Lmod* if modules are new for you.

---

# Bootstrapping EasyBuild

    !bash
    $ curl -O https://raw.github.../bootstrap_eb.py
    $ python bootstrap_eb.py $EASYBUILD_PREFIX
    $ module use $EASYBUILD_PREFIX/modules/all
    $ module load EasyBuild

---

# build something

Let's build something

---

# EasyBuild @FredHutch

---

# Fred Hutch goals

- centralized, shared packages
- non-root, multiple, individual builders in a group
- fast release of new versions
- "fat" R with up-to-date libraries optimized

---

# Executive decisions

Toolchains:

- **foss-n:** Free Open Source Software - GCC and friends
- **intel-n:** Intel C and Fortran

Deployment:

- universal read-only NFS export
- common EasyBuild hierarchy
 
---

# Engineering Decisions

- NFS mounted /app on all systems (nfs ro)
- EasyBuild PREFIX owned by POSIX group
- EASYBUILD ENV VARS to support group building
- EasyBuild is easy-built; everything in /app
- EasyBuilds done on a build host (nfs rw)

---

# EasyBuild @FredHutch Metrics

- 218/757 software packages/versions built
- 4 builders (simultaneous)
- 86% built in 4 months
- 87% of our old software stack re-built

---

# What's the catch?

- writing easyconfigs is not easy
- existing community accepts all easyconfigs
- = proliferation of matrix
- = inconsistent support libraries
- = conflict

---

# FredHutch Next Step Goals

- build EB life-sciences community
- create github-based workflow for easyconfigs
- provide new versions quickly/publicly
- share tools and code
- reduce complexity of writing new easyconfigs

---

# Next Step Details

- help implement expanded toolchains
- publish easyconfigs
  - upstream
  - life-sciences github repo
- published detailed implementation example
- take ownership of R easyconfig
- implement EasyBuild in a container
