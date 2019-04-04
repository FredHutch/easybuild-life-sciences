---
title: Scientific Software 
layout: single
permalink: /
toc: true
toc_label: "On This Page"
sidebar:
  nav: "docs"
---

Fred Hutch maintains open source scientific software for use with HPC resouces
at the center. This site provides an inventory of available software and
information about using scientific software.  The Hutch uses Easybuild 
to build and manage software. Easybuild is a software build and installation
framework for managing scientific software. 

## Life Science Software Inventory
Inventory of module class "bio" software modules.  This is a subset of all modules that are available. The full
list of modules encompesses low level system libraries, math libraries, programming languages, and visualization tools. 

 - [Life Science Software Inventory]({{ site.baseurl }}/bio-modules-14.04/)

## Programming Language Modules
Scientific Computing maintains custom builds for R and Python. The
custom R and Python modules contain hundreds of packages. The package
list is a compilation of user requests.  Click on the links
bellow to see a list of available builds with a list of modules.

 - [R Modules]({{ site.baseurl }}/r/)
 - [Python Modules]({{ site.baseurl }}/python/)
 - [foss-2016b]({{ site.baseurl }}/foss-2016b/) (Free and Open Source Software)
 - [foss-2018b]({{ site.baseurl }}/foss-2018b/) (Free and Open Source Software)
 - [Dependency Graphs]({{ site.baseurl }}/dot/)

### Scientific Software Environment
The user interface for using software is Modules.
The module command is used to instantiate a specific software package.
Easybuild and Modules provide the tools for using scientific software in a
repeatable and documented way.

### Software Packages
Easybuild modules are built from toolkits. The toolkit provides the foundation for system level and
math libraries. Easybuild software modules names contain the software name-version and the toolkit name
used to build the module. Custom reciepies used at the Hutch also have suffix <tt>fh</tt> to distigue them
from other published reciepeis. Once a reciepe is published to this repository it is not changed. If changes
have to be made to a package the suffix is versioned. (fh1, fh2, fh3 etc).

### Using Modules

Modules can be loaded, unloaded, listed and searched.  The
```module load``` command can be abrivated to ```ml```.

Load a Python module
```ml Python/3.6.5-foss-2016b-fh1```

List currently loaded modules
```ml list```

Unload all modules
```ml purge```

Show all available Python modules
```ml avail Python``` 

### Using Modules in Scripts
Place the following in your bash or sbatch scripts to load modules within
your scripts.

```
source /app/Lmod/lmod/lmod/init/bash
module use /app/easybuild/modules/all
module load R/3.5.0-foss-2016b-fh1
```

### Software Requests
Please send software requests to SciComp support.
Requests can be made for additional packages to be added for R,
Python and for new software.
Packaged modules of R and Python are made as soon as possible 
after new releases are announced.
