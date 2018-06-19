---
layout: default
title: Using Software 
permalink: /software
---

### Scientific Software Environment
This repository contains the tools and recipes used to build software for the Fred Hutch. The Hutch
uses Easybuild software to build and manage software. Easybuild is a software build
and installation framework for managing scientific software. The user interface for using software is Modules.
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

 - ```ml Python/3.6.5-foss-2016b-fh1``` Load a Python module 
 - ```ml list``` List currently loaded modules
 - ```ml purge``` Unload all modules
 - ```ml avail Python``` Show all available Python modules

To you modules from a script the modules environment needs to be setup. Place the following
in your bash or sbatch scripts.

```
source /app/Lmod/lmod/lmod/init/bash
module use /app/easybuild/modules/all
module load R/3.5.0-foss-2016b-fh1
```

### Software Requests
Please send software requests to SciComp support.  
Requests can be made for additional packages to be added for R, Python and for new software. 
Packaged modules of R and Python are made as soon as possible after new releases are announced. 
