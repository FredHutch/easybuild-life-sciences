---
title: "Scientific Software"
permalink: /
layout: single
toc: true
toc_label: "On This Page"
---

## Fred Hutch Scientific Software
This repository contains the tools and recipes used to build software for the Fred Hutch. The Hutch
uses Easybuild frame to buiand and manage software. Easybuild is a software build
and installation framework for managing scientific software. The user interface for using software is Modules.
The module command is used to instantiate a specific software package. Easybuild and Modules 
together provide a method to make repeatable and documented scientific software. 

### Software Packages
Easybuild modules are build from toolkits. The toolkit provides the foundation for system level and
math libraries. Easybuild software modules names contain the software name-version and the toolkit name
used to build the module. Custom reciepies used at the Hutch also have suffix ```fhx``` to distigue them
from other published reciepeis. Once a reciepe is published to this repository it is not changed. If changes
have to be made to a package the suffix is versioned. (fh1, fh2, fh3 etc). 

### Using Modules
Modules can be loaded, unloaded, listed and searched.  The ``module load``` can be abrivated to ```ml```.
 - ml Python/3.6.5-foss-2016b-fh1  \# load a module
 - ml  \# list currently loaded modules
 - ml purge \# Unload all modules
 - ml avail Python  \# Show all available Python modules
 
### R 
Link to R page

### Python  
Link to Python Page

### Scientific Software
Link to all ```bio``` modules
