---
title: Scientific Software at Fred Hutch
layout: single
permalink: /about/
toc: true
toc_label: "On This Page"
sidebar:
  nav: "docs"
---

### Scientific Software at Fred Hutch

This website documents open source scientific software at the Fred Hutch that is maintained
by HPC.  These are the software modules that are available on the SciComp compute resources. The purpose of this site is to provide an up to date catalog of scientific software. Posts are created as new software is added and updated.  Use the [Software Updates](releases.html) link on the main page to check for recent updates. The Fred Hutch maintains local builds of R and Python that contain hundreds of packages. Provide a list of R and Python versions along with a detailed list for R and Python packages. 

### EasyBuild
FredHutch Scientific Computing is using Easybuild to create and document open source software (OSS). EasyBuild is a software build and installation framework.
EasyBuild is a community project supported by research centers from around the world.
EasyBuild builds software from a specification scripts named easyconfigs. Each version of every
software tool has a specific easyconfig file. Easyconfig follow a strict nameing 
scheme ```<name>-<version>[-<toolchain>][-<versionsuffix>].```
Easyconfigs are managed with GitHub.
The Hutch's HPC environment provides 100's of OSS
packages to our Scientists and staff. All software is built to offer high reproducibility, it can be rebuilt exactly even 10 years from now.

### Roadmap 
The Gizmo Cluster will be updated from Ubuntu Trusty to Bionic update for the summer of 2020. The upgrade will
nessitage rebuilding all software packages. Bionic packages are being built with the foss-2019b toolchain.
EasyBuild packages are built from a common set of build tools referred to as a toolchain.
Toolchains are standardized and updated twice per year, using the naming convention YEAR[a/b].
Beginning with the foss-2019b toolchain dependent libraries are being standradized for each toolchain for 
consistency within the toolchain.
The foss-2016b toolchain was used by Trusty version of Gizmo, and foss-2019b will be used for Bionic.

R and Python modules at the Hutch are customized to support our local user base.
The Hutch will continue to support custom R and Python modules but they will be based on
the published modules from EasyBuild. The custom packages will have a versionsuffix of "-fhx".
The "x" will be a version number that changes as new modules are added without changing
the versions of any of the existing packages. Example ```Python-3.7.4-foss-2019b``` is a base package
 that used by the EasyBuild community which only contains 44 packages. 
Our local build of Python ```3.7.4-foss-2019b-fh1``` has hundreds of packages based on users requests.

### Site Maintenance
This site [(EasyBuild-Life-Sciences)](http://fredhutch.github.io/easybuild-life-sciences)
is a GitHub Pages site hosted from GitHub [Easybuild-life-sciences](https://github.com/FredHutch/easybuild-life-sciences). The contents of the site are contained in the ```docs``` directory.

### Document a Module
When a new module is created or updated that is of interest to our user community it should be documented. Our EasyBuild environment contains hundreds of modules.  Many of these modules are system libraries that users do not interactive with.  We do not document libraries like ```libXcb```.
To document and a significate module like BLAST+  a notice should be made created and the software inventory should be updated.  Only moudles with a moduleclass of ``bio``` will be documented in the software inventory. But a notification can be made for any module.
```
# Perform all work with the repository directory tree.
cd  ~/easybuild-life-sciences
# Create a Post
./scripts/create_posts.sh BLAST+-2.7.1-foss-2018b.eb
# Update Software Inventory
./scripts/create_module_list.sh
Collecting Inventory
Generating Markdown
Wrote inventory to /home/jfdey/easybuild-life-sciences/docs/bio-modules-18.04.md
```
Notice that the inventory name is based on the Ubuntu version.

### Document R and Python Packages
```
# Perform all work with the repository directory tree.
cd  ~/easybuild-life-sciences
./scripts/easy_annotate.py fh_easyconfigs/Python-3.6.7-foss-2018b.eb
Package: Python-3.6.7-foss-2018b
# easy_annotate will create an Markdown file based on the package name
# Move the Markdown to the approreate directory
#  docs/_R or docs/_Python
mv Python-3.6.7-foss-2018b.md docs/_Python/
```
### Update Site
After your package is built and tested and the documentation is updated the changes need to be pushed to GitHub.
```
git add -A
git commit -m "new software tool build"
git push
```

