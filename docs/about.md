---
layout: default
title: Scientific Software at Fred Hutch
permalink: /about
---

### Scientific Software at Fred Hutch

This website documents open source scientific software at the Fred Hutch that is maintained
by HPC.  These are the software modules that are available on the SciComp compute resources. The purpose of this site is to provide an up to date catalog of scientific software. Posts are created as new software is added and updated.  Use the [Software Updates](releases.html) link on the main page to check for recent updates. The Fred Hutch maintains local builds of R and Python that contain hundreds of packages. Provide a list of R and Python versions along with a detailed list for R and Python packages. 

### EasyBuild
FredHutch Scientific Computing is using Easybuild to create and document open source software (OSS). EasyBuild is a software build and installation framework.
EasyBuild is a community project supported by research centers from around the world.
EasyBuild builds software from a software from a specification file called and
easyconfig. Easyconfigs are managed with GitHub.
software. The Hutch's HPC environment provides 100's of OSS
packages to our Scientists and staff. All software is built to offer high reproducibility, it can be rebuilt exactly even 10 years from now.

### Roadmap 
EasyBuild packages are built from a common set of build tools referred to as a toolchain.
The foss-2016b toolchain is the current predominate toolchain used at the Hutch. Foss-2016b modules are built on Ubuntu 14.04 Linux.
Hutch SciComp is planning a cluster an upgrade
to Ubuntu 8.04 in the Spring of 2019. To support the cluster upgrade the toolchain will
be updated to foss-2018b. Most of the Hutch's 2016b packages are unique to the
Hutch to due issues with local libraries, liberal usage of OS development packages and remapping of libraries for consistency within the toolchain.

The EasyBuild community will be maintaining consistent libraries dependencies for the foss-2018b toolchain. This is good news for the Hutch and other research centers. As we migrate all software packages to foss-2018b the Hutch will use modules as they are published from EasyBuild. If new package development is
necessary the results will be pushed upstream to EasyBuild.

R and Python modules at the Hutch are customized to support our local user base.
The Hutch will continue to support custom R and Python modules but they will be based on
the published modules from EasyBuild. The custom packages will have a versionsuffix of "-fhx".
The "x" will be a version number that changes as new modules are added without changing the versions of any of the existing packages. Example *Python-3.6.7-foss-2018b* is a base package that used by the EasyBuild community which only contains 44 packages.  *Python-3.6.7-foss-2018b-fh1* is our local build which has 471 packages.

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
Wrote inventory to /home/jfdey/easybuild-life-sciences/docs/bio-modules-14.04.md
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

