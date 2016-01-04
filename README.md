# easybuild-life-sciences
Implementation and use of EasyBuild at FredHutch

# Overview
This repo will reflect our decisions and strategy during implementation of EasyBuild to replace our manual build procedure. Prior to Easybuild, we built software packages by hand, installed them onto an NFS share mounted read-only on our cluster, and documented the build procedures as best we could. These software packages were then "wrapped" using Environment Modules to allow users and scripts to use them without multiple manual PATH changes.

If you are not familiar with Easybuild, check out the [repo](https://github.com/hpcugent/easybuild) and [excellent documentation](http://easybuild.readthedocs.org/en/latest/). One of the benefits of Easybuild is that it will construct a proper modulefile for you during the build process.

# Goals
Before and during implementation, we kept the following goals in mind:

   * software packages will be reproducable
   * modules will be easily loaded by user in interactive sessions and in scripts
   * default versions of software packages will be easy to manage (ex: R-3.2.3 may be most recent, but `module load R` will load R-3.2.1)
   * packages will be built by any member of a given POSIX group
   * new packages will be easily implemented (new versions and software packages without existing easyconfigs)

# TODO
Some of our goals were not met in the initial implementation. Mostly due to unimplemented features in Environment Modules and/or EasyBuild itself.

*New package implementation* (i.e.: new easyconfigs) - Easybuild may be in transition with regard to easyconfig implementation. Initially the idea of reproducibility and separate easyconfig files went hand-in-hand. You build R version 3.2.1 with an Intel toolchain so you have an easyconfig detailing just that. At this time there are 738 software packages with easyconfigs in the public repo. If we figure only two toolchains and three versions for each package, we will need 4428 easyconfig files. The number of easyconfigs will become more difficult to manage, and certainly more difficult to choose. Easybuild supports several options in the "try" family like "--try-software-version=" and "--try-toolchain=" that can help re-build packages under a different toolchain or with a new version.

*Default Module Version Management* - In many cases, modules does a good job of picking the most recent modulefile to use. However in several cases (R, PYthon, Intel toolchain, etc.) we do not want the most recent package to be the default. Easybuild (AFAIK) does not provide a mechanism for managing the contents of ".version" in modulefiles directories. Nor should it. However, a written procedure like this does not leave me with a good feeling toward ending up with correctly managed default versions effortlessly. Perhaps a wrapper script. Or an easybuild option "--make-default" that would modify the ".version" file in the resulting modulefile directory.

# Implementation Notes

## Prerequisites

Clearly, some software is required to compile software. We started with a base of Ubuntu 14.04 LTS along with the `build-essentials` meta package. Also, an implementation of Modules is required. We were already using Environment Modules, but if you do not have a Module suite in use, please check out [Lmod](https://www.tacc.utexas.edu/research-development/tacc-projects/lmod) - it is under current development and offers a different and more advanced feature set.

## Environment

First, some notes about our environment:

   * read-only NFS mount on all systems at /app that currently has hand-build software packages and our existing modulefile hierarchy
   * continued use of existing modulefiles (using [Environment Modules](http://modules.sourceforge.net/)
   * pre-existing POSIX group of all users expected to execute builds
   * user base that is highly varied with regard to Unix knowledge - *keeping things simple encourages more widespread use*

Second, some notes about paths:

  * easybuild was bootstrapped into /app/easybuild
  * we created /app/easybuild/etc to hold additional centralized configuration files
  * we created /app/easybuild/fh_easyconfigs to hold our custom easyconfig files while we are developing them

# Step-By-Step Easybuild installation
  1. bootstrap easybuild - follow the [excellent documentation](https://easybuild.readthedocs.org/en/latest/Installation.html#bootstrapping-easybuild)
   Since we have an existing location for software packages, I chose that for the EasyBuild bootstrap (/app/easybuild in our case)
  2. Set EasyBuild variables
   Easybuild is very consistent in how it can be configured. A configuration file, command-line parameters, or environment variables are all recognized in the same and consistent way. Since Easybuild uses modules, I decided to set environment variables there. I made the following changes to files:

    * in the easybuild modulefile (/app/easybuild/modules/all/EasyBuild/2.3.0 at this time), I added the following:

```Tcl
setenv EASYBUILD_SOURCEPATH /app/easybuild/sources
setenv EASYBUILD_BUILDPATH /app/easybuild/build
setenv EASYBUILD_INSTALLPATH_SOFTWARE /app/easybuild/software
setenv EASYBUILD_INSTALLPATH_MODULES /app/easybuild/modules
setenv EASYBUILD_REPOSITORYPATH /app/easybuild/ebfiles_repo
setenv EASYBUILD_LOGFILE_FORMAT "/app/easybuild/logs,easybuild-%(name)s-%(version)s-%(date)s.%(time)s.log"
# keep group writable bit
setenv EASYBUILD_GROUP_WRITABLE_INSTALLDIR 1
# set umask to preserve group write permissions on modulefiles
setenv EASYBUILD_UMASK 002
# create module dependencies to recursively unload
setenv EASYBUILD_RECURSIVE_MODULE_UNLOAD 1
# add our normal modulefile footer
setenv EASYBUILD_MODULES_FOOTER /app/easybuild/etc/fredhutch_modulefile_footer
# add our own easyconfig directory to robot paths
setenv EASYBUILD_ROBOT_PATHS :/app/easybuild/fh_easyconfigs
# Our licenses
setenv LM_LICENSE_FILE /app/easybuild/etc/licenses/intel.lic
```

    * in /app/easybuild/etc/fredhutch_modulefile_footer, we added these lines to write module loads to syslog for syslog-driven metrics of software use:

```Tcl
# enable logging to syslog
set curMod [module-info name]
if { [module-info mode load] } {
system "logger \$USER module load $curMod "
}
```

# Step-By-Step Build a package

# Step-By-Step Adapt an EasyConfig
