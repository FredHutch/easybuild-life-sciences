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

    * in the easybuild modulefile (/app/easybuild/modules/all/EasyBuild/2.3.0 at this time), I added the following (and also saved this separately to a file I can use during easybuilding Easybuild to automatically include this in the modulefile using the `--modules-footer` parameter):

```Tcl
set ebDir "/app/easybuild"
setenv EASYBUILD_SOURCEPATH "$ebDir/sources"
setenv EASYBUILD_BUILDPATH "$ebDir/build"
setenv EASYBUILD_INSTALLPATH_SOFTWARE "$ebDir/software"
setenv EASYBUILD_INSTALLPATH_MODULES "$ebDir/modules"
setenv EASYBUILD_REPOSITORYPATH "$ebDir/ebfiles_repo"
setenv EASYBUILD_LOGFILE_FORMAT "$ebDir/logs,easybuild-%(name)s-%(version)s-%(date)s.%(time)s.log"
# keep group writable bit
setenv EASYBUILD_GROUP_WRITABLE_INSTALLDIR 1
# set umask to preserve group write permissions on modulefiles
setenv EASYBUILD_UMASK 002
# create module dependencies to recursively unload
setenv EASYBUILD_RECURSIVE_MODULE_UNLOAD 1
# add our normal modulefile footer
setenv EASYBUILD_MODULES_FOOTER "$ebDir/etc/fredhutch_modulefile_footer"
# add our own easyconfig directory to robot paths
setenv EASYBUILD_ROBOT_PATHS ":$ebDir/fh_easyconfigs"
# Our licenses
setenv LM_LICENSE_FILE "$ebDir/etc/licenses/intel.lic"
```

    * in /app/easybuild/etc/fredhutch_modulefile_footer, we added these lines to write module loads to syslog for syslog-driven metrics of software use:

```Tcl
# enable logging to syslog
set curMod [module-info name]
if { [module-info mode load] } {
system "logger \$USER module load $curMod "
}
```

   3. At this point, you should be able to load the EasyBuild module:

```
$ module use /app/easybuild/modules/all   # adds this path to MODULEPATH
$ module load Easybuild/2.3.0             # you should use the version you just bootstrapped - it should also tab out
```

# Step-By-Step Build a package

Once you have EasyBuild bootstrapped, you can search for and build a package:

Begin by searching:

```
$ eb -S PCRE
== temporary log file in case of crash /tmp/eb-lz7d_6/easybuild-dKc03x.log
== Searching (case-insensitive) for 'PCRE' in /app/easybuild/software/EasyBuild/2.3.0/lib/python2.7/site-packages/easybuild_easyconfigs-2.3.0-py2.7.egg/easybuild/easyconfigs 
== Searching (case-insensitive) for 'PCRE' in /app/easybuild/fh_easyconfigs 
CFGS1=/app/easybuild/software/EasyBuild/2.3.0/lib/python2.7/site-packages/easybuild_easyconfigs-2.3.0-py2.7.egg/easybuild/easyconfigs/p/PCRE
 * $CFGS1/PCRE-8.12-goalf-1.1.0-no-OFED.eb
 * $CFGS1/PCRE-8.12-goolf-1.4.10.eb
 * $CFGS1/PCRE-8.12-ictce-4.0.6.eb
 * $CFGS1/PCRE-8.12-ictce-5.3.0.eb
 * $CFGS1/PCRE-8.12-ictce-5.5.0.eb
 * $CFGS1/PCRE-8.35-intel-2014b.eb
 * $CFGS1/PCRE-8.36-foss-2015a.eb
 * $CFGS1/PCRE-8.36-intel-2015a.eb
 * $CFGS1/PCRE-8.37-intel-2015a.eb
== Tmporary log file(s) /tmp/eb-lz7d_6/easybuild-dKc03x.log* have been removed.
== Temporary directory /tmp/eb-lz7d_6 has been removed.
```

Here we have found 9 different easyconfigs for PCRE. The keyword after the version is the toolchain. You should research toolchains at some point, but we chose to focus on two toolchain families: `intel` and `foss`. One is a closed-source optimized compiler for Intel CPUs and the other is an open-source chain using free open source software.

Once we have decided what to build, you can do a dry-run like this:

```
$ eb -r -D PCRE-8.36-foss-2015a.eb
== temporary log file in case of crash /tmp/eb-08QTaF/easybuild-r5D8gf.log
Dry run: printing build status of easyconfigs and dependencies
CFGS=/app/easybuild/software/EasyBuild/2.3.0/lib/python2.7/site-packages/easybuild_easyconfigs-2.3.0-py2.7.egg/easybuild/easyconfigs
 * [x] $CFGS/g/GCC/GCC-4.9.2.eb (module: GCC/4.9.2)
 * [x] $CFGS/o/OpenBLAS/OpenBLAS-0.2.13-GCC-4.9.2-LAPACK-3.5.0.eb (module: OpenBLAS/0.2.13-GCC-4.9.2-LAPACK-3.5.0)
 * [x] $CFGS/l/libtool/libtool-2.4.2-GCC-4.9.2.eb (module: libtool/2.4.2-GCC-4.9.2)
 * [x] $CFGS/m/M4/M4-1.4.17-GCC-4.9.2.eb (module: M4/1.4.17-GCC-4.9.2)
 * [x] $CFGS/a/Autoconf/Autoconf-2.69-GCC-4.9.2.eb (module: Autoconf/2.69-GCC-4.9.2)
 * [x] $CFGS/a/Automake/Automake-1.15-GCC-4.9.2.eb (module: Automake/1.15-GCC-4.9.2)
 * [x] $CFGS/n/numactl/numactl-2.0.10-GCC-4.9.2.eb (module: numactl/2.0.10-GCC-4.9.2)
 * [x] $CFGS/h/hwloc/hwloc-1.10.0-GCC-4.9.2.eb (module: hwloc/1.10.0-GCC-4.9.2)
 * [x] $CFGS/o/OpenMPI/OpenMPI-1.8.4-GCC-4.9.2.eb (module: OpenMPI/1.8.4-GCC-4.9.2)
 * [x] $CFGS/g/gompi/gompi-2015a.eb (module: gompi/2015a)
 * [x] $CFGS/f/FFTW/FFTW-3.3.4-gompi-2015a.eb (module: FFTW/3.3.4-gompi-2015a)
 * [x] $CFGS/s/ScaLAPACK/ScaLAPACK-2.0.2-gompi-2015a-OpenBLAS-0.2.13-LAPACK-3.5.0.eb (module: ScaLAPACK/2.0.2-gompi-2015a-OpenBLAS-0.2.13-LAPACK-3.5.0)
 * [x] $CFGS/f/foss/foss-2015a.eb (module: foss/2015a)
 * [ ] $CFGS/p/PCRE/PCRE-8.36-foss-2015a.eb (module: PCRE/8.36-foss-2015a)
== Tmporary log file(s) /tmp/eb-08QTaF/easybuild-r5D8gf.log* have been removed.
== Temporary directory /tmp/eb-08QTaF has been removed.
```

In this case, EasyBuild was given the '-r' robot flag so it automatically included dependencies it was able to meet with existing easyconfigs, and it was told to do a dry-run with '-D'. You will note the all of the packages except PCRE itself have 'X' noting they are already built and installed.

And finally, you can remove the '-D' and build the software:

```
$ eb -r -f PCRE-8.36-foss-2015a.eb
== temporary log file in case of crash /tmp/eb-1TnpU8/easybuild-3J4ttj.log
== resolving dependencies ...
== processing EasyBuild easyconfig /app/easybuild/software/EasyBuild/2.3.0/lib/python2.7/site-packages/easybuild_easyconfigs-2.3.0-py2.7.egg/easybuild/easyconfigs/p/PCRE/PCRE-8.36-foss-2015a.eb
== building and installing PCRE/8.36-foss-2015a...
== fetching files...
== creating build dir, resetting environment...
== unpacking...
== patching...
== preparing...
== configuring...
== building...
== testing...
== installing...
== taking care of extensions...
== postprocessing...
== sanity checking...
== cleaning up...
== creating module...
== permissions...
== packaging...
== COMPLETED: Installation ended successfully
== Results of the build can be found in the log file /app/easybuild/logs/easybuild-PCRE-8.36-20160104.164159.log
== Build succeeded for 1 out of 1
== Tmporary log file(s) /tmp/eb-1TnpU8/easybuild-3J4ttj.log* have been removed.
== Temporary directory /tmp/eb-1TnpU8 has been removed.
```
# Step-By-Step Adapt an EasyConfig
