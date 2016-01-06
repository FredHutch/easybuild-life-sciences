# easybuild-life-sciences
Implementation and use of EasyBuild at FredHutch

---

# Overview
- We identified EasyBuild as piece of software that could help us manage software package builds
- We use Environment Modules today, and more importantly, our users use modules
- We have a small group of people who build software packages (admins, not users)

---

# Goals
Before and during implementation, we kept the following goals in mind:

- software packages will be reproducable
- modules will be easily loaded by user in interactive sessions and in scripts
- default versions of software packages will be easy to manage (ex: R-3.2.3 may be most recent, but `module load R` will load R-3.2.1)
- packages will be built by any member of a given POSIX group
- new packages will be easily implemented (new versions and software packages without existing easyconfigs)

---

# TODO
Some of our goals were not met in the initial implementation. Mostly due to unimplemented features in Environment Modules and/or EasyBuild itself.

- New package implementation or **new easyconfigs made easy**
- - versions can be changed with `--try` options, but a new easyconfig per pkg/ver/compiler balloons quickly
- - Easybuild may be moving away from this model to more generic easyconfigs, but this may affect reproducibility

- Default Module Version Management
- - Modules picks the 'newest' version if not specified with `module load` or in `.version`
- - I find no option in Easybuild to manage `.version` files

---

# Prerequisites

Clearly, some software is required to compile software. We started with a base of Ubuntu 14.04 LTS along with the `build-essentials` meta package. Also, an implementation of Modules is required. We were already using Environment Modules, but if you do not have a Module suite in use, please check out [Lmod](https://www.tacc.utexas.edu/research-development/tacc-projects/lmod) - it is under current development and offers a different and more advanced feature set.

---

# Our environment

- read-only NFS mount on all systems mounted at `/app`
- hand-built software packages
- hand-managed modulefile hierarchy
- pre-existing POSIX group of all users expected to execute builds
- user base that is highly varied with regard to Unix knowledge - *keeping things simple encourages more widespread use*

---

# Bootstrap

  * easybuild was bootstrapped into `/app/easybuild`
  * we created `/app/easybuild/etc` to hold additional centralized configuration files
  * we created `/app/easybuild/fh_easyconfigs` to hold our custom easyconfig files while we are developing them

---

# Step-By-Step Easybuild Bootstrap

## Step One - RTFM

- Follow the [Fine Manual](https://easybuild.readthedocs.org/en/latest/Installation.html#bootstrapping-easybuild)

---

# Step-By-Step Easybuild Bootstrap

## Step Two - Environment

EasyBuild configuration

- consistent across methods:
- - config file(s)
- - environment variables
- - command-line parameters
- in that order
- we picked environment variables as it fits with Modules

---

# Modulefile Templating

## Paths and Logs

In the easybuild modulefile, I added the following:

!Tcl
set ebDir "/app/easybuild"
setenv EASYBUILD_SOURCEPATH "$ebDir/sources"
setenv EASYBUILD_BUILDPATH "$ebDir/build"
setenv EASYBUILD_INSTALLPATH_SOFTWARE "$ebDir/software"
setenv EASYBUILD_INSTALLPATH_MODULES "$ebDir/modules"
setenv EASYBUILD_REPOSITORYPATH "$ebDir/ebfiles_repo"
setenv EASYBUILD_LOGFILE_FORMAT "$ebDir/logs,easybuild-%(name)s-%(version)s-%(date)s.%(time)s.log"

---

# Modulefile Templating 2

## Easybuild parameters

!Tcl
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

---

# Easybuild Modulefile

- Our Easybuild config loads with `module load Easybuild/2.3.0` everytime for everyone
- With those modulefile additions saved to a separate file, they can be automatically included in future builds of EasyBuild itself
- Custom items can be added to all modulefiles produced by Easybuild, like logging module use:

!Tcl
# enable logging to syslog
set curMod [module-info name]
if { [module-info mode load] } {
system "logger \$USER module load $curMod "
}

---

 # EasyBuilt

```
$ module use /app/easybuild/modules/all   # adds this path to $MODULEPATH
$ module load Easybuild/2.3.0             # you should use the version you just bootstrapped - it should also tab out
$ eb --version
This is EasyBuild 2.3.0 (framework: 2.3.0, easyblocks: 2.3.0) on host rhino-d.
```

---

# Step-By-Step Build a package

Once you have EasyBuild bootstrapped, you can search for and build a package:

Begin by searching:

---

!bash
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

---

# Found PCRE esayconfigs!

We found 9 different easyconfigs for PCRE. Let's build this one: `PCRE-8.36-foss-2015a.eb`

You probably figured out that `8.36` is the version of PCRE we will build, but what is `foss`?

---

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

---

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
