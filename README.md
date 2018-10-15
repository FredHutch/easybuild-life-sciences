# EasyBuild at FredHutch

---

## Overview
- FredHutch Scientific Computing uses Easybuild to provide 100s of OSS packages to our Scientists
- Scientists can load multiple versions of any software via Environemnt modules (LMOD)
- All software is built to offer high reproducibility, it can be rebuilt exactly even 10 years from now 

---

## New package requests (including Python and R libraries/modules)
Please open an issue against this repo to request new softwares!

---

## Quickstart 

please follow these simple steps:

- use a system with at least 8GB RAM, 8GB Disk and 4 cores
- create a useraccount you using for builds, for example 'eb' :

```
sudo adduser --disabled-password --gecos "" eb
sudo sh -c "echo 'eb ALL=(ALL:ALL) NOPASSWD:ALL' > /etc/sudoers.d/zz_eb"

```

- and then simply launch the bootstrap process:

```
curl -s https://raw.githubusercontent.com/FredHutch/easybuild-life-sciences/master/easybuild_bootstrap.sh | bash
```
- after easybuild is installed simply log  out and login again and as an example (installing R) execute this

```
module load EasyBuild
eb R-3.3.1-foss-2016b.eb --robot
```


---

## Presentation
- This readme is also a [presentation](http://fredhutch.github.io/easybuild-life-sciences)

---

## Goals
Before and during implementation, we kept the following goals in mind:

- software packages will be reproducable
- modules will be easily loaded by user in interactive sessions and in scripts
- default versions of software packages will be easy to manage (ex: R-3.2.3 may be most recent, but `module load R` will load R-3.2.1)
- packages will be built by any member of a given POSIX group
- new packages will be easily implemented (new versions and software packages without existing easyconfigs)

---

## Prerequisites

You need money to make money, and you need software to build software.

- Ubuntu 14.04
- `build-essentials`
- an implementation of Modules - we use Environment Modules, but check out [Lmod](https://www.tacc.utexas.edu/research-development/tacc-projects/lmod)

---

## Our environment

- read-only NFS mount on all systems mounted at `/app`
- hand-built software packages
- hand-managed modulefile hierarchy
- pre-existing POSIX group of all users expected to execute builds
- user base that is highly varied with regard to Unix knowledge - *keeping things simple encourages more widespread use*

---

## Bootstrap

- easybuild was bootstrapped into `/app/easybuild`
- we created `/app/easybuild/etc` to hold additional centralized configuration files
- we created `/app/easybuild/fh_easyconfigs` to hold our custom easyconfig files while we are developing them

---

## Bootstrap - Step One - RTFM

- Follow the [Fine Manual](https://easybuild.readthedocs.org/en/latest/Installation.html#bootstrapping-easybuild)

---

## Bootstrap - Step Two - Environment

EasyBuild configuration

Configuration is consistent across methods:
- config file(s)
- environment variables
- command-line parameters

Easybuild applies them in that order (meaning command-line overrides everything)

Since we use Modules, it made sense to use Environment Variables in our case

---

## Bootstrap - Paths and Logs

In the easybuild modulefile, I added the following:

    !Tcl
    set ebDir "/app/easybuild"
    setenv EASYBUILD_SOURCEPATH "$ebDir/sources"
    setenv EASYBUILD_BUILDPATH "$ebDir/build"
    setenv EASYBUILD_INSTALLPATH_SOFTWARE "$ebDir/software"
    setenv EASYBUILD_INSTALLPATH_MODULES "$ebDir/modules"
    setenv EASYBUILD_REPOSITORYPATH "$ebDir/ebfiles_repo"
    setenv EASYBUILD_LOGFILE_FORMAT "$ebDir/logs,easybuild-%(name)s-%(version)s-%(date)s.%(time)s.log"

The modulefile is a tcl snippet and this sets environment variables for us.

---

## Bootstrap Easybuild Parameters

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

These are more complex, and will be documented soon.

---

## Bootstrap - Ownership and Permissions

There are a number of manual steps that were performed that can best be described as messy, and also perhaps make up the bulk of the useful information here.

Since we decided to have building be performed by members of a POSIX group, and we want produced software and modules centrally located for the use of everyone, we have to tell Easybuild how to do that.

- `GROUP_WRITABLE_INSTALLDIR` - this lets our group write to easybuild-created directories
- `UMASK 002` - this tells Easybuild to set this umask when writing directories/files

Of course, some manual adjusting was needed:

- `chgrp -R <build group> /app/easybuild/*` - change to our build group from the default group of the installer account
- `chmod -R g+s /app/easybuild/*` - setgid bit for dirs created during bootstrap

---

## Bootstrap - Easybuild Parameters

- `RECURSIVE_MODULE_UNLOAD 1` - this causes Easybuild to create modulefiles that will auto-load *and* auto-unload dependent modules
- `ROBOT_PATHS <path>` - adds our local easyconfig development dir to robot paths
- `LM_LICENSE_FILE <file>` - where our licenses are located

## Bootstrap - Modulefile Manipulating

- Now, Easybuild configured for us loads with `module load Easybuild/2.3.0` everytime for everyone
- `MODULES_FOOTER` - code to include in every modulefile created by Easybuild

Ex:

    !Tcl
    set curMod [module-info name]
    if { [module-info mode load] } {
        system "logger \$USER module load $curMod "
    }

---
## User Documentation

[Easybuild packages available at the Fred Hutch](https://fredhutch.github.io/easybuild-life-sciences/)

Github pages documentation is maintained for HPC users at the Fred Hutch.
 The GitHub pages have user-facing documentation about how to use modules,
 posts of recent package builds, and an inventory of packages. Users can access
 an inventory of available of 'bio' modules at the Hutch. Detailed documentation
 is given for R and Python modules which list all the versions available and
 the libraries they were built with. Links to "homepage" are available for
 users to access additional documentation about the software. 

### Maintaining Github Pages 
The content of the GitHub pages is generated by scripts in this repository.
 Inventory, Posts and detailed package information about R and Python are
 created by different scripts. The scripts need to be run by hand after
 packages have been built. After the documentation and package have been
 built the commit your changes and push to github. 

What needs to documented?  The GitHub pages are intended for SciComp users at
 the Hutch. Users are not interested in knowing about the latest version of
 libXext or lib. Only update documentation for Bio packages, R and Python.

It is assumed the maintainers are working from a clone of this repository. The
 scripts have to be run from within the directory tree of the repository.  

#### Generate a post to document a new module
The `create_post.sh` script takes a easybuild module name as an argument. The
 easyconfig is assumed to be in the fh_easyconfig unless a full path is
 specified. 
```
../scripts/create_post.sh beagle-lib-3.0.2-foss-2018b.eb
CFGS1=/app/easybuild/software/EasyBuild/3.7.0/lib/python2.7/site-packages/easybuild_easyconfigs-3.7.0-py2.7.egg/easybuild/easyconfigs
./create_posts.sh $CFGS1/p/PHASE/PHASE-2.1.1.eb
```

#### Generate new software inventory
The `create_module_list.sh` uses module spider to create a list of installed
 modules.  `create_module_list.sh` checks OS distributions and creates seperate
 output for Ubuntu 14.04 and 16.04.
```
~/scripts/create_module_list.sh
```

#### Generate module list for R and python
The script `easy_annotate.py` is used to create a Markdown page containing all
 the modules for an R or Python easyconfig. `easy_annotate.py` requires
 greater than Python-2.7.12. The output is written with .md file extention
 to the local directory.  The output has to be moved to the
 ~/easybuild-life-sciences/docs/[_R | _Python] directory.

#### Create software depenancy graph
Depenancy graphs require EasyBuid and Graphviz.  Run easybuild from the module
 repository with robot and 'dot' in the path.

```
module load EasyBuild
cd ~/easybuild-life-sciences/fh_easyconfigs
eb --dep-graph=texlive-20180531.dot ../fh_easyconfigs/texlive-20180531-foss-2016b.eb --robot .
module load Graphviz
dot -Tsvg texlive-20180531.dot >texlive-20180531.html
mv  texlive-20180531.html ~/easybuild-lif-sciences/docs/_Dot
mv texlive-20180531.dot   ~/easybuild-lif-sciences/docs/_Dot
```

---

## EasyBuilt

To use:

- Add the Easybuild modules directory to your MODULEPATH environment variable:

`$ module use /app/easybuild/modules/all`

- Load the EasyBuild module (it should tab out, these are just files - use this to find newest ver '...EasyBuild/<tab>'):

`$ module load EasyBuild/2.3.0`

- Did it work?

`$ eb --version  `
`This is EasyBuild 2.3.0 (framework: 2.3.0, easyblocks: 2.3.0) on host rhino-d.  `

*Note - you should always use the newest version of EasyBuild that has been built as easyconfigs are distributed with EB.
---

## Step-By-Step Build a package

Once you have EasyBuild bootstrapped, you can search for and build a package:

Begin by searching:

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

## Found easyconfigs!

We found 9 different easyconfigs for PCRE. Let's build this one:

`PCRE-8.36-foss-2015a.eb`

You probably figured out that `8.36` is the version of PCRE we will build, but what is `foss`?

That is the Easybuild toolchain for this easyconfig. You can get a list of toolchains with:

`eb --list-toolchains`

I prefer to just browse the [repo](https://github.com/hpcugent/easybuild-easyconfigs) - toolchains are just another easyconfig to Easybuild.

---

## Dry Run

Once we have decided what to build, you can do a dry-run like this:

    $ eb PCRE-8.36-foss-2015a.eb --robot --dry-run
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

---

## Dependencies

By giving the `-r` flag to Easybuild:

- Easybuild automatically spidered all dependencies, including the specified toolchain
- Already built packages are indicated with an `X`
- In this case, only PCRE will be built

---

## Build

And finally, you can remove the '-D' and build the software:

    $ eb PCRE-8.36-foss-2015a.eb --robot --force
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

---

## Step-By-Step EasyConfigs

Again, [RTFM](http://easybuild.readthedocs.org/en/latest/Writing_easyconfig_files.html) - it is very good!

There are two reason rou might want to modify or build an easyconfig file:

- Update versions (of software, toolchain, or dependencies)
- No easyconfig exists

I'll demonstrate creating a new easyconfig as the procedure is the same, only generally easier for version updates.

---

## Create an easyconfig file

Easybuild logic is contained in easyblocks - these are what execute the build. You can get a list of easyblocks with: `eb --list-easyblocks`.

There are a number of mandatory parameters for each easyblock, which can be displayed with:

`eb -a -e <easyblock>`

For this explanation, we will use the `ConfigureMake` easyblock, which should be familiar to anyone who has manually built software: `./configure && make && make install`.

The naming convention is typically `<name>-<version>-<toolchain name>-<toolchain version>.eb`.

---

## ConfigureMake Mandatory Parameters

This is a skeleton ConfigureMake easyconfig with all mandatory parameters:

    !python
    easyblock = 'ConfigureMake'
    name = 
    version = 
    toolchain = 
    description = 
    homepage = 
    docurls = 
    software_license = 
    software_license_urls = 

Except for easyblock, these will all default to `None` if not supplied in the file (so I guess they are not really mandatory, huh?)

---

## Parameters: name, version

## `name`

This name is the name of the software package, will be the name of the modulefile, and will be in the path of the software install directory. It is sometimes referenced later in the easyconfig file.

Ex:

    !python
    name = 'zlib'

## `version`

This is the version of the software to build. It is referenced later in the easyconfig file.

Ex:

    !python
    version = '1.2.8'

---

## Parameters: toolchain

## `toolchain`

This is the toolchain (compilers, supplemental libraries, etc.) that easybuild will use to build the software. It must be specified in an existing easyconfig (though does not need to be pre-built - easybuild will take care of building it).

Ex:

    !python
    toolchain = {'name': 'foss', 'version': '2015b'}

This is a python dict specifying the name and version of the toolchain.

---

## Parameters: description, homepage

## `description`

This is a generally free-form description that will appear as metadata in the modulefile, and therefore be availabe to users through the `module` command.

Ex:

    !python
    description = """zlib is designed to be a free, general-purpose, legally
                     unencumbered -- that is, not covered by any patents --
                     lossless data-compression library for use on virtually any
                     computer hardware and operating system."""
## `homepage`

This is a URL also included in modulefile metadata. It should be the homepage of the software.

Ex:

    !python
    homepage = 'http://www.zlib.net/'

---

## Parameters: source, source_urls

We will need to specify a few more parameters for easybuild to handle things correctly:

    !python
    sources = 
    source_urls = 

These will specify where easybuild should find the sourcecode for the software package. There are some shortcuts:

- Use `[SOURCELOWER_TAR_GZ]` to produce `<name>-<version>.tar.gz`
- Use `%{version}s` and `%{name}s` to use `name` and `version` from the easyconfig

Ex:

    !python
    sources = [SOURCELOWER_TAR_GZ]
    source_urls = ['http://sourceforge.net/projects/libpng/files/zlib/%(version)s']

These are python lists.

---

## Build it!

That should be sufficient to build a basic package. Let's see what a failure looks like.

Here is my perfect easyconfig for rsync:

    !python
    easyblock = 'ConfigureMake'
    name = 'rsync'
    version = '3.1.2'
    toolchain = {'name': 'foss', 'version': '2015b'}
    description = """rsync is an open source utility that provides fast incremental file transfer"""
    homepage = 'https://rsync.samba.org'
    sources = [SOURCELOWER_TAR_GZ]
    source_urls = ['https://download.samba.org/pub/rsync/src/rsync-3.1.2.tar.gz']

I save this as `rsync-3.1.2-foss-2015b.eb` and it should build!

---

## Build Example

    $  eb rsync-3.1.2-foss-2015b.eb
    == temporary log file in case of crash /tmp/eb-j_sVge/easybuild-cY3SFZ.log
    == processing EasyBuild easyconfig /app/easybuild/fh_easyconfigs/rsync-3.1.2-foss-2015b.eb
    == building and installing rsync/3.1.2-foss-2015b...
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
    == FAILED: Installation ended unsuccessfully (build directory: /app/easybuild/build/rsync/3.1.2/foss-2015b): build failed (first 300 chars): Sanity check failed: no dir of ('lib', 'lib64') in /app/easybuild/software/rsync/3.1.2-foss-2015b
    == Results of the build can be found in the log file /tmp/eb-j_sVge/easybuild-rsync-3.1.2-20160107.110632.Lepgl.log
    ERROR: Build of /app/easybuild/fh_easyconfigs/rsync-3.1.2-foss-2015b.eb failed (err: "build failed (first 300 chars): Sanity check failed: no dir of ('lib', 'lib64') in /app/easybuild/software/rsync/3.1.2-foss-2015b")

---

## Troubleshooting

We can look into the logfile mentioned (`Results of the build can be found in the log file /tmp/eb-j_sVge/easybuild-rsync-3.1.2-20160107.110632.Lepgl.log`') but in this case, the error is shown:

`Sanity check failed: no dir of ('lib', 'lib64') in /app/easybuild/software/rsync/3.1.2-foss-2015b`

And a quick search of `sanity check` in the Easybuild docs reveals that by default `bin` and `lib` or `lib64` must not be empty after install. Rsync builds no `lib` directory, so we add the following to the easyconfig file:

`sanity_check_paths = {'dirs': ['bin','share'], 'files': ['bin/rsync']}`

And now it builds (trust me).

**An interesting note** I expected this build to fail as version `3.1.2` of rsync is distributed in `rsync-3.1.2.tar.gz` but is not actually compressed, only a tarball. Easybuild built it anyway!
=======
## ls2_python3


Please look at [ls2](https://github.com/FredHutch/ls2) for details on how to build these Dockerfiles and how to use them to deploy the same software to a local archive.

This container adds: Python3

>>>>>>> 34c2607c52ea2576a50a407f68333ff624148831
