# easybuild-life-sciences
Implementation and use of EasyBuild at FredHutch
================================================

# Overview
Over
----
This repo will reflect our decisions and strategy during implementation of EasyBuild to replace our manual build procedure

# Goals
Before and during implementation, we kept the following goals in mind:

   * packages will be reproducable
   * modules will be easily loaded by user in interactive sessions and in scripts
   * default versions of software packages will be easy to manage
   * packages will be built by members of a POSIX group
   * new packages will be easily implemented

# TODO
Some of our goals were not met in the initial implementation. Mostly due to unimplemented features in Environment Modules and/or EasyBuild itself.

New package implementation (i.e.: new easyconfigs) - Easybuild may be in transition with regard to easyconfig implementation. Initially the idea of reproducibility and separate easyconfig files went hand-in-hand. You build R version 3.2.1 with an Intel toolchain so you have an easyconfig detailing just that. At this time there are 738 software packages with easyconfigs in the public repo. If we figure only two toolchains and three versions for each package, we will need 4428 easyconfig files. The number of easyconfigs will become more difficult to manage, and certainly more difficult to choose. Easybuild supports several options in the "try" family like "--try-software-version=" and "--try-toolchain=" that can help re-build packages under a different toolchain or with a new version.

Default module selection - In many cases, modules does a good job of picking the most recent modulefile to use. However in several cases (R, PYthon, Intel toolchain, etc.) we do not want the most recent package to be the default. Easybuild (AFAIK) does not provide a mechanism for managing the contents of ".version" in modulefiles directories. Nor should it. However, a written procedure like this does not leave me with a good feeling toward ending up with correctly managed default versions effortlessly. Perhaps a wrapper script. Or an easybuild option "--make-default" that would modify the ".version" file in the resulting modulefile directory.

# Implementation Notes
First, some notes about 

# Step-By-Step Easybuild installation

# Step-By-Step Build a package

# Step-By-Step Adapt an EasyConfig
