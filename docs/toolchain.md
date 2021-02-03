---
title: toolchains
permalink: /toolchains/
layout: single
date:  2020-02-01
sidebar:
  nav: "docs"
---

Most software packages are built from a compiler toolchain. A toolchain consists of one or more compilers, put
together with some libraries for using an MPI stack and commonly used math libraries. Well-know packages 
include BLAS/LAPACK APIs for linear algebra routines.

There are several different toolchains supported by EasyBuild. The Hutch uses the foss toolchain, which
is an acronym for Free Open Source Software. Foss is based on the GNU compiler tools.  
Toolchains are currated twice a year and named for the year and labeled a or b, ie: foss-2019b.
When loading multible modules for a workflow always load tools from the same toolchain. Do not mix and
match tools from different toolchains. j

Begining in the fall of 2020, new software will be built with the foss-2020b toolchain at the Hutch. The
Foss-2019b was the primary toolchain to support Ubuntu Bionic 18.04. And foss-2016b was used for Ubuntu 14.04.

## GNU Toolchain Versions

| Toolchain | GCC Version | Notes |
| ----------|-------------| ---------|
| [foss-2020b]({{ site.baseurl }}/foss-2020b/) | GCC 10.2.0 | Currently supported package development |
| foss-2020a | GCC 9.3.0 | Mostly skipped | 
| [foss-2019b]({{ site.baseurl }}/foss-2019b/) | GCC 8.3.0 | Primay tool chain for 2019 cluster update |
| foss-2019a | GCC 8.2.0 | skipped |
| [foss-2018b]({{ site.baseurl }}/foss-2018b/) | GCC 7.3.0 | Many module were built for the new cluster, but are outdated by 2019 |
| foss-2016b | GCC 5.4.0 | Old cluster |
|------------|-----------|-------------|


