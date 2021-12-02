---
title: foss-2019b
date: 2019-08-01
---

Foss-2019b is the primary toolchain to support Ubuntu Bionic 18.04. Foss-2019b is based on GCC 8.3.0.

### Package List
 * [GCCcore-8.3.0](https://gcc.gnu.org/)The GNU Compiler Collection includes front ends for
 C, C++, Objective-C, Fortran, Java, and Ada, as well as libraries for these languages (libstdc++, libgcj,...).
 * [binutils-2.32](http://directory.fsf.org/project/binutils/) binutils: GNU binary utilities
 * [GCC-8.3.0](http://gcc.gnu.org/) The GNU Compiler Collection.
 * [numactl-2.0.12](http://oss.sgi.com/projects/libnuma/) The numactl program allows you to run your
   application program on specific cpus and memory nodes.
 It does this by supplying a NUMA memory policy to the operating system before running your program.
 The libnuma library provides convenient ways for you to add NUMA memory policies into your own program.
 * [hwloc-1.11.12](http://www.open-mpi.org/projects/hwloc/) The Portable Hardware Locality (hwloc) software package provides a portable abstraction
 (across OS, versions, architectures, ...) of the hierarchical topology of modern architectures, including
 NUMA memory nodes, sockets, shared caches, cores and simultaneous multithreading. It also gathers various
 system attributes such as cache and memory information as well as the locality of I/O devices such as
 network interfaces, InfiniBand HCAs or GPUs. It primarily aims at helping applications with gathering
 information about modern computing hardware so as to exploit it accordingly and efficiently.
 * [OpenMPI-3.1.4](http://www.open-mpi.org/) The Open MPI Project is an open source MPI-2 implementation.
 * [OpenBLAS-0.3.7](http://xianyi.github.com/OpenBLAS/) OpenBLAS is an optimized BLAS library based on GotoBLAS2 1.13 BSD version.
 including OpenMPI for MPI support.
 * [FFTW-3.3.8](http://www.fftw.org) FFTW is a C subroutine library for computing the discrete Fourier transform (DFT)
 in one or more dimensions, of arbitrary input size, and of both real and complex data.
 * [ScaLAPACK-2.0.2](http://www.netlib.org/scalapack/) The ScaLAPACK (or Scalable LAPACK) library includes a subset of LAPACK routines redesigned for distributed memory MIMD parallel computers.
 * [foss-2019b](https://raw.githubusercontent.com/easybuilders/easybuild-easyconfigs/master/easybuild/easyconfigs/f/foss/foss-2019b.eb)
 FOSS - Free Open Source Software. GNU Compiler Collection (GCC) based compiler toolchain, including
 OpenMPI for MPI support, OpenBLAS (BLAS and LAPACK support), FFTW and ScaLAPACK.
