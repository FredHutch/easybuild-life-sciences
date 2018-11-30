---
layout: post
title: foss-2016b
---

Foss 2016b is currently in production on all HPC resources.  The packages are built to s
support Ubunutu 14.04.

### Package List
 * [*GCCcore-5.4.0*](http://gcc.gnu.org/) The GNU Compiler Collection includes front ends for C, C++, Objective-C, Fortran, Java, and Ada,
 as well as libraries for these languages (libstdc++, libgcj,...).

 * [binutils-2.26](http://directory.fsf.org/project/binutils/)binutils: GNU binary utilities
 * [numactl-2.0.11](http://oss.sgi.com/projects/libnuma/) The numactl program allows you to run your application program on specific cpu's and memory nodes.
 It does this by supplying a NUMA memory policy to the operating system before running your program.
 The libnuma library provides convenient ways for you to add NUMA memory policies into your own program.'
 * [hwloc-1.11.3](http://www.open-mpi.org/projects/hwloc/) The Portable Hardware Locality (hwloc) software package provides a portable abstraction
 (across OS, versions, architectures, ...) of the hierarchical topology of modern architectures, including
 NUMA memory nodes, sockets, shared caches, cores and simultaneous multithreading. It also gathers various
 system attributes such as cache and memory information as well as the locality of I/O devices such as
 network interfaces, InfiniBand HCAs or GPUs. It primarily aims at helping applications with gathering
 information about modern computing hardware so as to exploit it accordingly and efficiently.
 * [*OpenMPI-1.10.3*](http://www.open-mpi.org/) The Open MPI Project is an open source MPI-2 implementation.
 * [OpenBLAS-0.2.18](http://xianyi.github.com/OpenBLAS/) OpenBLAS is an optimized BLAS library based on GotoBLAS2 1.13 BSD version.
 * [FFTW-3.3.4](http://www.fftw.org) FFTW is a C subroutine library for computing the discrete Fourier transform (DFT)
 in one or more dimensions, of arbitrary input size, and of both real and complex data.
 * [*ScaLAPACK-2.0.2 OpenBLAS-0.2.18 LAPACK-3.6.1*](http://www.netlib.org/scalapack/) The ScaLAPACK (or Scalable LAPACK) library includes a subset of LAPACK routines redesigned for distributed memory MIMD parallel computers.
 * [*foss-2016b*](https://github.com/easybuilders/easybuild-easyconfigs/blob/afb34474f4d35d4676b4e5d1cbd94c8e0c978afe/easybuild/easyconfigs/f/foss/foss-2016b.eb) Free Open Source Software. GNU Compiler Collection (GCC) based compiler toolchain, including
 OpenMPI for MPI support, OpenBLAS (BLAS and LAPACK support), FFTW and ScaLAPACK.
