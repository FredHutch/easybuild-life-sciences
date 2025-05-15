---
title: Bio Modules 24.04
layout: single
permalink: /sw_inventory/quokka-bio-modules-24.04/
created: 2025-05-15
toc: true
toc_label: "On This Page"
sidebar:
  nav: "docs"
---

 - [BEDTools/2.31.0-GCC-12.3.0](https://bedtools.readthedocs.io/)
BEDTools: a powerful toolset for genome arithmetic.
The BEDTools utilities allow one to address common genomics tasks such as finding feature overlaps and
computing coverage.
The utilities are largely based on four widely-used file formats: BED, GFF/GTF, VCF, and SAM/BAM.
 - [BamTools/2.5.2-GCC-12.3.0](https://github.com/pezmaster31/bamtools)
BamTools provides both a programmer's API and an end-user's toolkit for handling BAM files.
 - [Biopython/1.83-foss-2023a](https://www.biopython.org)
Biopython is a set of freely available tools for biological
 computation written in Python by an international team of developers. It is
 a distributed collaborative effort to develop Python libraries and
 applications which address the needs of current and future work in
 bioinformatics. 
 - [Eigen/3.4.0-GCCcore-13.3.0](https://eigen.tuxfamily.org)
Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers,
 and related algorithms.
 - [GEOS/3.12.0-GCC-12.3.0](https://trac.osgeo.org/geos)
GEOS (Geometry Engine
 - [GMP/6.2.1-GCCcore-12.2.0](https://gmplib.org/)
GMP is a free library for arbitrary precision arithmetic, operating on signed
 integers, rational numbers, and floating point numbers.

 - [MPC/1.3.1-GCCcore-12.3.0](http://www.multiprecision.org/)
Gnu Mpc is a C library for the arithmetic of
 complex numbers with arbitrarily high precision and correct
 rounding of the result. It extends the principles of the IEEE-754
 standard for fixed precision real floating point numbers to
 complex numbers, providing well-defined semantics for every
 operation. At the same time, speed of operation at high precision
 is a major design goal.
 - [MPFR/4.2.0-GCCcore-12.2.0](https://www.mpfr.org)
The MPFR library is a C library for multiple-precision floating-point
 computations with correct rounding.

 - [MUMmer/4.0.0rc1-GCCcore-12.3.0](https://mummer.sourceforge.net/)
MUMmer is a system for rapidly aligning entire genomes,
 whether in complete or draft form. AMOS makes use of it.

 - [Porechop/0.2.4-GCCcore-12.3.0](https://github.com/rrwick/Porechop)
[easyconfig](https://github.com/FredHutch/easybuild-life-sciences/blob/master/easyconfigs/p/Porechop/Porechop-0.2.4-GCCcore-12.3.0.eb)
Porechop is a tool for finding and removing adapters from Oxford Nanopore reads.
 Adapters on the ends of reads are trimmed off, and when a read has an adapter in its middle,
 it is treated as chimeric and chopped into separate reads. Porechop performs thorough alignments
 to effectively find adapters, even at low sequence identity
 - [Pysam/0.22.0-GCC-12.3.0](https://github.com/pysam-developers/pysam)
Pysam is a python module for reading and manipulating Samfiles.
 It's a lightweight wrapper of the samtools C-API. Pysam also includes an interface for tabix.
 - [Qhull/2020.2-GCCcore-12.2.0](http://www.qhull.org)
Qhull computes the convex hull, Delaunay triangulation, Voronoi diagram,
 halfspace intersection about a point, furthest-site Delaunay triangulation,
 and furthest-site Voronoi diagram. The source code runs in 2-d, 3-d, 4-d, and
 higher dimensions. Qhull implements the Quickhull algorithm for computing the
 convex hull.

 - [Seaborn/0.13.2-gfbf-2023a](https://seaborn.pydata.org/)
Seaborn is a Python visualization library based on matplotlib.
 It provides a high-level interface for drawing attractive statistical graphics. 
 - [Shapely/2.0.1-gfbf-2023a](https://github.com/Toblerity/Shapely)
Shapely is a BSD-licensed Python package for manipulation and analysis of planar geometric objects.
It is based on the widely deployed GEOS (the engine of PostGIS) and JTS (from which GEOS is ported) libraries.
 - [anndata/0.10.5.post1-foss-2023a](https://github.com/scverse/anndata)
anndata is a Python package for handling annotated data matrices in memory and on disk,
 positioned between pandas and xarray
 - [biom-format/2.1.15-foss-2023a](https://biom-format.org)
The BIOM file format (canonically pronounced biome) is designed to be
 a general-use format for representing biological sample by observation
 contingency tables. BIOM is a recognized standard for the Earth Microbiome
 Project and is a Genomics Standards Consortium supported project.

 - [gmpy2/2.1.5-GCC-12.3.0](https://github.com/aleaxit/gmpy)
GMP/MPIR, MPFR, and MPC interface to Python 2.6+ and 3.x
 - [libcerf/2.3-GCCcore-12.3.0](https://jugit.fz-juelich.de/mlz/libcerf)
libcerf is a self-contained numeric library that provides an efficient and
 accurate implementation of complex error functions, along with Dawson,
 Faddeeva, and Voigt functions.

 - [pybedtools/0.9.1-foss-2023a](https://daler.github.io/pybedtools)
pybedtools wraps and extends BEDTools and offers feature-level manipulations from within Python.
 - [scanpy/1.9.8-foss-2023a](https://scanpy.readthedocs.io/en/stable/)
[easyconfig](https://github.com/FredHutch/easybuild-life-sciences/blob/master/easyconfigs/s/scanpy/scanpy-1.9.8-foss-2023a.eb)
Scanpy is a scalable toolkit for analyzing single-cell gene expression data built
 jointly with anndata. It includes preprocessing, visualization, clustering, trajectory inference
 and differential expression testing. The Python-based implementation efficiently deals with
 datasets of more than one million cells.

 - [scikit-bio/0.6.0-foss-2023a](http://scikit-bio.org)
scikit-bio is an open-source, BSD-licensed Python 3 package providing data structures, algorithms
and educational resources for bioinformatics.
 - [statsmodels/0.14.1-gfbf-2023a](https://www.statsmodels.org/)
Statsmodels is a Python module that allows users to explore data, estimate statistical models,
and perform statistical tests.
 - [sympy/1.12-gfbf-2023a](https://sympy.org/)
SymPy is a Python library for symbolic mathematics. It aims to
 become a full-featured computer algebra system (CAS) while keeping the code as
 simple as possible in order to be comprehensible and easily extensible. SymPy
 is written entirely in Python and does not require any external libraries.
