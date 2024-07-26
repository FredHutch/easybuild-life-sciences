---
title: Scientific Software 
layout: single
permalink: /
toc: true
toc_label: "On This Page"
sidebar:
  nav: "docs"
---

Fred Hutch maintains open source scientific software for use with HPC resouces
at the center. This site provides an inventory of available software and
information about using scientific software.  The Hutch uses Easybuild 
to build and manage software. Easybuild is a software build and installation
framework for managing scientific software. 

## Life Science Software Inventory
Inventory of module class "bio" software modules.  This is a subset of all modules that are available. The full
list of modules encompesses low level system libraries, math libraries, programming languages, and visualization tools. 

 - [Life Science Software Inventory]({{ site.baseurl }}/bio-modules-18.04/)

## Programming Language Modules
Scientific Computing maintains custom builds for R and Python. The
custom R and Python modules contain hundreds of packages. The package
list is a compilation of user requests.  Click on the links
bellow to see a list of available builds with a list of modules.

 - [R Modules]({{ site.baseurl }}/r/)
 - [Python Modules]({{ site.baseurl }}/python/)
 - [Dependency Graphs]({{ site.baseurl }}/dot/)

### Scientific Software Environment
The user interface for using software is Modules.
The module command is used to instantiate a specific software package.
Easybuild and Modules provide the tools for using scientific software in a
repeatable and documented way.

### Software Packages
Easybuild modules are built from [toolchains]({{ site.baseurl }}/toolchains/).
 The toolchain provides the foundation for system level libraries and
math libraries. Easybuild software modules names contain the software name-version and the toolchain name
used to build the module. Custom reciepies used at the Hutch also have suffix <tt>fh</tt> to distigue them
from other published reciepeis. Once a reciepe is published to this repository it is not changed. If changes
have to be made to a package the suffix is versioned. (fh1, fh2, fh3 etc).

### Using Modules

Modules can be loaded, unloaded, listed and searched.  The
```module load``` command can be abrivated to ```ml```.

Load a Python module
```ml Python/3.7.4-foss-2019b-fh1```

List currently loaded modules
```ml list```

Unload all modules
```ml purge```

Show all available Python modules
```ml avail Python``` 

### Using Modules in Scripts
Place the following in your bash or sbatch scripts to load modules within
your scripts.

```
source /app/lmod/lmod/init/profile
module load R/4.0.0-foss-2019b
```

### Software Requests
Please send software requests to SciComp support.
Requests can be made for additional packages to be added for R,
Python and for new software.
Packaged modules of R and Python are made as soon as possible 
after new releases are announced.

### Module Usage Report
Most frequently loaded modules.  (Jan 1, 2024 -> July 1, 2024)

  | Module Name | Count |
  |-------------|-------|
  | R/4.4.0-gfbf-2023b | 515887 |
  | Boost/1.83.0-GCC-13.2.0 | 442110 |
  | Eigen/3.4.0-GCCcore-13.2.0 | 442038 |
  | R/4.3.3-gfbf-2023b | 410328 |
  | fhR/4.3.1-foss-2022b | 315485 |
  | Python/3.7.4-foss-2019b-fh1 | 290945 |
  | R/4.3.1-gfbf-2022b | 128351 |
  | GSL/2.7-GCCcore-12.2.0 | 83736 |
  | ImageMagick/7.1.0-53-GCCcore-12.2.0 | 83665 |
  | cuDNN/8.4.1.50-CUDA-11.7.0 | 83516 |
  | fhR/4.1.2-foss-2021b | 80743 |
  | fhR/4.0.4-foss-2020b | 73462 |
  | fhR/4.4.0-foss-2023b | 40147 |
  | fhR/4.3.3-foss-2023b | 35929 |
  | gdc-client/1.6.1-GCCcore-10.2.0 | 34472 |
  | R/4.3.0-foss-2022b | 26893 |
  | Bowtie2/2.4.2-GCC-10.2.0 | 24845 |
  | BEDTools/2.30.0-GCC-11.2.0 | 23939 |
  | R/4.2.0-foss-2021b | 23696 |
  | BCFtools/1.19-GCC-13.2.0 | 19256 |
  | GATK/4.4.0.0-GCCcore-12.2.0-Java-17 | 18423 |
  | R/4.2.2-foss-2022b | 17425 |
  | fhR/4.2.0-foss-2021b | 16071 |
  | MACS2/2.2.9.1-foss-2022b | 14896 |
  | PyTorch/1.11.0-foss-2021b-CUDA-11.4.1 | 14791 |
  | FastQC/0.11.9-Java-11 | 13491 |
  | torchvision/0.13.0-foss-2021b-CUDA-11.4.1 | 13433 |
  | Trinity/2.12.0-foss-2020b | 12633 |
  | IgBLAST/1.22.0-x64-linux | 12185 |
  | OpenSlide/3.4.1-GCCcore-12.3.0-largefiles | 11734 |
  | Java/17.0.6 | 8672 |
  | SAMtools/1.14-GCC-11.2.0 | 6238 |
  | Kent_tools/20201201-linux.x86_64 | 5732 |
  | BEDTools/2.30.0-GCC-12.2.0 | 5388 |
  | Python/3.8.2-GCCcore-9.3.0 | 5381 |
  | Bowtie2/2.4.4-GCC-11.2.0 | 4941 |
  | Miniconda3/4.9.2 | 4840 |
  | SAMtools/1.11-GCC-10.2.0 | 4519 |
  | BEDTools/2.29.2-GCC-9.3.0 | 4020 |
  | Anaconda3/2023.09-0 | 3821 |
  | cutadapt/4.1-GCCcore-11.2.0 | 3791 |
  | SAMtools/1.16.1-GCC-11.2.0 | 3583 |
  | Java/11.0.2 | 3542 |
  | deepTools/3.5.4.post1-gfbf-2022b | 3512 |
  | Apptainer/1.1.6 | 3478 |
  | SAMtools/1.10-GCCcore-8.3.0 | 2556 |
  | openslide-python/1.1.2-GCCcore-11.2.0 | 2527 |
  | picard/2.18.29-Java | 2358 |
  | fhPython/3.9.6-foss-2021b | 1994 |
  | picard/2.21.6-Java-11 | 1851 |
  | SAMtools/1.17-GCC-12.2.0 | 1733 |
  | STAR/2.7.7a-GCC-10.2.0 | 1711 |
  | GATK/4.2.6.1-GCCcore-11.2.0 | 1638 |
  | Anaconda3/2020.02 | 1621 |
  | GATK/4.1.8.1-GCCcore-8.3.0-Java-11 | 1466 |
  | SRA-Toolkit/3.0.0-ubuntu64 | 1458 |
  | scikit-learn/1.0.1-foss-2021b | 1404 |
  | deepTools/3.5.1-foss-2021b | 1331 |
  | picard/2.6.0-Java-11 | 1325 |
  | SAMtools/1.19.2-GCC-13.2.0 | 1316 |
  | PostgreSQL/10.3-foss-2019b | 1121 |
  | Java/21.0.2 | 1111 |
  | CUDA/12.1.1 | 1046 |
  | MACS2/2.2.6-foss-2019b-Python-3.7.4 | 1024 |
  | STAR/2.7.10b-GCC-12.2.0 | 1023 |
  | CellRanger/8.0.0 | 1004 |
  | PyTorch/1.12.1-foss-2022a-CUDA-11.7.0 | 1002 |
  | matplotlib/3.5.2-foss-2022a | 974 |
  | BLAST+/2.14.0-gompi-2022b | 894 |
  | Python/3.10.8-GCCcore-12.2.0 | 874 |
  | Nextflow/23.04.2 | 755 |
  | Singularity/3.5.3 | 750 |
  | picard/2.25.1-Java-11 | 727 |
  | Kraken2/2.1.3-gompi-2022b | 701 |
  | seqtk/1.3-GCC-10.2.0 | 655 |
  | GATK/4.2.5.0-GCCcore-11.2.0-Java-11 | 652 |
  | BWA/0.7.17-GCCcore-11.2.0 | 613 |
  | SlamDunk/0.4.3-foss-2021b | 606 |
