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
Most frequently loaded modules.  (Oct 3, 2022 - Apr 3, 2023)

 
  | Name | Count | Percent |
  | -----|-------|--------| |
  | R | 1,431,523 | 39.443534 |
  | fhR | 1,253,561 | 34.540050 |
  | BCFtools | 549,348 | 15.136485 |
  | Java | 71,586 | 1.972448 |
  | picard | 36,661 | 1.010141 |
  | Python | 32,381 | 0.892211 |
  | SAMtools | 24,188 | 0.666465 |
  | BEDTools | 20,939 | 0.576944 |
  | deepTools | 16,995 | 0.468273 |
  | Bowtie2 | 16,891 | 0.465407 |
  | GATK | 10,888 | 0.300003 |
  | Kent_tools | 10,266 | 0.282865 |
  | Kraken2 | 9,188 | 0.253162 |
  | beagle-lib | 8,806 | 0.242637 |
  | FastQC | 8,787 | 0.242113 |
  | STAR | 8,286 | 0.228309 |
  | MUMmer | 7,857 | 0.216488 |
  | cutadapt | 6,010 | 0.165597 |
  | Anaconda3 | 5,627 | 0.155044 |
  | Trinity | 5,123 | 0.141157 |
  | PLINK2 | 4,997 | 0.137685 |
  | Singularity | 4,994 | 0.137602 |
  | IgBLAST | 4,393 | 0.121043 |
  | fhPython | 3,903 | 0.107541 |
  | Pandoc | 3,835 | 0.105668 |
  | Arriba | 3,152 | 0.086849 |
  | BWA | 3,003 | 0.082743 |
  | SRA-Toolkit | 2.468 | 0.068002 |
  | QUAST | 2.300 | 0.063373 |
  | RStudio-Server | 2.290 | 0.063098 |
  | Beast | 2,084 | 0.057420 |
  | tmux | 2,057 | 0.056676 |
  | BLAST+ | 1,775 | 0.048906 |
  | MACS2 | 1,612 | 0.044415 |
  | HTSeq | 1,569 | 0.043230 |
  | libGLU | 1,554 | 0.042817 |
  | Nextflow | 1,469 | 0.040475 |
  | TensorFlow | 1,262 | 0.034772 |
  | Emacs | 1,251 | 0.034469 |
  | CellRanger | 1,200 | 0.033063 |
  | strelka | 1,153 | 0.031768 |
  | fgbio | 1,121 | 0.030887 |
  | rstudio-server | 1,113 | 0.030666 |
  | Subread | 1,089 | 0.030005 |
  | Trimmomatic | 1,087 | 0.029950 |
  | CNVkit | 1,053 | 0.029013 |
  | scanpy | 1,044 | 0.028765 |
  | XGBoost | 1,004 | 0.027663 |
  | PLINK | 910 | 0.025073 |
  | foss | 900 | 0.024798 |
  | kallisto | 894 | 0.024632 |
  | git | 874 | 0.024081 |
  | fhDev | 805 | 0.022180 |
  | X11 | 793 | 0.021849 |
  | GLPK | 750 | 0.020665 |
  | OptiType | 736 | 0.020279 |
  | PCRE2 | 659 | 0.018157 |
  | Perl | 644 | 0.017744 |
  | OpenBLAS | 642 | 0.017689 |
  | CMake | 616 | 0.016973 |
  | Automake | 607 | 0.016725 |
  | MariaDB | 596 | 0.016421 |
  | Pysam | 591 | 0.016284 |
  | ncurses | 585 | 0.016118 |
  | HISAT2 | 580 | 0.015981 |
  | libreadline | 574 | 0.015815 |
  | cURL | 572 | 0.015760 |
  | flex | 561 | 0.015457 |
  | Bison | 561 | 0.015457 |
  | pkg-config | 559 | 0.015402 |
  | texinfo | 558 | 0.015374 |
  | libtool | 558 | 0.015374 |
  | bsddb3 | 558 | 0.015374 |
  | PCRE | 558 | 0.015374 |
  | DB | 558 | 0.015374 |
  | parallel | 555 | 0.015292 |
  | tbb | 554 | 0.015264 |
  | itpp | 554 | 0.015264 |
  | Trim_Galore | 553 | 0.015237 |
  | datamash | 540 | 0.014879 |
  | RNA-SeQC | 531 | 0.014631 |
  | Bowtie | 514 | 0.014162 |
  | cuDNN | 475 | 0.013088 |
  | gnuplot | 452 | 0.012454 |
  | snakemake | 447 | 0.012316 |
  | DIAMOND | 414 | 0.011407 |
  | BBMap | 387 | 0.010663 |
  | IPython | 377 | 0.010387 |
  | nextflow | 370 | 0.010195 |
  | FLASH | 324 | 0.008927 |
  | gatk | 283 | 0.007797 |
  | bcl2fastq2 | 262 | 0.007219 |
  | GLib | 256 | 0.007054 |
  | JupyterLab | 249 | 0.006861 |
  | Mesa | 248 | 0.006833 |
  | GSL | 248 | 0.006833 |
  | HLA-HD | 247 | 0.006806 |
  | GTK3 | 247 | 0.006806 |
  | PEAR | 222 | 0.006117 |
  | MATLAB | 217 | 0.005979 |

### Modules ranked by most unique users
Splunk Query: `index = syslog LMOD: top=yes | stats dc(user) as Users by name | sort - Users`

  | Name | Users |
  |------|-------|
  | R | 175 |
  | RStudio-Server | 133 |
  | Python | 122 |
  | fhR | 121 |
  | SAMtools | 102 |
  | Singularity | 67 |
  | fhPython | 57 |
  | Anaconda3 | 54 |
  | FastQC | 46 |
  | BEDTools | 44 |
  | Java | 44 |
  | Bowtie2 | 43 |
  | picard | 40 |
  | CellRanger | 37 |
  | SRA-Toolkit | 34 |
  | STAR | 34 |
  | OpenBLAS | 31 |
  | BCFtools | 28 |
  | GATK | 28 |
  | awscli | 28 |
  | BWA | 25 |
  | Nextflow | 25 |
  | MariaDB | 24 |
  | cutadapt | 24 |
  | nextflow | 24 |
  | CMake | 23 |
  | deepTools | 23 |
  | JupyterLab | 21 |
  | Perl | 21 |
  | bcl2fastq2 | 21 |
  | MACS2 | 19 |
  | MultiQC | 18 |
  | foss | 18 |
  | snakemake | 18 |
  | tabix | 18 |
  | cURL | 17 |
  | tmux | 17 |
  | X11 | 16 |
  | cellranger | 16 |
  | rstudio | 16 |
  | "Kent_tools" | 15 |
  | Pysam | 15 |
  | "Trim_Galore" | 15 |
  | GCC | 14 |
  | cromwell | 14 |
  | "rstudio-server" | 14 |
  | TensorFlow | 13 |
  | "BLAST+" | 12 |
  | Automake | 11 |
  | IGV | 11 |
  | Miniconda3 | 11 |
  | PCRE2 | 11 |
  | libreadline | 11 |
  | ncurses | 11 |
  | plink | 11 |
  | "Aspera-Connect" | 10 |
  | BBMap | 10 |
  | Beast | 10 |
  | Bison | 10 |
  | GCCcore | 10 |
  | HISAT2 | 10 |
  | MATLAB | 10 |
  | OptiType | 10 |
  | Subread | 10 |
  | flex | 10 |
  | parallel | 10 |
  | "pkg-config" | 10 |
  | BEDOPS | 9 |
  | DB | 9 |
  | PCRE | 9 |
  | PubWeb | 9 |
  | "SciPy-bundle" | 9 |
  | bsddb3 | 9 |
  | libpciaccess | 9 |
  | libtool | 9 |
  | texinfo | 9 |
  | AlphaFold | 8 |
  | cuDNN | 8 |
  | itpp | 8 |
  | kallisto | 8 |
  | tbb | 8 |
  | BamTools | 7 |
  | "FASTX-Toolkit" | 7 |
  | FriBidi | 7 |
  | HarfBuzz | 7 |
  | Homer | 7 |
  | MEME | 7 |
  | OpenMPI | 7 |
  | PEAR | 7 |
  | PLINK2 | 7 |
  | Pandoc | 7 |
  | RSeQC | 7 |
  | Trimmomatic | 7 |
  | XZ | 7 |
  | "beagle-lib" | 7 |
  | git | 7 |
  | seqtk | 7 |
  | Apptainer | 6 |
  | Bowtie | 6 |
  | CUDA | 6 |
  | EMBOSS | 6 |
  | Emacs | 6 |
  | HTSeq | 6 |
  | IPython | 6 |
  | Ruby | 6 |
  | Trinity | 6 |
  | minimap2 | 6 |
  | ArrayFire | 5 |
  | CNVkit | 5 |
  | GLib | 5 |
  | GMP | 5 |
  | IGVTools | 5 |
  | "MAGeCK-VISPR" | 5 |
  | MiXCR | 5 |
  | PLINK | 5 |
  | "RNA-SeQC" | 5 |
  | bcl2fastq | 5 |
  | cirro | 5 |
  | fhCellRanger | 5 |
  | texlive | 5 |
