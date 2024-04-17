---
title: "Available R Modules"
layout: collection
classes: wide
permalink: /r/
collection: r
entries_layout: list # list or grid (default),
show_excerpts:  true #(default), false
sort_by: date # date or title (default)
sort_order:  reverse # (default),forward or reverse
sidebar:
  nav: "docs"
---

The Fred Hutch R modules have the prefix `fh`. fhR modules are specific to the Hutch and
contain libraries that Hutch users have requested. The Fred Hutch module inherits all
the libraries from EasyBuild R.  fhR module contains many Bioconductor libraries. Cluster
users should use the 'fh' versions of R.

### Requesting Modules 
The Fred Hutch R module has over 1,000 libraries. Adding every user request for libraries is becoming a challenge to support. Users are encouraged to install R libraries in their home directories. To have complex libraries installed which require additional system libraries or data open a ticket with Scicomp.

### Installing R Modules
Additional R libraries can be installed into your home directory.
Use the `install.packages()` function to install R libraries.
Since verion 4.0 of R the library paths in your home directory are
based on the Major.Minor version of R. Users should not set
R_LIB or any other environment variables that effect searching of
R library paths. After starting R you can verify library paths with
`> .libPaths()`. Your personal library path in your home directory
should app first in the list.

```
> .libPaths()
[1] "/home/jfdey/R/x86_64-pc-linux-gnu-library/4.2"
[2] "/app/software/fhR/4.2.2.1-foss-2021b"
[3] "/app/software/R/4.2.2-foss-2021b/lib/R/library"
```

The first path is in my "Home" directory and contains the major.minor
 version of `R-4.2` in the path. If your library path is not versioned you might have
defined `R_LIBS` or `R_LIBS_USER` in your .Rprofile configuation file.  Use `Sys.getenv()`
 to check your default user path.

```
Sys.getenv("R_LIBS_USER")
> Sys.getenv("R_LIBS_USER")
[1] "~/R/x86_64-pc-linux-gnu-library/4.2"
```

### Install BioConductor Package
BioConductor packages are released and updated in step with R releases. Each release of Bioconductor matches a specific version of R. BiocManager is installed in fhR. BiocManager will only install libraries that match the release version of BioConductor. When working with an older R version, you need to know the Bioconductor version that matches. Search the older releases of libraries that match the Bioconductor release. 

| R Version | Bioconductor Version |
|---|---|
| 4.3.3 | 3.18 |
| 4.3.1 | 3.17 |
| 4.3.0 | 3.17 |
| 4.2.2 | 3.16 |
| 4.2.0 | 3.15 |
| 4.1.2 | 3.14 |
| 4.1.1 | 3.13 |
| 4.1.0 | 3.13 |
| 4.0.5 | 3.12 |
| 4.0.4 | 3.12 |
| 4.0.3 | 3.12 |
| 4.0.2 | 3.11 |

### Issues with R Libraries
One of the most frequent issues user issues are with loading libraries.  There are
two major issues; A user installed R library that is out of date or a library 
requires a newer version of a user installed library.
Use the `packageVersion("snow")`
function to show the library version and `update.packages()` function to update
the out of date packages. Often a user installed library is already part of the
`fhR` and is conflicting with the module version. 
