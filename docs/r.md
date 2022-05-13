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

The Fred Hutch R modules have the suffix `-fh1`. The `-fh1` module is spefific
 to the Hutch, and contains libraries that have been requested by Hutch users.
 The Fred Hutch module inherits the modules from a base R module that 
 is maintained by the EasyBuild community. 
 Users should use the 'fh' versions of R.

### Requesting Modules ###
Adding every user request for libraries is becoming a challenge to support.
 The Fred Hutch R module has close to 1,000 libraries. Users are encouraged
 to install custom R libraries in their home directories. Users can submit install
 request for libraries that require system libraries.

### User Installed R Modules
The Fred Hutch R module has over 1,000 libraries, but it might not have the one
library that you need. R libraries into your home directory if they not otherwise
available.  Use the `install.packages()` function to install R libraries.

```
> install.packages("dsa")
Installing package into ‘/app/easybuild/software/R/3.6.0-foss-2016b-fh1’
(as ‘lib’ is unspecified)
Warning in install.packages("dsa") :
  'lib = "/app/easybuild/software/R/3.6.0-foss-2016b-fh1"' is not writable
Would you like to use a personal library instead? (yes/No/cancel) yes
```

The default path for libraries is the system path which is not writable. Chose the
option to create a personal library. Before adding a personal library you should
check the location where it will be written.  The `.libPaths()` fuction will show
the locations that libraries are searched.

```
> .libPaths()
[1] "/app/easybuild/software/R/3.6.0-foss-2016b-fh1"
[2] "/home/jfdey/R/x86_64-pc-linux-gnu-library/3.6"
[3] "/app/easybuild/software/R/3.6.0-foss-2016b/lib/R/library"
>
```

Notice one of the paths is in my "Home" directory and contains the major.minor
 version of R `3.6` in the path. If your library path is not versioned you might have
defined **R_LIBS_USER** in your .Rprofile configuation file.  Use `Sys.getenv()`
 to check your default user path.

```
Sys.getenv("R_LIBS_USER")
> Sys.getenv("R_LIBS_USER")
[1] "~/R/x86_64-pc-linux-gnu-library/3.6"
```
### Install BioConductor Package
BioConductor packages are released and updated with R releases. Each release of Bioconductor is developed to work with a
speicific version of R. BiocManager is part of fhR, and BiocManaager will only install packages which match the
version of BioCondutor. If you attempt to install a mismatch version module you will get an error. When working with
a version of R that is not current check the version BiocManager.

|---|---|
| R Version | Bioconductor Version | 
| 4.2.0 | 3.15 |
| 4.1.2 | 3.14 |
| 4.1.1 | 3.13 |
| 4.1.0 | 3.13 |
| 4.0.5 | 3.12 |
| 4.0.4 | 3.12 |
| 4.0.3 | 3.12 |
| 4.0.2 | 3.11 |
|---|---|

### Issues with R Libraries
One of the most frequent issues are errors with loading libraries.  There are
two major issues; the library is out of date or there is a newer version available
from the system that is not being loaded. Use the `packageVersion("snow")`
function to show the library version and `update.packages()` function to update
out of date packages.

### [R Issues and Solutions]({{ site.baseurl }}/rissues/)
