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

Beginning with R-3.5.1 there are two modules for each release of R. The two
 modules are a base package and a custom Fred Hutch module.
 The Fred Hutch R modules have the suffix `-fh1`. The `-fh1` module is specific
 to the Hutch, and contains libraries that have been requested by Hutch users beyond those in the base module.
 The Fred Hutch module inherits the modules from the base module. The base
 module is maintained by the EasyBuild community. 

### Requesting Additional Libraries
Adding every user request for libraries is becoming a challenge to support.
 The Fred Hutch R module has close to 1,000 libraries. Users are encouraged
 to install custom R libraries in their home directories. Users can submit install
 request for libraries that require system libraries by emailing `scicomp`.  

### User Installed Libraries
Use `install.packages("package-name")` to install packages in your home directory.
 Newer versions of R set the search path for user installed libraries. The
 search path is based on Major and Minor version numbers. Packages installed for
 R-3.4.x will not be loaded by R-3.5.x. If you are having install issues verify
 that your install library path is defined correctly. 
```
> Sys.getenv('R_LIBS_USER')
[1] "~/R/x86_64-pc-linux-gnu-library/3.6"
``` 
If an error arises during the installation that says the library path is not writable, you may need to specify the desired library path like this:
```
install.pacakages("package-name", lib = "~/R/x86_64-pc-linux-gnu-library/3.6")
```


