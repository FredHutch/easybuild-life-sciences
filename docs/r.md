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

Beginning with R-3.5.1 there are two modudles for each release of R. The two
 modules are a base package and a custom Fred Hutch module.
 The Fred Hutch R modules have the suffix `-fh1`. The `-fh1` module is spefific
 to the Hutch, and contains libraries that have been requested by Hutch users.
 The Fred Hutch module inherits the modules from the base module. The base
 module is maintained by the EasyBuild community. 

### Requesting Modules ###
Adding every user request for libraries is becoming a challenge to support.
 The Fred Hutch R module has close to 1,000 libraries. Users are encouraged
 to install custom R libraries in their home directories. Users can submit install
 request for libraries that require system libraries.

### User Installed Libraries ###
R looks for libraries based on two environment variables; a system path,
 and a user path. Check your environment after starting R with the command
 `Sys.getenv('R_LIBS_USER')`. `R_LIBS_USER` path starts with a Tilda
 to point to your directory. The output should look like this
```
> Sys.getenv('R_LIBS_USER')
[1] "~/R/x86_64-pc-linux-gnu-library/3.6"
``` 
The major and minor version numbers are part of the user search path.
 User installed libraries for 3.5.x R need to be re-installed for use with 3.6.x.
 **Note:** There is a bug since R-3.5.3 which does not set the default path for
 installing local libraries. When installing local libraries, explicitly set
 the install path:
```
# CRAN example
install.packages("ggplot2", lib=Sys.getenv("R_LIBS_USER"))
# Github Example
devtools::install_github('FredHutch/tgR', lib=Sys.getenv("R_LIBS_USER"))
```
