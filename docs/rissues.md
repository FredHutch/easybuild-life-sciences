---
title: "R Issues Solutions"
layout: collection
classes: wide
permalink: /rissues/
collection: rissues
entries_layout: list # list or grid (default),
show_excerpts:  true #(default), false
sort_by: date # date or title (default)
sort_order:  reverse # (default),forward or reverse
sidebar:
  nav: "docs"
---

### User Installed Libraries ###
Use `install.packages("package-name")` to install packages in your home directory.
 Newer versions of R set the search path for user installed libraries. The
 search path is based on Major and Minor version numbers. Packages installed for
 R-3.4.x will not be loaded by R-3.5.x. If you are having install issues verify
 that your install library path is defined correctly.

```
> Sys.getenv('R_LIBS_USER')
[1] "~/R/x86_64-pc-linux-gnu-library/3.6"
```
