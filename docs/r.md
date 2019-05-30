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

### [R Issues and Solutions]({{ site.baseurl }}/rissues/)
