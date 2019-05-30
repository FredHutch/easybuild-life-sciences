---
title: libtiff not found
date: 2017-01-01
---

You may encounter errors with some tools when you have this R module loaded.  You will get messages indicating version `LIBTIFF_4.0' not found .  This is because of differences in the way that tools are built for the operating system and the library paths set up when loading the module.  What is typically needed is to only use tools loaded by modules. 

Specific notes on tools are listed below.

emacs
Loading one of these environment modules will break the installed emacs like this:

$ ml R/3.3.3-foss-2016b-fh2
$ emacs
emacs: /app/easybuild/software/LibTIFF/4.0.6-foss-2016b/lib/libtiff.so.5: version `LIBTIFF_4.0' not found (required by emacs)
The fix is to load an Emacs module along with your R module. You will get a newer version of Emacs as well.

$ ml R/3.3.3-foss-2016b-fh2
$ ml Emacs/25.1-foss-2016b
As long as the toolchain string (in this case 'foss-2016b') matches between R and Emacs, they will work together.

ESS - Emacs Speaks Statistics
Gizmo nodes (and rhino) have the OS package 'ess' installed. Emacs Speaks Statistics provides a bridge between Emacs and R. We are currently working on a central install of ess, but in the meantime, you will need to install it after loading the R and Emacs modules if you need to use it. The instructions are here .

Latex/Tetex Document Generation
Some operations (e.g. building packages, vignettes, etc) may fail with an error similar to this:

pdflatex: /app/easybuild/software/LibTIFF/4.0.6-foss-2016b/lib/libtiff.so.5: version `LIBTIFF_4.0' not found (required by /usr/lib/x86_64-linux-gnu/libpoppler.so.44)

The tetex package is used to build and convert much of the documentation in R libraries and such.  In these cases, load the "tetex" module prior to building:

$ ml tetex 

