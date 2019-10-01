---
title: Reduced Graphics Capabilities
date: 2018-01-01
---

## R Graphics Capabilities

R defaults to Xlib graphics capabilities. If DISPLAY is not set the Xlib
graphics engine will issue a cryptic error.  

```
Error in .External2(C_X11, paste0("png::", filename), g$width, g$height, :
 unable to start device PNG
In addition: Warning message:
In png(tempfile()) : unable to open connection to X11 display ''
```

To fix either define an X11 display or
you can change the default grahics engine. If you mainly create graphs for
printing with png or PDF it is recomended that you use the **cairo** graphics.

#### How to change the default graphics?
Use ```getOption()``` and ```option()``` to query and set graphics capablities.
The default can also be changed by adding bitmapType to .Rprofile

```
> getOption('bitmapType')
[1] "Xlib"
> options(bitmapType = 'cairo')
> getOption('bitmapType')
[1] "cairo"
```

#### Graphics Issues
*This issue has beeen resolved since R-3.3.2.*

We have recently discovered that recent easybuilds of R (R modules containing
"foss" or "intel") were not compiled correctly, leading to missing graphics
capabilities if you are not running an X display:

```
> capabilities()
       jpeg         png        tiff       tcltk         X11        aqua 
      FALSE       FALSE       FALSE        TRUE       FALSE       FALSE 
   http/ftp     sockets      libxml        fifo      cledit       iconv 
       TRUE        TRUE        TRUE        TRUE        TRUE        TRUE 
        NLS     profmem       cairo         ICU long.double     libcurl 
       TRUE       FALSE       FALSE        TRUE        TRUE        TRUE 
```

note in the above that jpeg, png, tiff are not indicated as capabilities.
We will be rebuilding these R modules to include cairo which will enable
these capabilities.  To work around this, use an earlier version of
R (3.2.2 and earlier) or ensure that your session has an X display- the
easiest way to do this is to use a NoMachine client to connect to the
rhino/gizmo nodes you need.
