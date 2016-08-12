%title: EasyBuild _is_ easy
%author: Ben McGough
%date: 2016.08.10
 
-> EasyBuild <-
 
-> # EasyBuild is building software with ease <-
 
*Features*
 
- multiple versions and builds of many software packages
- automated dependency system
- re-producible software package builds
- automated environment modules managment
- curated toolchains
 

-------------------------------------------------

-> # EasyBuild Terms <-

- *easyconfig:* a file that describes a software package version
- *toolchain:* a defined collection of compilers and support libraries
- *environment modules:* a system of managing environment variables


-------------------------------------------------

-> # Example of an EasyBuild <-
=====

Demo a quick EasyBuild in shell


-------------------------------------------------
 
-> # Environment Modules <-
=====
 
Easybuild uses Environment Modules to automatically wrap software packages it builds.

This is how it makes software packages available, and isolated from eachother.

Environment Modules is a system that:
 
- sets and unsets environment variables
- avoids conflicts
- allows some administrative control

 
-------------------------------------------------
 
-> # Environment Variables: shell integration <-
 
Environment Variables integrates with the users shell

- shell function module()
- calls modulecmd with options and shell name
- modulecmd outputs commands for specified shell to set/unset variables
- I recommend the use of Lmod as it is the most widely used by EasyBuilders

 
-------------------------------------------------
 
-> # Example of module use <-
 
Here is a quick example of a module in use:
 
    $ echo $LD_LIBRARY_PATH

    $ module load zlib
    $ echo $LD_LIBRARY_PATH
    LD_LIBRARY_PATH="/app/zlib/1.2.8/lib:/app/lib:/app/lib64:/app/lib:/app/lib64";
    $ module unload zlib
    $ echo $LD_LIBRARY_PATH
    
    $

 
-------------------------------------------------
 
-> # EasyBuild Toolchains <-
 
The example showed a quick easyconfig using the "dummy" toolchain - whatever compilers your OS provides.

EasyBuild provides 54 toolchains that include multiple architecures, propietary compilers, and various features.

We decided to use:

- *foss-n:* Free Open Source Software - GCC and friends
- *intel-n:* Intel C and Fortran

 
-------------------------------------------------
-> # Why different toolchains? <-

Demo of R performance OS vs. foss cs. intel
_

-------------------------------------------------

-> # Bootstrapping EasyBuild <-

- bootstrap.py

-------------------------------------------------

-> # EasyBuild Environment <-

EasyBuild environment and options - paths, shared paths, Lmod hooks, logging


-------------------------------------------------

move toolchain difference earlier - toolchaina nd reproducibility
better transition between env modules and toolchains

bootstrap slide - quick steps... "getting easybuild going quickly" is python

easybuild environment: easybuid at FH - link to easybuild-lifesciences
 multi user
 shared paths
 no root

swap zlib for a binary with 'which' example

new slide: R dependency example

new slides: easybuild metrics/stats - totals, FH totals + packages for labs, life sciences

new slide: easyconfig simple example

new slide: easyconfig complex example

new slide: collaboration and community

new slide: new easyconfig example

new slide: easybuild github and easyconfigs search

new slide: contribution, mailing list, community, processes, testing, github-based

new slide: outlook/next steps - lifesciences easybuild community for collaboration using github,

new slide: containers, dockers, contrast and compare

new slide: reproducibility/citation lxc/docker containers
