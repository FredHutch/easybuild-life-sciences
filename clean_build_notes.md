# Notes from running EasyBuild on a 'clean' system

These notes come from my experience using LXD to run EasyBuild inside 'base' Linux containers.

My base containers are all from [linuxcontainers.org](https://images.linuxcontainers.org:8443/):

 * Ubuntu 14.04 amd64
 * CentOS 7 amd64
 * Opensuse 13.2 amd64

## Notes from trying builds:

| Build Attempted | OS | Notes | Remedy | Temporary Remedy |
| --- | --- | --- | --- | --- |
| EasyBuild | Ubuntu 14.04 | os pkgs needed: build-essential, python-pygraph, libreadline-dev | Script installation into container | N/A |
| EasyBuild | Ubuntu 16.04 | os pkgs needed: build-essentials, python-minimal, python-pygraph | Script installation into container | N/A |
| foss-2016a | Ubuntu 14.04 | os pkgs needed: libibverbs-dev | None, these are called out in osdependencies | N/A |
| intel-2016a | Ubuntu 14.04 | os pkgs needed: gcc-multilib, g++-multilib (may be an Ubuntu thing) | Add to .eb as os dependency? | Manual install |
| many pkgs | Ubuntu 14.04 | missed build dependency: pkgconfig | Should be in toolchain? | Manual install of os pkg :( |
| flex-2.5.39-foss-2016a.eb | Ubuntu 14.04 | missed build dependency: m4 | Add to .eb as build dependency? | manual install of os pkg :( |

## Base OS package inclusion table

| Pkg | Ubuntu 14.04 | Ubuntu 16.04 | CentOS 7 | Opensuse 13.2 |
| --- | --- | --- | --- | --- |
| pkgconfig | No | No | Yes | No |
| m4 | No | No | No | No | No |
