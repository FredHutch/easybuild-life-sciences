# ls2
Life Sciences Software

# Overview
The Life Sciences Software (LS2) project aims to normalize the build of software packages across multiple technologies.

## Components
LS2 is a collection of open source components.

*[EasyBuild](https://easybuilders.github.io/easybuild/)
  *We use EasyBuild to compile software packages and provide reproducibility.
*[Lmod](https://github.com/TACC/Lmod)
  *EasyBuild uses Environment Modules to manage software packages, and LS2 uses Lmod to provide Environment Modules.
*[Docker](https://www.docker.com/)
  *LS2 can produce Docker containers with one or more software packages.
*[Ubuntu](https://www.ubuntu.com/)
  *LS2 uses Ubuntu as its primary platform, and creates Docker containers based on Ubuntu containers. Note that EasyBuild uses CentOS as their primary platform, and extending LS2 to CentOS would likely be pretty easy.

## LS2 Architecture
This is the hierarchy of LS2 containers:

Name/Repo | FROM | Reason | Notes
--- | --- | ---
[ls2_ubuntu](https://github.com/FredHutch/ls2_ubuntu) | ubuntu | simple 'freeze' of the public ubuntu container | OS pkgs added: bash, curl, git
[ls2_easybuild](https://github.com/FredHutch/ls2_easybuild) | ls2_ubuntu | Adding EasyBuild and Lmod | OS pkgs added: python, lua
[ls2_easybuild_foss](https://github.com/FredHutch/ls2_easybuild_foss) | ls2_easybuild | Adding the 'foss' toolchain | OS pkgs added: libibverbs-dev, lib6c-dev, bzip2, unzip, make, xz-utils
[ls2](https://github.com/FredHutch/ls2) | ls2_easybuild_foss | This 'demo' repo | does not produce a container directly
[ls2_r](https://github.com/FredHutch/ls2_r) | ls2_easybuild_foss | Our 'R' build | OS pkgs added: awscli

## Container Architecture
* default user 'neo' (UID 500, GID 500) /home/neo
* Lmod and EasyBuild are installed into /app
* EasyConfigs are copied into /app/fh_easyconfigs
* sources are copied and downloaded into /app/sources
* EasyBuild is run in a single 'RUN' command to reduce container layer size:
  * installs specified OS packages
  * runs EasyBuild
  * uninstalls specified OS packages

## Use Cases

### Create a container
The initial reason for LS2 is to create Docker containers with EasyBuilt software packages to mirror those available on our HPC systems. We realize that containerizing common software packages will be key in leveraging many new technologies like AWS Batch.

The intention is to use multiple LS2 containers in step to achieve a pipeline. Having the same software packages compiled in the same ways as we have deployed to our traditional HPC cluster enables users to focus on pipeline building and not software troubleshooting when moving to different compute methods.

### EasyBuild testing
Building or updating an EasyConfig can be time-consuming. Many existing technologies help to automated the `docker build` process, so LS2 opens these up to EasyBuild. Even if you run EasyBuild traditionally to deploy built packages to a private directory or shared software archive, you can test those EasyConfigs in an LS2 build to ensure your production build will be successful.

As we are building in a container with minimal installed packages, it is easy to find OS dependencies that are unstate in EasyConfigs. Some examples of this range from pkg-config (a default package in CentOS but not Ubuntu) to utilities like unzip and bzip2. Some dependencies are intentional like OpenSSL (better to pull presumably-updated OS packages than possibly stale EasyBuild packages) and some are oversights easy to miss when you are building in a fully-installed OS (ex: make is not present in the foss-n toolchains).

### Manage a traditional archive
We use EasyBuild to install software packages onto an NFS volume. This volume is then shared to our HPC and other systems to enable software package use on those platforms. LS2 can still be used to deploy packages in this way by mounting the NFS volume into the container and performing a build. This process isolates the EasyConfig development process from your live package archive or volume.

# HOWTO
There are two sections here. First case covers building an existing or new EasyConfig, and the second covers using a built container to deploy a software package to an existing software archive or volume.

## Fork and Add
Steps to build a new LS2 container are pretty straight-forward, but assume some knowledge of EasyBuild, Lmod, and Docker. A more details explanation follows.

1. Fork this repo
1. Add required EasyConfig files that are not in the EasyBuild repo to /easyconfigs
1. Add sources to the sources/ folder of the repo (for sources <50MB in size that cannot easily be downloaded)
1. Add URLs to sources/download_sources.sh to download sources during `docker build` (for larger sources, perhaps placed in the cloud for easier download)
1. Edit the Dockerfile to adjust the following:
  * EASYCONFIG_NAME - this is the name of the package to be built
  * INSTALL_OS_PKGS - these packages will be installed (by root) prior to running EasyBuild
  * UNINSTALL_OS_PKGS - these packages will be uninstalled at the end of running EasyBuild
1. Run `docker build`

## Then add to /app
We keep our deployed software package on an NFS volume that we mount at /app on our systems (can you guess why LS2 builds into /app rather than .local in the container?). In order to use your recently build LS2 software package container to deploy the same package into our /app NFS volume, use these steps:

1. Complete 'Fork and Add' steps to produce a successful container with your software package
1. Run `docker build` once again, this time with `--target fh_deploy` - this will run quickly and build a second container
1. Run that container with our package deploy location mapped in to /app like this: `docker run ls2_r_fh_deploy -v /app:/app`

The steps above will produce a container with EasyBuild and all the pieces necessary, with the actual EasyBuild command set as the entrypoint. Running the container will trigger the EasyBuild run, and the resulting output will be placed into the /app volume outside the container.

Note that this overrides the Lmod in the container, so if version parity is important to you, you'll want to keep your Lmod in sync with the LS2 Lmod.
