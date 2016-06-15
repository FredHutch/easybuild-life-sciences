# Running EasyBuild in an LXD container

We found that our default installs of Ubuntu included too many *-dev packages. In a world where ./configure is greedy, this can cause include confusion.

## Pre-requisites
This is all running on an Ubuntu 16.04 LTS host system. The LXD packages must be installed (lxd and lxd-client). Sufficient space shoud be given to LXD. We are using ZFS on a local disk volume.

## lxc_eb_setup.sh
The script `lxc_eb_setup.sh` automates the creation of containers with EasyBuild pre-installed (including toolchains optionally). It takes several arguments:

 * `--image-name` (required) - image to be used as a base image for the container (string) (remote ex: ubuntu:16.04/amd64)
 * `--container-name` (required) - specify name for the container (string)
 * `--group` - group to manage permission in the easybuild directories (string)
 * `--user` - additional users (beyond the user running the script) to be included inside the container (string) (multiple)
 * `--ephemeral` - create the new container as ephemeral (yes/no - default yes)
 * `--luaver` - specify the version of Lua to install (string - default 5.3.3)
 * `--luarocksver` - specify the version of LuaRocks to install (string - default 2.3.0)
 * `--lmodver` - specify the version of Lmod to install (string - default 6.3)
 * `--ebver` - specify the version of EasyBuild to install (string - default 2.8.1)
 * `--lm-license-file` - value to set LM_LICENSE_FILE environment variable in EasyBuild modulefile
 * `--toolchain` - EasyBuild toolchain to pre-build in container (string - '.eb' will be appended) (multiple)
 
 Values are separated from argument with a space (Ex: `lxc_eb_setup.sh --image-name ubuntu:14.04/amd64 --container-name eb-14.04`)
 
 ## Caveats/Assumptions/Best Practices
 Odd things about us:
 * we use EasyBuild with our regular POSIX user accounts
 * we are members of a POSIX group, which is uesd to manage permissions
 
 
