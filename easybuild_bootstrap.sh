#!/bin/bash

#####
# # bootstrap script to demo easybuild
# # only tested on Ubuntu 16.04 
# # on fresh linux install create a new user eb
# # and add that user to sudoers 
# adduser --disabled-password --gecos "" eb
# echo "eb ALL=(ALL:ALL) NOPASSWD:ALL" >  /etc/sudoers.d/zz_eb
#  # then execute script as user eb
######

# variables

# Internal
LUA_BASE_URL="http://www.lua.org/ftp/lua-"
LUAROCKS_BASE_URL="http://luarocks.org/releases/luarocks-"
LMOD_BASE_URL="https://github.com/TACC/Lmod/archive/"
EB_BASE_URL="https://github.com/hpcugent/easybuild-framework/raw/easybuild-framework-v%s/easybuild/scripts/bootstrap_eb.py"

LUA_VER="5.3.3"   # verion of lua to install into the container
LUAROCKS_VER="2.3.0"   # version of luarocks package manager to install into the container
LMOD_VER="6.3"   # version of Lmod to install into the container
EB_VER="2.9.0"   # verison of EasyBuild to bootstrap in the container

EB_DIR="/easybuild"   # location for EasyBuild directory tree

# install lua, luarocks, luafilesystem, and luaposix
function lua_install {

  # install readline dev package
  sudo apt-get install -y libreadline-dev make unzip


  # get Lua
  lua_url="${LUA_BASE_URL}$LUA_VER.tar.gz"
  wget -O /tmp/lua-$LUA_VER.tar.gz $lua_url
  if [ "$?" != "0" ]; then
    echo "Oops, retrieving lua failed!" 1>&2
    exit 1
  fi

  # extract
  tar -xf /tmp/lua-$LUA_VER.tar.gz -C /tmp

  # build and install
  cd /tmp/lua-$LUA_VER && sudo make linux install

  # get luarocks
  luarocks_url="${LUAROCKS_BASE_URL}$LUAROCKS_VER.tar.gz"
  wget -O /tmp/luarocks-$LUAROCKS_VER.tar.gz $luarocks_url
  if [ "$?" != "0" ]; then
    echo "Oops, retrieving luarocks failed!" 1>&2
    exit 1
  fi

  # extract
  tar -xf /tmp/luarocks-$LUAROCKS_VER.tar.gz -C /tmp

  # build and install
  cd /tmp/luarocks-$LUAROCKS_VER && ./configure && make build && sudo make install

  # use luarocks to install luaposix and luafilesystem
  sudo luarocks install luaposix
  sudo luarocks install luafilesystem

  # uninstall libreadline-dev for a clean system
  sudo apt-get remove -y libreadline-dev
}

# install Lmod
function lmod_install {

  # install Tcl 
  sudo apt-get install -y tcl

  # get Lmod
  lmod_url="${LMOD_BASE_URL}$LMOD_VER.tar.gz"
  wget -O /tmp/Lmod-$LMOD_VER.tar.gz $lmod_url
  if [ "$?" != "0" ]; then
    echo "Oops, retrieving Lmod failed!" 1>&2
    exit 1
  fi

  # extract
  tar -xf /tmp/Lmod-$LMOD_VER.tar.gz -C /tmp

  # build and install
  cd /tmp/Lmod-$LMOD_VER && ./configure && sudo make install

  # link into /etc/profile so shells can use module function
  sudo ln -s /usr/local/lmod/lmod/init/profile /etc/profile.d/modules.sh
  source /etc/profile.d/modules.sh
}

# bootstrap easybuild
function eb_bootstrap {

  sudo apt-get install -y python-minimal python-pygraph git libibverbs-dev libssl-dev build-essential
  #sudo apt-get install -y environment-modules
  # get EB
  eb_url=$(printf "$EB_BASE_URL" "$EB_VER")
  wget -O /tmp/bootstrap_eb.py $eb_url
  if [ "$?" != "0" ]; then
    echo "Oops, retrieving bootstrap_eb.py failed!" 1>&2
    exit 1
  fi

  # create $EB_DIR
  myself=$(whoami)
  sudo mkdir -p $EB_DIR
  sudo chown $myself $EB_DIR

  # bootstrap it
  python /tmp/bootstrap_eb.py $EB_DIR

  # pop in some useful environment variables to our EasyBuild modulefile
  (cat <<EOF
set ebDir "$EB_DIR"
setenv EASYBUILD_SOURCEPATH "\$ebDir/sources"
setenv EASYBUILD_BUILDPATH "\$ebDir/build"
setenv EASYBUILD_INSTALLPATH_SOFTWARE "\$ebDir/software"
setenv EASYBUILD_INSTALLPATH_MODULES "\$ebDir/modules"
setenv EASYBUILD_REPOSITORYPATH "\$ebDir/ebfiles_repo"
setenv EASYBUILD_LOGFILE_FORMAT "\$ebDir/logs,easybuild-%(name)s-%(version)s-%(date)s.%(time)s.log"
setenv EASYBUILD_MODULES_TOOL "Lmod"
EOF
  ) >> $EB_DIR/modules/all/EasyBuild/$EB_VER

  sudo sh -c "echo export MODULEPATH=/easybuild/modules/all:'\$MODULEPATH' > /etc/profile.d/modules_eb.sh"
  sudo sh -c "export EASYBUILD_MODULES_TOOL=Lmod >> /etc/profile.d/modules_eb.sh"

  source /etc/profile.d/modules_eb.sh
}

# main

# checking requirements
mem=$(awk '( $1 == "MemTotal:" ) { printf "%.0f", $2/1024/1024 }' /proc/meminfo)
if [[ $mem -lt 8 ]]; then
    printf "This script requires a minimum of 8GB ram, 8 GB disk and 8 cores !"
    exit 1
fi 
printf "\nUpdating packages..."
apt-get update
printf "\nInstalling Lua...\n"
lua_install
printf "Installing Lmod...\n"
lmod_install
printf "Bootstrapping EasyBuild...\n"
eb_bootstrap
echo "please log out and log in again or source /etc/profile.d/modules.sh and /etc/profile.d/modules_eb.sh"
exit 0

