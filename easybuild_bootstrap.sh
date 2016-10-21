#!/bin/bash

# variables

# Internal
LUA_BASE_URL="http://www.lua.org/ftp/lua-"
LUAROCKS_BASE_URL="http://luarocks.org/releases/luarocks-"
LMOD_BASE_URL="https://github.com/TACC/Lmod/archive/"
EB_BASE_URL="https://github.com/hpcugent/easybuild-framework/raw/easybuild-framework-v%s/easybuild/scripts/bootstrap_eb.py"

LUA_VER="5.3.3"   # verion of lua to install into the container
LUAROCKS_VER="2.3.0"   # version of luarocks package manager to install into the container
LMOD_VER="6.3"   # version of Lmod to install into the container
EB_VER="2.8.2"   # verison of EasyBuild to bootstrap in the container

EB_DIR="/easybuild"   # location for EasyBuild directory tree

# install lua, luarocks, luafilesystem, and luaposix
function lua_install {

  # install readline dev package
  apt-get install ibreadline-dev

  # install unzip package - apparently luarocks silently requires this
  apt-get install unzip

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
  cd /tmp/lua-$LUA_VER && make linux install

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
  cd /tmp/luarocks-$LUAROCKS_VER && ./configure && make build && make install

  # use luarocks to install luaposix and luafilesystem
  luarocks install luaposix
  luarocks install luafilesystem

  # uninstall libreadline-dev for a clean system
  apt-get remove libreadline-dev
}

# install Lmod
function lmod_install {

  # install Tcl 
  apt-get install tcl

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
  cd /tmp/Lmod-$LMOD_VER && ./configure && make install

  # link into /etc/profile so shells can use module function
  ln -s /usr/local/lmod/lmod/init/profile /etc/profile.d/modules.sh
}

# bootstrap easybuild
function eb_bootstrap {

  # get EB
  eb_url=$(printf "$EB_BASE_URL" "$EB_VER")
  wget -O /tmp/bootstrap_eb.py $eb_url
  if [ "$?" != "0" ]; then
    echo "Oops, retrieving bootstrap_eb.py failed!" 1>&2
    exit 1
  fi

  # create $EB_DIR
  mkdir -p $EB_DIR

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
}

# main

printf "\nInstalling Lua...\n"
lua_install
printf "Installing Lmod...\n"
lmod_install
printf "Bootstrapping EasyBuild...\n"
eb_bootstrap
exit 0
