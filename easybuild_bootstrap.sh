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
LMOD_VER="7.0"   # version of Lmod to install into the container

EB_VER="3.0.0"   # verison of EasyBuild to bootstrap in the container
EB_DIR="/easybuild"   # location for EasyBuild directory tree

EB_CFG_DEVELOP="hpcugent FredHutch"
EB_TOOLCHAIN_ONLY="foss-2016b"

SOURCE_JAVA="http://ftp.osuosl.org/pub/funtoo/distfiles/oracle-java/jdk-8u92-linux-x64.tar.gz"

myself=$(whoami)
for clon in $EB_CFG_DEVELOP; do
  robot_paths=${EB_DIR}/github/develop/${clon}/easybuild/easyconfigs:${robot_paths}
done

# install lua, luarocks, luafilesystem, and luaposix
function lua_install {

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

}

# install Lmod
function lmod_install {

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

  # add our easybuild modulepath now
  echo "export MODULEPATH=\$MODULEPATH:/${EB_DIR}/modules/all" | sudo tee -a /etc/profile.d/modules.sh

  # source /etc/profile.d/modules.sh to enable modules in this shell
  # while this script is not intended to be sourced, it works to do so
  source /etc/profile.d/modules.sh
}

# install OS packages required by EasyBuild and foss toolchains
# encryption-related packages here on purpose as OS udpates should be more frequent
function install_EB_OS_pkgs {
  if hash apt-get 2>/dev/null; then
    sudo apt-get update
    sudo apt-get install -y python-minimal python-pygraph build-essential libibverbs-dev libssl-dev libffi-dev libreadline-dev unzip tcl git
  elif hash yum 2>/dev/null; then
    echo "redhat based install, not currently supported"
  else
    echo "unknown unix, not currently supported"
  fi
}

# install OS packages to satisfy unstated dependencies in common easyconfigs
# this is usually due to Ubuntu<->RedHat differences and these should be included in easyconfigs eventually
function install_missed_dependency_OS_pkgs {
  if hash apt-get 2>/dev/null; then
    sudo apt-get install -y pkg-config m4 libx11-dev
  elif hash yum 2>/dev/null; then
    echo "redhat based install, not currently supported"
  else
    echo "unknown unix, not currently supported"
  fi
}

# remove OS packages after dependency install to get clean OS
function remove_OS_pkgs {
  if hash apt-get 2>/dev/null; then
    sudo apt-get remove -y libreadline-dev
  elif hash yum 2>/dev/null; then
    echo "redhat based install, not currently supported"
  else
    echo "unknown unix, not currently supported"
  fi
}

# download some required stuff and to it into source 
function download_extra_sources {

  mkdir -p $EB_DIR/sources
  wget -P "${EB_DIR}/sources" "${SOURCE_JAVA}"

  mkdir -p $EB_DIR/github/develop
  for clon in $EB_CFG_DEVELOP; do 
    if ! [[ -d ${EB_DIR}/github/develop/${clon} ]]; then
      git clone -b develop --single-branch https://github.com/${clon}/easybuild-easyconfigs ${EB_DIR}/github/develop/${clon}
      if [[ -n $EB_TOOLCHAIN_ONLY ]]; then
        find ${EB_DIR}/github/develop/${clon} ! -name *-${EB_TOOLCHAIN_ONLY}-* -name *.eb -exec rm {} \;
      fi
    fi
  done

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
  sudo mkdir -p "${EB_DIR}"
  sudo chown -R ${myself} ${EB_DIR}
  homedir=~
  sudo chown -R ${myself} ${homedir}

  # bootstrap it
  python /tmp/bootstrap_eb.py $EB_DIR

  # pop in some useful environment variables to our EasyBuild modulefile
  (cat <<EOF
setenv ("EASYBUILD_SOURCEPATH", "${EB_DIR}/sources")
setenv ("EASYBUILD_BUILDPATH", "${EB_DIR}/build")
setenv ("EASYBUILD_INSTALLPATH_SOFTWARE", "${EB_DIR}/software")
setenv ("EASYBUILD_INSTALLPATH_MODULES", "${EB_DIR}/modules")
setenv ("EASYBUILD_REPOSITORYPATH", "${EB_DIR}/ebfiles_repo")
setenv ("EASYBUILD_ROBOT_PATHS", "${robot_paths}")
setenv ("EASYBUILD_LOGFILE_FORMAT", "${EB_DIR}/logs,easybuild-%(name)s-%(version)s-%(date)s.%(time)s.log")
setenv ("EASYBUILD_MODULES_TOOL", "Lmod")
EOF
  ) >> $EB_DIR/modules/all/EasyBuild/$EB_VER.lua

  if [[ -n $EB_TOOLCHAIN_ONLY ]]; then
    easyconfigs="${EB_DIR}/software/EasyBuild/${EB_VER}/lib/python2.7/site-packages/easybuild_easyconfigs-${EB_VER}-py2.7.egg/easybuild/easyconfigs"
    find ${easyconfigs} ! -name *-${EB_TOOLCHAIN_ONLY}-* -name *.eb -exec rm {} \;
  fi

}

# main

# checking requirements
mem=$(awk '( $1 == "MemTotal:" ) { printf "%.0f", $2/1024/1024 }' /proc/meminfo)
if [[ $mem -lt 8 ]]; then
  printf "This script requires a minimum of 8GB ram, 8 GB disk and 8 cores !\n"
  exit 1
fi
if [[ "$myself" == "root" ]]; then 
  printf "please do not run as root, run under a user with sudo access"
  exit 1
fi

printf "Installing required OS pkgs...\n"
install_EB_OS_pkgs
printf "Installing Lua...\n"
lua_install
printf "Installing Lmod...\n"
lmod_install
printf "Removing packages to clean system after dependency install...\n"
remove_OS_pkgs
printf "Installing OS pkgs to satify missing dependencies in easyconfigs...\n"
install_missed_dependency_OS_pkgs
printf "Bootstrapping EasyBuild...\n"
eb_bootstrap
printf "Downloading some extras...\n"
download_extra_sources
printf "\n\n"
printf "************It appears to have worked! *************\n"
printf "Login shells now have modules enabled, so all you need to do is log out and back in.\n"
printf "   ... or you can source /etc/profile.d/modules.sh .\n"
printf "Once you are back, you should have modules, and can run 'module load EasyBuild' to get started building!\n"
exit 0

