#!/bin/bash

#####
# # bootstrap script to demo easybuild
# # only tested on Ubuntu 14.04/16.04 
# # run as root inside fresh linux os 
# # or inside a container
######

# variables

# Internal

EB_DIR="/easybuild" # root folder to install easybuild into (can be a nfs mount)
EB_VER="3.7"      # version of EasyBuild to bootstrap in the container

# install these develop branches of easybuild-easyconfigs from github repos
EB_CFG_DEVELOP="hpcugent"  #EB_CFG_DEVELOP="hpcugent FredHutch"
# remove all older easyconfigs with these pattern
EB_OLDSTUFF=".*\(2014a\|2014b\|2015a\|2015b\|goolf\|ictce\|iimpi\|ifort\|icc-\|CrayGNU\|iomkl\|gimkl\).*.eb"

LUA_BASE_URL="http://www.lua.org/ftp/lua-"
LUAROCKS_BASE_URL="http://luarocks.org/releases/luarocks-"
LMOD_BASE_URL="https://github.com/TACC/Lmod/archive/"
EB_BASE_URL="https://raw.githubusercontent.com/easybuilders/easybuild-framework/master/easybuild/scripts/bootstrap_eb.py"

LUA_VER="5.4"   # verion of lua to install into the container
LUAROCKS_VER="3.0.3"   # version of luarocks package manager to install into the container
LMOD_VER="7.0"   # version of Lmod to install into the container

SOURCE_JAVA="http://ftp.osuosl.org/pub/funtoo/distfiles/oracle-java/jdk-8u92-linux-x64.tar.gz"

myself=$(whoami)
robot_paths=${EB_DIR}/github/easybuild-life-sciences/easybuild/easyconfigs
for clon in $EB_CFG_DEVELOP; do
  robot_paths=${robot_paths}:${EB_DIR}/github/develop/${clon}/easybuild/easyconfigs
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
  cd /tmp/Lmod-$LMOD_VER && ./configure && make install

  # link into /etc/profile so shells can use module function
  ln -s /usr/local/lmod/lmod/init/profile /etc/profile.d/modules.sh

  # add our easybuild modulepath now
  echo "export MODULEPATH=\$MODULEPATH:/${EB_DIR}/modules/all" | tee -a /etc/profile.d/modules.sh

  # source /etc/profile.d/modules.sh to enable modules in this shell
  # while this script is not intended to be sourced, it works to do so
  source /etc/profile.d/modules.sh
}

# install OS packages required by EasyBuild and foss toolchains
# encryption-related packages here on purpose as OS udpates should be more frequent
function install_EB_OS_pkgs {
  if hash apt-get 2>/dev/null; then
    apt-get update
    apt-get install -y wget python-minimal python-setuptools build-essential libibverbs-dev libssl-dev libffi-dev libreadline-dev unzip tcl git
    apt-get install -y python-pygraph
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
    wget -O /tmp/os-dependencies.apt https://raw.githubusercontent.com/FredHutch/easybuild-life-sciences/master/sandbox/os-dependencies.apt
    apt-get install -y pkg-config m4
    for mypkg in $(cat /tmp/os-dependencies.apt); do
      apt-get install -y $mypkg
    done
    # xorg-dev is bigger than libx11-dev and may not be needed.
    # libglu1-mesa-dev is needed for R rgl (R)
    # libcairo2-dev libxt-dev are needed for Cairo (R)
    # libpq-dev is needed for RPostgreSQL (R)
    # libnetcdf-dev is for netcdf4 (R)
    # libglpk-dev is for Rglpk (R)
    # RODBC rzmq
  elif hash yum 2>/dev/null; then
    #echo "redhat based install, not currently supported"
    wget -O /tmp/os-dependencies.yum https://raw.githubusercontent.com/FredHutch/easybuild-life-sciences/master/sandbox/os-dependencies.yum
    for mypkg in $(cat /tmp/os-dependencies.yum); do
      yum -y install $mypkg
    done
  else
    echo "unknown unix, not currently supported"
  fi
}

# remove OS packages after dependency install to get clean OS
function remove_OS_pkgs {
  if hash apt-get 2>/dev/null; then
    apt-get remove -y libreadline-dev
  elif hash yum 2>/dev/null; then
    echo "redhat based install, not currently supported"
  else
    echo "unknown unix, not currently supported"
  fi
}

# download some required stuff and to it into source 
function download_extra_sources {

  /bin/su - eb -c "mkdir -p $EB_DIR/sources"
  /bin/su - eb -c "wget -P \"${EB_DIR}/sources\" \"${SOURCE_JAVA}\""

  /bin/su - eb -c "mkdir -p $EB_DIR/github"
  /bin/su - eb -c "git clone https://github.com/FredHutch/easybuild-life-sciences ${EB_DIR}/github/easybuild-life-sciences"

  /bin/su - eb -c "mkdir -p $EB_DIR/github/develop"
  for clon in $EB_CFG_DEVELOP; do 
    if ! [[ -d ${EB_DIR}/github/develop/${clon} ]]; then
      /bin/su - eb -c "git clone -b develop --single-branch https://github.com/${clon}/easybuild-easyconfigs ${EB_DIR}/github/develop/${clon}"
      if [[ -n $EB_OLDSTUFF ]]; then
        printf "deleting old easyconfigs in ${clon} ...\n\n"
        find "${EB_DIR}/github/develop/${clon}" -regex "${EB_OLDSTUFF}" -exec rm "{}" \;
      fi
    fi
  done

}

# bootstrap easybuild
function eb_bootstrap {

  # create eb user if needed
  if id "eb" >/dev/null 2>&1; then
    echo "user eb exists"
  else
    adduser --disabled-password --gecos "" eb
  fi

  # get EB
  eb_url=$(printf "$EB_BASE_URL" "$EB_VER")
  /bin/su - eb -c "wget -O /tmp/bootstrap_eb.py $eb_url"
  if [ "$?" != "0" ]; then
    echo "Oops, retrieving bootstrap_eb.py failed!" 1>&2
    exit 1
  fi

  # create $EB_DIR
  mkdir -p "${EB_DIR}"
  chown -R eb ${EB_DIR}

  # bootstrap it
  /bin/su - eb -c "/usr/bin/python /tmp/bootstrap_eb.py $EB_DIR"

  # pop in some useful environment variables to our EasyBuild modulefile
  (/bin/su - eb -c cat <<EOF
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

  if [[ -n $EB_OLDSTUFF ]]; then    
    easyconfigs="${EB_DIR}/software/EasyBuild/${EB_VER}/lib/python2.7/site-packages/easybuild_easyconfigs-${EB_VER}-py2.7.egg/easybuild/easyconfigs"
    printf "deleting old easyconfigs in ${easyconfigs} ...\n\n"
    find "${easyconfigs}" -regex "${EB_OLDSTUFF}" -exec rm "{}" \;
  fi

}

# main

# checking requirements
mem=$(awk '( $1 == "MemTotal:" ) { printf "%.0f", $2/1024/1024 }' /proc/meminfo)
if [[ $mem -lt 8 ]]; then
  printf "This script requires a minimum of 8GB ram, 8 GB disk and 8 cores !\n"
  exit 1
fi
if ! [[ "$myself" == "root" ]]; then 
  printf "  This script needs to be run as root.\n"
  exit 1
fi
if [[ "$(ls -A $EB_DIR)" ]]; then
  echo "$EB_DIR is not empty. Please empty it or choose a different folder."
  exit 1
fi

printf "This script will install Lua & Lmod from source in /usr/local,\n"
printf "and add configuration /etc/profile.d/modules.sh\n"

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
printf "Login shells now have modules enabled, so all you need to do is login as user 'eb' or switch to the eb user 'su - eb'.\n"
printf "Once you are back, you should have modules, and can run 'ml EasyBuild' to get started building!\n"
exit 0
