#!/bin/bash

# Environment setup script for building RStudio-Server

version='2023.12.1+402'
toolchain='foss-2022b'

source_urls = ['https://github.com/rstudio/rstudio/archive']
sources = ['%(version)s.tar.gz']
patches = [
    # '%(name)s-2022.07.1+554_allow-disabling-quarto.patch',  patch failed
    # '%(name)s-2022.07.1+554_fix-libsoci-search.patch',
    # '%(name)s-2022.07.1+554_use-XDG_CACHE_HOME.patch',
    # '%(name)s-2023.09.0+463-yarn-bash.patch',
    # '%(name)s-2022.07.2+576-std-set.patch',
]
checksums = [
]

#  builddependencies
ml ant/1.10.11', '-Java-%(javaver)s', SYSTEM),
ml CMake/3.24.3'),
ml Ninja/1.11.1'),
ml pkgconf/1.9.3'),

# dependencies 
ml Boost', '1.81.0'),
ml Java', '11', '', SYSTEM),
ml R', '4.3.1'),
ml SOCI', '4.0.3'),
ml yaml-cpp/0.7.0

apt-get update -y
apt-get install libpam0g-dev

# rstudio-" + version.replace("+", "-")
fname_version=${version:+:-}
mkdir /app/software/RStudio-Server/${version}-${toolchain}/rstudio-" + version.replace("+", "-")

preconfigopts = " && ".join([
    # Install dependencies via scripts. Done in subshell to preserve PWD
    "(export RSTUDIO_TOOLS_ROOT='%(builddir)s'",
    "cd '%s/dependencies/common'" % local_start_dir,
    "./install-cef",
    "./install-dictionaries",
    "./install-mathjax",
    "./install-pandoc",
    "./install-packages",
    "./install-npm-dependencies)",
    ""
])

configopts = " ".join([
    "-DRSTUDIO_TOOLS_ROOT='%(builddir)s'",
    "-DRSTUDIO_TARGET=Server",
    "-DRSTUDIO_USE_SYSTEM_BOOST=ON",
    "-DRSTUDIO_USE_SYSTEM_SOCI=ON",
    "-DRSTUDIO_USE_SYSTEM_YAML_CPP=ON",
    "-DQUARTO_ENABLED=OFF",  # Not available on all archs, use pandoc fallback
    "-DRSTUDIO_GIT_REVISION_HASH=" + local_git_rev
])

sanity_check_commands = [
    # RSession requires environment variables R_HOME and R_DOC_DIR
    'R_HOME="$EBROOTR/lib64/R" R_DOC_DIR="$R_HOME/doc" rsession --verify-installation=1',
    # RServer requires a db conf (this may also be needed for live use)
    # Also create and set a soem dirs so it doesn't try to use $HOME
    'MYTMP=`mktemp -d`'
    ' && export RSTUDIO_CONFIG_DIR="$MYTMP"'
    ' && export XDG_DATA_HOME="$MYTMP/.data"'
    ' && export XDG_CACHE_HOME="$MYTMP/.cache"'
    ' && mkdir "$XDG_DATA_HOME" "$XDG_CACHE_HOME"'
    ' && export RS_LOG_DIR="$MYTMP/log"'
    ' && echo -e "provider=sqlite\\ndirectory=$MYTMP/db" >> "$MYTMP/db.conf"'
    ' && rserver ' + ' '.join([
        '--verify-installation=1',
        '--server-user="$USER"',
        '--database-config-file="$MYTMP/db.conf"',
        '--server-data-dir="$MYTMP/sdd"',
        '--secure-cookie-key-file="$MYTMP/secure-cookie-key"',
    ]),
]

sanity_check_paths = {
    'files': ['bin/rstudio-server'],
    'dirs': ['bin', 'extras', 'resources', 'www', 'www-symbolmaps', 'R'],
}

moduleclass = 'lang'
