!#/bin/bash

name=oarfish
version='0.6.5'

toolchain=GCCcore-13.2.0


ml Rust/1.76.0-GCCcore-13.2.0
ml binutils/2.40-GCCcore-13.2.0
ml pkgconf/2.0.3-GCCcore-13.2.0
cargo-download/v0.1.2-GCCcore-13.2.0

DIR="/app/software/${name}/${version}-$toolchain"

mkdir -p ${DIR}
cd ${DIR}
export CARGO_HOME=$DIR

cargo install ${name}

#  URL=https://crates.io/api/v1/crates/${name}/${version}/download
#  curl -L https://crates.io/api/v1/crates/oarfish/0.6.5/download --output ${name}-${version}.tar.gz 
#  tar -xf ${name}-${version}.tar.gz 
#  cd ${name}-${version}

#  Other Cargo Vars
# CARGO_INSTALL_ROOT

### Build Issues
add: #![feature(c_str_literals)]
to: registry/src/index.crates.io-6f17d22bba15001f/oarfish-0.6.5/src/main.rs

add: #![feature(iter_repeat_n)] 
to: registry/src/index.crates.io-6f17d22bba15001f/noodles-bam-0.70.0/src/lib/.rs

