list='ant
  build-essential
  clang
  cmake
  debsigs
  dpkg-sig
  expect
  fakeroot
  gnupg1
  libacl1-dev
  libattr1-dev
  libbz2-dev
  libcap-dev
  libclang-6.0-dev
  libclang-dev
  libcurl4-openssl-dev
  libegl1-mesa
  libfuse2
  libgl1-mesa-dev
  libgtk-3-0
  libpam-dev
  libpango1.0-dev
  libssl-dev
  libuser1-dev
  libxslt1-dev
  lsof
  patchelf
  pkg-config
  python
  rrdtool
  software-properties-common
  unzip
  uuid-dev
  wget
  zlib1g-dev'

dpkg -l >pkg-list
for pkg in $list; do 
    echo package: $pkg
    grep $pkg pkg-list
done
