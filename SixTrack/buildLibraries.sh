#!/usr/bin/env bash
# This script can be used to automatically build the "one time build"
# libraries boinc and libarchive.

set -e #Exit on error

#Make sure we have the right submodule versions
if [[ $(uname) != MINGW* ]]; then
    # Use MSYS on Windows, git on MINGW is buggy.
    git submodule init
    git submodule update
fi

### BOINC ###
cd boinc

./_autosetup
./configure --disable-client --disable-server --disable-manager --disable-boinczip

make -j 4 # Machines with low memory doesn't like an automatic -j

cd api
make boinc_api_fortran.o
cd ..

cd ..

### libArchive ###
rm -rf libarchive_build
mkdir libarchive_build
cd libarchive_build

if [[ $(uname)=!= MINGW* ]]; then
    #Windows build; use ICONV
    cmake -DENABLE_BZip2=ON -DENABLE_ZLIB=ON -DENABLE_CAT=OFF -DENABLE_CPIO=OFF -DENABLE_EXPAT=OFF -DENABLE_INSTALL=OFF -DENABLE_LIBXML2=OFF -DENABLE_LZMA=OFF -DENABLE_NETTLE=OFF -DENABLE_OPENSSL=OFF -DENABLE_TAR=OFF -DENABLE_CNG=OFF -DENABLE_ICONV=ON -DENABLE_TEST=OFF -G "Unix Makefiles" ../libarchive -LH
else
    #Normal (UNIX) build.
    cmake -DENABLE_BZip2=ON -DENABLE_ZLIB=ON -DENABLE_CAT=OFF -DENABLE_CPIO=OFF -DENABLE_EXPAT=OFF -DENABLE_INSTALL=OFF -DENABLE_LIBXML2=OFF -DENABLE_LZMA=OFF -DENABLE_NETTLE=OFF -DENABLE_OPENSSL=OFF -DENABLE_TAR=OFF -DENABLE_CNG=OFF -DENABLE_ICONV=OFF -DENABLE_TEST=OFF -G "Unix Makefiles" ../libarchive -LH
fi

make -j 4 # Machines with low memory doesn't like an automatic -j

cd ..
