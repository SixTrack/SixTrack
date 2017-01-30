#!/usr/bin/env bash
# This script can be used to automatically build the "one time build"
# libraries boinc and libarchive.

set -e #Exit on error

#Make sure we have the right submodule versions
if [[ $(uname) != MINGW* ]]; then     # Use MSYS on Windows, git on MINGW is buggy.
    #On LxPlus, you need to run git submodule init
    # from the toplevel of the working tree
    cd ..
    git submodule init
    git submodule update
    cd SixTrack
fi

####################################
### BOINC ##########################
####################################
cd boinc

if [[ $(uname) == FreeBSD* ]]; then
    MAKE=/usr/local/bin/gmake ./_autosetup -f
else
    ./_autosetup -f
fi

./configure --disable-client --disable-server --disable-manager --disable-boinczip

if [[ $(pwd) == /afs/* ]]; then
    #AFS doesn't like hardlinks between files in different directories and configure doesn't check for this corner case...
    sed -i 's/\/bin\/ln/cp/g' Makefile
    sed -i 's/\/bin\/ln/cp/g' api/Makefile
    sed -i 's/\/bin\/ln/cp/g' lib/Makefile
    #AFS doesn't like parallel make
    make
else
    # Machines with low memory doesn't like an automatic -j
    make -j 4
fi

#Need to build the boinc/api/boinc_api_fortran.o separately
cd api
make boinc_api_fortran.o
cd ..

cd ..
####################################
### libArchive #####################
####################################
rm -rf libarchive_build
mkdir libarchive_build
cd libarchive_build

if [[ $(uname) == MINGW* ]]; then
    #Windows build; use ICONV
    cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_BZip2=ON -DENABLE_ZLIB=ON -DENABLE_CAT=OFF -DENABLE_CPIO=OFF -DENABLE_EXPAT=OFF -DENABLE_INSTALL=OFF -DENABLE_LIBXML2=OFF -DENABLE_LZMA=OFF -DENABLE_NETTLE=OFF -DENABLE_OPENSSL=OFF -DENABLE_TAR=OFF -DENABLE_CNG=OFF -DENABLE_ICONV=ON -DENABLE_TEST=OFF -G "Unix Makefiles" ../libarchive -LH
else
    #Normal (UNIX) build.
    cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_BZip2=ON -DENABLE_ZLIB=ON -DENABLE_CAT=OFF -DENABLE_CPIO=OFF -DENABLE_EXPAT=OFF -DENABLE_INSTALL=OFF -DENABLE_LIBXML2=OFF -DENABLE_LZMA=OFF -DENABLE_NETTLE=OFF -DENABLE_OPENSSL=OFF -DENABLE_TAR=OFF -DENABLE_CNG=OFF -DENABLE_ICONV=OFF -DENABLE_TEST=OFF -G "Unix Makefiles" ../libarchive -LH
fi

if [[ $(pwd) == /afs/* ]]; then
    #AFS doesn't like parallel make
    make
else
    # Machines with low memory doesn't like an automatic -j
    make -j 4
fi

cd ..
