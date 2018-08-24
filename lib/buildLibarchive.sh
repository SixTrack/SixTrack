#!/usr/bin/env bash
# This script can be used to automatically build the "one time build"
# libraries boinc and libarchive.

set -e #Exit on error

echo ""
echo " Building libArchive and zlib"
echo "=============================="
echo
echo "**** ZLIB ****"
echo
#Support library zlib
cd zlib
ZLIB_BASE=$(pwd)
if [[ -d build ]]; then
   rm -rf build
fi
mkdir build
if [[ -d install ]]; then
   rm -rf install
fi
mkdir install

cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$ZLIB_BASE/install -DCMAKE_C_FLAGS=-fPIC -G "Unix Makefiles"
make
make install
cd ..

if [[ $(uname) == MINGW* ]]; then
    #Windows build creates a differently-named libz.a.
    cd install/lib
    cp libzlibstatic.a libz.a
    cd ../..
fi
ZLIB_PATH=$ZLIB_BASE/install/lib/libz.a

cd ..

echo
echo "**** libArchive ****"
echo

#Then: LibArchive itself
if [[ -d libarchive_build ]]; then
    rm -rf libarchive_build
fi
mkdir libarchive_build
cd libarchive_build

cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_BZip2=OFF -DENABLE_ZLIB=ON -DENABLE_CAT=OFF -DENABLE_CPIO=OFF -DENABLE_EXPAT=OFF -DENABLE_INSTALL=OFF -DENABLE_LIBXML2=OFF -DENABLE_LZMA=OFF -DENABLE_NETTLE=OFF -DENABLE_OPENSSL=OFF -DENABLE_TAR=OFF -DENABLE_CNG=OFF -DENABLE_ICONV=OFF -DENABLE_TEST=OFF -DZLIB_LIBRARY=$ZLIB_PATH -DZLIB_INCLUDE_DIR=$ZLIB_BASE/install/include -G "Unix Makefiles" ../libarchive -LH

if [[ $(pwd) == /afs/* ]]; then
    #AFS doesn't like parallel make
    make
else
    # Machines with low memory doesn't like an automatic -j
    make -j 4
fi
