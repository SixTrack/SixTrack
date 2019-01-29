#!/usr/bin/env bash
# This script can be used to automatically build the "one time build"
# libraries boinc and libarchive.

set -e #Exit on error

echo ""
echo " Building zlib"
echo "==============="
echo

cd zlib
export ZLIB_BASE=$(pwd)

if [[ -d build ]]; then
   rm -rf build
fi
mkdir -p build
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
