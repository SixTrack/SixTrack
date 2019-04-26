#!/usr/bin/env bash
# This script can be used to automatically build the "one time build"
# libraries boinc and libarchive.

set -e # Exit on error

echo ""
echo " Building libArchive"
echo "====================="
echo ""
echo " Linking to ZLib from: $ZLIB_BASE"
echo ""

if [[ -d libarchive_build ]]; then
    rm -rf libarchive_build
fi
mkdir libarchive_build
cd libarchive_build

cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_BZip2=OFF -DENABLE_ZLIB=ON -DENABLE_CAT=OFF -DENABLE_CPIO=OFF \
    -DENABLE_EXPAT=OFF -DENABLE_INSTALL=OFF -DENABLE_LIBXML2=OFF -DENABLE_LZMA=OFF -DENABLE_NETTLE=OFF \
    -DENABLE_OPENSSL=OFF -DENABLE_TAR=OFF -DENABLE_CNG=OFF -DENABLE_ICONV=OFF -DENABLE_TEST=OFF \
    -DZLIB_LIBRARY=$ZLIB_BASE/install/lib/libz.a -DZLIB_INCLUDE_DIR=$ZLIB_BASE/install/include \
    -G "Unix Makefiles" ../libarchive -LH

if [[ $(pwd) == /afs/* ]]; then
    # AFS doesn't like parallel make
    make
else
    # Machines with low memory doesn't like an automatic -j
    make -j6
fi
