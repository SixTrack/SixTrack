#!/usr/bin/env bash

set -e #Exit on error

echo ""
echo " Building HDF5"
echo "==============="
echo ""

VERSION=1.10
PATCH=1.10.1

mkdir -pv hdf5
cd hdf5

# Download links from https://support.hdfgroup.org/HDF5/release/cmakebuild.html
if [[ $(uname) == MINGW* ]]; then
  if [ ! -f CMake-hdf5-1.10.1.zip ]; then
    SOURCE=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$VERSION/hdf5-$PATCH/src/CMake-hdf5-$PATCH.zip
    wget $SOURCE
    unzip CMake-hdf5-$PATCH.zip
  fi
else
  if [ ! -f CMake-hdf5-1.10.1.tar.gz ]; then
    SOURCE=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$VERSION/hdf5-$PATCH/src/CMake-hdf5-$PATCH.tar.gz
    wget $SOURCE
    tar -xf CMake-hdf5-$PATCH.tar.gz
  fi
fi

if [ -d build ]; then
  rm -rf build
fi

mkdir build
cd build
cmake ../CMake-hdf5-$PATCH/hdf5-$PATCH -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DBUILD_SHARED_LIBS=OFF \
  -DHDF5_ENABLE_PARALLEL=OFF -DHDF5_BUILD_CPP_LIB=OFF -DBUILD_TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF \
  -DH5_ZLIB_HEADER=$(pwd)/../../zlib/install/include/zlib.h \
  -DZLIB_INCLUDE_DIRS=$(pwd)/../../zlib/install/include/ \
  -DZLIB_STATIC_LIBRARY=$(pwd)/../../zlib/install/lib/libz.a
make -j4
cd ..
cd ..
