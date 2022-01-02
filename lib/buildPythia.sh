#!/usr/bin/env bash

VERSION=pythia8243

echo ""
echo " Building PYTHIA support for SixTrack "
echo "======================================"
echo ""

mkdir -pv pythia
cd pythia
rm -rf pythia*/

if [ ! -e $VERSION.tgz ]; then
  echo ""
  echo "Downloading:"
  wget https://pythia.org/download/pythia82/$VERSION.tgz
fi

echo ""
echo "Extracting:"
tar -xf $VERSION.tgz --totals

echo ""
echo "Building:"
cd $VERSION
./configure --enable-64bit --enable-shared
make

echo ""
echo "Making Symlinks:"
cd ..
ln -sfv $VERSION/lib     lib
ln -sfv $VERSION/include include
cd ..

echo ""
echo "Done!"
echo ""
