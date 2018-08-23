#!/usr/bin/env bash
# Script for building SixTrack dependencies that do not eed to be re-built every time SixTrack is built.

set -e # Exit on error

echo ""
echo " Building SixTrack Library Dependecies"
echo "========================================"
echo ""

# Make sure we have the right submodule versions
echo "Updating submodules ..."
git submodule init
git submodule update
echo "Done"

cd lib

# Run individual build scripts
./buildBoinc.sh
./buildLibarchive.sh

cd ..