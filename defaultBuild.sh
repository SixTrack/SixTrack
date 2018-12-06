#!/usr/bin/env bash

set -e #Exit on error

# Build default external libraries
./buildLibraries.sh naff libarchive

echo ""
echo " Building SixTrack with Defaults"
echo "================================="
echo ""

./cmake_six gfortran release LIBARCHIVE
