#!/usr/bin/env bash

set -e #Exit on error

echo ""
echo " Building Standard SixTrack"
echo "============================"
echo ""

./buildLibraries.sh naff
./cmake_six gfortran release
