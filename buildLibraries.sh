#!/usr/bin/env bash
# Script for building SixTrack dependencies that do not eed to be re-built every time SixTrack is built.

set -e # Exit on error

echo ""
echo " Building SixTrack Library Dependecies"
echo "========================================"
echo ""

ALL=true
BOINC=false
LIBARCH=false
HDF5=false
NAFF=false

for ARG in "$@"; do
    if [[ $ARG == "boinc" ]]; then
        BOINC=true
        LIBARCH=true
        echo "Boinc depends on libarchive, libarchive enabled as well."
    elif [[ $ARG == "libarchive" ]]; then
        LIBARCH=true
    elif [[ $ARG == "hdf5" ]]; then
        HDF5=true
    elif [[ $ARG == "naff" ]]; then
        NAFF=true
    else
        echo "Unknown library $ARG requested."
        exit 1
    fi
    echo "Will build $ARG"
    ALL=false
done

if [ $BOINC = true ] || [ $ALL = true ]; then
    git submodule init lib/boinc
    git submodule update lib/boinc
    cd lib
    ./buildBoinc.sh
    cd ..
fi

if [ $LIBARCH = true ] || [ $ALL = true ]; then
    git submodule init lib/libarchive
    git submodule update lib/libarchive
    cd lib
    ./buildLibarchive.sh
    cd ..
fi

if [ $HDF5 = true ] || [ $ALL = true ]; then
    cd lib
    ./buildHDF5.sh
    cd ..
fi

if [ $NAFF = true ] || [ $ALL = true ]; then
    cd lib
    ./buildNAFF.sh
    cd ..
fi
