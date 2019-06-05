#!/usr/bin/env bash
# Script for building SixTrack dependencies that do not need to be re-built every time SixTrack is built.

echo ""
echo " Building SixTrack Library Dependecies"
echo "======================================="
echo ""

ALL=true
BOINC=false
LIBARCH=false
ZLIB=false
HDF5=false
PYTHIA=false
NAFF=false

for ARG in "$@"; do
    if [[ $ARG == "boinc" ]]; then
        BOINC=true
    elif [[ $ARG == "libarchive" ]]; then
        LIBARCH=true
        ZLIB=true
        echo "Libarchive depends on zlib, zlib enabled as well."
    elif [[ $ARG == "hdf5" ]]; then
        HDF5=true
        ZLIB=true
        echo "HDF5 depends on zlib, zlib enabled as well."
    elif [[ $ARG == "pythia" ]]; then
        PYTHIA=true
    elif [[ $ARG == "naff" ]]; then
        NAFF=true
    elif [[ $ARG == "zlib" ]]; then
        ZLIB=true
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

# If building libArchive or HDF5, ZLib must be built first!
if [ $ZLIB = true ] || [ $ALL = true ]; then
    cd lib
    source ./buildZlib.sh
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

if [ $PYTHIA = true ] || [ $ALL = true ]; then
    cd lib
    ./buildPythia.sh
    cd ..
fi

if [ $NAFF = true ] || [ $ALL = true ]; then
    cd lib
    ./buildNAFF.sh
    cd ..
fi
