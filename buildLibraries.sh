#!/usr/bin/env bash
# Script for building SixTrack dependencies that do not need to be re-built every time SixTrack is built.

echo ""
echo " Building SixTrack Library Dependecies"
echo "======================================="
echo ""

ALL=true
BOINC=false
HDF5=false
PYTHIA=false

for ARG in "$@"; do
    if [[ $ARG == "boinc" ]]; then
        BOINC=true
    elif [[ $ARG == "hdf5" ]]; then
        HDF5=true
    elif [[ $ARG == "pythia" ]]; then
        PYTHIA=true
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
