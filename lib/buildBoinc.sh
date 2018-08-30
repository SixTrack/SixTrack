#!/usr/bin/env bash
# This script can be used to automatically build the "one time build"
# libraries boinc and libarchive.

set -e #Exit on error

echo ""
echo " Building BOINC"
echo "================"
echo ""

cd boinc

if [[ $(uname) == FreeBSD* ]]; then
    MAKE=/usr/local/bin/gmake ./_autosetup -f
elif [[ $(uname) == OpenBSD* ]]; then
#These numbers will need updating in the future.
    AUTOCONF_VERSION=2.69 AUTOMAKE_VERSION=1.15 MAKE=/usr/local/bin/gmake ./_autosetup -f
elif [[ $(uname) == NetBSD* ]]; then
    MAKE=/usr/pkg/bin/gmake ./_autosetup -f
elif [[ $(uname) == Darwin* ]]; then
    LIBTOOLIZE=/usr/local/bin/glibtoolize ./_autosetup -f
else
    ./_autosetup -f
fi

./configure --disable-client --disable-server --disable-manager --disable-boinczip

# This is a terrible hack for building on MinGW, but it works.
# Line numbers may need to be updated if BOINC is updated.
if [[ $(uname) == MINGW* ]]; then
    cd lib

    if [ ! -f "boinc_win.h.bak" ]; then
        mv boinc_win.h boinc_win.h.bak
        cat boinc_win.h.bak | head -n27 > boinc_win.h
        echo "#include \"windows.h\"" >> boinc_win.h
        cat boinc_win.h.bak | tail -n+28 >> boinc_win.h
    fi

    if [ ! -f "util.cpp.bak" ]; then
        mv util.cpp util.cpp.bak
        cat util.cpp.bak | head -n631 > util.cpp
        echo "int get_real_executable_path(char* , size_t ) {return ERR_NOT_IMPLEMENTED;}" >> util.cpp
    fi

    cd ..
fi

if [[ $(pwd) == /afs/* ]]; then
    #AFS doesn't like hardlinks between files in different directories and configure doesn't check for this corner case...
    sed -i 's/\/bin\/ln/cp/g' Makefile
    sed -i 's/\/bin\/ln/cp/g' api/Makefile
    sed -i 's/\/bin\/ln/cp/g' lib/Makefile
    #AFS doesn't like parallel make
    make
else
    # Machines with low memory doesn't like an automatic -j
    make -j4
fi

# Need to build the boinc/api/boinc_api_fortran.o separately
cd api
if [[ -e boinc_api_fortran.o ]]; then
    # In case we already have such an .o file, it may be of the wrong arch, so recompile.
    rm boinc_api_fortran.o
fi
make boinc_api_fortran.o
