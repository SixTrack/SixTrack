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
    # These numbers will need updating in the future.
    AUTOCONF_VERSION=2.69 AUTOMAKE_VERSION=1.16 MAKE=/usr/local/bin/gmake ./_autosetup -f
elif [[ $(uname) == NetBSD* ]]; then
    MAKE=/usr/pkg/bin/gmake ./_autosetup -f
elif [[ $(uname) == Darwin* ]]; then
    LIBTOOLIZE=/usr/local/bin/glibtoolize ./_autosetup -f
elif [[ $(uname) == SunOS* ]]; then
    MAKE=/usr/bin/gmake ./_autosetup -f
else
    ./_autosetup -f
fi

./configure --disable-client --disable-server --disable-manager --disable-boinczip

if [[ $(pwd) == /afs/* ]]; then
    # AFS doesn't like hardlinks between files in different directories and configure doesn't check for this corner case...
    sed -i 's/\/bin\/ln/cp/g' Makefile
    sed -i 's/\/bin\/ln/cp/g' api/Makefile
    sed -i 's/\/bin\/ln/cp/g' lib/Makefile
    # AFS doesn't like parallel make
    make
else
    # Machines with low memory doesn't like an automatic -j
    # dmake on solaris needs a space between j and the number.
    make -j 4
fi

# Need to build the boinc/api/boinc_api_fortran.o separately
cd api
if [[ -e boinc_api_fortran.o ]]; then
    # In case we already have such an .o file, it may be of the wrong arch, so recompile.
    rm boinc_api_fortran.o
fi
make boinc_api_fortran.o
