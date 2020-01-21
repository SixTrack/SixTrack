#!/usr/bin/env bash
# Script to check that the singletrackfiles and the fort.90 files are matching

set -e

read90=$1
echo "Using read90 command='"$read90"'"
if [ "$read90" == "" ]; then
    echo "ERROR: No read90 command was specified"
    exit 1
fi
if [ ! -f $read90 ]; then
   echo "ERROR: read90 command '"$read90"' was not found"
   exit 1
fi
if [ ! -x $read90 ]; then
   echo "ERROR: read90 command '"$read90"' was not executable"
   exit 1
fi

if [ -d tmp ]; then
    echo "ERROR: Folder 'tmp' already exists in" $(pwd)
    exit 1
fi
mkdir tmp

for i in $(ls -d */); do
    if [ ! -f ${i}fort.3 ]; then
    #echo "Folder '"$i"' contains no tests; skipping."
    continue
    fi

    #echo "Checking folder '"${i%%/}"'";

    if [ ! -f ${i}fort.90.canonical ]; then
    echo "WARNING: No fort.90.canonical found in '"${i%%/}"'"
    continue
    fi

    if [ -f ${i}singletrackfile.dat.canonical ]; then
    $read90 --STF --SP 1 --fname ${i}singletrackfile.dat.canonical --ofname tmp/singletrackfile.dat.canonical.ascii
    $read90  --fname ${i}fort.90.canonical --ofname tmp/fort.90.canonical.ascii

    set +e
    diff -q tmp/singletrackfile.dat.canonical.ascii tmp/fort.90.canonical.ascii
    if [[ $? != 0 ]]; then
        echo "ERROR: Canonicals does not match in '"${i%%/}"'"
    fi
    set -e
    else
    echo "WARNING: No singletrackfile.dat.canonical found in '"${i%%/}"'"
    continue
    fi
done

rm -rf tmp
