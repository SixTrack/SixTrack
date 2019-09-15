#!/usr/bin/env bash
set -e
#set -v

NPAIRS=1000

gfortran mod_splitSingletrack.f90 splitSingletrack.f90 -o splitSingletrack

./splitSingletrack

for i in `seq 1 $NPAIRS`;
do
    i2=$(( 2*($i-1)+1 ))
    echo "pair#" $i " pIdx1=" $i2

    echo "./read90 --STF --fname singletrackfile.dat --SP $i"
    ./read90 --STF --fname singletrackfile.dat --SP $i2 > f1.txt

    echo "./read90 --fname singletrackfile.dat.`printf "%0*d" 6 $i`"
    ./read90 --fname singletrackfile.dat.`printf "%0*d" 6 $i`  > f2.txt

    echo "diff f1.txt f2.txt"
    diff f1.txt f2.txt
done
