#!/usr/bin/env bash
# Input arguments:
# $1 splitSingetrack binary
# $2 read90 binary
# $3 singletrackfile to use

set -e
#set -v

if [ "$#" -ne 3 ]; then
    echo "Expected some arguments!"
    exit 2
fi

cp $3 singletrackfile.dat

NPAIRS=`$1 --getNumPairs`
echo NPAIRS = $NPAIRS

$1 #Run the splitter, no arguments

for i in `seq 1 $NPAIRS`;
do
    i2=$(( 2*($i-1)+1 ))
    echo "pair#" $i " pIdx1=" $i2

    echo "read90 --STF --fname singletrackfile.dat --SP $i"
    $2  --STF --fname singletrackfile.dat --SP $i2 > f1.txt

    echo "./read90 --fname singletrackfile.dat.`printf "%0*d" 6 $i`"
    $2 --fname singletrackfile.dat.`printf "%0*d" 6 $i`  > f2.txt

    echo "diff f1.txt f2.txt"
    diff f1.txt f2.txt
done
