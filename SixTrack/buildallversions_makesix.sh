#!/bin/sh
# Eric McIntosh, 2016
# Script for building SixTrack with all known supportedg
# compiler- and flag combinations.

# ADD a build of bignblz and _da hugenblz
for C in ifort
do
  for A in bnlelens
  do
    for B in "" "boinc" "boinc api"
    do
      for SSE in "" sse2 sse3
      do
        make_six -cernlib $C $A $B $SSE </dev/null
      done
    done
  done
done
for C in gfortran
do
  for A in bnlelens
  do
# Add the api option when Windows gfortran libraries are ready.
    for B in "" "boinc"
    do
      make_six -cernlib $C $A $B </dev/null
    done
  done
done
for C in nagfor
do
  for A in bnlelens
  do
    for B in "" "boinc"
    do
      make_six -cernlib $C $A $B </dev/null
    done
  done
done 
