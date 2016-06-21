#!/usr/bin/env bash
# K. Sjobak, 2016
# Script for building SixTrack with all know supported compiler- and flag combinations.

#Sucess:
echo "Defaults:"
make                     # SixTrack_4534_cernlib_crlibm_fast_gfortran_m32_tilt_vvector_O3

#Adds:
echo
echo "m64 -> removes cernlib"
make add=m64             # SixTrack_4534_crlibm_fast_gfortran_m64_tilt_vvector_O3

echo
echo "collimat"
make add=collimat        # SixTrack_4534_cernlib_collimat_fast_gfortran_m32_tilt_vvector_O3

echo
echo "collimat m64"
make add="collimat m64"  # SixTrack_4534_cernlib_collimat_fast_gfortran_m32_tilt_vvector_O3

echo
echo "da"
make add=da              # SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_tilt_vvector_O3

echo
echo "da m64"
make add=da              # SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_tilt_vvector_O3

echo
echo "beamgas -> adds collimat"
make add="beamgas"

echo
echo "bpm"
make add="bpm"


echo
echo "hdf5"
make add=hdf5            # SixTrack_4534_cernlib_collimat_fast_gfortran_hdf5_m32_tilt_vvector_O3
                         # makefile:449: recipe for target 'SixTrack_4534_cernlib_collimat_fast_gfortran_hdf5_m32_tilt_vvector_O3/beamgas.o' failed

#Removes:
echo
echo "-cernlib"
make remove=cernlib

echo
echo "-cernlib -crlibm"
make remove="cernlib crlibm"




#Fails (and should fail):
echo
echo "da hdf5 (not actually compatible)"
make add="da hdf5"       

echo
echo "naglib"
make add="naglib"

#Fails (bug):
echo
echo "da m64 naglib"
make add="da m64 naglib" #Fails because naglib can't be linked. Works in make_six...

