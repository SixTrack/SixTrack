#!/usr/bin/env bash
# K. Sjobak, 2016
# Script for building SixTrack with all know supported compiler- and flag combinations.

#Sucess:
make                  # SixTrack_4534_cernlib_crlibm_fast_gfortran_m32_tilt_vvector_O3
make add=m64          # SixTrack_4534_crlibm_fast_gfortran_m64_tilt_vvector_O3
make add=collimat     # SixTrack_4534_cernlib_collimat_fast_gfortran_m32_tilt_vvector_O3
make add=da           # SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_tilt_vvector_O3


#Fails:
make add=hdf5         # SixTrack_4534_cernlib_collimat_fast_gfortran_hdf5_m32_tilt_vvector_O3
                      # makefile:449: recipe for target 'SixTrack_4534_cernlib_collimat_fast_gfortran_hdf5_m32_tilt_vvector_O3/beamgas.o' failed

make add="da hdf5"    # SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3
                      # SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixtrack
                      # SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.o: In function `coruglo_':
                      # /home/kyrsjo/code/SixTrack/SixTrack/SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.f:1507: undefined reference to `x04abf_'
                      # /home/kyrsjo/code/SixTrack/SixTrack/SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.f:1621: undefined reference to `e04uef_'
                      # /home/kyrsjo/code/SixTrack/SixTrack/SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.f:1625: undefined reference to `e04udm_'
                      # /home/kyrsjo/code/SixTrack/SixTrack/SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.f:1625: undefined reference to `e04ucf_'
                      # SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.o: In function `coruord_':
                      # /home/kyrsjo/code/SixTrack/SixTrack/SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.f:469: undefined reference to `x04abf_'
                      # /home/kyrsjo/code/SixTrack/SixTrack/SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.f:592: undefined reference to `e04uef_'
                      # /home/kyrsjo/code/SixTrack/SixTrack/SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.f:596: undefined reference to `e04udm_'
                      # /home/kyrsjo/code/SixTrack/SixTrack/SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixda.f:596: undefined reference to `e04ucf_'
                      # collect2: error: ld returned 1 exit status
                      # makefile:369: recipe for target 'SixTrack_4534_cernlib_crlibm_da_fast_gfortran_m32_naglib_tilt_vvector_O3/sixtrack' failed

