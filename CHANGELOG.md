# SixTrack Changelog

### Version 4.7.18 [21.12.2017]

* Fixing the RF multipole elements
  - RF multipoles now affect the energy directly, not the momentum
  - Fix missing RV update in RF multipole element with FOX (DA package); this improves the symplecticity check
  - Fix errorenous factor 10^-3 in FOX part for RF mulitpole
  - Fix broken 2nd order RF multipole in FOX part (missing `+ca alignf` call)
  - Take care of px vs. xp. coordinates
  - NOTE: Compilation with TILT flag OFF is broken.
* Fix the makefile for build_manual
* Fix building the testing wrapper on some older systems by explicitly enable C++11, even if the compiler does not use this by default
* Update DA version tests and use explicit write formats for the DA output

### Version 4.7.17 [11.11.2017]

* Fix CR for DUMP format 7 and 8.

### Version 4.7.16 [06.11.2017]

* Added option in collimation block for reading normalized distribution
* Fix an ifort compiler warning in SCATTER.

### Version 4.7.15 [02.11.2017]

* Fixing the collimat 10k turn limit
* SCATTER:
  - Adds offset to Gaussian beam profile.
  - Some minor changes to floating point values in calculations.

### Version 4.7.14 [31.10.2017]

* Remove boost dependency in NAFF by including the required header in the distribution + some small changes.

### Version 4.7.13 [27.10.2017]

* Round properly when reading fort.13 when CRLIBM is active.

### Version 4.7.12 [26.10.2017]

* Added descriptive(-ish) headers to subroutines like in the DYNK module
* Added the input parsing and calculation for Gaussian profile

### Version 4.7.11 [22.10.2017]

* Added detection and error message for wrong LINEAR OPTICS setup when using collimation
* Document the DIPEDGE element - now we know it is working, including conversion from MadX.
* Fix building of the DA version
* Fix linking to NAGLIB
* Fix building with NAGFOR
* Improve build speeds on LXPLUS when not doing BUILD_TESTING
* Let Astuce++ build with the same flags as the rest of the project (i.e. build it static for a static SixTrack build, dynamic for a dynamically linked SixTrack build)
* Fix building on recent versions of MSYS2
* Fix bulding on newer versions of OS X / Homebrew
* Update a Sixin.zip that was forgotten
* Clean up the build of ZLIB, fixing building on GNU HURD

### Version 4.7.10 [03.10.2017]

* DYNK: Rewrite from FixedForm to FreeForm FORTRAN and rename dynk.s to dynk.s90
* DYNK: Allocate memory dynamically
* DYNK: Write only what is needed to the CR files
* DYNK: Various tweaks in parseFUN
* FMA/PLATO: Fix compilation with gfortran 4.4.7

### Version 4.7.9 [20.09.2017]

* FMA: Adding support for setting first and last turn of analysis, independent of DUMP
* FMA: Write the NORM_ dump files once and only once
* PLATO is now a f90 module
* Resurrect the -CRLIBM version, it was broken in master due to a missing "then"
* Cleanup FMA in physical coordinates

### Version 4.7.8 [14.09.2017]

* Add experimental SCATTER block
* Build system improvements for handling multiple .s and .s90 files
* Split DYNK and BEAMGAS into separate files

### Version 4.7.7 [14.09.2017]

* Normalized DUMPs (normalized particle coordinates in ASCII & binary, normalized matrix)
* Read StartDUMP in FMA
* Possibility to do 4D FMA at StartDUMP
* Change constant "1e-3" to "c1m3=1d-3", used when setting the right units in the matrices used for normalization. Update NORM_DUMP canonicals & fma_sixtrack canonicals.

### Version 4.7.6 [01.09.2017]

* Integrate NAFF library for FMA analysis

### Version 4.7.5 [25.08.2017]

* Fix stringlength from stringzerotrim; This fixes a bug introduced in v4.7.3 (PR #310).

### Version 4.7.4 [18.08.2017]

* Allow FMA analysis to read binary DUMP files.
* Add NoNormDUMP flag to disable writing of normalized dump files from FMA.

### Version 4.7.3 [16.08.2017]

* Special DUMP at injection point.

### Version 4.7.2 [13.07.2017]

* Compare machine length in SYNC block (TLEN) with calculated value and print a warning if they are different.

### Version 4.7.1 [07.07.2017]

* Fix output formats to avoid ifort warnings like `"Ifort: remark #8291: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+7'"`
* Adding support for gprof with the GPROF flag to our CMakelists, various cleanup.

## Version 4.7 [03.07.2017]

* Split off the collimation code into a new file `collimation.s`
* Split off `+cd` blocks `crlibco` (`sin_rn` etc. and `crcoall` (`lout`) to new file `common_headers.s`
* Geant4 support in the collimation version
* Fixes for platform PPC64
* Documentation fixes

### Version 4.6.30 [09.06.2017]

* Fix sign for coupling in BB lens.
* Add new test `beam-HO_6D-simple-newstyle-coupling`, including fort.13 generator and analysis tool to extract the charge density distribution of the strong beam.
* Add new tests `beam-HO_6D-simple-oldstyle` and `beam-HO_6D-simple-newstyle`
* Update description of elens in physics manual
* Fix test javier_bignpart

### Version 4.6.29 [08.06.2017]

* Replace calls to abend() with write(lout,*) + call prror(-1); this will cause an error stop that leads to a non-zero exit status, and for CR (including BOINC), triggers the copying of the tail of fort.6 and fort.93 to STDOUT.

### Version 4.6.28 [08.06.2017]

* Fix the max amplitude to exactly amp0, do not calculate it. This makes it possible to compare top amplitude of one run with bottom amplitude of another.

### Version 4.6.27 [07.06.2017]

* Move the "BEAM-EXPERT" conversion output to file 'beam_expert.txt', write this file both for EXPERT and non-expert BEAM block.

### Version 4.6.26 [07.06.2017]

* Safe and correct repositioning of output files when doing CR
* Fix the logo in the README displaying on github.

### Version 4.6.25 [06.06.2017]

* Change STOP statements for pre-Fortran2008 compatibility

### Version 4.6.24 [05.06.2017]

* In BOINC version, call BOINCRF to convert filenames before calling open().

### Version 4.6.23 [04.06.2017]

* Always call prror(ierr) when something goes wrong
* In CR version (including BOINC), copy the last 40 lines of fort.6 and fort.93 to stderr.

### Version 4.6.22 [01.06.2017]

* In collimation version, check that number of turns is < 10k, if not then refuse to run.

### Version 4.6.21 [30.05.2017]

* Fix wrong positioning of dynksets.dat
* Check that SixTrack version and moddate matches before trying to restart
* Add dynksets.dat.canonical for test dynk_globalvars
* Some cleanup in SixTestWrapper
* Include ZLIB in the distribution; this fixes lxplus compilation problems & simplifies the scripts.

### Version 4.6.20 [26.05.2017]

* Fix beam beam EXPERT mode in 4D tracking (issue #269)

### Version 4.6.19 [09.05.2017]

* Increase nbb to 500 for big and hugenblz so that long-range beam-beam works with FCC; fixes issue #262.

### Version 4.6.18 [08.05.2017]

* Fixes bug in 6D beam beam lens (formula for the inverse of the boost, see issue #260)
* Documentation improvements
* Various small fixes to and build/test system
* Various small code cleanups

### Version 4.6.17 [20.04.2017]

* Added new FUN in DYNK "RANDON"
* Fixes and improvements to documentation

### Version 4.6.16 [18.04.2017]

* New EXPERT version of BEAM block
* Write fort.6 in SixTestWrapper under unix when no CR.
* Force static linking of libstdc++ if a static build is requested.
* Fix building old astuce with nagfor

### Version 4.6.15 [07.04.2017]

* New astuce
* Add missing input file for FCC tests with BOINC.

### Version 4.6.14 [04.04.2017]

* Allow case that wire separation is zero in one plane

    **Note:** The BOINC executables are built from commit fa708fac04c23fabf524276b14c22eb2fa7ee5d4 , not this one! The difference is that this includes Sixin.zip for the test `fcc` so that the BOINC version can be correctly tested.
    The BOINC executables were built with: `./cmake_six BUILD_TESTING BOINC API LIBARCHIVE CR BIGNBLZ`

### Version 4.6.13 [31.03.2017]

* Adding a FCC lattice test to SixTest; Requires big- or hugenblz. Also some restructuring of SixTest/CMakeLists.txt.
* Adding DYNK functions RANDU and ONOFF, modifying subroutine RANECU to also generate uniformly distributed numbers.

### Version 4.6.12 [24.03.2017]

* Remove duplicated write statements for CR

### Version 4.6.11 [21.03.2017]

* Remove output file `FirstImpacts_AcceleratorFrame.dat` from collimation version
* Fix building on aarch64
* Fix CheckTestInputForBOINC.sh
* Add Sixin.zip to tests dump4, dump5, dump6, dump_binary, 
* Fix headers of DUMP canonicals for test elensidealthin6d_DYNK_ZIPF

### Version 4.6.10 [14.03.2017]

* Add option to use the Merlin scattering routines inside SixTrack/collimation
* Fix compilation for non-x86_64 systems

### Version 4.6.9 [08.03.2017]

* Fix collimaton array bugs
* Reduce collimation console output
* Enable collimat+hugenpart. 
* Fix "tagging" of array for "other" halo - it was not correctly compacted when particles were lost, now this is fixed.
* Fix bug in DUNK allowing two dumps to ALL at the same time, with only the last having an effect
* Remove the binary dump.dat (copy of the canonical) from SixTest/dump_binary (which was masking that the binary dump was never written...), 
* Remove the ASCII version of the ALL dump from dump_binary, which was written and tested successfully, but which can't be used together with the binary dump which is what we are actually trying to test here...

### Version 4.6.8 [24.02.2017]

* Add format 3 (binary version of format 2) to DUMP

### Version 4.6.7 [22.02.2017]

* DUMP format 5, 6, and 7: Write out the beam means and matrix (angle and canonical coordinates)

### Version 4.6.6 [20.02.2017]

* Fixed windows linking problems
* Fixed test programs (checkf10 etc.) when using NAGFOR.

### Version 4.6.5 [09.02.2017]

* Added ZIPF block
* Fixed "ALL" option in DUMP
* Various fixes to CMake and CTest
* Fix potential segfault in libArchive_Fwrapper.c

### Version 4.6.4 [06.02.2017]

* Add wire element

### Version 4.6.3 [06.02.2017]

* CTEST wrapper works correctly in Windows with Checkpoint/Restart.
* In CTEST, check that Sixin.zip matches the uncompressed inputs.
* General consolidation of CTEST setup -- all active tests now pass!
* Big update of the build manual
* AVX / AVX2 build flags

### Version 4.6.2 [31.01.2017]

* In collimation version, refuse to run if no collimation block is present.

### Version 4.6.1 [31.01.2017]

* Stop if RIPP block was encountered in fort.3, instead of just skipping. RIPP was removed in version 4.5.20, replaced by DYNK. Also add the RIPP->DYNK conversion script.
* Now builds and tests on OS X, FreeBSD, NetBSD.

## Version 4.6 [27.01.2017]

* Initial version of CMAKE build system and CTEST test framework. It now builds on Linux and Windows.
* Add BOINC and LIBARCHIVE as sub-projects and link with them. Use LibArchive instead of boinc-zip to unzip Sixin.zip.
* Separate enable/disable_xp, various ascii<->double conversion routines from crlibm sources
* FMA: Work with last dump turn = -1 (i.e. until the end of the simulation)
* FMA: Detect 4D tracking and exit.
* SixTest: Delete all the testnameIN folders from SixTest.
* Various fixes for SixTest inputs and canonicals, especially for BOINC (Sixin.zip)
* Small update of user manual
* Initial add of build manual

### Version 4.5.45 [19.12.2016]

* Detect missing XSTR parameter in the BEAM block; if so print warning and set XSTR=XANG.

### Version 4.5.44 [16.12.2016]

* COLLIMATION: Modification of outputfile survival.dat: Print Nturns and Nparticles with 7 digits instead of 4.
* DYNK: Always trim zero-termination off strings before writing to STDOUT.
* Fix test SixTest check tools to not print out random memory contents, so that the basic consistency of the test tools can be checked by CTEST.

### Version 4.5.43 [12.12.2016]

* Support for BIGNPART/HUGENPART: Use 7 instead of 3 digits when writing out particle indices to STDOUT (fort.6).

### Version 4.5.42 [12.12.2016]

* If SixTrack is compiled without collimation active, STOP the program if a collimation block is present and switched on (do_coll = .TRUE.)
* Fix FMA result inconsistency between compilers by using double precision everywhere, careful conversion between complex and real numbers, fix various bugs

### Version 4.5.41 [08.12.2016]

* The COLLIMATION now block uses the emittance read from itself for placing collimators and for generating beam distributions. Also, these two emittances (distribution and collimator position) are now separated.
* DAFOR can now read the names of the input- and output-file from the command line.
* DAFOR returns a non-zero exit status in case of errors.
* DYNK can modify the reference energy
* READ90 can output to a user-selectable file name
* Bugfixes in PLATO
* Add back missing test lostnumxv

### Version 4.5.40 [22.11.2016]

* Remove comments `!hrXX <code without paranthesis>`, which is followed by the paranthesiz'ed code (which may also contain bugfixes...)

### Version 4.5.39 [22.11.2016]

* Change output format for tracks2.dat to work with FCC -- two more characters in the "s" column.
* DYNK for electron lens
* Fix "Warning: Rank mismatch in argument ‘x’ at (1) (rank-1 and scalar)" by calling daall with a size-one array instead of a scalar.

### Version 4.5.38 [08.09.2016]

* Adding compile flags "datablocks" and "hugenpart". Using these together allows tracking up to 65'536 particles in one go (thin tracking).
* Also adding flag bignpart for up to 2048 particles/run.
* Adding momentum dependence on elens kick.

### Version 4.5.37 [30.08.2016]

* Single Track File functionality.

### Version 4.5.36 [24.07.2016]

* Added new element ELEN, i.e. electron lens. This is a general infrastructure for different e-lenses (= different current distributions). Currently only the ideal hollow e-lens is implemented. The structure of the input has been chosen in order to later easily extend it also to other e-beam current distributions.
* Modifications for collimation halo studies:
  - Increased the range of the efficiency histogram vs. radial amplitude (file "efficiency.dat"), now ranging from 5-20.5 sigma.
  - Fixed radial amplitude definition used for halo studies to include also the alpha and angle.
  - Added to output files: "efficiency_dpop.dat": efficiency vs. Delta p/p to be used for off-momentum halo studies; "efficiency_2d.dat": 2D efficiency vs. radial amplitude and Delta p/p.
* Five new materials added to the collimator database: Molybdenum-graphite (MoGR), Copper-diamond (CuCD), Molybdenum (Mo), Glidcop (Glid) and Inermet-180 (Iner). Variables and indexes linked to collimator materials throughout the code also reviewed.
* Makefile now contains a proper option parser / dependency resolver. Some new checks of source consistency (syntax of astuce "if not"s) during compilation were added to both makefile and make_six . Makefile can now build 64-bit version, which works on Linux and OS X.
* Modifications of output files for collimation halo studies.
* Bugfix in DYNK for thick tracking and skew elements.
* Various fixes to SixTest, support for input file list `extra_inputs.txt` + new tests `thick6dblocks`, `thick6dsingles`, `thick6ddynk`
* Increased number of element per block nelb=280 for hugenblz 

### Version 4.5.35 [13.06.2018]

* Add FMA (Frequency Map Analysis) as a built-in postprocessing using the PLATO library for tune calculation.
* Add "check_extras" functionality to SixTest, for automatically checking any type of output file against a canonical.
* Refactorize subroutine dynk_stringzerotrim() -> stringzerotrim(). Move this and getfields_split into new deck "stringhandling".
* Fix various small bugs, misspellings, and missing comments throughout the code.
* Teach SixTest/check_10 to handle empty fort.10 files, as happens for the "exact" test.
* Small bugfixes in DYNK
* Fixes for naglib/nagdummy and DA
* Fix makefile; it now passes all tests:
* Rename makefiles used by make_six to avoid name crashes on case-insensitive operating systems (Windows and most OsX)
* Make_six simplifications: Delete twice checking for the same flag conflict, remove extraenous hardcoded -L paths for unrecognized gfortran versions, don't take externally defined shell variables for CERNLIB etc.
* Cleanup build process by removing various unused decks and flags from astuce and sixtrack.s

### Version 4.5.34 [29.03.2016]

* Single diffraction differential cross section typo fix
* Fixed string matching for TCTPH.
* Fix from Eric which makes the Differential Algebra version compile without complaining about linker error to COMNUL and various DYNK routines.
* In LINOPT: Treat crab cavity dipoles the same way as other RF multipoles.
* Improved error detection and bugfixes in COLLIMAT:
  - Detect if no particles will be tracked (such as if there is no collimation block)
  - Detect if initial distribution file is missing
  - Use an unused FORTRAN unit number when loading initial distribution
  - Detect if collimation database is missing and output appropriate error message

### Version 4.5.33 [05.02.2016]

* Recognize TCLX and TCLD collimators and make the sigma possible to set from fort.3
* Detect un-recognized collimator names, print a warning, and set nsig=1000.

### Version 4.5.32 [22.01.2016]

* Cleaned up nzfz
* Error checking: When preparing random numbers, confirm that nzfz is actually big enough.
* Input checking: Single element names (bez entries) should be unique. At least for some FCC lattices this was not the case; this will now be detected and produce an error.
* Add the unit [1] in the header for the DUMP dE/E output
* Comments and code formatting.

### Version 4.5.31 [19.01.2016]

**New features**
* Bignblz+collimat now works, fixes issue #31 .
* New make_six option "hugenblz" for nblz=400'000, requested in issue #31 .
* Force thin6d when collimation is on, such that FCC collimation studies can be done without RF, fixes #30 .
* Collimation input checking: Write and error and stop if more samples than what is possible is requested, fixes #6 .
* Clarifications in the manual
* Collimation: Detect and STOP if an unknown collimator material is encountered in the CollDB. This fixes #29 .
* General cleanup of LINOPT code.
* Remove unused arrays from collimation code.
* Detect and issue warning in case the user tries to built with crlibm + colimat; See issue #33 
* Add DYNK sanity checks during initialization.

**Bug fixes**
* Collimation: Invalid index when initializing the rvv particle array, causing potentially invalid tracking results. After fixing this, it is now possible to track up to 20'000 particles per collimation run. It also means that tracking results using collimation version may change slightly, as particle array initialization is now correct.
* LINOPT: For collimation, call WRITELIN on ALL elements, including disabled RF cavities and Crab Cavities. This fixes #8 .
* WRITELIN+CORRORB: Avoid invalid access of variable ixwl and kp() for BLOCs, which are anyway neither correctors or observation points.
* DUMP+COLLIMAT: Fix bug resulting in wrong particle indices. This fixes issue #27 .
* LINOPT: Get block names form BEZB, not BEZ.



















