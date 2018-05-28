# SixTrack Changelog

### List of Changes in Development Branch

**User Changes**

* DUMP now accepts -1 as file unit input, in which case a unit will be assigned dynamically.
* Added optional QUIET flag that will stop SixTrack from reporting initial and final values of particle pairs.
* Fixed the previously broken DA version of the code. Still may have issues with naglib when compiled with gfortran.
* Added HDF5 block for alternative ouptut for DUMP and Scatter

**Code Internal Changes**

* File units can now be dynamically allocated by the file_units module.
* Added a strings tools module to do string manipulation.
* Moved DA version core code to sixda.s90 and tshift_corr.s90.

### Version 5.0 RC1 [06.04.2018] - Release Candidate

Following is only a few key points. Full changelog will be included in final release.

**Code**

* SixTrack is now Fortran 2008 Free Form, and depcrecated syntax has been converted or removed.
* Significant portions of the code has been split out into Fortran modules, including many common blocks.
* SixTrack can now be compiled in single, double or quad precision.
* Partially completed:
  * Autoscaling of arrays. Will replace BIGNPART, HUGENPART, etc. flags eventually.
  * Dynamic file unit assignment.
* Fixed:
  * Closed orbit search.
* Removed:
  * thin/thck6dua. Replacd by DYNK. See documenttaion.
* Currently broken:
  * Electron lens.

**Documentation**

* Manual has been completely reformatted and restructured.
* Documentation still missing or incomplete for several modules.
### Version 4.7.18 [21.12.2017]

* Fixing the RF multipole elements
  * RF multipoles now affect the energy directly, not the momentum
  * Fix missing RV update in RF multipole element with FOX (DA package); this improves the symplecticity check
  * Fix errorenous factor 10^-3 in FOX part for RF mulitpole
  * Fix broken 2nd order RF multipole in FOX part (missing `+ca alignf` call)
  * Take care of px vs. xp. coordinates
  * NOTE: Compilation with TILT flag OFF is broken.
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
  * Adds offset to Gaussian beam profile.
  * Some minor changes to floating point values in calculations.

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

    **Note:** The BOINC executables are built from commit fa708fac04c23fabf524276b14c22eb2fa7ee5d4, not this one! The difference is that this includes Sixin.zip for the test `fcc` so that the BOINC version can be correctly tested.
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
  * Increased the range of the efficiency histogram vs. radial amplitude (file "efficiency.dat"), now ranging from 5-20.5 sigma.
  * Fixed radial amplitude definition used for halo studies to include also the alpha and angle.
  * Added to output files: "efficiency_dpop.dat": efficiency vs. Delta p/p to be used for off-momentum halo studies; "efficiency_2d.dat": 2D efficiency vs. radial amplitude and Delta p/p.
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
  * Detect if no particles will be tracked (such as if there is no collimation block)
  * Detect if initial distribution file is missing
  * Use an unused FORTRAN unit number when loading initial distribution
  * Detect if collimation database is missing and output appropriate error message

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

# Content from old README

### SixTrack Version: 4.5.26 

* Removed the two calls to boinc_zipitall

_McIntosh, 8th June, 2015_
 
### SixTrack Version: 4.5.25

* Give DYNK the capability of changing RF cavity settings. The new features are documented in the draft [SixTrack manual](https://github.com/kyrsjo/SixtrackTools/raw/master/user_manual/six.pdf) When implementing this, a bug affecting RF cavities which were present in the SINGLE ELEMENTS and STRUCTURE INPUT list of fort.2 but switched off (voltage = 0). These were *in some cases* still present in the lattice, slightly affecting the beam. This is the case for the LHC sequence (200 MHz RF system). This means that simulation results will from now on be *slightly* different.

_Kyrre, 3rd June, 2015_

### SixTrack Version: 4.5.24

* Reinstated the two call to boinc_zipitall

_McIntosh, 1st June, 2015_

### SixTrack Version: 4.5.23

**Eric's mods merged with Kyrre**

* Does NOT open fort.10 except for bnlelens with BOINC
* Opens and writes in postpr and leaves open for sumpos
* Subroutine abend checks fort.10 and if inexistent or empty it writes 0d0s with turn number and CPU time
* Replaced wzsubv with latest and greatest as developed for Frank
* Put close check under if boinc
* Added more CLOSE to closeUnits subroutine
* Removed tabs, ^I form source
* Fixed ampersands on calls to  h5 routines
* Commented out boinc_zipitall calls
* Added close 110/111 and added err= to close statements
* Changed iostat variable to be consistently ierro
* Added many iostat=ierro to many I/O statements
* Fixed nasty bug for nagfor, checkdist.dat, if boinc where the last character of a line in fort.10 was undefined

_McIntosh, 1st June, 2015_

### SixTrack Version: 4.5.20

* This is SVN version 201 onwards

_Kyrre, April, 2015_

### SixTrack Version: 4.5.19

* Uncommented the BUG fix? go to 360 in LIMI block
* Fixed undefine variable time in sixda
* Added SixTrack_da link for Frank

_McIntosh, 14th February, 2015_

### SixTrack Version: 4.5.18

* Collimation changes
* Line  29951: formatting in the output file 38 (tracks2.dat)
* Lines 57712 and 57715: input data for tracking (carbon density and free path length respectively)
* Lines 58583 and 58584 : correction of the calculation of momentum and energy in the subroutine CALC_ION_LOSS

_Adriana and McIntosh, 14th May, 2014_

### SixTrack Version: 4.5.17 Eric

* Removed boinc_zipitall calls
* added CLOSE 6 and 93 to routine abend

_McIntosh, 9th May, 2014_

### SixTrack Version: 4.5.16

* Going back to Version 4513 (best so far)
* Added Riccardo's physics from SVN Version 176 
* From Javier: increase bb 6d max number of slices, change logic for lhc option (strong beam optics):
  * lhc=0 symmetric optics; lhc=1 antysymmetric; lhc=2, like lhc=0 but sigma 11,33,55 overwriten

_McIntosh, Riccardo, Javier, 5th May, 2014_

### SixTrack Version: 4.5.15

* Rebuilding with absolutely NO zip stuff
* Modified the Makefile.boinc.template and removed boinc_unzip

_McIntosh, 3rd May, 2014_

### SixTrack Version: 4.5.14 Eric

* Rebuilding WITHOUT any zip calls but with the lib boinc_zip
* Just replaced the Makefile.boinc.template

_McIntosh, 3rd May, 2014_

### SixTrack Version: 4.5.14 Eric

* Rebuilding WITHOUT any zip stuff.

_McIntosh, 3rd May, 2014_

### SixTrack Version: 4.5.13 Adriana

**Modifications to Collimation Routines**
* Add output file Coll_Scatter.dat (within the Jaw subroutine) with all interactions at collimators
* Add output file FLUKA_impacts_all.dat with all absorptions at collimators
* Add Single Diffractive events in FLUKA_impacts.dat
* Add Pencil Beam type 3 (R.Bruce): The new pencil beam model uses the existing routines for annular halo generation and samples a matched halo at the face or end of any given collimator, with all particles hitting the jaws at a variable depth. There is no unphysical tilt of the jaws with the beam envelope as in the previous pencil beam
* Replacement of the scattering routine with an upgraded version (C.Tambasco and B.Salvachua). At the moment the old routine is commented out, to be seen if to delete it or put it as a compilation option. The new routine includes:
  * A different Carbon density that better fits the real collimator jaw material.
  * The addition of the logarithmic term in the in the rms angle equation of the multiple Coulomb scattering has been added.
  * Implementation of the Bethe-Bloch equation (in place of a fix value) to better estimate the ionization losses.
  * Computation of new fits of experimental data (TOTEM collaboration) to better extrapolate the proton-proton elastic and total cross sections.
  * Updating of the proton-proton single diffractive cross section considering a recent parametrization based on the renormalized pomeron flux exchange.
  * Updating of the proton-nucleus inelastic and total cross sections considering the new available data from the Particle Data Group.
* Further modification to the computation of ionisation losses taking into account the tale of the Landau distribution, since Bethe-Bloch is not valid for beta.gamma>1000 (D. Mirarchi).

_McIntosh/Adriana, 30th April, 2014_

### SixTrack Version: 4.5.12 Eric

* Just for testing new build by Xavier (not committed)

_McIntosh, 30th April, 2014_

### SixTrack Version: 4.5.11 Eric

* Put back call to unzip and zipitall
* Waiting for Riccardo's changes
* Commit for a build anyway to use on boinctest

_McIntosh, 25th April, 2014_

### SixTrack Version: 4.5.10 Eric

* NO call boinc_zipitall (missed one in 4.5.09)
* Re-instated call boinc_time_to checkpoint
* Changed CR ABEND Message

_McIntosh, 14th April, 2014_

### SixTrack Version: 4.5.09 Eric

* Based on 4508 but NO call boinc_unzip
* Based on 4508 but NO call boinc_zipitall
* Just to test ifort 2013_sp1.2.144 as well

_McIntosh, 26th March, 2014_

### SixTrack Version: 4.5.08 Eric

* Added parentheses to FOXY EXACT code to match the actual computation (and removed redundant parentheses)

_McIntosh, 7th March, 2014_

### SixTrack Version: 4.5.07 Eric

* Just trying a Windows boinc api WITHOUT call boinc_unzip

_McIntosh, 26th February, 2014_

### SixTrack Version: 4.5.06 Eric and Mattias

* Implemented "exact" drift from Mattias
* Eric changed to just use sqrt

_McIntosh 6th March, 2014_

### SixTrack Version: 4.5.05 Riccardo

* Added licence info

_Riccardo 24th February, 2014_

### SixTrack Version: 4.5.04 Eric

* Removed +if crlibm from the new vars posttime, pretime and tottime

_McIntosh 24th January, 2014_

### SixTrack Version: 4.5.03 Eric

* First attempt to fix reported CPU times. Finally used new names, removed crtime0/crtime1 and used time3 in crpoint (as well as post) for crtime2

_McIntosh 10th January, 2014_

### SixTrack Version: 4.5.02 Eric

* Now try and correct the turns reported in fort.6
* Should also try and fix CPU time reports
* Leaving debug code until I check LOST particles

_McIntosh 6th January, 2014_

### SixTrack Version: 4.5.01 Eric

* Now try and correct sumda(22/23) and Total Turn(s) Tricky, because of nnuml and multiple restarts Recomputing napxto in maincr and sumda(22/23) in postpr. Leaving debug code until I check LOST particles

_McIntosh 5th January, 2014_

### SixTrack Version: 4.4.99 Eric

* Removed the now hopefully unnecessary checkpoint of OIDPSV

_McIntosh 29th December, 2013_

### SixTrack Version: 4.4.98 Riccardo

* Fixes from Riccardo for the OIDPSV problem etc
* Based on Version 4.4.97

_McIntosh 28th December, 2013_

### SixTrack Version: 4.4.97 Eric

* Default numlcp to 1000
* Reading numlcp and numlmax as 10th/11th items in fort.3. Tracking Parameters Line 1 (was done already) but remember to set 8th/9th items niu(1 and 2) to 0 and re-create Sixin.zip
* Commented out call system in crpoint to check timing of C/R

_McIntosh 26th December, 2013_

### SixTrack Version: 4.4.96 Eric

* Added a fort.93 log message to clarify start up
* Fixed missing +ei at 611 for =if boinc

_McIntosh 23rd December, 2013_

### SixTrack Version: 4.4.95 Eric

* Put back and added checks "unzip of Sixin.zip"

_McIntosh 23rd December, 2013_

### SixTrack Version: 4.4.94 Eric

* Commented out the "unzip of Sixin.zip" for tests with Igor

_McIntosh 23rd December, 2013_

### SixTrack Version: 4.4.93 Eric

* FIXED the mywwerf problem by using a CONSTANT P
* Found that unzip Sixin.zip was over-writing fort.93 etc
* Decided to just re-position fort.93 again and empty fort.92 
* WRITEBIN not to write message if +bnlelens
* Do NOT enable_xp in dabnew.s for Lahey lf95
* Enable for other compilers if +crlibm, but I use the binary fort.111
* In sixtrack.s it is done only for +fio which isn't supported by lf95 (should carefully check when testing fio)
* Commented out some +if debug commands
* Riccardo will add fix for OIDPSV and then I will test
* See README for all changes since Version 4462

_McIntosh 19th December, 2013_

### SixTrack Version: 4.4.92 Eric

* Changing dabnew.s to handle better fort.111 etc
* Note fort.111 is different (just * format!) and fort.18 is OK.
* Problem is bnl and nagfor combination
* Don't write to binary files in WRITE6 +if bnlelens

_McIntosh 20th November, 2013_

### SixTrack Version: 4.4.91 Eric

* Merged version of sixtrack.s.v4490 with Yngve's
* collimations mods based on v4463

_McIntosh 13th November, 2013_

### SixTrack Version: 4.4.90 Eric

* Adding oidpsv(j) to C/R as a TEMPORARY??? fix for test javier.
* Still to run regression tests before testing new NUMLMAX NUMLCP
* Modified from Version 4.4.87.
* See README for full info.

_McIntosh 12th November, 2013_

### SixTrack Version: 4.4.89 Eric

* Added a dump/abend after restart before Tracking Turn 2.
* This shows that oidpsv(j) is DIFFERENT (dpsv(j) is changed in both CRABAMP elements and by multipoles).

_McIntosh 30th October, 2013_

### SixTrack Version: 4.4.88 Eric

* Adding debug for every element and a dump/abend in order to save C/R files after 1 Turn.

_McIntosh 30th October, 2013_

### SixTrack Version: 4.4.88xory Eric

* Just saved this with xory added to C/R. Doesn't help.

_McIntosh 28th October, 2013_

### SixTrack Version: 4.4.87 Eric

* Still a problem with test javier on "cr" version (+/- bnl???)
* Fixed problem on +if cr for binrec in WRITEBIN
* Changed "writing" message in WRITEBIN to match "written" (binrec+1)
* Fixed the over-write fort.10 for only if .not.boinc 
* Still to run regression tests before testing new NUMLMAX NUMLCP

_McIntosh 25th October, 2013_

### SixTrack Version: 4.4.86 Eric

* Fixed MAJOR design issue with RESTART by introducing a new START flag and optionally calling UNZIP Sixin.zip
* Running regression tests before testing new NUMLMAX NUMLCP

_McIntosh 22nd October, 2013_

### SixTrack Version: 4.4.85 Eric

* First version for Igor
* Changed fort.93 messages in WRITEBIN and CALLCRP to improve readability
* Added fort.93 messages for call POSTPR
* Changed DO until for clarity of tracking loop

_McIntosh 21st October, 2013_

### SixTrack Version: 4.4.84 Eric

* Various bug fixes in new NUMLMAX option

_McIntosh 18th October, 2013_

### SixTrack Version: 4.4.83 Eric

* Introduced an inetrface routine CALLCRP to acll CRPOINT

_McIntosh 18th October, 2013_

### SixTrack Version: 4.4.82 Eric

* Changed to use Sixin.zip instead of fort.zip
* ENDFILE and CLOSE 92 in ABEND to save time and space
* Fewer messages from writebin
* Added boinc_unzip_ to boinc_zipitall.cpp
* Split writebin and crpoint to writebin and callcrp

_McIntosh 18th October, 2013_

### SixTrack Version: 4.4.81 Eric

* Fixing a problem with numlcp=9 numlmax=100 when we never get past 200 turns!
* Still to implement Sixin.zip/Sixout.zip
* Added call to boinc_zipitall and boinc_zipitall.cpp
* Added crlibm option for fort.54 blrdis read dist

_McIntosh 14th October, 2013_

### SixTrack Version: 4.4.80 Eric

* Added new daten input numlcp and nummlmax
* All the changes for numlmax turns per job
* Use numlcp instead of defaulting to nwri (nwr(3))
* Need to verify a problem with BNL .dat output

_McIntosh 13th October, 2013_

### SixTrack Version: 4.4.71 Eric

* Removed two dangling +ei
* Fixed up bnlelens write(10|51|52|97) to use dtostr
* Merged Collimation mods from SVN Versions 149 and 150
* Added NAGLIBC for SixTrack_DA to make_six (messy as is the whole question of cernlib/gcc libs etc) Tested SixTrack_DA with naglib, gfortran and ifort.

_McIntosh 12th October, 2013_

### SixTrack Version: 4.4.70 Eric

* Just testing nagfor with -float-store???

_McIntosh 7th October, 2013_

### SixTrack Version: 4.4.69 Eric

* Removed a redundant integer int in function bran and changed a 10 to 10.0D0 both in dabnew.s
* Added a message when there is an error reading fort.23
* Fixed a nasty bug with [en|dis]able_xp for nagfor+crlibm only.

_McIntosh 7th October, 2013_

### SixTrack Version: 4.4.68 Eric

* Just testing ifort xe2013/composer_xe_2013_sp1.0.080 O0, O1 and O2 and ia32, sse2 and sse3.

_McIntosh 30th September, 2013_

### SixTrack Version: 4.4.67 Eric

* Just testing gfortran O4

_McIntosh 26th September, 2013_

### SixTrack Version: 4.4.66 Eric

* Just testing ifort 2013 with O1 (ia32, sse2 and sse3)

_McIntosh 20th September, 2013_

### SixTrack Version: 4.4.65 Eric

* Re-build with new boinc libs from server_stable (Riccardo)

_McIntosh 18th September, 2013_

### SixTrack Version: 4.4.64 Eric

* Just testing new ifort 2013 (usevnewifort)

_McIntosh 16th September, 2013_

### SixTrack Version: 4.4.63 Eric

* Just added a call to boincrf for fort.zip for Windows.

_McIntosh 31st August, 2013_

### SixTrack Version: 4.4.62 Eric

* Fixed nasty bug in daten iclr.eq.2 concerning exz in tracking input
* Fixing all "formatted read" for fio and _xp stuff.
* dabnew.s and sixtrack.s. if fio (but NOT Lahey lf95) then enable/disable _xp before/after the READ "NEAREST" fio overrules crlibm for formatted input

_McIntosh 29th September, 2013_

### SixTrack Version: 4.4.61 Eric

* Added code to daten for crlibm and read of xstr

_McIntosh 16th August, 2013_

### SixTrack Version: 4.4.60 Eric and Javier

* Included Laurent's new disable_xp, enable_xp for WINDOWS
* Removed "MY" crlibm/features.h (but kept fpu_control.h)
* Used Javier's latest features

_McIntosh 16th August, 2013_

### SixTrack Version: 4.4.56 eric and Igor and Laurent

* Same as 4.4.52 but pgf90 has -Mnoflushz (for SSE)

_McIntosh 28th July, 2013_

### SixTrack Version: 4.4.55 eric and Igor and Laurent

* Same as 4.4.52 but pgf90 traps underflow

_McIntosh 26th July, 2013_

### SixTrack Version: 4.4.54 eric and Igor and Laurent

* Same as 4.4.52 but ifort has -no-ftz for Linux

_McIntosh 21st July, 2013_

### SixTrack Version: 4.4.53 eric and Igor and Laurent

* Same as 4.4.52 but ifort has -ftz for Linux

_McIntosh 20th July, 2013_

### SixTrack Version: 4.4.52 eric and Igor and Laurent

* Just fixed the type definitions in myboinc.f which is for internal use only.

_McIntosh 12th August, 2013_

### SixTrack Version: 4.4.52 eric and Igor and Laurent

* Deleted a redundant +ca cseeds (came in v4443!)???

_McIntosh 4th July, 2013_

### SixTrack Version: 4.4.51 eric

* Really fixes checkpoint and adds log message to fort.93 for the last checkpoint
* Adds parameters 8,1 to boinc_zip call (string lengths)

_McIntosh 1st July, 2013_

### SixTrack Version: 4.4.50 eric

* Always checkpoint (do NOT call BOINC time to checkpoint)
* Do a checkpoint when all turns completed before post-processing (unless restarted from that checkpoint) in order to facilitate ten million turns etc
* Replaced g77 option by pgf90 (preliminary)

_McIntosh 27th June, 2013_

### SixTrack Version: 4.4.49 eric 

* Only disable_xp/enable_xp for Nagfor ifdef crlibm
* Do not use errno in round_near.c (for Cygwin)

_McIntosh 30th May, 2013_

### SixTrack Version: 4.4.48 eric 

* Fixed crlibm/round_near.c and function fround (nasty BUG)
* Used maxf = 30 in call to round_near

_McIntosh 29th May, 2013_

### SixTrack Version: 4.4.47 eric 

* Fixed read(23 for crlibm in routine readd1

_McIntosh November 2012_

### SixTrack Version: 4.4.46 eric and laurent

* Now call boinc_zip_ directly and use boinc_fraction_done instead of boinc_sixtrack_progress 
* Added boinc_fraction_done etc to myboinc.f

_McIntosh 26th September, 2012_

### SixTrack Version: 4.4.45 eric

* Just made sure MMAC is zero (when NMS=0) for lost particles

_McIntosh 26th September, 2012_

### SixTrack Version: 4.4.44 eric and laurent

* Removed the integer parameter from call to boinc_init

_McIntosh 26th September, 2012_

### SixTrack Version: 4.4.43 eric and javier

* Javier's mods for new elements called allMULT
* Fix for sumpos for mmac and nms
* make_six new SSE defaults and checks and removed set -x

_McIntosh 22nd September, 2012_

### SixTrack Version: 4.4.42 eric

* Fixed SUMDA 59/60 for all particles lost early

_McIntosh 22nd July, 2012_

### SixTrack Version: 4.4.41 eric

* Quickly changed to SSE3 (but still called SSE2) for Igor for BOINC
* Added total turns and CPU time to fort.10 words 59/60
* uname -m instead of fs sysname
* sse3 instead of sse2  and nothing for Macs

_McIntosh 3rd July, 2012_

### SixTrack Version: 4.4.40 eric

* Same SixTrack but new make_six for MacOS
* Defaults to ifort instead of lf95
* uname -m instead of fs sysname
* sse3 instead of sse2 for Macs
* enable_xp.c/disable_xp.c dummied on ifdef APPLE
* Comma separated map options for FCLMAP
* NO static linking for Mac
* Fixed bug with make_six options and added lf95

_McIntosh 22nd June, 2012_

### SixTrack Version: 4.4.40 eric

* Put write 210 (fort.10 HEX) under +if debug
* Added the crlibm .f routines (acos_rn.f asin_rn.f atan2_rn.f) to sixtrack.s
* make_six (for MacOS), boinc_api_fortran[.windows] and new make_windows with new templates and flags for +/-SSE2 

_McIntosh 18th June, 2012_

### SixTrack Version: 4.4.37 eric

* based on 4.4.35
* Changes for writing fort.10 using dtoa_c.c, 
* function dtostr_ calls function dtoaf_ calls dtoa. 
* Removed redundant error messages in dtostr

_McIntosh 29th May, 2012_

### SixTrack Version: 4.4.36 eric 

* Never committed (just to test g77) which does
* NOT work as it can't even open Unit 90!
* So I just deleted all the references to binary writes and built NO bnlelens NO C/R etc

_McIntosh 29th May, 2012_

### SixTrack Version: 4.4.35 eric and Harry

* based on 4.4.34
* Javier's mods for crabamp with Harry's fixes labelled hr13
* Removed  a double assignment to ferrno in round_near.c

_McIntosh 18th May, 2012_

### SixTrack Version: 4.4.34 eric

* Fixed assignment of data to sigma0 and qw0 in daten.
* Fixed a few problems with g77 and concatenation in new daten subs splitfld, spliterr, rounderr.
* Unfixed initialisation of DA variables in dafor.f.
* Not committed

_McIntosh 28th April, 2012_

### SixTrack Version: 4.4.33 eric

* based on 4.4.11 
* Uses strtod and sprintf and fround for all input variables and for formatted write of fort.10. This is the first 0 ULP version giving identical results on Linux with 4 different Fortran compilers. 

_McIntosh 11th April, 2012_

### SixTrack Version: 4.4.11 eric

* based on 4.4.01 (never committed)
* Uses round_ulp from Laurent roundnulp( ,1024) only for the same variables, extaux(k), alignx|z and tilt

_McIntosh 25th March, 2012_

### SixTrack Version: 4.4.01 eric

* Based on 4.4.00 (never committed) which merges all my latest ifort/portability mods including DBLE(SNGL)), exponentiation, DP lfitw, errf, postprocessing etc. with SVN Revision 119 to trunk/SixTrack
* Now uses new Fortran Function myround13??? to (hopefully) solve problems with different results on ifort formatted input between Linux and Windows (and Nagfor lf95).
* No longer calls "worker" (no graphics) at L. Deniau request
* Call disable_xp for nagfor only (will also be needed for gfortran)

_McIntosh 8th February, 2012_

### SixTrack Version: 4.2.11 eric

* Mainly fixing unitialised variables astuce: Removed tabs and initialised ILEV
* dafor:
  * Removed BLOCKDATA and used subroutine predata instead
  * Initialised many variables and moved the ILIN = (ILAST(AC)...
  * Generated code to initialise variables allocated by DAALL
* sixtrack:
  * Changed comnul for AAIV, BBIV, and TASAU
* New make_six with nagfor and FCF and added  `-L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/32 to NAGLIBC`
* New crlibm Makefile with -fno-strict-aliasing for the type-punned problem
* Fixed lfit/lfitw argument dimensions
* changed log10(real(xxaux)) to real(log10(xxaux))
* This version was checked out on SLC5 as well and with NAG nagfor C=all and -C=undefined and -float-store and Lahey for production Changed the NAGLIBC definition for SLC4/5 compatibility

_McIntosh 8th August, 2010_

### SixTrack Version: 4.2.10 FRS/eric

* Changed dimension of phase(3,npos+1) in common phasecomm Corrected an overwrite in postpr in rare cases

_FRS/McIntosh 1st July, 2010_

### SixTrack Version: 4.2.9 FS

* Proper stop for single elements with names longer then 16 characters which may happen for SixTrack input file fort.2 produced by MAD-X (could be 24). Before, the 17th column was checked for non-blank character which is a special case only. It is checked for 16 non-blank characters, leading blanks allowed. Moreover, the code is stopped in case the length of 80 characters is exceeded without minimal input, i.e. a mininum of 8 characters is needed for valid input.

### SixTrack Version: 4.2.8 eric

* Commented out the boinc graphic init, finish and progress calls
* Commented out the routines in myboinc.f

_McIntosh 7th March, 2010_

### SixTrack Version: 4.2.7 eric

* Just fixed NOT to close 51, 52,53 if bnlelens and BOINC.

_McIntosh February, 2010_

### SixTrack Version: 4.2.6 Revision  1.35 eric and frs

* Added Closed Orbit to sumda(53:58) for fort.10.
* Changed a couple of constants in "Wire" to use c1m7.

_Frank Schmidt/McIntosh 15th January, 2010_

### SixTrack Version: 4.2.5 Revision 1.33/1.34 frs

* Number of BB encounters up to 350

_Frank Schmidt 5th/6th October, 2009_

### SixTrack Version: 4.2.4 Revision 1.31/1.32 adriana

**Changes in the collimation part:**

* Added and updated alignment errors (tilt, offset, gap-size): all values of alignment errors can be kept identical for different runs, by applying the same seed. Therefore the random function myran_gauss and rndm5 have been added and the inputs fort.3 adapted.
* Alignment errors are also applied for deformed jaws.
* Added to the pencil beam section the possibility to generate different particle distributions (Gaussian and rectangular in x and y) on the selected collimator.
* Changes to the do_select option to get the multi-turn halo information for all particle packets.

_C. Bracco / A. Rossi / Th. Weiler_

### SixTrack Version: 4.2.3 Revision 1.30 frs

* Fix the missing phase advance in the solenoid. This thin element is special in the sense that the it creates a direct phase advance that has to be calculated at every element occurrence. Fix the missing phase advance in the solenoid. This thin element is special in the sense that the it creates a direct phase advance that has to be calculated at every element occurrence.

_Frank Schmidt 27th July, 2009_

### SixTrack Version: 4.2.2 Revision 1.29 yipeng

* First attempt at thin solenoid (work done by Yipeng Sun) The tune calculation in the traditional SixTrack part seems slightly off. DA part seems okay and in agreement with MAD-X.

_Yipeng Sun & Frank Schmidt 24th July, 2009_

### SixTrack Version: 4.2.1 Revision 1.28 frs

* Fixed computed go to dipedge element
* Got rid of some warnings in crlibm/csh_fast.c
* New make_six for 32/64-bit MACHTYPE with SLC4/SLC5

_Frank Schmidt 6th July, 2009_

### SixTrack Version: 4.2.0 Revision 1.27 frs

* Adding dipedge element which allows to track all relevent elements except the solenoids

_Frank Schmidt 14th May, 2009_

### SixTrack Version: 4.1.16 CVS Version 1.26 McIntosh

* Deleted the debug FCL for g77
* Fixed README and update make_six for SLC5 (cernlib and -static)
* Added lib32 with the cernlib .a libraries and X11

_McIntosh 20th March, 2009_

### SixTrack Version: 4.1.16 CVS Version 1.26 McIntosh

* Small fix to phase trombone from Guillaume and Yun

_McIntosh 11th November, 2008_

### SixTrack Version: 4.1.15 CVS Version 1.25 McIntosh

* Write CRPOINT messages 5 times maximum
* Set n_cut and n_nocut to 0 after printing

_McIntosh 21st October, 2008_

### SixTrack Version: 4.1.14 CVS Version 1.24 McIntosh

* Now endfile SixTwiss, checkdist, and beambeam-output and fort.10 fort.97 fort.51. 
* Debug code for SIGSEV on BNL in crstart.
* Implemented additional debug dump routines (needs testing)
* Remove n_cut=0 and n_nocut=0 (wrongly added by me)
* Re-instated the second C/R file fort.96 and added code to check extended checkpoint in crcheck. 
* Two versions of writelin (extra 7th argument for collimat and bnlelens)
* Moved REAL time2 to COMMON ttime and read95/read96 to COMMON as well.

_McIntosh 3rd October, 2008_

### SixTrack Version: 4.1.13 CVS Version 1.23 McIntosh

* If bnlelens AND lhc.eq.9 AND NOT boinc write one line to fort.10 containing the sixtit(1:60)
* if boinc AND NOT restart write and add one to bnlrecs Cleaned up all `read/write[](` to `read/write(` Fixed all calls to rndm4 in +collimat to be to `dble(rndm4())`

_McIntosh 8th September, 2008_

### SixTrack Version: 4.1.12 CVS Version 1.22 McIntosh

* Fixed the SixTwiss output to write one line F20.13 and most importantly changed the bnlrec count for C/R

_McIntosh 2nd September, 2008_

### SixTrack Version: 4.1.11 CVS Version 1.21 McIntosh

* make_six copies only necessary .ast and .f files Sorted .f files for crlibm for windows and makes a link from the new executable to SixTrack
* SixTrack for BNL only, use fort.54 for beambeamdist.dat for BOINC and CPSS. Use fort.52, fort.53, fort.51 and fort.97 for beambeam-output.dat, beambeam-lostID.dat, SixTwiss.dat, and checkdist.dat for CPSS but ONLY fort.10 for BOINC.

_McIntosh 31st August, 2008_

### SixTrack Version: 4.1.10 CVS Version 1.20 McIntosh

* Added SixTwiss output to all tracking routines if bnldata
* Use napx to write checkdist if bnlelens
* Do NOT write a second C/R file fort.96

_McIntosh 23rd August, 2008_

### SixTrack Version: 4.1.9 CVS Version 1.19 McIntosh

* Fixed a problem with binary output files when using the bnlelens option for normal DA runs.
* Removed redundant plotting initialisation.
* Forced -cernlib for +windows in make_six.
* Added !GRDRHIC/!GRD-042008 comments.

_McIntosh 23rd August, 2008_

### SixTrack Version: 4.1.8 CVS Version 1.18 McIntosh

* Fixed a problem with the collimat option (my fault) Need a regression test for this.
* Added C/R for synuthck, 6d thick lens and variable IL. Modified crpoint,crcheck, and crstart and the size of the C/R file.
* If bnlelens, use napx rather napx00/npart for number of pairs
* Most important is that bnlens option is now an addition to the normal tracking so that the SAME executable can be used for normal DA and bnlelens LHC=9 runs (important for BOINC/CPSS) but LHC=9 does NOT write binary files nor post-process.

_McIntosh 20th August, 2008_

### SixTrack Version: 4.1.7 CVS Version 1.17 McIntosh

* make_six has a new 'debug' option to aid SixTrack development.
* If debug is selected, a new flag to ASTUTE SixTrack makes  available Unit 99 for messages and dumps and a set of dump routines (of which only a full dump in this version).
* A couple of write(* in SUBRE are now handled by C/R.
* A couple of bugs with C/R and IDFOR have been fixed.
* All subsequent changes noted here are for the 'bnlelens' and other SixTrack functionality should be unchanged (except that fort.95 and 96 are used for Checkpoint/Restart, C/R, instead of 12 and 13.)
* The bnlelens option has been implemented in comdecks bnlin and bnlout principally so it works for all six tracking routines. This has been tested only for THCK6D so far. There are now no changes to trauthck and trauthin.
* The C/R option for 6d and thick lens is disabled unless debug is selected.
* C/R now handles additional variables and correctly positions the beambeam-output.dat and beambeam-lostID.dat files.
* The binary files 90-59 are not used; no post-processing is performed. All OPEN/CLOSE have been moved to the standard comdecks and the variables moved to COMMON initialised by COMNUL.
* The READDIS routine for bnlelens is renamed to BNLRDIS to avoid confusion with collimation and it uses UNIT 54 to  read beambeamdist.dat which may now have only one sample.
* The n_cut and n_nocut variables are set to zero before the j loop.

_McIntosh 18th August, 2008_

### SixTrack Version: 4.1.6 CVS Version 1.16 McIntosh

* and of course I forgot the SixTrack Version
* Version and moddate set to 4.1.6, I am skipping 4.1.4 and 4.1.5 and then CVS and SixTrack versions should correspond. Otherwise we get confused.

_McIntosh 5th August, 2008_

### SixTrack Version: 4.1.3 CVS Version 1.15 McIntosh

* Second interim update for BNL.
* synuthck replaced computed goto by IF's due to particularly nasty Lahey lf95 bug.
* Patched make_six mkwindows to sort out logsix.c and .h until we can sort out CVS
* Under +if bnlelens set ch1 to "" in subroutine intepr
* Under +if bnlelens replaced  rvv(j)= by rvv(i)= in trauthck
* Under +if bnlelens added k=0 before 1st linopt call to writelin

_McIntosh 5th August, 2008_

### SixTrack Version: 4.1.3 CVS Version 1.14 McIntosh

* First interim update to facilitate development with BNL.
* make_six creates a directory and executable name based on options.
* make_six now supports [+-]bignblz which sets nblz to 200,000!!! if selected.
* The Makefile creates a map and uses -g in FC for NAG.
* crlibm is cleaned up for Linux/Windows and logsix.dat and logsix.h are used instead of log.c and log.h and inlining is handled with an IFDEF.
* In bnlelens Unit 97 (not 98) is used for checkdist.dat.
* A problem with open(10... and open(99 is fixed for NAG.
* sigsecut2 comment is cleaned up.

_McIntosh 30th July, 2008_

### SixTrack Version: 4.1.2 CVS revision 1.13 rtomas

* Fixed "2d0*pi" bug in crab cavity kick

### SixTrack Version: 4.1.1 CVS revision 1.12 McIntosh

* New make_six with new flags bpm and bnlelens and the ast's
* sixtrack.s: moved GRDRHIC comments inside +if bnlelens
* Changed some writelin statements to be consistent; under investigation.
* Fixed FORMAT bug in bnlelens.

_McIntosh 11 July, 2008_

### SixTrack Version: 4.1.0 CVS revision 1.11 McIntosh

* The first commited version with Guillaume's preliminary changes for "bnlelens" using the deck "rhicelens" required for RHIC BEAM-BEAM studies.
* make_six now determines the version and modification date from sixtrack.s, and displays them
* Sets the BOINC variable correctly for the Makefile.
* The program maincr to call worker is added to myboinc.f.

_McIntosh 11 July, 2008_

### SixTrack Version: ?  CVS revision 1.10 frs

* bpmdata are only written to files units > 100 if bpm flag in track.ast.

_FRS 10 July, 2008_

### SixTrack Version: 4.0.10 CVS revision 1.9 McIntosh

* The version and last modification date for both SixTrack and SixTrack_da are now specified in a common deck at the very beginning of sixtrack.s. Version is an 8 character string [v]v.[v]v.[v]v and moddate a 10 character string dd.mm.yyyy. The version is stored in floating-point in sumda(52) for BOINC in particular. In this way fort.10 identifies the version which produced it.
* Small format changes to a couple of comments and the SIXDA date/time corrected.

_McIntosh 24 June, 2008_

### SixTrack Version: 4.0.9    CVS revision 1.8 frs

* sumda(52) set to 4009.0

_FRS 19 June, 2008_

### SixTrack Version: ? CVS Revision 1.7 frs

* replace CRAB /c1e3 by *c1m3

_FRS 15 January, 2007_

### SixTrack Version: ? CVS Revision 1.6

* ???  CRAB Cavity?

### SixTrack Version: ? CVS Revision 1.5

* ???  CRAB Cavity?

### SixTrack Version: ? CVS Revision 1.4 frs

* Increase number of elements nblz=20000 (actually 200,000)

_FRS 05 July, 2007_

### Sixtrack Version: ? CVS Revision 1.3 frs

* Program stops when single element names longer than 16 characters.
* Program stops when when order of multipoles in fort.3 exceed MMUL.
* MMUL is set to 20 for all non-collimation versions.

_FS 09 May 2007_

### Sixtrack Version: ? CVS Revision 1.2 robertde

**Latest version of the expanded version for collimation studies (collimat flag)**

* change of the collimator labels step to take into account the new Beam 2 lattice from MADX
* "STOP" flags have been checked and removed if necessary
* fort.XX units checked: no unit number above 99 shall be used

_GRD-SR 26/09/2006_

### SixTrack Version: 4.0.08 

**Latest version of the expanded version for collimation studies (collimat flag)**

* creation of 3 additive "makedis" subroutines, incl. one that allows to read any beam distribution from a given input file
* correction of "bug" from previous RHIC version that added a fake aperture limitation of 4 cm during tracking (particles over that amplitude were considered as absorbed)
* new treatment of collimator material: one can now simulate the flatness of any/both jaws with fit parameters to be specified in the fort.3 file (COLLIMATION block)

_GRD-SR 20/10/2005_

### SixTrack Version: 4.0.07 

**Version 1.3 make_six 05/09/2005 and Makefile updated**

* CERNLIB only for SixTrack (not SixTrack_da); CERNGRAF dropped.
* The -naglib option now works for SixTrack_da
* Fixed problem with missing +ei in crpoint for CPSS version 

### SixTrack Version: 4.0.06

**Version 1.2 make_six 29/08/2005**

* CERNLIB//GRAFLIB re-defined
* Forced lf95 for crlibm
* Added g95 option
* Redefined "our" isnan as myisnan and dropped it from crlibm (Should actually make a crlibm.s for crlibm with a Windows option)

### SixTrack Version: 4.0.05

**Latest version of the expanded version for collimation studies (collimat flag)**

* creation of impact, absorption and FLUKA files
* correction/upgrade of writing process for halo files: trajectory written at entrance and exit of each collimator
* new treatment of RHIC lattice with possibility of playing with the 2 different apertures of its primary collimator

_GRD 14/6/2005_

### SixTrack Version: 4.0.04 boinc bis

* Fixed a corrupted makefile and added Makefile.boinc

_Eric 2/4/2005_

### SixTrack Version: 4.0.04 boinc

* Added boinc option to sixtrack.s dabnew.s and to make_six using myboinc.f

_Eric 31/3/2005_

### SixTrack Version: 4.0.04 bis

* Didn't touch SixTrack at all.
* Fixed comment and bug in make_six re naglib default.
* Created directory windows for crlibm for Windows.
* Added some options for cpss/windows to generate the src code and crlibm using the script mkwindows and added checks for options.

_Eric. 22/3/2005._

### SixTrack Version: 4.0.04

* Fix the problem with common blocks in sixda.f.
* Fix a small bug in dafor.f 
* Minor changes to make_six and Makefile.

### SixTrack Version: 4.0.03

* Completely new make_six to make SixTrack,(default), SixTrack_da (da) with crlibm/crlibm.a,Cernlib and Naglib.
* Options are: tilt tracking fast crlibm windows cernlib naglib da collimat cpss boinc cr nag g77

**Defaults are:**
* tilt tracking fast crlibm cernlib naglib
* -option to delete it e.g -cernlib
* [+]option to activate e.g +cr for Checkpoint/Restart

Uses the ast files in ast_mask e.g. make_six (defaults to tracking tilt fast crlibm cernlib) make_six collimat (does the same but with collimation) Supported compilers are lf95 (default) or NAG f95 or g77. Fixed a bug with crbinrecs in writebin in sixtrack.s

### SixTrack Version: 4.0.0

**Comments:** Includes beam collimation, checkpoint/restart, and crlib options along with major code cleanup.

**Options:** The following sub-directories, ast_crcrlib ast_crlib ast_straight ast_tilt,  contain the .ast files for building as follows: The straight/tilt builds without/with tilt. The crcrlib builds tilt with checkpoint/restart and crlib. The crlib builds tilt with crlib.

**Default:** The .ast files in the current directory are for ast_crlib on SLC3 Linux at CERN. The crlibm routines in CRLIBOBJS.a are for Linux and have been compiled with lf95 6.2 and gcc 3.2.3.

The following options are available (after copying the desired ast files to the current directory if non-default):

* makex clean
* makex sixtrack_vector_new
* makex sixtrack_da_new
