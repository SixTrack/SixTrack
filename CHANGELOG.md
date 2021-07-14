# SixTrack Changelog

### merged in master

**Bugfixes**

* Fix the root build. PR # 1040 (J. Molson)
* More robust detection of lxplus at compilation. PR #1045 (J. Molson).
* Fix pencil beam type 3 - the optics function at the entrance of the collimator were always used for beam sampling, even when those at the exit should have been used (e.g. because the beam is divergent on the cleaning plane). PR #1046 (A. Mereghetti).
* Do not update the pair mapping for non-primary particles. PR #1050 (A. Mereghetti)
* Increased number of digits for particle ID in FirstImpacts.dat and in collimator length in coll_summary.dat (to properly display crystal collimators which are usually a few mm long). First impacts on crystal collimators are now correctly flagged and a missing check on the `dowrite_impact` flag when writing Coll_Scatter.dat has been added. PR #1053 (M. D'Andrea)
* When collimator settings are required to match those read from an old format CollDB, a separate subroutine reconstructs the family settings based on the most frequent setting in each family. PR #1053 (M. D'Andrea)
* Do not perform the pair mapping when geant4 collimation is enabled. PR #1059 (J. Molson)
* Enable single sided collimators with geant4 collimation. PR #1060 (J. Molson)
* Fix building with geant4 collimation with geant4 releases >= 10.06. PR #1060 (J. Molson)
* Fix a mass miss-match with geant4 when entering non-ground state ions into geant4. PR #1062 (J. Molson)
* Fix a crash with miss-matched format strings when writing the aperture losses file with geant4 enabled (and not FLUKA). PR #1062 (J. Molson)
* Fix building with gcc >= 10. PR #1076 (J. Molson)
* Fix 2 regressions in the K2 collimation cross section calculations from version 4 to 5. PR #1077 (J. Molson)
* Fix the time coordinate when using geant4 based collimation. PR #1078 (J. Molson)
* Fix an incorrect file header in mod_dist. (J. Molson)
* Fix a problem with the circular aperture check not working as expected (J. Molson)

**User Side Changes**

* When specifying `XP` and `YP` in the `FORMAT` statement of the `DIST` block, the units are parsed. Accepted values are [1], [1000], [MRAD], [RAD]. PR # 1054 (A. Mereghetti)
* Electron lenses have been inserted into FOX - PR #839 and #1056 (A. Mereghetti).
* Increased flexibility of e-lens module - PR #841 and 1056 (A. Mereghetti):
  * elens module fully dynamic allocatable;
  * give possibility to express R_1 and R_2 in sigma;
  * add any ion species to be defined as possible lens beam;
  * degenerate WIRE type of e-beam distribution is correctly handled;
  * other changes, including:
    * relativistic gamma of lens beam added to calculation of theta_R2;
    * removed remaining signs of chebyshev polynomials in elens module;
    * empty lines allowed in file describing the radial profile;
    * fixed bug in geometric normalisation factor of GAUSSIAN and RADIAL prpofiles;
  Documentation changed accordingly (user and physics manual).
* When sending particles to geant4, if the particle mass is within a tolerance of the geant4 value, update the mass to this value and re-scale the particle energy. PR # 1055 (J. Molson).
* If no collimator are found for a given family, the aperture of that family is set to zero. PR #1053 (M. D'Andrea)
* If a particle interacts with a crystal collimator after having previously interacted with another or the same crystal collimator, the process ID of the previous interaction is stored in cry_interaction.dat. PR #1058 (M. D'Andrea)
* Collimator material names are now case insensitive in geant4. PR #1062 (J. Molson)
* In the HION block the PDG ID can now be set as the 5th value. PR #1062 (J. Molson)
* Foxified lenses with Chebyshev maps - PR 849 (A. Mereghetti). In addition:
  * Chebyshev maps are actually used on a squared domain;
  * the check against R1 and R2 takes into account rounding issues;
  * max order of cheby polynomials is a parameter in the module header;
  * derivatives of Chebyshev polynomials, used to compute the kick, no longer go through a dangerous division by 1-u^2;
  * update chebyshev tests;
  * tabular method for inputting/outputting coefficients;
  * some house-keeping.
* Only disable generating fma instructions with gcc if CRLIBM is enabled.

**Code Improvements and Changes**

* Removed updating napxo variable in the context of the Fluka-SixTrack coupling. This allows not to screw-up pair mapping in the context of DA studies. PR # 1052 (A. Mereghetti)
* Removed the un-used fluka_init_brhono function. PR # 1055 (J. Molson).
* Print error codes from the fluka coupling. PR # 1055 (J. Molson)
* Update FLUKAIO reference. PR #1057 (J. Molson)
* Add Si as a possible collimator material for G4.  PR #1059 (J. Molson)
* Add particle ID and parent ID tracking with geant4. PR #1062 (J. Molson)
* Allow setting particle statistical weights in geant4 PR #1062 (J. Molson)
* Always enable the EMD physics process in geant4. PR #1062 (J. Molson)
* Use global id/parent/weight variables in the FLUKA coupling. PR #1062 (J. Molson)
* Start to enable the ability to use collimation with thick lens lattices. PR #1062 (J. Molson)
* Updated some physical constants to use their now fixed values. PR #1077 (J. Molson)
* Always output the collimator energy loss file fort.208 in all collimation modes. (J. Molson)
* Add the ability to use the dump module with root output. (J. Molson)

### Version 5.4.3 [19.12.2019] - Release

**Bugfixes**

* Fixed parsing of the `DIFF` block where only the first element name was parsed. The remaining elements were ignored. PR #1031 (F. Schmidt, V.K. Berglyd Olsen, J. Molson)
* Added a missing else statement in the DA code that was unintentionally removed when the exact drift code was added in 2014. PR #1028 (J. Molson, V.K. Berglyd Olsen)
* Fixed the conversion of zeta to sigma in `DIST` and the user manual. It was previously assuming the inverse definition of `rvv` used by SixTrackLib. PR #1027 (R. De Maria, V.K. Berglyd Olsen)
* Fixed a type definition inconsistency introduced by #878 that affected the ROOT interface. PR #1026 (J. Molson)
* Fixed the `CERNLIB` build system to support 64-bit. PR #1025 (F. Schmidt, J. Molson)
* Fixed a faulty loop in the FFT routine in post-processing. PR #1025 (J. Molson, F. Schmidt)
* Fixed a difference in sign in the thin combined function code. PR #1005 (T. Persson)

**User Side Changes**

* Minor changes to the formatting of `fort.18` to ensure correct column width. PR #1029 (J. Molson, F. Schmidt)
* Added better error reporting for the FLUKA interface. PR #1028 (J. Molson)
* Particles that pass through a collimator, but don't interact with it, no longer have their coordinates changed. Previously, these particles were shifted to the closed orbit, and had their units changed, for then to be changed back after the collimator. This added unnecessary numerical noise. PR #1023 (V.K. Berglyd Olsen)
* The scatter module has been rewritten, and the PYTHIA integration updated to work with PYTHIA 8.243. This version allows for sending tracked particles, in terms of their momentum vector, as well as a sampled colliding particle, to the event generator, This means we no longer have to project a PYTHIA event onto the SixTrack particle, but can instead send the tracked particles back and forth between the codes. The Scatter and Pythia modules have been extended to allow for this. In addition, more density profile options have been added: collision with a reference particle at a given probability density, or a model of Beam 2 given by a set of Twiss parameters. The latter can also be set up to mirror Beam 1 using the internal optic parameters calculated by SixTrack. PRs #1017 and #1038 (V.K. Berglyd Olsen)

**Test Suite**

* Added a number of SixDA test to increase coverage of the Differential Algebra version of SixTrack. PRs #1030, #1032 (F. Schmidt, V.K. Berglyd Olsen)

**BOINC Interface**

* Extended the BOINC interface code to produce a new validation file for BOINC jobs. The new file supports a wider range of job types, especially jobs that do not produce a `fort.10` file. PR #878 (V.K. Berglyd Olsen, A. Mereghetti)

**Code Improvements and Changes**

* Cleaned up a number of unused variables throughout the SixTrack source code. PR #1037 (J. Molson)

### Version 5.4.2 [22.11.2019] - Release

**Bugfixes**

* Fixed a missing use statement in collimation under the `ROOT` compiler flag. Root support should now build again. PR #1015 (V.K. Berglyd Olsen, J. Molson)
* Removed one sigma cuts on all normal random distributions in the crystal module. These were ported over from the old version of the crystal collimation code, but never worked properly there due to a datatype bug. There should not be any such cuts in the physics in the first place, so they have now been removed. PR #1016 (M. D'Andrea, V.K. Berglyd Olsen)
* Fixed an inconsistency in the splitting tool for `singletrackfile.dat` where the header was not padded with zeroes when files were split up into particle pair files. The padding was always in the post-processing code, but not in the read90 tool which the splitting tool was derived from. Since the test suite relies on the read90 tool, the split files still passed, but user analysis codes failed. The splitting tool is now consistent with SixTrack post-processing. PR #1013 (V.K. Berglyd Olsen)

**User Side Changes**

* The turn number in the particle state files is now the smallest of current turn and number of particle tracked turns. Previously, the default value was always the total number of turns, and only adjusted if a particle was lost before the end of tracking. This was misleading users to think that the `initial_state.dat` file wasn't really written before tracking, as the turn number was non-zero. PR #1013 (V.K. Berglyd Olsen)
* The main crystal interaction file is now written only if the `WRITE_CRYCOORD` flag is enabled in collimation. The entrance and exit files are, in addition, only written when the global `DEBUG` flag is enabled in the `SETTINGS` block. PR #1019 (V.K. Berglyd Olsen, M. D'Andrea, A. Mereghetti)
* Added a command line flag `--notrack` that disables tracking in SixTrack. This allows for the checking of input files and simulation initialisation without running the full job. PR #1020 (V.K. Berglyd Olsen)

**Documentation**

* Documentation section for collimation updated to include crystal collimators. PRs #1012 and #1019 (M. D'Andrea, V.K. Berglyd Olsen)

### Version 5.4.1 [01.11.2019] - Release

**Bugfixes**

* Fixed an issue when opening `fort.208` when building with both FLUKA and G4COLLIMATION compiler flags at the same time. PR #1008 (V.K. Berglyd Olsen, J. Molson)
* Fixed a bug in the efficiency calculation where the vertical sigma was not correctly calculated when filling the histograms due to a missing set of parentheses. PR #987 (V.K. Berglyd Olsen)

**User Side Changes**

* The collimation module now supports crystal collimation. The crystal collimation code that was written based on an older version of SixTrack has now been updated and added to SixTrack 5. Using crystal collimators requires the use of the new collimator database format. PR #1004 (M. D'Andrea, V.K. Berglyd Olsen)
* The collimation code is no longer hardcoded to assume LHC element naming convention. The collimator names are still processed assuming LHC convention when using the old input block and database format, but the new format does not rely on any naming convention. A collimator is a collimator if it's defined in the collimator database file. PR #970 (V.K. Berglyd Olsen)
* Since the collimators are no longer dependent on naming convention, the collimator stage (primary, secondary, tertiary, etc) has to be specified in the database if an on-line analysis based on these are intended. This has been added as an optional parameter in the collimator family definition. PR #986 (V.K. Berglyd Olsen)
* Collimator efficiency studies (histograps on normalised amplitudes) are now disabled by default, and need to be explicitly requested. A flag has been added to both the old and the new input block for this. The feature is by default disabled as it requires substantial amounts of memory relative to the total memory usage of SixTrack. PR #987 (V.K. Berglyd Olsen)
* Aperture tilt in the input is now specified in degrees, not radians. PR #1003 (T. Persson)

**Code Improvements and Changes**

* This release includes a major clean-up of the collimation code. The old collimation module has been split up into separate modules. The K2 physics routines have been split out into their own module, and connected to the crystal collimation module. The funlux code has also been moved to a new module. A lot of statistics calculations and other sections of code have been extracted and put into separate routines to make the code both more readable and to avoid duplication of code. Many unused or unneeded variables have been removed, and collimator materials have been split out into a separate module; so has a lot of shared collimator variables. A test for the `tracks2.dat` file has been added. PRs #975, #987, #997, #998 and #1002 (V.K. Berglyd Olsen, A. Mereghetti)
* Exact drifts are now computed without first converting the coordinates to m and rad and then back to mm and mrad for each drift. This adds numerical noise and unnecessary CPU cost. PR #999 (V.K. Berglyd Olsen)
* Cavities and phase trombone elements now uses the general routine for updating particle energy arrays, ensuring consistency between elements. PR #1001 (V.K. Berglyd Olsen)

### Version 5.4 [10.10.2019] - Release

**Bug Fixes**

* A series of long standing issues caused by uninitialised variables have been resolved. PRs #983 and #988 (V.K. Berglyd Olsen).<br>The following modules and settings were affected:
  * Beam--beam simulations with the `ibeco` flag set to 0. This caused `NaN` particle coordinates on some systems.
  * In the differential algebra version of SixTrack, the longitudinal part of the normalisation matrix was uninitialised when running in 4D, but the values still used in some calculations, causing `NaN` values on some systems.
  * When running SixTrack with the `ntwin` parameter set to 1, post-processing would still compute the amplitude for the second particle from values not being initialised. Now, if `ntwin` is not set to 2, these values are set to zero and the extra calculations skipped.
  * A number of variables were uninitialised in parts of the initialisation code when running 6D simulations (in subroutine `umlauda`), these have been cleared up, but are not known to have caused any issues.
* Fixed a bug in collimation where the collimator families were not generated when using the old database file format with the new `COLL` block format. PR #984 (V.K. Berglyd Olsen, M. D'Andrea)
* Fixed a bug in `aperture_losses.dat` where the header was missing the `#` char so that it can be identified as a comment line for analysis code. PR #995 (A. Mereghetti, V.K. Berglyd Olsen)
* Fixed an issue with writing the header of `singletrackfile.dat` where SixTrack would skip to after tracking if it failed to write to file. PR #996 (V.K. Berglyd Olsen, A. Mereghetti)
* Added a check of `iostat` when closing file units in the internal file units handler module. Any error is reported to stderr and to log file. PR #996 (V.K. Berglyd Olsen, E. Mcintosh)

**User Side Changes**

* The STF build flag has been removed. That means SixTrack now always produces a single track file `singletrackfile.dat` instead of the optional pair track files `fort.59` - `fort.90`. A tool for converting the full track file to a pair track file has been added. See the user manual for further details. PRs #967 and 989 (V.K. Berglyd Olsen, K.N. Sjobak)
* A random numbers module has been added, which introduces the `RAND` block in the input file. The block gives more control over how the internal random number generators are initialised and also provides a better framework for the internal management of the various random number sequences used by different modules of the code. It is also designed to be easier to manage when using checkpoint/restart. Many modules still use their own seed, so this will be implemented gradually. Currently, the new random numbers module is only used by the `DIST` block. PR #978 (V.K Berglyd Olsen)
* The main debug file `dynksets.dat` in `DYNK` is no longer written by default. Previously, this file could be disabled with the `NOFILE` flag in `fort.3`. Instead, it now has to be explicitly requested with the `DYNKSETS` flag. PRs #992 and #993, solving issue #982 (V.K. Berglyd Olsen)
* The `TILT` build flag has been removed. The feature is now always enabled. PR #985 (V.K. Berglyd Olsen)

### Version 5.3.4 [01.10.2019] - Release

**Bug Fixes**

* Fixing an error in the integration of the electron lens radial profiles from files. The error was twofold: (1) it had the wrong unit conversion from A/cm2 to A/mm2, and (2) integration was not deploying the rule of the trapezoid on the radius, but the lower square. PR #968 (A. Mereghetti)
* Minor fix where the collimation exit routine was not called if SixTrack was run in thin 4D mode. This only prevented some final summary files to be written, and had no effect on the simulation itself. PR #976 (V.K. Berglyd Olsen, M. D'Andrea)

**User Side Changes**

* Onesided collimators can now be set in the new collimation database. The LHC hard coded naming convention restrictions have been removed for the new database format, but are still in place for the old format. PRs #958 and #966 (V.K. Berglyd Olsen, A. Mereghetti)
* The option to slice the collimator jaw and apply a deformation has been added to the new collimator database format. The previous restrictions on LHC naming convention have been removed, and the slicing can now be applied to any collimator, with different fit models as needed. PR #969 (V.K Berglyd Olsen, A. Mereghetti)
* The number of beam--beam elements now scale dynamically with user input. It was previously restricted to 500 elements. PR #974 (V.K. Berglyd Olsen)

**Test Suite**

* A series of tests of the collimator jaw profiling has been added to ensure the code is stable through development stages. An additional output file for the jaw was added as well. PR #961 (A. Mereghetti)

**Code Improvements and Changes**

* The collimation jaw fit code has been cleaned up, improved for speed, and moved to a separate module. PRs #933 and #966 (V.K. Berglyd Olsen)
* The separate tracking initialisation routines for thick and thin tracking have now been merged into one routine in a dedicated tracking module. PR #977 (V.K. Berglyd Olsen)

### Version 5.3.3 [09.09.2019] - Release

**Bug Fixes**

* Fixed bug in the `DIST` module with conversion of longitudinal emittance from eVs to µm. The conversion was off by a factor `1e6` due to the energy variable being in `MeV` not `eV`. PR #950 (V.K. Berglyd Olsen)
* Fixed bug in the `DIST` module where emittance was sent to `DISTLIB` as normalised emittance, while `DISTLIB` expected geometric emittance. PR #952 (V.K. Berglyd Olsen)
* Fixed a bug in binary particle state files where the internal normalisation matrix was written to file instead of the ones which have all elements scaled to the same unit. PR #952 (V.K. Berglyd Olsen)

**Documentation**

* Minor changes to the documentation (LaTeX manual and markdown files) to correct outdated information. PR #948 (V.K. Berglyd Olsen)
* Updated the README to add a paper that can be cited by studies using SixTrack 5. PR #949 (V.K. Berglyd Olsen, R. De Maria)

**Code Improvements and Changes**

* The way SixTrack keeps track of particle pairs for DA studies has changed. Earlier, the pairing was preserved by a reverse map from original particle index to the current index. The index of a particle changes when it is lost in an aperture, collimator or interaction point. The reverse map was effectively a record of the particle ID. However, since the new `DIST` module makes it possible to set the particle ID to any value, the reverse map has now been replaced by a map containing a particle pairID as well as whether it is particle 1 or 2 of the pair. This map is not under the user's control, and therefore preserves the pair structure through tracking. Rewriting the code to use this map eliminates a potential memory access violation due to a corner case when particles are lost and initiated with a non-incremental particle ID. The rewrite is otherwise identical to old functionality, and shouldn't alter any results. The main benefit is cleaner code, and the particle ID now being entirely passthrough as far as SixTrack is concerned, making it easier to interface with external codes that inject new particles into SixTrack. PR #938 (V.K. Berglyd Olsen, A. Mereghetti)
* Added a check in the parsing of multi-column `STRUCT` input blocks that ensures that the element position of a given lattice element has a position larger or equal to the previous element. This prevents the accidental initialisation of negative length elements. PR #955 (V.K. Berglyd Olsen)
* Moved a number of subroutines related to initialisation of beam--beam elements and elements in general out of the main `sixtrack` source file and to more appropriate files. PRs #957 and #960 (V.K. Berglyd Olsen)

### Version 5.3.2 [23.08.2019] - Release

**Bug Fixes**

* Fixed a bug with `DUMP` format 101 when using HDF5 output. The memory map used was mixed up with the map for format 3. PR #937 (V.K. Berglyd Olsen)
* Fixed a bug in the `DIST` block module where reading less than all particles of a file would fail. PR #939 (V.K. Berglyd Olsen)
* Fixed an issue where calculating PDGID would overflow due to intermediate integer variables being 16 bit. PR #940 (J. Molson)
* Fixed a bug with saving int16 ion variables to HDF5 files. PR #942 (V.K. Berglyd Olsen)

**Documentation**

* Added full documentation of the new and improved `DIST` block. PR #941 (V.K. Berglyd Olsen)

**Code Improvements and Changes**

* General particle transport for FLUKA and changing FLUKAIO to a submodule pulled from the GitLab repository (requires CERN Kerberos access). PR #919 (J. Molson)
* The submodule for libArchive, and its interface and wrapper code, has been removed. PR #920 (V.K. Berglyd Olsen)
* Clean-up of the formatting of the header in a number of output files in the aperture module. PR #923 (J. Molson)
* The arrays for the per-element normalisation matrix is now no longer sparse, but instead using a compact array of structs. This reduces the memory usage of SixTrack by up to about 35%. PRs #934 and #935 (V.K. Berglyd Olsen)
* The linear optics subroutines have been moved the a new module together with the parsing of the `LINEAR OPTICS` input block. PR #936 (V.K. Berglyd Olsen)

**Test Suite**

* Removed unused and incomplete test `thick6dsingles` and duplicate test `thick4in_da` as well as a lot of comments in the code used as "version control". PR #944 (V.K. Berglyd Olsen)

### Version 5.3.1 [02.08.2019] - Release

**Bug Fixes**

* Fixed and issue when using the collimation module with a thin 4D simulation. In this setup, the module would not be properly initialised due to an erroneous if-condition for the initialisation call. PR #931 (V.K. Berglyd Olsen)
* Fixed a minor issue with the formatting of the tracking progress printout. PR #925 (V.K. Berglyd Olsen)
* Fixed and issue with missing labels in aperture losses file. PR #928, issue #926 (A. Gorzawski, A. Mereghetti)

**Documentation**

* Some inconsistencies and out-of-date information has been corrected in the user manual. PRs #921, #922 and #927 (V.K. Berglyd Olsen, R. De Maria)

**Code Improvements and Changes**

* The `DIST` block has been rewritten and a number of new parsing options added for integrating with a new external library for generating beam distributions. The library is not yet completed, so the new block format is not finalised or documented. However, the `DIST` block is backwards compatible with the old options, and should be working as before. PRs #905 and #930 (V.K. Berglyd Olsen, T. Persson)
* Particle spin arrays have been added to SixTrack intended for future code, but not yet in use. The arrays have been added nonetheless so they can be included in the new `DIST` block. PR #916 (J. Molson)

### Version 5.3.0 [11.07.2019] - Release

**Bug Fixes**

* The Fortran `.eqv.` operator has a lower precedence than `.eq.`. This was not accounted for in the `SETTINGS` block when the particle summary output after tracking is requested. PR #902 (V.K. Berglyd Olsen, A. Mereghetti)
* Fixed an issue where the wrong index value was written to the `fort.208` file used in the FLUKA coupling, making the file useless. PR #911 (A. Gorzawski, A. Mereghetti)

**User Side Changes**

* A new Collimation Database file format is now supported by the Collimation module. The new format is column-wise as opposed to the old single column file. The old database is converted for the user to the new format and written to a file with the same name as the old database but `.new` added to it. The new format is in preparation for adding new collimator types to the database. PR #903 (V.K. Berglyd Olsen, A. Mereghetti)
* A Geant4 block `GNT4` has been added in order to use Geant4 with the collimation module. The features are enabled by building with the `G4COLLIMATION` flag. PR #713 (J. Molson, A. Mereghetti, K. Sjobak)
* SixTrack can now track any particle with a charge. The charge setting that exists in the `SIMU` block, and has been added to the `HION` block, now properly separates the particle charge from the ion Z value. PR #713 (J. Molson, A. Mereghetti, K. Sjobak)

**Build System**

* Fixes to ARM and OSX builds. PR #901 (J. Molson)
* Removed the .exe in the middle of the file name when building on Windows. PR #900 (V.K. Berglyd Olsen)
* The `ZLIB` flag is now on by default, meaning the `ZIPF` block also works by default. PR #899 (V.K. Berglyd Olsen)
* The NAFFlib submodule is updated by default during build if the user has git installed. PR #899 (V.K. Berglyd Olsen)
* Some minor tweaks to the build system were made to support AVX-512. Seems to work best when building with Intel's Fortran compiler. PR #898 (V.K. Berglyd Olsen)

**Code Improvements and Changes**

* The `FAST` compiler flag and code has been removed. PR #914 (R. De Maria)
* The collimation beam distribution generator has been moved to a separate module named `coll_dist` and undergone a significant cleanup. A number of now redundant particle arrays have been removed in the process. PR #885 (V.K. Berglyd Olsen)

### Version 5.2.10 [13.06.2019] - Release

**Bug Fixes**

* Technically not a bug, but the Chromaticity Corrections (`CHRO`) block did not have proper checks on the validity of the settings in the input block. It never reported on the non-existence of the elements requested (although there are checks later on in SixTrack) and if the `ichrom` flag was set to an invalid value, it was simply set to a valid value without telling the user. These issues have now been corrected. It may cause previously accepted input files with faulty parameters to now fail. PR #894 (V.K. Berglyd Olsen)

**User Side Changes**

* Added `ZLIB` as a build option, which now includes a much simpler zipping library than libArchive. The `ZLIB` flag is off by default, but when enabled, will allow the user to use the `ZIPF` module to add simulation files to a zip file. The `ZIPF` module has also been improved a bit, and the output file can be specified as well as the compression level. In the near future, libArchive will be removed, and `ZLIB` enabled by default. PR #882 (V.K. Berglyd Olsen)
* Rotations (also known as patches) have been implemented as new elements. It is following the same logic as MAD-X, and the conversion from MAD-X is already in the master branch. PR #892 (T. Persson)
* The particle state files that can be dumped before and after tracking now include the charge column if ions are enabled. The file also includes additional values in the header covering the settings for the reference particle, the 4D and 6D closed orbit, the tunes, and the 6-by-6 TA matrix (eigenvector/normalisation matrix). PRs #893, #895 and #896 (V.K. Berglyd Olsen, R. De Maria)

**Build System**

* Removed the `API` build flag for building the BOINC API. Building without it would previously build a dummy API for BOINC. There is no longer any apparent reason to do this as the rewrite of checkpoint/restart ensures that the tests can run fine with the full API linked. The API is now always linked when the `BOINC` flag is enabled. PR #890 (V.K. Berglyd Olsen)

### Version 5.2.9 [05.06.2019] - Release

**Bug Fixes**

* Added a check that the `BEAM` block is present when beam--beam elements exist in the lattice. This would previously run, but produce NaN particle coordinates due to division by zero. Thanks to Sofia Kostoglou for the example job. PR #887 (V.K. Berglyd Olsen, R. De Maria)
* Fixed a minor issue with the MD5 interface when running on Windows where the return from `fwrite` was treated as an error value. This resulted in an error message being printed for each line of the hashed file for no reason at all. PR #886 (V.K. Berglyd Olsen)

**User Side Changes**

* The Collimation Module now allows for name/value format to be used in the `COLL` block in `fort.3`. The full description of the formatting of the block is available in the user manual. The old format is still supported, but they cannot be mixed. PR #796 (V.K. Berglyd Olsen)

**Code Improvements and Changes**

* The timing module `mod_time` is now fully checkpointed, meaning the timing data in the log file `sim_time.dat` contains the correct timing information for the full simulation across multiple checkpoints and restarts. PR #884 (V.K. Berglyd Olsen)

**Build System**

* Solaris build support in the BOINC build script hs been added. PR #888 (J. Molson)


### Version 5.2.8 [28.05.2019] - Release

**Bug Fixes**

* Fixed a bug in the routine that inserts structure elements into the lattice. The routine uses a circular shift to cycle a fresh element (added at the end of the lattice structure by memory allocation) into the new position, but was cycling the wrong way. The result was that if the shift was over more than 2 elements, the order of the lattice would be wrong. Bug experienced when inserting the aperture markers downstream of Fluka insertions or of the lattice edges. PR #879 (A. Mereghetti)
* Fixed an error message in `DUMP` that would segfault due to using the wrong index variable in an array lookup. PR #879 (A. Mereghetti)
* Fixed a bug in the Fluka-SixTrack coupling when SixTrack errors while the beam is on the Fluka side. Solved with a temporary call to `kernel_fluka_exit` at `prror`. This required to change the interface of the `kernel_fluka_exit` subroutine. PR #879 (A. Mereghetti)
* Fixed checkpointing of statistical scaling arrays in the Scatter Module. These were previously not checkpointed at all. PR #880 (V.K. Berglyd Olsen)
* Fixed a bug in the SixTestWrapper where not passing any checks at all would default the SixTestWrapper to pass the test as a whole. It now properly requires at least one check to be considered passed. Previously, one of the tests would always pass without verifying a single output file. PR #881 (V.K. Berglyd Olsen)

**BOINC Integration**

* The dummy API has been replaced by a BOINC service module. PR #875 (V.K. Berglyd Olsen, A. Mereghetti)
  * When running on BOINC, SixTrack no longer checkpoints at the turn interval specified in `fort.3`, but instead tries to checkpoint once a minute. The volunteer's setting on minimum checkpoint interval in their BOINC Manager is still obeyed.
  * When BOINC is run in standalone mode, that is in the test suite rather than by a volunteer, the checkpoint interval is 10 seconds. The settings for the `CRKILLSWITCH` are still obeyed.
  * For the time being, at least, the checkpointing decisions taken by the service module are logged to the file `cr_boinc.log`. This log file may be removed when the new scheme is properly tested.
* The BOINC library is now built from a point on the upstream master branch after some modifications to the Fortran API were merged. This means that we are running off-tag on the BOINC API, which now reports version 7.15. PR #877 (V.K. Berglyd Olsen)
* Tracking progress reported to BOINC now stops at 99%, and the remaining 1% is reserved for post-processing. PR #877 (V.K. Berglyd Olsen)

**Code Improvements and Changes**

* Cleanup in the FLUKA integration modules, merging them to a single module. At the same time, refactored code for checking integrity of user info concerning coupling. Also added state variables concerning the current fluka insertion and the last messages sent/received, in view of a thorough revision of the error handling. PR #879 (A. Mereghetti)

### Version 5.2.7 [17.05.2019] - Release

**Bug Fixes**

* Fixed an infinite loop corner case from previous release. It was only triggered if one tried to rum the Differential Algebra (DA) version of SixTrack with checkpoint/restarting, which isn't something that makes sense doing anyway. The loop was only triggered with ifort builds, and used all available memory when triggered, so it was a bit nasty. In addition, the DA executable is no longer built when building with checkpoint/restart. PR #871 (V.K. Berglyd Olsen)

**User Side Changes**

* A new block has been added to `fort.3`. It's named `SIMU`, and holds the main tracking variables for setting up a simulations. When used, it replaces the `TRAC`, `INIT`, and `HION` block. It is intended to work with the `DIST` block, and therefore does not support any of the settings for generating distributions. If those are needed, please use the old blocks. For further details, see the user manual. The block is considered experimental until it's thoroughly tested. PR #872 (V.K. Berglyd Olsen, R. De Maria, A. Mereghetti)

**Documentation**

* Fixed the RF-multipole definitions. PR #869 (R. De Maria)

**Code Improvements and Changes**

* Removed the thick arrays: `ekv`, `fokqv` and removed the checkpointing of all remaining thick arrays. These are recomputed by the subroutine `synuthck` anyway. PR #870 (V.K. Berglyd Olsen, K. Sjobak)
* Some cleanup in the main common module for particle arrays. Mostly added comments on what variables do, and removed a few unused or redundant arrays. The deleted arrays are: `xsv`, `zsv`, `dp0v`, `sigmv6`, `dpsv6`, `xlv`, `zlv`, `nms`, and `ekkv`. PR #874 (V.K. Berglyd Olsen)

### Version 5.2.6 [07.05.2019] - Release

**Bug Fixes**

* Fixed a bug with linking zlib and libarchive when building on Windows. PR #853 (V.K. Berglyd Olsen)
* Fixed a bug in postprocessing where the binary postprocessing summary file would append its data to another file in some cases. This was caused by the file having a fixed file unit in the range which is otherwise reserved for dynamically allocated file units. This file now also gets a unit assigned to it, avoiding this issue. Some checks have been added to the file units module to try and prevent similar bugs in the future. PR #855 (V.K. Berglyd Olsen)
* The recently added error tests failed on some operating systems due to the python wrapper not finding a symlinked include file. PR #862 (K.N Sjobak)
* Fixed a floating point exception and array bounds violation in the fringe field module. PR #866 (V.K. Berglyd Olsen)
* Fixed a bug where sometimes a segfault would be triggered when SixTrack exited with an error due to the error routine itself missing an explicit interface. PR #866 (V.K. Berglyd Olsen, K. Sjobak)

**User Side Changes**

* All the checkpoint/restart files have now been given more descriptive names with a `cr` prefix. The log file, formerly `fort.93`, has been renamed `cr_status.log`. The output to the log file has also been cleaned up significantly so that it is now more readable. PR #854 (V.K. Berglyd Olsen)
* The `numlmax` parameter in the `TRAC` block of `fort.3` has been removed. The feature was either not fully implemented as described, or has been broken at some point. The option to use checkpoint/restart files to extend earlier simulations will be reimplemented at a later point. PR #854 (V.K. Berglyd Olsen)

**Code Improvements and Changes**

* When running SixTrack with BOINC, checkpoint files are written on turn 1 without requesting permission from the BOINC api. This ensures that the tests that require a successful restart to pass also pass on executables with the BOINC api linked. PR #853 (V.K. Berglyd Olsen)
* Rewritten the checkpoint/restart module and cleaned out all code duplication. It should now be a lot easier to read the code and extend it. PR #854 (V.K. Berglyd Olsen)

### Version 5.2.5 [26.04.2019] - Release

**User Side Changes**

* It is now allowed to set `R1=0.0` for both elenses and chebyshev maps, even though the fox implementation is still on-going. This will allow for running full elens studies, even though the closed orbit does not take into account these elements yet. PR #850 (A. Mereghetti)
* The main SixTrack executable will now print version information and exit if provided with one of the arguments `-v`, `-V` or `-nv`. Representing two levels of detail, and the numerical version number, respectively. PR #846 (V.K. Berglyd Olsen)
* Alternative file names for the `fort.3` and `fort.2` input files can be given as first and second argument, respectively. PR #846 (V.K. Berglyd Olsen)

**Bug Fixes**

* Setting the `TRAC` block variable `imc` to anything other than 1 now properly triggers an error. This is now consistent with the manual, which states "Number of variations of the relative momentum deviation has been removed. This value must be 1." PRs #847 and #848 (V.K. Berglyd Olsen)
* The wrapper for CRLIBM had its own set of parameters `pi` and `pi2`. The wrapper now uses the new hex constants defined in `numerical_constants`. PR #845 (V.K. Berglyd Olsen, Eric Mcintosh)

**Test Suite**

* Added a new test category that checks error messages when SixTrack is made to fail by providing invalid input files. This bot increases code coverage, and ensures that faulty simulation input is caught and reported correctly. PRs #825, #847, #848 and #851 (V.K. Berglyd Olsen, K. Sjobak)

### Version 5.2.4 [23.04.2019] - Release

**User Side Changes**

* Error messages are now written to `stderr` rather than `stdout`. When building with checkpoint/restart, this is piped to `fort.91`. PR #834 (V.K. Berglyd Olsen)
* The tracking progress output to `stdout` now also contains info about how many particles are being tracked. This is useful in simulations with particle losses. This information was previously only printed by the collimation module. PR #836 (V.K. Berglyd Olsen)

**Bug Fixes**

* Minor inconsistency in the writing of the `initial_state.bin` and `final_state.bin` files between compilers as logical values are stored differently. This is now set explicitly for logical values. PR #832 (V.K. Berglyd Olsen)
* Electron lenses/chebyshev lenses flagged in `fort.2`, but not declared as such in `fort.3`, were caught only in case the respective blocks were active, The consistency check was applied in the wrong stage of parsing, and would therefore not necessarily catch the inconsistency, causing a segfault. Fixes Issue #826. PR #833 (A. Mereghetti)

**Code Improvements and Changes**

* Removed 30 arrays used as temporary arrays for thick tracking (but was always allocated). These arrays were allocated to the number of particles, but could be replaced by scalars. This change frees up `(NPART-1)*240` bytes of memory, where `NPART` is the number of tracked particles. PR #822 (V.K. Berglyd Olsen, K. Sjobak)
* Some minor reordering of the logic in handling RF cavities when below transition energy. The initialisation of the cavity elements is now in one place only. An integer array with the length of the single element list (`NELE`) was also removed. Tests have been added for below transition energy tracking. PR #828 (V.K. Berglyd Olsen, K. Sjobak)
* Commented out and no longer maintained code under the `DEBUG` flag has been deleted in coordination with Eric McIntosh. PR #835 (V.K. Berglyd Olsen, K. Sjobak)
* The constants `pi`, `pi/2`, `2*pi`, `sqrt(pi)`, `pi/180`, and `1/ln(2)` are now set with hex values to bypass conversion from decimal to binary. This is enabled for single and double precision, and for quad precision for the GNU compiler only. Intel and NAG does not support this for quad precision. PRs #840 and #842 (V.K. Berglyd Olsen)

### Version 5.2.3 [13.04.2019] - Release

**Bug Fixes**

* Fixed a bug in collimation where the delta_p and momentum of particle distributions generated by the collimation module were not properly recalculated from the energy. Only the first particle was handled, and the rest retained the values from the distribution generated from the `TRAC`/`INIT` blocks. The bug originated in version 5.0.3. PR #818 (V.K. Berglyd Olsen, A. Mereghetti)
* Fixed a bug where a negative kz of -12 (cavity element when below transition energy) was not converted to positive value if the print flag was off. Neither was the phase converted to radians. The conversion was embedded in a loop otherwise used for printing the lattice. PR #823 (V.K. Berglyd Olsen)
* Fixed a bug in the SixTestWrapper in release 5.2.2 where extra checks would always be marked as passed. PR #818 (V.K. Berglyd Olsen)

**Code Improvements and Changes**

* The collimator database has been moved to a new module `coll_db`. This is in preparation for moving to a new and more flexible collimator database format. PR #792 (V.K. Berglyd Olsen, A. Mereghetti)

### Version 5.2.2 [08.04.2019] - Release

**User Side Changes**

* Support for variations of momentum offset in the `INIT` block has been removed. This feature has not been maintained for a long time, and was incompatible with the several modules of SixTrack. It was also wasteful in terms of memory usage. Since SixTrack accepts input distributions from file, the feature is also redundant. We therefore decided to remove it rather than bringing it up to speed with the rest of the code. This change fixes the `imc` flag to a value of 1 (7th value on line 1 of the `INIT` block). Any other value will cause an error. This change frees up `(NPART-1)*680` bytes of memory, where `NPART` is the number of tracked particles. PR #804 (V.K. Berglyd Olsen)

**Test Suite**

* The output of the wrapper executable for the test suite has been made a bit more readable. This is visible when ctest is run with the `-V` flag. PR #807 (V.K. Berglyd Olsen)

### Version 5.2.1 [05.04.2019] - Release

**New Features**

* The Structure Block in `fort.2` is now also available in a multicolumn format. This mode is switched on by specifying `MULTICOL` on the first line. The column format consists of a minimum of 3 columns. The first being the proper element name as used in MadX, the second column is the Single Element name (corresponding to the single column format), and the third one is the element position as described in MadX (the s-coordinate of the centre of the element). A `multicol` flag has also been added to the SixTrack converter in MadX, and will be available in the next release of MadX following version 5.04.02. The new format is currently only implemented in the aperture losses file to list the proper element name for lost particles. PR #799 (V.K. Berglyd Olsen)
* A module for quadrupole fringe fields has been added to SixTrack, and added as an `FFIE` block in `fort.3`. This method allows for the usage of a longitudinal description of the quadrupole magnetic field, adapted for each magnet specifically selected for the study, without changing the reference optics of SixTrack. PR #776 (T. Pugnat, B. Dalena)

**User Side Changes**

* The `PRINT_DCUM` flag in the `SETTINGS` block that printed the full lattice with s-coordinates to stdout is now instead written to a file named `machine_length.dat`. The file now also prints the information from the multicolumn lattice description in the Structure Block, and compares the computed position `dcum` with the one read from MadX. The delta in nanometres is listed. PR #799 (V.K. Berglyd Olsen)

**Bug Fixes**

* Fixed a bug in `plato_seq.f90` where the comment states in `tuneffti` that it will not accept frequency indices of 0 or 1, but does anyway and returns a value NaN. This causes a segfault on Debug type builds. The routine now returns one instead in these cases. A similar fix was applied to `tunebt2`. PR #803 (V.K. Berglyd Olsen)
* Removed two `close(6)` calls in `abend` that caused an abort with BOINC when building with nagfor. PR #803 (V.K. Berglyd Olsen)

### Version 5.2 [27.03.2019] - Release

**New Features**

* The SCATTER module has been completely rewritten to better allow for adding multiple scattering processes (generators) to scatter elements. It does this by computing branching ratios either from provided cross sections, or by fixed values. PR #670 (V.K. Berglyd Olsen)
* SixTrack can now integrate with PYTHIA to generate scattering events for the SCATTER module. Currently only head on events at a single centre of mass energy are generated. This will be extended further in the future. It is, however, possible to extract elastic and diffractive events, with and without particle losses. PR #670 (V.K. Berglyd Olsen)
* The `FINALSTATE` and `INITIALSTATE` flags in the `SETTINGS` block now take `ions` as a second keyword, enabling the dumping of the ion columns in addition to the main particle arrays. PR #777 (V.K. Berglyd Olsen)

**User Side Changes**

* When reading `fort.13` the particle energy is used to set the momentum and delta_p arrays instead of the delta_p overwriting energy and momentum. This means the delta_p values in the file are ignored. The change was made due to the intended usage of `fort.13` being continuation of tracking, not initialisation of particles. During tracking only energy should be changed, and the other dependent variables calculated from the energy arrays. PR #766 (V.K. Berglyd Olsen)
* The trombone and beam--beam elements now use the updateEnergy routine from `mod_particles` as well to recompute the energy-dependant arrays. The change in beam--beam slightly alters output due to the previous code using a different way of calculating one of the values. PRs #786 and #795 (T. Persson)

**Bug Fixes**

* The module handling opening and closing of files wrote warnings to stderr. This is now written to stdout instead to avoid cluttering the stderr logs when running on BOINC. PR #774 (V.K. Berglyd Olsen)
* Fixed a bug with some some of the output files under the FLUKA flag. The bug was introduced in 5.1.2. PR #780 (V.K. Berglyd Olsen)
* Fixed some `write(unit,*)` of header files in the collimation module that produced different results on the `nagfor` compiler than the others. PR #790 (V.K. Berglyd Olsen)

**Code Improvements and Changes**

* Interface routines for adding attributes to HDF5 datasets have been added. PR #670 (V.K. Berglyd Olsen)
* Longer filenames are now allowed in the `mod_units` module. The new maximum is 255 characters (Windows limit). PR #783 (V.K. Berglyd Olsen)
* Moved the parameters derived from pi `twopi`, `pi2`, `rad`, and `pisqrt` to module `numerical_constants`. PR #785 (V.K. Berglyd Olsen)
* Added echo on STDOUT of calls to `contour_aperture_marker` and `contour_aperture_markers` subroutines and respective header, improved dump of aperture markers and header now fully aligned to dumped table, improved readability of messages on STDOUT from `contour_aperture_marker`, and added a function that checks that the lattice structure does not start inside a `FLUKA` insertion. PR #789 (A. Mereghetti)

**Documentation**

* Fixed an issue in the documentation of `DYNK` where a table was split in two. PR #782 (T. Persson)

**Build System and Test Suite**

* Added a `MOSYMLINK` flag to CMake that makes the test suite copy the input files to the build folder rather than symlink on Unix environments. PR #781 (V.K. Berglyd Olsen)
* Some minor cleanup CMake, including the way version numbers are handled. PR #788 (V.K. Berglyd Olsen)
* Disabled collimation tests when checkpoint/restart is active. PR #791 (V.K. Berglyd Olsen)

### Version 5.1.3 [25.01.2019] - BOINC Release

**User Side Changes**

* Misalignment of RF-multipoles is now possible. PR #763 (T. Persson)

**Bug Fixes**

* Aperture checking along a transition takes into account RACETRACK properly. PR #771 (A. Mereghetti)
* When built with checkpoint/restart, the `mod_units` would try to write to the `file_units.log` file after it was closed. This caused SixTrack to exit with an error code after tracking when building with the NAG compiler. PR #770 (V.K. Berglyd Olsen)

**Documentation**

* Build scripts for HTML version of the manual have been fixed. PR #767 (V.K. Berglyd Olsen)
* Documentation on aperture limitations have been updated. PR #771 (A. Mereghetti)

**Code Improvements and Changes**

* The old `comnul` subroutine used to zero out old common blocks has been reduced to almost nothing. The default values are now set in their respective modules, mainly in `common_modules.f90`. PRs #765, #768 and #771 (V.K. Berglyd Olsen, A. Mereghetti)
* Some general code clean-up of the aperture module. PR #771 (A. Mereghetti)
* Symplecticity deviation is now written to `sim_meta.dat` instead of its own file. PR #770 (V.K. Berglyd Olsen)
* Files are no generally opened where they're needed rather than at the start of `maincr`. This helps to reduce the number of zero size fort.* files in the simulation folder. PR #770 (V.K. Berglyd Olsen)

### Version 5.1.2 [14.01.2019] - Release

**User Side Changes**

* The memory allocation log is now only written when SixTrack is built with the `DEBUG` flag. This file can grow fairly large, and has caused issues when running on BOINC. PR #750 (V.K. Berglyd Olsen)
* Added the option to write a file with the initial (pre-tracking) coordinates of all particles. This is the same format as the final state files added earlier. The feature is enabled with the `INITIALSTATE` and `FINALSTATE` keywords in the `SETTINGS` block. PR #760 (V.K. Berglyd Olsen)
* Generalised RF-Multipoles have been added to SixTrack. See manual. PR #756 (T. Persson)

**Bug Fixes**

* An allocatable array was declared and passed to `mod_alloc` per particle/turn/element in subroutine `sbc` in `beam6d.f90`. This caused the memory log file to grow very large, and likely had a performance cost. The call to `mod_alloc` has been removed. PR #749 (R. De Maria)
* Fixed an issue where the extraction of `Sixin.zip` would not trigger properly. The decision whether to call the extraction is taken on the existence of a specific file in the run folder. The file exist check didn't trigger properly with the new file open wrapper in `mod_units`. PR #751 (V.K. Berglyd Olsen)
* When writing `fort.4` (a duplicate of `fort.2`), both files need to be opened in module `mod_fluc`. Previously `fort.2` was opened globally during init, but is now only opened when input parsing is run. An additional open call was missing in `mod_fluc`. This made the nagfor compiler very unhappy. PR #752 (V.K. Berglyd Olsen)
* Fixed a little bug in aperture which becomes evident if `tlen` starts to deviate from `dcum(iu)`. All the calculations in aperture (especially for the backtracking algorithm) are based on `dcum`, but if the code has to deal with the beginning/end of the ring, it was using `tlen`, whereas, to be consistent, it should use `dcum(iu)`. PR #755 (A. Mereghetti)
* Fixed a bug where lines after the `LOAD` keyword in the `LIMI` block were ignored. PR #754
* Fixed a potential bug with fixed file units set in the `elens` module where the unit was within the range of those automatically assigned by `mod_units`. PR #759  (V.K. Berglyd Olsen)
* In a number of places in tracking, a subroutine looping over all particles was called from within a particle loop. These were moved to the outside of the loop. PR #759 (V.K. Berglyd Olsen)

**Code Improvements and Changes**

* The two modules handling file units have been merged to a single module. Each previous module was handling dynamic allocation of file unit and wrapping the open statement for compiler options separately. The new module now also writes a log file listing all files opened and closed through the module. PR #747 (V.K. Berglyd Olsen)
* The parsing of the aperture limitations file in the `LIMI` block has been updated and made more consistent with how other secondary input files are handled elsewhere in SixTrack. This simplifies the code. PR #754 (A. Mereghetti)
* No-longer used particle array size variables in the collimation code were cleaned up. PR #665 (K. Sjobak)
* The cavities now call the `part_updatePartEnergy` routine after changing particle energy. The subroutine now also optionally computes the change in particle angle. PR #759 (V.K. Berglyd Olsen)

### Version 5.1.1 [13.12.2018] - BOINC Release

**BOINC Specific Changes**

* Updated to use BOINC lib and API 7.14.2.

**Known Issues**

* BOINC with API does not run properly when built on latest Ubuntu LTS and Debian when both the BOINC API and SixTrack is built with the Gnu compiler. Mixing gcc and ifort or nagfor runs fine. This does not seem to be an issue when building on Fedora and CentOS. The Linux executables provided for this release are built on CentOS 7.

**Test Suite**

* A test can now be configured to automaticall stop on certain turn numbers using the `CRKILLSWITCH` flag in the `SETTINGS` block. This will help ensuring that tests actually restart from checkpoint data.
* The test suite can now verify that specified tests actually do restart when building with checkpoint/restart support.

### Version 5.1.0 [11.12.2018] - Release

While this release includes regular bug fixes and changes, the primary focus is on making code improvements that allows for a wider range of studies to be run on BOINC.

**BOINC Specific Changes**

* The elens module is now compatible with BOINC, meaning that it is properly checkpointed and the output files are properly wrapped for use with BOINC.
* The aperture module is now compatible with BOINC, meaning that it is properly checkpointed and the output files are properly wrapped for use with BOINC.

**Input File Format Changes**

* The input parser now enforces the use of a `NEXT` flag after the `PRINT` block. That is, the block must be properly closed like all other input blocks. The fact that such a flag has not been required in the past is due to a loophole in the old input parsing code. Note that usage of the `PRINT` block is deprecated. The same functionality is achieved by specifying a `PRINT` command in the `SETTINGS` block. The `PRINT` block will be removed in a future release.
* The format of aperture input information (`LIMI` block) has been changed, breaking compatibility with previous versions. In particular, tilt angle is now set as last column, as implemented by MADX-to-SixTrack converter, instead of being last but two. Moreover, the aperture offset is not expressed in terms of local value of the survey but as actual offset; hence, offset values are subtracted from particle coordinates and not summed when the aperture is checked. All users of the Fluka-SixTrack coupling will have to update their aperture model - the change in the pre-processing script to do so will come shortly.

**Other User Side Changes**

* Added a `FINALSTATE` flag in the `SETTINGS` block in `fort.3` that writes a binary or text file (via roundctl) of all particles at the end of tracking (but before post-processing). The flag takes `binary` or `text` as an option, specifying the file format. The `final_state.dat` or `final_state.bin` file produced also contains the particles flagged as lost during tracking.
* Added a `HASH` module that can be used for computing the md5sum of output files. The primary purpose of this is for checking that the output is consistent in the test suite or when results are returned from BOINC.
* Increased information in the error message produced when an error is encountered in the parsing of the `fort.3` input file.
* E-lens module can now handle a radial profile read from text file.
* E-lens kick are now fully chromatic.
* E-lens current and kinetic energy can be modified during tracking via `DYNK`.

**Build System**

* The build system now requires CMake 3.2; up from 3.0.

**Bug Fixes**

* Minor bug in parsing of the `ELENS` input block where specifying a non-existent element would trigger the wrong error message to be returned.
* Fixes a bug in checkpoint/restarting where `DUMP` restart information was not written to the secondary checkpoint file.
* Fixed a bug in checkpoint/restarting where an infinite loop might occur when both primary and secondary checkpoint files were corrupt.
* Checkpoint/restarting did not work as expected when building with the nagfor compiler. This compiler is more strict than gfortran and ifort on how files are accessed, which caused a nagfor built executable to try to overwrite `fort.10` with a dummy file even if it existed.
* Fixed a bug in beam--beam in the case of the `ibeco` flag being set to 0. In this case the beam offset would be computed with uninitialised variables.

**Code Improvements and Changes**

* All open file units registered in the module `mod_units` and `file_units` are now flushed after post-processing, and before the `HASH` and `ZIPF` modules are called.

**Test Suite**

* The `Sixin.zip` files used for testing BOINC builds have been removed from the repository. These are now generated when needed by the test suite CMake.

### Version 5.0.3 [22.11.2018] - Release

**Bug Fixes**

* Due to inconsistent if-statements, FMA would sometimes write to a non-existent file unit without opening the file properly.
* DYNK would not build with nagfor due to two malformed strings.
* Fixed building with `ROUND_ZERO` flag. The `matlib_bouncer` was trying to call non-existent round-towards-zero routines for `exp_mb()` functions (which is handled by round down `exp_rd()`), and `log10_mb` which was missing for `ROUND_ZERO` and has been added to `CRLIBM`.
* A variable for DYNK FIR/IIR functions was written to memory twice, but with an invalid index for one of those writes. This write has been removed.
* Removed a remaining `COLLIMAT` flag in the beam-gas module that prevented the `BEAMGAS` flag from building.
* Fixed an infinite loop bug in Checkpoint/Restart where the restart would repeatedly try to open `fort.96` if the reading failed while reading DUMP checkpoint data. This would only occur if both `fort.95` and `fort.96` files were corrupt.
* Fixed a bug where ifort would not accept `180_fPrec` as a valid floating point number in the `aperture` module.
* Removed references to intrinsic `ieee_arithmetic` due to significant performance loss when building with gfortran.
* The particle array `dpsv1` was wrongly updated in various parts of SixTrack. This has now been corrected. The error only affected ion tracking, but not protons.
* Previous fix to trombone elements in v5.0.1 was only for thin tracking. The same fix has now been applied to thick tracking.
* Loading beam population with the `DIST` block did not update all energy/momentum arrays correctly on initialisation. This is now fixed. The change is very small.

**User Side Changes**

* Electron lenses have been re-enabled, and the code reverted/rewritten to produce the same results as for SixTrack 4.7.18. Further updates will be made.
* The FMA output files will no longer contain NaN values. The NaN values were deliberate return values for `atan2(0,0)` call. This behaviour has been changed to return `0d0` instead, following [IEEE Std 1003.1-2017 (Revision of IEEE Std 1003.1-2008)](http://pubs.opengroup.org/onlinepubs/9699919799/). The nagfor compiler does not comply with this standard, and the return value for the nagfor built executables is handled by an additional if statement.
* Added DYNK support for coupled 4D beam-beam elements.
* Introduced consistent settings of coupling in the strong beam for 4D.
* A new output file, `sim_meta.dat`, has been added. The file lists name value pairs of information about the last run simulation.
* A new output file, `sim_time.dat`, has been added. The file lists time stamps throughout the simulation in key points, as well as compute cpu time averages.
* Parsing of beam distribution file `fort.13` has been updated. The energy per particle value read from the file is ignored and computed from the delta_p (`dpsv`) value provided.
* The `DIST` block has been updated to use the `CRLIBM` rounding library.
* The Collimation module now accepts dist format 0, which bypasses the internal beam distribution generator and instead uses the one defined by the `INIT` or `DIST` block.

**Build System**

* Added a `defaultBuild.sh` script that builds SixTrack with NAFF and libArchive support.
* The `buildLibraries.sh` script has been split up and rewritten to support more libraries. The script can be called with no arguments to build all libraries, or with `boinc`, `libarchive` or `hdf5` to build the respective libraries only. Dependencies are handled automatically.
* Building with flags `G4COLLIMAT`, `BEAMGAS` and `MERLINSCATTER` now disables building of the Differential Algebra executable (SixDA).
* The first 7 characters of the git hash is now added to the executable name as part of the version number.
* Added symlink name `sixtrack` pointing to the SixTrack executable in the build directory.
* Changed the way lquadmath is linked on Mac to make it more robust.

**Other Changes**

* A set of developer tools for MAD-X/SixTrack output testing and comparison has been added in the `devtools` folder.
* The standard output from SixTrack has been cleaned up and tweaked a little.

**Code Improvements and Changes**

* A new subroutine for initialising the random number generator has been added. This routine has proper boundary checks on the seeds. It can also optionally accept one seed. instead of two. The second seed is then calculated by a fixed offset.
* The internal NAFF library has been replaced by a rewritten external library now included as a submodule.
* Remaining code for multiple machines (different random seeds) has bean cleaned out. The feature was already disabled.
* Aperture checks have been inlined to improve performance.
* Particle vectors for [x,y] and [xp,yp] has been split up into 4 separate vectors.

### Version 5.0.2 [23.08.2018] - Release

**Bug Fixes**

* Fixed bug where the Solenoid element would use the energy of the previous turn for computations.

### Version 5.0.1 [22.08.2018] - Release

**User Side Changes**

* DUMP now uses the roundctl library when writing float values (except for the s coordinate). This fixes the failing test(s) where some values are printed as -0.0E0. This ensures that the output is consitent across compilers and platforms.
* Added DUMP format 101 for as a debugging dump format.
* MULTIPOLES now consider the curvature effect when there is a quadrupolar field in the dipole (element 11). This is mainly useful to model combined function magnets.
* Added new "debugging" format to DUMP, format 101.
* Updated DUMP manual and headers for consistency.
* COLLIMAT compilation flag removed. The collimation code is always compiled, and is controlled by the `COLL` block in `fort.3`.
* Fixes to phase trombone implementation.
* FMA is now using `atan2` in stead of `atan` to get the phase.
* ABS function added to DYNK.
* Flag `ibidu` and file `fort.32` was removed.
* Improved error reporting and consistency throughout the code.
* SixTrack output header updated, adding including GIT hash and key build flags. Timing report at the end of execution also improved.

**Documentation**

* Collimation was added to manual.
* HDF5 module was documented.

**Code Improvements and Changes**

* The old fixed length getfields_split routine has been removed entirely from the source. It is replaced by `chr_split` from the `string_tools` module.
* All remaining usage of zero chars have has been removed from the source. Zero-terminated strings should now be handled by interface routines for c. The fixed length `stringzerotrim` routine has been removed, but the variable length chr_trimZero remains in case it is needed in the future.
* Collimation arrays only allocated when COLL block is present. When present, memory use is significantly increased (as before).
* Particles and energy update codes centralized in new module `mod_particles`.
* Some NaNs are now quiet.
* SixTrack no longer exits on underflow when compiled with build type `debug`.
* Libraries boinc, libarchive, and zlib were moved from `source` to `lib`.
* Consolidated `alloc`/`resize` calls and removed redundant interface `resize` from `mod_alloc`.
* Updated `DYNK`, `FMA`, `ZIPF`, and `ROOT` blocks to use the new string split routines. 
* New flexible wrapping system `chr_fromReal()` for converting floating point numbers to strings. Added to `DUMP`, making the DUMP outputs consistent between compilers.

**Tests**

* Added tests for nearly all of the DYNK FUN statements. Only PELP, PIPE, and FIR are currently untested.
* The old SixTest test harness was finally deleted, and the `test` folder was cleaned up.
* Updated canonicals for `prob1` / `prob3`, `fma*`, `dump*`. All tests except `*elens*` should now work.
* Added test that the user manual builds correctly.
* Continuous integration (CI) was added.

**Bug Fixes**

* MULT block with non-existent single element now causes an error instead of writing the coefficients to invalid memory.
* Fixed write-to-outside-an-array error in `track_thin`.

### Version 5.0 RC3 [12.07.2018] - Release Candidate

**User Side Changes**

* All flags related to scaling of particle and element arrays have been removed. All these arrays now scale dynamically. The related compiler flags have been removed. that is: BIGNPART, HUGENPART, BIGNBLZ, HUGENBLZ.
* Major updates to the build system. Both versions of SixTrack (DA and non-DA) executables are now built, and the test folders have been merged.
* Major changes to collimation and aperture losses. Fixes a number of issues. Results will be different due to improved precision.
* Completely rewritten input parsing. General changes are:
  * Allows for blank lines and in-line comments starting with `!`
  * Allows for indented block content. Block name and NEXT still has to start at the beginning of the line.
  * Full support for single and double quoted strings.
  * Note that the parser is now stricter on blocks being properly closed and not repeated (unless repetition of blocks is specifically suported, see the manual.)
* Solenoid code has been fixed.
* Fixes to HDF5 and ROOT file output.
* New SETTINGS block where the PRINT flag has been moved, and where a QUIET flag has been added. The QUIET flag takes an integer number from 0 to 3 that reduces the amount of output generated. A global DEBUG flag has also been added that will cause te input parser to echo back parsed values.
* Major cleanup of runtime output in general.

**Code Improvements and Changes**

* The custom astuce pre-processor has been removed and everything is now handled by the C pre-processor. This implies that:
  * Includes are now handled by `#include` statements and the files are located in the include folder. `+cd` blocks are gone.
  * Nearly all common blocks have been removed and the variables moved into modules.
* Renamed the folders `SixTrack` to `source` and `SixTest` and `SixTest_da` to `test`.
* Major cleanup of the source and test folders, and the build target directory is now outside the source folder.
* Major restructuring of the source code. Mainly the splitting up of old `sixtrack.s` file into sensible files and converted to modules when reasonable.
* Started using mod_alloc also for internal subroutine/module arrays. This helps tracking memory usage.
* Added a completely new string split routine, and started phasing out both the old ones. The new routine uses dynamic array allocation and therefore have no strict limits to size and number of elements.
* Completely rewritten subroutine daten. All blocks are now parsed line by line in a single fort.2/fort.3 loop. Lines are either parsed in the module corresponding to the given block, or in a separate subroutine in the new module sixtrack_input. Diagnostics and debug routines have been added.
* Opening of standard SixTrack I/O files are now handled by module mod_units. This is different from the dynamically assigned file units module as these are fixed units with compiler flags for FIO, CRLIBM and BOINC handled internally. This removes the need to specifiy open calls for each set of flags every time a file is opened.
* Wrapper for converting strings to real(32/64/128), presenting one user interface `chr_cast()` to the programmer, internally using different codes for CRLIBM, non-CRLIBM, and FIO. This cleanes up input reading a lot.

**Disabled or Removed Features**

* Tune-shift corrections have been removed. This also removes the need for NAGLIB.
* Electron lens is currently disabled, but is planned re-enabled before final release.
* BNLELENS and RHICELENS flags removed.

### Version 5.0 RC2 [11.06.2018] - Release Candidate

**Bug Fixes and Updates**

* Fixed checkpoint/restarting.
* Updated BOINC.
* Time compilation flag was removed.

### Version 5.0 RC1 [06.04.2018] - Release Candidate

**User Side Changes**

* DUMP now accepts -1 as file unit input, in which case a unit will be assigned dynamically.
* Added optional QUIET flag that will stop SixTrack from reporting initial and final values of particle pairs.
* Fixed the previously broken DA version of the code. Still may have issues with naglib when compiled with gfortran.
* Added HDF5 block for alternative ouptut for DUMP and Scatter.
* BDEX and FLUKA was merged.

**Code Improvements and Changes**

* SixTrack is now Fortran 2008 Free Form, and depcrecated syntax has been converted or removed.
* Significant portions of the code has been split out into Fortran modules, including many common blocks.
* Abstraction of math functions: One single routine to call for each math function (sin(x) etc) which then calls the appropriate crlibm / libm method.
* SixTrack can now be compiled in single, double or quad precision.
* Partially completed:
  * Autoscaling of arrays. Will replace BIGNPART, HUGENPART, etc. flags eventually.
  * Dynamic file unit assignment.
  * ROOT/XROOTD (EOS) support. Not documented.
* Fixed:
  * Closed orbit search.
* Removed:
  * thin/thck6dua. Replacd by DYNK, see documenttaion.
  * BPM flag. Replaced by DUMP, see documentation.
* Currently broken:
  * Electron lens.
* File units can now be dynamically allocated by the file_units module.
* Added a strings tools module to do string manipulation.

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
