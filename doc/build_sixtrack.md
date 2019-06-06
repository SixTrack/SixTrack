# Building SixTrack

## Dependencies

SixTrack is built with CMake. SixTrack is written in Fortran 2008, and can be built using either gfortran, ifort or nagfor. A number of wrappers and dependencies are written in c and c++, and can be built with gcc. Building with icc and clang generally works, but the test suite wrapper does not currently build with icpc.

### Versions

The following versions are known minimum restrictions, or compiler versions that are known to work. Certain versions have known bugs, like gfortran 7.3, however these are not critical.

* CMake: Requires at least 3.0.
* GNU Compiler: Requires at least 5.x. Recommended 8.0 or higher.
* Intel Compiler: Tested with 18.0 or higher.
* Nag Compiler: Tested with 6.1.

## Build Options

A build script `cmake_six` is provided to help build SixTrack. Each feature below can be enabled by specifying the flag, or disabled by specifying the flag with a minus before it.

Example:

    ./cmake_six gfortran release CR -STF

Builds SixTrack release with gfortran, with checkpoint/restaring support and without singletrackfile (see below). Compiler and build type are optional arguments.

### Build Types

* **release**:       Builds the standard version of SixTrack
* **debug**:         Builds a debug version of SixTrack with the `-g` flag and no compiler optimisation. It is significantly slower, and only to be used for code debugging.

### Options Enabled by Default

* **TILT**:          Allow elements to be tilted (by error
* **FAST**:          Which implementation of drifts to use in thin6d
* **CRLIBM**:        Use correctly rounded libmath instead of system libmath
* **SIXDA**:         Build differential algebra version (NOT dynamic aperture!)
* **STF**:           Single Track File. Write all tracks for postprocessing to singletrackfile.dat instead of 32 separate files (fort.59 - fort.90). This option is required for more than 64 particles.

### Options Disabled by Default

* **BUILD_TESTING**: Enable the test suite.
* **BOINC**:         Builds BOINC version of SixTrack.
* **CR**             Enables checkpoint/restart support. Required for BOINC.
* **LIBARCHIVE**:    Link with LIBARCHIVE, and is required for BOINC and ZIPF block if ZIPF is off.
* **ZLIB**:          Link with zlib and minizip, and is required for BOINC and the ZIPF block if LIBARCHIVE is off.
* **BEAMGAS**:       Beam-gas scattering.
* **FIO**:           Use FortranIO from Fortran2003 to correctly round ASCII input/output. this option overrides CRLIBM when reading/writing.
* **CERNLIB**:       Link to external CERNLIB library for PAW plots. Otherwise use internally defined dummy functions.

### Linking of External Tools (Disabled by Default)

* **HDF5**:          Adds support for the HDF5 block which enables writing output to a single HDF5 file.
* **ROOT**:          Adds support for writing to ROOT files. Experimental and undocumented.
* **MERLINSCATTER**: Interaction physics for collimation from Merlin.
* **G4COLLIMAT**:    Interaction physics for collimation from Geant4.
* **FLUKA**:         Couple to FLUKA for beam collimation.


### Binary Type

* **32BIT**:         Create a 32bit binary. Default OFF.
* **64BIT**:         Create a 64bit binary. Default ON.
* **STATIC**:        Create a statically linked binary. Default ON.
* **NATIVE**:        Enable optimizations for the current machine. Default OFF.
* **AVX**:           Enable use of the Advanced Vector Extensions (AVX) instruction set. Sandy Bridge and later. Default OFF.
* **AVX2**:          Enable use of the Advanced Vector Extensions 2 (AVX2) instruction set. Haswell and later. Default OFF:
* **AVX-512**:       Enable use of the Advanced Vector Extensions 512 (AVX-512) instruction set. Currently targeting Skylake-X. Default OFF.

### Numerical Precision (Mutually Exclusive)

* **32BITM**:        Floats are 32bit (single precision. Default OFF.
* **64BITM**:        Floats are 64bit (double precision. Default ON.
* **128BITM**        Floats are 128bit (quad precision). Default OFF.

### Rounding (Mutually Exclusive)

* **ROUND_NEAR**:    Always round to the nearest number after floating point operations. Default ON.
* **ROUND_UP**:      Always round up after floating point operations. Default OFF.
* **ROUND_DOWN**:    Always round down after floating point operations. Default OFF.
* **ROUND_ZERO**:    Always round towards zero after floating point operations. Default OFF.

### Debugging Options

* **COVERAGE**:      Enable build flags for testing code coverage with gcov. Only works with GNU compilers. Default OFF.
* **GPROF**:         Enable build flags for code profiling with gprof. Only works with GNU compilers. Default OFF.
* **WARN**:          Enable build flags for additional compiler warnings. Default OFF.
