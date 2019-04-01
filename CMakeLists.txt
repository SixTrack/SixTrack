#
#  SixTrack CMake Build File
# ===========================
#  Making use of https://cmake.org/Wiki/CMakeForFortranExample
#  and the file created for building Merlin.
#

###################################################################################################
# General configuration
###################################################################################################

cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

# Force an out of source build to stop clutter
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(STATUS "SixTrack requires an out-of-source build.")
  message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
  message(STATUS "  CMakeCache.txt")
  message(STATUS "  CMakeFiles")
  message(STATUS "Once these files are removed, create a separate directory and run CMake from there")
  message(FATAL_ERROR "In-source build detected")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
if(NOT CMAKE_GENERATOR)
  set(CMAKE_GENERATOR "Unix Makefiles")
endif()

# Extract SixTrack version info
file(READ ${CMAKE_SOURCE_DIR}/source/version.f90 IFILE)
string(REGEX MATCH ".*version = \"([0-9.]+)\"" VERSION_TMP ${IFILE})
set(STRING_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH ".*numvers = ([0-9]+)"      VERSION_TMP ${IFILE})
set(NUM_VERSION    "${CMAKE_MATCH_1}")

# Set the project
project(SixTrack VERSION ${STRING_VERSION})
enable_language(Fortran)

# Get the current git commit hash
find_package(Git)
if(GIT_FOUND)
  execute_process(COMMAND git log -n1 --pretty=format:%H WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE GIT_REVISION)
  string(STRIP ${GIT_REVISION} GIT_REVISION)
  string(SUBSTRING ${GIT_REVISION} 0 7 GIT_BUILD)
  set(GIT_BUILD "-${GIT_BUILD}")
else()
  set(GIT_REVISION "Unknown")
  set(GIT_BUILD    "")
endif()

# Print some build info
message(STATUS "SixTrack Version: ${${PROJECT_NAME}_VERSION}")
message(STATUS "Version Major:    ${${PROJECT_NAME}_VERSION_MAJOR}")
message(STATUS "Version Minor:    ${${PROJECT_NAME}_VERSION_MINOR}")
message(STATUS "Version Patch:    ${${PROJECT_NAME}_VERSION_PATCH}")
message(STATUS "Numeric Version:  ${NUM_VERSION}")
message(STATUS "Git Revision:     ${GIT_REVISION}")
message(STATUS "System Name:      ${CMAKE_SYSTEM_NAME}")
message(STATUS "System Processor: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "CMake Generator:  ${CMAKE_GENERATOR}")
message(STATUS "CMake Source:     ${CMAKE_SOURCE_DIR}")
message(STATUS "CMake Binary:     ${CMAKE_BINARY_DIR}")
message(STATUS "Project Source:   ${PROJECT_SOURCE_DIR}")
message(STATUS "Project Binary:   ${PROJECT_BINARY_DIR}")

###################################################################################################
# SixTrack configuration options
###################################################################################################

# Default ON
option(TILT "Allow elements to be tilted (by error)" ON)
option(FAST "Which implementation of drifts to use in thin6d" ON)
option(CRLIBM "Use correctly rounded libmath instead of system libmath" ON)
option(SIXDA "Build differential algebra version (NOT dynamic aperture!)" ON)
option(STF "Single Track File: Write all tracks for postprocessing to fort.90 instead of 32 separate files fort.59 to fort.90. Required for >64 particles" ON)
option(NAFF "Link to external NAFFlib for FMA." ON)
option(HASHLIB "Build the md5 hash library" ON)

# Default OFF
OPTION(CERNLIB "Link to external CERNLIB library for PAW plots. Otherwise use internally defined dummy functions from deck plotdumy" OFF)
OPTION(API "Link to external BOINC API (if BOINC but not API, just use internal dummy library)" OFF)
OPTION(BOINC "BOINC distributed computing support" OFF)
OPTION(CR "Checkpoint / restarting" OFF)
OPTION(BEAMGAS "Beam-gas scattering" OFF)
OPTION(HDF5 "Adds support for the HDF5 block which enables writing output to a single HDF5 file." OFF)
OPTION(FIO "Use FortranIO from Fortran2003 to correctly round ASCII input/output" OFF)
OPTION(MERLINSCATTER "Interaction physics for collimation from Merlin" OFF)
OPTION(G4COLLIMAT "Interaction physics for collimation from Geant4" OFF)
OPTION(FLUKA "Couple to FLUKA for beam collimation" OFF)
OPTION(FFIELD "Fringe Field for Quadrupole" OFF)
OPTION(PYTHIA "Link with Pythia8 for SCATTER module" OFF)
OPTION(LIBARCHIVE "Link with LIBARCHIVE, needed for BOINC API and ZIPF block" OFF)
OPTION(ROOT "Link to the root libraries" OFF)

# Build Options
option(32BIT "Create a 32bit binary" OFF)
option(64BIT "Create a 64bit binary" ON)
option(STATIC "Create a statically linked binary" ON)
option(NATIVE "Enable optimisations for the current machine" OFF)
option(AVX "Enable use of the Advanced Vector Extensions (AVX) instruction set - Sandy Bridge and later" OFF)
option(AVX2 "Enable use of the Advanced Vector Extensions 2 (AVX2) instruction set - Haswell and later" OFF)
option(AVX-512 "Enable use of the Advanced Vector Extensions 512 (AVX-512) instruction set - (Currently targeting Skylake-X)" OFF)
option(LTO "Enable link time optimisations" OFF)

# Float Precision (Default: 64BITM)
option(32BITM "Floats are 32bit (single precision)" OFF)
option(64BITM "Floats are 64bit (double precision)" ON)
option(128BITM "Floats are 128bit (quad precision)" OFF)

# Rounding (Default: ROUND_NEAR)
option(ROUND_NEAR "Always round to the nearest number after floating point operations" ON)
option(ROUND_UP "Always round up after floating point operations" OFF)
option(ROUND_DOWN "Always round down after floating point operations" OFF)
option(ROUND_ZERO "Always round towards zero after floating point operations" OFF)

# Testing and Debug Options
option(BUILD_TESTING "Enable the test suite" OFF)
option(COVERAGE "Enable build flags for testing code coverage with gcov (only works with GNU compilers)" OFF)
option(GPROF "Enable build flags for code profiling with gprof (only works with GNU compilers)" OFF)
option(WARN "Enable build flags for additional compiler warnings" OFF)
option(DEBUG "Debugging output" OFF)
option(NOSYMLINK "Don't use symlinks for test inputs and canonicals" OFF)

set(BOINC_DIR ${CMAKE_SOURCE_DIR}/lib/boinc CACHE PATH "Path to the boinc build folder")
set(LIBARCHIVE_BUILD_DIR ${CMAKE_SOURCE_DIR}/lib/libarchive_build CACHE PATH "Path to the libarchive build folder")
set(LIBARCHIVE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/lib/libarchive CACHE PATH "Path to the libarchive source folder")
set(PYTHIA_DIR ${CMAKE_SOURCE_DIR}/lib/pythia CACHE PATH "Path to the Pythia build folder")

###################################################################################################
# Sixtrack binary name construction
###################################################################################################

# FFLAGS depend on the compiler
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set(SIXTRACK_FEAT_CR "${CMAKE_SYSTEM_NAME}_${Fortran_COMPILER_NAME}")
set(SIXTRACK_FEAT_DA "${CMAKE_SYSTEM_NAME}_${Fortran_COMPILER_NAME}")

if(NOT TILT)
  set(SIXTRACK_FEAT_CR "notilt_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "notilt_${SIXTRACK_FEAT_DA}")
endif()
if(NOT FAST)
  set(SIXTRACK_FEAT_CR "nofast_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "nofast_${SIXTRACK_FEAT_DA}")
endif()
if(NOT STF)
  set(SIXTRACK_FEAT_CR "nostf_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "nostf_${SIXTRACK_FEAT_DA}")
endif()
if(CRLIBM AND ROUND_NEAR)
  set(SIXTRACK_FEAT_CR "crlibm_rn_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "crlibm_rn_${SIXTRACK_FEAT_DA}")
endif(CRLIBM AND ROUND_NEAR)
if(CRLIBM AND ROUND_UP)
  set(SIXTRACK_FEAT_CR "crlibm_ru_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "crlibm_ru_${SIXTRACK_FEAT_DA}")
endif(CRLIBM AND ROUND_UP)
if(CRLIBM AND ROUND_DOWN)
  set(SIXTRACK_FEAT_CR "crlibm_rd_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "crlibm_rd_${SIXTRACK_FEAT_DA}")
endif(CRLIBM AND ROUND_DOWN)
if(CRLIBM AND ROUND_ZERO)
  set(SIXTRACK_FEAT_CR "crlibm_rz_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "crlibm_rz_${SIXTRACK_FEAT_DA}")
endif(CRLIBM AND ROUND_ZERO)
if(CERNLIB)
  set(SIXTRACK_FEAT_CR "cernlib_${SIXTRACK_FEAT_CR}")
endif()
if(API)
  set(SIXTRACK_FEAT_CR "api_${SIXTRACK_FEAT_CR}")
endif()
if(BOINC)
  SET(SIXTRACK_FEAT_CR "boinc_${SIXTRACK_FEAT_CR}")
endif()
if(CR)
  set(SIXTRACK_FEAT_CR "cr_${SIXTRACK_FEAT_CR}")
endif()
if(BEAMGAS)
  set(SIXTRACK_FEAT_CR "beamgas_${SIXTRACK_FEAT_CR}")
endif()
if(DEBUG)
  set(SIXTRACK_FEAT_CR "debug_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "debug_${SIXTRACK_FEAT_DA}")
endif()
if(HDF5)
  set(SIXTRACK_FEAT_CR "hdf5_${SIXTRACK_FEAT_CR}")
endif()
if(FIO)
  set(SIXTRACK_FEAT_CR "fio_${SIXTRACK_FEAT_CR}")
  set(SIXTRACK_FEAT_DA "fio_${SIXTRACK_FEAT_DA}")
endif()
if(MERLINSCATTER)
  set(SIXTRACK_FEAT_CR "merlinscatter_${SIXTRACK_FEAT_CR}")
endif()
if(G4COLLIMAT)
  set(SIXTRACK_FEAT_CR "g4collimat_${SIXTRACK_FEAT_CR}")
endif()
if(FLUKA)
  set(SIXTRACK_FEAT_CR "FLUKA_${SIXTRACK_FEAT_CR}")
endif(FLUKA)
if(FFIELD)
	SET(SIXTRACK_FEAT_CR "FFIELD_${SIXTRACK_FEAT_CR}")
endif(FFIELD)
if(PYTHIA)
  set(SIXTRACK_FEAT_CR "pythia_${SIXTRACK_FEAT_CR}")
endif()
if(LIBARCHIVE)
  set(SIXTRACK_FEAT_CR "libarchive_${SIXTRACK_FEAT_CR}")
endif()

if(ROOT)
  set(SIXTRACK_FEAT_CR "root_${SIXTRACK_FEAT_CR}")
endif(ROOT)

if(STATIC)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_static")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_static")
endif()
if(AVX)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_avx")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_avx")
endif(AVX)
if(AVX2)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_avx2")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_avx2")
endif(AVX2)
if(AVX-512)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_avx512")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_avx512")
endif(AVX-512)
if(NATIVE)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_native")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_native")
endif()

set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_${CMAKE_SYSTEM_PROCESSOR}")
set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_${CMAKE_SYSTEM_PROCESSOR}")

if(32BIT)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_32bit")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_32bit")
endif()
if(64BIT)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_64bit")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_64bit")
endif()

if(32BITM)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_single")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_single")
endif()
if(64BITM)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_double")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_double")
endif()
if(128BITM)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_quad")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_quad")
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_DEBUG")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_DEBUG")
endif()

# Create the sixtrack binary name string
set(SIXTRACK_BINARY_NAME    "${PROJECT_NAME}_${NUM_VERSION}${GIT_BUILD}_${SIXTRACK_FEAT_CR}")
set(SIXTRACK_BINARY_NAME_DA "SixDA_${NUM_VERSION}${GIT_BUILD}_${SIXTRACK_FEAT_DA}")

###################################################################################################
# Setting the c pre-processor definitons
###################################################################################################

list(APPEND PREPRO_FLAGS_TOOLS TILT FAST BOINC BEAMGAS STF MERLINSCATTER G4COLLIMAT FLUKA CR PYTHIA NAFF FFIELD)
list(APPEND PREPRO_FLAGS_LIB CRLIBM ROUND_NEAR ROUND_UP ROUND_DOWN ROUND_ZERO CERNLIB DEBUG HDF5 FIO LIBARCHIVE ROOT HASHLIB)
list(APPEND PREPRO_FLAGS_COMPILE GFORTRAN IFORT NAGFOR WIN32)

###################################################################################################
# Select which fortran files to generate and compile
###################################################################################################

list(APPEND FORT90_LIB
  sixtrack version core_tools constants sixtrack_input mod_units
  common_modules close strings string_tools mainarrays mod_meta
  mod_time mod_alloc zipf matrix_inv lielib beam6d
  bouncy_castle plato_seq bdex fma dabnew end_sixtrack
  ranecu wire hions elens mod_dist scatter dump dynk
  aperture mod_fluc beam6d_fox six_fox postprocessing
  read_write ranlux collimation mod_particles cheby utils
)

# Files unique to DA or CR version
list(APPEND FORT90_DA main_da sixda)
list(APPEND FORT90_CR main_cr sixve track_thin track_thick)

# List of files that must be treated by DAFOR
list(APPEND DAFOR_F90 beam6d_fox six_fox sixda)

if(CR)
  list(APPEND FORT90_LIB checkpoint_restart)
endif(CR)

# if(DEBUG)
#   list(APPEND FORT90_LIB debug_dump)
# endif(DEBUG)

if(HDF5)
  list(APPEND FORT90_LIB hdf5_output hdf5_sixtrack)
endif(HDF5)

if(HASHLIB)
  list(APPEND FORT90_LIB mod_hash)
endif(HASHLIB)

if(BEAMGAS)
  LIST(APPEND FORT90_CR beamgas)
endif(BEAMGAS)

if(FLUKA)
  list(APPEND FORT90_LIB mod_fluka fluka)
endif(FLUKA)

if(FFIELD)
	list(APPEND FORT90_LIB  mod_ffield mod_ffield_Table_n_Track)
endif(FFIELD)


if(ROOT)
  LIST(APPEND FORT90_LIB root_output)
endif(ROOT)

if(BOINC AND NOT API)
  # We have a BOINC build but no API, load the fake API substitute for debugging.
  list(APPEND FORT90_LIB boinc_dummy)
endif(BOINC AND NOT API)

if(PYTHIA)
  add_library(pythia_wrapper STATIC pythia_wrapper.h pythia_wrapper.cpp)
  target_compile_options(pythia_wrapper PRIVATE)
  list(APPEND FORT90_LIB mod_pythia)
endif(PYTHIA)

foreach(loop IN LISTS FORT90_LIB)
  list(APPEND SIXTRACK_SOURCES ${loop}.f90)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${loop}.f90)
endforeach()
foreach(loop IN LISTS FORT90_CR)
  list(APPEND SIXTRACK_SOURCES_CR ${loop}.f90)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${loop}.f90)
endforeach()
foreach(loop IN LISTS FORT90_DA)
  list(APPEND SIXTRACK_SOURCES_DA ${loop}.f90)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${loop}.f90)
endforeach()
foreach(loop IN LISTS DAFOR_F90)
  set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${loop}a.f90)
endforeach()

if(LIBARCHIVE)
  add_library(libArchive_wrapper STATIC libArchive_Fwrapper.c libArchive_wrapper.c libArchive_wrapper.h)
  target_compile_options(libArchive_wrapper PRIVATE -std=c99)
endif(LIBARCHIVE)

if(LIBARCHIVE OR HDF5)
  add_library(z STATIC IMPORTED)
  set_target_properties(z PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/zlib/install/lib/libz.a)
endif(LIBARCHIVE OR HDF5)

###################################################################################################
# Checking if some requested packages exist
###################################################################################################

if(G4COLLIMAT)
  if(STATIC)
    find_package(Geant4 REQUIRED static)
  else()
    find_package(Geant4 REQUIRED)
  endif(STATIC)

  #no error checking yet
  message(STATUS "Geant4 include folders: ${Geant4_INCLUDE_DIRS}")
  message(STATUS "Geant4 use file: ${Geant4_USE_FILE}")
  message(STATUS "Geant4 libraries: ${Geant4_LIBRARIES}")
  message(STATUS "Geant4 definitions: ${Geant4_DEFINITIONS}")
# include(${Geant4_USE_FILE})
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/g4collimat)
endif(G4COLLIMAT)

##Possible check for Root
if(ROOT)
  #/usr/share/root/cmake/FindROOT.cmake
  SET(CMAKE_MODULE_PATH /usr/share/root/cmake/ ${CMAKE_MODULE_PATH})
  if(STATIC)
    find_package(ROOT REQUIRED static)
  else()
    find_package(ROOT REQUIRED)
  endif(STATIC)
  if(NOT ROOT_FOUND)
    MESSAGE(FATAL_ERROR "Root enabled build requested but root was not found!")
  endif(NOT ROOT_FOUND)
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/root_output)
endif(ROOT)

if(BUILD_TESTING OR BOINC OR LIBARCHIVE OR G4COLLIMAT OR CERNLIB)

  find_package(Threads REQUIRED)
  if(NOT Threads_FOUND)
    message(FATAL_ERROR "Could not find a thread library! This is required for your selected build configuration.")
  endif(NOT Threads_FOUND)

endif(BUILD_TESTING OR BOINC OR LIBARCHIVE OR G4COLLIMAT OR CERNLIB)

###################################################################################################
# Checking if requested options are compatible
###################################################################################################

#Check we do not have both 32 and 64 bit options enabled
if(32BIT)
  if(64BIT)
    message(FATAL_ERROR "Cannot build with both 64 and 32 bit options enabled. Please pick one only.")
  endif()
endif()

#Check we have at least one 32 or 64 bit option enabled
if(NOT 32BIT)
  if(NOT 64BIT)
    message(FATAL_ERROR "Must pick either a 64 or 32 bit build. Currently neither are selected")
  endif()
endif()

if(HDF5)
  if(CR)
    message(FATAL_ERROR "Cannot build HDF5 with CR enabled")
  endif()
endif()
if(BOINC)
  if(NOT CRLIBM)
    message(FATAL_ERROR "Building with BOINC requires CRLIBM enabled")
  endif(NOT CRLIBM)
  if(NOT CR)
    message(FATAL_ERROR "Building with BOINC requires CR enabled")
  endif()
  if(NOT LIBARCHIVE)
    message(FATAL_ERROR "Building with BOINC requires LIBARCHIVE enabled")
  endif()
  if(SIXDA)
    message(STATUS "Will not build SixDA when BOINC is enabled.")
    set(SIXDA OFF)
  endif()
endif(BOINC)
if(API)
  if(NOT BOINC)
    message(FATAL_ERROR "Building with API requires BOINC enabled")
  endif(NOT BOINC)
endif(API)
if(FLUKA)
  if(SIXDA)
    message(STATUS "Will not build SixDA when FLUKA is enabled.")
    set(SIXDA OFF)
  endif()
endif()
if(G4COLLIMAT)
  if(SIXDA)
    message(STATUS "Will not build SixDA when G4COLLIMAT is enabled.")
    set(SIXDA OFF)
  endif()
endif()
if(MERLINSCATTER)
  if(SIXDA)
    message(STATUS "Will not build SixDA when MERLINSCATTER is enabled.")
    set(SIXDA OFF)
  endif()
endif()
if(BEAMGAS)
  if(SIXDA)
    message(STATUS "Will not build SixDA when BEAMGAS is enabled.")
    set(SIXDA OFF)
  endif()
endif()

if(64BIT)
  if(CERNLIB)
    message(FATAL_ERROR "Cannot build CERNLIB with 64BIT enabled.")
  endif()
endif()

if(32BITM)
  if(64BITM)
    message(FATAL_ERROR "Cannot use both single and double precision math.")
  endif(64BITM)
  if(128BITM)
    message(FATAL_ERROR "Cannot use both single and quad precision math.")
  endif(128BITM)
  if(FLUKA)
    message(FATAL_ERROR "Cannot use FLUKA with single precision math.")
  endif(FLUKA)
  if(CRLIBM)
    message(FATAL_ERROR "Cannot use CRLIBM with single precision math.")
  endif(CRLIBM)
endif(32BITM)

if(64BITM)
  if(32BITM)
    message(FATAL_ERROR "Cannot use both double and single precision math.")
  endif(32BITM)
  if(128BITM)
    message(FATAL_ERROR "Cannot use both double and quad precision math.")
  endif(128BITM)
endif(64BITM)

if(128BITM)
  if(32BITM)
    message(FATAL_ERROR "Cannot use both quad and single precision math.")
  endif(32BITM)
  if(64BITM)
    message(FATAL_ERROR "Cannot use both quad and double precision math.")
  endif(64BITM)
  if(FLUKA)
    message(FATAL_ERROR "Cannot use FLUKA with quad precision math.")
  endif(FLUKA)
  if(CRLIBM)
    message(FATAL_ERROR "Cannot use CRLIBM with quad precision math.")
  endif(CRLIBM)
endif(128BITM)

if(NOT SIXDA)
  message(STATUS "Not building SixDA version.")
endif()

###################################################################################################
# Settings for compiler build flags
###################################################################################################

if(STATIC)
  #This just sets the default build setting for shared libraries.
  SET(BUILD_SHARED_LIBS OFF)
endif()

#FLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

#lots of games to play here.... start with fortran flags
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  ##############
  #  gfortran  #
  ##############

  set(GFORTRAN ON)
  set(CMAKE_Fortran_FLAGS "-frecord-marker=4 -fno-second-underscore -funroll-loops -std=f2008")
  set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
  set(PREPRO_FLAGS -cpp -E -P)

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES aarch64 )
    set (CMAKE_Fortran_FLAGS "-ffp-contract=off ${CMAKE_Fortran_FLAGS}")
  endif()

  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -ffpe-trap=invalid,zero,overflow -fcheck=all")

  if(WARN)
    set(CMAKE_Fortran_FLAGS_RELEASE "-Wuninitialized ${CMAKE_Fortran_FLAGS_RELEASE}")
  endif()

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES ppc64le )
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fno-expensive-optimizations -fno-inline-functions")
  endif()
  if(LTO)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -flto")
  endif()

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mfpmath=sse -msse2")

    if(32BIT)
      set(CMAKE_Fortran_FLAGS "-m32 ${CMAKE_Fortran_FLAGS}")
    endif(32BIT)

    if(64BIT)
      set(CMAKE_Fortran_FLAGS "-m64 ${CMAKE_Fortran_FLAGS}")
    endif(64BIT)

    if(AVX)
      set(CMAKE_Fortran_FLAGS "-mavx ${CMAKE_Fortran_FLAGS}")
    endif(AVX)

    if(AVX2)
      set(CMAKE_Fortran_FLAGS "-mavx -mavx2 -mno-fma4 -mno-fma ${CMAKE_Fortran_FLAGS}")
    endif(AVX2)

    if(AVX-512)
      set(CMAKE_Fortran_FLAGS "-mavx -mavx2 -mavx512f -mavx512cd -mavx512vl -mavx512bw -mavx512dq -mno-fma4 -mno-fma ${CMAKE_Fortran_FLAGS}")
    endif(AVX-512)

    if(NATIVE)
      set(CMAKE_Fortran_FLAGS "-march=native -mno-fma4 -mno-fma ${CMAKE_Fortran_FLAGS}")
    endif(NATIVE)
  endif()

  if(STATIC)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      SET(CMAKE_Fortran_LINKER_FLAGS "-Bstatic -static-libgcc -static-libgfortran ${CMAKE_Fortran_LINKER_FLAGS}")
    else()
      SET(CMAKE_Fortran_LINKER_FLAGS "-static ${CMAKE_Fortran_LINKER_FLAGS}")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  endif(STATIC)

  if(COVERAGE)
    SET(CMAKE_Fortran_FLAGS "-fprofile-arcs -ftest-coverage ${CMAKE_Fortran_FLAGS}")
    SET(CMAKE_Fortran_LINKER_FLAGS "-fprofile-arcs ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(COVERAGE)

  if(GPROF)
    SET(CMAKE_Fortran_FLAGS "-pg ${CMAKE_Fortran_FLAGS}")
    SET(CMAKE_Fortran_LINKER_FLAGS "-pg ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(GPROF)

elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  ###########
  #  ifort  #
  ###########

  set(IFORT ON)
  set(CMAKE_Fortran_FLAGS_RELEASE "-fp-model source -fp-model strict -fp-model no-except -O3 -no-fma -std15 -assume realloc-lhs")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-fp-model source -fp-model strict -fp-model no-except -O0 -g -no-fma -std15 -assume realloc-lhs")
  set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
  set(PREPRO_FLAGS -I${CMAKE_SOURCE_DIR}/source/include -fpp -E -P)

  if(LTO)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -ffat-lto-objects -ipo")
  endif()

  if(WARN)
    set(CMAKE_Fortran_FLAGS "-warn unused ${CMAKE_Fortran_FLAGS}")
  endif(WARN)

  if(32BIT)
    set(CMAKE_Fortran_FLAGS_RELEASE "-m32 ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG "-m32 ${CMAKE_Fortran_FLAGS_DEBUG}")
    set(CMAKE_Fortran_LINKER_FLAGS "-m32 ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(32BIT)

  if(64BIT)
    set(CMAKE_Fortran_FLAGS_RELEASE "-m64 ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG "-m64 ${CMAKE_Fortran_FLAGS_DEBUG}")
    set(CMAKE_Fortran_LINKER_FLAGS "-m64 ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(64BIT)

  if(AVX)
    set(CMAKE_Fortran_FLAGS "-mavx ${CMAKE_Fortran_FLAGS}")
  endif(AVX)

  if(AVX2)
    set(CMAKE_Fortran_FLAGS "-mavx -mavx2 -no-fma ${CMAKE_Fortran_FLAGS}")
  endif(AVX2)

  if(AVX-512)
    set(CMAKE_Fortran_FLAGS "-axcommon-avx512 -no-fma ${CMAKE_Fortran_FLAGS}")
  endif(AVX-512)

  if(NATIVE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-xHost ${CMAKE_Fortran_FLAGS_RELEASE}")
  endif(NATIVE)

  if(STATIC)
    set(CMAKE_Fortran_LINKER_FLAGS "-static ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(STATIC)

elseif (Fortran_COMPILER_NAME MATCHES "nagfor.*")
  ############
  #  nagfor  #
  ############

  set(NAGFOR ON)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O4 -dusty -ieee=full -dcfuns")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -dusty -ieee=full -dcfuns -gline")
  set(PREPRO_EXEC  ${CMAKE_CXX_COMPILER})
  set(PREPRO_FLAGS -cpp -E -P)
  # set(PREPRO_EXEC  ${CMAKE_Fotran_COMPILER})
  # set(PREPRO_FLAGS -F)

  if(WARN)
    set(CMAKE_Fortran_FLAGS "-w=ulv -w=uparam ${CMAKE_Fortran_FLAGS}")
  else()
    set(CMAKE_Fortran_FLAGS "-w=all ${CMAKE_Fortran_FLAGS}")
  endif(WARN)

  if(32BIT)
    set(CMAKE_Fortran_FLAGS_RELEASE "-abi=32 ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG "-abi=32 ${CMAKE_Fortran_FLAGS_DEBUG}")
  endif(32BIT)

  if(64BIT)
    #See the nagfor man page for the different 64bit -abi options.
    set(CMAKE_Fortran_FLAGS_RELEASE "-abi=64t ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG "-abi=64t ${CMAKE_Fortran_FLAGS_DEBUG}")
  endif(64BIT)

  if(NATIVE)
    SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")
  endif(NATIVE)

  if(STATIC)
    SET(CMAKE_Fortran_LINKER_FLAGS "-Bstatic ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(STATIC)

  if(NOT STATIC)
    SET(CMAKE_Fortran_LINKER_FLAGS "-Bdynamic ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(NOT STATIC)

#generic catch option
else()
  message(STATUS "CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message(STATUS "Fortran compiler: " ${Fortran_COMPILER_NAME})
  message(STATUS "No optimized Fortran compiler flags are known, we just try -O2 or -O0 -g...")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 ${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g ${CMAKE_Fortran_FLAGS_DEBUG}")
  set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
  set(PREPRO_FLAGS -cpp -E -P)
endif ()

#Now for C flags (crlibm for now)
#${CMAKE_C_COMPILER_ID} : From the docs: "one of "Clang", "GNU", "Intel", or "MSVC". This works even if a compiler wrapper like ccache is used."
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES aarch64 )
    set (CMAKE_C_FLAGS "-ffp-contract=off ${CMAKE_C_FLAGS}")
  endif()
  set (CMAKE_C_FLAGS_RELEASE "-O3")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -g3")

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
    if(32BIT)
      set(CMAKE_C_FLAGS "-m32 ${CMAKE_C_FLAGS}")
    endif(32BIT)

    if(64BIT)
      set(CMAKE_C_FLAGS "-m64 ${CMAKE_C_FLAGS}")
    endif(64BIT)

    if(AVX)
      set(CMAKE_C_FLAGS "-mavx ${CMAKE_C_FLAGS}")
    endif(AVX)

    if(AVX2)
      set(CMAKE_C_FLAGS "-mavx -mavx2 -mno-fma4 -mno-fma ${CMAKE_C_FLAGS}")
    endif(AVX2)

    if(AVX-512)
      set(CMAKE_C_FLAGS "-mavx -mavx2 -mavx512f -mavx512cd -mavx512vl -mavx512bw -mavx512dq -mno-fma4 -mno-fma ${CMAKE_C_FLAGS}")
    endif(AVX-512)

    if(NATIVE)
      set(CMAKE_C_FLAGS "-march=native ${CMAKE_C_FLAGS}")
    endif(NATIVE)
  endif()

  if(STATIC)
    SET(CMAKE_C_LINKER_FLAGS "-static ${CMAKE_C_LINKER_FLAGS}")
  endif(STATIC)

  if(COVERAGE)
    SET(CMAKE_C_FLAGS "-fprofile-arcs -ftest-coverage ${CMAKE_C_FLAGS}")
    SET(CMAKE_C_LINKER_FLAGS "-fprofile-arcs ${CMAKE_C_LINKER_FLAGS}")
  endif(COVERAGE)

  if(GPROF)
    SET(CMAKE_C_FLAGS "-pg ${CMAKE_C_FLAGS}")
    SET(CMAKE_C_LINKER_FLAGS "-pg ${CMAKE_C_LINKER_FLAGS}")
  endif(GPROF)
endif()

if(${CMAKE_C_COMPILER_ID} MATCHES "Intel")

  if(32BIT)
    set(CMAKE_C_FLAGS "-m32 ${CMAKE_C_FLAGS}")
  endif(32BIT)

  if(64BIT)
    set(CMAKE_C_FLAGS "-m64 ${CMAKE_C_FLAGS}")
  endif(64BIT)

  if(AVX)
    set(CMAKE_C_FLAGS "-mavx ${CMAKE_C_FLAGS}")
  endif(AVX)

  if(AVX2)
    set(CMAKE_C_FLAGS "-mavx -mavx2 -no-fma ${CMAKE_C_FLAGS}")
  endif(AVX2)

  if(AVX-512)
    set(CMAKE_C_FLAGS "-axcommon-avx512 -no-fma ${CMAKE_C_FLAGS}")
  endif(AVX-512)

  if(NATIVE)
    set(CMAKE_C_FLAGS "-xHost ${CMAKE_C_FLAGS}")
  endif(NATIVE)

  set (CMAKE_C_FLAGS_RELEASE "-O3")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -g")

  if(STATIC)
    SET(CMAKE_C_LINKER_FLAGS "-static ${CMAKE_C_LINKER_FLAGS}")
  endif(STATIC)
endif()

#Now for CXX flags (merlinScatter for now)
#${CMAKE_CXX_COMPILER_ID} : From the docs: "one of "Clang", "GNU", "Intel", or "MSVC". This works even if a compiler wrapper like ccache is used."
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES aarch64 )
    set (CMAKE_CXX_FLAGS "-ffp-contract=off ${CMAKE_CXX_FLAGS}")
  endif()
  set (CMAKE_CXX_FLAGS_RELEASE "-O3")
  set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
    if(32BIT)
      set(CMAKE_CXX_FLAGS "-m32 ${CMAKE_CXX_FLAGS}")
    endif(32BIT)

    if(64BIT)
      set(CMAKE_CXX_FLAGS "-m64 ${CMAKE_CXX_FLAGS}")
    endif(64BIT)

    if(AVX)
      set(CMAKE_CXX_FLAGS "-mavx ${CMAKE_CXX_FLAGS}")
    endif(AVX)

    if(AVX2)
      set(CMAKE_CXX_FLAGS "-mavx -mavx2 -mno-fma4 -mno-fma ${CMAKE_CXX_FLAGS}")
    endif(AVX2)

    if(AVX-512)
      set(CMAKE_CXX_FLAGS "-mavx -mavx2 -mavx512f -mavx512cd -mavx512vl -mavx512bw -mavx512dq -mno-fma4 -mno-fma ${CMAKE_CXX_FLAGS}")
    endif(AVX-512)

    if(NATIVE)
      set(CMAKE_CXX_FLAGS "-march=native ${CMAKE_CXX_FLAGS}")
    endif(NATIVE)
  endif()

  if(STATIC)
    SET(CMAKE_CXX_LINKER_FLAGS "-static -static-libstdc++ ${CMAKE_CXX_LINKER_FLAGS}")
  endif(STATIC)

  if(COVERAGE)
    SET(CMAKE_CXX_FLAGS "-fprofile-arcs -ftest-coverage ${CMAKE_CXX_FLAGS}")
    SET(CMAKE_CXX_LINKER_FLAGS "-fprofile-arcs ${CMAKE_CXX_LINKER_FLAGS}")
  endif(COVERAGE)

  IF(GPROF)
    SET(CMAKE_CXX_FLAGS "-pg ${CMAKE_CXX_FLAGS}")
    SET(CMAKE_CXX_LINKER_FLAGS "-pg ${CMAKE_CXX_LINKER_FLAGS}")
  ENDIF(GPROF)

endif()

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")

  if(32BIT)
    set(CMAKE_CXX_FLAGS "-m32 ${CMAKE_CXX_FLAGS}")
  endif(32BIT)

  if(64BIT)
    set(CMAKE_CXX_FLAGS "-m64 ${CMAKE_CXX_FLAGS}")
  endif(64BIT)

  if(AVX)
    set(CMAKE_CXX_FLAGS "-mavx ${CMAKE_CXX_FLAGS}")
  endif(AVX)

  if(AVX2)
    set(CMAKE_CXX_FLAGS "-mavx -mavx2 -no-fma ${CMAKE_CXX_FLAGS}")
  endif(AVX2)

  if(AVX-512)
    set(CMAKE_CXX_FLAGS "-axcommon-avx512 -no-fma ${CMAKE_CXX_FLAGS}")
  endif(AVX-512)

  if(NATIVE)
    set(CMAKE_CXX_FLAGS "-march=xHost ${CMAKE_CXX_FLAGS}")
  endif(NATIVE)

  set (CMAKE_CXX_FLAGS_RELEASE "-O3")
  set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")

  #When using icc remember to link to libstdc++
  if(STATIC)
    find_library(stdc++ STATIC stdc++)
    SET(CMAKE_CXX_LINKER_FLAGS "-static ${CMAKE_CXX_LINKER_FLAGS}")
  else()
    find_library(stdc++ stdc++)
  endif(STATIC)

endif()

#check that the compiler has the features we need!
# Disabling this test until its cmake version is in sync with the main build script as this breaks on Ubuntu 16.04
# See https://github.com/SixTrack/SixTrack/issues/594
#include(${CMAKE_SOURCE_DIR}/fortran_support.cmake)

###################################################################################################
# Settings for what exactly to build
###################################################################################################
# Add the source files that will be used to build SixTrack
add_library(SixTrackCore STATIC ${SIXTRACK_SOURCES})
set_target_properties(SixTrackCore PROPERTIES OUTPUT_NAME "SixTrackCore")
set_target_properties(SixTrackCore PROPERTIES PREFIX "")
add_executable(SixTrackCR ${SIXTRACK_SOURCES_CR})
target_link_libraries(SixTrackCR SixTrackCore)
set_target_properties(SixTrackCR PROPERTIES OUTPUT_NAME "${SIXTRACK_BINARY_NAME}" LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}" LINKER_LANGUAGE Fortran)
if(SIXDA)
  add_executable(SixTrackDA ${SIXTRACK_SOURCES_DA})
  target_link_libraries(SixTrackDA SixTrackCore)
  set_target_properties(SixTrackDA PROPERTIES OUTPUT_NAME "${SIXTRACK_BINARY_NAME_DA}" LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}" LINKER_LANGUAGE Fortran)
endif()

# Pre-processor flags
foreach(loop IN LISTS PREPRO_FLAGS_TOOLS PREPRO_FLAGS_LIB PREPRO_FLAGS_COMPILE)
  if(${loop})
    message(STATUS "Adding preprocessor definition -D${loop} to target SixTrack")
    list(APPEND PREPRO_FLAGS -D${loop})
  endif()
endforeach()
if(32BITM)
  message(STATUS "Adding preprocessor macro -DSIXTRACK_REAL=real32 to target SixTrack")
  list(APPEND PREPRO_FLAGS -DSIXTRACK_REAL=real32 -DSINGLE_MATH)
endif()
if(64BITM)
  message(STATUS "Adding preprocessor macro -DSIXTRACK_REAL=real64 to target SixTrack")
  list(APPEND PREPRO_FLAGS -DSIXTRACK_REAL=real64 -DDOUBLE_MATH)
endif()
if(128BITM)
  message(STATUS "Adding preprocessor macro -DSIXTRACK_REAL=real128 to target SixTrack")
  list(APPEND PREPRO_FLAGS -DSIXTRACK_REAL=real128 -DQUAD_MATH)
endif()

  list(APPEND PREPRO_FLAGS -DSIXTRACK_GIT_REVISION="'${GIT_REVISION}'")

#build dafor
add_executable(dafor ${CMAKE_SOURCE_DIR}/source/tools/dafor/dafor.f90)
#we do not need any fancy flags for dafor but nagfor does not like the following
if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  SET_TARGET_PROPERTIES(dafor PROPERTIES COMPILE_FLAGS "-mno-avx -mno-avx2 -mno-avx512f -mno-avx512cd -mno-avx512vl -mno-avx512bw -mno-avx512dq")
  SET_TARGET_PROPERTIES(dafor PROPERTIES LINK_FLAGS "-mno-avx -mno-avx2 -mno-avx512f -mno-avx512cd -mno-avx512vl -mno-avx512bw -mno-avx512dq")
  #target_compile_definitions(dafor PUBLIC "-mno-avx -mno-avx2 -mno-avx512f -mno-avx512cd -mno-avx512vl -mno-avx512bw -mno-avx512dq")
endif(Fortran_COMPILER_NAME MATCHES "gfortran.*")

if(CRLIBM)
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/crlibm)
endif()
if(CRLIBM OR BUILD_TESTING)
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/roundctl) #Also needed by read90
endif()

if(HDF5)
  find_library(dl STATIC dl) #HDF5 has a capability to dynamically load filters
  add_library(hdf5_local_f STATIC IMPORTED)
  add_library(hdf5_local_s STATIC IMPORTED)
  add_library(hdf5_local_c STATIC IMPORTED)
  set_target_properties(hdf5_local_f PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/hdf5/build/bin/libhdf5_fortran-static.a)
  set_target_properties(hdf5_local_s PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/hdf5/build/bin/libhdf5_f90cstub-static.a)
  set_target_properties(hdf5_local_c PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/hdf5/build/bin/libhdf5-static.a)
  target_link_libraries(SixTrackCore hdf5_local_f hdf5_local_s hdf5_local_c z dl)
  target_include_directories(SixTrackCore PRIVATE ${CMAKE_SOURCE_DIR}/lib/hdf5/build/bin/static)
  target_include_directories(SixTrackCore PRIVATE ${CMAKE_SOURCE_DIR}/lib/hdf5/build/fortran)
endif(HDF5)

if(HASHLIB)
  set(MD5_SOURCES
    ${CMAKE_SOURCE_DIR}/source/md5/md5.c
    ${CMAKE_SOURCE_DIR}/source/md5/md5interface.c
  )
  add_library(md5_lib STATIC  ${MD5_SOURCES})
  target_include_directories(md5_lib PRIVATE ${CMAKE_SOURCE_DIR}/source/md5/)
  set_property(TARGET md5_lib PROPERTY C_STANDARD 99)
  if(32BIT)
    set_target_properties(md5_lib PROPERTIES LINK_FLAGS "-m32")
  endif(32BIT)
  if(64BIT)
    set_target_properties(md5_lib PROPERTIES LINK_FLAGS "-m64")
  endif(64BIT)
  target_link_libraries(SixTrackCore md5_lib)
endif(HASHLIB)

if(MERLINSCATTER)
  file(GLOB MERLINSCATTER_SOURCES ${CMAKE_SOURCE_DIR}/source/merlinscatter/*.cpp ${CMAKE_SOURCE_DIR}/source/merlinscatter/RandomNG/*.cpp)
  add_library(merlinscatter STATIC ${MERLINSCATTER_SOURCES})

  #Enable c++11
  set_property(TARGET merlinscatter PROPERTY CXX_STANDARD 11)

  if(32BIT)
    SET_TARGET_PROPERTIES(merlinscatter PROPERTIES LINK_FLAGS "-m32")
  endif(32BIT)
  if(64BIT)
    SET_TARGET_PROPERTIES(merlinscatter PROPERTIES LINK_FLAGS "-m64")
  endif(64BIT)
endif()

# Build NAFFlib
if(NAFF)
  if(NOT EXISTS ${CMAKE_SOURCE_DIR}/lib/NAFFlib/NAFFlib/interface/naff_sixtrack_interface.c)
    message(FATAL_ERROR "The git submodules NAFFlib is not available. Please run:
    git submodule update --init lib/NAFFlib")
  endif()

  file(GLOB NAFF_SOURCES ${CMAKE_SOURCE_DIR}/lib/NAFFlib/NAFFlib/source/*.c ${CMAKE_SOURCE_DIR}/lib/NAFFlib/NAFFlib/interface/naff_sixtrack_interface.c)
  add_library(naff STATIC ${NAFF_SOURCES})
  target_include_directories(naff PRIVATE ${CMAKE_SOURCE_DIR}/lib/NAFFlib/NAFFlib/include)
  target_compile_options(naff PRIVATE -DCOMPILE_WITHOUT_PYTHON)
  if(CRLIBM)
    target_include_directories(naff PRIVATE ${CMAKE_SOURCE_DIR}/source/crlibm)
    target_compile_options(naff PRIVATE -DCOMPILE_WITH_CRLIBM)
    target_link_libraries(naff crlibm)
  endif()

  # Enable c99
  set_property(TARGET naff PROPERTY C_STANDARD 99)

  if(32BIT)
    SET_TARGET_PROPERTIES(naff PROPERTIES LINK_FLAGS "-m32")
  endif(32BIT)
  if(64BIT)
    SET_TARGET_PROPERTIES(naff PROPERTIES LINK_FLAGS "-m64")
  endif(64BIT)
endif(NAFF)

if(BOINC AND API)
  find_library(stdc++ STATIC stdc++)
  include_directories(${BOINC_DIR}/api ${BOINC_DIR}/lib ${BOINC_DIR})

  add_library(boinc       STATIC IMPORTED)
  add_library(boinc_api   STATIC IMPORTED)
  add_library(boinc_api_f STATIC IMPORTED)

  set_target_properties(boinc       PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/lib/libboinc.a)
  set_target_properties(boinc_api   PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/api/libboinc_api.a)
  set_target_properties(boinc_api_f PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/api/boinc_api_fortran.o)

  target_link_libraries(SixTrackCR boinc_api_f boinc_api boinc stdc++)
  # target_link_libraries(SixTrackCR Threads::Threads)

  # if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  #   find_library(nsl STATIC nsl)
  #   target_link_libraries(SixTrackCR nsl)
  # endif()

  # Fix OSX linking to the correct c++ stdlib
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    find_library(c++ STATIC c++)
    target_link_libraries(SixTrackCore c++)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif()

if(FLUKA)
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/FlukaIO)
  target_link_libraries(SixTrackCore flukaIO)
endif(FLUKA)

if(PYTHIA)

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pythia_wrapper.h POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/source/pythia/pythia_wrapper.h ${CMAKE_BINARY_DIR}/pythia_wrapper.h
    DEPENDS copy_input_files
  )
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pythia_wrapper.cpp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/source/pythia/pythia_wrapper.cpp ${CMAKE_BINARY_DIR}/pythia_wrapper.cpp
    DEPENDS copy_input_files
  )

  include_directories(${PYTHIA_DIR}/include)
  include_directories(${CMAKE_CURRENT_BINARY_DIR})

  target_link_libraries(SixTrackCore pythia_wrapper)
  add_library(pythia STATIC IMPORTED)
  set_target_properties(pythia PROPERTIES IMPORTED_LOCATION ${PYTHIA_DIR}/lib/libpythia8.a)
  target_link_libraries(SixTrackCore pythia dl)

endif(PYTHIA)

add_custom_target(copy_input_files)
if(LIBARCHIVE)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libArchive_Fwrapper.c POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/source/libarchive_wr/libArchive_Fwrapper.c ${CMAKE_BINARY_DIR}/libArchive_Fwrapper.c
    DEPENDS copy_input_files
  )
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libArchive_wrapper.c POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/source/libarchive_wr/libArchive_wrapper.c ${CMAKE_BINARY_DIR}/libArchive_wrapper.c
    DEPENDS copy_input_files
  )
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libArchive_wrapper.h POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/source/libarchive_wr/libArchive_wrapper.h ${CMAKE_BINARY_DIR}/libArchive_wrapper.h
    DEPENDS copy_input_files
  )

  include_directories(${LIBARCHIVE_SOURCE_DIR}/libarchive)
  include_directories(${CMAKE_CURRENT_BINARY_DIR})

  #A very small but convenient test code for debugging libArchive.
  add_executable(libArchiveTester ${CMAKE_SOURCE_DIR}/source/libarchive_wr/libArchiveTester.cpp)
  target_link_libraries(libArchiveTester libArchive_wrapper)

  #See also further down for SixTestWrapper linking
  target_link_libraries(SixTrackCore libArchive_wrapper)
  add_library(libarchive STATIC IMPORTED)
  if(WIN32)
    set_target_properties(libarchive PROPERTIES IMPORTED_LOCATION ${LIBARCHIVE_BUILD_DIR}/libarchive/libarchive_static.a)
    add_definitions( -DLIBARCHIVE_STATIC )
  else()
    set_target_properties(libarchive PROPERTIES IMPORTED_LOCATION ${LIBARCHIVE_BUILD_DIR}/libarchive/libarchive.a)
  endif()

  target_link_libraries(SixTrackCore libarchive z Threads::Threads)
  target_link_libraries(libArchiveTester libarchive z Threads::Threads)

  add_definitions( -DLIBARCHIVE )
endif()

#We need to build dafor before we build SixTrack, but crlibm and so on can be built at the same time.
add_dependencies(SixTrackCore copy_input_files generate_input)
add_dependencies(SixTrackCR SixTrackCore)
if(SIXDA)
  add_dependencies(SixTrackDA SixTrackCore)
endif()

###################################################################################################
# Running Pre-Processors
###################################################################################################

add_custom_target(generate_input)
add_dependencies(generate_input dafor)

## Run the preprocessor on all .f90 files and output them to the build directory
foreach(loop IN ITEMS ${FORT90_LIB} ${FORT90_CR} ${FORT90_DA})
  list(FIND DAFOR_F90 ${loop} DAFOR_IDX)
  if(${DAFOR_IDX} GREATER -1)
    message(STATUS "Will run DAFOR to generate ${loop}.f90")
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f90
      COMMAND ${PREPRO_EXEC} ${PREPRO_FLAGS} ${CMAKE_SOURCE_DIR}/source/${loop}.f90 > ${CMAKE_CURRENT_BINARY_DIR}/${loop}a.f90
      COMMAND $<TARGET_FILE:dafor> ${loop}a.f90 ${loop}.f90 > ${loop}_dafor.log
      DEPENDS generate_input ${CMAKE_SOURCE_DIR}/source/${loop}.f90
    )
  else()
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f90
      COMMAND ${PREPRO_EXEC} ${PREPRO_FLAGS} ${CMAKE_SOURCE_DIR}/source/${loop}.f90 > ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f90
      DEPENDS generate_input ${CMAKE_SOURCE_DIR}/source/${loop}.f90
    )
  endif()
endforeach()

###################################################################################################
#linking
###################################################################################################
set(MSYS_QUADMATH "/usr/lib/gcc/x86_64-pc-msys/6.3.0/libquadmath.a")
#just a work around to build on MSYS2 for now
if(MSYS AND STATIC)
  target_link_libraries(dafor ${MSYS_QUADMATH})
  target_link_libraries(SixTrackCore ${MSYS_QUADMATH})
endif()

#same with OSX, using gfortran from homebrew
if(STATIC)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    get_filename_component(COMPILERDIR ${CMAKE_Fortran_COMPILER} DIRECTORY)
    string(SUBSTRING ${CMAKE_Fortran_COMPILER_VERSION} 0 1 COMPILER_MAJOR_VERSION)
    # search relative to the compiler executable
    find_library(lquadmath NAMES libquadmath.a PATHS "${COMPILERDIR}/../lib/gcc${COMPILER_MAJOR_VERSION}/" "${COMPILERDIR}/../lib/gcc/${COMPILER_MAJOR_VERSION}/")
    if (lquadmath)
      message(STATUS "Looking for lquadmath... - found")
    else()
      message(FATAL_ERROR "lquadmath not found!")
    endif()
    target_link_libraries(SixTrackCore ${lquadmath})
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(STATIC)

if(CERNLIB)
  # Use precompiled CERNLIB from lib32 folder
  add_library(g2c STATIC IMPORTED)
  add_library(graflib STATIC IMPORTED)
  add_library(grafX11 STATIC IMPORTED)
  add_library(kernlib_noshift STATIC IMPORTED)
  add_library(packlib_noshift STATIC IMPORTED)
  add_library(X11 STATIC IMPORTED)
  add_library(Xau STATIC IMPORTED)
  add_library(xcb STATIC IMPORTED)
  add_library(Xdmcp STATIC IMPORTED)

  set_target_properties(g2c PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libg2c.a)
  set_target_properties(graflib PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libgraflib.a)
  set_target_properties(grafX11 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libgrafX11.a)
  set_target_properties(kernlib_noshift PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libkernlib_noshift.a)
  set_target_properties(packlib_noshift PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libpacklib_noshift.a)
  set_target_properties(X11 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libX11.a)
  set_target_properties(Xau PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libXau.a)
  set_target_properties(xcb PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libxcb.a)
  set_target_properties(Xdmcp PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libXdmcp.a)

  find_library(dl STATIC dl)
  find_library(gcc_eh STATIC gcc_eh)

  target_link_libraries(SixTrackCore graflib grafX11 packlib_noshift kernlib_noshift g2c X11 xcb Xau Xdmcp dl Threads::Threads gcc_eh)
endif(CERNLIB)

#if we are using crlibm, link to the crlibm library
if(CRLIBM)
  target_link_libraries(SixTrackCore crlibm roundctl)
endif(CRLIBM)

#if we are using merlinscatter, link to the merlinscatter library
if(MERLINSCATTER)
  target_link_libraries(SixTrackCR merlinscatter)
endif(MERLINSCATTER)

# if we are using naff, link to the naff library
if(NAFF)
  target_link_libraries(SixTrackCore naff)
endif(NAFF)

#if we are using root, link to all the root libraries
if(ROOT)
  target_link_libraries(SixTrackCore root_output)
endif(ROOT)

if(HDF5)
  target_link_libraries(SixTrackCore ${HDF5_LIBRARIES})
endif(HDF5)

if(HASHLIB)
  target_link_libraries(SixTrackCore md5_lib)
endif(HASHLIB)

if(G4COLLIMAT)
  target_link_libraries(SixTrackCR g4collimat ${Geant4_LIBRARIES} Threads::Threads)
endif(G4COLLIMAT)

if(COVERAGE)
  target_compile_options(SixTrackCore PRIVATE -fprofile-arcs -ftest-coverage)
  target_compile_options(SixTrackCR   PRIVATE -fprofile-arcs -ftest-coverage)
  set_target_properties(SixTrackCore  PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  set_target_properties(SixTrackCR    PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  if(SIXDA)
    target_compile_options(SixTrackDA PRIVATE -fprofile-arcs -ftest-coverage)
    set_target_properties(SixTrackDA  PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  endif()
endif(COVERAGE)

set(BUILDNAME "${SIXTRACK_BINARY_NAME}" CACHE STRING "Name of build on the dashboard")
mark_as_advanced(BUILDNAME)

# For testing
IF(BUILD_TESTING)

  enable_testing()
  include(CTest)

  add_subdirectory(${CMAKE_SOURCE_DIR}/test)
  configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake  ${CMAKE_BINARY_DIR})

  #need to build the testing tools
  add_executable(verify10 ${CMAKE_SOURCE_DIR}/test/bin/verify10.f90)
  add_executable(checkf10 ${CMAKE_SOURCE_DIR}/test/bin/checkf10.f90)
  add_executable(compf10 ${CMAKE_SOURCE_DIR}/test/bin/compf10.f90)
  add_executable(checkf110 ${CMAKE_SOURCE_DIR}/test/bin/checkf110.f90)
  add_executable(checkf1014 ${CMAKE_SOURCE_DIR}/test/bin/checkf1014.f90)
  add_executable(read90 ${CMAKE_SOURCE_DIR}/test/read90/read90.f90)
  add_executable(readDump3 ${CMAKE_SOURCE_DIR}/test/readDump3/readDump3.f90)

  target_link_libraries(read90 roundctl)

  set_target_properties(verify10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  set_target_properties(checkf10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  set_target_properties(compf10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  set_target_properties(checkf110 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  set_target_properties(checkf1014 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  set_target_properties(read90 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
  set_target_properties(readDump3 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")

  #The testing wrapper
  add_executable(SixTestWrapper ${CMAKE_SOURCE_DIR}/test/bin/SixTestWrapper.cpp)
  if (LIBARCHIVE)
    target_link_libraries(SixTestWrapper libArchive_wrapper libarchive z Threads::Threads)
  endif()
  set_target_properties(SixTestWrapper PROPERTIES LINK_FLAGS "${CMAKE_CXX_LINKER_FLAGS}")

  target_link_libraries(SixTestWrapper Threads::Threads)

  if(COVERAGE)
    target_compile_options(SixTestWrapper PRIVATE -fprofile-arcs -ftest-coverage)
    SET_TARGET_PROPERTIES(SixTestWrapper PROPERTIES LINK_FLAGS "-fprofile-arcs ${CMAKE_CXX_LINKER_FLAGS}")
  endif(COVERAGE)

  if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    target_link_libraries(SixTestWrapper stdc++)
  endif()
  set_property(TARGET SixTestWrapper PROPERTY CXX_STANDARD 11)

  #just a work around to build on MSYS2 for now
  if(MSYS AND STATIC)
    #set(CMAKE_EXECUTABLE_SUFFIX_Fortran "")
    target_link_libraries(verify10 ${MSYS_QUADMATH})
    target_link_libraries(checkf10 ${MSYS_QUADMATH})
    target_link_libraries(checkf110 ${MSYS_QUADMATH})
    target_link_libraries(checkf1014 ${MSYS_QUADMATH})
    target_link_libraries(compf10 ${MSYS_QUADMATH})
    target_link_libraries(read90 ${MSYS_QUADMATH})
    target_link_libraries(SixTestWrapper ${MSYS_QUADMATH})
    target_link_libraries(readDump3 ${MSYS_QUADMATH})
  endif()
endif()

#To make packages, run "make package" in the build directory
#set(CPACK_GENERATOR RPM DEB STGZ TGZ TBZ2 TXZ ZIP 7Z)
set(CPACK_GENERATOR RPM DEB TGZ ZIP) #Take only the most common targets
set(CPACK_PACKAGE_CONTACT "sixtrack.contact@cern.ch")
set(CPACK_PACKAGE_VERSION ${SixTrack_VERSION})
set(CPACK_PACKAGE_VENDOR "CERN")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SixTrack symplectic particle tracking code for particle accelerators.")
set(CPACK_PACKAGE_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${${PROJECT_NAME}_VERSION_PATCH})
include(CPack)

install(TARGETS SixTrackCR RUNTIME DESTINATION bin )

if(NOT WIN32)
  add_custom_command(TARGET SixTrackCR POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${SIXTRACK_BINARY_NAME}" "sixtrack")
endif()
