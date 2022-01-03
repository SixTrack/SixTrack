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
option(CRLIBM "Use correctly rounded libmath instead of system libmath" ON)
option(SIXDA "Build differential algebra version (NOT dynamic aperture!)" ON)
option(NAFF "Link to external NAFFlib for FMA" ON)
option(DISTLIB "Link to external DISTlib for the beam distribution" OFF)
option(HASHLIB "Build the md5 hash library" ON)
option(ZLIB "Build and link zlib. Needed for BOINC and ZIPF block" ON)

# Default OFF
option(CERNLIB "Link to external CERNLIB library for PAW plots. Otherwise use internally defined dummy functions from deck plotdumy" OFF)
option(BOINC "BOINC distributed computing support" OFF)
option(CR "Checkpoint/restarting" OFF)
option(BEAMGAS "Beam-gas scattering" OFF)
option(HDF5 "Adds support for the HDF5 block which enables writing output to a single HDF5 file." OFF)
option(FIO "Use FortranIO from Fortran2003 to correctly round ASCII input/output" OFF)
option(MERLINSCATTER "Interaction physics for collimation from Merlin" OFF)
option(G4COLLIMATION "Interaction physics for collimation from Geant4" OFF)
option(FLUKA "Couple to FLUKA for beam collimation" OFF)
option(PYTHIA "Link with Pythia8 for SCATTER module" OFF)
option(ROOT "Link to the root libraries" OFF)

# Build Options
option(32BIT "Create a 32bit binary" OFF)
option(64BIT "Create a 64bit binary" ON)
option(STATIC "Create a statically linked binary" ON)
option(NATIVE "Enable optimisations for the current machine" OFF)
option(AVX "Enable use of the Advanced Vector Extensions (AVX) instruction set - Sandy Bridge and later" OFF)
option(AVX2 "Enable use of the Advanced Vector Extensions 2 (AVX2) instruction set - Haswell and later" OFF)
option(AVX512 "Enable use of the Advanced Vector Extensions 512 (AVX512) instruction set" OFF)
option(AVX512_SKYLAKE "Enable use of the Advanced Vector Extensions 512 (AVX512) instruction set targeting Skylake" OFF)
option(LTO "Enable link time optimisations" OFF)
OPTION(PREFETCH "Enables the -fprefetch-loop-arrays compiler flag for gcc" OFF)

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
option(INITNAN "Enable build flags for initialising reals to signaling NaN" OFF)
option(DEBUG "Debugging output" OFF)
option(NOSYMLINK "Don't use symlinks for test inputs and canonicals" OFF)
option(MEMUSAGE "Writes PID and actual memory usage to sim_meta.dat (GNU/UNIX only)" OFF)

set(BOINC_DIR ${CMAKE_SOURCE_DIR}/lib/boinc CACHE PATH "Path to the boinc build folder")
set(PYTHIA_DIR ${CMAKE_SOURCE_DIR}/lib/pythia CACHE PATH "Path to the Pythia build folder")

###################################################################################################
# Sixtrack binary name construction
###################################################################################################

# FFLAGS depend on the compiler
get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# Remove .exe from compiler name on Windows
string(REPLACE ".exe" "" Fortran_COMPILER_NAME ${Fortran_COMPILER_NAME})

set(SIXTRACK_FEAT_CR "${CMAKE_SYSTEM_NAME}_${Fortran_COMPILER_NAME}")
set(SIXTRACK_FEAT_DA "${CMAKE_SYSTEM_NAME}_${Fortran_COMPILER_NAME}")

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
if(BOINC)
  set(SIXTRACK_FEAT_CR "boinc_${SIXTRACK_FEAT_CR}")
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
if(G4COLLIMATION)
  set(SIXTRACK_FEAT_CR "g4collimation_${SIXTRACK_FEAT_CR}")
endif(G4COLLIMATION)
if(FLUKA)
  set(SIXTRACK_FEAT_CR "FLUKA_${SIXTRACK_FEAT_CR}")
endif(FLUKA)
if(PYTHIA)
  set(SIXTRACK_FEAT_CR "pythia_${SIXTRACK_FEAT_CR}")
endif()
if(ZLIB)
  set(SIXTRACK_FEAT_CR "zlib_${SIXTRACK_FEAT_CR}")
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
if(AVX512)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_avx512")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_avx512")
endif(AVX512)
if(AVX512_SKYLAKE)
  set(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_avx512skylake")
  set(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_avx512skylake")
endif(AVX512_SKYLAKE)
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

list(APPEND PREPRO_FLAGS_TOOLS BOINC BEAMGAS MERLINSCATTER G4COLLIMATION FLUKA CR PYTHIA NAFF)
list(APPEND PREPRO_FLAGS_LIB CRLIBM ROUND_NEAR ROUND_UP ROUND_DOWN ROUND_ZERO CERNLIB HDF5 FIO ROOT HASHLIB ZLIB DISTLIB)
list(APPEND PREPRO_FLAGS_COMPILE GFORTRAN IFORT NAGFOR WIN32 DEBUG MEMUSAGE)

###################################################################################################
# Select which fortran files to generate and compile
###################################################################################################

list(APPEND FORT90_LIB
  aperture
  bdex
  beam6d
  beam6d_fox
  bouncy_castle
  cheby
  close
  coll_common
  coll_crystal
  coll_db
  coll_dist
  coll_jawfit
  coll_k2
  collimation
  common_modules
  constants
  core_tools
  dabnew
  dump
  dynk
  elements
  elens
  end_sixtrack
  fma
  lielib
  linear_optics
  mainarrays
  matrix_inv
  mod_alloc
  mod_dist
  mod_ffield
  mod_ffield_aux
  mod_find
  mod_fluc
  mod_funlux
  mod_geometry
  mod_meta
  mod_particles
  mod_random
  mod_ranecu
  mod_ranlux
  mod_time
  mod_units
  mod_utils
  pdgid
  plato_seq
  postprocessing
  read_write
  scatter
  six_fox
  sixtrack
  sixtrack_input
  string_tools
  strings
  version
  wire
  zipf
)

# Files unique to DA or CR version
list(APPEND FORT90_DA main_da sixda)
list(APPEND FORT90_CR main_cr sixve tracking track_thin track_thick)

# List of files that must be treated by DAFOR
list(APPEND DAFOR_F90 beam6d_fox cheby elens six_fox sixda)

if(CR)
  list(APPEND FORT90_LIB checkpoint_restart)
endif(CR)

if(HDF5)
  list(APPEND FORT90_LIB hdf5_output hdf5_sixtrack)
endif(HDF5)

if(HASHLIB)
  list(APPEND FORT90_LIB mod_hash)
endif(HASHLIB)

if(BEAMGAS)
  list(APPEND FORT90_CR beamgas)
endif(BEAMGAS)

if(FLUKA)
  list(APPEND FORT90_LIB mod_fluka)
endif(FLUKA)

if(ROOT)
  list(APPEND FORT90_LIB root_output)
endif(ROOT)

if(G4COLLIMATION)
  list(APPEND FORT90_LIB geant4)
endif(G4COLLIMATION)

if(BOINC)
  list(APPEND FORT90_LIB mod_boinc)
endif(BOINC)

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

if(ZLIB)
  set(ZLIB_DIR ${CMAKE_SOURCE_DIR}/lib/zlib)
  include_directories(${ZLIB_DIR} ${CMAKE_CURRENT_BINARY_DIR}/zlib)
  add_subdirectory(${ZLIB_DIR} ${CMAKE_CURRENT_BINARY_DIR}/zlib)
endif(ZLIB)

###################################################################################################
# Checking if some requested packages exist
###################################################################################################

##Possible check for Root
if(ROOT)
  set(CMAKE_MODULE_PATH /usr/share/root/cmake/ ${CMAKE_MODULE_PATH})
  if(STATIC)
    find_package(ROOT REQUIRED static)
  else()
    find_package(ROOT REQUIRED)
  endif(STATIC)
  if(NOT ROOT_FOUND)
    message(FATAL_ERROR "Root enabled build requested but root was not found!")
  endif(NOT ROOT_FOUND)
  message(STATUS "ROOT include folders: ${ROOT_INCLUDE_DIRS}")
  message(STATUS "ROOT definitions:${ROOT_DEFINITIONS}")
  message(STATUS "ROOT CXX flags:  ${ROOT_CXX_FLAGS}")
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/root_output)
endif(ROOT)

if(G4COLLIMATION)
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
  message(STATUS "Geant4 CXX flags: ${Geant4_CXX_FLAGS}")
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/g4collimation)
endif(G4COLLIMATION)

if(BUILD_TESTING OR BOINC OR G4COLLIMATION OR CERNLIB)

  find_package(Threads REQUIRED)
  if(NOT Threads_FOUND)
    message(FATAL_ERROR "Could not find a thread library! This is required for your selected build configuration.")
  endif(NOT Threads_FOUND)

endif(BUILD_TESTING OR BOINC OR G4COLLIMATION OR CERNLIB)

###################################################################################################
# Checking if requested options are compatible
###################################################################################################

# Check we do not have both 32 and 64 bit options enabled
if(32BIT)
  if(64BIT)
    message(FATAL_ERROR "Cannot build with both 64 and 32 bit options enabled. Please pick one only.")
  endif()
endif()

# Check we have at least one 32 or 64 bit option enabled
if(NOT 32BIT)
  if(NOT 64BIT)
    message(FATAL_ERROR "Must pick either a 64 or 32 bit build. Currently neither are selected.")
  endif()
endif()

if(HDF5)
  if(CR)
    message(FATAL_ERROR "Cannot build HDF5 with CR enabled.")
  endif()
  if(NOT ZLIB)
    message(FATAL_ERROR "Building with HDF5 requires ZLIB enabled.")
  endif()
endif()
if(SIXDA)
  if(CR)
    message(STATUS "Will not build SixDA when checkpoint/restart is enabled.")
    set(SIXDA OFF)
  endif()
endif()
if(BOINC)
  if(NOT CRLIBM)
    message(FATAL_ERROR "Building with BOINC requires CRLIBM enabled.")
  endif(NOT CRLIBM)
  if(NOT CR)
    message(FATAL_ERROR "Building with BOINC requires CR enabled.")
  endif()
  if(NOT ZLIB)
    message(FATAL_ERROR "Building with BOINC requires ZLIB enabled.")
  endif()
  if(NOT HASHLIB)
    message(FATAL_ERROR "Building with BOINC requires HASHLIB enabled.")
  endif()
endif(BOINC)
if(FLUKA)
  if(SIXDA)
    message(STATUS "Will not build SixDA when FLUKA is enabled.")
    set(SIXDA OFF)
  endif()
endif()
if(G4COLLIMATION)
  if(SIXDA)
    message(STATUS "Will not build SixDA when G4COLLIMATION is enabled.")
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
  if(DISTLIB)
    message(FATAL_ERROR "Cannot use DISTLIB with single precision math.")
  endif(DISTLIB)
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
  if(DISTLIB)
    message(FATAL_ERROR "Cannot use DISTLIB with quad precision math.")
  endif(DISTLIB)
endif(128BITM)

if(NOT SIXDA)
  message(STATUS "Not building SixDA version.")
endif(NOT SIXDA)

# Since each AVX option includes the lower ones, just turn the lower flags off to avoid duplication
if(AVX512_SKYLAKE)
  set(AVX OFF)
  set(AVX2 OFF)
  set(AVX512 OFF)
endif(AVX512_SKYLAKE)
if(AVX512)
  set(AVX OFF)
  set(AVX2 OFF)
endif(AVX512)
if(AVX2)
  set(AVX OFF)
endif(AVX2)

###################################################################################################
# Settings for compiler build flags
###################################################################################################

if(STATIC)
  # This just sets the default build setting for shared libraries.
  set(BUILD_SHARED_LIBS OFF)
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(MIN_OSX_VERSION 10.10.5)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# FLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
  ##############
  #  gfortran  #
  ##############

  set(GFORTRAN ON)
  set(CMAKE_Fortran_FLAGS "-frecord-marker=4 -fno-second-underscore -funroll-loops -std=f2008")
  set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
  set(PREPRO_FLAGS -cpp -E -P)

  #get the compiler version
  string(REGEX MATCH "^([0-9]*)" VERSION_TMP ${CMAKE_Fortran_COMPILER_VERSION})
  set(GFORTRAN_VERSION "${CMAKE_MATCH_1}")
  #fix for building with gfortran >= 10
  if(${GFORTRAN_VERSION} GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif(${GFORTRAN_VERSION} GREATER_EQUAL 10)

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES aarch64 )
    set (CMAKE_Fortran_FLAGS "-ffp-contract=off ${CMAKE_Fortran_FLAGS}")
  endif()

  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -fcheck=all")

  if(WARN)
    set(CMAKE_Fortran_FLAGS_RELEASE "-Wuninitialized ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-Wuninitialized -ffpe-trap=invalid,zero,overflow,underflow ${CMAKE_Fortran_FLAGS_DEBUG}")
  else()
    set(CMAKE_Fortran_FLAGS_DEBUG   "-ffpe-trap=invalid,zero,overflow ${CMAKE_Fortran_FLAGS_DEBUG}")
  endif()

  if(INITNAN)
    set(CMAKE_Fortran_FLAGS_RELEASE "-finit-real=snan -finit-character=0 ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-finit-real=snan -finit-character=0 ${CMAKE_Fortran_FLAGS_DEBUG}")
  endif()

  if(MEMUSAGE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-fall-intrinsics ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-fall-intrinsics ${CMAKE_Fortran_FLAGS_DEBUG}")
  endif()

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES ppc64le )
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fno-expensive-optimizations -fno-inline-functions")
  endif()
  if(LTO)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -flto -Wl,-flto")
  endif(LTO)

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)
    set(CMAKE_Fortran_FLAGS "-mmacosx-version-min=${MIN_OSX_VERSION} ${CMAKE_Fortran_FLAGS}")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mfpmath=sse -msse2")

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
      set(CMAKE_Fortran_FLAGS "-mavx -mavx2 ${CMAKE_Fortran_FLAGS}")
    endif(AVX2)

    if(AVX512)
      set(CMAKE_Fortran_FLAGS "-mavx -mavx2 -mavx512f ${CMAKE_Fortran_FLAGS}")
    endif(AVX512)

    if(AVX512_SKYLAKE)
      set(CMAKE_Fortran_FLAGS "-mavx -mavx2 -mavx512f -mavx512cd -mavx512bw -mavx512vl ${CMAKE_Fortran_FLAGS}")
    endif(AVX512_SKYLAKE)

    if(NATIVE)
      set(CMAKE_Fortran_FLAGS "-march=native ${CMAKE_Fortran_FLAGS}")
    endif(NATIVE)

    if(CRLIBM)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mno-fma4 -mno-fma")
    endif(CRLIBM)

    if(PREFETCH)
      set(CMAKE_Fortran_FLAGS "-fprefetch-loop-arrays ${CMAKE_Fortran_FLAGS}")
    endif(PREFETCH)
  endif()

  if(STATIC)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set(CMAKE_Fortran_LINKER_FLAGS "-Bstatic -static-libgcc -static-libgfortran ${CMAKE_Fortran_LINKER_FLAGS}")
    else()
      set(CMAKE_Fortran_LINKER_FLAGS "-static ${CMAKE_Fortran_LINKER_FLAGS}")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  endif(STATIC)

  if(COVERAGE)
    set(CMAKE_Fortran_FLAGS "-fprofile-arcs -ftest-coverage ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_LINKER_FLAGS "-fprofile-arcs ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(COVERAGE)

  if(GPROF)
    set(CMAKE_Fortran_FLAGS "-pg ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_LINKER_FLAGS "-pg ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(GPROF)

elseif(Fortran_COMPILER_NAME MATCHES "ifort.*")
  ###########
  #  ifort  #
  ###########

  set(IFORT ON)
  set(CMAKE_Fortran_FLAGS_RELEASE "-fp-model source -fp-model strict -fp-model no-except -O3 -no-fma -std18 -assume realloc-lhs")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-fp-model source -fp-model strict -fp-model no-except -O0 -g -traceback -no-fma -std18 -assume realloc-lhs")
  set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
  set(PREPRO_FLAGS -I${CMAKE_SOURCE_DIR}/source/include -fpp -E -P)

  if(LTO)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -ffat-lto-objects -ipo")
  endif(LTO)

  if(WARN)
    set(CMAKE_Fortran_FLAGS "-warn unused -check bounds ${CMAKE_Fortran_FLAGS}")
  endif(WARN)

  if(INITNAN)
    set(CMAKE_Fortran_FLAGS_RELEASE "-init=snan,arrays ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-init=snan,arrays ${CMAKE_Fortran_FLAGS_DEBUG}")
  endif()

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

  if(AVX512)
    set(CMAKE_Fortran_FLAGS "-axcommon-avx512 -no-fma ${CMAKE_Fortran_FLAGS}")
  endif(AVX512)

  if(AVX512_SKYLAKE)
    set(CMAKE_Fortran_FLAGS "-axskylake-avx512 -no-fma ${CMAKE_Fortran_FLAGS}")
  endif(AVX512_SKYLAKE)

  if(NATIVE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-xHost ${CMAKE_Fortran_FLAGS_RELEASE}")
  endif(NATIVE)

  if(STATIC)
    set(CMAKE_Fortran_LINKER_FLAGS "-static ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(STATIC)

elseif(Fortran_COMPILER_NAME MATCHES "nagfor.*")
  ############
  #  nagfor  #
  ############

  set(NAGFOR ON)
  set(CMAKE_Fortran_FLAGS_RELEASE "-O4 -dusty -ieee=full -dcfuns")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -dusty -ieee=full -dcfuns -gline")
  set(PREPRO_EXEC  ${CMAKE_CXX_COMPILER})
  set(PREPRO_FLAGS -cpp -E -P)

  if(WARN)
    set(CMAKE_Fortran_FLAGS "-w=ulv -w=uparam ${CMAKE_Fortran_FLAGS}")
  else()
    set(CMAKE_Fortran_FLAGS "-w=all ${CMAKE_Fortran_FLAGS}")
  endif(WARN)

  if(INITNAN)
    set(CMAKE_Fortran_FLAGS_RELEASE "-nan ${CMAKE_Fortran_FLAGS_RELEASE}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-nan ${CMAKE_Fortran_FLAGS_DEBUG}")
  endif()

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
    set(CMAKE_Fortran_LINKER_FLAGS "-Bstatic ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(STATIC)

  if(NOT STATIC)
    set(CMAKE_Fortran_LINKER_FLAGS "-Bdynamic ${CMAKE_Fortran_LINKER_FLAGS}")
  endif(NOT STATIC)

else()
  # Generic catch option
  message(STATUS "CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message(STATUS "Fortran compiler: " ${Fortran_COMPILER_NAME})
  message(STATUS "No optimized Fortran compiler flags are known, we just try -O2 or -O0 -g...")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 ${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g ${CMAKE_Fortran_FLAGS_DEBUG}")
  set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
  set(PREPRO_FLAGS -cpp -E -P)
endif()

# C Flags
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES aarch64 )
    set (CMAKE_C_FLAGS "-ffp-contract=off ${CMAKE_C_FLAGS}")
  endif()
  set (CMAKE_C_FLAGS_RELEASE "-O3")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -g3")

  if(LTO)
    set(CMAKE_Fortran_C_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto -Wl,-flto")
  endif(LTO)

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
      set(CMAKE_C_FLAGS "-mavx -mavx2 ${CMAKE_C_FLAGS}")
    endif(AVX2)

    if(AVX512)
      set(CMAKE_C_FLAGS "-mavx -mavx2 -mavx512f ${CMAKE_C_FLAGS}")
    endif(AVX512)

    if(AVX512_SKYLAKE)
      set(CMAKE_C_FLAGS "-mavx -mavx2 -mavx512f -mavx512cd -mavx512bw -mavx512vl ${CMAKE_C_FLAGS}")
    endif(AVX512_SKYLAKE)

    if(NATIVE)
      set(CMAKE_C_FLAGS "-march=native ${CMAKE_C_FLAGS}")
    endif(NATIVE)

    if(CRLIBM)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mno-fma4 -mno-fma")
    endif(CRLIBM)

    if(PREFETCH)
      set(CMAKE_C_FLAGS "-fprefetch-loop-arrays ${CMAKE_C_FLAGS}")
    endif(PREFETCH)
  endif()

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)
    set(CMAKE_C_FLAGS "-mmacosx-version-min=${MIN_OSX_VERSION} ${CMAKE_C_FLAGS}")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)

  if(STATIC)
    set(CMAKE_C_LINKER_FLAGS "-static ${CMAKE_C_LINKER_FLAGS}")
  endif(STATIC)

  if(COVERAGE)
    set(CMAKE_C_FLAGS "-fprofile-arcs -ftest-coverage ${CMAKE_C_FLAGS}")
    set(CMAKE_C_LINKER_FLAGS "-fprofile-arcs ${CMAKE_C_LINKER_FLAGS}")
  endif(COVERAGE)

  if(GPROF)
    set(CMAKE_C_FLAGS "-pg ${CMAKE_C_FLAGS}")
    set(CMAKE_C_LINKER_FLAGS "-pg ${CMAKE_C_LINKER_FLAGS}")
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

  if(AVX512)
    set(CMAKE_C_FLAGS "-axcommon-avx512 -no-fma ${CMAKE_C_FLAGS}")
  endif(AVX512)

  if(AVX512_SKYLAKE)
    set(CMAKE_C_FLAGS "-axskylake-avx512 -no-fma ${CMAKE_C_FLAGS}")
  endif(AVX512_SKYLAKE)

  if(NATIVE)
    set(CMAKE_C_FLAGS "-xHost ${CMAKE_C_FLAGS}")
  endif(NATIVE)

  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_DEBUG   "-O0 -g")

  if(STATIC)
    set(CMAKE_C_LINKER_FLAGS "-static ${CMAKE_C_LINKER_FLAGS}")
  endif(STATIC)

endif()

# AppleClang
if(${CMAKE_C_COMPILER_ID} MATCHES "AppleClang")

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

    if(AVX512)
      set(CMAKE_C_FLAGS "-mavx -mavx2 -mavx512f -mno-fma4 -mno-fma ${CMAKE_C_FLAGS}")
    endif(AVX512)

    if(AVX512_SKYLAKE)
      set(CMAKE_C_FLAGS "-mavx -mavx2 -mavx512f -mavx512cd -mavx512bw -mavx512vl -mno-fma4 -mno-fma ${CMAKE_C_FLAGS}")
    endif(AVX512_SKYLAKE)

    if(NATIVE)
      set(CMAKE_C_FLAGS "-march=native ${CMAKE_C_FLAGS}")
    endif(NATIVE)
  endif()

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)
    set(CMAKE_C_FLAGS "-mmacosx-version-min=${MIN_OSX_VERSION} ${CMAKE_C_FLAGS}")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)

  if(STATIC)
    set(CMAKE_C_LINKER_FLAGS "-static ${CMAKE_C_LINKER_FLAGS}")
  endif(STATIC)

endif()

# CXX Flags
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")

  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES aarch64 )
    set (CMAKE_CXX_FLAGS "-ffp-contract=off ${CMAKE_CXX_FLAGS}")
  endif()
  set (CMAKE_CXX_FLAGS_RELEASE "-O3")
  set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")

  if(LTO)
    # set(CMAKE_CXX_LINK_EXECUTABLE "")
    set(CMAKE_Fortran_CXX_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto -Wl,-flto")
  endif(LTO)

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
      set(CMAKE_CXX_FLAGS "-mavx -mavx2 ${CMAKE_CXX_FLAGS}")
    endif(AVX2)

    if(AVX512)
      set(CMAKE_CXX_FLAGS "-mavx -mavx2 -mavx512f ${CMAKE_CXX_FLAGS}")
    endif(AVX512)

    if(AVX512_SKYLAKE)
      set(CMAKE_CXX_FLAGS "-mavx -mavx2 -mavx512f -mavx512cd -mavx512bw -mavx512vl ${CMAKE_CXX_FLAGS}")
    endif(AVX512_SKYLAKE)

    if(NATIVE)
      set(CMAKE_CXX_FLAGS "-march=native ${CMAKE_CXX_FLAGS}")
    endif(NATIVE)

    if(CRLIBM)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-fma4 -mno-fma")
    endif(CRLIBM)

    if(PREFETCH)
      set(CMAKE_CXX_FLAGS "-fprefetch-loop-arrays ${CMAKE_CXX_FLAGS}")
    endif(PREFETCH)
  endif()

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)
    set(CMAKE_CXX_FLAGS "-mmacosx-version-min=${MIN_OSX_VERSION} ${CMAKE_CXX_FLAGS}")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)

  if(STATIC)
    set(CMAKE_CXX_LINKER_FLAGS "-static -static-libstdc++ ${CMAKE_CXX_LINKER_FLAGS}")
  endif(STATIC)

  # Disable COVERAGE for C++ source as it breaks SixTestWrapper on some systems
  # if(COVERAGE)
  #   set(CMAKE_CXX_FLAGS "-fprofile-arcs -ftest-coverage ${CMAKE_CXX_FLAGS}")
  #   set(CMAKE_CXX_LINKER_FLAGS "-fprofile-arcs ${CMAKE_CXX_LINKER_FLAGS}")
  # endif(COVERAGE)

  if(GPROF)
    set(CMAKE_CXX_FLAGS "-pg ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_LINKER_FLAGS "-pg ${CMAKE_CXX_LINKER_FLAGS}")
  endif(GPROF)

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

  if(AVX512)
    set(CMAKE_CXX_FLAGS "-axcommon-avx512 -no-fma ${CMAKE_CXX_FLAGS}")
  endif()

  if(AVX512_SKYLAKE)
    set(CMAKE_CXX_FLAGS "-axskylake-avx512 -no-fma ${CMAKE_CXX_FLAGS}")
  endif()

  if(NATIVE)
    set(CMAKE_CXX_FLAGS "-march=xHost ${CMAKE_CXX_FLAGS}")
  endif(NATIVE)

  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")

  # When using icc remember to link to libstdc++
  if(STATIC)
    find_library(stdc++ STATIC stdc++)
    set(CMAKE_CXX_LINKER_FLAGS "-static ${CMAKE_CXX_LINKER_FLAGS}")
  else()
    find_library(stdc++ stdc++)
  endif(STATIC)

endif()

#AppleClang
if(${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang")

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

    if(AVX512)
      set(CMAKE_CXX_FLAGS "-mavx -mavx2 -mavx512f -mno-fma4 -mno-fma ${CMAKE_CXX_FLAGS}")
    endif(AVX512)

    if(AVX512_SKYLAKE)
      set(CMAKE_CXX_FLAGS "-mavx -mavx2 -mavx512f -mavx512cd -mavx512bw -mavx512vl -mno-fma4 -mno-fma ${CMAKE_CXX_FLAGS}")
    endif(AVX512_SKYLAKE)

    if(NATIVE)
      set(CMAKE_CXX_FLAGS "-march=native ${CMAKE_CXX_FLAGS}")
    endif(NATIVE)
  endif()

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)
    set(CMAKE_CXX_FLAGS "-mmacosx-version-min=${MIN_OSX_VERSION} ${CMAKE_CXX_FLAGS}")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)

  if(STATIC)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)
      #do nothing here - just to be clear c++ does not get fully statically linked
      set(CMAKE_CXX_LINKER_FLAGS "${CMAKE_CXX_LINKER_FLAGS}")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND BOINC)
  endif(STATIC)

endif()

# Check that the compiler has the features we need
if(${CMAKE_VERSION} VERSION_GREATER "3.6.0")
  include(${CMAKE_SOURCE_DIR}/cmake/fortran_support.cmake)
endif()

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
if(ZLIB)
  target_link_libraries(SixTrackCore zlibstatic)
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
  list(APPEND PREPRO_FLAGS -DSINGLE_MATH)
endif()
if(64BITM)
  message(STATUS "Adding preprocessor macro -DSIXTRACK_REAL=real64 to target SixTrack")
  list(APPEND PREPRO_FLAGS -DDOUBLE_MATH)
endif()
if(128BITM)
  message(STATUS "Adding preprocessor macro -DSIXTRACK_REAL=real128 to target SixTrack")
  list(APPEND PREPRO_FLAGS -DQUAD_MATH)
endif()

  list(APPEND PREPRO_FLAGS -DSIXTRACK_GIT_REVISION="'${GIT_REVISION}'")

#build dafor
add_executable(dafor ${CMAKE_SOURCE_DIR}/source/tools/dafor/dafor.f90)
#we do not need any fancy flags for dafor but nagfor does not like the following
#only do this on x86 + friends
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
  if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
    set_target_properties(dafor PROPERTIES COMPILE_FLAGS "-mno-avx -mno-avx2")
    set_target_properties(dafor PROPERTIES LINK_FLAGS "-mno-avx -mno-avx2")
  endif(Fortran_COMPILER_NAME MATCHES "gfortran.*")
endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)

if(CRLIBM)
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/crlibm)
endif()
if(CRLIBM OR BUILD_TESTING)
  add_subdirectory(${CMAKE_SOURCE_DIR}/source/roundctl) #Also needed by read90
endif()

if(HDF5)
  find_library(dl STATIC dl) # HDF5 has a capability to dynamically load filters
  add_library(hdf5_local_f STATIC IMPORTED)
  add_library(hdf5_local_s STATIC IMPORTED)
  add_library(hdf5_local_c STATIC IMPORTED)
  set_target_properties(hdf5_local_f PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/hdf5/build/bin/libhdf5_fortran-static.a)
  set_target_properties(hdf5_local_s PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/hdf5/build/bin/libhdf5_f90cstub-static.a)
  set_target_properties(hdf5_local_c PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/hdf5/build/bin/libhdf5-static.a)
  target_link_libraries(SixTrackCore hdf5_local_f hdf5_local_s hdf5_local_c dl)
  target_include_directories(SixTrackCore PRIVATE ${CMAKE_SOURCE_DIR}/lib/hdf5/build/bin/static)
  target_include_directories(SixTrackCore PRIVATE ${CMAKE_SOURCE_DIR}/lib/hdf5/build/fortran)
endif(HDF5)

if(HASHLIB)
  set(MD5_SOURCES
    ${CMAKE_SOURCE_DIR}/source/md5/md5.c
    ${CMAKE_SOURCE_DIR}/source/md5/md5interface.c
  )
  add_library(md5 STATIC  ${MD5_SOURCES})
  target_include_directories(md5 PRIVATE ${CMAKE_SOURCE_DIR}/source/md5/)
  set_property(TARGET md5 PROPERTY C_STANDARD 99)
  if(32BIT)
    set_target_properties(md5 PROPERTIES LINK_FLAGS "-m32")
  endif(32BIT)
  if(64BIT)
    set_target_properties(md5 PROPERTIES LINK_FLAGS "-m64")
  endif(64BIT)
  target_link_libraries(SixTrackCore md5)
endif(HASHLIB)

if(MERLINSCATTER)
  file(GLOB MERLINSCATTER_SOURCES ${CMAKE_SOURCE_DIR}/source/merlinscatter/*.cpp ${CMAKE_SOURCE_DIR}/source/merlinscatter/RandomNG/*.cpp)
  add_library(merlinscatter STATIC ${MERLINSCATTER_SOURCES})

  #Enable c++11
  set_property(TARGET merlinscatter PROPERTY CXX_STANDARD 11)

  if(32BIT)
    set_target_properties(merlinscatter PROPERTIES LINK_FLAGS "-m32")
  endif(32BIT)
  if(64BIT)
    set_target_properties(merlinscatter PROPERTIES LINK_FLAGS "-m64")
  endif(64BIT)
endif()

# Build DISTlib
if(DISTLIB)
  if(NOT EXISTS ${CMAKE_SOURCE_DIR}/lib/DISTlib/source/distinterface.c)
    if(GIT_FOUND)
      message(STATUS "Checking out DISTlib submodule")
      execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init lib/DISTlib
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT
        ERROR_QUIET
      )
      if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init lib/DISTlib failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif()
    else()
      message(FATAL_ERROR "The git submodules DISTlib is not available. Please download it into lib/DISTlib")
    endif()
  endif()

  file(GLOB DIST_SOURCES ${CMAKE_SOURCE_DIR}/lib/DISTlib/source/*.c)
  add_library(distlib STATIC ${DIST_SOURCES})
  target_include_directories(distlib PRIVATE ${CMAKE_SOURCE_DIR}/lib/DISTlib/source)
  if(CRLIBM)
    target_include_directories(distlib PRIVATE ${CMAKE_SOURCE_DIR}/source/crlibm)
    target_include_directories(distlib PRIVATE ${CMAKE_SOURCE_DIR}/source/roundctl)
    target_compile_options(distlib PRIVATE -DCRLIBM -fcommon) # GCC 10 defaults to -fno-common https://gcc.gnu.org/gcc-10/porting_to.html
    target_link_libraries(distlib crlibm roundctl)
  endif()

  if(32BIT)
    set_target_properties(distlib PROPERTIES LINK_FLAGS "-m32")
  endif(32BIT)
  if(64BIT)
    set_target_properties(distlib PROPERTIES LINK_FLAGS "-m64")
  endif(64BIT)

endif(DISTLIB)

# Build NAFFlib
if(NAFF)
  if(NOT EXISTS ${CMAKE_SOURCE_DIR}/lib/NAFFlib/NAFFlib/interface/naff_sixtrack_interface.c)
    if(GIT_FOUND)
      message(STATUS "Checking out NAFFlib submodule")
      execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init lib/NAFFlib
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT
        ERROR_QUIET
      )
      if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init lib/NAFFlib failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif()
    else()
      message(FATAL_ERROR "The git submodules NAFFlib is not available. Please download it into lib/NAFFlib")
    endif()
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
    set_target_properties(naff PROPERTIES LINK_FLAGS "-m32")
  endif(32BIT)
  if(64BIT)
    set_target_properties(naff PROPERTIES LINK_FLAGS "-m64")
  endif(64BIT)
endif(NAFF)

if(BOINC)
  find_library(stdc++ STATIC stdc++)
  include_directories(${BOINC_DIR}/api ${BOINC_DIR}/lib ${BOINC_DIR})

  add_library(boinc       STATIC IMPORTED)
  add_library(boinc_api   STATIC IMPORTED)
  add_library(boinc_api_f STATIC IMPORTED)

  set_target_properties(boinc       PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/lib/libboinc.a)
  set_target_properties(boinc_api   PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/api/libboinc_api.a)
  set_target_properties(boinc_api_f PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/api/boinc_api_fortran.o)

  target_link_libraries(SixTrackCR boinc_api_f boinc_api boinc stdc++ Threads::Threads)

  # if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  #   find_library(nsl STATIC nsl)
  #   target_link_libraries(SixTrackCR nsl)
  # endif()

  # Fix OSX linking to the correct c++ stdlib
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    find_library(c++ STATIC c++)
    target_link_libraries(SixTrackCore c++)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(BOINC)

if(FLUKA)
  if(GIT_FOUND)
    message(STATUS "Checking out FlukaIO submodule")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init lib/FlukaIO
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT
      ERROR_QUIET
    )
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init lib/FlukaIO failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  else()
    message(FATAL_ERROR "The git submodule FlukaIO is not available. Please download it into lib/FlukaIO")
  endif()

  add_subdirectory(${CMAKE_SOURCE_DIR}/lib/FlukaIO)
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

if(ZLIB)
  set(MINIZIP_DIR  ${ZLIB_DIR}/contrib/minizip)
  set(MINIZIP_WRAP ${CMAKE_SOURCE_DIR}/source/minizip)
  set(MINIZIP_SOURCES ${MINIZIP_DIR}/zip.c ${MINIZIP_DIR}/unzip.c ${MINIZIP_DIR}/ioapi.c ${MINIZIP_WRAP}/utils.c ${MINIZIP_WRAP}/libminizip.c)
  if(WIN32)
    set(MINIZIP_SOURCES ${MINIZIP_SOURCES} ${MINIZIP_DIR}/iowin32.c)
  endif(WIN32)
  include_directories(${MINIZIP_DIR})
  add_library(minizip STATIC ${MINIZIP_SOURCES})
  target_link_libraries(minizip zlibstatic)
  add_dependencies(SixTrackCore minizip)
  target_link_libraries(SixTrackCore minizip)
endif(ZLIB)

# We need to build dafor before we build SixTrack, but crlibm and so on can be built at the same time.
add_custom_target(copy_input_files)
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
    #Get the compiler version
    string(REGEX MATCH "^([0-9]*)" VERSION_TMP ${CMAKE_Fortran_COMPILER_VERSION})
    set(COMPILER_MAJOR_VERSION "${CMAKE_MATCH_1}")

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
  add_library(graflib STATIC IMPORTED)
  add_library(grafX11 STATIC IMPORTED)
  add_library(kernlib_noshift STATIC IMPORTED)
  add_library(packlib_noshift STATIC IMPORTED)
  add_library(X11 STATIC IMPORTED)
  add_library(Xau STATIC IMPORTED)
  add_library(xcb STATIC IMPORTED)
  add_library(Xdmcp STATIC IMPORTED)
  add_library(g2c STATIC IMPORTED)
  if(32BIT)
    # Use precompiled CERNLIB from lib32 folder
    set_target_properties(g2c PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libg2c.a)
    set_target_properties(graflib PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libgraflib.a)
    set_target_properties(grafX11 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libgrafX11.a)
    set_target_properties(kernlib_noshift PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libkernlib_noshift.a)
    set_target_properties(packlib_noshift PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libpacklib_noshift.a)
    set_target_properties(X11 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libX11.a)
    set_target_properties(Xau PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libXau.a)
    set_target_properties(xcb PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libxcb.a)
    set_target_properties(Xdmcp PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib32/libXdmcp.a)
  endif(32BIT)
  if(64BIT)
    # Use precompiled CERNLIB from lib64 folder
    set_target_properties(g2c PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libg2c.a)
    set_target_properties(graflib PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libgraflib.a)
    set_target_properties(grafX11 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libgrafX11.a)
    set_target_properties(kernlib_noshift PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libkernlib_noshift.a)
    set_target_properties(packlib_noshift PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libpacklib_noshift.a)
    set_target_properties(X11 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libX11.a)
    set_target_properties(Xau PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libXau.a)
    set_target_properties(xcb PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libxcb.a)
    set_target_properties(Xdmcp PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/lib64/libXdmcp.a)
  endif(64BIT)
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

if(DISTLIB)
  target_link_libraries(SixTrackCore distlib)
endif(DISTLIB)

#if we are using root, link to all the root libraries
if(ROOT)
  target_link_libraries(SixTrackCore root_output)
endif(ROOT)

if(HDF5)
  target_link_libraries(SixTrackCore ${HDF5_LIBRARIES})
endif(HDF5)

if(G4COLLIMATION)
  target_link_libraries(SixTrackCR g4collimation ${Geant4_LIBRARIES} Threads::Threads)
endif(G4COLLIMATION)

if(HASHLIB)
  target_link_libraries(SixTrackCore md5)
endif(HASHLIB)

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

#Tools
set(MOD_SPLITSINGLETRACK_SOURCES
  ${CMAKE_SOURCE_DIR}/source/splitSingletrack/mod_splitSingletrack.f90
)
add_library          (mod_splitSingletrack STATIC ${MOD_SPLITSINGLETRACK_SOURCES})
add_executable       (splitSingletrack ${CMAKE_SOURCE_DIR}/source/splitSingletrack/splitSingletrack.f90)
target_link_libraries(splitSingletrack mod_splitSingletrack)
set_target_properties(splitSingletrack PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")

# For testing
if(BUILD_TESTING)

  enable_testing()
  include(CTest)

  add_subdirectory(${CMAKE_SOURCE_DIR}/test)
  configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake  ${CMAKE_BINARY_DIR})

  #need to build the testing tools
  add_executable(verify10 ${CMAKE_SOURCE_DIR}/test/tools/verify10.f90)
  add_executable(checkf10 ${CMAKE_SOURCE_DIR}/test/tools/checkf10.f90)
  add_executable(compf10 ${CMAKE_SOURCE_DIR}/test/tools/compf10.f90)
  add_executable(checkf110 ${CMAKE_SOURCE_DIR}/test/tools/checkf110.f90)
  add_executable(checkf1014 ${CMAKE_SOURCE_DIR}/test/tools/checkf1014.f90)
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

  # The Test Wrapper
  add_executable(SixTestWrapper ${CMAKE_SOURCE_DIR}/test/tools/SixTestWrapper.cpp)
  set_target_properties(SixTestWrapper PROPERTIES LINK_FLAGS "${CMAKE_CXX_LINKER_FLAGS}")

  target_link_libraries(SixTestWrapper Threads::Threads)

  if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    target_link_libraries(SixTestWrapper stdc++)
  endif()
  set_property(TARGET SixTestWrapper PROPERTY CXX_STANDARD 11)

  #just a work around to build on MSYS2 for now
  if(MSYS AND STATIC)
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

install(TARGETS SixTrackCR RUNTIME DESTINATION bin)

if(NOT WIN32)
  add_custom_command(
    TARGET SixTrackCR POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${SIXTRACK_BINARY_NAME}" "sixtrack"
  )
  if(SIXDA)
    add_custom_command(
      TARGET SixTrackDA POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E create_symlink "${SIXTRACK_BINARY_NAME_DA}" "sixtrack_da"
    )
  endif()
endif()
