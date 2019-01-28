#Sixtrack CMAKE build file
#Making use of https://cmake.org/Wiki/CMakeForFortranExample
#and the file created for building Merlin.

###################################################################################################
# General configuration
###################################################################################################

#Set the cmake version required
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
#set(CMAKE_LEGACY_CYGWIN_WIN32 0)

# Set the project name and options
PROJECT(SixTrack)
ENABLE_LANGUAGE(Fortran)

# force an out of source build to stop clutter
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
	message(STATUS "Sixtrack requires an out-of-source build.")
	message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
	message(STATUS "CMakeCache.txt")
	message(STATUS "CMakeFiles")
	message(STATUS "Once these files are removed, create a separate directory")
	message(STATUS "and run CMake from there")
	message(FATAL_ERROR "in-source build detected")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"# MPI MPIDEBUG."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Set Sixtrack version info
file(READ ${CMAKE_SOURCE_DIR}/source/version.f90 IFILE)
STRING(REGEX MATCH ".*data version /'([0-9]+)[.]([0-9]+)[.]?([0-9]+)?[-]?(rc[0-9])?'/" VERSION_STRING ${IFILE})
SET(VERSION_STRING "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}")
SET(${PROJECT_NAME}_VERSION_MAJOR "${CMAKE_MATCH_1}")
SET(${PROJECT_NAME}_VERSION_MINOR "${CMAKE_MATCH_2}")
if ("${CMAKE_MATCH_3}" STREQUAL "")
	SET(${PROJECT_NAME}_VERSION_PATCH "0")
else()
	SET(VERSION_STRING "${VERSION_STRING}.${CMAKE_MATCH_3}")
	SET(${PROJECT_NAME}_VERSION_PATCH "${CMAKE_MATCH_3}")
endif()
string(REGEX MATCH "..$" NUM_MINOR "00${${PROJECT_NAME}_VERSION_MINOR}")
string(REGEX MATCH "..$" NUM_PATCH "00${${PROJECT_NAME}_VERSION_PATCH}")
set(NUM_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}${NUM_MINOR}${NUM_PATCH}")

if ("${CMAKE_MATCH_4}" STREQUAL "")
	SET(${PROJECT_NAME}_VERSION_RC "")
else()
	SET(VERSION_STRING "${VERSION_STRING}-${CMAKE_MATCH_4}")
	SET(${PROJECT_NAME}_VERSION_RC ${CMAKE_MATCH_4})
endif()
message(STATUS "SixTrack Version: ${VERSION_STRING}")
message(STATUS "Numeric Version:  ${NUM_VERSION}")
SET(${PROJECT_NAME}_VERSION "${VERSION_STRING}")
MESSAGE(STATUS "Running on: ${CMAKE_SYSTEM_NAME} with ${CMAKE_SYSTEM_PROCESSOR} and generating for ${CMAKE_GENERATOR}")

#get the current git commit hash
EXECUTE_PROCESS(COMMAND git log -n1 --pretty=format:%H WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE GIT_REVISION)
STRING(STRIP ${GIT_REVISION} GIT_REVISION)
message(STATUS "GIT_REVISION: ${GIT_REVISION}")

#replace the "GIT_HASH" in version.f90
#add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/source/version.f90 ${CMAKE_BINARY_DIR})
#STRING(REPLACE "GIT_HASH" "\"${GIT_REVISION}\"" VERSION_FILE_OUT ${IFILE})
#file(WRITE ${CMAKE_BINARY_DIR}/version.f90 ${VERSION_FILE_OUT})

###################################################################################################
# Sixtrack configuration options
###################################################################################################

# Now we have the default settings
# we have ON
OPTION(TILT "Allow elements to be tilted (by error)" ON)
OPTION(FAST "Which implementation of drifts to use in thin6d" ON)
OPTION(CRLIBM "Use correctly rounded libmath instead of system libmath" ON)
OPTION(SIXDA "Build differential algebra version (NOT dynamic aperture!)" ON)
OPTION(STF "Single Track File - write all tracks for postprocessing to fort.90, instead of 32 separate files fort.90 to fort.(90-32). Required for increased particle numbers" ON)
# and we have OFF
OPTION(CERNLIB "Link to external CERNLIB library for PAW plots. Otherwise use internally defined dummy functions from deck plotdumy" OFF)
OPTION(API "Link to external BOINC API (if BOINC but not API, just use internal 'myboinc.f' dummy library)" OFF)
OPTION(COLLIMAT "Collimation code" OFF)
#OPTION(CPSS "CERN physics screensaver support" OFF)
OPTION(BOINC "BOINC distributed computing support" OFF)
OPTION(CR "Checkpoint / restarting" OFF)
OPTION(BEAMGAS "Beam-gas scattering" OFF)
OPTION(DEBUG "Debugging output" OFF)
OPTION(HDF5 "Adds support for the HDF5 block which enables writing output to a single HDF5 file." OFF)
OPTION(FIO "Use FortranIO from Fortran2003 to correctly round ASCII input/output" OFF)
OPTION(MERLINSCATTER "Interaction physics for collimation from Merlin" OFF)
OPTION(NAFF "Use NAFF library for FMA analysis" OFF)
OPTION(G4COLLIMAT "Interaction physics for collimation from Geant4" OFF)
OPTION(FLUKA "Couple to FLUKA for beam collimation" OFF)
OPTION(FFIELD "Fringe Field for Quadrupole" OFF)
OPTION(COVERAGE "Enable build flags for testing code coverage with gcov (only works with GNU compilers)" OFF)
OPTION(LIBARCHIVE "Link with LIBARCHIVE, needed for BOINC API and ZIPIT block" OFF)
OPTION(GPROF "Enable build flags for code profiling with gprof (only works with GNU compilers)" OFF)
OPTION(WARN "Enable build flags for additional compiler warnings" OFF)

#Compiler stuff - not needed any more with cmake
#OPTION(IFORT "Use the intel fortran compiler" OFF)
#OPTION(LF95 "Use the Lahey-Fujitsu lf95 fortran compiler" OFF)
#OPTION(PGF90 "Use the pgf90 fortran compiler" OFF)
#OPTION(F95 "Use the g95 fortran compiler" OFF)
#OPTION(GFORTRAN "Use the gfortran fortran compiler" OFF)
#OPTION(NAGFOR "Use NAGFOR compiler" OFF)

#Build options - some are commented out for future support
#OPTION(OPENMP "Enable OpenMP support" OFF)
#OPTION(MPI "Enable MPI support" OFF)
#OPTION(CUDA "Enable CUDA support" OFF)
OPTION(32BIT "Create a 32bit binary" OFF)
OPTION(64BIT "Create a 64bit binary" ON)
OPTION(STATIC "Create a statically linked binary" ON)
OPTION(NATIVE "Enable optimizations for the current machine" OFF)
OPTION(AVX "Enable use of the Advanced Vector Extensions (AVX) instruction set - Sandy Bridge and later" OFF)
OPTION(AVX2 "Enable use of the Advanced Vector Extensions 2 (AVX2) instruction set - Haswell and later" OFF)
OPTION(AVX-512 "Enable use of the Advanced Vector Extensions 512 (AVX-512) instruction set - (Currently targeting Skylake-X)" OFF)

OPTION(BUILD_TESTING "Enable the test suite" OFF)

OPTION(32BITM "Floats are 32bit (single precision)" OFF)
OPTION(64BITM "Floats are 64bit (double precision)" ON)
OPTION(128BITM "Floats are 128bit (quad precision)" OFF)

OPTION(ROUND_NEAR "Always round to the nearest number after floating point operations" ON)
OPTION(ROUND_UP "Always round up after floating point operations" OFF)
OPTION(ROUND_DOWN "Always round down after floating point operations" OFF)
OPTION(ROUND_ZERO "Always round towards zero after floating point operations" OFF)

OPTION(ROOT "Link to the root libraries" OFF)

SET(BOINC_DIR ${CMAKE_SOURCE_DIR}/source/boinc CACHE PATH "Path to the boinc build folder")
SET(LIBARCHIVE_BUILD_DIR ${CMAKE_SOURCE_DIR}/source/libarchive_build CACHE PATH "Path to the libarchive build folder")
SET(LIBARCHIVE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/source/libarchive CACHE PATH "Path to the libarchive source folder")

###################################################################################################
# Sixtrack binary name construction
###################################################################################################

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

SET(SIXTRACK_FEAT_CR "${CMAKE_SYSTEM_NAME}_${Fortran_COMPILER_NAME}")
SET(SIXTRACK_FEAT_DA "${CMAKE_SYSTEM_NAME}_${Fortran_COMPILER_NAME}")

if(NOT TILT)
	SET(SIXTRACK_FEAT_CR "notilt_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "notilt_${SIXTRACK_FEAT_DA}")
endif()
if(NOT FAST)
	SET(SIXTRACK_FEAT_CR "nofast_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "nofast_${SIXTRACK_FEAT_DA}")
endif()
if(CRLIBM AND ROUND_NEAR)
	SET(SIXTRACK_FEAT_CR "crlibm_rn_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "crlibm_rn_${SIXTRACK_FEAT_DA}")
endif(CRLIBM AND ROUND_NEAR)
if(CRLIBM AND ROUND_UP)
	SET(SIXTRACK_FEAT_CR "crlibm_ru_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "crlibm_ru_${SIXTRACK_FEAT_DA}")
endif(CRLIBM AND ROUND_UP)
if(CRLIBM AND ROUND_DOWN)
	SET(SIXTRACK_FEAT_CR "crlibm_rd_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "crlibm_rd_${SIXTRACK_FEAT_DA}")
endif(CRLIBM AND ROUND_DOWN)
if(CRLIBM AND ROUND_ZERO)
	SET(SIXTRACK_FEAT_CR "crlibm_rz_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "crlibm_rz_${SIXTRACK_FEAT_DA}")
endif(CRLIBM AND ROUND_ZERO)
if(CERNLIB)
	SET(SIXTRACK_FEAT_CR "cernlib_${SIXTRACK_FEAT_CR}")
endif()
if(API)
	SET(SIXTRACK_FEAT_CR "api_${SIXTRACK_FEAT_CR}")
endif()
if(COLLIMAT)
	SET(SIXTRACK_FEAT_CR "collimat_${SIXTRACK_FEAT_CR}")
endif()
if(BOINC)
	SET(SIXTRACK_FEAT_CR "boinc_${SIXTRACK_FEAT_CR}")
endif()
if(CR)
	SET(SIXTRACK_FEAT_CR "cr_${SIXTRACK_FEAT_CR}")
endif()
if(BEAMGAS)
	SET(SIXTRACK_FEAT_CR "beamgas_${SIXTRACK_FEAT_CR}")
endif()
if(STF)
	SET(SIXTRACK_FEAT_CR "stf_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "stf_${SIXTRACK_FEAT_DA}")
endif()
if(DEBUG)
	SET(SIXTRACK_FEAT_CR "debug_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "debug_${SIXTRACK_FEAT_DA}")
endif()
if(HDF5)
	SET(SIXTRACK_FEAT_CR "hdf5_${SIXTRACK_FEAT_CR}")
endif()
if(FIO)
	SET(SIXTRACK_FEAT_CR "fio_${SIXTRACK_FEAT_CR}")
	SET(SIXTRACK_FEAT_DA "fio_${SIXTRACK_FEAT_DA}")
endif()
if(MERLINSCATTER)
	SET(SIXTRACK_FEAT_CR "merlinscatter_${SIXTRACK_FEAT_CR}")
endif()
if(NAFF)
	SET(SIXTRACK_FEAT_CR "naff_${SIXTRACK_FEAT_CR}")
endif()
if(G4COLLIMAT)
	SET(SIXTRACK_FEAT_CR "g4collimat_${SIXTRACK_FEAT_CR}")
endif()
if(FLUKA)
	SET(SIXTRACK_FEAT_CR "FLUKA_${SIXTRACK_FEAT_CR}")
endif(FLUKA)
if(FFIELD)
	SET(SIXTRACK_FEAT_CR "FFIELD_${SIXTRACK_FEAT_CR}")
endif(FFIELD)
if(LIBARCHIVE)
	SET(SIXTRACK_FEAT_CR "libarchive_${SIXTRACK_FEAT_CR}")
endif()

if(ROOT)
	SET(SIXTRACK_FEAT_CR "root_${SIXTRACK_FEAT_CR}")
endif(ROOT)

#if(OPENMP)
#	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_openmp")
#endif()
#if(MPI)
#	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_mpi")
#endif()
#if(CUDA)
#	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_cuda")
#endif()
if(STATIC)
	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_static")
	SET(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_static")
endif()
if(AVX)
	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_avx")
	SET(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_avx")
endif(AVX)
if(AVX2)
	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_avx2")
	SET(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_avx2")
endif(AVX2)
if(AVX-512)
	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_avx512")
	SET(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_avx512")
endif(AVX-512)
if(NATIVE)
	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_native")
	SET(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_native")
endif()

SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_${CMAKE_SYSTEM_PROCESSOR}")
SET(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_${CMAKE_SYSTEM_PROCESSOR}")

if(32BIT)
	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_32bit")
	SET(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_32bit")
endif()
if(64BIT)
	SET(SIXTRACK_FEAT_CR "${SIXTRACK_FEAT_CR}_64bit")
	SET(SIXTRACK_FEAT_DA "${SIXTRACK_FEAT_DA}_64bit")
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
endif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")

# Create the sixtrack binary name string
set(SIXTRACK_BINARY_NAME    "${PROJECT_NAME}_${NUM_VERSION}_${SIXTRACK_FEAT_CR}")
set(SIXTRACK_BINARY_NAME_DA "SixDA_${NUM_VERSION}_${SIXTRACK_FEAT_DA}")

###################################################################################################
# Setting the c pre-processor definitons
###################################################################################################

LIST(APPEND PREPRO_FLAGS_TOOLS TILT FAST COLLIMAT BOINC BEAMGAS STF MERLINSCATTER NAFF G4COLLIMAT FLUKA CR FFIELD)
LIST(APPEND PREPRO_FLAGS_LIB CRLIBM ROUND_NEAR ROUND_UP ROUND_DOWN ROUND_ZERO CERNLIB DEBUG HDF5 FIO LIBARCHIVE ROOT)
LIST(APPEND PREPRO_FLAGS_COMPILE NAGFOR)

###################################################################################################
# Select which fortran files to generate and compile
###################################################################################################

list(APPEND FORT90_LIB
	sixtrack core_tools constants sixtrack_input mod_units
	common_modules close strings string_tools mainarrays
	file_units mod_alloc zipf matrix_inv lielib beam6d
	bouncy_castle plato_seq bdex fma dabnew end_sixtrack
	ranecu wire hions elens mod_dist scatter dump dynk
	aperture mod_fluc beam6d_fox six_fox postprocessing
	read_input
)

# Files unique to DA or CR version
list(APPEND FORT90_DA main_da sixda)
list(APPEND FORT90_CR main_cr sixve track_thin track_thick)

# List of files that must be treated by DAFOR
list(APPEND DAFOR_F90 beam6d_fox six_fox sixda)

if(CR)
	list(APPEND FORT90_LIB checkpoint_restart)
endif(CR)

if(DEBUG)
	list(APPEND FORT90_LIB debug_dumps)
endif(DEBUG)

if(COLLIMAT)
	list(APPEND FORT90_LIB ranlux collimation)
endif(COLLIMAT)

if(HDF5)
	list(APPEND FORT90_LIB hdf5_output hdf5_sixtrack)
endif(HDF5)

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

###################################################################################################
# Checking if some requested packages exist
###################################################################################################
#Check for HDF5
if(HDF5)
	find_package(HDF5 COMPONENTS Fortran REQUIRED)
	if(STATIC)
		SET(HDF5_USE_STATIC_LIBRARIES "True")
	endif(STATIC)
endif(HDF5)

if(HDF5_FOUND)
	#uncomment if you want more information
	message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
	message(STATUS "HDF5_INCLUDE_DIR: ${HDF5_INCLUDE_DIR}")
	#message(STATUS "HDF5_DEFINITIONS: ${HDF5_DEFINITIONS}")
	message(STATUS "HDF5_Fortran_LIBRARIES: ${HDF5_Fortran_LIBRARIES}")
	message(STATUS "HDF5_Fortran_HL_LIBRARIES: ${HDF5_Fortran_HL_LIBRARIES}")
	#message(STATUS "HDF5_HL_LIBRARIES: ${HDF5_HL_LIBRARIES}")
	#message(STATUS "HDF5_C_LIBRARIES: ${HDF5_C_LIBRARIES}")
	#message(STATUS "HDF5_CXX_LIBRARIES: ${HDF5_CXX_LIBRARIES}")
	message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
	message(STATUS "HDF5_FOUND: ${HDF5_FOUND}")
	message(STATUS "HDF5_VERSION: ${HDF5_VERSION}")
	#message(STATUS "HDF5_LIBRARY_DIRS: ${HDF5_LIBRARY_DIRS}")
	#message(STATUS "HDF5_IS_PARALLEL: ${HDF5_IS_PARALLEL}")
	message(STATUS "HDF5_Fortran_COMPILER_EXECUTABLE: ${HDF5_Fortran_COMPILER_EXECUTABLE}")
endif(HDF5_FOUND)

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
#	include(${Geant4_USE_FILE})
	add_subdirectory(${CMAKE_SOURCE_DIR}/source/g4collimat)
endif(G4COLLIMAT)

#Check for OpenMP
#if(OPENMP)
#	find_package(OpenMP REQUIRED)
#endif(OPENMP)
#
##If we have OpenMP enabled, set the compiler flags needed
#if(OPENMP_FOUND)
#	SET_TARGET_PROPERTIES(SixTrack PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
#	SET_TARGET_PROPERTIES(SixTrack PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
#endif()
#
##Possible future  MPI test
#if(CMAKE_BUILD_TYPE STREQUAL MPI OR CMAKE_BUILD_TYPE STREQUAL MPIDEBUG)
#	FIND_PACKAGE(MPI REQUIRED)
#	if(NOT MPI_CXX_FOUND)
#		MESSAGE(FATAL_ERROR "MPI build requested but no MPI libraries found!")
#	endif()
#	if(MPI_CXX_FOUND)
#		message(STATUS "Using MPI C++ compiler: ${MPI_CXX_COMPILER}")
#		message(STATUS "Using MPI runtime: ${MPIEXEC}")
#		message(STATUS "Using MPI process number switch: ${MPIEXEC_NUMPROC_FLAG}")
#		message(STATUS "Using MPI pre-flags: ${MPIEXEC_PREFLAGS}")
#		message(STATUS "Using MPI post-flags: ${MPIEXEC_POSTFLAGS}")
#		message(STATUS "Using MPI includes: ${MPI_CXX_INCLUDE_PATH}")
#		message(STATUS "Using MPI libraries: ${MPI_CXX_LIBRARIES}")
#		message(STATUS "Using MPI compile flags: ${MPI_CXX_COMPILE_FLAGS}")
#		message(STATUS "Using MPI link flags: ${MPI_CXX_LINK_FLAGS}")
#
#		#We need the MPI headers to be included!
#		include_directories(${MPI_CXX_INCLUDE_PATH})
#	endif()
#endif()
#
##possible future CUDA test.
#if(CUDA)
#	find_package(CUDA REQUIRED)
#	if(NOT CUDA_FOUND)
#		MESSAGE(FATAL_ERROR "CUDA enabled build requested but CUDA was not found!")
#	endif()
#	if(CUDA_FOUND)
#		message(STATUS "Found CUDA version ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")
#		message(STATUS "Includes at ${CUDA_INCLUDE_DIRS}")
#		ADD_DEFINITIONS("-I${CUDA_INCLUDE_DIRS}")
#	endif()
#endif(CUDA)
#
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
if(BEAMGAS)
	if(NOT COLLIMAT)
		message(FATAL_ERROR "Building with BEAMGAS requires COLLIMAT enabled")
	endif(NOT COLLIMAT)
endif(BEAMGAS)
if(COLLIMAT)
	if(CR)
		message(FATAL_ERROR "Cannot build COLLIMAT with CR enabled.")
	endif()
	if(SIXDA)
		message(STATUS "Will not build SixDA when COLLIMAT is enabled.")
		set(SIXDA OFF)
	endif()
endif()
if(FLUKA)
	if(SIXDA)
		message(STATUS "Will not build SixDA when FLUKA is enabled.")
		set(SIXDA OFF)
	endif()
endif()
if(MERLINSCATTER)
	if(NOT COLLIMAT)
		message(FATAL_ERROR "Building with MERLINSCATTER requires COLLIMAT enabled.")
	endif()
endif()
if(G4COLLIMAT)
	if(NOT COLLIMAT)
		message(FATAL_ERROR "Building with G4COLLIMAT requires COLLIMAT enabled.")
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
	
	SET(CMAKE_Fortran_FLAGS "-frecord-marker=4 -fno-second-underscore -funroll-loops -std=f2008")
	set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
	set(PREPRO_FLAGS -cpp -E -P)

	if(${CMAKE_SYSTEM_PROCESSOR} MATCHES aarch64 )
		set (CMAKE_Fortran_FLAGS "-ffp-contract=off ${CMAKE_Fortran_FLAGS}")
	endif()

	SET(CMAKE_Fortran_FLAGS_RELEASE "-O3")
	SET(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all")

	if(${CMAKE_SYSTEM_PROCESSOR} MATCHES ppc64le )
		SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fno-expensive-optimizations -fno-inline-functions")
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
			SET(CMAKE_Fortran_FLAGS "-march=native -mno-fma4 -mno-fma ${CMAKE_Fortran_FLAGS}")
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
	
	set(CMAKE_Fortran_FLAGS_RELEASE "-fp-model source -fp-model strict -fp-model no-except -O3 -no-fma -std15 -assume realloc-lhs")
	set(CMAKE_Fortran_FLAGS_DEBUG   "-fp-model source -fp-model strict -fp-model no-except -O0 -g -no-fma -std15 -assume realloc-lhs")
	set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
	set(PREPRO_FLAGS -fpp -E -P)

	if(WARN)
		set(CMAKE_Fortran_FLAGS "-warn unused ${CMAKE_Fortran_FLAGS}")
	endif(WARN)
	
	if(32BIT)
		set(CMAKE_Fortran_FLAGS_RELEASE "-m32 ${CMAKE_Fortran_FLAGS_RELEASE}")
		set(CMAKE_Fortran_FLAGS_DEBUG "-m32 ${CMAKE_Fortran_FLAGS_DEBUG}")
		SET(CMAKE_Fortran_LINKER_FLAGS "-m32 ${CMAKE_Fortran_LINKER_FLAGS}")
	endif(32BIT)

	if(64BIT)
		set(CMAKE_Fortran_FLAGS_RELEASE "-m64 ${CMAKE_Fortran_FLAGS_RELEASE}")
		set(CMAKE_Fortran_FLAGS_DEBUG "-m64 ${CMAKE_Fortran_FLAGS_DEBUG}")
		SET(CMAKE_Fortran_LINKER_FLAGS "-m64 ${CMAKE_Fortran_LINKER_FLAGS}")
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
		SET(CMAKE_Fortran_FLAGS_RELEASE "-xHost ${CMAKE_Fortran_FLAGS_RELEASE}")
		#Nothing to do for DEBUG
	endif(NATIVE)

	if(STATIC)
		SET(CMAKE_Fortran_LINKER_FLAGS "-static ${CMAKE_Fortran_LINKER_FLAGS}")
	endif(STATIC)

elseif (Fortran_COMPILER_NAME MATCHES "nagfor.*")
	############
	#  nagfor  #
	############
	
	SET(NAGFOR ON)
	#Here I am going to assume that since we are using nagfor we are on lxplus for now
	set(CMAKE_Fortran_FLAGS_RELEASE "-O4 -dusty -ieee=full -dcfuns -w=ulv -w=uparam")
	set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -dusty -ieee=full -dcfuns -w=ulv -w=uparam -gline")
	set(PREPRO_EXEC  fpp)
	set(PREPRO_FLAGS -P -free)

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
	message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})

	message ("No optimized Fortran compiler flags are known, we just try -O2 or -O0 -g...")
	set (CMAKE_Fortran_FLAGS_RELEASE "-O2 ${CMAKE_Fortran_FLAGS_RELEASE}")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g ${CMAKE_Fortran_FLAGS_DEBUG}")
	set(PREPRO_EXEC  ${CMAKE_Fortran_COMPILER})
	set(PREPRO_FLAGS -cpp -E -P)
endif ()

#Now for C flags (crlibm for now)
#${CMAKE_C_COMPILER_ID} : From the docs: "one of "Clang", "GNU", "Intel", or "MSVC". This works even if a compiler wrapper like ccache is used."
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")

	if(${CMAKE_SYSTEM_PROCESSOR} MATCHES aarch64 )
		set (CMAKE_C_FLAGS "-ffp-contract=off ${CMAKE_C_FLAGS}")
	endif()
	set (CMAKE_C_FLAGS_RELEASE "-O2")
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
    include(${CMAKE_SOURCE_DIR}/fortran_support.cmake)

#MPI default build options
#SET( CMAKE_CXX_FLAGS_MPI "-DNDEBUG -Wall -O3 -pedantic -DENABLE_MPI -march=native" CACHE STRING
#	"Flags used by the C++ compiler during MPI builds."
#	FORCE )
#SET( CMAKE_EXE_LINKER_FLAGS_MPI
#	"-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
#	"Flags used for linking binaries during MPI builds."
#	FORCE )
#SET( CMAKE_SHARED_LINKER_FLAGS_MPI
#	"-Wl,--warn-unresolved-symbols,--warn-once -fPIC -flto" CACHE STRING
#	"Flags used by the shared libraries linker during MPI builds."
#	FORCE )
#MARK_AS_ADVANCED(
#	CMAKE_CXX_FLAGS_MPI
#	#CMAKE_C_FLAGS_MPI
#	CMAKE_EXE_LINKER_FLAGS_MPI
#	CMAKE_SHARED_LINKER_FLAGS_MPI)

#MPI debug default build options
#SET( CMAKE_CXX_FLAGS_MPIDEBUG "-Wall -g3 -ggdb -pedantic -DENABLE_MPI -march=native" CACHE STRING
#	"Flags used by the C++ compiler during MPI Debug builds."
#	FORCE )
#SET( CMAKE_EXE_LINKER_FLAGS_MPIDEBUG
#	"-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
#	"Flags used for linking binaries during MPI Debug builds."
#	FORCE )
#SET( CMAKE_SHARED_LINKER_FLAGS_MPIDEBUG
#	"-Wl,--warn-unresolved-symbols,--warn-once -fPIC" CACHE STRING
#	"Flags used by the shared libraries linker during MPI Debug builds."
#	FORCE )
#MARK_AS_ADVANCED(
#	CMAKE_CXX_FLAGS_MPIDEBUG
#	#CMAKE_C_FLAGS_MPIDEBUG
#	CMAKE_EXE_LINKER_FLAGS_MPIDEBUG
#	CMAKE_SHARED_LINKER_FLAGS_MPIDEBUG )

###################################################################################################
# Settings for what exactly to build
###################################################################################################
# Add the source files that will be used to build SixTrack
add_library(SixTrackLib STATIC ${SIXTRACK_SOURCES})
set_target_properties(SixTrackLib PROPERTIES OUTPUT_NAME "SixTrack")
add_executable(SixTrackCR ${SIXTRACK_SOURCES_CR})
target_link_libraries(SixTrackCR SixTrackLib)
set_target_properties(SixTrackCR PROPERTIES OUTPUT_NAME "${SIXTRACK_BINARY_NAME}" LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}" LINKER_LANGUAGE Fortran)
if(SIXDA)
	add_executable(SixTrackDA ${SIXTRACK_SOURCES_DA})
	target_link_libraries(SixTrackDA SixTrackLib)
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
	list(APPEND PREPRO_FLAGS -DSIXTRACK_REAL=real32)
endif()
if(64BITM)
	message(STATUS "Adding preprocessor macro -DSIXTRACK_REAL=real64 to target SixTrack")
	list(APPEND PREPRO_FLAGS -DSIXTRACK_REAL=real64)
endif()
if(128BITM)
	message(STATUS "Adding preprocessor macro -DSIXTRACK_REAL=real128 to target SixTrack")
	list(APPEND PREPRO_FLAGS -DSIXTRACK_REAL=real128)
endif()

	list(APPEND PREPRO_FLAGS -DSIXTRACK_GIT_REVISION="'${GIT_REVISION}'")

#build dafor
add_executable(dafor ${CMAKE_SOURCE_DIR}/source/tools/dafor/dafor.f90)
#we do not need any fancy flags for dafor
#but nagfor does not like the following
if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
	SET_TARGET_PROPERTIES(dafor PROPERTIES COMPILE_FLAGS "-mno-avx -mno-avx2 -mno-avx512f -mno-avx512cd -mno-avx512vl -mno-avx512bw -mno-avx512dq")
	SET_TARGET_PROPERTIES(dafor PROPERTIES LINK_FLAGS "-mno-avx -mno-avx2 -mno-avx512f -mno-avx512cd -mno-avx512vl -mno-avx512bw -mno-avx512dq")
	#target_compile_definitions(dafor PUBLIC "-mno-avx -mno-avx2 -mno-avx512f -mno-avx512cd -mno-avx512vl -mno-avx512bw -mno-avx512dq")
endif(Fortran_COMPILER_NAME MATCHES "gfortran.*")

if(CRLIBM)
	add_subdirectory(${CMAKE_SOURCE_DIR}/source/crlibm)
endif()
if(CRLIBM OR BUILD_TESTING )
	add_subdirectory(${CMAKE_SOURCE_DIR}/source/roundctl) #Also needed by read90
endif()

if(HDF5)
	if(HDF5_FOUND)
		SET_TARGET_PROPERTIES(SixTrackLib PROPERTIES INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}")
		SET_TARGET_PROPERTIES(SixTrackLib PROPERTIES COMPILE_DEFINITIONS "${HDF5_DEFINITIONS}")
	else()
		MESSAGE(FATAL_ERROR "HDF5 build requested, but could not find the HDF5 package.")
	endif(HDF5_FOUND)
endif(HDF5)

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

if(NAFF)
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/source/cmake_modules/")
	if (STATIC)
		SET(FFTW_USE_STATIC_LIBS 1)
	endif()
	find_package(FFTW REQUIRED)
	include_directories( ${FFTW_INCLUDES} )
	message(STATUS "Include dirs of FFTW: " ${FFTW_INCLUDES} )
	message(STATUS "Libs of of FFTW:      " ${FFTW_LIBRARIES} )
	
	
	file(GLOB NAFF_SOURCES ${CMAKE_SOURCE_DIR}/source/naff/*.cpp ${CMAKE_SOURCE_DIR}/source/naff/*.cc)
	add_library(naff STATIC ${NAFF_SOURCES})
	target_link_libraries(naff ${FFTW_LIBRARIES})
	
	#Enable c++11
	set_property(TARGET naff PROPERTY CXX_STANDARD 11)

	if(32BIT)
		SET_TARGET_PROPERTIES(naff PROPERTIES LINK_FLAGS "-m32")
	endif(32BIT)
	if(64BIT)
		SET_TARGET_PROPERTIES(naff PROPERTIES LINK_FLAGS "-m64")
	endif(64BIT)
endif()


add_custom_target(copy_input_files)
add_dependencies(copy_input_files include_mkdir)

add_custom_target(include_mkdir)
file(GLOB SIXTRACK_INCL ${CMAKE_SOURCE_DIR}/source/include/*.f90)
add_custom_command(TARGET include_mkdir POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/include/)
foreach(loop IN ITEMS ${SIXTRACK_INCL})
	add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${loop} ${CMAKE_BINARY_DIR}/include/)
endforeach()
#add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/source/version.f90 ${CMAKE_BINARY_DIR})

if(BOINC AND API)
	include_directories(${BOINC_DIR}/api ${BOINC_DIR}/lib ${BOINC_DIR})
	
	#The fortran interface needs boinc_api_fortran.o building, which does not seem to happen by default
	#go into the boinc source folder, cd api; make boinc_api_fortran.o
	add_library(boinc_api_fortran STATIC IMPORTED)
	set_target_properties(boinc_api_fortran PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/api/boinc_api_fortran.o)
	# add_library(boinc_api_fortran STATIC ${BOINC_DIR}/api/boinc_api_fortran.cpp) #Alternative way of building and linking boinc_api_fortran.cpp/.o
	target_link_libraries(SixTrackCR boinc_api_fortran)
	
	#add the libs
	add_library(boinc STATIC IMPORTED)
	set_target_properties(boinc PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/lib/libboinc.a)
	
	add_library(boinc_api STATIC IMPORTED)
	set_target_properties(boinc_api PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/api/libboinc_api.a)
	target_link_libraries(SixTrackCR boinc_api)
	
	target_link_libraries(SixTrackCR boinc)
	
	find_library(stdc++ STATIC stdc++)
	target_link_libraries(SixTrackLib stdc++ Threads::Threads)

	if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
		find_library(nsl STATIC nsl)
		target_link_libraries(SixTrackCR nsl)
	endif()

	#fix osx linking to the correct c++ stdlib
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		find_library(c++ STATIC c++)
		target_link_libraries(SixTrackCR c++)
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif()

if(FLUKA)
	add_subdirectory(${CMAKE_SOURCE_DIR}/source/FlukaIO)
	target_link_libraries(SixTrackLib flukaIO)
endif(FLUKA)

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
	target_link_libraries(SixTrackCR libArchive_wrapper)
	add_library(libarchive STATIC IMPORTED)
	if(WIN32)
		set_target_properties(libarchive PROPERTIES IMPORTED_LOCATION ${LIBARCHIVE_BUILD_DIR}/libarchive/libarchive_static.a)
		add_definitions( -DLIBARCHIVE_STATIC )
	else()
		set_target_properties(libarchive PROPERTIES IMPORTED_LOCATION ${LIBARCHIVE_BUILD_DIR}/libarchive/libarchive.a)
	endif()
	
	#find_library(z STATIC z)
	add_library(z STATIC IMPORTED)
	set_target_properties(z PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/source/zlib/install/lib/libz.a)
	
	target_link_libraries(SixTrackCR libarchive z Threads::Threads)
	target_link_libraries(libArchiveTester libarchive z Threads::Threads)
	
	add_definitions( -DLIBARCHIVE )
endif()

#We need to build dafor before we build SixTrack, but crlibm and so on can be built at the same time.
add_dependencies(SixTrackLib copy_input_files generate_input)
add_dependencies(SixTrackCR SixTrackLib)
if(SIXDA)
	add_dependencies(SixTrackDA SixTrackLib)
endif()

###################################################################################################
# Running dafor
###################################################################################################

add_custom_target(generate_input)
add_custom_target(run_prepro)
add_dependencies(run_prepro dafor)
add_dependencies(generate_input run_prepro)

## Run the preprocessor on all .f90 files and output them to the build directory
foreach(loop IN ITEMS ${FORT90_LIB} ${FORT90_CR} ${FORT90_DA})
	list(FIND DAFOR_F90 ${loop} DAFOR_IDX)
	if(${DAFOR_IDX} GREATER -1)
		message(STATUS "Will run DAFOR to generate ${loop}.f90")
		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f90
			COMMAND ${PREPRO_EXEC} ${PREPRO_FLAGS} ${CMAKE_SOURCE_DIR}/source/${loop}.f90 > ${CMAKE_CURRENT_BINARY_DIR}/${loop}a.f90
			COMMAND $<TARGET_FILE:dafor> ${loop}a.f90 ${loop}.f90 > ${loop}_dafor.log
			DEPENDS run_prepro
		)
	else()
		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f90
			COMMAND ${PREPRO_EXEC} ${PREPRO_FLAGS} ${CMAKE_SOURCE_DIR}/source/${loop}.f90 > ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f90
			DEPENDS run_prepro
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
	target_link_libraries(SixTrackLib ${MSYS_QUADMATH})
endif()

#same with OSX, using gfortran from homebrew
if(STATIC)
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		find_library(lquadmath NAMES libquadmath.a PATHS /usr/local/opt/gcc/lib/gcc/6/ /usr/local/opt/gcc/lib/gcc/7/ /usr/local/opt/gcc/lib/gcc/8/)
		target_link_libraries(SixTrackLib ${lquadmath})
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(STATIC)

if(CERNLIB)
#	if(EXISTS "/cern/pro/bin/cernlib")
		#running on lxplus and friends
#		SET_TARGET_PROPERTIES(SixTrack PROPERTIES LINK_FLAGS "`cernlib -u graflib packlib_noshift kernlib_noshift` -lg2c -lpthread")
#	else()
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

		target_link_libraries(SixTrackLib graflib grafX11 packlib_noshift kernlib_noshift g2c X11 xcb Xau Xdmcp dl Threads::Threads gcc_eh)
#	endif()
endif(CERNLIB)

#if we are using crlibm, link to the crlibm library
if(CRLIBM)
	target_link_libraries(SixTrackLib crlibm roundctl)
endif(CRLIBM)

#if we are using merlinscatter, link to the merlinscatter library
if(MERLINSCATTER)
	target_link_libraries(SixTrackCR merlinscatter)
endif(MERLINSCATTER)
#if we are using naff, link to the naff library
if(NAFF)
	target_link_libraries(SixTrackLib naff)
endif(NAFF)

#if we are using root, link to all the root libraries
if(ROOT)
	target_link_libraries(SixTrackLib root_output)
endif(ROOT)

if(HDF5)
	target_link_libraries(SixTrackLib ${HDF5_LIBRARIES})
endif(HDF5)

if(G4COLLIMAT)
	target_link_libraries(SixTrackCR g4collimat ${Geant4_LIBRARIES} Threads::Threads)
endif(G4COLLIMAT)

#Again for naglib, I assume we are on lxplus or have CERN afs access.
# if(NAGLIB)
# 	if(Fortran_COMPILER_NAME MATCHES "gfortran.*")
# 		message(WARNING "Compiling with NAGLIB and gfortran may cause SixTrack to segfault.")
# 	endif()
	
# 	if(64BIT)
# 		if(STATIC)
# 			if(NOT Fortran_COMPILER_NAME MATCHES "ifort.*")
# 				find_library(nglibifcoremt NAMES libifcoremt.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			endif()
# 			find_library(nglibsvml NAMES libsvml.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			find_library(nglibirng NAMES libirng.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			find_library(nglibimf NAMES libimf.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			find_library(nglibirc NAMES libirc.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			find_library(nglib NAMES libnag_nag.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/lib)
# 		else()
# 			if(NOT Fortran_COMPILER_NAME MATCHES "ifort.*")
# 				find_library(nglibifcoremt NAMES libifcoremt.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			endif()
# 			find_library(nglibsvml NAMES libsvml.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			find_library(nglibirng NAMES libirng.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			find_library(nglibimf NAMES libimf.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			find_library(nglibirc NAMES libirc.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
# 			find_library(nglib NAMES libnag_nag.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/lib)
# 		endif(STATIC)
# 	endif(64BIT)

# 	if(32BIT)
# 		if(STATIC)
# 			find_library(nglibifcoremt NAMES libifcoremt.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglibsvml NAMES libsvml.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglibirng NAMES libirng.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglibimf NAMES libimf.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglibirc NAMES libirc.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglib NAMES libnag_nag.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/lib)
# 		else()
# 			find_library(nglibifcoremt NAMES libifcoremt.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglibsvml NAMES libsvml.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglibirng NAMES libirng.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglibimf NAMES libimf.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglibirc NAMES libirc.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
# 			find_library(nglib NAMES libnag_nag.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/lib)
# 		endif(STATIC)
# 	endif(32BIT)

# 	target_link_libraries(SixTrack ${nglib} ${nglibifcoremt} ${nglibsvml} ${nglibirng} ${nglibimf} ${nglibirc} dl pthread )
# endif(NAGLIB)

IF(COVERAGE)
	target_compile_options(SixTrackLib PRIVATE -fprofile-arcs -ftest-coverage)
	target_compile_options(SixTrackCR  PRIVATE -fprofile-arcs -ftest-coverage)
	SET_TARGET_PROPERTIES(SixTrackLib  PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(SixTrackCR   PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	if(SIXDA)
		target_compile_options(SixTrackDA PRIVATE -fprofile-arcs -ftest-coverage)
		SET_TARGET_PROPERTIES(SixTrackDA  PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	endif()
ENDIF(COVERAGE)

SET(BUILDNAME "${SIXTRACK_BINARY_NAME}" CACHE STRING "Name of build on the dashboard")
MARK_AS_ADVANCED(BUILDNAME)

#For testing
IF(BUILD_TESTING)

	ENABLE_TESTING()
	INCLUDE(CTest)

	add_subdirectory(${CMAKE_SOURCE_DIR}/test)
	configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake  ${CMAKE_BINARY_DIR})

	#Checks of specialized input/canoncials
	# Anyway hardcoded to use them from the SixTest source folder...
	#add_executable(CheckTestInputForBOINC.sh IMPORTED)
	#SET_TARGET_PROPERTIES(CheckTestInputForBOINC.sh PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/test/CheckTestInputForBOINC.sh)

	#add_executable(CheckTestInputForSTF.sh IMPORTED)
	#SET_TARGET_PROPERTIES(CheckTestInputForSTF.sh PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/test/CheckTestInputForBOINC.sh)

	#need to build the testing tools
	add_executable(verify10 ${CMAKE_SOURCE_DIR}/test/bin/verify10.f90)
	add_executable(checkf10 ${CMAKE_SOURCE_DIR}/test/bin/checkf10.f90)
	add_executable(compf10 ${CMAKE_SOURCE_DIR}/test/bin/compf10.f90)
	add_executable(checkf110 ${CMAKE_SOURCE_DIR}/test/bin/checkf110.f90)
	add_executable(checkf1014 ${CMAKE_SOURCE_DIR}/test/bin/checkf1014.f90)
	add_executable(read90 ${CMAKE_SOURCE_DIR}/test/read90/read90.f90)
	add_executable(readDump3 ${CMAKE_SOURCE_DIR}/test/readDump3/readDump3.f90)

	target_link_libraries(read90 roundctl)

	SET_TARGET_PROPERTIES(verify10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(checkf10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(compf10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(checkf110 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(checkf1014 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(read90 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(readDump3 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")

	#The testing wrapper
	add_executable(SixTestWrapper ${CMAKE_SOURCE_DIR}/test/bin/SixTestWrapper.cpp)
	if (LIBARCHIVE)
		target_link_libraries(SixTestWrapper libArchive_wrapper libarchive z Threads::Threads)
	endif()
	SET_TARGET_PROPERTIES(SixTestWrapper PROPERTIES LINK_FLAGS "${CMAKE_CXX_LINKER_FLAGS}")

	target_link_libraries(SixTestWrapper Threads::Threads)

	IF(COVERAGE)
		target_compile_options(SixTestWrapper PRIVATE -fprofile-arcs -ftest-coverage)
		SET_TARGET_PROPERTIES(SixTestWrapper PROPERTIES LINK_FLAGS "-fprofile-arcs ${CMAKE_CXX_LINKER_FLAGS}")
	ENDIF(COVERAGE)

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
ENDIF()

#To make packages, run "make package" in the build directory
#SET(CPACK_GENERATOR RPM DEB STGZ TGZ TBZ2 TXZ ZIP 7Z)
SET(CPACK_GENERATOR RPM DEB TGZ ZIP) #Take only the most common targets
#SET(CPACK_PACKAGE_CONTACT "sixtrack.support@cern.ch") # non-existent email address
SET(CPACK_PACKAGE_VERSION ${SixTrack_VERSION})
SET(CPACK_PACKAGE_VENDOR "CERN")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SixTrack symplectic particle tracking code for particle accelerators.")
SET(CPACK_PACKAGE_VERSION_MAJOR ${SixTrack_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${SixTrack_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${SixTrack_VERSION_PATCH})
INCLUDE(CPack)

install(TARGETS SixTrackCR RUNTIME DESTINATION bin )
