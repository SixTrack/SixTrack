#Sixtrack CMAKE build file
#Making use of https://cmake.org/Wiki/CMakeForFortranExample
#and the file created for building Merlin.

###################################################################################################
# General configuration
###################################################################################################

#Set the cmake version required
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
#set(CMAKE_LEGACY_CYGWIN_WIN32 0)

# Set the project name and options
PROJECT(SixTrack)
ENABLE_LANGUAGE(Fortran)

file(READ ${CMAKE_SOURCE_DIR}/SixTrack/sixtrack.s IFILE)
STRING(REGEX MATCH ".*data version /'([0-9]+)[.]([0-9]+)[.]([0-9]+)'/" VERSION_STRING ${IFILE} )
set(VERSION_STRING "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
MESSAGE(STATUS "Sixtrack version: ${VERSION_STRING}")
MESSAGE(STATUS "Running on: ${CMAKE_SYSTEM_NAME} with ${CMAKE_SYSTEM_PROCESSOR} and generating for ${CMAKE_GENERATOR}")

#Set Sixtrack version info
set(${PROJECT_NAME}_VERSION "${VERSION_STRING}")
set(${PROJECT_NAME}_VERSION_MAJOR "${CMAKE_MATCH_1}")
set(${PROJECT_NAME}_VERSION_MINOR "${CMAKE_MATCH_2}")
set(${PROJECT_NAME}_VERSION_PATCH "${CMAKE_MATCH_3}")

# force an out of source build to stop clutter
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
	message(STATUS "Sixtrack requires an out-of-source build.")
	message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
	message(STATUS "CMakeCache.txt")
	message(STATUS "CMakeFiles")
	message(STATUS "Once these files are removed, create a separate directory")
	message(STATUS "and run CMake from there")
	message(FATAL_ERROR "in-source build detected")
endif()

if(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"# MPI MPIDEBUG."
		FORCE)
endif()

###################################################################################################
# Sixtrack configuration options
###################################################################################################

# Now we have the default settings
# we have ON
OPTION(TILT "Allow elements to be tilted (by error)" ON)
OPTION(FAST "Which implementation of drifts to use in thin6d" ON)
OPTION(CRLIBM "Use correctly rounded libmath instead of system libmath" ON)
OPTION(CERNLIB "Link to external CERNLIB library for PAW plots. Otherwise use internally defined dummy functions from deck plotdumy" OFF)
OPTION(NAGLIB "Link to external NAGLIB library, used in subroutine coruglo and coruord for DA version" OFF)
# used for ifort only SSE
#SSE="NULL"
# and we have OFF
OPTION(API "Link to external BOINC API (if BOINC but not API, just use internal 'myboinc.f' dummy library)" OFF)
OPTION(DA "Build differential algebra version (NOT dynamic aperture!)" OFF)
OPTION(COLLIMAT "Collimation code" OFF)
OPTION(CPSS "CERN physics screensaver support" OFF)
OPTION(BOINC "BOINC distributed computing support" OFF)
OPTION(CR "Checkpoint / restarting" OFF)
OPTION(BPM "Activates the collection of tracking data for sussix (deprecated)" OFF)
OPTION(BEAMGAS "Beam-gas scattering" OFF)
OPTION(BNLELENS "Support for the BNL electron lens (deprecated)" OFF)
OPTION(BIGNBLZ "Allow up to 200'000 structure elements (normally 20'000 (15k for collimat)) - useful for FCC etc." OFF)
OPTION(HUGENBLZ "Allow up to 400'000 structure elements (normally 20'000 (15k for collimat)) - useful for FCC etc." OFF)
OPTION(STF "Single Track File - write all tracks for postprocessing to fort.90, instead of 32 separate files fort.90 to fort.(90-32). Required for increased particle numbers" OFF)
OPTION(BIGNPART "Allow tracking up to 2048 particles (normally 64)." OFF)
OPTION(HUGENPART "Allow tracking up to 65536 particles (normally 64)" OFF)
OPTION(DATAMODS "Allocate big matrices for thick tracking on the heap and only if required" OFF)
OPTION(DEBUG "Debugging output" OFF)
OPTION(HDF5 "Link to external HDF5 library, write tracks2.dat data (for postprocessing with BeamLossPattern) there instead of writing it to ASCII." OFF)
OPTION(FIO "Use FortranIO from Fortran2003 to correctly round ASCII input/output" OFF)
OPTION(MERLINSCATTER "Interaction physics for collimation from Merlin" OFF)
OPTION(NAFF "Use NAFF library for FMA analysis" OFF)
OPTION(G4COLLIMAT "Interaction physics for collimation from Geant4" OFF)
OPTION(COVERAGE "Enable build flags for testing code coverage with gcov (only works with GNU compilers)" OFF)
OPTION(LIBARCHIVE "Link with LIBARCHIVE, needed for BOINC API and ZIPIT block" OFF)
OPTION(GPROF "Enable build flags for code profiling with gprof (only works with GNU compilers)" OFF)

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

OPTION(BUILD_TESTING "Enable the SixTest suite" OFF)

#OPTION(ROOT "Link to the root libraries" OFF)

SET(BOINC_DIR ${CMAKE_SOURCE_DIR}/SixTrack/boinc CACHE PATH "Path to the boinc build folder")
SET(LIBARCHIVE_BUILD_DIR ${CMAKE_SOURCE_DIR}/SixTrack/libarchive_build CACHE PATH "Path to the libarchive build folder")
SET(LIBARCHIVE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/SixTrack/libarchive CACHE PATH "Path to the libarchive source folder")


###################################################################################################
# Sixtrack binary name construction
###################################################################################################

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

SET(SIXTRACK_FEATURES "cmake")
SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_${CMAKE_SYSTEM_NAME}_${Fortran_COMPILER_NAME}")

if(TILT)
	SET(SIXTRACK_FEATURES "tilt_${SIXTRACK_FEATURES}")
endif()
if(FAST)
	SET(SIXTRACK_FEATURES "fast_${SIXTRACK_FEATURES}")
endif()
if(CRLIBM)
	SET(SIXTRACK_FEATURES "crlibm_${SIXTRACK_FEATURES}")
endif()
if(CERNLIB)
	SET(SIXTRACK_FEATURES "cernlib_${SIXTRACK_FEATURES}")
endif()
if(NAGLIB)
	SET(SIXTRACK_FEATURES "naglib_${SIXTRACK_FEATURES}")
endif()
if(API)
	SET(SIXTRACK_FEATURES "api_${SIXTRACK_FEATURES}")
endif()
if(DA)
	SET(SIXTRACK_FEATURES "da_${SIXTRACK_FEATURES}")
endif()
if(COLLIMAT)
	SET(SIXTRACK_FEATURES "collimat_${SIXTRACK_FEATURES}")
endif()
if(CPSS)
	SET(SIXTRACK_FEATURES "cpss_${SIXTRACK_FEATURES}")
endif()
if(BOINC)
	SET(SIXTRACK_FEATURES "boinc_${SIXTRACK_FEATURES}")
endif()
if(CR)
	SET(SIXTRACK_FEATURES "cr_${SIXTRACK_FEATURES}")
endif()
if(BPM)
	SET(SIXTRACK_FEATURES "bpm_${SIXTRACK_FEATURES}")
endif()
if(BEAMGAS)
	SET(SIXTRACK_FEATURES "beamgas_${SIXTRACK_FEATURES}")
endif()
if(BNLELENS)
	SET(SIXTRACK_FEATURES "bnlelens_${SIXTRACK_FEATURES}")
endif()
if(BIGNBLZ)
	SET(SIXTRACK_FEATURES "bignblz_${SIXTRACK_FEATURES}")
endif()
if(HUGENBLZ)
	SET(SIXTRACK_FEATURES "hugenblz_${SIXTRACK_FEATURES}")
endif()
if(STF)
	SET(SIXTRACK_FEATURES "stf_${SIXTRACK_FEATURES}")
endif()
if(BIGNPART)
	SET(SIXTRACK_FEATURES "bignpart_${SIXTRACK_FEATURES}")
endif()
if(DATAMODS)
	SET(SIXTRACK_FEATURES "datamods_${SIXTRACK_FEATURES}")
endif()
if(HUGENPART)
	SET(SIXTRACK_FEATURES "hugenpart_${SIXTRACK_FEATURES}")
endif()
if(DEBUG)
	SET(SIXTRACK_FEATURES "debug_${SIXTRACK_FEATURES}")
endif()
if(HDF5)
	SET(SIXTRACK_FEATURES "hdf5_${SIXTRACK_FEATURES}")
endif()
if(FIO)
	SET(SIXTRACK_FEATURES "fio_${SIXTRACK_FEATURES}")
endif()
if(MERLINSCATTER)
	SET(SIXTRACK_FEATURES "merlinscatter_${SIXTRACK_FEATURES}")
endif()
if(NAFF)
	SET(SIXTRACK_FEATURES "naff_${SIXTRACK_FEATURES}")
endif()
if(G4COLLIMAT)
	SET(SIXTRACK_FEATURES "g4collimat_${SIXTRACK_FEATURES}")
endif()
if(LIBARCHIVE)
	SET(SIXTRACK_FEATURES "libarchive_${SIXTRACK_FEATURES}")
endif()

#if(ROOT)
#	SET(SIXTRACK_FEATURES "root_${SIXTRACK_FEATURES}")
#endif()
#
#if(OPENMP)
#	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_openmp")
#endif()
#if(MPI)
#	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_mpi")
#endif()
#if(CUDA)
#	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_cuda")
#endif()
if(STATIC)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_static")
endif()
if(AVX)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_avx")
endif(AVX)
if(AVX2)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_avx2")
endif(AVX2)
if(NATIVE)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_native")
endif()

SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_${CMAKE_SYSTEM_PROCESSOR}")

if(32BIT)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_32bit")
endif()
if(64BIT)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_64bit")
endif()

#Create the sixtrack binary name string
SET(SIXTRACK_BINARY_NAME "${PROJECT_NAME}_${${PROJECT_NAME}_VERSION_MAJOR}${${PROJECT_NAME}_VERSION_MINOR}${${PROJECT_NAME}_VERSION_PATCH}_${SIXTRACK_FEATURES}" )

# Select which fortran files to generate and compile
if(DA)
	LIST(APPEND FORS_FIXFORM sixda sixsc sixscfox dabnew lielib)
	LIST(APPEND FORS_FREEFORM scatter dynk)
else()
	LIST(APPEND FORS_FIXFORM track sixve sixvefox dabnews lielib plato_seq)
	LIST(APPEND FORS_FREEFORM beamgas scatter dynk)
endif()

if(NOT DA AND DATAMODS)
	LIST(APPEND FORS_FIXFORM datamods)
endif()

if(COLLIMAT)
	LIST(APPEND FORS_FIXFORM collimation)
endif()

foreach(loop IN LISTS FORS_FIXFORM)
	set(SIXTRACK_SOURCES ${SIXTRACK_SOURCES} ${loop}.f)
endforeach()
foreach(loop IN LISTS FORS_FREEFORM)
	set(SIXTRACK_SOURCES ${SIXTRACK_SOURCES} ${loop}.f90)
endforeach()

#List of files that must be treated by DAFOR
LIST(APPEND FORS_DAFOR sixvefox sixscfox sixda)

if(BOINC AND NOT API)
	#We have a BOINC build but no API, load the fake API substitute for debugging.
	set(SIXTRACK_SOURCES ${SIXTRACK_SOURCES} myboinc.f)
endif()

if(LIBARCHIVE)
	add_library(libArchive_wrapper STATIC  libArchive_Fwrapper.c libArchive_wrapper.c libArchive_wrapper.h)
	target_compile_options(libArchive_wrapper PRIVATE -std=c99)
endif()

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
	##message(STATUS "HDF5_INCLUDE_DIR: ${HDF5_INCLUDE_DIR}")
	message(STATUS "HDF5_DEFINITIONS: ${HDF5_DEFINITIONS}")
	message(STATUS "HDF5_Fortran_LIBRARIES: ${HDF5_Fortran_LIBRARIES}")
	##message(STATUS "HDF5_Fortran_HL_LIBRARIES: ${HDF5_Fortran_HL_LIBRARIES}")
	##message(STATUS "HDF5_HL_LIBRARIES: ${HDF5_HL_LIBRARIES}")
	##message(STATUS "HDF5_C_LIBRARIES: ${HDF5_C_LIBRARIES}")
	##message(STATUS "HDF5_CXX_LIBRARIES: ${HDF5_CXX_LIBRARIES}")
	message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
	##message(STATUS "HDF5_FOUND: ${HDF5_FOUND}")
	##message(STATUS "HDF5_VERSION: ${HDF5_VERSION}")
	message(STATUS "HDF5_LIBRARY_DIRS: ${HDF5_LIBRARY_DIRS}")
	##message(STATUS "HDF5_IS_PARALLEL: ${HDF5_IS_PARALLEL}")
	##message(STATUS "HDF5_Fortran_COMPILER_EXECUTABLE: ${HDF5_Fortran_COMPILER_EXECUTABLE}")
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
	add_subdirectory(${CMAKE_SOURCE_DIR}/SixTrack/g4collimat)
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
#if(ROOT)
#	#/usr/share/root/cmake/FindROOT.cmake
#	SET(CMAKE_MODULE_PATH /usr/share/root/cmake/ ${CMAKE_MODULE_PATH})
#	find_package(ROOT REQUIRED)
#	include_directories(${ROOT_INCLUDE_DIR})
#
#	if(NOT ROOT_FOUND)
#		MESSAGE(FATAL_ERROR "Root enabled build requested but root was not found!")
#	endif()
#endif(ROOT)

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

if(CRLIBM)
	if(COLLIMAT)
		message(FATAL_ERROR "Cannot build CRLIBM with COLLIMAT enabled..")
	endif()
endif()
if(HDF5)
	if(NOT COLLIMAT)
		message(FATAL_ERROR "Building with HDF5 requires COLLIMAT enabled")
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
	if(BIGNBLZ)
		message(FATAL_ERROR "Cannot build BEAMGAS with BIGNBLZ enabled.")
	endif(BIGNBLZ)
	if(HUGENBLZ)
		message(FATAL_ERROR "Cannot build BEAMGAS with HUGENBLZ enabled.")
	endif(HUGENBLZ)
endif(BEAMGAS)
if(BIGNBLZ)
	if(HUGENBLZ)
		message(FATAL_ERROR "Cannot combine BIGNBLZ with HUGENBLZ.")
	endif()
endif()
if(HUGENPART)
	if(BIGNPART)
		message(FATAL_ERROR "Cannot combine BIGNPART with HUGENPART.")
	endif()
	if(BNLELENS)
		message(FATAL_ERROR "Cannot build HUGENPART with BNLELENS enabled.")
	endif()
	if(NOT STF)
		message(FATAL_ERROR "Building with HUGENPART requires STF enabled.")
	endif()
	if(NOT DATAMODS)
		message(FATAL_ERROR "Building with HUGENPART requires DATAMODS enabled.")
	endif()
endif()
if(DA)
	if(COLLIMAT)
		message(FATAL_ERROR "Cannot build DA with COLLIMAT enabled.")
	endif()
	if(CPSS)
		message(FATAL_ERROR "Cannot build DA with CPSS enabled.")
	endif()
	if(BPM)
		message(FATAL_ERROR "Cannot build DA with BPM enabled.")
	endif()
	if(NOT NAGLIB)
		message(FATAL_ERROR "Building with DA requires NAGLIB enabled.")
	endif()
endif()
if(COLLIMAT)
	if(DA)
		message(FATAL_ERROR "Cannot build COLLIMAT with DA enabled.")
	endif()
	if(CPSS)
		message(FATAL_ERROR "Cannot build COLLIMAT with CPSS enabled.")
	endif()
	if(BPM)
		message(FATAL_ERROR "Cannot build COLLIMAT with BPM enabled.")
	endif()
	if(CR)
		message(FATAL_ERROR "Cannot build COLLIMAT with CR enabled.")
	endif()
	if(CRLIBM)
		message(FATAL_ERROR "Cannot build COLLIMAT with CRLIBM enabled.")
	endif()
endif()
if(CPSS)
	if(CERNLIB)
		message(FATAL_ERROR "Cannot build CPSS with CERNLIB enabled.")
	endif()
	if(NOT CRLIBM)
		message(FATAL_ERROR "Building with CPSS requires CRLIBM enabled.")
	endif()
	if(NOT CR)
		message(FATAL_ERROR "Building with CPSS requires CR enabled.")
	endif()
endif()
if(STF)
	if(BNLELENS)
		message(FATAL_ERROR "Cannot build STF with BNLELENS enabled.")
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
if(BIGNPART)
	if(COLLIMAT)
		message(FATAL_ERROR "Cannot build BIGNPART with COLLIMAT enabled.")
	endif()
	if(BNLELENS)
		message(FATAL_ERROR "Cannot build BIGNPART with BNLELENS enabled.")
	endif()
	if(NOT STF)
		message(FATAL_ERROR "Building with BIGNPART requires STF enabled.")
	endif()
endif()
if(64BIT)
	if(CERNLIB)
		message(FATAL_ERROR "Cannot build CERNLIB with 64BIT enabled.")
	endif()
endif()
if(AVX2)
	if(NOT AVX)
		message(FATAL_ERROR "Building with AVX2 requires AVX enabled")
	endif(NOT AVX)
endif(AVX2)

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
	# gfortran
	SET(CMAKE_Fortran_FLAGS "-frecord-marker=4 -fno-second-underscore -funroll-loops")

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
			set(CMAKE_Fortran_FLAGS "-mavx2 -mno-fma4 -mno-fma ${CMAKE_Fortran_FLAGS}")
		endif(AVX2)

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
	# ifort
	set (CMAKE_Fortran_FLAGS_RELEASE "-fp-model source -fp-model strict -fp-model no-except -O3 -no-fma")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-fp-model source -fp-model strict -fp-model no-except -O0 -g -no-fma")

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
		set(CMAKE_Fortran_FLAGS "-mavx2 -no-fma ${CMAKE_Fortran_FLAGS}")
	endif(AVX2)

	if(NATIVE)
		SET(CMAKE_Fortran_FLAGS_RELEASE "-xHost ${CMAKE_Fortran_FLAGS_RELEASE}")
		#Nothing to do for DEBUG
	endif(NATIVE)

	if(STATIC)
		SET(CMAKE_Fortran_LINKER_FLAGS "-static ${CMAKE_Fortran_LINKER_FLAGS}")
	endif(STATIC)

elseif (Fortran_COMPILER_NAME MATCHES "nagfor.*")
	# nagfor
	SET(NAGFOR ON)
	#Here I am going to assume that since we are using nagfor we are on lxplus for now
	set (CMAKE_Fortran_FLAGS_RELEASE "-O4 -dusty -ieee=full -dcfuns -w=ulv -w=uparam")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -dusty -ieee=full -dcfuns -w=ulv -w=uparam")

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
			set(CMAKE_C_FLAGS "-mavx2 -mno-fma4 -mno-fma ${CMAKE_C_FLAGS}")
		endif(AVX2)

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
		set(CMAKE_C_FLAGS "-mavx2 -no-fma ${CMAKE_C_FLAGS}")
	endif(AVX2)

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
			set(CMAKE_CXX_FLAGS "-mavx2 -mno-fma4 -mno-fma ${CMAKE_CXX_FLAGS}")
		endif(AVX2)

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
		set(CMAKE_CXX_FLAGS "-mavx2 -no-fma ${CMAKE_CXX_FLAGS}")
	endif(AVX2)

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
#add the source files that will be used to build sixtrack
add_executable(SixTrack ${SIXTRACK_SOURCES})
SET_TARGET_PROPERTIES(SixTrack PROPERTIES OUTPUT_NAME "${SIXTRACK_BINARY_NAME}" LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}" LINKER_LANGUAGE Fortran)

##build astuce
add_executable(astuce_old ${CMAKE_SOURCE_DIR}/SixTrack/tools/astuce/astuce.f)
SET_TARGET_PROPERTIES(astuce_old PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}" OUTPUT_NAME "astuce")

add_library(astuce++ STATIC ${CMAKE_SOURCE_DIR}/SixTrack/tools/astuce/astuce.cpp)
add_executable(astuce ${CMAKE_SOURCE_DIR}/SixTrack/tools/astuce/astuce++.cpp)
SET_TARGET_PROPERTIES(astuce PROPERTIES LINK_FLAGS "${CMAKE_CXX_LINKER_FLAGS}" OUTPUT_NAME "astuce++")
target_link_libraries(astuce astuce++)

#build dafor
add_executable(dafor ${CMAKE_SOURCE_DIR}/SixTrack/tools/dafor/dafor.f)
SET_TARGET_PROPERTIES(dafor PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")

if(CRLIBM)
  add_subdirectory(${CMAKE_SOURCE_DIR}/SixTrack/crlibm)
endif()
if(CRLIBM OR BUILD_TESTING )
  add_subdirectory(${CMAKE_SOURCE_DIR}/SixTrack/roundctl) #Also needed by read90
endif()

if(HDF5)
	if(HDF5_FOUND)
		SET_TARGET_PROPERTIES(SixTrack PROPERTIES COMPILE_FLAGS "-I${HDF5_INCLUDE_DIRS} ${HDF5_DEFINITIONS}")
	else()
		MESSAGE(FATAL_ERROR "HDF5 build requested, but could not find the HDF5 package.")
	endif(HDF5_FOUND)
endif(HDF5)

if(MERLINSCATTER)
	file(GLOB MERLINSCATTER_SOURCES ${CMAKE_SOURCE_DIR}/SixTrack/merlinscatter/*.cpp ${CMAKE_SOURCE_DIR}/SixTrack/merlinscatter/RandomNG/*.cpp)
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
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/SixTrack/cmake_modules/")
	if (STATIC)
		SET(FFTW_USE_STATIC_LIBS 1)
	endif()
	find_package(FFTW REQUIRED)
	include_directories( ${FFTW_INCLUDES} )
	message(STATUS "Include dirs of FFTW: " ${FFTW_INCLUDES} )
	message(STATUS "Libs of of FFTW:      " ${FFTW_LIBRARIES} )
	
	
	file(GLOB NAFF_SOURCES ${CMAKE_SOURCE_DIR}/SixTrack/naff/*.cpp ${CMAKE_SOURCE_DIR}/SixTrack/naff/*.cc)
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
add_dependencies(copy_input_files ast_mkdir)
file(GLOB SIXTRACK_INPUT ${CMAKE_SOURCE_DIR}/SixTrack/*.s ${CMAKE_SOURCE_DIR}/SixTrack/*.s90)

#copy all the input files *.s to our build folder
foreach(loop IN ITEMS ${SIXTRACK_INPUT})
	add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${loop} ${CMAKE_BINARY_DIR})
endforeach()

if(BOINC)
	if(API)
		include_directories(${BOINC_DIR}/api ${BOINC_DIR}/lib ${BOINC_DIR})
		
		#The fortran interface needs boinc_api_fortran.o building, which does not seem to happen by default
		#go into the boinc source folder, cd api; make boinc_api_fortran.o
		add_library(boinc_api_fortran STATIC IMPORTED)
		set_target_properties(boinc_api_fortran PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/api/boinc_api_fortran.o)
		# add_library(boinc_api_fortran STATIC ${BOINC_DIR}/api/boinc_api_fortran.cpp) #Alternative way of building and linking boinc_api_fortran.cpp/.o
		target_link_libraries(SixTrack boinc_api_fortran)
		
		#add the libs
		add_library(boinc STATIC IMPORTED)
		set_target_properties(boinc PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/lib/libboinc.a)
		
		add_library(boinc_api STATIC IMPORTED)
		set_target_properties(boinc_api PROPERTIES IMPORTED_LOCATION ${BOINC_DIR}/api/libboinc_api.a)
		target_link_libraries(SixTrack boinc_api)
		
		target_link_libraries(SixTrack boinc)
		
		find_library(stdc++ STATIC stdc++)
		find_library(pthread STATIC pthread)
		target_link_libraries(SixTrack stdc++ pthread)

		if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
			find_library(nsl STATIC nsl)
			target_link_libraries(SixTrack nsl)
		endif()

		#fix osx linking to the correct c++ stdlib
		if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
			find_library(c++ STATIC c++)
			target_link_libraries(SixTrack c++)
		endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

	else()
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/myboinc.f POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/SixTrack/myboinc.f ${CMAKE_BINARY_DIR}/myboinc.f DEPENDS copy_input_files)
	endif()
endif()

if(LIBARCHIVE)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libArchive_Fwrapper.c POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/SixTrack/libArchive_Fwrapper.c ${CMAKE_BINARY_DIR}/libArchive_Fwrapper.c DEPENDS copy_input_files)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libArchive_wrapper.c POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/SixTrack/libArchive_wrapper.c ${CMAKE_BINARY_DIR}/libArchive_wrapper.c DEPENDS copy_input_files)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libArchive_wrapper.h POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/SixTrack/libArchive_wrapper.h ${CMAKE_BINARY_DIR}/libArchive_wrapper.h DEPENDS copy_input_files)
	
	include_directories(${LIBARCHIVE_SOURCE_DIR}/libarchive)
	include_directories(${CMAKE_CURRENT_BINARY_DIR})

	#A very small but convenient test code for debugging libArchive.
	add_executable(libArchiveTester ${CMAKE_SOURCE_DIR}/SixTrack/libArchiveTester.cpp)
	target_link_libraries(libArchiveTester libArchive_wrapper)

	#See also further down for SixTestWrapper linking
	target_link_libraries(SixTrack libArchive_wrapper)
	add_library(libarchive STATIC IMPORTED)
	if(WIN32)
		set_target_properties(libarchive PROPERTIES IMPORTED_LOCATION ${LIBARCHIVE_BUILD_DIR}/libarchive/libarchive_static.a)
		add_definitions( -DLIBARCHIVE_STATIC )
	else()
		set_target_properties(libarchive PROPERTIES IMPORTED_LOCATION ${LIBARCHIVE_BUILD_DIR}/libarchive/libarchive.a)
	endif()
	
	find_library(pthread STATIC pthread)
	#find_library(z STATIC z)
	add_library(z STATIC IMPORTED)
	set_target_properties(z PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/zlib/install/lib/libz.a)
	
	target_link_libraries(SixTrack libarchive z pthread)
	target_link_libraries(libArchiveTester libarchive z pthread)
	
	add_definitions( -DLIBARCHIVE )
endif()

#make a folder for the astuce masks and copy the input files to this folder
add_custom_target(ast_mkdir)
file(GLOB ASTUCE_MASKS ${CMAKE_SOURCE_DIR}/SixTrack/ast_mask/*.ast)
add_custom_command(TARGET ast_mkdir POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/ast_mask/)
foreach(loop IN ITEMS ${ASTUCE_MASKS})
	add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${loop} ${CMAKE_BINARY_DIR}/ast_mask/ DEPENDS ast_mkdir)
endforeach()

#We need to build asatuce and dafor before we build SixTrack, but crlibm and so on can be built at the same time.
add_dependencies(SixTrack copy_input_files generate_input)

###################################################################################################
#running astuce and dafor
###################################################################################################
add_custom_target(run_astuce)

#make running astuce depend on astuce being built and the input files being in place
add_dependencies(run_astuce astuce astuce_clean_define_flags astuce_clean_decks astuce_old)

#This target cleans the define flags
add_custom_target(astuce_clean_define_flags)
#this requires the input files to be in place
add_dependencies(astuce_clean_define_flags copy_input_files)

#this target cleans the decks
add_custom_target(astuce_clean_decks)
#and requires the define flags to already be cleaned (this is just to stop parallel writes on the same file)
add_dependencies(astuce_clean_decks astuce_clean_define_flags)

#This macro calls the regex in masks_define_flags.cmake to remove the requested flags
macro(ast_masks_flag OPT FILE_NAME)
	add_custom_command(TARGET astuce_clean_define_flags POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -DFILE_NAME:FILEPATH=${CMAKE_BINARY_DIR}/ast_mask/${FILE_NAME}.ast -DOPT:STRING=${OPT} -P ${CMAKE_SOURCE_DIR}/SixTrack/masks_define_flags.cmake)
endmacro()

#This macro calls the regex in masks_decks.cmake to remove the requested decks
macro(ast_masks_deck OPT FILE_NAME)
	add_custom_command(TARGET astuce_clean_decks POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -DFILE_NAME:FILEPATH=${CMAKE_BINARY_DIR}/ast_mask/${FILE_NAME}.ast -DOPT:STRING=${OPT} -P ${CMAKE_SOURCE_DIR}/SixTrack/masks_decks.cmake)
endmacro()

foreach(loop IN LISTS FORS_FIXFORM FORS_FREEFORM)
	#Clean up the astuce flags for only what is enabled
	if(NOT CRLIBM)
		ast_masks_flag(crlibm ${loop})
	endif()
	if(NOT TILT)
		ast_masks_flag(tilt ${loop})
	endif()
	if(NOT FAST)
		ast_masks_flag(fast ${loop})
	endif()
	if(NOT BPM)
		ast_masks_flag(bpm ${loop})
	endif()
	if(NOT BEAMGAS)
		ast_masks_flag(beamgas ${loop})
	endif()
	if(NOT COLLIMAT)
		ast_masks_flag(collimat ${loop})
	endif()
	if(NOT BNLELENS)
		ast_masks_flag(bnlelens ${loop})
	endif()
	if(NOT CR)
		ast_masks_flag(cr ${loop})
	endif()
	if(NOT BOINC)
		ast_masks_flag(boinc ${loop})
	endif()
	if(NOT NAGFOR)
		ast_masks_flag(nagfor ${loop})
	endif()
	if(NOT BIGNBLZ)
		ast_masks_flag(bignblz ${loop})
	endif()
	if(NOT DEBUG)
		ast_masks_flag(debug ${loop})
	endif()
	if(NOT HDF5)
		ast_masks_flag(hdf5 ${loop})
	endif()
	if(NOT FIO)
		ast_masks_flag(fio ${loop})
	endif()
	if(NOT LF95)
		ast_masks_flag(lf95 ${loop})
	endif()
	if(NOT HUGENBLZ)
		ast_masks_flag(hugenblz ${loop})
	endif()
	if(NOT STF)
		ast_masks_flag(stf ${loop})
	endif()
	if(NOT BIGNPART)
		ast_masks_flag(bignpart ${loop})
	endif()
	if(NOT DATAMODS)
		ast_masks_flag(datamods ${loop})
	endif()
	if(NOT HUGENPART)
		ast_masks_flag(hugenpart ${loop})
	endif()
	if(NOT MERLINSCATTER)
		ast_masks_flag(merlinscatter ${loop})
	endif()
	if(NOT NAFF)
		ast_masks_flag(naff ${loop})
	endif()

	if(NOT LIBARCHIVE)
		ast_masks_flag(libarchive ${loop})
	endif()
	if(NOT G4COLLIMAT)
		ast_masks_flag(g4collimat ${loop})
	endif()

	#Clean up the astuce decks for only what is enabled
	if(NOT COLLIMAT)
		ast_masks_deck(nwrtcoll ${loop})
	endif()
	if(NOT DATAMODS)
		ast_masks_deck(datamods ${loop})
	endif()
	if(NOT DEBUG)
		ast_masks_deck(dumps ${loop})
	endif()
	if(NOT CR)
		ast_masks_deck(checkpt ${loop})
	endif()
	if(NOT BNLELENS)
		ast_masks_deck(nwrtbnl ${loop})
	endif()

	#From make_six:
	#NOTA BENE! Definitions of cernlib and naglib are reversed
	# in the sense that while a "" null value implies option is ON
	# the reverse is true for cernlib and naglib
	if(CERNLIB)
		ast_masks_deck(plotdumy ${loop})
	endif()
	if(NAGLIB)
		ast_masks_deck(nagdumy ${loop})
	endif()
endforeach(loop)

add_custom_target(generate_input)
add_dependencies(generate_input run_astuce dafor)

## Run ASTUCE and DAFOR to generate Fortran source files ##
foreach(loop IN LISTS FORS_FIXFORM FORS_FREEFORM)
	message(STATUS "Will run astuce to generate ${loop}n.f using ${loop}.ast")
	add_custom_command(TARGET run_astuce POST_BUILD COMMAND $<TARGET_FILE:astuce> ${CMAKE_CURRENT_BINARY_DIR}/ast_mask/${loop}.ast)
	
	LIST(FIND FORS_DAFOR   ${loop}  FORS_DAFOR_LISTINDEX )
	LIST(FIND FORS_FIXFORM ${loop}  FORS_FIXFORM_LISTINDEX )
	LIST(FIND FORS_FREEFORM ${loop} FORS_FREEFORM_LISTINDEX )
	
	if(${FORS_DAFOR_LISTINDEX} GREATER -1)
		message(STATUS "Will run dafor to generate ${loop}.f from ${loop}n.f")
		if(NOT ${FORS_FIXFORM_LISTINDEX} GREATER -1)
			message(FATAL_ERROR "Was configured to run dafor on ${loop}, but this is not a fixform file!")
		endif()
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f COMMAND $<TARGET_FILE:dafor> ${loop}n.f ${loop}.f > ${loop}_dafor.log DEPENDS run_astuce)
		
	elseif(${FORS_FIXFORM_LISTINDEX} GREATER -1)
		message(STATUS "Will copy ${loop}n.f to ${loop}.f")
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f COMMAND ${CMAKE_COMMAND} -E copy ${loop}n.f ${loop}.f DEPENDS run_astuce)
		
	elseif(${FORS_FREEFORM_LISTINDEX} GREATER -1)
		message(STATUS "Will copy ${loop}n.f to ${loop}.f90")
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f90 COMMAND ${CMAKE_COMMAND} -E copy ${loop}n.f ${loop}.f90 DEPENDS run_astuce)
		
	else()
		message(FATAL_ERROR "Error when preprocessing ${loop}: Logic error in CMakeLists.txt.")
	endif()
endforeach()

###################################################################################################
#linking
###################################################################################################
set(MSYS_QUADMATH "/usr/lib/gcc/x86_64-pc-msys/6.3.0/libquadmath.a")
#just a work around to build on MSYS2 for now
if(MSYS AND STATIC)
	target_link_libraries(astuce_old ${MSYS_QUADMATH})
	target_link_libraries(dafor ${MSYS_QUADMATH})
	target_link_libraries(SixTrack ${MSYS_QUADMATH})
endif()

#same with OSX, using gfortran from homebrew
if(STATIC)
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		find_library(lquadmath NAMES libquadmath.a PATHS /usr/local/opt/gcc/lib/gcc/6/ /usr/local/opt/gcc/lib/gcc/7/ /usr/local/opt/gcc/lib/gcc/8/)
		target_link_libraries(SixTrack ${lquadmath})
	endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
endif(STATIC)

if(CERNLIB)
#	if(EXISTS "/cern/pro/bin/cernlib")
		#running on lxplus and friends
#		SET_TARGET_PROPERTIES(SixTrack PROPERTIES LINK_FLAGS "`cernlib -u graflib packlib_noshift kernlib_noshift` -lg2c -lpthread")
#	else()
                #Use precompiled CERNLIB from lib32 folder
		add_library(g2c STATIC IMPORTED)
		add_library(graflib STATIC IMPORTED)
		add_library(grafX11 STATIC IMPORTED)
		add_library(kernlib_noshift STATIC IMPORTED)
		add_library(packlib_noshift STATIC IMPORTED)
		add_library(X11 STATIC IMPORTED)
		add_library(Xau STATIC IMPORTED)
		add_library(xcb STATIC IMPORTED)
		add_library(Xdmcp STATIC IMPORTED)

		set_target_properties(g2c PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libg2c.a)
		set_target_properties(graflib PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libgraflib.a)
		set_target_properties(grafX11 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libgrafX11.a)
		set_target_properties(kernlib_noshift PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libkernlib_noshift.a)
		set_target_properties(packlib_noshift PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libpacklib_noshift.a)
		set_target_properties(X11 PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libX11.a)
		set_target_properties(Xau PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libXau.a)
		set_target_properties(xcb PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libxcb.a)
		set_target_properties(Xdmcp PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTrack/lib32/libXdmcp.a)

		find_library(pthread STATIC pthread)
		find_library(dl STATIC dl)
		find_library(gcc_eh STATIC gcc_eh)

		target_link_libraries(SixTrack graflib grafX11 packlib_noshift kernlib_noshift g2c X11 xcb Xau Xdmcp dl pthread gcc_eh)
#	endif()
endif(CERNLIB)

#if we are using crlibm, link to the crlibm library
if(CRLIBM)
	target_link_libraries(SixTrack crlibm roundctl)
endif(CRLIBM)

#if we are using merlinscatter, link to the merlinscatter library
if(MERLINSCATTER)
	target_link_libraries(SixTrack merlinscatter)
endif(MERLINSCATTER)
#if we are using naff, link to the naff library
if(NAFF)
	target_link_libraries(SixTrack naff)
endif(NAFF)


#if we are using root, link to all the root libraries
#if(ROOT)
#	target_link_libraries(SixTrack ${ROOT_LIBRARIES})
#endif(ROOT)

if(HDF5)
	target_link_libraries(SixTrack ${HDF5_LIBRARIES})
endif(HDF5)

if(G4COLLIMAT)
	find_library(pthread STATIC pthread)
	target_link_libraries(SixTrack g4collimat ${Geant4_LIBRARIES} pthread)
endif(G4COLLIMAT)

#Again for naglib, I assume we are on lxplus or have CERN afs access.
if(NAGLIB)
	if(64BIT)
		if(STATIC)
			find_library(nglibifcoremt NAMES libifcoremt.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglibsvml NAMES libsvml.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglibirng NAMES libirng.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglibimf NAMES libimf.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglibirc NAMES libirc.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglib NAMES libnag_nag.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/lib)
		else()
			find_library(nglibifcoremt NAMES libifcoremt.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglibsvml NAMES libsvml.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglibirng NAMES libirng.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglibimf NAMES libimf.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglibirc NAMES libirc.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/rtl)
			find_library(nglib NAMES libnag_nag.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fll6i24dcl/lib)
		endif(STATIC)
	endif(64BIT)

	if(32BIT)
		if(STATIC)
			find_library(nglibifcoremt NAMES libifcoremt.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglibsvml NAMES libsvml.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglibirng NAMES libirng.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglibimf NAMES libimf.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglibirc NAMES libirc.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglib NAMES libnag_nag.a PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/lib)
		else()
			find_library(nglibifcoremt NAMES libifcoremt.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglibsvml NAMES libsvml.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglibirng NAMES libirng.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglibimf NAMES libimf.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglibirc NAMES libirc.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/rtl)
			find_library(nglib NAMES libnag_nag.so PATHS /afs/cern.ch/sw/nag/mark24/lnx/fllux24dcl/lib)
		endif(STATIC)
	endif(32BIT)

	target_link_libraries(SixTrack ${nglib} ${nglibifcoremt} ${nglibsvml} ${nglibirng} ${nglibimf} ${nglibirc} dl pthread )
endif(NAGLIB)

IF(COVERAGE)
	target_compile_options(SixTrack PRIVATE -fprofile-arcs -ftest-coverage)
	SET_TARGET_PROPERTIES(SixTrack PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
ENDIF(COVERAGE)

SET(BUILDNAME "${SIXTRACK_BINARY_NAME}" CACHE STRING "Name of build on the dashboard")
MARK_AS_ADVANCED(BUILDNAME)

#For testing
IF(BUILD_TESTING)

	ENABLE_TESTING()
	INCLUDE(CTest)

	if(DA)
		add_subdirectory(${CMAKE_SOURCE_DIR}/SixTest_da)
	else()
		add_subdirectory(${CMAKE_SOURCE_DIR}/SixTest)
	endif()

	configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake  ${CMAKE_BINARY_DIR})

	#Checks of specialized input/canoncials
	# Anyway hardcoded to use them from the SixTest source folder...
	#add_executable(CheckTestInputForBOINC.sh IMPORTED)
	#SET_TARGET_PROPERTIES(CheckTestInputForBOINC.sh PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTest/CheckTestInputForBOINC.sh)

	#add_executable(CheckTestInputForSTF.sh IMPORTED)
	#SET_TARGET_PROPERTIES(CheckTestInputForSTF.sh PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/SixTest/CheckTestInputForBOINC.sh)

	#need to build the testing tools
	add_executable(verify10 ${CMAKE_SOURCE_DIR}/SixTest/bin/verify10.f)
	add_executable(checkf10 ${CMAKE_SOURCE_DIR}/SixTest/bin/checkf10.f)
	add_executable(compf10 ${CMAKE_SOURCE_DIR}/SixTest/bin/compf10.f)
	add_executable(checkf110 ${CMAKE_SOURCE_DIR}/SixTest/bin/checkf110.f)
	add_executable(checkf1014 ${CMAKE_SOURCE_DIR}/SixTest/bin/checkf1014.f)
	add_executable(read90 ${CMAKE_SOURCE_DIR}/SixTest/read90/read90.f90)
	add_executable(readDump3 ${CMAKE_SOURCE_DIR}/SixTest/readDump3/readDump3.f90)

	target_link_libraries(read90 roundctl)

	SET_TARGET_PROPERTIES(verify10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(checkf10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(compf10 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(checkf110 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(checkf1014 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(read90 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")
	SET_TARGET_PROPERTIES(readDump3 PROPERTIES LINK_FLAGS "${CMAKE_Fortran_LINKER_FLAGS}")

	#The testing wrapper
	add_executable(SixTestWrapper ${CMAKE_SOURCE_DIR}/SixTest/bin/SixTestWrapper.cpp)
	if (LIBARCHIVE)
		target_link_libraries(SixTestWrapper libArchive_wrapper libarchive z pthread)
	endif()
	SET_TARGET_PROPERTIES(SixTestWrapper PROPERTIES LINK_FLAGS "${CMAKE_CXX_LINKER_FLAGS}")

	find_library(pthread STATIC pthread)
	target_link_libraries(SixTestWrapper pthread)

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

install(TARGETS SixTrack RUNTIME DESTINATION bin )
