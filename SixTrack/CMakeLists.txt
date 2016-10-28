#Sixtrack CMAKE build file
#Making use of https://cmake.org/Wiki/CMakeForFortranExample
#and the file created for building Merlin.

#!/bin/bash
#
# Version 1.2 29/08/2005
# CERNLIB/GRAFLIB re-defined
# Forces lf95 for crlibm
# Adds g95 option
# Version 1.3 05/09/2005
# CERNLIB only for SixTrack (not SixTrack_da); CERNGRAF dropped.
# The -naglib option now works for SixTrack_da
# Version 1.4 1/8/2008
# Gets the version and modification date from sixtrack.s
# Corrects the BOINC variable in the Makefile
# make_six clean now exits after the clean just like a make clean
# Version 1.5 11/07/2008
# Adds the deck nwrtbnl and the flags bpm and bnlelens
# Version 2.0 Major change in philosophy; to clarify and document
# all types of make.
# Version 2.1 Adds bignblz option for Frank.
#
# July, 2009. Updated all libraries etc to use Frank's cernlib and X11 libraries
# if they are not present (SLC5!). Also verified correct linking on 32/64 bit
# systems, SLC4 and SLC5.    A NIGHTMARE!!!
# Default: tilt tracking fast crlibm cernlib naglib -bpm -api -da -collimat -bnlelens -cpss -boinc -cr
# and uses the Lahey-Fujitsu lf95 compiler 
# 12th August, 2009. Just changed print_map to print-map. Eric.
# Options:
# da => -tracking -cpss -boinc -collimat -cr +naglib
# collimat => -da -cr -boinc -cpss
# cpss => +cr -collimat -da -cernlib -naglib -boinc
# boinc => +cr +crlibm -collimat -da -cernlib -naglib -cpss
# bpm activates the collection of tracking data for sussix => basic SixTrack
# bnlelens activates the bnl electron lens options possibly => boinc, cr, crlibm
# nag => use f95 instead of LAHEY
# pgf90 => use pgf90
# g95 => use g95
# gfortran => use gfortran
# Version 2.2 August 2010
# Changed to nagfor, changed NAGLIBC for SLC4/5 compatibility
# Version 2.3 27th June, 2011 (including ylevinse collimation mods
# creating a stand alone Makefile and with an hdf5 flag).
# Replaced windows flag by api flag for BOINC
# Added Intel ifort, fio (rounded formatted I/O) flags,
# and an optional optimisation parameter [O|o][0-4].
# A new deck dump for debug. Changed nag to nagfor everywhere.
# Allowed crlibm for NAG nagfor and Intel ifort.
# Version 2.4 January, 2012. merged latest SixTrack production Rev 119
# with all my changes for ifort and FIO.
# Version 2.5 June, 2012.
# Added option SSE2 to try and fix HORRIBLE Intel/AMD performance issues.
# and many other changes for Windows and BOINC including templates.
# Version 2.6. 4th July
# uname -s for Darwin SSE as sse2 or sse3
# and added SSE4.1 and SSE4.2 (just in case)
# Note that SSE$.1 and SSE4.2 are specified as sse41 sse42 (aesthetics)

#Possible debug options, uncomment to use
#set -e #Stop on error - This is sometimes problematic due to the use of $? in the script.
#set -v #Verbose: Print input lines as they are read from the shell script file
#set -x #XTrace:  Print command traces before executing

###################################################################################################
# General configuration
###################################################################################################

#Set the cmake version required
cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)

# Set the project name and options
PROJECT(SixTrack)
ENABLE_LANGUAGE(Fortran)

file(READ ${CMAKE_SOURCE_DIR}/sixtrack.s IFILE)
STRING(REGEX MATCH ".*data version /'([0-9]+)[.]([0-9]+)[.]([0-9]+)'/" VERSION_STRING ${IFILE} )
set(VERSION_STRING "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
MESSAGE(STATUS "Sixtrack version: ${VERSION_STRING}")

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

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs to include the MPI options
#SET(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
SET(CMAKE_BUILD_TYPE Release CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel MPI MPIDEBUG."
	FORCE)

###################################################################################################
# Sixtrack configuration options
###################################################################################################

# Now we have the default settings
# we have ON
OPTION(TILT "??" ON)
OPTION(FAST "??" ON)
OPTION(CRLIBM "Constant rounding library support" ON)
OPTION(CERNLIB "??" ON)
OPTION(NAGLIB "??" ON)
# used for ifort only SSE
#SSE="NULL"
# and we have OFF
OPTION(API "API" OFF)
OPTION(DA "Dynamic aperture code" OFF)
OPTION(COLLIMAT "Collimation code" OFF)
OPTION(CPSS "???" OFF)
OPTION(BOINC "boinc distributed computing support" OFF)
OPTION(CR "??" OFF)
OPTION(NAGFOR "??" OFF)
OPTION(BPM "Activates the collection of tracking data for sussix" OFF)
OPTION(BEAMGAS "Beamgas" OFF)
OPTION(BNLELENS "Support for the BNL electron lens" OFF)
OPTION(BIGNBLZ "??" OFF)
OPTION(HUGENBLZ "??" OFF)
OPTION(STF "??" OFF)
OPTION(BIGNPART "Allow more particles" OFF)
OPTION(DATAMODS "??" OFF)
OPTION(HUGENPART "Allow lots of particles" OFF)
OPTION(DEBUG "Debugging" OFF)
OPTION(HDF5 "??" OFF)
OPTION(FIO "rounded formatted I/O flags" OFF)

#Compiler stuff
OPTION(IFORT "Use the intel fortran compiler" OFF)
OPTION(LF95 "Use the Lahey-Fujitsu lf95 fortran compiler" OFF)
OPTION(PGF90 "Use the pgf90 fortran compiler" OFF)
OPTION(F95 "Use the g95 fortran compiler" OFF)
OPTION(GFORTRAN "Use the gfortran fortran compiler" OFF)

#Build options
OPTION(OPENMP "Enable OpenMP support" OFF)
OPTION(MPI "Enable MPI support" OFF)
OPTION(CUDA "Enable CUDA support" OFF)
OPTION(32BIT "Create a 32bit binary" ON)
OPTION(64BIT "Create a 64bit binary" OFF)

OPTION(ROOT "Link to the root libraries" OFF)

###################################################################################################
# Sixtrack binary name construction
###################################################################################################
SET(SIXTRACK_FEATURES "cmake")
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
if(NAGFOR)
	SET(SIXTRACK_FEATURES "nagfor_${SIXTRACK_FEATURES}")
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

if(ROOT)
	SET(SIXTRACK_FEATURES "root_${SIXTRACK_FEATURES}")
endif()

if(OPENMP)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_openmp")
endif()
if(MPI)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_mpi")
endif()
if(CUDA)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_cuda")
endif()
if(32BIT)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_32")
endif()
if(64BIT)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_64")
endif()

#Ugly, but does the job...
SET(SIXTRACK_BINARY_NAME "${PROJECT_NAME}_${${PROJECT_NAME}_VERSION_MAJOR}${${PROJECT_NAME}_VERSION_MINOR}${${PROJECT_NAME}_VERSION_PATCH}_${SIXTRACK_FEATURES}" )
if(DA)
    SET(FORS sixda sixsc sixscfox dabnew lielib beamgas)
else()
    SET(FORS track sixve sixvefox dabnews lielib beamgas plato_seq datamods)
endif()

foreach(loop IN ITEMS ${FORS})
	set(SIXTRACK_SOURCES ${SIXTRACK_SOURCES} ${loop}.f)
endforeach()

#add the source files that will be used to build sixtrack
add_executable(${SIXTRACK_BINARY_NAME} ${SIXTRACK_SOURCES})

###################################################################################################
# Checking if some requested packages exist
###################################################################################################

#Check for OpenMP
if(OPENMP)
	find_package(OpenMP REQUIRED)
endif(OPENMP)

#If we have OpenMP enabled, set the compiler flags needed
if(OPENMP_FOUND)
	SET_TARGET_PROPERTIES(${SIXTRACK_BINARY_NAME} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
	SET_TARGET_PROPERTIES(${SIXTRACK_BINARY_NAME} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif()

#Possible future  MPI test
if(CMAKE_BUILD_TYPE STREQUAL MPI OR CMAKE_BUILD_TYPE STREQUAL MPIDEBUG)
	FIND_PACKAGE(MPI REQUIRED)
	if(NOT MPI_CXX_FOUND)
		MESSAGE(FATAL_ERROR "MPI build requested but no MPI libraries found!")
	endif()
	if(MPI_CXX_FOUND)
		message(STATUS "Using MPI C++ compiler: ${MPI_CXX_COMPILER}")
		message(STATUS "Using MPI runtime: ${MPIEXEC}")
		message(STATUS "Using MPI process number switch: ${MPIEXEC_NUMPROC_FLAG}")
		message(STATUS "Using MPI pre-flags: ${MPIEXEC_PREFLAGS}")
		message(STATUS "Using MPI post-flags: ${MPIEXEC_POSTFLAGS}")
		message(STATUS "Using MPI includes: ${MPI_CXX_INCLUDE_PATH}")
		message(STATUS "Using MPI libraries: ${MPI_CXX_LIBRARIES}")
		message(STATUS "Using MPI compile flags: ${MPI_CXX_COMPILE_FLAGS}")
		message(STATUS "Using MPI link flags: ${MPI_CXX_LINK_FLAGS}")

		#We need the MPI headers to be included!
		include_directories(${MPI_CXX_INCLUDE_PATH})
	endif()
endif()

#possible future CUDA test.
if(CUDA)
	find_package(CUDA REQUIRED)
	if(NOT CUDA_FOUND)
		MESSAGE(FATAL_ERROR "CUDA enabled build requested but CUDA was not found!")
	endif()
	if(CUDA_FOUND)
		message(STATUS "Found CUDA version ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")
		message(STATUS "Includes at ${CUDA_INCLUDE_DIRS}")
		ADD_DEFINITIONS("-I${CUDA_INCLUDE_DIRS}")
	endif()
endif(CUDA)

#Possible check for Root
if(ROOT)
	#/usr/share/root/cmake/FindROOT.cmake
	SET(CMAKE_MODULE_PATH /usr/share/root/cmake/)
	find_package(ROOT REQUIRED)
	include_directories(${ROOT_INCLUDE_DIR})

	if(NOT ROOT_FOUND)
		MESSAGE(FATAL_ERROR "Root enabled build requested but root was not found!")
	endif()
endif(ROOT)

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
	if(NOT COLLIMATE)
		message(FATAL_ERROR "Building with HDF5 requires COLLIMAT enabled")
	endif()
endif()
if(BONIC)
	if(NOT CPSS)
		message(FATAL_ERROR "Building with HDF5 requires CPSS enabled")
	endif()
	if(NOT CRLIBM)
		message(FATAL_ERROR "Building with BONIC requires CRLIBM enabled")
	endif()
endif()
if(API)
	if(NOT BONIC)
		message(FATAL_ERROR "Building with API requires BONIC enabled")
	endif()
endif()
if(BEAMGAS)
	if(NOT COLLIMAT)
		message(FATAL_ERROR "Building with BEAMGAS requires COLLIMAT enabled")
	endif()
	if(BIGNBLZ)
		message(FATAL_ERROR "Cannot build BEAMGAS with BIGNBLZ enabled.")
	endif()
	if(HUGENBLZ)
		message(FATAL_ERROR "Cannot build BEAMGAS with HUGENBLZ enabled.")
	endif()
endif()
if(BIGNBLZ)
	if(HUGENBLZ)
		message(FATAL_ERROR "Cannot combine BIGNBLZ with HUGENBLZ.")
	endif()
endif()
if(HUGENPART)
	if(BIGNPART)
		message(FATAL_ERROR "Cannot combine BIGNPART with HUGENPART.")
	endif()
	if(COLLIMAT)
		message(FATAL_ERROR "Cannot build HUGENPART with COLLIMAT enabled.")
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
	if(NAGLIB)
		message(FATAL_ERROR "Cannot build NAGLIB with 64BIT enabled.")
	endif()
	if(CERNLIB)
		message(FATAL_ERROR "Cannot build CERNLIB with 64BIT enabled.")
	endif()
	if(NAGFOR)
		message(FATAL_ERROR "Cannot build NAGFOR with 64BIT enabled.")
	endif()
	if(BONIC)
		message(FATAL_ERROR "Cannot build BONIC with 64BIT enabled.")
	endif()
endif()

###################################################################################################
# Settings for compiler build flags
###################################################################################################
#Additional Compiler flags (for all options)
#ADD_DEFINITIONS("-Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -pedantic -march=native")

if(32BIT)
	ADD_DEFINITIONS("-m32")
endif()

if(64BIT)
	ADD_DEFINITIONS("-m64")
endif()

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	# gfortran
	set (CMAKE_Fortran_FLAGS_RELEASE "-frecord-marker=4 -fno-second-underscore -funroll-loops -mfpmath=sse -msse2 -O3")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-frecord-marker=4 -fno-second-underscore -funroll-loops -mfpmath=sse -msse2 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
	# ifort (untested)
	set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
	# g77
	set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
	message ("No optimized Fortran compiler flags are known, we just try -O2...")
	set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

#MPI default build options
SET( CMAKE_CXX_FLAGS_MPI "-DNDEBUG -Wall -O3 -pedantic -DENABLE_MPI -march=native" CACHE STRING
	"Flags used by the C++ compiler during MPI builds."
	FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_MPI
	"-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
	"Flags used for linking binaries during MPI builds."
	FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_MPI
	"-Wl,--warn-unresolved-symbols,--warn-once -fPIC -flto" CACHE STRING
	"Flags used by the shared libraries linker during MPI builds."
	FORCE )
MARK_AS_ADVANCED(
	CMAKE_CXX_FLAGS_MPI
	#CMAKE_C_FLAGS_MPI
	CMAKE_EXE_LINKER_FLAGS_MPI
	CMAKE_SHARED_LINKER_FLAGS_MPI)

#MPI debug default build options
SET( CMAKE_CXX_FLAGS_MPIDEBUG "-Wall -g3 -ggdb -pedantic -DENABLE_MPI -march=native" CACHE STRING
	"Flags used by the C++ compiler during MPI Debug builds."
	FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_MPIDEBUG
	"-Wl,--warn-unresolved-symbols,--warn-once" CACHE STRING
	"Flags used for linking binaries during MPI Debug builds."
	FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_MPIDEBUG
	"-Wl,--warn-unresolved-symbols,--warn-once -fPIC" CACHE STRING
	"Flags used by the shared libraries linker during MPI Debug builds."
	FORCE )
MARK_AS_ADVANCED(
	CMAKE_CXX_FLAGS_MPIDEBUG
	#CMAKE_C_FLAGS_MPIDEBUG
	CMAKE_EXE_LINKER_FLAGS_MPIDEBUG
	CMAKE_SHARED_LINKER_FLAGS_MPIDEBUG )

###################################################################################################
# Settings for what exactly to build
###################################################################################################
if(CRLIBM)
	file(GLOB CRLIBM_SOURCES crlibm/*.c)
	add_library(crlibm STATIC ${CRLIBM_SOURCES})
	if(32BIT)
		SET_TARGET_PROPERTIES(crlibm PROPERTIES COMPILE_FLAGS "-DLINUX_INLINE -DHAVE_CONFIG_H -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith  -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing -m32")
		SET_TARGET_PROPERTIES(crlibm PROPERTIES LINK_FLAGS "-DLINUX_INLINE -DHAVE_CONFIG_H -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith  -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing -m32")
	endif()
	if(64BIT)
		SET_TARGET_PROPERTIES(crlibm PROPERTIES COMPILE_FLAGS "-DLINUX_INLINE -DHAVE_CONFIG_H -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith  -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing -m64")
		SET_TARGET_PROPERTIES(crlibm PROPERTIES LINK_FLAGS "-DLINUX_INLINE -DHAVE_CONFIG_H -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith  -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing -m64")
	endif()
endif()

#build astuce
add_executable(astuce astuce.f)

#build dafor
add_executable(dafor dafor.f)

if(32BIT)
	SET_TARGET_PROPERTIES(astuce PROPERTIES COMPILE_FLAGS "-m32")
	SET_TARGET_PROPERTIES(dafor PROPERTIES COMPILE_FLAGS "-m32")
	SET_TARGET_PROPERTIES(astuce PROPERTIES LINK_FLAGS "-m32")
	SET_TARGET_PROPERTIES(dafor PROPERTIES LINK_FLAGS "-m32")
endif()

if(32BIT)
	SET_TARGET_PROPERTIES(${SIXTRACK_BINARY_NAME} PROPERTIES LINK_FLAGS "-m32")
endif()

#if we have a boinc build, add myboinc.f
if(BOINC)
	add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/myboinc.f ${CMAKE_BINARY_DIR})
	set(SIXTRACK_SOURCES ${SIXTRACK_SOURCES} myboinc.f)
endif()

add_custom_target(copy_input_files)
add_dependencies(copy_input_files ast_mkdir)
file(GLOB SIXTRACK_INPUT ${CMAKE_SOURCE_DIR}/*.s)
foreach(loop IN ITEMS ${SIXTRACK_INPUT})
	add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${loop} ${CMAKE_BINARY_DIR})
endforeach()

add_custom_target(ast_mkdir)
file(GLOB ASTUCE_MASKS ${CMAKE_SOURCE_DIR}/ast_mask/*.ast)
add_custom_command(TARGET ast_mkdir POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/ast_mask/)
foreach(loop IN ITEMS ${ASTUCE_MASKS})
	add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${loop} ${CMAKE_BINARY_DIR}/ast_mask/ DEPENDS ast_mkdir)
endforeach()

#We need to build asatuce and dafor before we build SixTrack, but crlibm and so on can be built at the same time.
add_dependencies(${SIXTRACK_BINARY_NAME} copy_input_files generate_input)

###################################################################################################
#running astuce and dafor
###################################################################################################
add_custom_target(run_astuce)

#make running astuce depend on astuce being built
add_dependencies(run_astuce astuce astuce_clean_a astuce_clean_b)
add_custom_target(astuce_clean_a)
add_custom_target(astuce_clean_b)
add_dependencies(astuce_clean_a copy_input_files)
add_dependencies(astuce_clean_b astuce_clean_a)

macro(ast_masks OPT FILE_NAME)
	add_custom_command(TARGET astuce_clean_a POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -DFILE_NAME:FILEPATH=${CMAKE_BINARY_DIR}/ast_mask/${FILE_NAME}.ast -DOPT:STRING=${OPT} -P ${CMAKE_SOURCE_DIR}/masks.cmake) 
endmacro()

macro(ast_masks2 OPT)
	add_custom_command(TARGET astuce_clean_b POST_BUILD COMMAND sed -i -e '/^e  *'${OPT}'/d' ${CMAKE_BINARY_DIR}/ast_mask/${loop}.ast)
endmacro()

#First clean up the astuce decks for only what is enabled
#set(OPT crlibm tilt fast bpm beamgas collimat bnlelens cr boinc nagfor bignblz debug hdf5 fio lf95 hugenblz stf bignpart datamods hugenpart)
foreach(loop IN ITEMS ${FORS})
	if(NOT CRLIBM)
		ast_masks(crlibm ${loop})
	endif()
	if(NOT TILT)
		ast_masks(tilt ${loop})
	endif()
	if(NOT FAST)
		ast_masks(fast ${loop})
	endif()
	if(NOT BPM)
		ast_masks(bpm ${loop})
	endif()
	if(NOT BEAMGAS)
		ast_masks(beamgas ${loop})
	endif()
	if(NOT COLLIMAT)
		ast_masks(collimat ${loop})
	endif()
	if(NOT BNLELENS)
		ast_masks(bnlelens ${loop})
	endif()
	if(NOT CR)
		ast_masks(cr ${loop})
	endif()
	if(NOT BOINC)
		ast_masks(boinc ${loop})
	endif()
	if(NOT NAGFOR)
		ast_masks(nagfor ${loop})
	endif()
	if(NOT BIGNBLZ)
		ast_masks(bignblz ${loop})
	endif()
	if(NOT DEBUG)
		ast_masks(debug ${loop})
	endif()
	if(NOT HDF5)
		ast_masks(hdf5 ${loop})
	endif()
	if(NOT FIO)
		ast_masks(fio ${loop})
	endif()
	if(NOT LF95)
		ast_masks(lf95 ${loop})
	endif()
	if(NOT HUGENBLZ)
		ast_masks(hugenblz ${loop})
	endif()
	if(NOT STF)
		ast_masks(stf ${loop})
	endif()
	if(NOT BIGNPART)
		ast_masks(bignpart ${loop})
	endif()
	if(NOT DATAMODS)
		ast_masks(datamods ${loop})
	endif()
	if(NOT HUGENPART)
		ast_masks(hugenpart ${loop})
	endif()

	if(NOT COLLIMAT)
		ast_masks2(nwrtcoll)
	endif()
	if(NOT DATAMODS)
		ast_masks2(datamods)
	endif()
	if(NOT DEBUG)
		ast_masks2(dumps)
	endif()
	if(NOT CR)
		ast_masks2(checkpt)
	endif()
	if(NOT BNLELENS)
		ast_masks2(nwrtbnl)
	endif()

	#From make_six:
	#NOTA BENE! Definitions of cernlib and naglib are reversed
	# in the sense that while a "" null value implies option is ON
	# the reverse is true for cernlib and naglib
	if(CERNLIB)
		ast_masks2(plotdumy)
	endif()
	if(NAGLIB)
		ast_masks2(nagdumy)
	endif()
endforeach(loop)

add_custom_target(generate_input)
add_dependencies(generate_input run_astuce dafor)
foreach(loop IN ITEMS ${FORS})
	#This runs astuce
	add_custom_command(TARGET run_astuce POST_BUILD COMMAND astuce < ${CMAKE_BINARY_DIR}/ast_mask/${loop}.ast )

	#We make custom targets/commands and deps for each dafor step
	#make a subdir so we can build in parallel
	add_custom_target(dafor_mkdir_${loop})

	#do the symlink/copy
	add_custom_target(dafor_copy_${loop})

	#do the symlink/copy
	add_custom_target(dafor_copy_binary_${loop})

	#run dafor
	add_custom_target(dafor_run_${loop})

	#remove the link and other files
	add_custom_target(dafor_remove_${loop})

	#copy the output we need
	add_custom_target(dafor_move_${loop})

	#chain together the dependencies, each with a custom command
	add_dependencies(dafor_remove_${loop} dafor_move_${loop})
	add_dependencies(dafor_move_${loop} dafor_run_${loop})
	add_dependencies(dafor_run_${loop} dafor_copy_${loop} dafor_copy_binary_${loop})
	add_dependencies(dafor_copy_${loop} dafor_mkdir_${loop} run_astuce)
	add_dependencies(dafor_copy_binary_${loop} dafor_mkdir_${loop} dafor)

	if(${loop} STREQUAL sixvefox OR ${loop} STREQUAL sixscfox OR ${loop} STREQUAL sixda)
		message(STATUS "Will run dafor on ${loop}.f")
		add_dependencies(generate_input dafor_move_${loop})

		#make a subdir so we can build in parallel
		add_custom_command(TARGET dafor_mkdir_${loop} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${loop})

		#do the copy (not symlink so it will work on windows too)
		add_custom_command(TARGET dafor_copy_${loop} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${loop}n.f ${loop}/fort.1)

		#do the copy of dafor (not symlink so it will work on windows too)
		add_custom_command(TARGET dafor_copy_binary_${loop} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy dafor ${loop}/dafor)

		#run dafor
		add_custom_command(TARGET dafor_run_${loop} COMMAND ${CMAKE_COMMAND} -E chdir ${loop} ./dafor)

		#copy the output we need
		add_custom_command(TARGET dafor_move_${loop} COMMAND ${CMAKE_COMMAND} -E copy ${loop}/fort.2 ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f)

		#remove the link and other files
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/${loop} DEPENDS dafor_remove_${loop})
	else()
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f COMMAND ${CMAKE_COMMAND} -E copy ${loop}n.f ${loop}.f DEPENDS run_astuce)
	endif()
endforeach()

###################################################################################################
#linking
###################################################################################################

#if we are using crlibm, link to the crlibm library
if(CRLIBM)
	target_link_libraries(${SIXTRACK_BINARY_NAME} crlibm)
endif(CRLIBM)

#if we are using root, link to all the root libraries
if(ROOT)
	target_link_libraries(${SIXTRACK_BINARY_NAME} ${ROOT_LIBRARIES})
endif(ROOT)

