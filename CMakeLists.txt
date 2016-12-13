#Sixtrack CMAKE build file
#Making use of https://cmake.org/Wiki/CMakeForFortranExample
#and the file created for building Merlin.

###################################################################################################
# General configuration
###################################################################################################

#Set the cmake version required
cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)

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

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs to include the MPI options
SET(CMAKE_BUILD_TYPE Release CACHE STRING
#SET(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel MPI MPIDEBUG."
	FORCE)

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
OPTION(API "Link to external BOINC API" OFF)
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

#Compiler stuff
OPTION(IFORT "Use the intel fortran compiler" OFF)
OPTION(LF95 "Use the Lahey-Fujitsu lf95 fortran compiler" OFF)
OPTION(PGF90 "Use the pgf90 fortran compiler" OFF)
OPTION(F95 "Use the g95 fortran compiler" OFF)
OPTION(GFORTRAN "Use the gfortran fortran compiler" OFF)
OPTION(NAGFOR "Use NAGFOR compiler" OFF)

#Build options
OPTION(OPENMP "Enable OpenMP support" OFF)
OPTION(MPI "Enable MPI support" OFF)
OPTION(CUDA "Enable CUDA support" OFF)
OPTION(32BIT "Create a 32bit binary" OFF)
OPTION(64BIT "Create a 64bit binary" ON)
OPTION(STATIC "Create a statically linked binary" ON)
OPTION(NATIVE "Enable optimizations for the current machine" OFF)

OPTION(BUILD_TESTING "Enable the SixTest suite" OFF)

OPTION(ROOT "Link to the root libraries" OFF)

###################################################################################################
# Sixtrack binary name construction
###################################################################################################

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

SET(SIXTRACK_FEATURES "cmake")
SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_${Fortran_COMPILER_NAME}")

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
if(MERLINSCATTER)
	SET(SIXTRACK_FEATURES "merlinscatter_${SIXTRACK_FEATURES}")
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
if(STATIC)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_static")
endif()
if(NATIVE)
	SET(SIXTRACK_FEATURES "${SIXTRACK_FEATURES}_native")
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
    SET(FORS sixda sixsc sixscfox dabnew lielib)
else()
    SET(FORS track sixve sixvefox dabnews lielib plato_seq)
endif()

if(NOT DA AND DATAMODS)
	SET(FORS datamods ${FORS})
endif()

foreach(loop IN ITEMS ${FORS})
	set(SIXTRACK_SOURCES ${SIXTRACK_SOURCES} ${loop}.f)
endforeach()

#add back in beamgas, which is in free format
SET(SIXTRACK_SOURCES ${SIXTRACK_SOURCES} beamgas.f90)
SET(FORS ${FORS} beamgas)

###################################################################################################
# Checking if some requested packages exist
###################################################################################################
#Check for HDF5
if(HDF5)
	if(STATIC)
		SET(HDF5_USE_STATIC_LIBRARIES "True")
	endif(STATIC)

	find_package(HDF5 COMPONENTS Fortran REQUIRED)
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

#Check for OpenMP
if(OPENMP)
	find_package(OpenMP REQUIRED)
endif(OPENMP)

#If we have OpenMP enabled, set the compiler flags needed
if(OPENMP_FOUND)
	SET_TARGET_PROPERTIES(SixTrack PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
	SET_TARGET_PROPERTIES(SixTrack PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
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
	if(NOT COLLIMAT)
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

if(STATIC)
	#This just sets the default build setting for shared libraries.
	SET(BUILD_SHARED_LIBS OFF)
endif()

#FLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

#lots of games to play here.... start with fortran flags
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	# gfortran
	SET(CMAKE_Fortran_FLAGS_RELEASE "-frecord-marker=4 -fno-second-underscore -funroll-loops -mfpmath=sse -msse2 -O3")
	SET(CMAKE_Fortran_FLAGS_DEBUG   "-frecord-marker=4 -fno-second-underscore -funroll-loops -mfpmath=sse -msse2 -g -ffpe-trap=invalid,zero,overflow,underflow -fcheck=all")

	if(32BIT)
		set(CMAKE_Fortran_FLAGS_RELEASE "-m32 ${CMAKE_Fortran_FLAGS_RELEASE}")
		set(CMAKE_Fortran_FLAGS_DEBUG "-m32 ${CMAKE_Fortran_FLAGS_DEBUG}")
		SET(CMAKE_Fortran_LINKER_FLAGS "-m32")
	endif(32BIT)

	if(64BIT)
		set(CMAKE_Fortran_FLAGS_RELEASE "-m64 ${CMAKE_Fortran_FLAGS_RELEASE}")
		set(CMAKE_Fortran_FLAGS_DEBUG "-m64 ${CMAKE_Fortran_FLAGS_DEBUG}")
		SET(CMAKE_Fortran_LINKER_FLAGS "-m64")
	endif(64BIT)

	if(NATIVE)
		SET(CMAKE_Fortran_FLAGS_RELEASE "-march=native -mno-fma4 -mno-fma ${CMAKE_Fortran_FLAGS_RELEASE}")
	endif(NATIVE)

	if(STATIC)
		SET(CMAKE_Fortran_LINKER_FLAGS "-static")
	endif(STATIC)

elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
	# ifort
	set (CMAKE_Fortran_FLAGS_RELEASE "-fp-model source -fp-model strict -fp-model no-except -O3 -no-fma")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-fp-model source -fp-model strict -fp-model no-except -O0 -g -no-fma")

	if(NATIVE)
		SET(CMAKE_Fortran_FLAGS_RELEASE "-xHost ${CMAKE_Fortran_FLAGS_RELEASE}")
	endif(NATIVE)

	if(STATIC)
		SET(CMAKE_Fortran_LINKER_FLAGS "-Bstatic")
	endif(STATIC)

	if(NOT STATIC)
		SET(CMAKE_Fortran_LINKER_FLAGS "-Bdynamic")
	endif(NOT STATIC)

elseif (Fortran_COMPILER_NAME MATCHES "nagfor.*")
	# nagfor
	#Here I am going to assume that since we are using nagfor we are on lxplus for now
	set (CMAKE_Fortran_FLAGS_RELEASE "-O4 -dusty -ieee=full -dcfuns")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -dusty -ieee=full -dcfuns")

	if(32BIT)
		set(CMAKE_Fortran_FLAGS_RELEASE "-abi=32 ${CMAKE_Fortran_FLAGS_RELEASE}")
		set(CMAKE_Fortran_FLAGS_DEBUG "-abi=32 ${CMAKE_Fortran_FLAGS_DEBUG}")
		SET(CMAKE_Fortran_LINKER_FLAGS "-abi=32 ${CMAKE_Fortran_LINKER_FLAGS}")
		find_library(g2cLink NAMES g2c PATHS /usr/lib /usr/lib/gcc/x86_64-redhat-linux/3.4.6/32 )
	endif(32BIT)

	if(64BIT)
		set(CMAKE_Fortran_FLAGS_RELEASE "-abi=64t ${CMAKE_Fortran_FLAGS_RELEASE}")
		set(CMAKE_Fortran_FLAGS_DEBUG "-abi=64t ${CMAKE_Fortran_FLAGS_DEBUG}")
		SET(CMAKE_Fortran_LINKER_FLAGS "-abi=64t ${CMAKE_Fortran_LINKER_FLAGS}")
		find_library(g2cLink NAMES g2c PATHS /usr/lib64 /usr/lib/gcc/x86_64-redhat-linux/3.4.6 )
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

elseif (Fortran_COMPILER_NAME MATCHES "g77")
	# g77
	set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")

#generic catch option
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
	message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
	message ("No optimized Fortran compiler flags are known, we just try -O2...")
	set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
	set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

#Now for C flags (crlibm for now)
#${CMAKE_C_COMPILER_ID} : From the docs: "one of "Clang", "GNU", "Intel", or "MSVC". This works even if a compiler wrapper like ccache is used."
if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")

	if(32BIT)
		set(CMAKE_C_FLAGS "-m32")
	endif(32BIT)

	if(64BIT)
		set(CMAKE_C_FLAGS "-m64")
	endif(64BIT)

	if(NATIVE)
		set(CMAKE_C_FLAGS "-march=native")
	endif(NATIVE)

	set (CMAKE_C_FLAGS_RELEASE "-O3 ${CMAKE_C_FLAGS}")
	set (CMAKE_C_FLAGS_DEBUG   "-O0 -g ${CMAKE_C_FLAGS}")

	if(STATIC)
		SET(CMAKE_C_LINKER_FLAGS "-static")
	endif(STATIC)

endif()

if(${CMAKE_C_COMPILER_ID} MATCHES "Intel")

	if(32BIT)
		set(CMAKE_C_FLAGS "-m32")
	endif(32BIT)

	if(64BIT)
		set(CMAKE_C_FLAGS "-m64")
	endif(64BIT)

	if(NATIVE)
		set(CMAKE_C_FLAGS "-xHost")
	endif(NATIVE)

	set (CMAKE_C_FLAGS_RELEASE "-O3 ${CMAKE_C_FLAGS}")
	set (CMAKE_C_FLAGS_DEBUG   "-O0 -g ${CMAKE_C_FLAGS}")

	if(STATIC)
		SET(CMAKE_C_LINKER_FLAGS "-Bstatic")
	endif(STATIC)

	if(NOT STATIC)
		SET(CMAKE_C_LINKER_FLAGS "-Bdynamic")
	endif(NOT STATIC)

endif()

#Now for CXX flags (merlinScatter for now)
#${CMAKE_CXX_COMPILER_ID} : From the docs: "one of "Clang", "GNU", "Intel", or "MSVC". This works even if a compiler wrapper like ccache is used."
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")

	if(32BIT)
		set(CMAKE_CXX_FLAGS "-m32 ${CMAKE_CXX_FLAGS}")
	endif(32BIT)

	if(64BIT)
		set(CMAKE_CXX_FLAGS "-m64 ${CMAKE_CXX_FLAGS}")
	endif(64BIT)

	if(NATIVE)
		set(CMAKE_CXX_FLAGS "-march=native ${CMAKE_CXX_FLAGS}")
	endif(NATIVE)

	set (CMAKE_CXX_FLAGS_RELEASE "-O3 ${CMAKE_CXX_FLAGS}")
	set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g ${CMAKE_CXX_FLAGS}")

	if(STATIC)
		SET(CMAKE_CXX_LINKER_FLAGS "-static")
	endif(STATIC)

endif()

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")

	if(32BIT)
		set(CMAKE_CXX_FLAGS "-m32 ${CMAKE_CXX_FLAGS}")
	endif(32BIT)

	if(64BIT)
		set(CMAKE_CXX_FLAGS "-m64 ${CMAKE_CXX_FLAGS}")
	endif(64BIT)

	if(NATIVE)
		set(CMAKE_CXX_FLAGS "-march=xHost ${CMAKE_CXX_FLAGS}")
	endif(NATIVE)

	set (CMAKE_CXX_FLAGS_RELEASE "-O3 ${CMAKE_CXX_FLAGS}")
	set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g ${CMAKE_CXX_FLAGS}")

	if(STATIC)
		SET(CMAKE_CXX_LINKER_FLAGS "-Bstatic")
	endif(STATIC)

	if(NOT STATIC)
		SET(CMAKE_CXX_LINKER_FLAGS "-Bdynamic")
	endif(NOT STATIC)
endif()

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
#add the source files that will be used to build sixtrack
add_executable(SixTrack ${SIXTRACK_SOURCES})
SET_TARGET_PROPERTIES(SixTrack PROPERTIES OUTPUT_NAME "${SIXTRACK_BINARY_NAME}")

#build astuce
add_executable(astuce ${CMAKE_SOURCE_DIR}/SixTrack/astuce.f)
if(Fortran_COMPILER_NAME MATCHES "nagfor.*")
	target_link_libraries(astuce ${g2cLink})
endif()

#build dafor
add_executable(dafor ${CMAKE_SOURCE_DIR}/SixTrack/dafor.f)

if(CRLIBM)
	file(GLOB CRLIBM_SOURCES ${CMAKE_SOURCE_DIR}/SixTrack/crlibm/*.c)
	add_library(crlibm STATIC ${CRLIBM_SOURCES})

	target_compile_definitions(crlibm PRIVATE -DLINUX_INLINE -DHAVE_CONFIG_H)

	if(32BIT)
		target_compile_options(crlibm PRIVATE -m32)
	endif(32BIT)
	if(64BIT)
		target_compile_options(crlibm PRIVATE -m64)
	endif(64BIT)

	if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")
		target_compile_options(crlibm PRIVATE -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing)
	elseif(${CMAKE_C_COMPILER_ID} MATCHES "Intel")
		target_compile_options(crlibm PRIVATE -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wconversion -Wstrict-prototypes -Winline -pedantic -fno-strict-aliasing)
	else()
		target_compile_options(crlibm PRIVATE -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing)
	endif()
endif()

if(HDF5_FOUND)
	SET_TARGET_PROPERTIES(SixTrack PROPERTIES COMPILE_FLAGS "-I${HDF5_INCLUDE_DIRS} ${HDF5_DEFINITIONS}")
endif(HDF5_FOUND)

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

add_custom_target(copy_input_files)
add_dependencies(copy_input_files ast_mkdir)
file(GLOB SIXTRACK_INPUT ${CMAKE_SOURCE_DIR}/SixTrack/*.s)

#copy all the input files *.s to our build folder
foreach(loop IN ITEMS ${SIXTRACK_INPUT})
	add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${loop} ${CMAKE_BINARY_DIR})
endforeach()

#if we have a boinc build, add myboinc.f
if(BOINC)
	add_custom_command(TARGET copy_input_files POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/SixTrack/myboinc.f ${CMAKE_BINARY_DIR})
	set(SIXTRACK_SOURCES ${SIXTRACK_SOURCES} myboinc.f)
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
add_dependencies(run_astuce astuce astuce_clean_define_flags astuce_clean_decks)

#This target cleans the define flags
add_custom_target(astuce_clean_define_flags)
#this requires the input files to be in place
add_dependencies(astuce_clean_define_flags copy_input_files)

#this target cleans the decks
add_custom_target(astuce_clean_decks)
#and requires the define flags to already be cleaned (this is just to stop parallel writes on the same file)
add_dependencies(astuce_clean_decks astuce_clean_define_flags)

#This macro calls the regex in masks_define_flags.cmake to remove the requested flags
macro(ast_masks OPT FILE_NAME)
	add_custom_command(TARGET astuce_clean_define_flags POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -DFILE_NAME:FILEPATH=${CMAKE_BINARY_DIR}/ast_mask/${FILE_NAME}.ast -DOPT:STRING=${OPT} -P ${CMAKE_SOURCE_DIR}/SixTrack/masks_define_flags.cmake)
endmacro()

#This macro calls the regex in masks_decks.cmake to remove the requested decks
macro(ast_masks2 OPT FILE_NAME)
	add_custom_command(TARGET astuce_clean_decks POST_BUILD COMMAND ${CMAKE_COMMAND} ARGS -DFILE_NAME:FILEPATH=${CMAKE_BINARY_DIR}/ast_mask/${FILE_NAME}.ast -DOPT:STRING=${OPT} -P ${CMAKE_SOURCE_DIR}/SixTrack/masks_decks.cmake)
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
	if(NOT MERLINSCATTER)
		ast_masks(merlinscatter ${loop})
	endif()

	if(NOT COLLIMAT)
		ast_masks2(nwrtcoll ${loop})
	endif()
	if(NOT DATAMODS)
		ast_masks2(datamods ${loop})
	endif()
	if(NOT DEBUG)
		ast_masks2(dumps ${loop})
	endif()
	if(NOT CR)
		ast_masks2(checkpt ${loop})
	endif()
	if(NOT BNLELENS)
		ast_masks2(nwrtbnl ${loop})
	endif()

	#From make_six:
	#NOTA BENE! Definitions of cernlib and naglib are reversed
	# in the sense that while a "" null value implies option is ON
	# the reverse is true for cernlib and naglib
	if(CERNLIB)
		ast_masks2(plotdumy ${loop})
	endif()
	if(NAGLIB)
		ast_masks2(nagdumy ${loop})
	endif()
endforeach(loop)

add_custom_target(generate_input)
add_dependencies(generate_input run_astuce dafor rename_beamgas)

#we rename beamgas.f to beamgas.f90 since it should be in free format fortran
add_custom_target(rename_beamgas)
add_dependencies(rename_beamgas run_astuce)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/beamgas.f90 COMMAND ${CMAKE_COMMAND} -E rename beamgasn.f beamgas.f90 DEPENDS rename_beamgas)

foreach(loop IN ITEMS ${FORS})
	#This runs astuce
	add_custom_command(TARGET run_astuce POST_BUILD COMMAND astuce < ${CMAKE_CURRENT_BINARY_DIR}/ast_mask/${loop}.ast  > ${loop}_ast.log)

	if(${loop} STREQUAL sixvefox OR ${loop} STREQUAL sixscfox OR ${loop} STREQUAL sixda)
		message(STATUS "Will run dafor on ${loop}.f")
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f COMMAND ${CMAKE_CURRENT_BINARY_DIR}/dafor ${loop}n.f ${loop}.f > ${loop}_dafor.log DEPENDS run_astuce)
	else()
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${loop}.f COMMAND ${CMAKE_COMMAND} -E copy ${loop}n.f ${loop}.f DEPENDS run_astuce)
	endif()
endforeach()

###################################################################################################
#linking
###################################################################################################

if(CERNLIB)
#	if(EXISTS "/cern/pro/bin/cernlib")
		#running on lxplus and friends
#		SET_TARGET_PROPERTIES(SixTrack PROPERTIES LINK_FLAGS "`cernlib -u graflib packlib_noshift kernlib_noshift` -lg2c -lpthread")
#	else()
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
	target_link_libraries(SixTrack crlibm)
endif(CRLIBM)

#if we are using merlinscatter, link to the crlibm library
if(MERLINSCATTER)
	target_link_libraries(SixTrack merlinscatter)
endif(MERLINSCATTER)

#if we are using root, link to all the root libraries
if(ROOT)
	target_link_libraries(SixTrack ${ROOT_LIBRARIES})
endif(ROOT)

if(HDF5)
	target_link_libraries(SixTrack ${HDF5_LIBRARIES})
endif(HDF5)

SET(BUILDNAME "${SIXTRACK_BINARY_NAME}" CACHE STRING "Name of build on the dashboard")
MARK_AS_ADVANCED(BUILDNAME)

#For testing
IF(BUILD_TESTING)

	ENABLE_TESTING()
	INCLUDE(CTest)

	add_subdirectory(${CMAKE_SOURCE_DIR}/SixTest)

	configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake  ${CMAKE_BINARY_DIR})

	#need to build the testing tools
	add_executable(verify10 ${CMAKE_SOURCE_DIR}/SixTest/bin/verify10.f)
	add_executable(checkf10 ${CMAKE_SOURCE_DIR}/SixTest/bin/checkf10.f)
	add_executable(compf10 ${CMAKE_SOURCE_DIR}/SixTest/bin/compf10.f)
	add_executable(checkf110 ${CMAKE_SOURCE_DIR}/SixTest/bin/checkf110.f)
	add_executable(checkf1014 ${CMAKE_SOURCE_DIR}/SixTest/bin/checkf1014.f)
	add_executable(read90 ${CMAKE_SOURCE_DIR}/SixTest/read90/read90.f90 ${CMAKE_SOURCE_DIR}/SixTrack/crlibm/dtoa_c.c ${CMAKE_SOURCE_DIR}/SixTrack/crlibm/dtoaf.c)

	if(Fortran_COMPILER_NAME MATCHES "nagfor.*")
		target_link_libraries(verify10 ${g2cLink})
		target_link_libraries(checkf10 ${g2cLink})
		target_link_libraries(compf10 ${g2cLink})
	endif()

	#The testing wrapper
	add_executable(SixTestWrapper ${CMAKE_SOURCE_DIR}/SixTest/bin/SixTestWrapper.cpp)
	#set_property(TARGET SixTestWrapper PROPERTY CXX_STANDARD 11)

ENDIF()
