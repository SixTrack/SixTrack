#CHECK_Fortran_SOURCE_COMPILES(<code> <var> [FAIL_REGEX <fail-regex>] [SRC_EXT <ext>])
include(CheckFortranSourceCompiles)

SET(FORTRAN_TEST_SOURCE
"
      program ftest
      end program
")
check_fortran_source_compiles(${FORTRAN_TEST_SOURCE} FORTRAN_TEST_BUILD)
if(NOT FORTRAN_TEST_BUILD)
    message(SEND_ERROR "ERROR: Your compiler failed to build a trivial executable. This is very bad!")
    SET(FAILURE 1)
endif(NOT FORTRAN_TEST_BUILD)

#test for fortran 90 freeform 
SET(FORTRAN_FREEFORM_TEST_SOURCE
"
program ftest
end program
")
#SRC_EXT needs cmake 3.7
check_fortran_source_compiles(${FORTRAN_FREEFORM_TEST_SOURCE} FORTRAN_FREEFORM_TEST_BUILD SRC_EXT ".f90")
if(NOT FORTRAN_FREEFORM_TEST_BUILD)
    message(SEND_ERROR "ERROR: Your compiler failed to build a trivial executable as freeform. This is very bad!")
    SET(FAILURE 1)
endif(NOT FORTRAN_FREEFORM_TEST_BUILD)

#test for fortran 2003 features
SET(FORTRAN_ISO_FORTRAN_ENV_STDOUT_TEST_SOURCE
"
program ftest
  use, intrinsic :: iso_fortran_env, only : output_unit
end program
")
check_fortran_source_compiles(${FORTRAN_ISO_FORTRAN_ENV_STDOUT_TEST_SOURCE} FORTRAN_ISO_FORTRAN_ENV_STDOUT_TEST_BUILD SRC_EXT ".f90")
if(NOT FORTRAN_ISO_FORTRAN_ENV_STDOUT_TEST_BUILD)
    message(SEND_ERROR 
"ERROR: Your compiler does not support the iso_fortran_env intrinsic (output_unit, etc)! This is required to build SixTrack! Please update your compiler to one that supports fortran 2008 (e.g. a recent gfortran build).")
    SET(FAILURE 1)
endif(NOT FORTRAN_ISO_FORTRAN_ENV_STDOUT_TEST_BUILD)


#test for fortran 2008 features
SET(FORTRAN_ISO_FORTRAN_ENV_REAL_TEST_SOURCE
"
program ftest
  use, intrinsic :: iso_fortran_env, only : real32, real64, real128
end program
")
check_fortran_source_compiles(${FORTRAN_ISO_FORTRAN_ENV_REAL_TEST_SOURCE} FORTRAN_ISO_FORTRAN_ENV_REAL_TEST_BUILD SRC_EXT ".f90")
if(NOT FORTRAN_ISO_FORTRAN_ENV_REAL_TEST_BUILD)
    message(SEND_ERROR 
"ERROR: Your compiler does not support the iso_fortran_env intrinsic (real64 types)! This is required to build SixTrack! Please update your compiler to one that supports fortran 2008 (e.g. a recent gfortran build).")
    SET(FAILURE 1)
endif(NOT FORTRAN_ISO_FORTRAN_ENV_REAL_TEST_BUILD)

#special check for lxplus
cmake_host_system_information(RESULT THIS_HOST QUERY HOSTNAME)
message(STATUS "Building on: ${THIS_HOST}")

#does this have lxplus at the start of the hostname
string(FIND ${THIS_HOST} "lxplus" LXPLUS)

if(FAILURE)
    if(NOT LXPLUS)
        message(STATUS "!")
        message(STATUS "!")
        message(STATUS "!")
        message(STATUS "!")
        message(STATUS "!")
        message(STATUS "!")
        message(STATUS "If you have hit this, SixTrack is not going to build!")
        message(STATUS "Your compiler is too old (or maybe just broken).")
        message(STATUS "CMake thinks you are on lxplus")
        message(STATUS "Newer compilers can be easily enabled via CVMFS")
        message(STATUS "!")
        message(STATUS "!")
        message(STATUS "On lxplus, add the following (for example) to your .bash_profile:")
        message(STATUS "source /cvmfs/sft.cern.ch/lcg/contrib/gcc/9/x86_64-centos7-gcc9-opt/setup.sh")
        message(STATUS "!")
        message(STATUS "Or add \"lxplus\" entry (lower case) to your build options if you are using cmake_six")
        message(STATUS "!")
        message(STATUS "!")
        message(STATUS "!")
        message(STATUS "!")
    endif(NOT LXPLUS)

    message(FATAL_ERROR "ERROR: Your compiler lacks the required features from fortran 2003 and 2008 to build SixTrack! You will need to update to a newer version (or a better compiler).")

endif(FAILURE)
