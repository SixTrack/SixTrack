! =================================================================================================
!  STANDARD OUTPUT MODULE
!  Last modified: 2019-05-08
!  Default to standard output and error, but can be set to other units for CR version
! =================================================================================================
module crcoall
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
  implicit none
  integer, public, save :: lout  = output_unit  ! Should be unit 6
  integer, public, save :: lerr  = error_unit   ! Should be unit 0
  integer, public, save :: crlog = -1           ! Writing to this outside of an #ifdef CR will segfault
end module crcoall

! =================================================================================================
!  FLOAT PRECISION MODULE
!  Last modified: 2019-05-08
! =================================================================================================
module floatPrecision
#ifdef SINGLE_MATH
  use, intrinsic :: iso_fortran_env, only : real32
  implicit none
  integer, parameter :: fPrec = real32
#endif
#ifdef DOUBLE_MATH
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  integer, parameter :: fPrec = real64
#endif
#ifdef QUAD_MATH
  use, intrinsic :: iso_fortran_env, only : real128
  implicit none
  integer, parameter :: fPrec = real128
#endif
end module floatPrecision
