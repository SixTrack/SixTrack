! =================================================================================================
!  STANDARD OUTPUT MODULE
!  Last modified: 2018-03-22
!  For CR version, this is the "buffer file" fort.92;
!  Otherwise write directly to "*" aka iso_fortran_env : output_unit (usually unit 6)
! =================================================================================================
module crcoall
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
  implicit none
  integer, public, save :: lout = output_unit
  integer, public, save :: lerr = error_unit
end module crcoall

! =================================================================================================
!  FLOAT PRECISION MODULE
!  Last modified: 2018-05-25
! =================================================================================================
module floatPrecision
  use, intrinsic :: iso_fortran_env, only : real32, real64, real128
  implicit none
  integer, parameter :: fPrec = SIXTRACK_REAL
end module floatPrecision
