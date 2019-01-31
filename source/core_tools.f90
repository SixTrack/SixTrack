! =================================================================================================
!  STANDARD OUTPUT MODULE
!  Last modified: 2018-03-22
!  For CR version, this is the "buffer file" fort.92;
!  Otherwise write directly to "*" aka iso_fortran_env : output_unit (usually unit 6)
! =================================================================================================
module crcoall
  use, intrinsic :: iso_fortran_env, only : output_unit
  implicit none
  integer, public, save :: lout = output_unit
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

! ================================================================================================ !
!  NUMERICAL PRECISION CHECKS
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-09
! ================================================================================================ !
logical function checkMultUnderflow(varA, varB)
  use, intrinsic :: iso_fortran_env, only : real32, real64, real128
  real(kind=SIXTRACK_REAL) varA, varB
  checkMultUnderflow = .false.
  if(varA /= 0.0 .and. varB /= 0.0) then
    if(log10(abs(varA)) + log10(abs(varB)) < -305) then
      checkMultUnderflow = .true.
    end if
  end if
end function checkMultUnderflow
