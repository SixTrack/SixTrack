! =================================================================================================
!  STANDARD OUTPUT MODULE
!  Last modified: 2018-03-22
!  For CR version, this is the "buffer file" fort.92;
!  Otherwise write directly to "*" aka iso_fortran_env::output_unit (usually unit 6)
! =================================================================================================
module crcoall
  
  implicit none
  
  integer lout
  save lout
  
end module crcoall

! =================================================================================================
!  FLOAT PRECISION MODULE
!  Last modified: 2018-05-25
! =================================================================================================
module floatPrecision
  use, intrinsic :: iso_fortran_env, only : real32, real64, real128
  implicit none
#ifdef P32BITM
  integer, parameter :: fPrec = real32
#endif
#ifdef P64BITM
  integer, parameter :: fPrec = real64
#endif
#ifdef P128BITM
  integer, parameter :: fPrec = real128
#endif
end module floatPrecision
