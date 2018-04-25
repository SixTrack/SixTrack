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
