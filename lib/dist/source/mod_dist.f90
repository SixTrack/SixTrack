module mod_dist

implicit none

! C Interface
interface

subroutine dist_getrefpara(energy0,mass0,a0,z0) bind(C, name="getrefpara")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE), intent(inout) :: energy0, mass0
  integer(kind=C_INT), intent(inout) :: a0, z0
end subroutine dist_getrefpara


subroutine dist_writefile(fileName, strlen) bind(C, name="writefile_f")
  use, intrinsic :: iso_c_binding
  character(kind=C_CHAR,len=1), intent(in) :: fileName
  integer(kind=C_INT), value, intent(in) :: strlen
end subroutine dist_writefile

subroutine dist_readfile(fileName, strlen) bind(C, name="readfile_f")
  use, intrinsic :: iso_c_binding
  character(kind=C_CHAR,len=1), intent(in) :: fileName
  integer(kind=C_INT), value, intent(in) :: strlen
end subroutine dist_readfile

subroutine dist_initializedistribution(numdist) bind(C, name="initializedistribution")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: numdist
end subroutine dist_initializedistribution

subroutine dist_setphysicalcut(variable, min, max) bind(C, name="setphysicalcut")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: variable
  real(kind=C_DOUBLE), value, intent(in) :: min, max
end subroutine dist_setphysicalcut

subroutine dist_setnormalizedcut(variable, min, max) bind(C, name="setnormalizedcut")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: variable
  real(kind=C_DOUBLE), value, intent(in) :: min, max
end subroutine dist_setnormalizedcut

subroutine dist_setdistribution(distcount) bind(C, name="setdistribution")
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), value, intent(in) :: distcount
end subroutine dist_setdistribution

subroutine dist_addclosedorbit(orbit) bind(C, name="addclosedorbit")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),dimension(6), intent(in) :: orbit
end subroutine dist_addclosedorbit

subroutine dist_get6trackcoord(x,xp,y,yp,sigma,deltap,npart) bind(C, name="get6trackcoord")
  use, intrinsic :: iso_c_binding
  real(kind=C_DOUBLE),dimension(*), intent(inout) :: x,xp,y,yp,sigma,deltap
  integer(kind=C_INT), intent(inout) :: npart
end subroutine dist_get6trackcoord

end interface 
end module