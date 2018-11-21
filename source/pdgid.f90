
module mod_pdgid

use, intrinsic :: iso_fortran_env

contains

! Calculates the  PDG id number given A + Z
subroutine CalculatePDGid(id, a, z)

implicit none

integer(kind=int32), intent(out) :: id
integer(kind=int16), intent(in) :: a
integer(kind=int16), intent(in) :: z


!Nuclear codes are given as 10-digit numbers ±10LZZZAAAI.

  if(a .eq. 1 .and. z .eq. 1) then
    id = 2212
    goto 10
  end if

  id = 1000000000 + (a * 10) + (z * 10000)

  10 continue

end subroutine CalculatePDGid

! Calculates the A + Z given a PDG id number
subroutine CalculateAZ(id, a, z)

implicit none

integer(kind=int32), intent(in) :: id
integer(kind=int16), intent(out) :: a
integer(kind=int16), intent(out) :: z


end subroutine CalculateAZ


end module mod_pdgid
