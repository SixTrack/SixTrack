
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

integer(kind=int16) :: z1,z2,z3
integer(kind=int16) :: a1,a2,a3
integer(kind=int16) :: tmpid 

  tmpid = id - 1000000000;

  tmpid = tmpid / 10;
  a1 = mod(tmpid,10_int16);

  tmpid = tmpid / 10;
  a2 = mod(tmpid,10_int16);

  tmpid = tmpid / 10;
  a3 = mod(tmpid,10_int16);

  a = a1 + 10*a2 + 100*a3;

  tmpid = tmpid / 10;
  z1 = mod(tmpid,10_int16);

  tmpid = tmpid / 10;
  z2 = mod(tmpid,10_int16);

  tmpid = tmpid / 10;
  z3 = mod(tmpid,10_int16);

  z = z1 + 10*z2 + 100*z3;

end subroutine CalculateAZ

!Learn to speak "FLUKA"
subroutine GetFLUKAid_fromPDG(pdgid, fluka_id)

  use crcoall

implicit none

integer(kind=int32), intent(in)  :: pdgid
integer(kind=int32), intent(out) :: fluka_id

integer(kind=int16) :: A
integer(kind=int16) :: Z

!Only way is to do an A->B translation using the numbers in the FLUKA manual
!FM.pdf "Table 5.1: Fluka particle names and code numbers"

!proton
if(pdgid .eq. 2212) then
  fluka_id = 1

!anti proton
else if(pdgid .eq. -2212) then
  fluka_id = 2

!electron
else if(pdgid .eq. 11) then
  fluka_id = 3

!positron
else if(pdgid .eq. -11) then
  fluka_id = 4

!5 to 9 are neutrals

!mu+
else if(pdgid .eq. -13) then
  fluka_id = 10

!mu-
else if(pdgid .eq. 13) then
  fluka_id = 11

!pi+
else if(pdgid .eq. 211) then
  fluka_id = 13

!pi-
else if(pdgid .eq. -211) then
  fluka_id = 14

!k+
else if(pdgid .eq. 321) then
  fluka_id = 15

!k-
else if(pdgid .eq. -321) then
  fluka_id = 16

!lambda
else if(pdgid .eq. 3122) then
  fluka_id = 17

!anti lambda
else if(pdgid .eq. -3122) then
  fluka_id = 18

!sigma-
else if(pdgid .eq. 3112) then
  fluka_id = 20

!sigma+
else if(pdgid .eq. 3222) then
  fluka_id = 21

!anti-sigma-
else if(pdgid .eq. -3222) then
  fluka_id = 31

!anti-sigma+
else if(pdgid .eq. -3112) then
  fluka_id = 33

!Xi-
else if(pdgid .eq. 3312) then
  fluka_id = 36

!Xi+
else if(pdgid .eq. -3312) then
  fluka_id = 37

!Omega-
else if(pdgid .eq. 3334) then
  fluka_id = 38

!anti omega
else if(pdgid .eq. -3334) then
  fluka_id = 39


!tau+
else if(pdgid .eq. -15) then
  fluka_id = 41

!tau-
else if(pdgid .eq. 15) then
  fluka_id = 42

!ion
else if(pdgid .gt. 1000000000) then
!get A,Z
  call CalculateAZ(pdgid, A, Z)

  if(Z .eq. 1) then

!Protons, etc
    if(A .eq. 1) then
      fluka_id = 1
    else if(A .eq. 2) then
      fluka_id = -3
    else if(A .eq. 3) then
      fluka_id = -4
    else
      fluka_id = -2
    end if

!Helium, etc
  else if(Z .eq. 2) then
    if(A .eq. 3) then
      fluka_id = -5
    else if(A .eq. 4) then
      fluka_id = -6
    else 
      fluka_id = -2
    end if

  else
!Generic heavy ion
    fluka_id = -2
  end if
!something else
else
  write(lout,'(a,i12,a)') 'Unknown particle ID in PDG id to FLUKA conversion: ', pdgid, '  - exiting!'
  call prror(-1)
end if

return

end subroutine GetFLUKAid_fromPDG


!FLUKA -> PDGid conversion
subroutine GetPDGid_fromFLUKA(fluka_id, pdg_id, A, Z)

  use crcoall

implicit none

integer(kind=int32), intent(in)  :: fluka_id
integer(kind=int32), intent(out) :: pdg_id

integer(kind=int16), intent(in) :: A
integer(kind=int16), intent(in) :: Z

!Only way is to do an A->B translation using the numbers in the FLUKA manual
!FM.pdf "Table 5.1: Fluka particle names and code numbers"

!proton
if(fluka_id .eq. 1) then
  pdg_id = 2212

!anti proton
else if(fluka_id .eq. 2) then
  pdg_id = -2212

!electron
else if(fluka_id .eq. 3) then
  pdg_id = 11

!positron
else if(fluka_id .eq. 4) then
  pdg_id = -11

!5 to 9 are neutrals

!mu+
else if(fluka_id .eq. 10) then
  pdg_id = -13

!mu-
else if(fluka_id .eq. 11) then
  pdg_id = 13

!pi+
else if(fluka_id .eq. 13) then
  pdg_id = 211

!pi-
else if(fluka_id .eq. 14) then
  pdg_id = -211

!k+
else if(fluka_id .eq. 15) then
  pdg_id = 321

!k-
else if(fluka_id .eq. 16) then
  pdg_id = -321

!lambda
else if(fluka_id .eq. 17) then
  pdg_id = 3122

!anti lambda
else if(fluka_id .eq. 18) then
  pdg_id = -3122

!sigma-
else if(fluka_id .eq. 20) then
  pdg_id = 3112

!sigma+
else if(fluka_id .eq. 21) then
  pdg_id = 3222

!anti-sigma-
else if(fluka_id .eq. 31) then
  pdg_id = -3222

!anti-sigma+
else if(fluka_id .eq. 33) then
  pdg_id = -3112

!Xi-
else if(fluka_id .eq. 36) then
  pdg_id = 3312

!Xi+
else if(fluka_id .eq. 37) then
  pdg_id = -3312

!Omega-
else if(fluka_id .eq. 38) then
  pdg_id = 3334

!anti omega
else if(fluka_id .eq. 39) then
  pdg_id = -3334

!tau+
else if(fluka_id .eq. 41) then
  pdg_id = -15

!tau-
else if(fluka_id .eq. 42) then
  pdg_id = 15

!ion
else if(fluka_id .eq. -2 .or. fluka_id .eq. -3 .or. fluka_id .eq. -4 .or. fluka_id .eq. -5 .or. fluka_id .eq. -6) then
  call CalculatePDGid(pdg_id, A, Z)
!something else
else
  write(lout,'(a,i12,a)') 'Unknown particle ID in FLUKA to PDG id conversion: ', fluka_id, '  - exiting!'
  call prror(-1)
end if

return
end subroutine GetPDGid_fromFLUKA

end module mod_pdgid
