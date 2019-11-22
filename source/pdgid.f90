
module mod_pdgid

  implicit none

contains

! Calculates the  PDG id number given A + Z
subroutine CalculatePDGid(id, a, z)

  use, intrinsic :: iso_fortran_env

  integer(kind=int32), intent(out) :: id
  integer(kind=int16), intent(in) :: a
  integer(kind=int16), intent(in) :: z

  !Nuclear codes are given as 10-digit numbers ±10LZZZAAAI.

  if(a == 1 .and. z == 1) then
    id = 2212
    return
  end if
  id = 1000000000 + int(a,int32)*10 + int(z,int32)*10000

end subroutine CalculatePDGid

! Calculates the A + Z given a PDG id number
subroutine CalculateAZ(id, a, z)

  use, intrinsic :: iso_fortran_env

  integer(kind=int32), intent(in) :: id
  integer(kind=int16), intent(out) :: a
  integer(kind=int16), intent(out) :: z

  integer(kind=int16) :: z1,z2,z3
  integer(kind=int16) :: a1,a2,a3
  integer(kind=int32) :: tmpid

  if(id == 2212) then
    a = 1
    z = 1
    return
  end if

  if(id < 1000000000) then
    a = 0
    z = 0
    return
  end if

  tmpid = id - 1000000000

  tmpid = tmpid / 10
  a1 = mod(tmpid,10_int32)

  tmpid = tmpid / 10
  a2 = mod(tmpid,10_int32)

  tmpid = tmpid / 10
  a3 = mod(tmpid,10_int32)

  a = a1 + 10*a2 + 100*a3

  tmpid = tmpid / 10
  z1 = mod(tmpid,10_int32)

  tmpid = tmpid / 10
  z2 = mod(tmpid,10_int32)

  tmpid = tmpid / 10
  z3 = mod(tmpid,10_int32)

  z = z1 + 10*z2 + 100*z3

end subroutine CalculateAZ

!Learn to speak "FLUKA"
subroutine GetFLUKAid_fromPDG(pdgid, fluka_id)

  use, intrinsic :: iso_fortran_env

  use crcoall

  integer(kind=int32), intent(in)  :: pdgid
  integer(kind=int32), intent(out) :: fluka_id

  integer(kind=int16) :: A
  integer(kind=int16) :: Z

  !Only way is to do an A->B translation using the numbers in the FLUKA manual
  !FM.pdf "Table 5.1: Fluka particle names and code numbers"

  if(pdgid == 2212) then
    !proton
    fluka_id = 1

  else if(pdgid == -2212) then
    !anti proton
    fluka_id = 2

  else if(pdgid == 11) then
    !electron
    fluka_id = 3

  else if(pdgid == -11) then
    !positron
    fluka_id = 4

  !5 to 9 are neutrals

  else if(pdgid == -13) then
    !mu+
    fluka_id = 10

  else if(pdgid == 13) then
    !mu-
    fluka_id = 11

  else if(pdgid == 211) then
    !pi+
    fluka_id = 13

  else if(pdgid == -211) then
    !pi-
    fluka_id = 14

  else if(pdgid == 321) then
    !k+
    fluka_id = 15

  else if(pdgid == -321) then
    !k-
    fluka_id = 16

  else if(pdgid == 3122) then
    !lambda
    fluka_id = 17

  else if(pdgid == -3122) then
    !anti lambda
    fluka_id = 18

  else if(pdgid == 3112) then
    !sigma-
    fluka_id = 20

  else if(pdgid == 3222) then
    !sigma+
    fluka_id = 21

  else if(pdgid == -3222) then
    !anti-sigma-
    fluka_id = 31

  else if(pdgid == -3112) then
    !anti-sigma+
    fluka_id = 33

  else if(pdgid == 3312) then
    !Xi-
    fluka_id = 36

  else if(pdgid == -3312) then
    !Xi+
    fluka_id = 37

  else if(pdgid == 3334) then
    !Omega-
    fluka_id = 38

  else if(pdgid == -3334) then
    !anti omega
    fluka_id = 39

  else if(pdgid == -15) then
    !tau+
    fluka_id = 41

  else if(pdgid == 15) then
    !tau-
    fluka_id = 42

  else if(pdgid > 1000000000) then
    !ion
    call CalculateAZ(pdgid, A, Z)

    if(Z == 1) then

      if(A == 1) then
        !Protons, etc
        fluka_id = 1
      else if(A == 2) then
        fluka_id = -3
      else if(A == 3) then
        fluka_id = -4
      else
        fluka_id = -2
      end if

    else if(Z == 2) then
      !Helium, etc
      if(A == 3) then
        fluka_id = -5
      else if(A == 4) then
        fluka_id = -6
      else
        fluka_id = -2
      end if

    else
      !Generic heavy ion
      fluka_id = -2
    end if
  else
    !something else
    write(lerr,"(a,i0,a)") "PDGID> ERROR Unknown particle ID ",pdgid," in FLUKA conversion"
    call prror
  end if

end subroutine GetFLUKAid_fromPDG


!FLUKA -> PDGid conversion
subroutine GetPDGid_fromFLUKA(fluka_id, pdg_id, A, Z)

  use, intrinsic :: iso_fortran_env

  use crcoall

  integer(kind=int32), intent(in)  :: fluka_id
  integer(kind=int32), intent(out) :: pdg_id

  integer(kind=int16), intent(in) :: A
  integer(kind=int16), intent(in) :: Z

  !Only way is to do an A->B translation using the numbers in the FLUKA manual
  !FM.pdf "Table 5.1: Fluka particle names and code numbers"

  if(fluka_id == 1) then
    !proton
    pdg_id = 2212

  else if(fluka_id == 2) then
    !anti proton
    pdg_id = -2212

  else if(fluka_id == 3) then
    !electron
    pdg_id = 11

  else if(fluka_id == 4) then
    !positron
    pdg_id = -11

  !5 to 9 are neutrals

  else if(fluka_id == 10) then
    !mu+
    pdg_id = -13

  else if(fluka_id == 11) then
    !mu-
    pdg_id = 13

  else if(fluka_id == 13) then
    !pi+
    pdg_id = 211

  else if(fluka_id == 14) then
    !pi-
    pdg_id = -211

  else if(fluka_id == 15) then
    !k+
    pdg_id = 321

  else if(fluka_id == 16) then
    !k-
    pdg_id = -321

  else if(fluka_id == 17) then
    !lambda
    pdg_id = 3122

  else if(fluka_id == 18) then
    !anti lambda
    pdg_id = -3122

  else if(fluka_id == 20) then
    !sigma-
    pdg_id = 3112

  else if(fluka_id == 21) then
    !sigma+
    pdg_id = 3222

  else if(fluka_id == 31) then
    !anti-sigma-
    pdg_id = -3222

  else if(fluka_id == 33) then
    !anti-sigma+
    pdg_id = -3112

  else if(fluka_id == 36) then
    !Xi-
    pdg_id = 3312

  else if(fluka_id == 37) then
    !Xi+
    pdg_id = -3312

  else if(fluka_id == 38) then
    !Omega-
    pdg_id = 3334

  else if(fluka_id == 39) then
    !anti omega
    pdg_id = -3334

  else if(fluka_id == 41) then
    !tau+
    pdg_id = -15

  else if(fluka_id == 42) then
    !tau-
    pdg_id = 15

  else if(fluka_id == -2 .or. fluka_id == -3 .or. fluka_id == -4 .or. fluka_id == -5 .or. fluka_id == -6) then
    !ion
    call CalculatePDGid(pdg_id, A, Z)

  else
    ! something else
    write(lerr,"(a,i0,a)") "PDGID> ERROR Unknown particle ID ",fluka_id," in PDG ID conversion"
    call prror
  end if

end subroutine GetPDGid_fromFLUKA

end module mod_pdgid
