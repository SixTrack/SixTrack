! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-30
!  Read a beam distribution
!
!  Format of the input file:
!    id:         unique identifier of the particle (integer)
!    gen:        parent ID (integer)
!    weight:     statistical weight of the particle (double: >0.0)
!    x,  y,  s:  particle position  [m]
!    xp, yp, zp: particle direction (tangents) []
!    aa, zz:     mass and atomic number
!    m:          rest mass [GeV/c2]
!    pc:         particle momentum [GeV/c]
!    dt:         time delay with respect to the reference particle [s]
!
!    aa,zz and m are now taken into account for hisix!
!
!  NOTA BENE:
!  - id, gen and weight are assigned by the fluka_mod_init subroutine;
!  - z and zp are actually useless (but we never know);
!  - in case the file contains less particle than napx, napx is re-assigned
!
! ================================================================================================ !
module mod_dist

  use crcoall
  use floatPrecision

  implicit none

  logical,            public,  save :: dist_enable     ! DIST input block given
  logical,            public,  save :: dist_echo       ! Echo the read distribution?
  character(len=256), public,  save :: dist_readFile   ! File name for reading the distribution
  character(len=256), public,  save :: dist_echoFile   ! File name for echoing the distribution
  integer,            private, save :: dist_readUnit   ! Unit for reading the distribution
  integer,            private, save :: dist_echoUnit   ! Unit for echoing the distribution

contains

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-30
! ================================================================================================ !
subroutine dist_parseInputLine(inLine, iLine, iErr)

  use string_tools
  use file_units

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "DIST> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) return

  select case(lnSplit(1))

  case("READ")
    if(nSplit < 2) then
      write(lout,"(a)") "DIST> ERROR READ must be followed by one file name only."
      iErr = .true.
      return
    end if
    dist_readFile = trim(lnSplit(2))
    call funit_requestUnit(dist_readFile, dist_readUnit)
    if(.not.dist_enable) dist_enable = .true.

  case("ECHO")
    if(nSplit >= 2) then
      dist_echoFile = trim(lnSplit(2))
    else
      dist_echoFile = "echo_distribution.dat"
    end if
    dist_echo = .true.
    call funit_requestUnit(dist_echoFile, dist_echoUnit)

  case default
    write(lout,"(a)") "DIST> ERROR Unknown keyword '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

end subroutine dist_parseInputLine

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-30
! ================================================================================================ !
subroutine dist_readDist

  use parpro
  use mod_hions
  use mod_common
  use mod_commonmn
  use string_tools
  use mod_particles
  use physical_constants
  use numerical_constants

  implicit none

  integer                 id, gen, j, ln, nSplit
  real(kind=fPrec)        weight, z, zp, dt(npart)
  logical                 spErr, cErr
  character(len=mInputLn) inLine
  character(len=:), allocatable :: lnSplit(:)

  write(lout,"(a)") "DIST> Reading particles from '"//trim(dist_readFile)//"'"

  xv1(:)  = zero
  yv1(:)  = zero
  xv2(:)  = zero
  yv2(:)  = zero
  sigmv(:) = zero
  ejfv(:)  = zero
  naa(:)   = 0
  nzz(:)   = 0
  nucm(:)  = zero
  dt(:)    = zero

  j    = 0
  ln   = 0
  cErr = .false.

  open(unit=dist_readUnit, file=dist_readFile)

10 continue
  read(dist_readUnit,"(a)",end=30,err=20) inLine
  ln = ln+1

  if(inLine(1:1) == "*") goto 10
  if(inLine(1:1) == "#") goto 10
  if(inLine(1:1) == "!") goto 10
  j = j+1

  if(j > napx) then
    write(lout,"(a,i0,a)") "DIST> Stopping reading file as ",napx," particles have been read, as requested in fort.3"
    j = napx
    goto 30
  end if

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) goto 20
  if(nSplit > 0)  call chr_cast(lnSplit(1),  id,      cErr)
  if(nSplit > 1)  call chr_cast(lnSplit(2),  gen,     cErr)
  if(nSplit > 2)  call chr_cast(lnSplit(3),  weight,  cErr)
  if(nSplit > 3)  call chr_cast(lnSplit(4),  xv1(j),  cErr)
  if(nSplit > 4)  call chr_cast(lnSplit(5),  xv2(j),  cErr)
  if(nSplit > 5)  call chr_cast(lnSplit(6),  z,       cErr)
  if(nSplit > 6)  call chr_cast(lnSplit(7),  yv1(j),  cErr)
  if(nSplit > 7)  call chr_cast(lnSplit(8),  yv2(j),  cErr)
  if(nSplit > 8)  call chr_cast(lnSplit(9),  zp,      cErr)
  if(nSplit > 9)  call chr_cast(lnSplit(10), naa(j),  cErr)
  if(nSplit > 10) call chr_cast(lnSplit(11), nzz(j),  cErr)
  if(nSplit > 11) call chr_cast(lnSplit(12), nucm(j), cErr)
  if(nSplit > 12) call chr_cast(lnSplit(13), ejfv(j), cErr)
  if(nSplit > 13) call chr_cast(lnSplit(14), dt(j),   cErr)
  if(cErr) goto 20

  xv1(j)    = xv1(j)*c1e3
  xv2(j)    = xv2(j)*c1e3
  yv1(j)    = yv1(j)*c1e3
  yv2(j)    = yv2(j)*c1e3
  ejfv(j)   = ejfv(j)*c1e3
  nucm(j)   = nucm(j)*c1e3
  sigmv(j)  = -(e0f/e0)*((dt(j)*clight)*c1e3)
  mtc(j)    = (nzz(j)*nucm0)/(zz0*nucm(j))
  nlostp(j) = j
  pstop(j)  = .false.

  goto 10

20 continue
  write(lout,"(a,i0)") "DIST> ERROR Reading particles from line ",ln
  call prror(-1)
  return

30 continue
  if(j == 0) then
    write(lout,"(a)") "DIST> ERROR Reading particles. No particles read from file."
    call prror(-1)
    return
  end if

  close(dist_readUnit)
  write(lout,"(a,i0,a)") "DIST> Read ",j," particles from file '"//trim(dist_readFile)//"'"

  ! Update longitudinal particle arrays from read momentum
  call part_updatePartEnergy(2)

  if(j < napx) then
    write(lout,"(a,i0)") "DIST> WARNING Read a number of particles LOWER than requested: ",napx
    napx = j
  end if

end subroutine dist_readDist

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-31
! ================================================================================================ !
subroutine dist_finaliseDist

  use parpro
  use mod_hions
  use mod_common
  use mod_commont
  use mod_commonmn
  use numerical_constants

  implicit none

  integer          :: j
  real(kind=fPrec) :: chkP, chkE

  ! Check existence of on-momentum particles in the distribution
  do j=1, napx
    chkP = (ejfv(j)/nucm(j))/(e0f/nucm0)-one
    chkE = (ejv(j)/nucm(j))/(e0/nucm0)-one
    if(abs(chkP) < c1m15 .or. abs(chkE) < c1m15) then
      write(lout,"(a)")                "DIST> WARNING Encountered on-momentum particle."
      write(lout,"(a,4(1x,a25))")      "DIST>           ","momentum [MeV/c]","total energy [MeV]","Dp/p","1/(1+Dp/p)"
      write(lout,"(a,4(1x,1pe25.18))") "DIST> ORIGINAL: ", ejfv(j), ejv(j), dpsv(j), oidpsv(j)

      ejfv(j)     = e0f*(nucm(j)/nucm0)
      ejv(j)      = sqrt(ejfv(j)**2+nucm(j)**2)
      dpsv1(j)    = zero
      dpsv(j)     = zero
      oidpsv(j)   = one
      moidpsv(j)  = mtc(j)
      omoidpsv(j) = c1e3*(one-mtc(j))

      if(abs(nucm(j)/nucm0-one) < c1m15) then
        nucm(j) = nucm0
        if(nzz(j) == zz0 .or. naa(j) == aa0) then
          naa(j) = aa0
          nzz(j) = zz0
          mtc(j) = one
        else
          write(lout,"(a)") "DIST> ERROR Mass and/or charge mismatch with relation to sync particle"
          call prror(-1)
        end if
      end if

      write(lout,"(a,4(1x,1pe25.18))") "DIST> CORRECTED:", ejfv(j), ejv(j), dpsv(j), oidpsv(j)
    end if
  end do

  write(lout,"(a,2(1x,i0),1x,f15.7)") "DIST> Reference particle species [A,Z,M]:", aa0, zz0, nucm0
  write(lout,"(a,1x,f15.7)")       "DIST> Reference energy [Z TeV]:", c1m6*e0/zz0

  do j=napx+1,npart
    nlostp(j)   = j
    pstop(j)    = .true.
    ejv(j)      = zero
    dpsv(j)     = zero
    oidpsv(j)   = one
    mtc(j)      = one
    naa(j)      = aa0
    nzz(j)      = zz0
    nucm(j)     = nucm0
    moidpsv(j)  = one
    omoidpsv(j) = zero
  end do

end subroutine dist_finaliseDist

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-30
! ================================================================================================ !
subroutine dist_echoDist

  use mod_common
  use mod_commonmn

  integer j

  open(unit=dist_echoUnit, file=dist_echoFile)
  rewind(dist_echoUnit)
  write(dist_echoUnit,"(a,1pe25.18)") "# Total energy of synch part [MeV]: ",e0
  write(dist_echoUnit,"(a,1pe25.18)") "# Momentum of synch part [MeV/c]:   ",e0f
  write(dist_echoUnit,"(a)")          "#"
  write(dist_echoUnit,"(a)")          "# x[mm], y[mm], xp[mrad], yp[mrad], sigmv[mm], ejfv[MeV/c]"
  do j=1, napx
    write(dist_echoUnit,"(6(1x,1pe25.18))") xv1(j), yv1(j), xv2(j), yv2(j), sigmv(j), ejfv(j)
  end do
  close(dist_echoUnit)

end subroutine dist_echoDist

end module mod_dist
