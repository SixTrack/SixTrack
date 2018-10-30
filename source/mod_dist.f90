! ================================================================================================ !
!  A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!  last modified: 2018-06-02
!  read a beam distribution
!  always in main code
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
!  - aa,zz and m are now taken into account for hisix
!
!  NOTA BENE:
!  - id, gen and weight are assigned by the fluka_mod_init subroutine;
!  - z and zp are actually useless (but we never know);
!  - aa, zz and m are not stored at the moment (safer decision from the code
!    point of view, until a decision about ion tracking is taken);
!    the subroutine fluka_send is then responsible for using the corresponding
!    values for protons through the interface, whereas the subroutine fluka_receive
!    simply ignores the values passed through the FlukaIO interface;
!
!  variables in input to routine:
!  - napx: number of protons to be tracked (from fort.3 file);
!  - npart: max number of protons that can be tracked (array dimensioning);
!  - enom: nominal total energy of the beam (ie of synch particle) [MeV];
!  - pnom: nominal linear momentum of the beam (ie of synch particle) [MeV/c];
!  - clight: speed of light [m/s];
!  NB: in case the file contains less particle than napx, napx is
!      re-assigned
!
!  output variables:
!    all other variables in the interface (6D tracking variables);
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
      write(lout,"(a)") "DIST> ERROR READ must be followed by a one file name only."
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

  end select

end subroutine dist_parseInputLine

subroutine dist_readDist()

  use numerical_constants, only : zero, c1e3
  use physical_constants,  only : clight
  use parpro,              only : mInputLn, npart
  use mod_hions,           only : naa, nzz, nucm
  use mod_common,          only : napx, e0
  use mod_commonmn,        only : e0f, xv, yv, ejfv, sigmv
  use string_tools

  implicit none

  integer                 id, gen, jj, ln, nSplit
  real(kind=fPrec)        weight, z, zp, dt(npart)
  logical                 spErr, cErr
  character(len=mInputLn) inLine
  character(len=:), allocatable :: lnSplit(:)

  write(lout,"(a)") "DIST> Reading particles from '"//trim(dist_readFile)//"'"

  xv(:,:)  = zero
  yv(:,:)  = zero
  sigmv(:) = zero
  ejfv(:)  = zero
  naa(:)   = 0
  nzz(:)   = 0
  nucm(:)  = zero
  dt(:)    = zero

  jj   = 0
  ln   = 0
  cErr = .false.

  open(unit=dist_readUnit, file=dist_readFile)

10 continue
  read(dist_readUnit,"(a)",end=30,err=20) inLine
  ln = ln+1

  if(inLine(1:1) == "*") goto 10
  if(inLine(1:1) == "#") goto 10
  if(inLine(1:1) == "!") goto 10
  jj = jj+1

  if(jj > napx) then
    write(lout,"(a,i0,a)") "DIST> Stopping reading file as ",napx," particles have been read, as requested in fort.3"
    jj = napx
    goto 30
  end if

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) goto 20
  if(nSplit > 0)  call chr_cast(lnSplit(1),  id,       cErr)
  if(nSplit > 1)  call chr_cast(lnSplit(2),  gen,      cErr)
  if(nSplit > 2)  call chr_cast(lnSplit(3),  weight,   cErr)
  if(nSplit > 3)  call chr_cast(lnSplit(4),  xv(1,jj), cErr)
  if(nSplit > 4)  call chr_cast(lnSplit(5),  xv(2,jj), cErr)
  if(nSplit > 5)  call chr_cast(lnSplit(6),  z,        cErr)
  if(nSplit > 6)  call chr_cast(lnSplit(7),  yv(1,jj), cErr)
  if(nSplit > 7)  call chr_cast(lnSplit(8),  yv(2,jj), cErr)
  if(nSplit > 8)  call chr_cast(lnSplit(9),  zp,       cErr)
  if(nSplit > 9)  call chr_cast(lnSplit(10), naa(jj),  cErr)
  if(nSplit > 10) call chr_cast(lnSplit(11), nzz(jj),  cErr)
  if(nSplit > 11) call chr_cast(lnSplit(12), nucm(jj), cErr)
  if(nSplit > 12) call chr_cast(lnSplit(13), ejfv(jj), cErr)
  if(nSplit > 13) call chr_cast(lnSplit(14), dt(jj),   cErr)
  if(cErr) goto 20

  xv(1,jj)  = xv(1,jj)*c1e3
  xv(2,jj)  = xv(2,jj)*c1e3
  yv(1,jj)  = yv(1,jj)*c1e3
  yv(2,jj)  = yv(2,jj)*c1e3
  ejfv(jj)  = ejfv(jj)*c1e3
  nucm(jj)  = nucm(jj)*c1e3
  sigmv(jj) = -(e0f/e0)*((dt(jj)*clight)*c1e3)

  goto 10

20 continue
  write(lout,"(a,i0)") "DIST> ERROR Reading particles from line ",ln
  call prror(-1)
  return

30 continue
  if(jj == 0) then
    write(lout,"(a)") "DIST> ERROR Reading particles. No particles read from file."
    call prror(-1)
    return
  end if

  close(dist_readUnit)
  write(lout,"(a,i0,a)") "DIST> Read ",jj," particles from file '"//trim(dist_readFile)//"'"

  if(jj < napx) then
    write(lout,"(a,i0)") "DIST> WARNING Read a number of particles LOWER than requested: ",napx
    napx = jj
  end if

end subroutine dist_readDist

subroutine dist_echoDist()

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
    write(dist_echoUnit,"(6(1x,1pe25.18))") xv(1,j), yv(1,j), xv(2,j), yv(2,j), sigmv(j), ejfv(j)
  end do
  close(dist_echoUnit)

end subroutine dist_echoDist

end module mod_dist
