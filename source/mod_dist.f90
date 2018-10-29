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

  use, intrinsic :: iso_fortran_env, only : int16

  use numerical_constants, only : zero
  use physical_constants,  only : clight
  use parpro,              only : mInputLn, npart
  use mod_hions,           only : naa, nzz, nucm
  use mod_common,          only : napx, e0
  use mod_commonmn,        only : e0f, xv, yv, ejfv, sigmv

  implicit none

! interface variables:
  real(kind=fPrec) :: s(npart)  !(npart)
  real(kind=fPrec) :: pc(npart) !(npart)
  real(kind=fPrec) m0

! temporary variables:
  integer id, gen
  real(kind=fPrec) :: weight, z, zp
  real(kind=fPrec) :: dt(npart)

  integer jj
  character(len=mInputLn) tmp_line

  character, parameter :: comment_char = '*'

  write(lout,"(a)") "DIST> Reading particles from '"//trim(dist_readFile)//"'"
  
  ! Zero the arrays
  xv(:,:)  = zero
  yv(:,:)  = zero
  sigmv(:) = zero
  ejfv(:)  = zero
  naa(:)   = 0
  nzz(:)   = 0
  nucm(:)  = zero

! initialise particle counter
  jj = 0

  open( unit=dist_readUnit, file=dist_readFile )

! cycle on lines in file:
1981 continue
  read(dist_readUnit,"(a)",end=1983,err=1982) tmp_line
  if( tmp_line(1:1).eq.comment_char ) goto 1981
  jj = jj+1

  if( jj.gt.napx ) then
    write(lout,"(a,i0,a)") "DIST> Stopping reading file as ",napx," particles have been as requested in fort.3"
    jj = napx
    goto 1983
  end if

  read(tmp_line, *, err=1982) id,gen,weight,xv(1,jj),xv(2,jj),z,yv(1,jj),yv(2,jj),zp,ejfv(jj),dt(jj),naa(jj),nzz(jj),nucm(jj)
  sigmv(jj) = -e0f/e0 * dt(jj)*clight

  goto 1981

! error while parsing file:
1982 continue
  write(lout,"(a)") "DIST> ERROR Reading particles from line: '"//trim(tmp_line)//"'"
  call prror(-1)
  return

1983 continue
  if( jj.eq.0 ) then
    write(lout,"(a)") "DIST> ERROR Reading particles. No particles read from file."
    call prror(-1)
    return
  end if

  close(dist_readUnit)
  write(lout,"(a,i0)") "DIST> Number of particles read = ",jj

  if( jj.lt.napx ) then
    write(lout,"(a,i0)") "DIST> WARNING Read a number of particles LOWER than the one requested for tracking. Requested: ",napx
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
