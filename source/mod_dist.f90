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
  logical,            public,  save :: dist_hasFormat  ! Whether the format flag is set
  character(len=256), public,  save :: dist_readFile   ! File name for reading the distribution
  character(len=256), public,  save :: dist_echoFile   ! File name for echoing the distribution
  integer,            private, save :: dist_readUnit   ! Unit for reading the distribution
  integer,            private, save :: dist_echoUnit   ! Unit for echoing the distribution

  integer, allocatable, private, save :: dist_colFormat(:) ! The format of the file columns
  integer,              private, save :: dist_nColumns = 0 ! The number of columns in the file

  !  Column formats
  ! ================
  integer, parameter :: dist_fmtNONE        = 0
  integer, parameter :: dist_fmtPartID      = 1
  integer, parameter :: dist_fmtParentID    = 2

  ! Physical Coordinates
  integer, parameter :: dist_fmtX           = 11 ! Horizontal position
  integer, parameter :: dist_fmtY           = 12 ! Vertical positiom
  integer, parameter :: dist_fmtXP          = 13 ! Horizontal angle
  integer, parameter :: dist_fmtYP          = 14 ! Vertical angle
  integer, parameter :: dist_fmtPX          = 15 ! Horizontal momentum
  integer, parameter :: dist_fmtPY          = 16 ! Vertical momentum
  integer, parameter :: dist_fmtSIGMA       = 17 ! Longitudinal relative position
  integer, parameter :: dist_fmtDT          = 18 ! Time delay
  integer, parameter :: dist_fmtE           = 19 ! Particle energy
  integer, parameter :: dist_fmtP           = 20 ! Particle momentum
  integer, parameter :: dist_fmtDEE0        = 21 ! Relative particle energy (to reference particle)
  integer, parameter :: dist_fmtDPP0        = 22 ! Relative particle momentum (to reference particle)

  ! Normalised Coordinates
  integer, parameter :: dist_fmtX_NORM      = 31 ! Normalised horizontal position
  integer, parameter :: dist_fmtY_NORM      = 32 ! Normalised vertical positiom
  integer, parameter :: dist_fmtXP_NORM     = 33 ! Normalised horizontal angle
  integer, parameter :: dist_fmtYP_NORM     = 34 ! Normalised vertical angle
  integer, parameter :: dist_fmtPX_NORM     = 35 ! Normalised horizontal momentum
  integer, parameter :: dist_fmtPY_NORM     = 36 ! Normalised vertical momentum
  integer, parameter :: dist_fmtSIGMA_NORM  = 37 ! Normalised longitudinal relative position
  integer, parameter :: dist_fmtDT_NORM     = 38 ! Normalised time delay
  integer, parameter :: dist_fmtE_NORM      = 39 ! Normalised particle energy
  integer, parameter :: dist_fmtP_NORM      = 40 ! Normalised particle momentum
  integer, parameter :: dist_fmtDEE0_NORM   = 41 ! Normalised relative particle energy (to reference particle)
  integer, parameter :: dist_fmtDPP0_NORM   = 42 ! Normalised relative particle momentum (to reference particle)

  ! Ion Columns
  integer, parameter :: dist_fmtMASS        = 51 ! Particle mass
  integer, parameter :: dist_fmtCHARGE      = 52 ! Particle Charge
  integer, parameter :: dist_fmtIonA        = 53 ! Ion atomic number
  integer, parameter :: dist_fmtIonZ        = 54 ! Ion atomic charge
  integer, parameter :: dist_fmtPDGID       = 55 ! Particle PDG ID

contains

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-30
! ================================================================================================ !
subroutine dist_parseInputLine(inLine, iLine, iErr)

  use string_tools
  use mod_units

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, i
  logical spErr, cErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "DIST> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(lnSplit(1))

  case("FORMAT")
    if(nSplit < 2) then
      write(lerr,"(a,i0)") "DIST> ERROR FORMAT takes at least 1 argument, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       FORMAT [list of columns]"
      iErr = .true.
      return
    end if
    do i=2,nSplit
      call dist_setColumnFormat(lnSplit(i),cErr)
      if(cErr) then
        iErr = .true.
        return
      end if
    end do

  case("READ")
    if(nSplit < 2) then
      write(lerr,"(a,i0)") "DIST> ERROR READ takes 1 argument, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       READ filename"
      iErr = .true.
      return
    end if
    dist_readFile = trim(lnSplit(2))
    call f_requestUnit(dist_readFile, dist_readUnit)
    if(.not.dist_enable) dist_enable = .true.

  case("ECHO")
    if(nSplit >= 2) then
      dist_echoFile = trim(lnSplit(2))
    else
      dist_echoFile = "echo_distribution.dat"
    end if
    dist_echo = .true.
    call f_requestUnit(dist_echoFile, dist_echoUnit)

  case default
    write(lerr,"(a)") "DIST> ERROR Unknown keyword '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

end subroutine dist_parseInputLine

subroutine dist_setColumnFormat(fmtName, fErr)

  use crcoall
  use string_tools

  character(len=*), intent(in)  :: fmtName
  logical,          intent(out) :: fErr

  fErr = .false.

  select case(chr_toUpper(fmtName))

  case("NONE")
    call dist_appendFormat(dist_fmtNONE)        ! Ignored
  case("ID")
    call dist_appendFormat(dist_fmtPartID)      ! Particle ID
  case("PARENT")
    call dist_appendFormat(dist_fmtParentID)    ! Parent ID

  case("X")
    call dist_appendFormat(dist_fmtX)           ! Horizontal position
  case("Y")
    call dist_appendFormat(dist_fmtY)           ! Vertical positiom
  case("XP")
    call dist_appendFormat(dist_fmtXP)          ! Horizontal angle
  case("YP")
    call dist_appendFormat(dist_fmtYP)          ! Vertical angle
  case("PX")
    call dist_appendFormat(dist_fmtPX)          ! Horizontal momentum
  case("PY")
    call dist_appendFormat(dist_fmtPY)          ! Vertical momentum
  case("SIGMA","DS")
    call dist_appendFormat(dist_fmtSIGMA)       ! Longitudinal relative position
  case("DT")
    call dist_appendFormat(dist_fmtDT)          ! Time delay
  case("E")
    call dist_appendFormat(dist_fmtE)           ! Particle energy
  case("P")
    call dist_appendFormat(dist_fmtP)           ! Particle momentum
  case("DE/E0")
    call dist_appendFormat(dist_fmtDEE0)        ! Relative particle energy (to reference particle)
  case("DP/P0")
    call dist_appendFormat(dist_fmtDPP0)        ! Relative particle momentum (to reference particle)

  case("X_NORM")
    call dist_appendFormat(dist_fmtX_NORM)      ! Normalised horizontal position
  case("Y_NORM")
    call dist_appendFormat(dist_fmtY_NORM)      ! Normalised vertical positiom
  case("XP_NORM")
    call dist_appendFormat(dist_fmtXP_NORM)     ! Normalised horizontal angle
  case("YP_NORM")
    call dist_appendFormat(dist_fmtYP_NORM)     ! Normalised vertical angle
  case("PX_NORM")
    call dist_appendFormat(dist_fmtPX_NORM)     ! Normalised horizontal momentum
  case("PY_NORM")
    call dist_appendFormat(dist_fmtPY_NORM)     ! Normalised vertical momentum
  case("SIGMA_NORM","DS_NORM")
    call dist_appendFormat(dist_fmtSIGMA_NORM)  ! Normalised longitudinal relative position
  case("DT_NORM")
    call dist_appendFormat(dist_fmtDT_NORM)     ! Normalised time delay
  case("E_NORM")
    call dist_appendFormat(dist_fmtE_NORM)      ! Normalised particle energy
  case("P_NORM")
    call dist_appendFormat(dist_fmtP_NORM)      ! Normalised particle momentum
  case("DE/E0_NORM")
    call dist_appendFormat(dist_fmtDEE0_NORM)   ! Normalised relative particle energy (to reference particle)
  case("DP/P0_NORM")
    call dist_appendFormat(dist_fmtDPP0_NORM)   ! Normalised relative particle momentum (to reference particle)

  case("MASS","M")
    call dist_appendFormat(dist_fmtMASS)        ! Particle mass
  case("CHARGE","Q")
    call dist_appendFormat(dist_fmtCHARGE)      ! Particle charge
  case("ION_A")
    call dist_appendFormat(dist_fmtIonA)        ! Ion atomic number
  case("ION_Z")
    call dist_appendFormat(dist_fmtIonZ)        ! Ion atomic charge
  case("PDGID")
    call dist_appendFormat(dist_fmtPDGID)       ! Particle PDG ID

  case("DEFAULT_4D") ! 4D default coordinates
    call dist_appendFormat(dist_fmtX)           ! Horizontal position
    call dist_appendFormat(dist_fmtY)           ! Vertical positiom
    call dist_appendFormat(dist_fmtXP)          ! Horizontal angle
    call dist_appendFormat(dist_fmtYP)          ! Vertical angle

  case("DEFAULT_6D") ! 6D default coordinates
    call dist_appendFormat(dist_fmtX)           ! Horizontal position
    call dist_appendFormat(dist_fmtY)           ! Vertical positiom
    call dist_appendFormat(dist_fmtXP)          ! Horizontal angle
    call dist_appendFormat(dist_fmtYP)          ! Vertical angle
    call dist_appendFormat(dist_fmtSIGMA)       ! Longitudinal relative position
    call dist_appendFormat(dist_fmtDPP0)        ! Relative particle momentum (to reference particle)

  case("DEFAULT_4D_NORM") ! 4D normalised coordinates
    call dist_appendFormat(dist_fmtX_NORM)      ! Normalised horizontal position
    call dist_appendFormat(dist_fmtY_NORM)      ! Normalised vertical positiom
    call dist_appendFormat(dist_fmtXP_NORM)     ! Normalised horizontal angle
    call dist_appendFormat(dist_fmtYP_NORM)     ! Normalised vertical angle

  case("DEFAULT_6D_NORM") ! 6D normalised coordinates
    call dist_appendFormat(dist_fmtX_NORM)      ! Normalised horizontal position
    call dist_appendFormat(dist_fmtY_NORM)      ! Normalised vertical positiom
    call dist_appendFormat(dist_fmtXP_NORM)     ! Normalised horizontal angle
    call dist_appendFormat(dist_fmtYP_NORM)     ! Normalised vertical angle
    call dist_appendFormat(dist_fmtSIGMA_NORM)  ! Normalised longitudinal relative position
    call dist_appendFormat(dist_fmtDPP0_NORM)   ! Normalised relative particle momentum (to reference particle)

  case("DEFAULT_OLD") ! The old DIST block file format
    call dist_appendFormat(dist_fmtPartID)      ! Particle ID
    call dist_appendFormat(dist_fmtNONE)        ! Ignored
    call dist_appendFormat(dist_fmtNONE)        ! Ignored
    call dist_appendFormat(dist_fmtX)           ! Horizontal position
    call dist_appendFormat(dist_fmtY)           ! Vertical positiom
    call dist_appendFormat(dist_fmtNONE)        ! Ignored
    call dist_appendFormat(dist_fmtXP)          ! Horizontal angle
    call dist_appendFormat(dist_fmtYP)          ! Vertical angle
    call dist_appendFormat(dist_fmtNONE)        ! Ignored
    call dist_appendFormat(dist_fmtIonA)        ! Ion atomic number
    call dist_appendFormat(dist_fmtIonZ)        ! Ion atomic charge
    call dist_appendFormat(dist_fmtMASS)        ! Particle mass
    call dist_appendFormat(dist_fmtP)           ! Particle momentum
    call dist_appendFormat(dist_fmtDT)          ! Time delay
  ! call dist_appendFormat(dist_fmtCHARGE)      ! Particle charge
  ! call dist_appendFormat(dist_fmtPDGID)       ! Particle PDG ID

  case default
    write(lout,"(a)") "DIST> ERROR Unknown column format '"//trim(fmtName)//"'"
    fErr = .true.

  end select

end subroutine dist_setColumnFormat

subroutine dist_appendFormat(fmtID)

  use mod_alloc

  integer, intent(in) :: fmtID

  dist_nColumns = dist_nColumns + 1
  call alloc(dist_colFormat,dist_nColumns,dist_fmtNONE,"dist_colFormat")
  dist_colFormat(dist_nColumns) = fmtID

end subroutine dist_appendFormat

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-30
! ================================================================================================ !
subroutine dist_readDist

  use parpro
  use mod_common
  use mod_common_main
  use string_tools
  use mod_particles
  use physical_constants
  use numerical_constants
  use mod_units, only : f_open, f_close

  implicit none

  integer                 id, gen, j, ln, nSplit
  real(kind=fPrec)        weight, z, zp, dt(npart)
  logical                 spErr, cErr
  character(len=mInputLn) inLine
  character(len=:), allocatable :: lnSplit(:)

  write(lout,"(a)") "DIST> Reading particles from '"//trim(dist_readFile)//"'"

  xv1(:)   = zero
  yv1(:)   = zero
  xv2(:)   = zero
  yv2(:)   = zero
  sigmv(:) = zero
  ejfv(:)  = zero
  ejf0v(:) = zero
  naa(:)   = 0
  nzz(:)   = 0
  nucm(:)  = zero
  dt(:)    = zero

  j    = 0
  ln   = 0
  cErr = .false.

  call f_open(unit=dist_readUnit,file=dist_readFile,mode='r',err=cErr,formatted=.true.,status="old")
  if(cErr) goto 19

10 continue
  read(dist_readUnit,"(a)",end=30,err=20) inLine
  ln = ln+1

  if(inLine(1:1) == "*") goto 10
  if(inLine(1:1) == "#") goto 10
  if(inLine(1:1) == "!") goto 10
  j = j+1

  if(j > napx) then
    write(lout,"(a,i0,a)") "DIST> Stopping reading file as ",napx," particles have been read, as requested in "//trim(fort3)
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

  xv1(j)      = xv1(j)*c1e3
  xv2(j)      = xv2(j)*c1e3
  yv1(j)      = yv1(j)*c1e3
  yv2(j)      = yv2(j)*c1e3
  ejfv(j)     = ejfv(j)*c1e3
  nucm(j)     = nucm(j)*c1e3
  sigmv(j)    = -(e0f/e0)*((dt(j)*clight)*c1e3)
  mtc(j)      = (nzz(j)*nucm0)/(zz0*nucm(j))
  partID(j)   = j
  parentID(j) = j
  pstop(j)    = .false.
  ejf0v(j)    = ejfv(j)

  goto 10

19 continue
  write(lerr,"(a)") "DIST> ERROR Opening file '"//trim(dist_readFile)//"'"
  call prror
  return

20 continue
  write(lerr,"(a,i0)") "DIST> ERROR Reading particles from line ",ln
  call prror
  return

30 continue
  if(j == 0) then
    write(lerr,"(a)") "DIST> ERROR Reading particles. No particles read from file."
    call prror
    return
  end if

  call f_close(dist_readUnit)
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
  use mod_common
  use mod_common_track
  use mod_common_main
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
          write(lerr,"(a)") "DIST> ERROR Mass and/or charge mismatch with relation to sync particle"
          call prror
        end if
      end if

      write(lout,"(a,4(1x,1pe25.18))") "DIST> CORRECTED:", ejfv(j), ejv(j), dpsv(j), oidpsv(j)
    end if
  end do

  write(lout,"(a,2(1x,i0),1x,f15.7)") "DIST> Reference particle species [A,Z,M]:", aa0, zz0, nucm0
  write(lout,"(a,1x,f15.7)")       "DIST> Reference energy [Z TeV]:", c1m6*e0/zz0

  do j=napx+1,npart
    partID(j)   = j
    parentID(j) = j
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
  use mod_common_main
  use mod_units, only : f_open, f_close

  integer j
  logical cErr

  call f_open(unit=dist_echoUnit,file=dist_echoFile,mode='w',err=cErr,formatted=.true.)
  if(cErr) goto 19

  rewind(dist_echoUnit)
  write(dist_echoUnit,"(a,1pe25.18)") "# Total energy of synch part [MeV]: ",e0
  write(dist_echoUnit,"(a,1pe25.18)") "# Momentum of synch part [MeV/c]:   ",e0f
  write(dist_echoUnit,"(a)")          "#"
  write(dist_echoUnit,"(a)")          "# x[mm], y[mm], xp[mrad], yp[mrad], sigmv[mm], ejfv[MeV/c]"
  do j=1, napx
    write(dist_echoUnit,"(6(1x,1pe25.18))") xv1(j), yv1(j), xv2(j), yv2(j), sigmv(j), ejfv(j)
  end do
  call f_close(dist_echoUnit)

  return

19 continue
  write(lerr,"(a)") "DIST> ERROR Opening file '"//trim(dist_echoFile)//"'"
  call prror
  return

end subroutine dist_echoDist

end module mod_dist
