! ================================================================================================ !
!
!  Beam Distribution Block
! ~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-07-11
!
!  This module was completely rewritten to support DISTlib in July 2019
!
! ================================================================================================ !
module mod_dist

  use crcoall
  use floatPrecision
  use parpro, only : mFileName
  use numerical_constants, only : zero

  implicit none

  logical,                  public,  save :: dist_enable      = .false. ! DIST input block given
  logical,                  private, save :: dist_echo        = .false. ! Echo the read distribution?
  logical,                  private, save :: dist_hasFormat   = .false. ! Whether the FORMAT keyword exists in block or not
  logical,                  private, save :: dist_libRead     = .false. ! Read file with dist library instead of internal reader
  logical,                  private, save :: dist_distLib     = .false. ! DISTlib is needed to generate the distribution requested
  character(len=mFileName), private, save :: dist_distFile    = " "     ! File name for reading the distribution
  character(len=mFileName), private, save :: dist_echoFile    = " "     ! File name for echoing the distribution
  real(kind=fPrec),         public,  save :: dist_beamEmit(3) = zero    ! Beam emittance for generator and normalisation

  integer,          allocatable, private, save :: dist_colFormat(:)   ! The column types in the FORMAT
  real(kind=fPrec), allocatable, private, save :: dist_colScale(:)    ! Scaling factor for the columns
  integer,                       private, save :: dist_nColumns  = 0  ! The number of columns in the FORMAT

  character(len=:), allocatable, private, save :: dist_partLine(:)    ! PARTICLE definitions in the block
  integer,                       private, save :: dist_nParticle = 0  ! Number of PARTICLE keywords in block

  integer,                       private, save :: dist_numPart   = 0  ! Number of actual particles generated or read

  !  Column formats
  ! ================
  integer, parameter :: dist_fmtNONE        = 0  ! Column ignored
  integer, parameter :: dist_fmtPartID      = 1  ! Paricle ID
  integer, parameter :: dist_fmtParentID    = 2  ! Particle parent ID (for secondary particles, otherwise equal particle ID)

  ! Physical Coordinates
  integer, parameter :: dist_fmtX           = 11 ! Horizontal position
  integer, parameter :: dist_fmtY           = 12 ! Vertical positiom
  integer, parameter :: dist_fmtXP          = 13 ! Horizontal momentum ratio Px/|P|
  integer, parameter :: dist_fmtYP          = 14 ! Vertical momentum ratio Py/|P|
  integer, parameter :: dist_fmtPX          = 15 ! Horizontal momentum
  integer, parameter :: dist_fmtPY          = 16 ! Vertical momentum
  integer, parameter :: dist_fmtPXP0        = 17 ! Horizontal momentum Px/P0
  integer, parameter :: dist_fmtPYP0        = 18 ! Vertical momentum Py/P0
  integer, parameter :: dist_fmtZETA        = 19 ! Longitudinal relative position (canonical)
  integer, parameter :: dist_fmtSIGMA       = 20 ! Longitudinal relative position
  integer, parameter :: dist_fmtDT          = 21 ! Time delay
  integer, parameter :: dist_fmtE           = 22 ! Particle energy
  integer, parameter :: dist_fmtP           = 23 ! Particle momentum
  integer, parameter :: dist_fmtDEE0        = 24 ! Relative particle energy (to reference particle)
  integer, parameter :: dist_fmtDPP0        = 25 ! Relative particle momentum (to reference particle)
  integer, parameter :: dist_fmtPT          = 26 ! Delta energy over reference momentum (Pt)

  ! Normalised Coordinates
  integer, parameter :: dist_fmtXN          = 31 ! Normalised horizontal position
  integer, parameter :: dist_fmtYN          = 32 ! Normalised vertical position
  integer, parameter :: dist_fmtZN          = 33 ! Normalised longitudinal position
  integer, parameter :: dist_fmtPXN         = 34 ! Normalised horizontal momentum
  integer, parameter :: dist_fmtPYN         = 35 ! Normalised vertical momentum
  integer, parameter :: dist_fmtPZN         = 36 ! Normalised longitudinal momentum

  integer, parameter :: dist_fmtJX          = 41 ! Horizontal action
  integer, parameter :: dist_fmtJY          = 42 ! Vertical action
  integer, parameter :: dist_fmtJZ          = 43 ! Longitudinal action
  integer, parameter :: dist_fmtPhiX        = 44 ! Horizontal action angle
  integer, parameter :: dist_fmtPhiY        = 45 ! Vertical action angle
  integer, parameter :: dist_fmtPhiZ        = 46 ! Longitudinal action angle

  ! Ion Columns
  integer, parameter :: dist_fmtMASS        = 51 ! Particle mass
  integer, parameter :: dist_fmtCHARGE      = 52 ! Particle Charge
  integer, parameter :: dist_fmtIonA        = 53 ! Ion atomic mass
  integer, parameter :: dist_fmtIonZ        = 54 ! Ion atomic number
  integer, parameter :: dist_fmtPDGID       = 55 ! Particle PDG ID

  ! Flags for columns we've set that we need to track for later checks
  logical, private, save :: dist_readMass   = .false.
  logical, private, save :: dist_readIonZ   = .false.
  logical, private, save :: dist_readIonA   = .false.
  logical, private, save :: dist_readCharge = .false.
  logical, private, save :: dist_readPDGID  = .false.

  ! Temporary particle arrays
  real(kind=fPrec), allocatable, private, save :: dist_partCol1(:) ! Ends up in array xv1
  real(kind=fPrec), allocatable, private, save :: dist_partCol2(:) ! Ends up in array yv1
  real(kind=fPrec), allocatable, private, save :: dist_partCol3(:) ! Ends up in array xv2
  real(kind=fPrec), allocatable, private, save :: dist_partCol4(:) ! Ends up in array yv2
  real(kind=fPrec), allocatable, private, save :: dist_partCol5(:) ! Ends up in array sigmv
  real(kind=fPrec), allocatable, private, save :: dist_partCol6(:) ! Ends up in array dpsv, evj, and ejfv
  integer,                       private, save :: dist_partFmt(6) = dist_fmtNONE ! The format used for each column

#ifdef DISTLIB
  interface

    subroutine distlib_init(numDist) bind(C, name="initializedistribution")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value, intent(in) :: numDist
    end subroutine distlib_init

    subroutine distlib_readFile(fileName, strLen) bind(C, name="readfile_f")
      use, intrinsic :: iso_c_binding
      character(kind=C_CHAR,len=1), intent(in) :: fileName
      integer(kind=C_INT),   value, intent(in) :: strLen
    end subroutine distlib_readFile

    subroutine distlib_setEnergyMass(energy0, mass0) bind(C, name="sete0andmass0")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), value, intent(in) :: energy0
      real(kind=C_DOUBLE), value, intent(in) :: mass0
    end subroutine distlib_setEnergyMass

    subroutine distlib_setEmittance12(emit1, emit2) bind(C, name="setemitt12")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), value, intent(in) :: emit1
      real(kind=C_DOUBLE), value, intent(in) :: emit2
    end subroutine distlib_setEmittance12

    subroutine distlib_setEmittance3(emit3) bind(C, name="setemitt3")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), value, intent(in) :: emit3
    end subroutine distlib_setEmittance3

    subroutine distlib_setTasMatrix(flatTas) bind(C, name="settasmatrix")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), dimension(36), intent(in) :: flatTas
    end subroutine distlib_setTasMatrix

    subroutine distlib_setNormalised(arr1, arr2, arr3, arr4, arr5, arr6, arrLen) bind(C, name="setnormalizedcoords")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), dimension(*), intent(in) :: arr1
      real(kind=C_DOUBLE), dimension(*), intent(in) :: arr2
      real(kind=C_DOUBLE), dimension(*), intent(in) :: arr3
      real(kind=C_DOUBLE), dimension(*), intent(in) :: arr4
      real(kind=C_DOUBLE), dimension(*), intent(in) :: arr5
      real(kind=C_DOUBLE), dimension(*), intent(in) :: arr6
      integer(kind=C_INT),               intent(in) :: arrLen
    end subroutine distlib_setNormalised

    subroutine distlib_getPartCoords(arr1, arr2, arr3, arr4, arr5, arr6, arrLen) bind(C, name="get6trackcoord")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), dimension(*), intent(inout) :: arr1
      real(kind=C_DOUBLE), dimension(*), intent(inout) :: arr2
      real(kind=C_DOUBLE), dimension(*), intent(inout) :: arr3
      real(kind=C_DOUBLE), dimension(*), intent(inout) :: arr4
      real(kind=C_DOUBLE), dimension(*), intent(inout) :: arr5
      real(kind=C_DOUBLE), dimension(*), intent(inout) :: arr6
      integer(kind=C_INT),               intent(inout) :: arrLen
    end subroutine distlib_getPartCoords

  end interface
#endif

contains

! ================================================================================================ !
!  Master Module Subrotuine
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-08
!  Updated: 2019-07-09
!  Called from main_cr if dist_enable is true. This should handle everything needed.
! ================================================================================================ !
subroutine dist_generateDist

  use crcoall
  use mod_alloc
  use mod_common
  use mod_particles
  use mod_common_main
  use numerical_constants

  integer distLibPart

  call alloc(dist_partCol1, napx, zero, "dist_partCol1")
  call alloc(dist_partCol2, napx, zero, "dist_partCol2")
  call alloc(dist_partCol3, napx, zero, "dist_partCol3")
  call alloc(dist_partCol4, napx, zero, "dist_partCol4")
  call alloc(dist_partCol5, napx, zero, "dist_partCol5")
  call alloc(dist_partCol6, napx, zero, "dist_partCol6")

#ifdef DISTLIB
  if(dist_distLib) then
    call distlib_init(1)
    call distlib_setEnergyMass(e0, nucm0)
    call distlib_setEmittance12(dist_beamEmit(1), dist_beamEmit(2))
    call distlib_setEmittance3(dist_beamEmit(3))
    block
      real(kind=fPrec) tmpTas(6,6), flatTas(36)
      tmpTas(1:6,1:6) = tas(1:6,1:6)
      tmpTas(1:5,6)   = tmpTas(1:5,6)*c1m3
      tmpTas(6,1:5)   = tmpTas(6,1:5)*c1e3
      flatTas = pack(transpose(tmpTas),.true.)
      call distlib_setTasMatrix(flatTas)
    end block
  end if
#endif

  if(dist_nParticle > 0) then
    ! Particles are read from DIST block in fort.3
    call dist_parseParticles
    call dist_postprParticles
  elseif(dist_distFile /= " ") then
    if(dist_distLib) then
      ! Particles are read entirely in DISTlib
#ifdef DISTLIB
      call distlib_readFile(trim(dist_distFile), len_trim(dist_distFile))
#endif
    else
      call dist_readDist
      call dist_postprParticles
    end if
  end if

#ifdef DISTLIB
  if(dist_distLib) then
    ! Our final coordinates are taken from DISTlib
    distLibPart = napx
    call distlib_getPartCoords(                    &
      dist_partCol1, dist_partCol2, dist_partCol3, &
      dist_partCol4, dist_partCol5, dist_partCol6, &
      distLibPart                                  &
    )
    if(distLibPart == napx) then
      xv1(1:napx)   = dist_partCol1(1:napx)
      yv1(1:napx)   = dist_partCol2(1:napx)
      xv2(1:napx)   = dist_partCol3(1:napx)
      yv2(1:napx)   = dist_partCol4(1:napx)
      sigmv(1:napx) = dist_partCol5(1:napx)
      dpsv(1:napx)  = dist_partCol6(1:napx)
      call part_updatePartEnergy(3,.false.)
    else
      write(lerr,"(2(a,i0))") "DIST> ERROR DISTlib returned ",distLibPart," particles, while SixTrack expected ",napx
      call prror
    end if
  end if
#endif

  if(dist_numPart > 0) then
    call dist_finaliseDist
    call part_applyClosedOrbit
    if(dist_echo) then
      call dist_echoDist
    end if
  end if

  call dealloc(dist_partCol1, "dist_partCol1")
  call dealloc(dist_partCol2, "dist_partCol2")
  call dealloc(dist_partCol3, "dist_partCol3")
  call dealloc(dist_partCol4, "dist_partCol4")
  call dealloc(dist_partCol5, "dist_partCol5")
  call dealloc(dist_partCol6, "dist_partCol6")

end subroutine dist_generateDist

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-07-09
! ================================================================================================ !
subroutine dist_parseInputLine(inLine, iLine, iErr)

  use parpro
  use mod_alloc
  use mod_units
  use string_tools

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

  cErr = .false.

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
    dist_hasFormat = .true.

  case("READ")
    if(nSplit /= 2 .and. nSplit /= 3) then
      write(lerr,"(a,i0)") "DIST> ERROR READ takes 1 or 2 arguments, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       READ filename [LIBDIST]"
      iErr = .true.
      return
    end if
    dist_distFile = trim(lnSplit(2))
    if(nSplit > 2) call chr_cast(lnSplit(3), dist_libRead, cErr)
    if(dist_libRead) then
      dist_distLib = .true.
    end if

  case("PARTICLE")
    if(dist_hasFormat .eqv. .false.) then
      write(lerr,"(a,i0)") "DIST> ERROR PARTICLE keyword requires a FORMAT to be defined first"
      iErr = .true.
      return
    end if
    if(nSplit /= dist_nColumns + 1) then
      write(lerr,"(a)")       "DIST> ERROR PARTICLE values must match the number of definitions in FORMAT"
      write(lerr,"(2(a,i0))") "DIST>       Got ",nsplit-1," values, but FORMAT defines ",dist_nColumns
      iErr = .true.
      return
    end if
    dist_nParticle = dist_nParticle + 1
    call alloc(dist_partLine, mInputLn, dist_nParticle, " ", "dist_partLine")
    dist_partLine(dist_nParticle) = trim(inLine)

  case("ECHO")
    if(nSplit >= 2) then
      dist_echoFile = trim(lnSplit(2))
    else
      dist_echoFile = "echo_distribution.dat"
    end if
    dist_echo = .true.

  case default
    write(lerr,"(a)") "DIST> ERROR Unknown keyword '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  case("EMIT","EMITTANCE")
    if(nSplit /= 3 .and. nSplit /= 4) then
      write(lerr,"(a,i0)") "DIST> ERROR EMITTANCE takes 2 or 3 arguments, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       EMITTANCE emit1 emit2 [emit3]"
      iErr = .true.
      return
    end if
    if(nSplit > 1) call chr_cast(lnSplit(2), dist_beamEmit(1), cErr)
    if(nSplit > 2) call chr_cast(lnSplit(3), dist_beamEmit(2), cErr)
    if(nSplit > 3) call chr_cast(lnSplit(4), dist_beamEmit(3), cErr)

  end select

  if(cErr) then
    iErr = .true.
    return
  end if

end subroutine dist_parseInputLine

! ================================================================================================ !
!  Parse File Column Formats
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-05
!  Updated: 2019-07-10
!
!    This routine splits the unit from the format description and adds the format and the scaling
!  factor to the list of format columns for later file parsing.
! ================================================================================================ !
subroutine dist_setColumnFormat(fmtName, fErr)

  use crcoall
  use string_tools
  use numerical_constants

  character(len=*), intent(in)  :: fmtName
  logical,          intent(out) :: fErr

  character(len=20) fmtBase
  character(len=10) fmtUnit
  real(kind=fPrec)  uFac
  integer           c, unitPos, fmtLen

  fErr = .false.

  fmtLen = len_trim(fmtName)
  if(fmtLen < 1) then
    write(lerr,"(a)") "DIST> ERROR Unknown column format '"//trim(fmtName)//"'"
    fErr = .true.
  end if

  unitPos = -1
  do c=1,fmtLen
    if(fmtName(c:c) == "[") then
      unitPos = c
      exit
    end if
  end do

  if(unitPos > 1 .and. fmtLen > unitPos+1) then
    fmtBase = trim(fmtName(1:unitPos-1))
    fmtUnit = trim(fmtName(unitPos:fmtLen))
  else
    fmtBase = trim(fmtName)
    fmtUnit = "[1]"
  end if

  select case(chr_toUpper(fmtBase))

  case("OFF","SKIP")
    call dist_appendFormat(dist_fmtNONE,     one,  0)
  case("ID")
    call dist_appendFormat(dist_fmtPartID,   one,  0)
  case("PARENT")
    call dist_appendFormat(dist_fmtParentID, one,  0)

  case("X")
    call dist_unitScale(fmtName, fmtUnit, 1, uFac, fErr)
    call dist_appendFormat(dist_fmtX,        uFac, 1)
  case("Y")
    call dist_unitScale(fmtName, fmtUnit, 1, uFac, fErr)
    call dist_appendFormat(dist_fmtY,        uFac, 3)

  case("XP")
    call dist_appendFormat(dist_fmtXP,       uFac, 2)
  case("YP")
    call dist_appendFormat(dist_fmtYP,       uFac, 4)

  case("PX")
    call dist_unitScale(fmtName, fmtUnit, 3, uFac, fErr)
    call dist_appendFormat(dist_fmtPX,       uFac, 2)
  case("PY")
    call dist_unitScale(fmtName, fmtUnit, 3, uFac, fErr)
    call dist_appendFormat(dist_fmtPY,       uFac, 4)

  case("PX/P0","PXP0")
    call dist_appendFormat(dist_fmtPX,       one,  2)
  case("PY/P0","PYP0")
    call dist_appendFormat(dist_fmtPY,       one,  4)

  case("SIGMA")
    call dist_unitScale(fmtName, fmtUnit, 1, uFac, fErr)
    call dist_appendFormat(dist_fmtSIGMA,    uFac, 5)
  case("ZETA")
    call dist_unitScale(fmtName, fmtUnit, 1, uFac, fErr)
    call dist_appendFormat(dist_fmtZETA,     uFac, 5)
  case("DT")
    call dist_unitScale(fmtName, fmtUnit, 4, uFac, fErr)
    call dist_appendFormat(dist_fmtDT,       uFac, 5)

  case("E")
    call dist_unitScale(fmtName, fmtUnit, 3, uFac, fErr)
    call dist_appendFormat(dist_fmtE,        uFac, 6)
  case("P")
    call dist_unitScale(fmtName, fmtUnit, 3, uFac, fErr)
    call dist_appendFormat(dist_fmtP,        uFac, 6)
  case("DE/E0","DEE0")
    call dist_appendFormat(dist_fmtDEE0,     one,  6)
  case("DP/P0","DPP0","DELTA")
    call dist_appendFormat(dist_fmtDPP0,     one,  6)
  case("DE/P0","DEP0","PT")
    call dist_appendFormat(dist_fmtPT,       one,  6)

  case("XN")
    call dist_appendFormat(dist_fmtXN,       one,  1)
    dist_distLib = .true.
  case("YN")
    call dist_appendFormat(dist_fmtYN,       one,  3)
    dist_distLib = .true.
  case("ZN")
    call dist_appendFormat(dist_fmtZN,       one,  5)
    dist_distLib = .true.
  case("PXN")
    call dist_appendFormat(dist_fmtPXN,      one,  2)
    dist_distLib = .true.
  case("PYN")
    call dist_appendFormat(dist_fmtPYN,      one,  4)
    dist_distLib = .true.
  case("PZN")
    call dist_appendFormat(dist_fmtPZN,      one,  6)
    dist_distLib = .true.

  case("JX")
    call dist_appendFormat(dist_fmtJX,       one,  1)
    dist_distLib = .true.
  case("JY")
    call dist_appendFormat(dist_fmtJY,       one,  3)
    dist_distLib = .true.
  case("JZ")
    call dist_appendFormat(dist_fmtJZ,       one,  5)
    dist_distLib = .true.
  case("PHIX")
    call dist_unitScale(fmtName, fmtUnit, 2, uFac, fErr)
    call dist_appendFormat(dist_fmtPhiX,     uFac, 2)
    dist_distLib = .true.
  case("PHIY")
    call dist_unitScale(fmtName, fmtUnit, 2, uFac, fErr)
    call dist_appendFormat(dist_fmtPhiY,     uFac, 4)
    dist_distLib = .true.
  case("PHIZ")
    call dist_unitScale(fmtName, fmtUnit, 2, uFac, fErr)
    call dist_appendFormat(dist_fmtPhiZ,     uFac, 6)
    dist_distLib = .true.

  case("MASS","M")
    call dist_unitScale(fmtName, fmtUnit, 3, uFac, fErr)
    call dist_appendFormat(dist_fmtMASS,     uFac, 0)
  case("CHARGE","Q")
    call dist_appendFormat(dist_fmtCHARGE,   one,  0)
  case("ION_A")
    call dist_appendFormat(dist_fmtIonA,     one,  0)
  case("ION_Z")
    call dist_appendFormat(dist_fmtIonZ,     one,  0)
  case("PDGID")
    call dist_appendFormat(dist_fmtPDGID,    one,  0)

  case("4D") ! 4D default coordinates
    call dist_appendFormat(dist_fmtX,        one,  1)
    call dist_appendFormat(dist_fmtPX,       one,  2)
    call dist_appendFormat(dist_fmtY,        one,  3)
    call dist_appendFormat(dist_fmtPY,       one,  4)

  case("6D") ! 6D default coordinates
    call dist_appendFormat(dist_fmtX,        one,  1)
    call dist_appendFormat(dist_fmtPX,       one,  2)
    call dist_appendFormat(dist_fmtY,        one,  3)
    call dist_appendFormat(dist_fmtPY,       one,  4)
    call dist_appendFormat(dist_fmtZETA,     one,  5)
    call dist_appendFormat(dist_fmtDPP0,     one,  6)

  case("NORM") ! 6D normalised coordinates
    call dist_appendFormat(dist_fmtXN,       one,  1)
    call dist_appendFormat(dist_fmtPXN,      one,  2)
    call dist_appendFormat(dist_fmtYN,       one,  3)
    call dist_appendFormat(dist_fmtPYN,      one,  4)
    call dist_appendFormat(dist_fmtZN,       one,  5)
    call dist_appendFormat(dist_fmtPZN,      one,  6)
    dist_distLib = .true.

  case("ACTION") ! 6D action
    call dist_appendFormat(dist_fmtJX,       one,  1)
    call dist_appendFormat(dist_fmtPhiX,     one,  2)
    call dist_appendFormat(dist_fmtJY,       one,  3)
    call dist_appendFormat(dist_fmtPhiY,     one,  4)
    call dist_appendFormat(dist_fmtJZ,       one,  5)
    call dist_appendFormat(dist_fmtPhiZ,     one,  6)
    dist_distLib = .true.

  case("IONS") ! The ion columns
    call dist_appendFormat(dist_fmtMASS,     c1e3, 0)
    call dist_appendFormat(dist_fmtCHARGE,   one,  0)
    call dist_appendFormat(dist_fmtIonA,     one,  0)
    call dist_appendFormat(dist_fmtIonZ,     one,  0)
    call dist_appendFormat(dist_fmtPDGID,    one,  0)

  case("OLD_DIST") ! The old DIST block file format
    ! This is added for the block to be compatible with the old, fixed column DIST file format.
    ! The scaling is hardcoded, and the file format has a few columns that are not used.
    ! Specifically the weight (column 3), and the z and pz coordinates (columns 6 and 9).
    call dist_appendFormat(dist_fmtPartID,   one,  0)
    call dist_appendFormat(dist_fmtParentID, one,  0)
    call dist_appendFormat(dist_fmtNONE,     one,  0)
    call dist_appendFormat(dist_fmtX,        c1e3, 1)
    call dist_appendFormat(dist_fmtY,        c1e3, 3)
    call dist_appendFormat(dist_fmtNONE,     one,  0)
    call dist_appendFormat(dist_fmtXP,       c1e3, 2)
    call dist_appendFormat(dist_fmtYP,       c1e3, 4)
    call dist_appendFormat(dist_fmtNONE,     one,  0)
    call dist_appendFormat(dist_fmtIonA,     one,  0)
    call dist_appendFormat(dist_fmtIonZ,     one,  0)
    call dist_appendFormat(dist_fmtMASS,     c1e3, 0)
    call dist_appendFormat(dist_fmtP,        c1e3, 6)
    call dist_appendFormat(dist_fmtDT,       one,  5)

  case default
    write(lerr,"(a)") "DIST> ERROR Unknown column format '"//trim(fmtName)//"'"
    fErr = .true.
    return

  end select

#ifndef DISTLIB
  if(dist_distLib) then
    write(lerr,"(a)") "DIST> ERROR Format '"//trim(fmtName)//"' requires SixTrack to be built with DISTLIB enabled"
    fErr = .true.
    return
  end if
#endif

end subroutine dist_setColumnFormat

! ================================================================================================ !
!  Parse File Column Units
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-05
!  Updated: 2019-07-10
! ================================================================================================ !
subroutine dist_unitScale(fmtName, fmtUnit, unitType, unitScale, uErr)

  use crcoall
  use string_tools
  use numerical_constants

  character(len=*), intent(in)    :: fmtName
  character(len=*), intent(in)    :: fmtUnit
  integer,          intent(in)    :: unitType
  real(kind=fPrec), intent(out)   :: unitScale
  logical,          intent(inout) :: uErr

  unitScale = zero

  if(unitType == 1) then ! Positions
    select case(chr_toLower(fmtUnit))
    case("[m]")
      unitScale = c1e3
    case("[mm]")
      unitScale = one
    case("[1]")
      unitScale = one
    case default
      goto 100
    end select
  elseif(unitType == 2) then ! Angle
    select case(chr_toLower(fmtUnit))
    case("[rad]")
      unitScale = c1e3
    case("[mrad]")
      unitScale = one
    case("[1]")
      unitScale = one
    case default
      goto 100
    end select
  elseif(unitType == 3) then ! Energy/Momentum/Mass
    select case(chr_toLower(fmtUnit))
    case("[ev]", "[ev/c]", "[ev/c^2]", "[ev/c2]", "[ev/c**2]")
      unitScale = c1m6
    case("[kev]","[kev/c]","[kev/c^2]","[kev/c2]","[kev/c**2]")
      unitScale = c1m3
    case("[mev]","[mev/c]","[mev/c^2]","[mev/c2]","[mev/c**2]")
      unitScale = one
    case("[gev]","[gev/c]","[gev/c^2]","[gev/c2]","[gev/c**2]")
      unitScale = c1e3
    case("[tev]","[tev/c]","[tev/c^2]","[tev/c2]","[tev/c**2]")
      unitScale = c1e6
    case("[1]")
      unitScale = one
    case default
      goto 100
    end select
  elseif(unitType == 4) then ! Time
    select case(chr_toLower(fmtUnit))
    case("[s]")
      unitScale = one
    case("[ms]")
      unitScale = c1e3
    case("[us]","[µs]")
      unitScale = c1e6
    case("[ns]")
      unitScale = c1e9
    case("[ps]")
      unitScale = c1e12
    case("[1]")
      unitScale = one
    case default
      goto 100
    end select
  end if

  return

100 continue
  write(lerr,"(a)") "DIST> ERROR Unrecognised or invalid unit for format identifier '"//trim(fmtName)//"'"
  uErr = .true.

end subroutine dist_unitScale

! ================================================================================================ !
!  Save File Column Format
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-05
!  Updated: 2019-07-05
! ================================================================================================ !
subroutine dist_appendFormat(fmtID, colScale, partCol)

  use crcoall
  use mod_alloc
  use numerical_constants

  integer,          intent(in) :: fmtID
  real(kind=fPrec), intent(in) :: colScale
  integer,          intent(in) :: partCol

  dist_nColumns = dist_nColumns + 1

  call alloc(dist_colFormat, dist_nColumns, dist_fmtNONE, "dist_colFormat")
  call alloc(dist_colScale,  dist_nColumns, one,          "dist_colScale")

  dist_colFormat(dist_nColumns) = fmtID
  dist_colScale(dist_nColumns)  = colScale

  if(partCol > 0) then
    if(dist_partFmt(partCol) == 0) then
      dist_partFmt(partCol) = fmtID
    else
      write(lerr,"(a,i0)") "DIST> ERROR Multiple formats selected for particle coordinate ",partCol
      select case(partCol)
      case(1)
        write(lerr,"(a)") "DIST>      Choose one of: X, XN, JX"
      case(2)
        write(lerr,"(a)") "DIST>      Choose one of: XP, PX, PX/P0, PXN, PHIX"
      case(3)
        write(lerr,"(a)") "DIST>      Choose one of: Y, YN, JY"
      case(4)
        write(lerr,"(a)") "DIST>      Choose one of: YP, PY, PY/P0, PYN, PHIY"
      case(5)
        write(lerr,"(a)") "DIST>      Choose one of: SIGMA, ZETA, DT, ZN, JZ"
      case(6)
        write(lerr,"(a)") "DIST>      Choose one of: E, P, DE/E0, DP/P0, PT, PZN, PHIZ"
      end select
      call prror
    end if
  end if

end subroutine dist_appendFormat

! ================================================================================================ !
!  Parse PARTICLE Lines from DIST Block
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-08
!  Updated: 2019-07-09
! ================================================================================================ !
subroutine dist_parseParticles

  use mod_common
  use string_tools

  character(len=:), allocatable :: lnSplit(:)
  integer i, j, nSplit
  logical spErr, cErr

  if(dist_nParticle > 0) then
    do j=1,dist_nParticle
      call chr_split(dist_partLine(j), lnSplit, nSplit, spErr)
      if(spErr) goto 20
      do i=2,nSplit
        call dist_saveParticle(j, i-1, lnSplit(i), cErr)
      end do
      if(cErr) goto 20
      dist_numPart = dist_numPart + 1
    end do
  end if

  return

20 continue
  write(lout,"(a,i0,a)") "DIST> ERROR Could not parse PARTICLE definition number ",j," from "//trim(fort3)

end subroutine dist_parseParticles

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-07-09
! ================================================================================================ !
subroutine dist_readDist

  use parpro
  use mod_common
  use mod_common_main
  use string_tools
  use numerical_constants
  use mod_units

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn) inLine
  integer i, nSplit, lineNo, fUnit, nRead
  logical spErr, cErr

  write(lout,"(a)") "DIST> Reading particles from '"//trim(dist_distFile)//"'"

  nRead  = 0
  lineNo = 0
  cErr   = .false.

  if(dist_hasFormat .eqv. .false.) then
    call dist_setColumnFormat("OLD_DIST",cErr)
  end if

  call f_requestUnit(dist_distFile, fUnit)
  call f_open(unit=fUnit,file=dist_distFile,mode='r',err=cErr,formatted=.true.,status="old")
  if(cErr) goto 20

10 continue
  read(fUnit,"(a)",end=40,err=30) inLine
  lineNo = lineNo + 1

  if(inLine(1:1) == "*") goto 10
  if(inLine(1:1) == "#") goto 10
  if(inLine(1:1) == "!") goto 10

  dist_numPart = dist_numPart + 1
  nRead = nRead + 1

  if(dist_numPart > napx) then
    write(lout,"(a,i0,a)") "DIST> Stopping reading file as ",napx," particles have been read, as requested in '"//trim(fort3)//"'"
    goto 40
  end if

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) goto 30
  if(nSplit == 0) goto 10

  if(nSplit /= dist_nColumns) then
    write(lerr,"(3(a,i0),a)") "DIST> ERROR Number of columns in file on line ",lineNo," is ",nSplit,&
      " but FORMAT defines ",dist_nColumns," columns"
    call prror
  end if
  do i=1,nSplit
    call dist_saveParticle(dist_numPart, i, lnSplit(i), cErr)
  end do
  if(cErr) goto 30

  goto 10

20 continue
  write(lerr,"(a)") "DIST> ERROR Opening file '"//trim(dist_distFile)//"'"
  call prror
  return

30 continue
  write(lerr,"(a,i0)") "DIST> ERROR Reading particles from line ",lineNo
  call prror
  return

40 continue
  call f_close(fUnit)
  write(lout,"(a,i0,a)") "DIST> Read ",nRead," particles from file '"//trim(dist_distFile)//"'"

end subroutine dist_readDist

! ================================================================================================ !
!  Save Particle Data to Arrays
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-08
!  Updated: 2019-07-09
! ================================================================================================ !
subroutine dist_saveParticle(partNo, colNo, inVal, sErr)

  use mod_common
  use mod_common_main
  use string_tools

  integer,          intent(in)    :: partNo
  integer,          intent(in)    :: colNo
  character(len=*), intent(in)    :: inVal
  logical,          intent(inout) :: sErr

  select case(dist_colFormat(colNo))

  case(dist_fmtNONE)
    return

  case(dist_fmtPartID)
    call chr_cast(inVal, partID(partNo), sErr)

  case(dist_fmtParentID)
    call chr_cast(inVal, parentID(partNo), sErr)

  !  Horizontal Coordinates
  ! ========================

  case(dist_fmtX, dist_fmtXN, dist_fmtJX)
    call chr_cast(inVal, dist_partCol1(partNo), sErr)
    dist_partCol1(partNo) = dist_partCol1(partNo) * dist_colScale(colNo)

  case(dist_fmtPX, dist_fmtXP, dist_fmtPXP0, dist_fmtPXN, dist_fmtPhiX)
    call chr_cast(inVal, dist_partCol2(partNo), sErr)
    dist_partCol2(partNo) = dist_partCol2(partNo) * dist_colScale(colNo)

  !  Vertical Coordinates
  ! ========================

  case(dist_fmtY, dist_fmtYN, dist_fmtJY)
    call chr_cast(inVal, dist_partCol3(partNo), sErr)
    dist_partCol3(partNo) = dist_partCol3(partNo) * dist_colScale(colNo)

  case(dist_fmtPY, dist_fmtYP, dist_fmtPYP0, dist_fmtPYN, dist_fmtPhiY)
    call chr_cast(inVal, dist_partCol4(partNo), sErr)
    dist_partCol4(partNo) = dist_partCol4(partNo) * dist_colScale(colNo)

  !  Longitudinal Coordinates
  ! ==========================

  case(dist_fmtSIGMA, dist_fmtZETA, dist_fmtDT, dist_fmtZN, dist_fmtJZ)
    call chr_cast(inVal, dist_partCol5(partNo), sErr)
    dist_partCol5(partNo) = dist_partCol5(partNo) * dist_colScale(colNo)

  case(dist_fmtE, dist_fmtDEE0, dist_fmtPT)
    call chr_cast(inVal, dist_partCol6(partNo), sErr)
    dist_partCol6(partNo) = dist_partCol6(partNo) * dist_colScale(colNo)

  case(dist_fmtP)
    call chr_cast(inVal, dist_partCol6(partNo), sErr)
    dist_partCol6(partNo) = dist_partCol6(partNo) * dist_colScale(colNo)

  case(dist_fmtDPP0, dist_fmtPZN, dist_fmtPhiZ)
    call chr_cast(inVal, dist_partCol6(partNo), sErr)
    dist_partCol6(partNo) = dist_partCol6(partNo) * dist_colScale(colNo)

  !  Ion Columns
  ! =============

  case(dist_fmtMASS)
    call chr_cast(inVal, nucm(partNo), sErr)
    nucm(partNo)  = nucm(partNo) * dist_colScale(colNo)
    dist_readMass = .true.

  case(dist_fmtCHARGE)
    call chr_cast(inVal, nqq(partNo), sErr)
    dist_readCharge = .true.

  case(dist_fmtIonA)
    call chr_cast(inVal, naa(partNo), sErr)
    dist_readIonA = .true.

  case(dist_fmtIonZ)
    call chr_cast(inVal, nzz(partNo), sErr)
    dist_readIonZ = .true.

  case(dist_fmtPDGID)
    call chr_cast(inVal, pdgid(partNo), sErr)
    dist_readPDGID = .true.

  end select

end subroutine dist_saveParticle

! ================================================================================================ !
!  Post-Processing of Particle Arrays
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-05
!  Updated: 2019-07-09
! ================================================================================================ !
subroutine dist_postprParticles

  use crcoall
  use mod_common
  use mod_particles
  use mod_common_main
  use physical_constants
  use numerical_constants

  logical doAction, doNormal
  real(kind=fPrec) beta0

  doAction = .false.
  doNormal = .false.

  beta0 = e0f/e0

  ! Forward energy/momentum must be calculated first
  select case(dist_partFmt(6))
  case(dist_fmtNONE)
    ejv(1:napx) = e0
    call part_updatePartEnergy(1,.false.)
  case(dist_fmtE)
    ejv(1:napx) = dist_partCol6(1:napx)
    call part_updatePartEnergy(1,.false.)
  case(dist_fmtP)
    ejfv(1:napx) = dist_partCol6(1:napx)
    call part_updatePartEnergy(2,.false.)
  case(dist_fmtDEE0)
    ejv(1:napx) = (one + dist_partCol6(1:napx))*e0
    call part_updatePartEnergy(1,.false.)
  case(dist_fmtDPP0)
    dpsv(1:napx) = dist_partCol6(1:napx)
    call part_updatePartEnergy(3,.false.)
  case(dist_fmtPT)
    ejv(1:napx) = (dist_partCol6(1:napx)*e0f)*beta0 + e0
    call part_updatePartEnergy(1,.false.)
  case(dist_fmtPZN)
    doNormal = .true.
  case(dist_fmtPhiZ)
    doAction = .true.
  end select

  select case(dist_partFmt(5))
  case(dist_fmtNONE)
    sigmv(1:napx) = zero
  case(dist_fmtZETA)
    sigmv(1:napx) = dist_partCol5(1:napx)/rvv(1:napx)
  case(dist_fmtSIGMA)
    sigmv(1:napx) = dist_partCol5(1:napx)
  case(dist_fmtDT)
    sigmv(1:napx) = -beta0*(dist_partCol5(1:napx)*clight)
  case(dist_fmtZN)
    doNormal = .true.
  case(dist_fmtJZ)
    doAction = .true.
  end select

  select case(dist_partFmt(1))
  case(dist_fmtNONE)
    xv1(1:napx) = zero
  case(dist_fmtX)
    xv1(1:napx) = dist_partCol1(1:napx)
  case(dist_fmtXN)
    doNormal = .true.
  case(dist_fmtJX)
    doAction = .true.
  end select

  select case(dist_partFmt(2))
  case(dist_fmtNONE)
    yv1(1:napx) = zero
  case(dist_fmtXP)
    yv1(1:napx) = dist_partCol2(1:napx)
  case(dist_fmtPX)
    yv1(1:napx) = (dist_partCol2(1:napx)/ejfv(1:napx))*c1e3
  case(dist_fmtPXP0)
    yv1(1:napx) = ((dist_partCol2(1:napx)*e0f)/ejfv(1:napx))*c1e3
  case(dist_fmtPXN)
    doNormal = .true.
  case(dist_fmtPhiX)
    doAction = .true.
  end select

  select case(dist_partFmt(3))
  case(dist_fmtNONE)
    xv2(1:napx) = zero
  case(dist_fmtY)
    xv2(1:napx) = dist_partCol3(1:napx)
  case(dist_fmtYN)
    doNormal = .true.
  case(dist_fmtJY)
    doAction = .true.
  end select

  select case(dist_partFmt(4))
  case(dist_fmtNONE)
    yv2(1:napx) = zero
  case(dist_fmtYP)
    yv2(1:napx) = dist_partCol4(1:napx)
  case(dist_fmtPY)
    yv2(1:napx) = (dist_partCol4(1:napx)/ejfv(1:napx))*c1e3
  case(dist_fmtPYP0)
    yv2(1:napx) = ((dist_partCol4(1:napx)*e0f)/ejfv(1:napx))*c1e3
  case(dist_fmtPYN)
    doNormal = .true.
  case(dist_fmtPhiY)
    doAction = .true.
  end select

#ifdef DISTLIB
  if(doNormal .and. doAction) then
    write(lerr,"(a)") "DIST> ERROR Cannot mix normalised and action coordinates"
    call prror
  end if

  if(doNormal) then
    call distlib_setNormalised(                    &
      dist_partCol1, dist_partCol2, dist_partCol3, &
      dist_partCol4, dist_partCol5, dist_partCol6, &
      napx                                         &
    )
  end if

  if(doAction) then
    ! Call DISTlib
    continue
  end if
#endif

end subroutine dist_postprParticles

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-07-09
! ================================================================================================ !
subroutine dist_finaliseDist

  use parpro
  use mod_pdgid
  use mod_common
  use mod_common_main
  use mod_common_track
  use numerical_constants

  real(kind=fPrec) chkP, chkE
  integer j

  ! Finish the ION bits
  if(dist_numPart /= napx) then
    write(lerr,"(2(a,i0),a)") "DIST> ERROR Number of particles read or generated is ",dist_numPart,&
      " but the simulation setup requests ",napx," particles"
    call prror
  end if

  if(dist_readIonA .neqv. dist_readIonZ) then
    write(lerr,"(a)") "DIST> ERROR ION_A and ION_Z columns have to both be present if one of them is"
    call prror
  end if

  if(dist_readIonZ .and. .not. dist_readCharge) then
    nqq(1:napx) = nzz(1:napx)
  end if
  mtc(1:napx) = (nqq(1:napx)*nucm0)/(qq0*nucm(1:napx))

  if(dist_readIonA .and. dist_readIonZ .and. .not. dist_readPDGID) then
    do j=1,napx
      call CalculatePDGid(pdgid(j), naa(j), nzz(j))
    end do
  end if

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
        if(nzz(j) == zz0 .or. naa(j) == aa0 .or. nqq(j) == qq0 .or. pdgid(j) == pdgid0) then
          naa(j)   = aa0
          nzz(j)   = zz0
          nqq(j)   = qq0
          pdgid(j) = pdgid0
          mtc(j)   = one
        else
          write(lerr,"(a)") "DIST> ERROR Mass and/or charge mismatch with relation to sync particle"
          call prror
        end if
      end if

      write(lout,"(a,4(1x,1pe25.18))") "DIST> CORRECTED:", ejfv(j), ejv(j), dpsv(j), oidpsv(j)
    end if
  end do

  write(lout,"(a,2(1x,i0),1x,f15.7,2(1x,i0))") "DIST> Reference particle species [A,Z,M,Q,ID]:", aa0, zz0, nucm0, qq0, pdgid0
  write(lout,"(a,1x,f15.7)")       "DIST> Reference energy [Z TeV]:", c1m6*e0/qq0

end subroutine dist_finaliseDist

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-07-08
! ================================================================================================ !
subroutine dist_echoDist

  use mod_common
  use mod_common_main
  use mod_units

  integer j, fUnit
  logical cErr

  call f_requestUnit(dist_echoFile, fUnit)
  call f_open(unit=fUnit,file=dist_echoFile,mode='w',err=cErr,formatted=.true.)
  if(cErr) goto 19

  rewind(fUnit)
  write(fUnit,"(a,1pe25.18)") "# Total energy of synch part [MeV]: ",e0
  write(fUnit,"(a,1pe25.18)") "# Momentum of synch part [MeV/c]:   ",e0f
  write(fUnit,"(a)")          "#"
  write(fUnit,"(a)")          "# x[mm], y[mm], xp[mrad], yp[mrad], sigmv[mm], ejfv[MeV/c]"
  do j=1, napx
    write(fUnit,"(6(1x,1pe25.18))") xv1(j), yv1(j), xv2(j), yv2(j), sigmv(j), ejfv(j)
  end do
  call f_close(fUnit)

  return

19 continue
  write(lerr,"(a)") "DIST> ERROR Opening file '"//trim(dist_echoFile)//"'"
  call prror

end subroutine dist_echoDist

end module mod_dist
