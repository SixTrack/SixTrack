! ================================================================================================ !
!
!  Beam Distribution Block
! ~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-07-12
!
!  This module was completely rewritten to support DISTlib in July 2019
!
! ================================================================================================ !
module mod_dist

  use crcoall
  use floatPrecision
  use parpro, only : mFileName
  use numerical_constants, only : zero, one

  implicit none

  logical,                  public,  save :: dist_enable      = .false. ! DIST input block given
  logical,                  private, save :: dist_echo        = .false. ! Echo the read distribution?
  logical,                  private, save :: dist_hasFormat   = .false. ! Whether the FORMAT keyword exists in block or not
  logical,                  private, save :: dist_libRead     = .false. ! Read file with dist library instead of internal reader
  logical,                  private, save :: dist_distLib     = .false. ! DISTlib is needed to generate the distribution requested
  logical,                  private, save :: dist_distLibNorm = .true.  ! DISTlib is used to convert normalised coordinates
  logical,                  private, save :: dist_needsRnd    = .false. ! If the DIST block needs random numbers
  character(len=mFileName), private, save :: dist_distFile    = " "     ! File name for reading the distribution
  character(len=mFileName), private, save :: dist_echoFile    = " "     ! File name for echoing the distribution

  integer,          allocatable, private, save :: dist_colFormat(:)   ! The column types in the FORMAT
  real(kind=fPrec), allocatable, private, save :: dist_colScale(:)    ! Scaling factor for the columns
  integer,                       private, save :: dist_nColumns  = 0  ! The number of columns in the FORMAT

  character(len=:), allocatable, private, save :: dist_partLine(:)    ! PARTICLE definitions in the block
  integer,                       private, save :: dist_nParticle = 0  ! Number of PARTICLE keywords in block

  integer,                       private, save :: dist_numPart   = 0  ! Number of actual particles generated or read

  type, private :: dist_fillType
    character(len=:), allocatable :: fillName
    integer                       :: fillTarget = 0
    integer                       :: fillMethod = 0
    integer,          allocatable :: iParams(:)
    real(kind=fPrec), allocatable :: fParams(:)
  end type dist_fillType

  type(dist_fillType), allocatable, private, save :: dist_fillList(:)
  integer,                          private, save :: dist_nFill = 0

  !  Fill Methods
  ! ===============
  integer, parameter :: dist_fillNONE       = 0  ! No fill
  integer, parameter :: dist_fillINT        = 1  ! Fixed integer value
  integer, parameter :: dist_fillFLOAT      = 2  ! Fixed float value
  integer, parameter :: dist_fillGAUSS      = 3  ! Gaussian distribution
  integer, parameter :: dist_fillRAYLEIGH   = 4  ! Rayleigh distribution
  integer, parameter :: dist_fillUNIFORM    = 5  ! Uniform distribution
  integer, parameter :: dist_fillLINEAR     = 6  ! Linear fill
  integer, parameter :: dist_fillCOUNT      = 7  ! Integer range

  !  Normalisation Methods
  ! =======================
  integer, parameter :: dist_normNONE       = 0  ! No normalisation
  integer, parameter :: dist_normSIXTRACK   = 1  ! Use SixTrack tas matrix
  integer, parameter :: dist_normINPUT      = 2  ! Use matrix in DIST block
  integer, parameter :: dist_normTWISS      = 3  ! Use twiss in DIST block

  integer,          private, save :: dist_normMethod  = dist_normSIXTRACK ! Flag for source of normalisation matrix
  integer,          private, save :: dist_tMatRow     = 0                 ! How many T-matrix rows we've read
  real(kind=fPrec), private, save :: dist_tMat(6,6)   = zero              ! T-matrix to use for normalisation
  real(kind=fPrec), private, save :: dist_twBeta(2)   = one               ! Twiss beta
  real(kind=fPrec), private, save :: dist_twAlpha(2)  = zero              ! Twiss alpha
  real(kind=fPrec), private, save :: dist_bDisp(4)    = zero              ! Dispersion
  real(kind=fPrec), private, save :: dist_beamEmit(3) = zero              ! Beam emittance
  integer,          private, save :: dist_emit3Unit   = 0                 ! Flag for converting longitudinal emittance

  !  Column Formats
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
  integer, parameter :: dist_fmtPSIGMA      = 27 ! Delta energy over reference momentum, but without beta0

  ! Normalised Coordinates
  integer, parameter :: dist_fmtXN          = 41 ! Normalised horizontal position
  integer, parameter :: dist_fmtYN          = 42 ! Normalised vertical position
  integer, parameter :: dist_fmtZN          = 43 ! Normalised longitudinal position
  integer, parameter :: dist_fmtPXN         = 44 ! Normalised horizontal momentum
  integer, parameter :: dist_fmtPYN         = 45 ! Normalised vertical momentum
  integer, parameter :: dist_fmtPZN         = 46 ! Normalised longitudinal momentum

  integer, parameter :: dist_fmtJX          = 51 ! Horizontal action
  integer, parameter :: dist_fmtJY          = 52 ! Vertical action
  integer, parameter :: dist_fmtJZ          = 53 ! Longitudinal action
  integer, parameter :: dist_fmtPhiX        = 54 ! Horizontal action angle
  integer, parameter :: dist_fmtPhiY        = 55 ! Vertical action angle
  integer, parameter :: dist_fmtPhiZ        = 56 ! Longitudinal action angle

  ! Ion Columns
  integer, parameter :: dist_fmtMASS        = 61 ! Particle mass
  integer, parameter :: dist_fmtCHARGE      = 62 ! Particle Charge
  integer, parameter :: dist_fmtIonA        = 63 ! Ion atomic mass
  integer, parameter :: dist_fmtIonZ        = 64 ! Ion atomic number
  integer, parameter :: dist_fmtPDGID       = 65 ! Particle PDG ID

  ! Spin
  integer, parameter :: dist_fmtSPINX       = 71 ! Spin vector x component
  integer, parameter :: dist_fmtSPINY       = 72 ! Spin vector y component
  integer, parameter :: dist_fmtSPINZ       = 73 ! Spin vector z component

  ! Flags for columns we've set that we need to track for later checks
  logical, private, save :: dist_readMass     = .false.
  logical, private, save :: dist_readIonZ     = .false.
  logical, private, save :: dist_readIonA     = .false.
  logical, private, save :: dist_readCharge   = .false.
  logical, private, save :: dist_readPDGID    = .false.
  logical, private, save :: dist_readPartID   = .false.
  logical, private, save :: dist_readParentID = .false.

  ! Temporary particle arrays
  real(kind=fPrec), allocatable, private, save :: dist_partCol1(:) ! Ends up in array xv1
  real(kind=fPrec), allocatable, private, save :: dist_partCol2(:) ! Ends up in array yv1
  real(kind=fPrec), allocatable, private, save :: dist_partCol3(:) ! Ends up in array xv2
  real(kind=fPrec), allocatable, private, save :: dist_partCol4(:) ! Ends up in array yv2
  real(kind=fPrec), allocatable, private, save :: dist_partCol5(:) ! Ends up in array sigmv
  real(kind=fPrec), allocatable, private, save :: dist_partCol6(:) ! Ends up in array dpsv, evj, and ejfv
  integer,                       private, save :: dist_partFmt(6) = dist_fmtNONE ! The format used for each column

#ifdef DISTLIB

  ! Coordiante Types
  ! Keep in sync with lib/DISTlib/distinterface.c
  integer, parameter :: dist_coordTypeAction   = 0
  integer, parameter :: dist_coordTypeNormal   = 1
  integer, parameter :: dist_coordTypePhysical = 2
  integer, parameter :: dist_coordTypeMixed    = 3

  ! C interface
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

    subroutine distlib_setTwiss(betaX, alphaX, betaY, alphaY) bind(C, name="settwisstas")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), value, intent(in) :: betaX
      real(kind=C_DOUBLE), value, intent(in) :: alphaX
      real(kind=C_DOUBLE), value, intent(in) :: betaY
      real(kind=C_DOUBLE), value, intent(in) :: alphaY
    end subroutine distlib_setTwiss

    subroutine distlib_setDispersion(dX, dPx, dY, dPy) bind(C, name="setdisptas")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), value, intent(in) :: dX
      real(kind=C_DOUBLE), value, intent(in) :: dPx
      real(kind=C_DOUBLE), value, intent(in) :: dY
      real(kind=C_DOUBLE), value, intent(in) :: dPy
    end subroutine distlib_setDispersion

    subroutine distlib_setTasMatrix(flatTas) bind(C, name="settasmatrix")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), intent(in) :: flatTas(36)
    end subroutine distlib_setTasMatrix

    subroutine distlib_setCoords(arr1, arr2, arr3, arr4, arr5, arr6, arrLen, cType) bind(C, name="setcoords")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE),        intent(in) :: arr1(*)
      real(kind=C_DOUBLE),        intent(in) :: arr2(*)
      real(kind=C_DOUBLE),        intent(in) :: arr3(*)
      real(kind=C_DOUBLE),        intent(in) :: arr4(*)
      real(kind=C_DOUBLE),        intent(in) :: arr5(*)
      real(kind=C_DOUBLE),        intent(in) :: arr6(*)
      integer(kind=C_INT), value, intent(in) :: arrLen
      integer(kind=C_INT), value, intent(in) :: cType
    end subroutine distlib_setCoords

    subroutine distlib_getCoords(arr1, arr2, arr3, arr4, arr5, arr6, arrLen) bind(C, name="get6trackcoord")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE), intent(inout) :: arr1(*)
      real(kind=C_DOUBLE), intent(inout) :: arr2(*)
      real(kind=C_DOUBLE), intent(inout) :: arr3(*)
      real(kind=C_DOUBLE), intent(inout) :: arr4(*)
      real(kind=C_DOUBLE), intent(inout) :: arr5(*)
      real(kind=C_DOUBLE), intent(inout) :: arr6(*)
      integer(kind=C_INT), intent(inout) :: arrLen
    end subroutine distlib_getCoords

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
  use mod_random
  use mod_particles
  use mod_common_main
  use physical_constants
  use numerical_constants

  integer i, distLibPart
  real(kind=fPrec) flatTas(36)

  call alloc(dist_partCol1, napx, zero, "dist_partCol1")
  call alloc(dist_partCol2, napx, zero, "dist_partCol2")
  call alloc(dist_partCol3, napx, zero, "dist_partCol3")
  call alloc(dist_partCol4, napx, zero, "dist_partCol4")
  call alloc(dist_partCol5, napx, zero, "dist_partCol5")
  call alloc(dist_partCol6, napx, zero, "dist_partCol6")

  if(dist_needsRnd .and. .not. rnd_initOK) then
    write(lerr,"(a,i0,a)") "DIST> ERROR DIST module requires random numbers, but random numbers generator not initialised"
    call prror
  end if

  ! If the T matrix hasn't been set from input, we set it here
  if(dist_normMethod == dist_normSIXTRACK) then
    dist_tMat(1:6,1:6) = tas(1:6,1:6)
    dist_tMat(1:5,6)   = dist_tMat(1:5,6)*c1m3
    dist_tMat(6,1:5)   = dist_tMat(6,1:5)*c1e3
  end if
  if(dist_normMethod == dist_normINPUT .and. dist_tMatRow /= 6) then
    write(lerr,"(a,i0,a)") "DIST> ERROR The T-matrix must be 6 rows, but only ",dist_tMatRow," rows were given"
    call prror
  end if

  ! Scale longitudinal emittance
  select case(dist_emit3Unit)
  case(1) ! Micrometre
    dist_beamEmit(3) = dist_beamEmit(3)*c1m6
  case(2) ! eV s
    dist_beamEmit(3) = dist_beamEmit(3)*((clight*beta0)/((four*pi)*(e0*c1e6)))
  end select
  write(lout,"(a,1pe13.6,a)") "DIST> EMITTANCE X = ",dist_beamEmit(1)," m"
  write(lout,"(a,1pe13.6,a)") "DIST> EMITTANCE Y = ",dist_beamEmit(2)," m"
  write(lout,"(a,1pe13.6,a)") "DIST> EMITTANCE Z = ",dist_beamEmit(3)," m"

#ifdef DISTLIB
  if(dist_distLib) then
    call distlib_init(1)
    call distlib_setEnergyMass(e0, nucm0)
    call distlib_setEmittance12(dist_beamEmit(1)/gamma0, dist_beamEmit(2)/gamma0) ! DISTlib expects geometric emittance
    call distlib_setEmittance3(dist_beamEmit(3))
    select case(dist_normMethod)
    case(dist_normSIXTRACK, dist_normINPUT)
      do i=1,6
        write(lout,"(a,i0,a,6(1x,1pe13.6))") "DIST> TMATRIX(",i,",1:6) =",dist_tMat(i,1:6)
      end do
      flatTas = pack(transpose(dist_tMat),.true.)
      call distlib_setTasMatrix(flatTas)
    case(dist_normTWISS)
      call distlib_setTwiss(dist_twBeta(1), dist_twAlpha(1), dist_twBeta(2), dist_twAlpha(2))
      call distlib_setDispersion(dist_bDisp(1), dist_bDisp(2), dist_bDisp(3), dist_bDisp(4))
      write(lout,"(a,2(1x,1pe13.6))") "DIST> BETA  X/Y =",dist_twBeta
      write(lout,"(a,2(1x,1pe13.6))") "DIST> ALPHA X/Y =",dist_twAlpha
    end select
  end if
#endif

  if(dist_nParticle > 0) then
    ! Particles are read from DIST block in fort.3
    call dist_parseParticles
  elseif(dist_distFile /= " ") then
    if(dist_libRead) then
      ! Particles are read entirely in DISTlib
#ifdef DISTLIB
      call distlib_readFile(trim(dist_distFile), len_trim(dist_distFile))
#endif
    else
      call dist_readDist
    end if
  end if
  call dist_doFills         ! Do fills, and if requested, overwrite what is read from file
  call dist_postprParticles ! Copy 6D temp arrays to their correct ones, or send to DISTlib

#ifdef DISTLIB
  if(dist_distLib) then
    ! Our final coordinates are taken from DISTlib
    ! The 6D temp arrays were sent to DISTlib in dist_postprParticles
    distLibPart = napx
    call distlib_getCoords(                        &
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
  use mod_ranecu
  use string_tools
  use numerical_constants

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  real(kind=fPrec) fmtFac
  integer nSplit, i, fmtID, fmtCol, tmpOne, tmpTwo
  logical spErr, cErr, isValid

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
      call dist_parseColumn(lnSplit(i),cErr,fmtID,fmtFac,fmtCol,isValid)
      if(isValid) then
        call dist_appendFormat(fmtID,fmtFac,fmtCol)
        if(cErr) then
          iErr = .true.
          return
        end if
      else
        call dist_setMultiColFormat(lnSplit(i), isValid)
        if(isValid .eqv. .false.) then
          write(lerr,"(a)") "DIST> ERROR Unknown column format '"//trim(lnSplit(i))//"'"
          iErr = .true.
          return
        end if
      end if
    end do
    dist_hasFormat = .true.

  case("READ")
    if(nSplit /= 2 .and. nSplit /= 3) then
      write(lerr,"(a,i0)") "DIST> ERROR READ takes 1 or 2 arguments, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       READ filename [use_distlib]"
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

  case("EMIT","EMITTANCE")
    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "DIST> ERROR EMITTANCE takes 2 arguments, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       EMITTANCE emit1[mm mrad] emit2[mm mrad]"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dist_beamEmit(1), cErr)
    call chr_cast(lnSplit(3), dist_beamEmit(2), cErr)
    dist_beamEmit(1) = dist_beamEmit(1) * c1m6
    dist_beamEmit(2) = dist_beamEmit(2) * c1m6

  case("LONGEMIT")
    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "DIST> ERROR LONGEMIT takes 2 arguments, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       LONGEMIT emit3 unit[eVs|um]"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dist_beamEmit(3), cErr)
    select case(chr_toLower(lnSplit(3)))
    case("um","µm")
      dist_emit3Unit = 1
    case("evs")
      dist_emit3Unit = 2
    case default
      write(lerr,"(a)") "DIST> ERROR Unknown or unsupported unit '"//trim(lnSplit(3))//"' for LONGEMIT"
      iErr = .true.
      return
    end select

  case("TWISS")
    if(nSplit /= 5) then
      write(lerr,"(a,i0)") "DIST> ERROR TWISS takes 4 arguments, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       TWISS betaX[m] alphaX[1] betaY[m] alphaY[1]"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dist_twBeta(1),  cErr)
    call chr_cast(lnSplit(3), dist_twAlpha(1), cErr)
    call chr_cast(lnSplit(4), dist_twBeta(2),  cErr)
    call chr_cast(lnSplit(5), dist_twAlpha(2), cErr)
    dist_normMethod = dist_normTWISS

  case("DISPERSION")
    if(nSplit /= 5) then
      write(lerr,"(a,i0)") "DIST> ERROR DISPERSION takes 4 arguments, got ",nSplit-1
      write(lerr,"(a)")    "DIST>       DISPERSION dx dpx dy dpy"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dist_bDisp(1), cErr)
    call chr_cast(lnSplit(3), dist_bDisp(2), cErr)
    call chr_cast(lnSplit(4), dist_bDisp(3), cErr)
    call chr_cast(lnSplit(5), dist_bDisp(4), cErr)

  case("FILL")
    if(nSplit < 3) then
      write(lerr,"(a,i0)") "DIST> ERROR FILL takes at least 3 argument, got ",nSplit-1
      iErr = .true.
      return
    end if
    call dist_parseFill(lnSplit, nSplit, cErr)

  case("SEED")
    write(lerr,"(a,i0)") "DIST> ERROR SEED keyword has been removed and replaced by the RAND block"
    iErr = .true.
    return

  case("TMATRIX")
    if(nSplit /= 7) then
      write(lerr,"(a,i0)") "DIST> ERROR TMATRIX takes 6 arguments, got ",nSplit-1
      iErr = .true.
      return
    end if
    dist_tMatRow = dist_tMatRow + 1
    if(dist_tMatRow > 6) then
      write(lerr,"(a)") "DIST> ERROR Only 6 rows for the TMATRIX can be given"
      iErr = .true.
      return
    end if
    do i=1,6
      call chr_cast(lnSplit(i+1), dist_tMat(dist_tMatRow,i), cErr)
    end do
    dist_normMethod = dist_normINPUT

  case default
    write(lerr,"(a)") "DIST> ERROR Unknown keyword '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

  if(cErr) then
    iErr = .true.
    return
  end if

end subroutine dist_parseInputLine

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-30
!  Updated: 2019-07-30
! ================================================================================================ !
subroutine dist_parseFill(lnSplit, nSplit, iErr)

  use crcoall
  use mod_alloc
  use string_tools
  use numerical_constants

  character(len=:), allocatable, intent(in)    :: lnSplit(:)
  integer,                       intent(in)    :: nSplit
  logical,                       intent(inout) :: iErr

  type(dist_fillType), allocatable :: tmpFill(:)
  real(kind=fPrec),    allocatable :: fParams(:)
  integer,             allocatable :: iParams(:)

  real(kind=fPrec) fmtFac
  integer fmtCol, fillMethod, fillTarget, firstPart, lastPart
  logical fErr, cErr, isValid

  fErr = .false.
  cErr = .false.

  if(allocated(dist_fillList)) then
    allocate(tmpFill(dist_nFill+1))
    tmpFill(1:dist_nFill) = dist_fillList(1:dist_nFill)
    call move_alloc(tmpFill, dist_fillList)
    dist_nFill = dist_nFill + 1
  else
    allocate(dist_fillList(1))
    dist_nFill = 1
  end if

  call dist_parseColumn(lnSplit(2), fErr, fillTarget, fmtFac, fmtCol, isValid)
  if(fmtCol > 0) then
    if(dist_partFmt(fmtCol) == 0) then
      dist_partFmt(fmtCol) = fillTarget
    end if
    if(dist_partFmt(fmtCol) /= fillTarget) then
      call dist_multFormatError(fmtCol)
      iErr = .true.
      return
    end if
  end if

  select case(chr_toUpper(lnSplit(3)))

  case("NONE")

    allocate(fParams(1))
    allocate(iParams(2))

    fillMethod = dist_fillNONE
    fParams(1) =  zero
    iParams(1) =  1
    iParams(2) = -1

  case("INT")
    if(nSplit /= 4 .and. nSplit /= 6) then
      write(lerr,"(a,i0)") "DIST> ERROR FILL format INT expected 1 or 3 arguments, got ",nSplit-3
      write(lerr,"(a)")    "DIST>       FILL column INT value [first last]"
      iErr = .true.
      return
    end if

    allocate(fParams(1))
    allocate(iParams(3))

    fillMethod = dist_fillINT
    fParams(1) =  zero
    iParams(1) =  1
    iParams(2) = -1
    iParams(3) =  0

    if(nSplit > 3) call chr_cast(lnSplit(4), iParams(3), cErr)
    if(nSplit > 4) call chr_cast(lnSplit(5), iParams(1), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), iParams(2), cErr)

  case("FLOAT")
    if(nSplit /= 4 .and. nSplit /= 6) then
      write(lerr,"(a,i0)") "DIST> ERROR FILL format FLOAT expected 1 or 3 arguments, got ",nSplit-3
      write(lerr,"(a)")    "DIST>       FILL column FLOAT value [first last]"
      iErr = .true.
      return
    end if

    allocate(fParams(1))
    allocate(iParams(2))

    fillMethod = dist_fillFLOAT
    fParams(1) = zero
    iParams(1) =  1
    iParams(2) = -1

    if(nSplit > 3) call chr_cast(lnSplit(4), fParams(1), cErr)
    if(nSplit > 4) call chr_cast(lnSplit(5), iParams(1), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), iParams(2), cErr)

  case("GAUSS")
    if(nSplit /= 5 .and. nSplit /= 6 .and. nSplit /= 8) then
      write(lerr,"(a,i0)") "DIST> ERROR FILL format GAUSS expected 2, 3 or 5 arguments, got ",nSplit-3
      write(lerr,"(a)")    "DIST>       FILL column GAUSS sigma mu [cut] [first last]"
      iErr = .true.
      return
    end if

    allocate(fParams(3))
    allocate(iParams(2))

    fillMethod = dist_fillGAUSS
    fParams(1) = one
    fParams(2) = zero
    fParams(3) = zero
    iParams(1) =  1
    iParams(2) = -1

    if(nSplit > 3) call chr_cast(lnSplit(4), fParams(1), cErr)
    if(nSplit > 4) call chr_cast(lnSplit(5), fParams(2), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), fParams(3), cErr)
    if(nSplit > 6) call chr_cast(lnSplit(7), iParams(1), cErr)
    if(nSplit > 7) call chr_cast(lnSplit(8), iParams(2), cErr)

    dist_needsRnd = .true.

  case("RAYLEIGH")
    if(nSplit /= 4 .and. nSplit /= 5 .and. nSplit /= 6 .and. nSplit /= 8) then
      write(lerr,"(a,i0)") "DIST> ERROR FILL format RAYLEIGH expected 1, 2, 3 or 5 arguments, got ",nSplit-3
      write(lerr,"(a)")    "DIST>       FILL column RAYLEIGH sigma [maxcut] [mincut] [first last]"
      iErr = .true.
      return
    end if

    allocate(fParams(3))
    allocate(iParams(2))

    fillMethod = dist_fillRAYLEIGH
    fParams(1) = one
    fParams(2) = zero
    fParams(3) = zero
    iParams(1) =  1
    iParams(2) = -1

    if(nSplit > 3) call chr_cast(lnSplit(4), fParams(1), cErr)
    if(nSplit > 4) call chr_cast(lnSplit(5), fParams(2), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), fParams(3), cErr)
    if(nSplit > 6) call chr_cast(lnSplit(7), iParams(1), cErr)
    if(nSplit > 7) call chr_cast(lnSplit(8), iParams(2), cErr)

    dist_needsRnd = .true.

  case("UNIFORM")
    if(nSplit /= 5 .and. nSplit /= 7) then
      write(lerr,"(a,i0)") "DIST> ERROR FILL format UNIFORM expected 2 or 4 arguments, got ",nSplit-3
      write(lerr,"(a)")    "DIST>       FILL column UNIFORM lower upper [first last]"
      iErr = .true.
      return
    end if

    allocate(fParams(2))
    allocate(iParams(2))

    fillMethod = dist_fillUNIFORM
    fParams(1) = zero
    fParams(2) = one
    iParams(1) =  1
    iParams(2) = -1

    if(nSplit > 3) call chr_cast(lnSplit(4), fParams(1), cErr)
    if(nSplit > 4) call chr_cast(lnSplit(5), fParams(2), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), iParams(1), cErr)
    if(nSplit > 6) call chr_cast(lnSplit(7), iParams(2), cErr)

    dist_needsRnd = .true.

  case("LINEAR")
    if(nSplit /= 5 .and. nSplit /= 7) then
      write(lerr,"(a,i0)") "DIST> ERROR FILL format LINEAR expected 2 or 4 arguments, got ",nSplit-3
      write(lerr,"(a)")    "DIST>       FILL column LINEAR lower upper [first last]"
      iErr = .true.
      return
    end if

    allocate(fParams(2))
    allocate(iParams(2))

    fillMethod = dist_fillLINEAR
    fParams(1) = zero
    fParams(2) = one
    iParams(1) =  1
    iParams(2) = -1

    if(nSplit > 3) call chr_cast(lnSplit(4), fParams(1), cErr)
    if(nSplit > 4) call chr_cast(lnSplit(5), fParams(2), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), iParams(1), cErr)
    if(nSplit > 6) call chr_cast(lnSplit(7), iParams(2), cErr)

  case("COUNT")
    if(nSplit /= 5 .and. nSplit /= 7) then
      write(lerr,"(a,i0)") "DIST> ERROR FILL format COUNT expected 2 or 4 arguments, got ",nSplit-3
      write(lerr,"(a)")    "DIST>       FILL column COUNT start step [first last]"
      iErr = .true.
      return
    end if

    allocate(fParams(1))
    allocate(iParams(4))

    fillMethod = dist_fillCOUNT
    fParams(1) =  zero
    iParams(1) =  1
    iParams(2) = -1
    iParams(3) =  1
    iParams(4) =  1

    if(nSplit > 3) call chr_cast(lnSplit(4), iParams(3), cErr)
    if(nSplit > 4) call chr_cast(lnSplit(5), iParams(4), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), iParams(1), cErr)
    if(nSplit > 6) call chr_cast(lnSplit(7), iParams(2), cErr)

  case default
    write(lerr,"(a)") "DIST> ERROR Unknown FILL method '"//trim(lnSplit(3))//"'"
    iErr = .true.
    return

  end select

  if(iParams(1) < 1) then
    write(lerr,"(a)") "DIST> ERROR First particle cannot be smaller than 1"
    iErr = .true.
    return
  end if
  if(iParams(2) /= -1 .and. iParams(1) > iParams(2)) then
    write(lerr,"(a)") "DIST> ERROR First particle to fill cannot be after last particle"
    iErr = .true.
    return
  end if

  dist_fillList(dist_nFill)%fillName   = trim(lnSplit(2))
  dist_fillList(dist_nFill)%fillTarget = fillTarget
  dist_fillList(dist_nFill)%fillMethod = fillMethod
  dist_fillList(dist_nFill)%iParams    = iParams
  dist_fillList(dist_nFill)%fParams    = fParams

end subroutine dist_parseFill

! ================================================================================================ !
!  Parse File Column Formats
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-05
!  Updated: 2019-07-10
!
!    This routine splits the unit from the format description and adds the format and the scaling
!  factor to the list of format columns for later file parsing.
! ================================================================================================ !
subroutine dist_parseColumn(fmtName, fErr, fmtID, fmtFac, fmtCol, isValid)

  use crcoall
  use string_tools
  use numerical_constants

  character(len=*), intent(in)  :: fmtName
  logical,          intent(out) :: fErr
  integer,          intent(out) :: fmtID
  real(kind=fPrec), intent(out) :: fmtFac
  integer,          intent(out) :: fmtCol
  logical,          intent(out) :: isValid

  character(len=20) fmtBase
  character(len=10) fmtUnit
  real(kind=fPrec)  uFac
  integer           c, unitPos, fmtLen

  fErr    = .false.
  fmtID   = dist_fmtNONE
  fmtFac  = one
  fmtCol  = 0
  isValid = .true.

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

  case("SKIP")
    fmtID = dist_fmtNONE
  case("ID")
    fmtID = dist_fmtPartID
  case("PARENT")
    fmtID  = dist_fmtParentID

  case("X")
    call dist_unitScale(fmtName, fmtUnit, 1, fmtFac, fErr)
    fmtID  = dist_fmtX
    fmtCol = 1
  case("Y")
    call dist_unitScale(fmtName, fmtUnit, 1, fmtFac, fErr)
    fmtID  = dist_fmtY
    fmtCol = 3

  case("XP")
    call dist_unitScale(fmtName, fmtUnit, 2, fmtFac, fErr)
    fmtID  = dist_fmtXP
    fmtCol = 2
  case("YP")
    call dist_unitScale(fmtName, fmtUnit, 2, fmtFac, fErr)
    fmtID  = dist_fmtYP
    fmtCol = 4

  case("PX")
    call dist_unitScale(fmtName, fmtUnit, 3, fmtFac, fErr)
    fmtID  = dist_fmtPX
    fmtCol = 2
  case("PY")
    call dist_unitScale(fmtName, fmtUnit, 3, fmtFac, fErr)
    fmtID  = dist_fmtPY
    fmtCol = 4

  case("PX/P0","PXP0")
    fmtID  = dist_fmtPX
    fmtCol = 2
  case("PY/P0","PYP0")
    fmtID  = dist_fmtPY
    fmtCol = 4

  case("SIGMA")
    call dist_unitScale(fmtName, fmtUnit, 1, fmtFac, fErr)
    fmtID  = dist_fmtSIGMA
    fmtCol = 5
  case("ZETA")
    call dist_unitScale(fmtName, fmtUnit, 1, fmtFac, fErr)
    fmtID  = dist_fmtZETA
    fmtCol = 5
  case("DT")
    call dist_unitScale(fmtName, fmtUnit, 4, fmtFac, fErr)
    fmtID  = dist_fmtDT
    fmtCol = 5

  case("E")
    call dist_unitScale(fmtName, fmtUnit, 3, fmtFac, fErr)
    fmtID  = dist_fmtE
    fmtCol = 6
  case("P")
    call dist_unitScale(fmtName, fmtUnit, 3, fmtFac, fErr)
    fmtID  = dist_fmtP
    fmtCol = 6
  case("DE/E0","DEE0")
    fmtID  = dist_fmtDEE0
    fmtCol = 6
  case("DP/P0","DPP0","DELTA")
    fmtID  = dist_fmtDPP0
    fmtCol = 6
  case("PT")
    fmtID  = dist_fmtPT
    fmtCol = 6
  case("PSIGMA")
    fmtID  = dist_fmtPSIGMA
    fmtCol = 6

  case("XN")
    fmtID  = dist_fmtXN
    fmtCol = 1
    dist_distLib = dist_distLibNorm
  case("YN")
    fmtID  = dist_fmtYN
    fmtCol = 3
    dist_distLib = dist_distLibNorm
  case("ZN")
    fmtID  = dist_fmtZN
    fmtCol = 5
    dist_distLib = dist_distLibNorm
  case("PXN")
    fmtID  = dist_fmtPXN
    fmtCol = 2
    dist_distLib = dist_distLibNorm
  case("PYN")
    fmtID  = dist_fmtPYN
    fmtCol = 4
    dist_distLib = dist_distLibNorm
  case("PZN")
    fmtID  = dist_fmtPZN
    fmtCol = 6
    dist_distLib = dist_distLibNorm

  case("JX")
    fmtID  = dist_fmtJX
    fmtCol = 1
    dist_distLib = .true.
  case("JY")
    fmtID  = dist_fmtJY
    fmtCol = 3
    dist_distLib = .true.
  case("JZ")
    fmtID  = dist_fmtJZ
    fmtCol = 5
    dist_distLib = .true.
  case("PHIX")
    call dist_unitScale(fmtName, fmtUnit, 2, fmtFac, fErr)
    fmtID  = dist_fmtPhiX
    fmtCol = 2
    dist_distLib = .true.
  case("PHIY")
    call dist_unitScale(fmtName, fmtUnit, 2, fmtFac, fErr)
    fmtID  = dist_fmtPhiY
    fmtCol = 4
    dist_distLib = .true.
  case("PHIZ")
    call dist_unitScale(fmtName, fmtUnit, 2, fmtFac, fErr)
    fmtID  = dist_fmtPhiZ
    fmtCol = 6
    dist_distLib = .true.

  case("MASS","M")
    call dist_unitScale(fmtName, fmtUnit, 3, fmtFac, fErr)
    fmtID  = dist_fmtMASS
  case("CHARGE","Q")
    fmtID  = dist_fmtCHARGE
  case("ION_A")
    fmtID  = dist_fmtIonA
  case("ION_Z")
    fmtID  = dist_fmtIonZ
  case("PDGID")
    fmtID  = dist_fmtPDGID

  case("SX")
    fmtID  = dist_fmtSPINX
  case("SY")
    fmtID  = dist_fmtSPINY
  case("SZ")
    fmtID  = dist_fmtSPINZ

  case default
    isValid = .false.
    return

  end select

#ifndef DISTLIB
  if(dist_distLib) then
    write(lerr,"(a)") "DIST> ERROR Format '"//trim(fmtName)//"' requires SixTrack to be built with DISTLIB enabled"
    fErr = .true.
    return
  end if
#endif

end subroutine dist_parseColumn

! ================================================================================================ !
!  Parse Multi-Column Formats
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-05
!  Updated: 2019-08-13
! ================================================================================================ !
subroutine dist_setMultiColFormat(fmtName, isValid)

  use string_tools
  use numerical_constants

  character(len=*), intent(in)  :: fmtName
  logical,          intent(out) :: isValid

  isValid = .false.

  select case(chr_toUpper(fmtName))

  case("4D") ! 4D default coordinates
    call dist_appendFormat(dist_fmtX,        one,  1)
    call dist_appendFormat(dist_fmtPX,       one,  2)
    call dist_appendFormat(dist_fmtY,        one,  3)
    call dist_appendFormat(dist_fmtPY,       one,  4)
    isValid = .true.

  case("6D") ! 6D default coordinates
    call dist_appendFormat(dist_fmtX,        one,  1)
    call dist_appendFormat(dist_fmtPX,       one,  2)
    call dist_appendFormat(dist_fmtY,        one,  3)
    call dist_appendFormat(dist_fmtPY,       one,  4)
    call dist_appendFormat(dist_fmtZETA,     one,  5)
    call dist_appendFormat(dist_fmtDPP0,     one,  6)
    isValid = .true.

  case("NORM") ! 6D normalised coordinates
    call dist_appendFormat(dist_fmtXN,       one,  1)
    call dist_appendFormat(dist_fmtPXN,      one,  2)
    call dist_appendFormat(dist_fmtYN,       one,  3)
    call dist_appendFormat(dist_fmtPYN,      one,  4)
    call dist_appendFormat(dist_fmtZN,       one,  5)
    call dist_appendFormat(dist_fmtPZN,      one,  6)
    dist_distLib = dist_distLibNorm
    isValid = .true.

  case("ACTION") ! 6D action
    call dist_appendFormat(dist_fmtJX,       one,  1)
    call dist_appendFormat(dist_fmtPhiX,     one,  2)
    call dist_appendFormat(dist_fmtJY,       one,  3)
    call dist_appendFormat(dist_fmtPhiY,     one,  4)
    call dist_appendFormat(dist_fmtJZ,       one,  5)
    call dist_appendFormat(dist_fmtPhiZ,     one,  6)
    dist_distLib = .true.
    isValid = .true.

  case("IONS") ! The ion columns
    call dist_appendFormat(dist_fmtMASS,     c1e3, 0)
    call dist_appendFormat(dist_fmtCHARGE,   one,  0)
    call dist_appendFormat(dist_fmtIonA,     one,  0)
    call dist_appendFormat(dist_fmtIonZ,     one,  0)
    call dist_appendFormat(dist_fmtPDGID,    one,  0)
    isValid = .true.

  case("SPIN") ! The spin columns
    call dist_appendFormat(dist_fmtSPINX,    one,  0)
    call dist_appendFormat(dist_fmtSPINY,    one,  0)
    call dist_appendFormat(dist_fmtSPINZ,    one,  0)
    isValid = .true.

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
    isValid = .true.
  end select

end subroutine dist_setMultiColFormat

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
    case("[1000]")
      unitScale = c1e3
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
    case("[1000]")
      unitScale = c1e3
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
      call dist_multFormatError(partCol)
      call prror
    end if
  end if

end subroutine dist_appendFormat

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-05
!  Updated: 2019-07-05
! ================================================================================================ !
subroutine dist_multFormatError(partCol)

  use crcoall

  integer, intent(in) :: partCol

  write(lerr,"(a,i0)") "DIST> ERROR Multiple formats selected for particle coordinate ",partCol
  select case(partCol)
  case(1)
    write(lerr,"(a)") "DIST>      Choose only one of: X, XN, JX"
  case(2)
    write(lerr,"(a)") "DIST>      Choose only one of: XP, PX, PX/P0, PXN, PHIX"
  case(3)
    write(lerr,"(a)") "DIST>      Choose only one of: Y, YN, JY"
  case(4)
    write(lerr,"(a)") "DIST>      Choose only one of: YP, PY, PY/P0, PYN, PHIY"
  case(5)
    write(lerr,"(a)") "DIST>      Choose only one of: SIGMA, ZETA, DT, ZN, JZ"
  case(6)
    write(lerr,"(a)") "DIST>      Choose only one of: E, P, DE/E0, DP/P0, PT, PZN, PHIZ, PSIGMA"
  end select

end subroutine dist_multFormatError

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

  spErr = .false.
  cErr  = .false.

  if(dist_nParticle > 0) then
    do j=1,dist_nParticle
      call chr_split(dist_partLine(j), lnSplit, nSplit, spErr)
      if(spErr) then
        write(lerr,"(a,i0,a)") "DIST> ERROR Could not split PARTICLE definition number ",j," from "//trim(fort3)
        call prror
      end if
      do i=2,nSplit
        call dist_saveParticle(j, i-1, lnSplit(i), cErr)
        if(cErr) then
          write(lerr,"(a,2(i0,a))") "DIST> ERROR Could not parse PARTICLE definition number ",j,", column ",i," from "//trim(fort3)
          call prror
        end if
      end do
      dist_numPart = dist_numPart + 1
    end do
  end if

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
  real(kind=fPrec) fmtFac
  integer i, nSplit, lineNo, fUnit, nRead, fmtID, fmtCol
  logical spErr, cErr, isValid

  write(lout,"(a)") "DIST> Reading particles from '"//trim(dist_distFile)//"'"

  nRead  = 0
  lineNo = 0
  cErr   = .false.

  if(dist_hasFormat .eqv. .false.) then
    call dist_setMultiColFormat("OLD_DIST",isValid)
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

  nRead = nRead + 1
  if(nRead > napx) then
    write(lout,"(a,i0,a)") "DIST> Stopping reading file as ",napx," particles have been read, as requested in '"//trim(fort3)//"'"
    goto 40
  end if

  dist_numPart = dist_numPart + 1
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
  write(lout,"(a,i0,a)") "DIST> Read ",dist_numPart," particles from file '"//trim(dist_distFile)//"'"

end subroutine dist_readDist

! ================================================================================================ !
!  Save Particle Data to Arrays
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-08
!  Updated: 2019-07-13
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
    dist_readPartID = .true.

  case(dist_fmtParentID)
    call chr_cast(inVal, parentID(partNo), sErr)
    dist_readParentID = .true.

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

  case(dist_fmtE, dist_fmtDEE0, dist_fmtPT, dist_fmtPSIGMA, dist_fmtP, dist_fmtDPP0, dist_fmtPZN, dist_fmtPhiZ)
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

  !  Spin Columns
  ! ==============

  case(dist_fmtSPINX)
    call chr_cast(inVal, spin_x(partNo), sErr)

  case(dist_fmtSPINY)
    call chr_cast(inVal, spin_y(partNo), sErr)

  case(dist_fmtSPINZ)
    call chr_cast(inVal, spin_z(partNo), sErr)

  end select

end subroutine dist_saveParticle

! ================================================================================================ !
!  Run the FILLs
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-30
!  Updated: 2019-07-30
! ================================================================================================ !
subroutine dist_doFills

  use crcoall
  use mod_common
  use mod_common_main
  use, intrinsic :: iso_fortran_env, only : int16

  integer iFill, j, iA, iB, iVal, fM, fT, intVal

  if(dist_nFill == 0) then
    return
  end if

  ! If we have fills, we consider all particles generated
  dist_numPart = napx

  do iFill=1,dist_nFill

    iA = dist_fillList(iFill)%iParams(1)
    iB = dist_fillList(iFill)%iParams(2)
    if(iB == -1) then
      iB = napx
    end if

    fM = dist_fillList(iFill)%fillMethod
    fT = dist_fillList(iFill)%fillTarget

    select case(fT)

    case(dist_fmtNONE)
      return

    case(dist_fmtPartID)
      dist_readPartID = .true.
      if(fM == dist_fillCOUNT) then
        iVal = dist_fillList(iFill)%iParams(3)
        do j=iA,iB
          partID(j) = iVal
          iVal = iVal + dist_fillList(iFill)%iParams(4)
        end do
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be COUNT"
        call prror
      end if

    case(dist_fmtParentID)
      write(lerr,"(a)") "DIST> ERROR Variable "//trim(dist_fillList(iFill)%fillName)//" cannot be filled automatically"
      call prror

    !  Horizontal Coordinates
    ! ========================

    case(dist_fmtX, dist_fmtXN)
      if(fM == dist_fillFLOAT .or. fM == dist_fillGAUSS .or. fM == dist_fillUNIFORM .or. fM == dist_fillLINEAR) then
        call dist_fillThis(dist_partCol1, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT, GAUSS, UNIFORM or LINEAR"
        call prror
      end if

    case(dist_fmtPX, dist_fmtXP, dist_fmtPXP0, dist_fmtPXN)
      if(fM == dist_fillFLOAT .or. fM == dist_fillGAUSS .or. fM == dist_fillUNIFORM .or. fM == dist_fillLINEAR) then
        call dist_fillThis(dist_partCol2, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT, GAUSS, UNIFORM or LINEAR"
        call prror
      end if

    !  Vertical Coordinates
    ! ========================

    case(dist_fmtY, dist_fmtYN)
      if(fM == dist_fillFLOAT .or. fM == dist_fillGAUSS .or. fM == dist_fillUNIFORM .or. fM == dist_fillLINEAR) then
        call dist_fillThis(dist_partCol3, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT, GAUSS, UNIFORM or LINEAR"
        call prror
      end if

    case(dist_fmtPY, dist_fmtYP, dist_fmtPYP0, dist_fmtPYN)
      if(fM == dist_fillFLOAT .or. fM == dist_fillGAUSS .or. fM == dist_fillUNIFORM .or. fM == dist_fillLINEAR) then
        call dist_fillThis(dist_partCol4, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT, GAUSS, UNIFORM or LINEAR"
        call prror
      end if

    !  Longitudinal Coordinates
    ! ==========================

    case(dist_fmtSIGMA, dist_fmtZETA, dist_fmtDT, dist_fmtZN)
      if(fM == dist_fillFLOAT .or. fM == dist_fillGAUSS .or. fM == dist_fillUNIFORM .or. fM == dist_fillLINEAR) then
        call dist_fillThis(dist_partCol5, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT, GAUSS, UNIFORM or LINEAR"
        call prror
      end if

    case(dist_fmtE, dist_fmtDEE0, dist_fmtPT, dist_fmtPSIGMA, dist_fmtP, dist_fmtDPP0, dist_fmtPZN)
      if(fM == dist_fillFLOAT .or. fM == dist_fillGAUSS .or. fM == dist_fillUNIFORM .or. fM == dist_fillLINEAR) then
        call dist_fillThis(dist_partCol6, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT, GAUSS, UNIFORM or LINEAR"
        call prror
      end if

    !  Action Angle Special Case
    ! ===========================

    case(dist_fmtJX, dist_fmtJY, dist_fmtJZ)
      if(fM == dist_fillFLOAT .or. fM == dist_fillRAYLEIGH .or. fM == dist_fillUNIFORM .or. fM == dist_fillLINEAR) then
        call dist_fillThis(dist_partCol1, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT, RAYLEIGH, UNIFORM or LINEAR"
        call prror
      end if

    case(dist_fmtPhiX, dist_fmtPhiY, dist_fmtPhiZ)
      if(fM == dist_fillFLOAT .or. fM == dist_fillUNIFORM .or. fM == dist_fillLINEAR) then
        call dist_fillThis(dist_partCol1, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT, UNIFORM or LINEAR"
        call prror
      end if

    !  Ion Columns
    ! =============

    case(dist_fmtMASS)
      if(fM == dist_fillFLOAT) then
        call dist_fillThis(nucm, iA, iB, fM, dist_fillList(iFill)%fParams)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be FLOAT"
        call prror
      end if
      dist_readMass = .true.

    case(dist_fmtCHARGE)
      if(fM == dist_fillINT) then
        nqq(iA:iB) = int(dist_fillList(iFill)%iParams(3), kind=int16)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be INT"
        call prror
      end if
      dist_readCharge = .true.

    case(dist_fmtIonA)
      if(fM == dist_fillINT) then
        naa(iA:iB) = int(dist_fillList(iFill)%iParams(3), kind=int16)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be INT"
        call prror
      end if
      dist_readIonA = .true.

    case(dist_fmtIonZ)
      if(fM == dist_fillINT) then
        nzz(iA:iB) = int(dist_fillList(iFill)%iParams(3), kind=int16)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be INT"
        call prror
      end if
      dist_readIonZ = .true.

    case(dist_fmtPDGID)
      if(fM == dist_fillINT) then
        pdgid(iA:iB) = dist_fillList(iFill)%iParams(3)
      else
        write(lerr,"(a)") "DIST> ERROR FILL "//trim(dist_fillList(iFill)%fillName)//" must be INT"
        call prror
      end if
      dist_readPDGID = .true.

    !  Spin Columns
    ! ==============

    case(dist_fmtSPINX, dist_fmtSPINY, dist_fmtSPINZ)
      write(lerr,"(a)") "DIST> ERROR Variable "//trim(dist_fillList(iFill)%fillName)//" cannot be filled automatically"
      call prror

    end select
  end do

end subroutine dist_doFills

! ================================================================================================ !
!  Fill Arrays
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-30
!  Updated: 2019-07-30
! ================================================================================================ !
subroutine dist_fillThis(fillArray, iA, iB, fillMethod, fParams)

  use mod_random

  real(kind=fPrec), intent(inout) :: fillArray(*)
  integer,          intent(in)    :: iA, iB
  integer,          intent(in)    :: fillMethod
  real(kind=fPrec), intent(in)    :: fParams(*)

  integer nRnd, j
  real(kind=fPrec) rndVals(iB-iA+1), dStep

  nRnd = iB-iA+1

  select case(fillMethod)

  case(dist_fillFLOAT)
    fillArray(iA:iB) = fParams(1)

  case(dist_fillGAUSS)
    if(fParams(3) > zero) then
      call rnd_normal(rndser_distGen, rndVals, nRnd, fParams(3))
    else
      call rnd_normal(rndser_distGen, rndVals, nRnd)
    end if
    fillArray(iA:iB) = rndVals*fParams(1) + fParams(2)

  case(dist_fillRAYLEIGH)
    if(fParams(2) > zero) then
      if(fParams(3) > zero) then
        call rnd_rayleigh(rndser_distGen, rndVals, nRnd, fParams(2), fParams(3))
      else
        call rnd_rayleigh(rndser_distGen, rndVals, nRnd, fParams(2))
      end if
    else
      call rnd_rayleigh(rndser_distGen, rndVals, nRnd)
    end if
    fillArray(iA:iB) = rndVals*fParams(1)

  case(dist_fillUNIFORM)
    call rnd_uniform(rndser_distGen, rndVals, nRnd)
    fillArray(iA:iB) = rndVals*(fParams(2)-fParams(1)) + fParams(1)

  case(dist_fillLINEAR)
    dStep = (fParams(2)-fParams(1))/real(iB-iA,kind=fPrec)
    do j=iA,iB
      fillArray(j) = fParams(1) + dStep*real(j-iA,kind=fPrec)
    end do

  end select

end subroutine dist_fillThis

! ================================================================================================ !
!  Post-Processing of Particle Arrays
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-05
!  Updated: 2019-07-13
! ================================================================================================ !
subroutine dist_postprParticles

  use crcoall
  use mod_common
  use mod_particles
  use mod_common_main
  use physical_constants
  use numerical_constants

  logical doAction, doNormal

  doAction = .false.
  doNormal = .false.

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
    ejv(1:napx) = dist_partCol6(1:napx)*e0f + e0
    call part_updatePartEnergy(1,.false.)
  case(dist_fmtPSIGMA)
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
    sigmv(1:napx) = dist_partCol5(1:napx)*rvv(1:napx)
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

  if(doNormal .and. doAction) then
    write(lerr,"(a)") "DIST> ERROR Cannot mix normalised and action coordinates"
    call prror
  end if

#ifdef DISTLIB
  if(doNormal .and. dist_distLibNorm) then
    call distlib_setCoords(                        &
      dist_partCol1, dist_partCol2, dist_partCol3, &
      dist_partCol4, dist_partCol5, dist_partCol6, &
      napx, dist_coordTypeNormal                   &
    )
  end if

  if(doAction .and. dist_distLibNorm) then
    call distlib_setCoords(                        &
      dist_partCol1, dist_partCol2, dist_partCol3, &
      dist_partCol4, dist_partCol5, dist_partCol6, &
      napx, dist_coordTypeAction                   &
    )
  end if
#endif

  if(doNormal .and. .not. dist_distLibNorm) then
    call dist_normToPhysical
  end if

end subroutine dist_postprParticles

! ================================================================================================ !
!  Apply normalisation matrix
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-07-30
!  Updated: 2019-07-30
! ================================================================================================ !
subroutine dist_normToPhysical

  use mod_common
  use mod_particles
  use mod_common_main
  use numerical_constants

  real(kind=fPrec) sqEmitX, sqEmitY, sqEmitZ

  sqEmitX = sqrt(dist_beamEmit(1))
  sqEmitY = sqrt(dist_beamEmit(2))
  sqEmitZ = sqrt(dist_beamEmit(3))

  xv1(1:napx) = ( &
    dist_partCol1(1:napx) * (sqEmitX * dist_tMat(1,1))  + &
    dist_partCol2(1:napx) * (sqEmitX * dist_tMat(1,2))  + &
    dist_partCol3(1:napx) * (sqEmitY * dist_tMat(1,3))  + &
    dist_partCol4(1:napx) * (sqEmitY * dist_tMat(1,4))  + &
    dist_partCol5(1:napx) * (sqEmitZ * dist_tMat(1,5))  + &
    dist_partCol6(1:napx) * (sqEmitZ * dist_tMat(1,6))) * c1e3

  yv1(1:napx) = ( &
    dist_partCol1(1:napx) * (sqEmitX * dist_tMat(2,1))  + &
    dist_partCol2(1:napx) * (sqEmitX * dist_tMat(2,2))  + &
    dist_partCol3(1:napx) * (sqEmitY * dist_tMat(2,3))  + &
    dist_partCol4(1:napx) * (sqEmitY * dist_tMat(2,4))  + &
    dist_partCol5(1:napx) * (sqEmitZ * dist_tMat(2,5))  + &
    dist_partCol6(1:napx) * (sqEmitZ * dist_tMat(2,6))) * (c1e3*(one+dpsv(1:napx)))

  xv2(1:napx) = ( &
    dist_partCol1(1:napx) * (sqEmitX * dist_tMat(3,1))  + &
    dist_partCol2(1:napx) * (sqEmitX * dist_tMat(3,2))  + &
    dist_partCol3(1:napx) * (sqEmitY * dist_tMat(3,3))  + &
    dist_partCol4(1:napx) * (sqEmitY * dist_tMat(3,4))  + &
    dist_partCol5(1:napx) * (sqEmitZ * dist_tMat(3,5))  + &
    dist_partCol6(1:napx) * (sqEmitZ * dist_tMat(3,6))) * c1e3

  yv2(1:napx) = ( &
    dist_partCol1(1:napx) * (sqEmitX * dist_tMat(4,1))  + &
    dist_partCol2(1:napx) * (sqEmitX * dist_tMat(4,2))  + &
    dist_partCol3(1:napx) * (sqEmitY * dist_tMat(4,3))  + &
    dist_partCol4(1:napx) * (sqEmitY * dist_tMat(4,4))  + &
    dist_partCol5(1:napx) * (sqEmitZ * dist_tMat(4,5))  + &
    dist_partCol6(1:napx) * (sqEmitZ * dist_tMat(4,6))) * (c1e3*(one+dpsv(1:napx)))

  sigmv(1:napx) = ( &
    dist_partCol1(1:napx) * (sqEmitX * dist_tMat(5,1))  + &
    dist_partCol2(1:napx) * (sqEmitX * dist_tMat(5,2))  + &
    dist_partCol3(1:napx) * (sqEmitY * dist_tMat(5,3))  + &
    dist_partCol4(1:napx) * (sqEmitY * dist_tMat(5,4))  + &
    dist_partCol5(1:napx) * (sqEmitZ * dist_tMat(5,5))  + &
    dist_partCol6(1:napx) * (sqEmitZ * dist_tMat(5,6))) / rvv(1:napx)

  dpsv(1:napx) = ( &
    dist_partCol1(1:napx) * (sqEmitX * dist_tMat(6,1))  + &
    dist_partCol2(1:napx) * (sqEmitX * dist_tMat(6,2))  + &
    dist_partCol3(1:napx) * (sqEmitY * dist_tMat(6,3))  + &
    dist_partCol4(1:napx) * (sqEmitY * dist_tMat(6,4))  + &
    dist_partCol5(1:napx) * (sqEmitZ * dist_tMat(6,5))  + &
    dist_partCol6(1:napx) * (sqEmitZ * dist_tMat(6,6)))

  call part_updatePartEnergy(3,.true.)

end subroutine dist_normToPhysical

! ================================================================================================ !
!  A. Mereghetti and D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Rewritten: 2019-07-09
!  Updated:   2019-08-13
! ================================================================================================ !
subroutine dist_finaliseDist

  use parpro
  use mod_pdgid
  use mod_common
  use mod_particles
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

  if(dist_readParentID .and. .not. dist_readPartID) then
    write(lerr,"(a)") "DIST> ERROR If you set particle parent ID you must also set particle ID"
    call prror
  end if

  if(dist_readPartID) then
    if(dist_readParentID .eqv. .false.) then
      ! Particle parentID is not set, so we set its parentID to itself
      parentID(1:napx) = partID(1:napx)
    end if
    call part_setPairID ! Set the pairID only
  else
    ! No IDs provided, so we set them to the default range 1:napx and compute their corresponding pairID
    call part_setParticleID ! Set partID, parentID and pairID
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
  write(fUnit,"(a)")          "# x[mm], xp[mrad], y[mm], yp[mrad], sigmv[mm], ejfv[MeV/c]"
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
