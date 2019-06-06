! ================================================================================================ !
!  SixTrack SCATTER Module
!  V.K. Berglyd Olsen, K.N. Sjobak, H. Burkhardt, BE-ABP-HSS
!  Last modified: 2018-11-12
!
!  References:
!  "Elastic pp scattering estimates and simulation relevant for burn-off"
!    > https://indico.cern.ch/event/625576
!  "Elastic Scattering in SixTrack"
!    > https://indico.cern.ch/event/737429
! ================================================================================================ !
module scatter

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_ranecu
  use strings
  use string_tools
#ifdef HDF5
  use hdf5_output
#endif
#ifdef PYTHIA
  use mod_pythia
#endif

  implicit none

  ! Common variables for the SCATTER routines
  logical, public, save :: scatter_active      = .false.
  logical, public, save :: scatter_debug       = .false.
  logical, public, save :: scatter_allowLosses = .false.

  ! Scatter Parameters
  integer,          parameter :: scatter_nProc            = 9
  integer,          parameter :: scatter_idAbsorb         = 1
  integer,          parameter :: scatter_idNonDiff        = 2
  integer,          parameter :: scatter_idElastic        = 3
  integer,          parameter :: scatter_idSingleDiffXB   = 4
  integer,          parameter :: scatter_idSingleDiffAX   = 5
  integer,          parameter :: scatter_idDoubleDiff     = 6
  integer,          parameter :: scatter_idCentralDiff    = 7
  integer,          parameter :: scatter_idUnknown        = 8
  integer,          parameter :: scatter_idError          = 9
  character(len=8), parameter :: scatter_procNames(scatter_nProc) = (/  &
    "Absorbed","NonDiff ","Elastic ","SingD_XB","SingD_AX",             &
    "DoubD_XX","CentDiff","Unknown ","Error   "                         &
  /)

  ! Storage Structs
  type, private :: scatter_elemStore
    character(len=mNameLen)       :: bezName
    integer                       :: bezID
    real(kind=fPrec)              :: elemScale
    real(kind=fPrec)              :: sigmaTot
    real(kind=fPrec)              :: ratioTot
    logical                       :: autoRatio
    integer                       :: profileID
    integer,          allocatable :: generatorID(:)
    real(kind=fPrec), allocatable :: brRatio(:)
  end type scatter_elemStore

  type, private :: scatter_proStore
    character(len=:), allocatable :: proName
    integer                       :: proType
    real(kind=fPrec), allocatable :: fParams(:)
  end type scatter_proStore

  type, private :: scatter_genStore
    character(len=:), allocatable :: genName
    integer                       :: genType
    real(kind=fPrec)              :: crossSection
    real(kind=fPrec), allocatable :: fParams(:)
  end type scatter_genStore

  ! Storage Arrays
  type(scatter_elemStore), allocatable, private, save :: scatter_elemList(:)
  type(scatter_proStore),  allocatable, private, save :: scatter_proList(:)
  type(scatter_genStore),  allocatable, private, save :: scatter_genList(:)
  integer,                              private, save :: scatter_nElem = 0
  integer,                              private, save :: scatter_nPro  = 0
  integer,                              private, save :: scatter_nGen  = 0
  integer,                 allocatable, public,  save :: scatter_elemPointer(:) ! (nele)
  real(kind=fPrec),        allocatable, private, save :: scatter_statScale(:)   ! (npart)

  ! Random generator seeds
  integer, public,  save :: scatter_seed1      = -1
  integer, public,  save :: scatter_seed2      = -1

  ! Variable for file output
  integer, private, save :: scatter_logFile    = -1
  integer, private, save :: scatter_sumFile    = -1
#ifdef HDF5
  integer, private, save :: scatter_logDataSet = 0
  integer, private, save :: scatter_logFormat  = 0
  integer, private, save :: scatter_sumDataSet = 0
  integer, private, save :: scatter_sumFormat  = 0
#endif

#ifdef CR
  integer, public,  save :: scatter_logFilePos     = -1
  integer, public,  save :: scatter_logFilePos_CR  =  0
  integer, public,  save :: scatter_sumFilePos     = -1
  integer, public,  save :: scatter_sumFilePos_CR  =  0

  integer, private, save :: scatter_seed1_CR       = -1
  integer, private, save :: scatter_seed2_CR       = -1

  real(kind=fPrec), allocatable, private, save :: scatter_statScale_CR(:)
#endif

  public :: scatter_getScaling

contains

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-12
! =================================================================================================
subroutine scatter_expand_arrays(nele_new, npart_new)
  use mod_alloc
  implicit none
  integer, intent(in) :: nele_new, npart_new
  call alloc(scatter_elemPointer, nele_new,  0,   "scatter_elemPointer")
  call alloc(scatter_statScale,   npart_new, one, "scatter_statScale")
end subroutine scatter_expand_arrays

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-11-12
! =================================================================================================
subroutine scatter_init

  use crcoall
  use parpro
  use mod_units

  implicit none

  integer iError
  logical fErr
  real(kind=fPrec) sigmaTot
#ifdef HDF5
  type(h5_dataField), allocatable :: setFields(:)
#endif

  if(scatter_active .eqv. .false.) return

  ! Initialise data output
#ifdef HDF5
  if(h5_useForSCAT) then

    call h5_initForScatter

    ! Scatter Log
    allocate(setFields(14))

    setFields(1)  = h5_dataField(name="ID",      type=h5_typeInt)
    setFields(2)  = h5_dataField(name="TURN",    type=h5_typeInt)
    setFields(3)  = h5_dataField(name="BEZ",     type=h5_typeChar, size=mNameLen)
    setFields(4)  = h5_dataField(name="GEN",     type=h5_typeChar, size=mNameLen)
    setFields(5)  = h5_dataField(name="PROCESS", type=h5_typeChar, size=8)
    setFields(6)  = h5_dataField(name="LOST",    type=h5_typeInt)
    setFields(7)  = h5_dataField(name="T",       type=h5_typeReal)
    setFields(8)  = h5_dataField(name="THETA",   type=h5_typeReal)
    setFields(9)  = h5_dataField(name="PHI",     type=h5_typeReal)
    setFields(10) = h5_dataField(name="dE/E",    type=h5_typeReal)
    setFields(11) = h5_dataField(name="dP/P",    type=h5_typeReal)
    setFields(12) = h5_dataField(name="DENSITY", type=h5_typeReal)
    setFields(13) = h5_dataField(name="PROB",    type=h5_typeReal)
    setFields(14) = h5_dataField(name="STATCORR",type=h5_typeReal)

    call h5_createFormat("scatter_log_fmt", setFields, scatter_logFormat)
    call h5_createDataSet("scatter_log", h5_scatID, scatter_logFormat, scatter_logDataSet)
    block
      character(len=:), allocatable :: colNames(:)
      character(len=:), allocatable :: colUnits(:)
      logical spErr
      integer nSplit
      call chr_split("ID turn bez generator process lost t dE/E dP/P theta phi density prob stat_corr",colNames,nSplit,spErr)
      call chr_split("1 1 text text text 1 MeV^2 1 1 mrad rad mb^-1 1 1",colUnits,nSplit,spErr)
      call h5_writeDataSetAttr(scatter_logDataSet,"colNames",colNames)
      call h5_writeDataSetAttr(scatter_logDataSet,"colUnits",colUnits)
    end block

    deallocate(setFields)

    ! Scatter Summary
    allocate(setFields(10))

    setFields(1)  = h5_dataField(name="TURN",      type=h5_typeInt)
    setFields(2)  = h5_dataField(name="NAPX",      type=h5_typeInt)
    setFields(3)  = h5_dataField(name="BEZ",       type=h5_typeChar, size=mNameLen)
    setFields(4)  = h5_dataField(name="GEN",       type=h5_typeChar, size=mNameLen)
    setFields(5)  = h5_dataField(name="PROCESS",   type=h5_typeChar, size=8)
    setFields(6)  = h5_dataField(name="NSCAT",     type=h5_typeInt)
    setFields(7)  = h5_dataField(name="NLOST",     type=h5_typeInt)
    setFields(8)  = h5_dataField(name="SCATRATIO", type=h5_typeReal)
    setFields(9)  = h5_dataField(name="CROSSSEC",  type=h5_typeReal)
    setFields(10) = h5_dataField(name="SCALING",   type=h5_typeReal)

    call h5_createFormat("scatter_summary_fmt", setFields, scatter_sumFormat)
    call h5_createDataSet("summary", h5_scatID, scatter_sumFormat, scatter_sumDataSet)
    block
      character(len=:), allocatable :: colNames(:)
      character(len=:), allocatable :: colUnits(:)
      logical spErr
      integer nSplit
      call chr_split("turn napx bez generator process nScat nLost scatRat scaling",colNames,nSplit,spErr)
      call chr_split("1 1 text text text 1 1 1 1",colUnits,nSplit,spErr)
      call h5_writeDataSetAttr(scatter_sumDataSet,"colNames",colNames)
      call h5_writeDataSetAttr(scatter_sumDataSet,"colUnits",colUnits)
    end block

    deallocate(setFields)

    ! Write Attributes
    call h5_writeAttr(h5_scatID,"SEED",scatter_seed1)
  else
#endif

  ! Open scatter_log.dat
  call f_requestUnit("scatter_log.dat", scatter_logFile)
#ifdef CR
  if(scatter_logFilePos == -1) then
    write(crlog,"(a)") "CR_CHECK> SCATTER INIT opening new file 'scatter_log.dat'"
#endif
    call f_open(unit=scatter_logFile,file="scatter_log.dat",formatted=.true.,mode="w",err=fErr,status="replace")
    write(scatter_logFile,"(a)") "# scatter_log"
    write(scatter_logFile,"(a1,a8,1x,a8,2(1x,a20),1x,a8,1x,a4,2(1x,a12),1x,a9,5(1x,a16))") &
      "#","ID","turn",chr_rPad("bez",20),chr_rPad("generator",20),chr_rPad("process",8),&
      "lost","t[MeV^2]","theta[mrad]","phi[rad]","dE/E","dP/P","density","prob","stat_corr"
    flush(scatter_logFile)
#ifdef CR
    scatter_logFilePos = 2
    flush(scatter_logFile)
  else
    write(crlog,"(a)") "CR_CHECK> SCATTER kept already opened file 'scatter_log.dat'"
  end if
#endif
#ifdef HDF5
  end if
#endif

  ! Open scatter_summary.dat
  call f_requestUnit("scatter_summary.dat",scatter_sumFile)
#ifdef CR
  if(scatter_sumFilePos == -1) then
    write(crlog,"(a)") "CR_CHECK> SCATTER INIT opening new file 'scatter_summary.dat'"
#endif
    call f_open(unit=scatter_sumFile,file="scatter_summary.dat",formatted=.true.,mode="w",err=fErr,status="replace")
    call scatter_writeReport
    write(scatter_sumFile,"(a)") "#  Summary Log"
    write(scatter_sumFile,"(a)") "# ============="
    write(scatter_sumFile,"(a1,a8,1x,a8,2(1x,a20),3(1x,a8),1x,a9,2(1x,a13))") &
      "#","turn","napx",chr_rPad("bez",20),chr_rPad("generator",20),chr_rPad("process",8), &
      "nScat","nLost","scatRat","crossSec[mb]","scaling"
    flush(scatter_sumFile)
#ifdef CR
    scatter_sumFilePos = scatter_sumFilePos + 3
    flush(scatter_sumFile)
  else
    write(crlog,"(a)") "CR_CHECK> SCATTER kept already opened file 'scatter_summary.dat'"
  end if
#endif

end subroutine scatter_init

! =================================================================================================
!  BEGIN Input Parser Functions
! =================================================================================================

! =================================================================================================
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-20
! =================================================================================================
subroutine scatter_parseInputLine(inLine, iErr)

  use crcoall

  implicit none

  type(string), intent(in)    :: inLine
  logical,      intent(inout) :: iErr

  type(string), allocatable   :: lnSplit(:)
  type(string) keyWord
  integer      nSplit, i
  logical      spErr

  ! Split the input line
  call str_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lerr,"(a)") "SCATTER> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) then
    if(scatter_debug) then
      write (lout,"(a,i3,a)") "SCATTER> DEBUG Input line len=",len(inLine),": '"//trim(inLine)//"'"
      write (lout,"(a)")      "SCATTER> DEBUG  * No fields found."
    end if
    return
  end if

  keyWord = lnSplit(1)%upper()

  select case(keyWord%get())
  case("DEBUG")
    scatter_debug = .true.
    write(lout,"(a)") "SCATTER> Scatter block debugging is ON."

  case("LOSSES")
    scatter_allowLosses = .true.
    write(lout,"(a)") "SCATTER> Particle losses is ALLOWED."
#ifdef PYTHIA
    ! Pythia also needs to know that losses are allowed to determine which processes can be used
    pythia_allowLosses = .true.
#endif

  case("SEED")
    if(nSplit /= 2) then
      write(lerr,"(a)") "SCATTER> ERROR SEED expected 2 arguments:"
      write(lerr,"(a)") "SCATTER>       SEED seed"
      iErr = .true.
      return
    end if
    call str_cast(lnSplit(2), scatter_seed1, iErr)
    call recuinit(scatter_seed1)
    call recuut(scatter_seed1, scatter_seed2)

  case("ELEM")
    call scatter_parseElem(lnSplit, nSplit, iErr)

  case("PRO")
    call scatter_parseProfile(lnSplit, nSplit, iErr)

  case("GEN")
    call scatter_parseGenerator(lnSplit, nSplit, iErr)

  case default
    write(lerr,"(a)") "SCATTER> ERROR Keyword not recognised: '"//keyWord//"'"
    iErr = .true.
    return

  end select

end subroutine scatter_parseInputLine

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-10
! =================================================================================================
subroutine scatter_parseElem(lnSplit, nSplit, iErr)

  use crcoall
  use mod_alloc
  use mod_common
  use mod_settings
  use mod_common_main

  implicit none

  type(string), allocatable, intent(in)    :: lnSplit(:)
  integer,                   intent(in)    :: nSplit
  logical,                   intent(inout) :: iErr

  ! Temporary Variables
  type(scatter_elemStore), allocatable :: tmpElem(:)
  real(kind=fPrec),        allocatable :: brRatio(:)
  integer,                 allocatable :: generatorID(:)
  character(len=mNameLen)              :: bezName
  integer                              :: bezID, profileID
  real(kind=fPrec)                     :: elemScale, sigmaTot, ratioTot
  logical                              :: autoRatio

  integer i, j

  ! Check number of arguments
  if(nSplit < 5) then
    write(lerr,"(a)") "SCATTER> ERROR ELEM expected at least 5 arguments:"
    write(lerr,"(a)") "SCATTER>       ELEM elemname profile scaling gen1 (... genN)"
    iErr = .true.
    return
  end if

  if(allocated(scatter_elemList)) then
    allocate(tmpElem(scatter_nElem+1))
    tmpElem(1:scatter_nElem) = scatter_elemList(1:scatter_nElem)
    call move_alloc(tmpElem, scatter_elemList)
    scatter_nElem = scatter_nElem + 1
  else
    allocate(scatter_elemList(1))
    scatter_nElem = 1
  end if

  allocate(generatorID(nSplit-5))
  allocate(brRatio(nSplit-5))

  bezName        = trim(lnSplit(2))
  bezID          = -1
  elemScale      = zero
  sigmaTot       = zero
  ratioTot       = one
  autoRatio      = .false.
  profileID      = -1
  generatorID(:) = -1
  brRatio(:)     = zero

  ! Find the single element referenced
  do i=1,il
    if(bez(i) == lnSplit(2)) then
      if(scatter_elemPointer(i) /= 0) then
        write(lerr,"(a)") "SCATTER> ERROR Tried to define element '"//trim(lnSplit(2))//"' twice."
        iErr = .true.
        return
      end if

      if(kz(i) /= 40) then
        write(lerr,"(a)")    "SCATTER> ERROR SCATTER can only work on SINGLE ELEMENTs of type 40."
        write(lerr,"(a,i0)") "SCATTER>       The referenced element '"//trim(lnSplit(2))//"'is of type ", kz(i)
        iErr = .true.
        return
      end if

      if(el(i) /= 0 .or. ek(i) /= 0 .or. ed(i) /= 0) then
        write(lerr,"(a)") "SCATTER> ERROR Please check your input in the single element "//&
          "definition of your SCATTER. All values except for the type must be zero."
        iErr = .true.
        return
      end if

      bezID = i
      exit
    end if
  end do
  if(bezID == -1) then
    write(lerr,"(a)") "SCATTER> ERROR Could not find element '"//trim(lnSplit(2))//"'"
    iErr = .true.
    return
  end if

  ! Find the profile name referenced
  do i=1,scatter_nPro
    if(scatter_proList(i)%proName == lnSplit(3)) then
      profileID = i
      exit
    end if
  end do
  if(profileID == -1) then
    write(lerr,"(a)") "SCATTER> ERROR Could not find profile '"//trim(lnSplit(3))//"'"
    iErr = .true.
    return
  end if

  ! Store the ratio
  if(lnSplit(4)%lower() == "auto") then
    autoRatio = .true.
  else
    autoRatio = .false.
    call str_cast(lnSplit(4),ratioTot,iErr)
  end if

  ! Store the scaling
  call str_cast(lnSplit(5),elemScale,iErr)

  do i=6,nSplit
    ! Search for the generator with the right name
    do j=1, scatter_nGen
      if(scatter_genList(j)%genName == lnSplit(i)) then
        generatorID(i-5) = j
        brRatio(i-5)     = scatter_genList(j)%crossSection
        sigmaTot         = sigmaTot + scatter_genList(j)%crossSection
        exit
      end if
    end do

    ! If it is still -1, it wasn't found
    if(generatorID(i-5) == -1) then
      write(lerr,"(a)") "SCATTER> ERROR Parsing ELEM, generator '"//trim(lnSplit(i))//"' not found."
      iErr = .true.
      return
    end if

    ! Loop over those GENerators we've filled before
    ! (i.e. up to but not including column ii-4+2)
    ! to check for duplicates
    do j=1,i-6
      if(generatorID(i-5) == generatorID(j)) then
        write(lerr,"(a)") "SCATTER> ERROR Parsing ELEM, generator '"//trim(lnSplit(i))//"' used twice."
        iErr = .true.
        return
      end if
    end do
  end do

  do i=1,size(brRatio,1)
    brRatio(i) = brRatio(i) / sigmaTot
  end do

  scatter_elemPointer(bezID)                  = scatter_nElem
  scatter_elemList(scatter_nElem)%bezName     = bezName
  scatter_elemList(scatter_nElem)%bezID       = bezID
  scatter_elemList(scatter_nElem)%elemScale   = elemScale
  scatter_elemList(scatter_nElem)%sigmaTot    = sigmaTot
  scatter_elemList(scatter_nElem)%ratioTot    = ratioTot
  scatter_elemList(scatter_nElem)%autoRatio   = autoRatio
  scatter_elemList(scatter_nElem)%profileID   = profileID
  scatter_elemList(scatter_nElem)%generatorID = generatorID
  scatter_elemList(scatter_nElem)%brRatio     = brRatio

  if(scatter_debug .or. st_debug) then
    write(lout,"(a,i0,a)")    "SCATTER> DEBUG Element ",scatter_nElem,":"
    write(lout,"(a)")         "SCATTER> DEBUG  * bezName        = '"//trim(bezName)//"'"
    write(lout,"(a,i4)")      "SCATTER> DEBUG  * bezID          = ",bezID
    write(lout,"(a,f11.6)")   "SCATTER> DEBUG  * ratioTot       = ",ratioTot
    write(lout,"(a,3x,l1)")   "SCATTER> DEBUG  * autoRatio      = ",autoRatio
    write(lout,"(a,1pe15.6)") "SCATTER> DEBUG  * elemScale      = ",elemScale
    write(lout,"(a,f11.6)")   "SCATTER> DEBUG  * sigmaTot       = ",sigmaTot*c1e27
    write(lout,"(a,i4)")      "SCATTER> DEBUG  * profileID      = ",profileID
    do i=1,size(generatorID,1)
      write(lout,"(a,i0,a,i4)")    "SCATTER> DEBUG  * generatorID(",i,") = ",generatorID(i)
    end do
    do i=1,size(brRatio,1)
      write(lout,"(a,i0,a,f11.6)") "SCATTER> DEBUG  * brRatio(",i,")     = ",brRatio(i)
    end do
  end if

end subroutine scatter_parseElem

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-10
! =================================================================================================
subroutine scatter_parseProfile(lnSplit, nSplit, iErr)

  use crcoall
  use mod_alloc
  use mod_settings

  implicit none

  type(string), allocatable, intent(in)    :: lnSplit(:)
  integer,                   intent(in)    :: nSplit
  logical,                   intent(inout) :: iErr

  ! Temporary Variables
  type(scatter_proStore), allocatable :: tmpPro(:)
  character(len=:),       allocatable :: proName
  real(kind=fPrec),       allocatable :: fParams(:)
  integer proType
  integer i, ii, tmpIdx

  ! Check number of arguments
  if(nSplit < 3) then
    write(lerr,"(a)") "SCATTER> ERROR PRO expected at least 3 arguments:"
    write(lerr,"(a)") "SCATTER>       PRO name type (arguments...)"
    iErr = .true.
    return
  end if

  if(allocated(scatter_proList)) then
    allocate(tmpPro(scatter_nPro+1))
    tmpPro(1:scatter_nPro) = scatter_proList(1:scatter_nPro)
    call move_alloc(tmpPro, scatter_proList)
    scatter_nPro = scatter_nPro + 1
  else
    allocate(scatter_proList(1))
    scatter_nPro = 1
  end if

  allocate(fParams(nSplit-3))
  proName    = trim(lnSplit(2))
  proType    = -1
  fParams(:) = zero

  ! Check that the profile name is unique
  do i=1,scatter_nPro-1
    if(scatter_proList(i)%proName == lnSplit(2)) then
      write(lerr,"(a)") "SCATTER> ERROR Profile name '"//trim(lnSplit(2))//"' is not unique."
      iErr = .true.
      return
    end if
  end do

  ! Profile type dependent code
  select case (lnSplit(3)%get())
  case("FLAT")
    proType = 1 ! Integer code for FLAT
    if(nSplit /= 6) then
      write(lerr,"(a)") "SCATTER> ERROR PROfile type FLAT expected 6 arguments:"
      write(lerr,"(a)") "SCATTER>       PRO name FLAT density[targets/cm^2] mass[MeV/c^2] momentum[MeV/c]"
      iErr = .true.
      return
    end if

    call str_cast(lnSplit(4),fParams(1),iErr) ! Density
    call str_cast(lnSplit(5),fParams(2),iErr) ! Mass
    call str_cast(lnSplit(6),fParams(3),iErr) ! Momentum

  case("FIXED")
    proType = 2 ! Integer code for FIXED
    if(nSplit /= 4) then
      write(lerr,"(a)") "SCATTER> ERROR PROfile type FIXED expected 4 arguments:"
      write(lerr,"(a)") "SCATTER>       PRO name FIXED density[targets/m^2]"
      iErr = .true.
      return
    end if

    call str_cast(lnSplit(4),fParams(1),iErr) ! Density

  case("GAUSS1")
    proType = 10  ! Integer code for BEAM_GAUSS1
    if(nSplit /= 8) then
      write(lerr,"(a)") "SCATTER> ERROR PROfile type GAUSS1 expected 8 arguments:"
      write(lerr,"(a)") "SCATTER        PRO name GAUSS1 beamtot[particles] sigma_x[mm] sigma_y[mm] offset_x[mm] offset_y[mm]"
      iErr = .true.
      return
    end if

    call str_cast(lnSplit(4),fParams(1),iErr) ! Beam Charge
    call str_cast(lnSplit(5),fParams(2),iErr) ! Sigma X
    call str_cast(lnSplit(6),fParams(3),iErr) ! Sigma Y
    call str_cast(lnSplit(7),fParams(4),iErr) ! Offset X
    call str_cast(lnSplit(8),fParams(5),iErr) ! Offset Y

  case default
    write(lerr,"(a)") "SCATTER> ERROR PRO name '"//trim(lnSplit(3))//"' not recognized."
    iErr = .true.
    return

  end select

  scatter_proList(scatter_nPro)%proName = proName
  scatter_proList(scatter_nPro)%proType = proType
  scatter_proList(scatter_nPro)%fParams = fParams

  if(scatter_debug .or. st_debug) then
    write(lout,"(a,i0,a)")    "SCATTER> DEBUG Profile ",scatter_nPro,":"
    write(lout,"(a)")         "SCATTER> DEBUG  * proName    = '"//trim(proName)//"'"
    write(lout,"(a,i0)")      "SCATTER> DEBUG  * proType    = ",proType
    do i=1,size(fParams,1)
      write(lout,"(a,i0,a,e22.15)") "SCATTER> DEBUG  * fParams(",i,") = ",fParams(i)
    end do
  end if

end subroutine scatter_parseProfile

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-12
! =================================================================================================
subroutine scatter_parseGenerator(lnSplit, nSplit, iErr)

  use crcoall
  use mod_alloc
  use strings
  use string_tools
  use mod_settings

  implicit none

  type(string), allocatable, intent(in)    :: lnSplit(:)
  integer,                   intent(in)    :: nSplit
  logical,                   intent(inout) :: iErr

  ! Temporary Variables
  type(scatter_genStore), allocatable :: tmpGen(:)
  character(len=:),       allocatable :: genName
  real(kind=fPrec),       allocatable :: fParams(:)

  real(kind=fPrec) crossSection
  integer genType
  integer i, ii, tmpIdx

  ! Check number of arguments
  if(nSplit < 3) then
    write(lerr,"(a)") "SCATTER> ERROR GEN expected at least 3 arguments:"
    write(lerr,"(a)") "SCATTER>       GEN name type (arguments...)"
    iErr = .true.
    return
  end if

  if(allocated(scatter_genList)) then
    allocate(tmpGen(scatter_nGen+1))
    tmpGen(1:scatter_nGen) = scatter_genList(1:scatter_nGen)
    call move_alloc(tmpGen, scatter_genList)
    scatter_nGen = scatter_nGen + 1
  else
    allocate(scatter_genList(1))
    scatter_nGen = 1
  end if

  allocate(fParams(nSplit-3))
  genName      = trim(lnSplit(2))
  genType      = -1
  crossSection = zero
  fParams(:)   = zero

  ! Check that the generator name is unique
  do i=1,scatter_nGen-1
    if(scatter_genList(i)%genName == genName) then
      write(lerr,"(a)") "SCATTER> ERROR Generator name '"//trim(genName)//"' is not unique."
      iErr = .true.
      return
    end if
  end do

  ! Generator type-dependent code
  select case (lnSplit(3)%get())
  case("ABSORBER")

    genType = 1

  case("PPBEAMELASTIC")

    genType = 10
    if(nSplit /= 9) then
      write(lerr,"(a)") "SCATTER> ERROR GEN PPBEAMELASTIC expected 9 arguments:"
      write(lerr,"(a)") "SCATTER>       GEN name PPBEAMELASTIC a b1 b2 phi tmin crossSection"
      call prror
      iErr = .true.
      return
    end if

    call str_cast(lnSplit(4),fParams(1),iErr) ! a
    call str_cast(lnSplit(5),fParams(2),iErr) ! b1
    call str_cast(lnSplit(6),fParams(3),iErr) ! b2
    call str_cast(lnSplit(7),fParams(4),iErr) ! phi
    call str_cast(lnSplit(8),fParams(5),iErr) ! tmin
    call str_cast(lnSplit(9),fParams(6),iErr) ! crossSection
    crossSection = fParams(6) * c1m27         ! Set crossSection explicitly in mb

    ! Check sanity of input values
    if(fParams(2) <= zero) then
      write(lerr,"(a)") "SCATTER> ERROR GEN PPBEAMELASTIC 5th input (b1) must be larger than zero"
      iErr = .true.
      return
    end if
    if(fParams(3) <= zero) then
      write(lerr,"(a)") "SCATTER> ERROR GEN PPBEAMELASTIC 6th input (b2) must be larger than zero"
      iErr = .true.
      return
    end if

  case("PYTHIASIMPLE")

#ifndef PYTHIA
    write(lerr,"(a)") "SCATTER> ERROR GEN PYTHIA requested, but PYTHIA not compiled into SixTrack"
    iErr = .true.
    return
#endif
    if(nSplit /= 4) then
      write(lerr,"(a)") "SCATTER> ERROR GEN PYTHIASIMPLE expected 4 arguments:"
      write(lerr,"(a)") "SCATTER>       GEN name PYTHIASIMPLE crossSection"
      call prror
      iErr = .true.
      return
    end if

    genType = 20

    call str_cast(lnSplit(4),fParams(1),iErr) ! crossSection
    crossSection = fParams(1) * c1m27         ! Set crossSection explicitly in mb

  case default

    write(lerr,"(a)") "SCATTER> ERROR GEN name '"//trim(lnSplit(3))//"' not recognised."
    iErr = .true.
    return

  end select

  scatter_genList(scatter_nGen)%genName      = genName
  scatter_genList(scatter_nGen)%genType      = genType
  scatter_genList(scatter_nGen)%crossSection = crossSection
  scatter_genList(scatter_nGen)%fParams      = fParams

  if(scatter_debug .or. st_debug) then
    write(lout,"(a,i0,a)")    "SCATTER> DEBUG Generator ",scatter_nGen,":"
    write(lout,"(a)")         "SCATTER> DEBUG  * genName      = '"//trim(genName)//"'"
    write(lout,"(a,i0)")      "SCATTER> DEBUG  * genType      = ",genType
    write(lout,"(a,e22.15)")  "SCATTER> DEBUG  * crossSection = ",crossSection
    do i=1,size(fParams,1)
      write(lout,"(a,i0,a,e22.15)") "SCATTER> DEBUG  * fParams(",i,")   = ",fParams(i)
    end do
  end if

end subroutine scatter_parseGenerator

! =================================================================================================
! END Input Parser Functions
! =================================================================================================

! =================================================================================================
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-12
! =================================================================================================
subroutine scatter_setScaling(iElem, scaleVal)
  integer,          intent(in) :: iElem
  real(kind=fPrec), intent(in) :: scaleVal
  scatter_elemList(scatter_elemPointer(iElem))%elemScale = scaleVal
end subroutine scatter_setScaling

! =================================================================================================
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-12
! =================================================================================================
pure function scatter_getScaling(iElem) result(scaleVal)
  integer, intent(in) :: iElem
  real(kind=fPrec)    :: scaleVal
  scaleVal = scatter_elemList(scatter_elemPointer(iElem))%elemScale
end function scatter_getScaling

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 02-11-2017
! =================================================================================================
subroutine scatter_thin(iElem, ix, turn)

  use string_tools
  use crcoall
  use mod_time
  use mod_alloc
  use mod_common
  use mod_particles
  use mod_common_main
  use numerical_constants, only : pi
#ifdef HDF5
  use hdf5_output
#endif
  use collimation, only : do_coll, scatterhit, part_hit_pos, part_hit_turn

  implicit none

  integer, intent(in) :: iElem, ix, turn

  integer          idElem, idPro, iGen, nGen, idGen, iError
  integer          i, j, k
  integer          tmpSeed1, tmpSeed2
  logical          updateE, autoRatio, isDiff
  integer          iLost, procID
  real(kind=fPrec) t, dEE, dPP, theta, phi
  real(kind=fPrec) elemScale, sigmaTot, ratioTot, crossSection, scatterProb, targetDensity, scRatio, brRatio

  real(kind=fPrec), allocatable :: brThreshold(:)
  real(kind=fPrec), allocatable :: rndVals(:)
  logical,          allocatable :: pScattered(:)
  logical,          allocatable :: hasProc(:,:)
  integer,          allocatable :: nLost(:,:)
  integer,          allocatable :: nScattered(:,:)

#ifdef HDF5
  ! For HDF5 it is best to write in chuncks, so we will make arrays of size napx to cache everything
  integer,          allocatable :: iRecords(:,:)
  real(kind=fPrec), allocatable :: rRecords(:,:)
  character(len=:), allocatable :: cRecords(:,:)
  integer nRecords
#endif

  idElem    = scatter_elemPointer(ix)
  idPro     = scatter_elemList(idElem)%profileID
  elemScale = scatter_elemList(idElem)%elemScale
  sigmaTot  = scatter_elemList(idElem)%sigmaTot
  ratioTot  = scatter_elemList(idElem)%ratioTot
  autoRatio = scatter_elemList(idElem)%autoRatio
  nGen      = size(scatter_elemList(idElem)%generatorID,1)

  if(elemScale <= pieni) then
    ! Skip the whole thing if the scaling is ~zero
    return
  end if

  ! If not, we're doing something, so start the stop watch
  call time_startClock(time_clockSCAT)

  ! Store the seeds in the random number generator, and set ours
  call recuut(tmpSeed1,tmpSeed2)
  call recuin(scatter_seed1,scatter_seed2)

  call alloc(rndVals,   napx*3,             zero,   "rndVals")
  call alloc(pScattered,napx,               .false.,"pScattered")
  call alloc(hasProc,   nGen,scatter_nProc, .false.,"hasProc")
  call alloc(nLost,     nGen,scatter_nProc, 0,      "nLost")
  call alloc(nScattered,nGen,scatter_nProc, 0,      "nScattered")
#ifdef HDF5
  if(h5_useForSCAT) then
    call alloc(iRecords,         3,napx,    0,      "iRecords")
    call alloc(rRecords,         8,napx,    zero,   "rRecords")
    call alloc(cRecords,mNameLen,3,napx,    " ",    "cRecords")
  end if
  nRecords = 0
#endif

  t        = zero
  theta    = zero
  dEE      = zero
  dPP      = zero
  procID   = 0
  iLost    = 0
  isDiff   = .false.
  updateE  = .false.

  ! Compute Thresholds
  allocate(brThreshold(nGen))
  do i=1,nGen
    brThreshold(i) = sum(scatter_elemList(idElem)%brRatio(1:i))
    write(lout,"(a,i0,a,f13.6)") "SCATTER> Threshold ",i," = ",brThreshold(i)
  end do

  ! Generate random numbers for probability, branching ratio and phi angle
  call ranecu(rndVals, napx*3, -1)

  ! Loop over particles
  do j=1,napx

    k = 3*j-2 ! Indices in the random number array

    ! Compute Scattering Probability
    targetDensity = scatter_profile_getDensity(idPro,xv1(j),xv2(j))
    if(autoRatio) then
      scatterProb = (targetDensity*sigmaTot)*elemScale
    else
      scatterProb = ratioTot*elemScale
    end if
    if(rndVals(k) > scatterProb) then
      cycle
    else
      pScattered(j) = .true.
    end if

    ! Select Generator
    idGen = 0
    iGen  = 0
    do i=1,nGen
      if(rndVals(k+1) <= brThreshold(i)) then
        idGen = scatter_elemList(idElem)%generatorID(i)
        iGen  = i
        exit
      end if
    end do
    if(idGen == 0) then
      write(lout,"(a,i0,a)") "SCATTER> WARNING Scattering for particle #",j," occured, but no generator was selected."
      pScattered(j) = .false.
      cycle
    end if

    phi = (2*pi)*rndVals(k+2)

    ! If we're scaling the probability with DYNK, update the statistical weight
    scatter_statScale(partID(j)) = scatter_statScale(partID(j)) / elemScale

    ! Get event
    call scatter_generator_getEvent(idGen,j,t,theta,dEE,dPP,procID,iLost,isDiff)
    hasProc(iGen,procID)    = .true.
    nScattered(iGen,procID) = nScattered(iGen,procID) + 1

    if(scatter_allowLosses .and. iLost == 1) then
      ! If lost, flag it, count it, but no need to update energy and angle
      nLost(iGen,procID) = nLost(iGen,procID) + 1
      llostp(j) = .true.
      phi       = zero
    else
      ! Calculate new energy, and scattering angle
      if(abs(dPP) >= pieni) then ! Update the energy
        ejv(j)  = sqrt(((one+dPP)*ejfv(j))**2 + nucm(j)**2)
        updateE = .true.
      end if

      ! Update particle trajectory
      yv1(j) = theta*cos_mb(phi) + yv1(j)
      yv2(j) = theta*sin_mb(phi) + yv2(j)
    end if

    ! Output to file
#ifdef HDF5
    if(h5_useForSCAT) then
      nRecords = nRecords + 1
      iRecords(1,nRecords) = partID(j)
      iRecords(2,nRecords) = turn
      cRecords(1,nRecords) = bez(ix)
      cRecords(2,nRecords) = chr_rPadCut(scatter_genList(idGen)%genName,mNameLen)
      cRecords(3,nRecords) = scatter_procNames(procID)
      iRecords(3,nRecords) = iLost
      rRecords(1,nRecords) = t
      rRecords(2,nRecords) = theta
      rRecords(3,nRecords) = phi
      rRecords(4,nRecords) = dEE
      rRecords(5,nRecords) = dPP
      rRecords(6,nRecords) = targetDensity
      rRecords(7,nRecords) = scatterProb
      rRecords(8,nRecords) = scatter_statScale(partID(j))
    else
#endif
      write(scatter_logFile,"(2(1x,i8),2(1x,a20),1x,a8,1x,i4,1x,f12.3,1x,f12.6,1x,f9.6,5(1x,1pe16.9))") &
        partID(j), turn, chr_rPadCut(bez(ix),20), chr_rPadCut(scatter_genList(idGen)%genName,20), scatter_procNames(procID), &
        iLost, t, theta, phi, dEE, dPP, targetDensity, scatterProb, scatter_statScale(partID(j))
#ifdef CR
      scatter_logFilePos = scatter_logFilePos + 1
#endif
#ifdef HDF5
    end if
#endif

    if(do_coll) then
      scatterhit(j)    = 8
      part_hit_pos(j)  = iElem
      part_hit_turn(j) = turn
    endif

  end do

#ifdef HDF5
  if(h5_useForSCAT) then
    call h5_prepareWrite(scatter_logDataSet, nRecords)
    call h5_writeData(scatter_logDataSet, 1,  nRecords, iRecords(1,1:nRecords))
    call h5_writeData(scatter_logDataSet, 2,  nRecords, iRecords(2,1:nRecords))
    call h5_writeData(scatter_logDataSet, 3,  nRecords, cRecords(1,1:nRecords)(1:mNameLen))
    call h5_writeData(scatter_logDataSet, 4,  nRecords, cRecords(2,1:nRecords)(1:mNameLen))
    call h5_writeData(scatter_logDataSet, 5,  nRecords, cRecords(3,1:nRecords)(1:8))
    call h5_writeData(scatter_logDataSet, 6,  nRecords, iRecords(3,1:nRecords))
    call h5_writeData(scatter_logDataSet, 7,  nRecords, rRecords(1,1:nRecords))
    call h5_writeData(scatter_logDataSet, 8,  nRecords, rRecords(2,1:nRecords))
    call h5_writeData(scatter_logDataSet, 9,  nRecords, rRecords(3,1:nRecords))
    call h5_writeData(scatter_logDataSet, 10, nRecords, rRecords(4,1:nRecords))
    call h5_writeData(scatter_logDataSet, 11, nRecords, rRecords(5,1:nRecords))
    call h5_writeData(scatter_logDataSet, 12, nRecords, rRecords(6,1:nRecords))
    call h5_writeData(scatter_logDataSet, 13, nRecords, rRecords(7,1:nRecords))
    call h5_writeData(scatter_logDataSet, 14, nRecords, rRecords(8,1:nRecords))
    call h5_finaliseWrite(scatter_logDataSet)

    do i=1,nGen
      do k=1,scatter_nProc
        if(hasProc(i,k)) then
          idGen        = scatter_elemList(idElem)%generatorID(i)
          brRatio      = scatter_elemList(idElem)%brRatio(i)
          scRatio      = real(nScattered(i,k),fPrec)/napx
          crossSection = scatter_genList(idGen)%crossSection
          crossSection = crossSection*(scRatio/brRatio)
          call h5_prepareWrite(scatter_sumDataSet, 1)
          call h5_writeData(scatter_sumDataSet, 1,  1, turn)
          call h5_writeData(scatter_sumDataSet, 2,  1, napx)
          call h5_writeData(scatter_sumDataSet, 3,  1, bez(ix))
          call h5_writeData(scatter_sumDataSet, 4,  1, trim(scatter_genList(idGen)%genName))
          call h5_writeData(scatter_sumDataSet, 5,  1, scatter_procNames(k))
          call h5_writeData(scatter_sumDataSet, 6,  1, nScattered(i,k))
          call h5_writeData(scatter_sumDataSet, 7,  1, nLost(i,k))
          call h5_writeData(scatter_sumDataSet, 8,  1, scRatio)
          call h5_writeData(scatter_sumDataSet, 9,  1, crossSection*c1e27)
          call h5_writeData(scatter_sumDataSet, 10, 1, elemScale)
          call h5_finaliseWrite(scatter_sumDataSet)
        end if
      end do
    end do
  end if
#else
  flush(scatter_logFile)
#endif

  do i=1,nGen
    do k=1,scatter_nProc
      if(hasProc(i,k)) then
        idGen        = scatter_elemList(idElem)%generatorID(i)
        brRatio      = scatter_elemList(idElem)%brRatio(i)
        scRatio      = real(nScattered(i,k),fPrec)/napx
        crossSection = scatter_genList(idGen)%crossSection
        crossSection = crossSection*(scRatio/brRatio)
        write(scatter_sumFile,"(2(1x,i8),2(1x,a20),1x,a8,2(1x,i8),1x,f9.6,1x,f13.6,1x,1pe13.6)") &
          turn, napx, bez(ix)(1:20), chr_rPad(trim(scatter_genList(idGen)%genName),20), scatter_procNames(k), &
          nScattered(i,k), nLost(i,k), real(nScattered(i,k),fPrec)/napx, crossSection*c1e27, elemScale
#ifdef CR
        scatter_sumFilePos = scatter_sumFilePos + 1
#endif
      end if
    end do
  end do
  flush(scatter_sumFile)

  if(scatter_allowLosses) then
    call shuffleLostParticles
  end if

  if(updateE) then
    call part_updatePartEnergy(1,.false.)
  end if

  ! Restore seeds in random generator
  call recuut(scatter_seed1,scatter_seed2)
  call recuin(tmpSeed1,tmpSeed2)

  call dealloc(rndVals,   "rndVals")
  call dealloc(hasProc,   "hasProc")
  call dealloc(pScattered,"pScattered")
  call dealloc(nLost,     "nLost")
  call dealloc(nScattered,"nScattered")
#ifdef HDF5
  if(h5_useForSCAT) then
    call dealloc(iRecords,           "iRecords")
    call dealloc(rRecords,           "rRecords")
    call dealloc(cRecords, mNameLen, "cRecords")
  end if
#endif

  call time_stopClock(time_clockSCAT)

end subroutine scatter_thin

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 02-11-2017
! =================================================================================================
function scatter_profile_getDensity(idPro, x, y) result(retval)

  use string_tools
  use crcoall
  use mod_common
  use numerical_constants, only : pi, zero

  implicit none

  integer,          intent(in) :: idPro
  real(kind=fPrec), intent(in) :: x, y

  real(kind=fPrec) beamtot, sigmaX, sigmaY, offsetX, offsetY, retVal
  integer tmpIdx

  select case(scatter_proList(idPro)%proType)
  case (1)  ! FLAT
    retval  = scatter_proList(idPro)%fParams(1)

  case (2)  ! FIXED
    retval  = scatter_proList(idPro)%fParams(1)

  case (10) ! GAUSS1
    beamtot = scatter_proList(idPro)%fParams(1)
    sigmaX  = scatter_proList(idPro)%fParams(2)
    sigmaY  = scatter_proList(idPro)%fParams(3)
    offsetX = scatter_proList(idPro)%fParams(4)
    offsetY = scatter_proList(idPro)%fParams(5)
    retval  = ((beamtot/(two*(pi*(sigmaX*sigmaY)))) &
      * exp_mb(-half*((x-offsetX)/sigmaX)**2))      &
      * exp_mb(-half*((y-offsetY)/sigmaY)**2)

  case default
    write(lerr,"(a,i0,a)") "SCATTER> ERROR Type ",scatter_proList(idPro)%proType," for profile '"//&
      trim(scatter_proList(idPro)%proName)//"' not understood."
    call prror
  end select

end function scatter_profile_getDensity

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 09-2017
! =================================================================================================
subroutine scatter_profile_getParticle(idPro, x, y, xp, yp, E)

  implicit none

  integer,          intent(in)  :: idPro
  real(kind=fPrec), intent(in)  :: x, y
  real(kind=fPrec), intent(out) :: xp, yp, E

  ! Return a particle to collide with
  ! Add dummy variables for now, which stops ifort from complaining
  xp = zero
  yp = zero
  E  = zero

end subroutine scatter_profile_getParticle

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 09-2017
! =================================================================================================
function scatter_generator_getCrossSection(idPro, idGen, x, y, xp, yp, E) result(retVal)

  use string_tools
  use crcoall
  use numerical_constants, only : zero

  implicit none

  integer,          intent(in) :: idPro, idGen
  real(kind=fPrec), intent(in) :: x, y, xp, yp, E

  ! Temporary variables
  integer          tmpIdx
  real(kind=fPrec) xpTarget, ypTarget, ETarget, retVal

  ! Calculate S
  call scatter_profile_getParticle(idPro, x, y, xpTarget, ypTarget, ETarget)

  retVal = zero

  ! Calculate the cross section as function of S
  ! This is currently a fixed value
  select case(scatter_genList(idGen)%genType)
  case(1)  ! ABSORBER
    retVal = scatter_genList(idGen)%crossSection

  case(10) ! PPBEAMELASTIC
    retVal = scatter_genList(idGen)%crossSection

  case(20) ! PYTHIASIMPLE
    retVal = scatter_genList(idGen)%crossSection

  case default
    write(lerr,"(a,i0,a)") "SCATTER> ERROR Type ",scatter_genList(idGen)%genType," for generator '"//&
      trim(scatter_genList(idGen)%genName)//"' not understood."
    call prror

  end select

end function scatter_generator_getCrossSection

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2017-11-02
!  Updated: 2019-04-24
! =================================================================================================
subroutine scatter_generator_getEvent(genID, pID, t, theta, dEE, dPP, procID, iLost, isDiff)

  use crcoall
  use mod_common_main
  use mod_common, only : fort3
  use, intrinsic :: iso_c_binding

  implicit none

  integer,          intent(in)  :: genID
  integer,          intent(in)  :: pID
  real(kind=fPrec), intent(out) :: t
  real(kind=fPrec), intent(out) :: theta
  real(kind=fPrec), intent(out) :: dEE
  real(kind=fPrec), intent(out) :: dPP
  integer,          intent(out) :: procID
  integer,          intent(out) :: iLost
  logical,          intent(out) :: isDiff

  ! Temporary variables
  logical(kind=C_BOOL) evStat
  integer              tmpIdx, evType, nRetry
  real(kind=fPrec)     a, b1, b2, phi, tmin

  dEE    = zero
  dPP    = zero
  t      = zero
  theta  = zero
  iLost  = 0
  nRetry = 0
  isDiff = .false.

  select case(scatter_genList(genID)%genType)
  case(1)  ! ABSORBER

    procID = scatter_idAbsorb

  case(10) ! PPBEAMELASTIC

    a      = scatter_genList(genID)%fParams(1)
    b1     = scatter_genList(genID)%fParams(2)
    b2     = scatter_genList(genID)%fParams(3)
    phi    = scatter_genList(genID)%fParams(4)
    tmin   = scatter_genList(genID)%fParams(5)

    t      = scatter_generator_getPPElastic(a, b1, b2, phi, tmin)
    t      = t*c1e6                                      ! Scale return variable to MeV^2
    theta  = acos_mb(one - (t/(2*ejfv(pID)**2)))*c1e3 ! Get angle from t
    procID = scatter_idElastic

  case(20) ! PYTHIA

#ifdef PYTHIA
10  continue
    call pythia_getEvent(evStat, evType, t, theta, dEE, dPP)
    nRetry = nRetry + 1
    if(nRetry > 100) then
      write(lout,"(a)") "SCATTER> WARNING Pythia failed to generate event. Skipping Particle."
      procID = scatter_idError
      iLost  = 0
      return
    end if
    if(evStat) then
      t     = abs(t)*c1e6 ! Scale t to MeV^2
      theta = theta*c1e3  ! Scale angle to mrad
    ! theta = acos_mb(one - (t/((2*ejfv(j)**2)*(one+dPP))))*c1e3 ! Calculated from t
    select case(evType)
      case(pythia_idNonDiff)
        if(scatter_allowLosses) then
          procID = scatter_idNonDiff
          iLost  = 1
          isDiff = .false.
        else
          write(lerr,"(a)") "SCATTER> ERROR Particle lost, but losses not explicitly allowed in "//trim(fort3)
          call prror
        end if
      case(pythia_idElastic)
        procID = scatter_idElastic
        iLost  = 0
        isDiff = .false.
      case(pythia_idSingleDiffXB)
        if(scatter_allowLosses) then
          procID = scatter_idSingleDiffXB
          iLost  = 1
          isDiff = .true.
        else
          goto 10
        end if
      case(pythia_idSingleDiffAX)
        procID = scatter_idSingleDiffAX
        iLost  = 0
        isDiff = .true.
      case(pythia_idDoubleDiff)
        if(scatter_allowLosses) then
          procID = scatter_idDoubleDiff
          iLost  = 1
          isDiff = .true.
        else
          write(lerr,"(a)") "SCATTER> ERROR Particle lost, but losses not explicitly allowed in "//trim(fort3)
          call prror
        end if
      case(pythia_idCentralDiff)
        procID = scatter_idCentralDiff
        iLost  = 0
        isDiff = .true.
      case default
        procID = scatter_idUnknown
        iLost  = 0
        isDiff = .false.
      end select
    else
      write(lout,"(a)") "SCATTER> WARNING Pythia failed to generate event. Pythia error."
      goto 10
    end if
#else
    write(lerr,"(a,i0,a)") "SCATTER> ERROR This version of SixTrack was built without PYTHIA support,"
    call prror
#endif

  case default
    write(lerr,"(a,i0,a)") "SCATTER> ERROR Generator type ",scatter_genList(genID)%genType," not understood"
    call prror

  end select

end subroutine scatter_generator_getEvent

! =================================================================================================
!  H. Burkhardt, K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2017-11-02
!  Converted from C++ code from:
!  "Elastic pp scattering estimates and simulation relevant for burn-off"
!  https://indico.cern.ch/event/625576/
! =================================================================================================
function scatter_generator_getPPElastic(a, b1, b2, phi, tmin) result(t)

  use crcoall

  implicit none

  real(kind=fPrec), intent(in) :: a, b1, b2, phi, tmin

  ! Temp Variables
  integer          nItt, maxItt
  real(kind=fPrec) g1, g2, g3, gg, prob3, invB1, invB2, rndArr(3), t

  ! Approximate distribution
  g1    =           exp_mb(-b1*tmin)/b1  ! Soft scatter term
  g3    = (a**2) * (exp_mb(-b2*tmin)/b2) ! Hard scatter term
  prob3 = g3/(g1+g3)                     ! Probability of g3

  ! Pre-calculate inverses
  invB1 = one/b1
  invB2 = one/b2

  nItt   = 0
  maxItt = 1000000
  do
    nItt = nItt + 1
    call ranecu(rndArr, 3, -1)

    ! Randomly switch between g1 and g3 according to probability
    if(rndArr(1) > prob3) then
      t = tmin - invB1*log_mb(rndArr(2))
    else
      t = tmin - invB2*log_mb(rndArr(2))
    end if

    ! Exact distribution
    g1 =             exp_mb(-b1*t)                              ! Soft scatter term
    g2 =  ((two*a) * exp_mb(((-half)*(b1+b2))*t)) * cos_mb(phi) ! Interference
    g3 =    (a**2) * exp_mb(-b2*t)                              ! Hard scatter term
    gg = ((g1+g2)+g3)/(g1+g3)

    ! Check hit/miss on exact distribution
    ! if miss, get new t value
    ! if too many attempts, exit
    if(rndArr(3) < gg .or. nItt > maxItt) exit
  end do

  if(nItt > maxItt) then
    write(lerr,"(a)")      "SCATTER> ERROR in generator PPBEAMELASTIC"
    write(lerr,"(a,i0,a)") "SCATTER>       Limit of ",maxItt," misses in generator loop reached."
    call prror
  end if

end function scatter_generator_getPPElastic

! =================================================================================================
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-12
!  Write a report of all the calculated cross sections and branching ratios
! =================================================================================================
subroutine scatter_writeReport

  use crcoall

  implicit none

  integer i, iElem, iGen, nLine

  nLine = 0

  ! Calculate branching ratios and write summary
  write(scatter_sumFile,"(a)")      "#"
  write(scatter_sumFile,"(a)")      "#  SCATTER SUMMARY"
  write(scatter_sumFile,"(a)")      "# ================="
  write(scatter_sumFile,"(a)")      "#"
  write(scatter_sumFile,"(a,i0,a)") "#  Read ",scatter_nElem," element(s)"
  write(scatter_sumFile,"(a,i0,a)") "#  Read ",scatter_nPro," profile(s)"
  write(scatter_sumFile,"(a,i0,a)") "#  Read ",scatter_nGen," generator(s)"
  write(scatter_sumFile,"(a)")      "# "
  nLine = nLine + 8

  write(scatter_sumFile,"(a)")      "#  Generators"
  write(scatter_sumFile,"(a)")      "# ============="
  do iGen=1,scatter_nGen
    write(scatter_sumFile,"(a,i0,a)") "#  Generator(",iGen,"): '"//trim(scatter_genList(iGen)%genName)//"'"
    nLine = nLine + 1
  end do
  write(scatter_sumFile,"(a)")      "#"
  nLine = nLine + 3

  write(scatter_sumFile,"(a)")      "#  Cross Sections"
  write(scatter_sumFile,"(a)")      "# ================"
  do iElem=1,scatter_nElem
    if(scatter_elemList(iElem)%autoRatio) then
      write(scatter_sumFile,"(a,i0,a)") "#  Element(",iElem,"): '"//trim(scatter_elemList(iElem)%bezName)//"' "//&
        "[Probability: Auto]"
    else
      write(scatter_sumFile,"(a,i0,a,f8.6,a)") "#  Element(",iElem,"): '"//trim(scatter_elemList(iElem)%bezName)//"' "//&
        "[Probability: ",scatter_elemList(iElem)%ratioTot,"]"
    end if
    do i=1,size(scatter_elemList(iElem)%generatorID)
      iGen = scatter_elemList(iElem)%generatorID(i)
      write(scatter_sumFile,"(a,i0,a,f15.6,a,f8.6,a)") "#   + Generator(",iGen,"): ",&
        (scatter_genList(iGen)%crossSection*c1e27)," mb [BR: ",scatter_elemList(iElem)%brRatio(i),"]"
      nLine = nLine + 1
    end do
    write(scatter_sumFile,"(a,f15.6,a,f8.6,a)") "#   = Sigma Total:  ",(scatter_elemList(iElem)%sigmaTot*c1e27)," mb [BR: ",one,"]"
    nLine = nLine + 2
  end do
  write(scatter_sumFile,"(a)")      "#"
  nLine = nLine + 3

#ifdef CR
  scatter_sumFilePos = nLine
#endif

end subroutine scatter_writeReport

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-05-23
!  Called from CRCHECK; reads the _CR arrays for scatter from file
!  Sets readErr to true if something goes wrong while reading.
! =================================================================================================
#ifdef CR
subroutine scatter_crcheck_readdata(fileUnit, readErr)

  use parpro
  use crcoall
  use mod_alloc
  use numerical_constants

  implicit none

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: readErr

  integer j

  call alloc(scatter_statScale_CR, npart, zero, "scatter_statScale_CR")

  read(fileUnit, err=10, end=10) scatter_logFilePos_CR, scatter_sumFilePos_CR
  read(fileUnit, err=10, end=10) scatter_seed1_CR, scatter_seed2_CR
  read(fileUnit, err=10, end=10) scatter_statScale_CR

  readErr = .false.
  return

10 continue
  readErr = .true.
  write(lout,"(a,i0,a)") "SIXTRACR> ERROR Reading C/R file fort.",fileUnit," in SCATTER"
  write(crlog,  "(a,i0,a)") "SIXTRACR> ERROR Reading C/R file fort.",fileUnit," in SCATTER"
  flush(crlog)

end subroutine scatter_crcheck_readdata

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-10
!  Called from CRCHECK; resets the position of scatter_log.dat
! =================================================================================================
subroutine scatter_crcheck_positionFiles

  use parpro
  use crcoall
  use mod_units

  logical isOpen, fErr
  integer iError
  integer j
  character(len=mInputLn) aRecord

  call f_requestUnit("scatter_log.dat",scatter_logFile)
  inquire(unit=scatter_logFile, opened=isOpen)
  if(isOpen) then
    write(crlog,"(a)")      "CR_CHECK> ERROR Failed while repsositioning 'scatter_log.dat'"
    write(crlog,"(a,i0,a)") "CR_CHECK>       unit ",scatter_logFile," already in use!"
    flush(crlog)
    write(lerr,"(a)") "CR_CHECK> ERROR Failed positioning 'scatter_log.dat'"
    call prror
  end if

  if(scatter_logFilePos_CR /= -1) then
    call f_open(unit=scatter_logFile,file="scatter_log.dat",formatted=.true.,mode="rw",err=fErr,status="old")
    if(fErr) goto 10
    scatter_logFilePos = 0
    do j=1, scatter_logFilePos_CR
      read(scatter_logFile,"(a1024)",end=10,err=10,iostat=iError) aRecord
      scatter_logFilePos = scatter_logFilePos + 1
    end do
    endfile(scatter_logFile,iostat=iError)
    call f_close(scatter_logFile)

    call f_open(unit=scatter_logFile,file="scatter_log.dat",formatted=.true.,mode="w+",err=fErr,status="old")
    if(fErr) goto 10
    write(97,"(2(a,i0))") "CR_CHECK> Sucessfully repositioned 'scatter_log.dat': "//&
      "Position: ",scatter_logFilePos,", Position C/R: ",scatter_logFilePos_CR
    flush(crlog)

  else
    write(crlog,"(a,i0)") "CR_CHECK> Did not attempt repositioning 'scatter_log.dat' at line ",scatter_logFilePos_CR
    write(crlog,"(a)")    "CR_CHECK> If anything has been written to the file, it will be correctly truncated in Scatter INIT."
    flush(crlog)
  end if

  call f_requestUnit("scatter_summary.dat",scatter_sumFile)
  inquire(unit=scatter_sumFile, opened=isOpen)
  if(isOpen) then
    write(crlog,"(a)")      "CR_CHECK> ERROR Failed while repsositioning 'scatter_summary.dat'"
    write(crlog,"(a,i0,a)") "CR_CHECK>       Unit ",scatter_sumFile," already in use!"
    flush(crlog)
    write(lout,"(a)") "CR_CHECK> Failed positioning 'scatter_summary.dat'"
    call prror
  end if

  if(scatter_sumFilePos_CR /= -1) then
    call f_open(unit=scatter_sumFile,file="scatter_summary.dat",formatted=.true.,mode="rw",err=fErr,status="old")
    if(fErr) goto 10
    scatter_sumFilePos = 0
    do j=1, scatter_sumFilePos_CR
      read(scatter_sumFile,"(a1024)",end=10,err=10,iostat=iError) aRecord
      scatter_sumFilePos = scatter_sumFilePos + 1
    end do
    endfile(scatter_sumFile,iostat=iError)
    close(scatter_sumFile)

    call f_open(unit=scatter_sumFile,file="scatter_summary.dat",formatted=.true.,mode="w+",err=fErr,status="old")
    if(fErr) goto 10
    write(97,"(2(a,i0))") "CR_CHECK> Sucessfully repositioned 'scatter_summary.dat': "//&
      "Position: ",scatter_sumFilePos,", Position C/R: ",scatter_sumFilePos_CR
    flush(crlog)

  else
    write(crlog,"(a,i0)") "CR_CHECK> Did not attempt repositioning 'scatter_summary.dat' at line ",scatter_sumFilePos_CR
    write(crlog,"(a)")    "CR_CHECK> If anything has been written to the file, it will be correctly truncated in Scatter INIT."
    flush(crlog)
  end if

  return

10 continue
  write(crlog,"(a,i0)")    "CR_CHECK> ERROR While reading 'scatter_log.dat' or 'scatter_summary.dat', iostat = ",iError
  write(crlog,"(2(a,i0))") "CR_CHECK>       Position: ",scatter_logFilePos,", Position C/R: ",scatter_logFilePos_CR
  write(crlog,"(2(a,i0))") "CR_CHECK>       Position: ",scatter_sumFilePos,", Position C/R: ",scatter_sumFilePos_CR
  flush(crlog)
  write(lerr,"(a)")"CR_CHECK> ERROR CRCHECK failure positioning 'scatter_log.dat' or 'scatter_summary.dat'."
  call prror

end subroutine scatter_crcheck_positionFiles

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-05-23
!  Called from CRPOINT; write checkpoint data to fort.95/96
! =================================================================================================
subroutine scatter_crpoint(fileUnit, writeErr)

  use crcoall

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: writeErr

  writeErr = .false.

  write(fileunit,err=10) scatter_logFilePos, scatter_sumFilePos
  write(fileunit,err=10) scatter_seed1, scatter_seed2
  write(fileunit,err=10) scatter_statScale
  flush(fileUnit)

  return

10 continue
  writeErr = .true.
  write(lout, "(a,i0,a)") "CR_POINT> ERROR Writing C/R file fort.",fileUnit," in SCATTER"
  write(crlog,"(a,i0,a)") "CR_POINT> ERROR Writing C/R file fort.",fileUnit," in SCATTER"
  flush(crlog)

end subroutine scatter_crpoint

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-05-23
!  Called from CRSTART
! =================================================================================================
subroutine scatter_crstart

  use mod_alloc

  scatter_seed1 = scatter_seed1_CR
  scatter_seed2 = scatter_seed2_CR
  scatter_statScale(:) = scatter_statScale_CR(:)

  call dealloc(scatter_statScale_CR,"scatter_statScale_CR")

end subroutine scatter_crstart
#endif
! End of CR

end module scatter
