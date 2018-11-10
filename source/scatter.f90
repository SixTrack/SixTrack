! ================================================================================================ !
!  SixTrack SCATTER Module
!  V.K. Berglyd Olsen, K.N. Sjobak, H. Burkhardt, BE-ABP-HSS
!  Last modified: 2018-06-18
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
  use numerical_constants, only : zero, half, one, two, c1e3, c1e6, c1m27, pieni
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
  integer,          parameter :: scatter_idAbsorb         = 1
  integer,          parameter :: scatter_idNonDiff        = 2
  integer,          parameter :: scatter_idElastic        = 3
  integer,          parameter :: scatter_idSingleDiffXB   = 4
  integer,          parameter :: scatter_idSingleDiffAX   = 5
  integer,          parameter :: scatter_idDoubleDiff     = 6
  integer,          parameter :: scatter_idCentralDiff    = 7
  integer,          parameter :: scatter_idUnknown        = 8
  integer,          parameter :: scatter_idError          = 9
  character(len=8), parameter :: scatter_procNames(9)     = &
    (/"Absorbed","NonDiff ","Elastic ","SingD_XB","SingD_AX","DoubD_XX","CentDiff","Unknown ","Error   "/)

  ! Total cross section
  real(kind=fPrec), private, save :: scatter_sigmaTot     = one

  ! Pointer from an element back to an ELEM statement (0 => not used)
  integer,          allocatable, public, save :: scatter_elemPointer(:)

  ! Statistical correction factor for a specific particle
  real(kind=fPrec), allocatable, public, save :: scatter_statScale(:)

  ! Configuration for an ELEM, columns are:
  ! (1)   : pointer to the SingleElement
  ! (2)   : pointer to PROFILE
  ! (3-5) : pointer to GENERATORs
  integer,          allocatable, public, save :: scatter_ELEM(:,:)
  real(kind=fPrec), allocatable, public, save :: scatter_ELEM_scale(:)

  ! Configuration for PROFILE and GENERATOR
  ! Columns of scatter_PROFILE:
  ! (1)   : Profile name in fort.3 (points within scatter_cData)
  ! (2)   : Profile type
  ! (3-5) : Arguments (often pointing within scatter_{i|c|f}Data)
  integer, allocatable, public, save :: scatter_PROFILE(:,:)
  integer, allocatable, public, save :: scatter_GENERATOR(:,:)

  integer,          allocatable, private, save :: scatter_iData(:)
  real(kind=fPrec), allocatable, private, save :: scatter_fData(:)
  character(len=:), allocatable, private, save :: scatter_cData(:)

  ! Number of currently used positions in arrays
  integer, public,  save :: scatter_nELEM      = 0
  integer, public,  save :: scatter_nPROFILE   = 0
  integer, public,  save :: scatter_nGENERATOR = 0
  integer, private, save :: scatter_niData     = 0
  integer, private, save :: scatter_nfData     = 0
  integer, private, save :: scatter_ncData     = 0

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
#endif

contains

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 02-11-2017
! =================================================================================================
subroutine scatter_allocate

  use mod_alloc
  implicit none

  call alloc(scatter_ELEM,              1,5, 0,          "scatter_ELEM")
  call alloc(scatter_ELEM_scale,        1,   zero,       "scatter_ELEM_scale")
  call alloc(scatter_PROFILE,           1,5, 0,          "scatter_PROFILE")
  call alloc(scatter_GENERATOR,         1,5, 0,          "scatter_GENERATOR")

  call alloc(scatter_iData,             1,   0,          "scatter_iData")
  call alloc(scatter_fData,             1,   zero,       "scatter_fData")
  call alloc(scatter_cData,    mStrLen, 1,   str_dSpace, "scatter_cData")

end subroutine scatter_allocate

! =================================================================================================
!  Used for changing the allocation of arrays that scale with global parameters like NELE
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
!  Last Modified: 2018-04-16
! =================================================================================================
subroutine scatter_initialise

  use crcoall
  use parpro
  use file_units
  use mod_units

  implicit none

  integer iError
  logical fErr

#ifdef HDF5

  type(h5_dataField), allocatable :: setFields(:)

  if(h5_useForSCAT) then

    call h5_initForScatter()

    ! Scatter Log
    allocate(setFields(14))

    setFields(1)  = h5_dataField(name="ID",      type=h5_typeInt)
    setFields(2)  = h5_dataField(name="TURN",    type=h5_typeInt)
    setFields(3)  = h5_dataField(name="BEZ",     type=h5_typeChar, size=mNameLen)
    setFields(4)  = h5_dataField(name="GEN",     type=h5_typeChar, size=mNameLen)
    setFields(5)  = h5_dataField(name="PROCESS", type=h5_typeChar, size=8)
    setFields(6)  = h5_dataField(name="LOST",    type=h5_typeInt)
    setFields(7)  = h5_dataField(name="T",       type=h5_typeReal)
    setFields(8)  = h5_dataField(name="dE/E",    type=h5_typeReal)
    setFields(9)  = h5_dataField(name="dP/P",    type=h5_typeReal)
    setFields(10) = h5_dataField(name="THETA",   type=h5_typeReal)
    setFields(11) = h5_dataField(name="PHI",     type=h5_typeReal)
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
    allocate(setFields(8))

    setFields(1) = h5_dataField(name="TURN",     type=h5_typeInt)
    setFields(2) = h5_dataField(name="BEZ",      type=h5_typeChar, size=mNameLen)
    setFields(3) = h5_dataField(name="GEN",      type=h5_typeChar, size=mNameLen)
    setFields(4) = h5_dataField(name="PROCESS",  type=h5_typeChar, size=8)
    setFields(5) = h5_dataField(name="NSCATT",   type=h5_typeInt)
    setFields(6) = h5_dataField(name="NLOST",    type=h5_typeInt)
    setFields(7) = h5_dataField(name="CROSSSEC", type=h5_typeReal)
    setFields(8) = h5_dataField(name="SCALING",  type=h5_typeReal)

    call h5_createFormat("scatter_summary_fmt", setFields, scatter_sumFormat)
    call h5_createDataSet("summary", h5_scatID, scatter_sumFormat, scatter_sumDataSet)
    block
      character(len=:), allocatable :: colNames(:)
      character(len=:), allocatable :: colUnits(:)
      logical spErr
      integer nSplit
      call chr_split("turn bez generator process nScatt nLost crossSec scaling",colNames,nSplit,spErr)
      call chr_split("1 text text text 1 1 mb 1",colUnits,nSplit,spErr)
      call h5_writeDataSetAttr(scatter_sumDataSet,"colNames",colNames)
      call h5_writeDataSetAttr(scatter_sumDataSet,"colUnits",colUnits)
    end block

    deallocate(setFields)

    ! Write Attributes
    call h5_writeAttr(h5_scatID,"SEED",scatter_seed1)

    return ! No need to open files
  end if
#endif

  ! Open scatter_log.dat
  if(scatter_logFile == -1) call funit_requestUnit("scatter_log.dat", scatter_logFile)
#ifdef CR
  if(scatter_logFilePos == -1) then
    write(93,"(a)") "SCATTER> scatter_initialise opening new file 'scatter_log.dat'"
#endif
    call units_openUnit(unit=scatter_logFile,fileName="scatter_log.dat",formatted=.true.,mode="w",err=fErr,status="replace")
    write(scatter_logFile,"(a)") "# scatter_log"
    write(scatter_logFile,"(a1,a8,1x,a8,2(1x,a20),1x,a8,1x,a4,1x,a13,7(1x,a16))") &
      "#","ID","turn",chr_rPad("bez",20),chr_rPad("generator",20),chr_rPad("process",8),&
      "lost","t[MeV^2]","dE/E","dP/P","theta[mrad]","phi[rad]","density","prob","stat_corr"
    flush(scatter_logFile)
#ifdef CR
    scatter_logFilePos = 2
    endfile(scatter_logFile,iostat=iError)
    backspace(scatter_logFile,iostat=iError)
  else
    write(93,"(a)") "SCATTER> scatter_initialise kept already opened file 'scatter_log.dat'"
  end if
#endif

  ! Open scatter_summary.dat
  if(scatter_sumFile == -1) call funit_requestUnit("scatter_summary.dat",scatter_sumFile)
#ifdef CR
  if(scatter_sumFilePos == -1) then
    write(93,"(a)") "SCATTER> scatter_initialise opening new file 'scatter_summary.dat'"
#endif
    call units_openUnit(unit=scatter_sumFile,fileName="scatter_summary.dat",formatted=.true.,mode="w",err=fErr,status="replace")
    write(scatter_sumFile,"(a)") "# scatter_summary"
    write(scatter_sumFile,"(a1,a8,2(1x,a20),1x,a8,2(1x,a8),2(1x,a13))") &
      "#","turn",chr_rPad("bez",20),chr_rPad("generator",20),chr_rPad("process",8), &
      "nScatt","nLost","crossSec[mb]","scaling"
    flush(scatter_sumFile)
#ifdef CR
    scatter_sumFilePos = 2
    endfile(scatter_sumFile,iostat=iError)
    backspace(scatter_sumFile,iostat=iError)
  else
    write(93,"(a)") "SCATTER> scatter_initialise kept already opened file 'scatter_summary.dat'"
  end if
#endif

end subroutine scatter_initialise

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
    write(lout,"(a)") "SCATTER> ERROR Failed to parse input line."
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

  if(scatter_debug) then
    write (lout,"(a,i3,a)")  "SCATTER> DEBUG Input line len=",len(inLine),": '"//trim(inLine)//"'"
    write (lout,"(a,i2,a)") ("SCATTER> DEBUG  * Field(",i,") = '"//lnSplit(i)//"'",i=1,nSplit)
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
      write(lout,"(a)") "SCATTER> ERROR SEED expected 1 arguments:"
      write(lout,"(a)") "SCATTER>       SEED seed"
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
    write(lout,"(a)") "SCATTER> ERROR Keyword not recognised: '"//keyWord//"'"
    iErr = .true.
    return

  end select

end subroutine scatter_parseInputLine

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-20
! =================================================================================================
subroutine scatter_parseElem(lnSplit, nSplit, iErr)

  use crcoall
  use mod_alloc
  use mod_common
  use mod_commonmn

  implicit none

  type(string), allocatable, intent(in)    :: lnSplit(:)
  integer,                   intent(in)    :: nSplit
  logical,                   intent(inout) :: iErr

  ! Temporary Variables
  integer ii, j

  ! Check number of arguments
  if(nSplit < 5) then
    write(lout,"(a)") "SCATTER> ERROR ELEM expected at least 5 arguments:"
    write(lout,"(a)") "SCATTER>       ELEM elemname profile scaling gen1 (gen2 gen3)"
    iErr = .true.
    return
  end if

  ! Add the element to the list
  scatter_nELEM = scatter_nELEM + 1
  call alloc(scatter_ELEM,       scatter_nELEM, 5, 0, "scatter_ELEM")
  call alloc(scatter_ELEM_scale, scatter_nELEM, zero, "scatter_ELEM_scale")

  ! Find the single element referenced
  ii = -1
  do j=1,il
    if(bez(j) == lnSplit(2)) then
      if(ii /= -1) then
        write(lout,"(a)") "SCATTER> ERROR Found element '"//lnSplit(2)//"' twice in SINGLE ELEMENTS list."
        iErr = .true.
        return
      end if
      ii = j

      if(scatter_elemPointer(j) /= 0) then
        write(lout,"(a)") "SCATTER> ERROR Tried to define element '"//lnSplit(2)//"' twice."
        iErr = .true.
        return
      end if

      if(kz(j) /= 40) then
        write(lout,"(a)")    "SCATTER> ERROR SCATTER can only work on SINGLE ELEMENTs of type 40."
        write(lout,"(a,i0)") "SCATTER>       The referenced element '"//lnSplit(2)//"'is of type ", kz(j)
        iErr = .true.
        return
      end if

      if(el(j) /= 0 .or. ek(j) /= 0 .or. ed(j) /= 0) then
        write(lout,"(6(a,i0))") "SCATTER> ERROR Length el(j) (SCATTER is treated as thin element), "//&
          " and first and second field have to be zero: el(j)=ed(j)=ek(j)=0; "//&
          "but el(",j,")=",el(j),", ed(",j,")=",ed(j),", ek(",j,")=",ek(j),"."
        write(lout,"(a)") "SCATTER>       Please check your input in the single element "//&
          "definition of your SCATTER. All values except for the type must be zero."
        iErr = .true.
        return
      end if

      scatter_elemPointer(j) = scatter_nELEM
      scatter_ELEM(scatter_nELEM,1) = j
    end if
  end do

  if(scatter_ELEM(scatter_nELEM,1) == 0) then
    write(lout,"(a)") "SCATTER> ERROR Could not find element '"//lnSplit(2)//"'"
    iErr = .true.
    return
  end if

  ! Find the profile name referenced
  do j=1,scatter_nPROFILE
    if(trim(scatter_cData(scatter_PROFILE(j,1))) == lnSplit(3)) then
      scatter_ELEM(scatter_nELEM,2) = j
    end if
  end do

  if(scatter_ELEM(scatter_nELEM,2) == 0) then
    write(lout,"(a)") "SCATTER> ERROR Could not find profile '"//lnSplit(3)//"'"
    iErr = .true.
    return
  end if

  ! Store the scaling
  call str_cast(lnSplit(4),scatter_ELEM_scale(scatter_nELEM),iErr)

  ! Find the generator(s) referenced
  if(nSplit-4 > 3) then
    write(lout,"(a,i0,a)") "SCATTER> ERROR Parsing ELEM, ",nSplit-4," generators specified, max is 3"
    iErr = .true.
    return
  end if

  do ii=5,nSplit
    ! In case we won't find the generator name
    scatter_ELEM(scatter_nELEM,ii-4+2) = -1

    ! Search for the generator with the right name
    do j=1, scatter_nGENERATOR
      if(trim(scatter_cData(scatter_GENERATOR(j,1))) == lnSplit(ii)) then
        ! Found it
        scatter_ELEM(scatter_nELEM,ii-4+2) = j
      end if
    end do

    ! If it is still -1, it wasn't found
    if(scatter_ELEM(scatter_nELEM,ii-4+2) == -1) then
      write(lout,"(a)") "SCATTER> ERROR Parsing ELEM, generator '"//lnSplit(ii)//"' not found."
      iErr = .true.
      return
    end if

    ! Loop over those GENerators we've filled before
    ! (i.e. up to but not including column ii-4+2)
    ! to check for duplicates
    do j=3, ii-4+2-1
      if(scatter_ELEM(scatter_nELEM,j) == scatter_ELEM(scatter_nELEM,ii-4+2)) then
        write(lout,"(a)") "SCATTER> ERROR Parsing ELEM, generator '"//lnSplit(ii)//"' used twice."
        iErr = .true.
        return
      end if
    end do
  end do

end subroutine scatter_parseElem

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-20
! =================================================================================================
subroutine scatter_parseProfile(lnSplit, nSplit, iErr)

  use mod_alloc
  use crcoall

  implicit none

  type(string), allocatable, intent(in)    :: lnSplit(:)
  integer,                   intent(in)    :: nSplit
  logical,                   intent(inout) :: iErr

  ! Temporary Variables
  integer ii, tmpIdx

  ! Check number of arguments
  if(nSplit < 3) then
    write(lout,"(a)") "SCATTER> ERROR PRO expected at least 3 arguments:"
    write(lout,"(a)") "SCATTER>       PRO name type (arguments...)"
    iErr = .true.
    return
  end if

  ! Add a profile to the list
  scatter_nPROFILE = scatter_nPROFILE + 1
  call alloc(scatter_PROFILE, scatter_nPROFILE, 5, 0, "scatter_PROFILE")

  ! Store the profile name
  scatter_ncData = scatter_ncData + 1
  call alloc(scatter_cData, mStrLen, scatter_ncData, str_dSpace, "scatter_cData")

  scatter_cData(scatter_ncData)       = lnSplit(2)
  scatter_PROFILE(scatter_nPROFILE,1) = scatter_ncData

  ! Check that the profile name is unique
  do ii=1,scatter_nPROFILE-1
    if(trim(scatter_cData(scatter_PROFILE(ii,1))) == lnSplit(2)) then
      write(lout,"(a)") "SCATTER> ERROR Profile name '"//lnSplit(2)//"' is not unique."
      iErr = .true.
      return
    end if
  end do

  ! Profile type dependent code
  select case (lnSplit(3)%get())
  case("FLAT")
    scatter_PROFILE(scatter_nPROFILE,2) = 1 ! Integer code for FLAT
    if(nSplit /= 6) then
      write(lout,"(a)") "SCATTER> ERROR PROfile type FLAT expected 6 arguments:"
      write(lout,"(a)") "SCATTER>       PRO name FLAT density[targets/cm^2] mass[MeV/c^2] momentum[MeV/c]"
      iErr = .true.
      return
    end if

    ! Request space to store the density
    tmpIdx = scatter_nfData + 1
    scatter_PROFILE(scatter_nPROFILE,3) = tmpIdx
    scatter_nfData = scatter_nfData + 3
    call alloc(scatter_fData, scatter_nfData, zero, "scatter_fData")

    call str_cast(lnSplit(4),scatter_fData(tmpIdx),  iErr) ! Density
    call str_cast(lnSplit(5),scatter_fData(tmpIdx+1),iErr) ! Mass
    call str_cast(lnSplit(6),scatter_fData(tmpIdx+2),iErr) ! Momentum

  case("FIXED")
    scatter_PROFILE(scatter_nPROFILE,2) = 2 ! Integer code for FIXED
    if(nSplit /= 4) then
      write(lout,"(a)") "SCATTER> ERROR PROfile type FIXED expected 4 arguments:"
      write(lout,"(a)") "SCATTER>       PRO name FIXED density[targets/m^2]"
      iErr = .true.
      return
    end if

    ! Request space to store the density
    tmpIdx = scatter_nfData + 1
    scatter_PROFILE(scatter_nPROFILE,3) = tmpIdx
    scatter_nfData = scatter_nfData + 1
    call alloc(scatter_fData, scatter_nfData, zero, "scatter_fData")

    call str_cast(lnSplit(4),scatter_fData(tmpIdx),iErr) ! Density

  case("GAUSS1")
    scatter_PROFILE(scatter_nPROFILE,2) = 10  ! Integer code for BEAM_GAUSS1
    if(nSplit /= 8) then
      write(lout,"(a)") "SCATTER> ERROR PROfile type GAUSS1 expected 8 arguments:"
      write(lout,"(a)") "SCATTER        PRO name GAUSS1 beamtot[particles] sigma_x[mm] sigma_y[mm] offset_x[mm] offset_y[mm]"
      iErr = .true.
      return
    end if

    ! Request space to store the density
    tmpIdx = scatter_nfData + 1
    scatter_PROFILE(scatter_nPROFILE,3) = tmpIdx
    scatter_nfData = scatter_nfData + 5
    call alloc(scatter_fData, scatter_nfData, zero, "scatter_fData")

    call str_cast(lnSplit(4),scatter_fData(tmpIdx),iErr)   ! Beam Charge
    call str_cast(lnSplit(5),scatter_fData(tmpIdx+1),iErr) ! Sigma X
    call str_cast(lnSplit(6),scatter_fData(tmpIdx+2),iErr) ! Sigma Y
    call str_cast(lnSplit(7),scatter_fData(tmpIdx+3),iErr) ! Offset X
    call str_cast(lnSplit(8),scatter_fData(tmpIdx+4),iErr) ! Offset Y

  case default
    write(lout,"(a)") "SCATTER> ERROR PRO name '"//lnSplit(3)//"' not recognized."
    iErr = .true.
    return

  end select

end subroutine scatter_parseProfile

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 09-2017
! =================================================================================================
subroutine scatter_parseGenerator(lnSplit, nSplit, iErr)

  use crcoall
  use mod_alloc
  use strings
  use string_tools

  implicit none

  type(string), allocatable, intent(in)    :: lnSplit(:)
  integer,                   intent(in)    :: nSplit
  logical,                   intent(inout) :: iErr

  ! Temporary Variables
  integer ii, tmpIdx

  ! Check number of arguments
  if(nSplit < 3) then
    write(lout,"(a)") "SCATTER> ERROR GEN expected at least 3 arguments:"
    write(lout,"(a)") "SCATTER>       GEN name type (arguments...)"
    iErr = .true.
    return
  end if

  ! Add a generator to the list
  scatter_nGENERATOR = scatter_nGENERATOR + 1
  call alloc(scatter_GENERATOR, scatter_nGENERATOR, 5, 0, "scatter_GENERATOR")

  ! Store the generator name
  scatter_ncData = scatter_ncData + 1
  call alloc(scatter_cData, mStrLen, scatter_ncData, str_dSpace, "scatter_cData")

  scatter_cData(scatter_ncData)           = lnSplit(2)
  scatter_GENERATOR(scatter_nGENERATOR,1) = scatter_ncData

  ! Check that the generator name is unique
  do ii=1,scatter_nGENERATOR-1
    if(trim(scatter_cData(scatter_GENERATOR(ii,1))) == lnSplit(2)) then
      write(lout,"(a)") "SCATTER> ERROR Generator name '"//lnSplit(2)//"' is not unique."
      iErr = .true.
      return
    end if
  end do

  ! Generator type-dependent code
  select case (lnSplit(3)%get())
  case("ABSORBER")

    scatter_GENERATOR(scatter_nGENERATOR,2) = 1  ! Code for ABSORBER

  case("PPBEAMELASTIC")

    scatter_GENERATOR(scatter_nGENERATOR,2) = 10 ! Code for PPBEAMELASTIC
    if(nSplit < 8 .or. nSplit > 9) then
      write(lout,"(a)") "SCATTER> ERROR GEN PPBEAMELASTIC expected 8 or 9 arguments:"
      write(lout,"(a)") "SCATTER>       GEN name PPBEAMELASTIC a b1 b2 phi tmin (crossSection)"
      call prror(-1)
      iErr = .true.
      return
    end if

    ! Request space to store the arguments
    tmpIdx = scatter_nfData + 1
    scatter_GENERATOR(scatter_nGENERATOR,3) = tmpIdx ! Parameters
    scatter_GENERATOR(scatter_nGENERATOR,4) = 0      ! CrossSection
    scatter_nfData = scatter_nfData + nSplit - 3
    call alloc(scatter_fData, scatter_nfData, zero, "scatter_fData")

    call str_cast(lnSplit(4),scatter_fData(tmpIdx),iErr)   ! PPBEAMELASTIC a
    call str_cast(lnSplit(5),scatter_fData(tmpIdx+1),iErr) ! PPBEAMELASTIC b1
    call str_cast(lnSplit(6),scatter_fData(tmpIdx+2),iErr) ! PPBEAMELASTIC b2
    call str_cast(lnSplit(7),scatter_fData(tmpIdx+3),iErr) ! PPBEAMELASTIC phi
    call str_cast(lnSplit(8),scatter_fData(tmpIdx+4),iErr) ! PPBEAMELASTIC tmin

    if(nSplit == 9) then
      call str_cast(lnSplit(9),scatter_fData(tmpIdx+5),iErr)  ! crossSection
      scatter_fData(tmpIdx+5) = scatter_fData(tmpIdx+5)*c1m27 ! Scale to mb
      scatter_GENERATOR(scatter_nGENERATOR,4) = tmpIdx+5
    end if

    ! Check sanity of input values
    if(scatter_fData(tmpIdx+1) < pieni) then
      write(lout,"(a)") "SCATTER> ERROR GEN PPBEAMELASTIC 5th input (b1) must be larger than zero"
      iErr = .true.
      return
    end if
    if(scatter_fData(tmpIdx+2) < pieni) then
      write(lout,"(a)") "SCATTER> ERROR GEN PPBEAMELASTIC 6th input (b2) must be larger than zero"
      iErr = .true.
      return
    end if

  case("PYTHIASIMPLE")

#ifndef PYTHIA
    write(lout,"(a)") "SCATTER> ERROR GEN PYTHIA requested, but PYTHIA not compiled into SixTrack"
    iErr = .true.
    return
#endif
    scatter_GENERATOR(scatter_nGENERATOR,2) = 20 ! Code for PYTHIASIMPLE
    scatter_GENERATOR(scatter_nGENERATOR,3) = 0  ! Parameters
    scatter_GENERATOR(scatter_nGENERATOR,4) = 0  ! CrossSection

    if(nSplit == 4) then
      ! Save specified crossSection
      tmpIdx = scatter_nfData + 1
      scatter_GENERATOR(scatter_nGENERATOR,4) = tmpIdx
      scatter_nfData = scatter_nfData + 1
      call alloc(scatter_fData, scatter_nfData, zero, "scatter_fData")
      call str_cast(lnSplit(4),scatter_fData(tmpIdx),iErr)
      scatter_fData(tmpIdx) = scatter_fData(tmpIdx)*c1m27 ! Scale to mb
    end if

  case default

    write(lout,"(a)") "SCATTER> ERROR GEN name '"//lnSplit(3)//"' not recognized."
    iErr = .true.
    return

  end select

end subroutine scatter_parseGenerator

! =================================================================================================
! END Input Parser Functions
! =================================================================================================

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 02-11-2017
! =================================================================================================
subroutine scatter_thin(iElem, ix, turn)

  use string_tools
  use crcoall
  use mod_hions
  use mod_alloc
  use mod_common
  use mod_commonmn
  use mod_particles
  use numerical_constants
#ifdef HDF5
  use hdf5_output
#endif

  use collimation, only : do_coll, scatterhit, part_hit_pos, part_hit_turn

  implicit none

  integer, intent(in) :: iElem, ix, turn

  ! Temp variables
  integer          idElem, idPro, idGen, iError
  integer          i, j, k, iLost, nLost(9), nScatter(9), procID
  integer          tmpSeed1, tmpSeed2
  logical          isDiff, updateE, hasProc(9)
  real(kind=fPrec) s, t, dEE, dPP, theta, theta2
  real(kind=fPrec) crossSection, N, P
  real(kind=fPrec) rndPhi(napx), rndP(napx)
  real(kind=fPrec) scaling

  logical, allocatable :: pLost(:)
  logical, allocatable :: pScattered(:)

#ifdef HDF5
  ! For HDF5 it is best to write in chuncks, so we will make arrays of size napx
  integer                 :: iRecords(3,napx)
  real(kind=fPrec)        :: rRecords(8,napx)
  character(len=mNameLen) :: cRecords(3,napx)
  integer                 :: nRecords
#endif

  idElem  = scatter_elemPointer(ix)
  idPro   = scatter_ELEM(idElem,2)
  scaling = scatter_ELEM_scale(idElem)

  ! if(scatter_debug) then
  !   write(lout,"(a)")       "SCATTER> DEBUG In scatter_thin"
  !   write(lout,"(a,i0)")    "SCATTER> DEBUG  * ix      = ",ix
  !   write(lout,"(a)")       "SCATTER> DEBUG  * bez     = '"//trim(bez(ix))//"'"
  !   write(lout,"(a,i0)")    "SCATTER> DEBUG  * napx    = ",napx
  !   write(lout,"(a,i0)")    "SCATTER> DEBUG  * turn    = ",turn
  !   write(lout,"(a,e13.6)") "SCATTER> DEBUG  * scaling = ",scaling
  ! end if

  if(scaling <= pieni) then
    ! Skip the whole thing if the scaling is zero
    return
  end if

  ! Store the seeds in the randum number generator, and set ours
  call recuut(tmpSeed1,tmpSeed2)
  call recuin(scatter_seed1,scatter_seed2)

  if(scatter_allowLosses) then
    call alloc(pLost,napx,.false.,"pLost")
  end if
  call alloc(pScattered,napx,.false.,"pScattered")

  updateE = .false.

  ! Loop over generators
  do i=3,5

    nLost(:)    = 0
    nScatter(:) = 0
    hasProc(:)  = .false.

    t      = zero
    theta  = zero
    theta2 = zero
    dEE    = zero
    dPP    = zero

    idGen = scatter_ELEM(idElem,i)
    if(idGen == 0) exit ! No generator

    ! Generate a random phi
    call ranecu(rndPhi, napx, -1)
    call ranecu(rndP,   napx, -1)
    rndPhi = rndPhi*(two*pi)
#ifdef HDF5
    ! Reset counter
    nRecords = 0
#endif

    do j=1, napx

      ! Do not scatter the same particle twice in the same turn
      if(pScattered(j)) cycle

      ! Compute the cross section at this s
      ! (in most cases roughly equal for all particles; use mean x,y,xp,yp,E)
      crossSection = scatter_generator_getCrossSection(idPro,idGen,xv1(j),xv2(j),yv1(j),yv2(j),ejv(j))

      ! Ask profile for density at x,y
      N = scatter_profile_getDensity(idPro,xv1(j),xv2(j))

      ! Compute probability P
      P = (N*crossSection)*scaling

      ! If RNG > P -> go to next particle, else scatter
      if(rndP(j) > P) then
        cycle
      else
        pScattered(j) = .true.
      end if
      scatter_statScale(nlostp(j)) = scatter_statScale(nlostp(j)) / scaling

      ! Get scattering event
      call scatter_generator_getEvent(idGen,j,t,theta,dEE,dPP,procID,iLost,isDiff)
      hasProc(procID)  = .true.
      nScatter(procID) = nScatter(procID) + 1

      if(scatter_allowLosses .and. iLost == 1) then
        ! If lost, flag it, put no need to update energy and angle
        pLost(j)      = .true.
        nLost(procID) = nLost(procID) + 1
        theta         = zero
        rndPhi(j)     = zero
      else
        ! Calculate new energy, and scattering angle assuming energy >> mass
        if(abs(dPP) >= pieni) then
          ejfv(j) = (one+dPP)*ejfv(j) ! Momentum [MeV/c]
          updateE = .true.            ! Re-calculate energy-dependent vectors
        end if

        ! Update particle trajectory
        yv1(j) = theta*cos_mb(rndPhi(j)) + yv1(j)
        yv2(j) = theta*sin_mb(rndPhi(j)) + yv2(j)
      end if

      ! Output to file
#ifdef HDF5
      if(h5_useForSCAT) then
        nRecords = nRecords + 1
        iRecords(1,nRecords) = nlostp(j)
        iRecords(2,nRecords) = turn
        cRecords(1,nRecords) = bez(ix)
        cRecords(2,nRecords) = trim(scatter_cData(scatter_GENERATOR(idGen,1)))
        cRecords(3,nRecords) = scatter_procNames(procID)
        iRecords(3,nRecords) = iLost
        rRecords(1,nRecords) = t
        rRecords(2,nRecords) = dEE
        rRecords(3,nRecords) = dPP
        rRecords(4,nRecords) = theta
        rRecords(5,nRecords) = rndPhi(j)
        rRecords(6,nRecords) = N
        rRecords(7,nRecords) = P
        rRecords(8,nRecords) = scatter_statScale(nlostp(j))
      else
#endif
        write(scatter_logFile,"(2(1x,i8),2(1x,a20),1x,a8,1x,i4,1x,f13.3,7(1x,1pe16.9))") &
          nlostp(j), turn, bez(ix)(1:20), chr_rPad(trim(scatter_cData(scatter_GENERATOR(idGen,1))),20), &
          scatter_procNames(procID), iLost, t, dEE, dPP, theta, rndPhi(j), N, P, scatter_statScale(nlostp(j))
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
    end do ! END Loop over particles

#ifdef HDF5
    if(h5_useForSCAT) then
      call h5_prepareWrite(scatter_logDataSet, nRecords)
      call h5_writeData(scatter_logDataSet, 1,  nRecords, iRecords(1,1:nRecords))
      call h5_writeData(scatter_logDataSet, 2,  nRecords, iRecords(2,1:nRecords))
      call h5_writeData(scatter_logDataSet, 3,  nRecords, cRecords(1,1:nRecords))
      call h5_writeData(scatter_logDataSet, 4,  nRecords, cRecords(2,1:nRecords))
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

      do k=1,9
        if(hasProc(k)) then
          call h5_prepareWrite(scatter_sumDataSet, 1)
          call h5_writeData(scatter_sumDataSet, 1, 1, turn)
          call h5_writeData(scatter_sumDataSet, 2, 1, bez(ix))
          call h5_writeData(scatter_sumDataSet, 3, 1, scatter_cData(scatter_GENERATOR(idGen,1)))
          call h5_writeData(scatter_sumDataSet, 4, 1, scatter_procNames(k))
          call h5_writeData(scatter_sumDataSet, 5, 1, nScatter(k))
          call h5_writeData(scatter_sumDataSet, 6, 1, nLost(k))
          call h5_writeData(scatter_sumDataSet, 7, 1, crossSection*c1e27)
          call h5_writeData(scatter_sumDataSet, 8, 1, scaling)
          call h5_finaliseWrite(scatter_sumDataSet)
        end if
      end do
    else
#endif

      do k=1,9
        if(hasProc(k)) then
          write(scatter_sumFile,"(1x,i8,2(1x,a20),1x,a8,2(1x,i8),2(1x,f13.6))") &
            turn, bez(ix)(1:20), chr_rPad(trim(scatter_cData(scatter_GENERATOR(idGen,1))),20), &
            scatter_procNames(k), nScatter(k), nLost(k), crossSection*c1e27, scaling
#ifdef CR
            scatter_sumFilePos = scatter_sumFilePos + 1
#endif
        end if
      end do
      flush(scatter_logFile)
      flush(scatter_sumFile)

#ifdef HDF5
    end if
#endif
  end do ! END Loop over generators

  if(scatter_allowLosses) then
    call compactArrays(pLost)
  end if

  if(updateE) then
    call part_updatePartEnergy(2)
  end if

#ifdef CR
  endfile(scatter_logFile,iostat=iError)
  backspace(scatter_logFile,iostat=iError)
  endfile(scatter_sumFile,iostat=iError)
  backspace(scatter_sumFile,iostat=iError)
#endif

  ! Restore seeds in random generator
  call recuut(scatter_seed1,scatter_seed2)
  call recuin(tmpSeed1,tmpSeed2)

  if(scatter_allowLosses) then
    call dealloc(pLost,"pLost")
  end if
  call dealloc(pScattered,"pScattered")

end subroutine scatter_thin

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 02-11-2017
! =================================================================================================
real(kind=fPrec) function scatter_profile_getDensity(idPro, x, y) result(retval)

  use string_tools
  use crcoall
  use mod_common
  use numerical_constants, only : pi

  implicit none

  integer,          intent(in) :: idPro
  real(kind=fPrec), intent(in) :: x, y

  real(kind=fPrec) beamtot, sigmaX, sigmaY, offsetX, offsetY
  integer tmpIdx

  tmpIdx = scatter_PROFILE(idPro,3)

  select case(scatter_PROFILE(idPro,2))
  case (1)  ! FLAT
    retval  = scatter_fData(tmpIdx)

  case (2)  ! FIXED
    retval  = scatter_fData(tmpIdx)

  case (10) ! GAUSS1
    beamtot = scatter_fData(tmpIdx)
    sigmaX  = scatter_fData(tmpIdx + 1)
    sigmaY  = scatter_fData(tmpIdx + 2)
    offsetX = scatter_fData(tmpIdx + 3)
    offsetY = scatter_fData(tmpIdx + 4)
    retval  = ((beamtot/(two*(pi*(sigmaX*sigmaY)))) &
      * exp_mb(-half*((x-offsetX)/sigmaX)**2))      &
      * exp_mb(-half*((y-offsetY)/sigmaY)**2)

  case default
    write(lout,"(a,i0,a)") "SCATTER> ERROR Type ", scatter_PROFILE(idPro,2)," for profile '"//&
      trim(scatter_cData(scatter_PROFILE(idPro,1)))//"' not understood."
    call prror(-1)
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
real(kind=fPrec) function scatter_generator_getCrossSection(idPro, genID, x, y, xp, yp, E)

  use string_tools
  use crcoall

  implicit none

  integer,          intent(in) :: idPro, genID
  real(kind=fPrec), intent(in) :: x, y, xp, yp, E

  ! Temporary variables
  integer          tmpIdx
  real(kind=fPrec) xpTarget, ypTarget, ETarget

  ! Calculate S
  call scatter_profile_getParticle(idPro, x, y, xpTarget, ypTarget, ETarget)

  ! Calculate the cross section as function of S
  select case(scatter_GENERATOR(genID,2))
  case(1)  ! ABSORBER
    !...

  case(10) ! PPBEAMELASTIC
    tmpIdx = scatter_GENERATOR(genID,4)
    if(tmpIdx == 0) then
      scatter_generator_getCrossSection = 30d-27
    else
      scatter_generator_getCrossSection = scatter_fData(tmpIdx)
    end if

  case(20) ! PYTHIASIMPLE
    tmpIdx = scatter_GENERATOR(genID,4)
    if(tmpIdx == 0) then
      scatter_generator_getCrossSection = 30d-27
    else
      scatter_generator_getCrossSection = scatter_fData(tmpIdx)
    end if

  case default
    write(lout,"(a,i0,a)") "SCATTER> ERROR Type ",scatter_PROFILE(idPro,2)," for profile '"//&
      trim(scatter_cData(scatter_PROFILE(idPro,1)))//"' not understood."
    call prror(-1)

  end select

end function scatter_generator_getCrossSection

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 02-11-2017
! =================================================================================================
subroutine scatter_generator_getEvent(genID, partID, t, theta, dEE, dPP, procID, iLost, isDiff)

  use crcoall
  use mod_commonmn
  use, intrinsic :: iso_c_binding

  implicit none

  integer,          intent(in)  :: genID
  integer,          intent(in)  :: partID
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

  select case(scatter_GENERATOR(genID,2))
  case(1)  ! ABSORBER

    procID = scatter_idAbsorb

  case(10) ! PPBEAMELASTIC

    tmpIdx = scatter_GENERATOR(genID,3)
    a      = scatter_fData(tmpIdx)
    b1     = scatter_fData(tmpIdx+1)
    b2     = scatter_fData(tmpIdx+2)
    phi    = scatter_fData(tmpIdx+3)
    tmin   = scatter_fData(tmpIdx+4)

    t      = scatter_generator_getPPElastic(a, b1, b2, phi, tmin)
    t      = t*c1e6                                      ! Scale return variable to MeV^2
    theta  = acos_mb(one - (t/(2*ejfv(partID)**2)))*c1e3 ! Get angle from t
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
          write(lout,"(a)") "SCATTER> ERROR Particle lost, but losses not explicitly allowed in fort.3"
          call prror(-1)
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
          write(lout,"(a)") "SCATTER> ERROR Particle lost, but losses not explicitly allowed in fort.3"
          call prror(-1)
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
    write(lout,"(a,i0,a)") "SCATTER> ERROR This version of SixTrack was built without PYTHIA support,"
    call prror(-1)
#endif

  case default
    write(lout,"(a,i0,a)") "SCATTER> ERROR Generator type ",scatter_GENERATOR(genID,2)," not understood"
    call prror(-1)

  end select

end subroutine scatter_generator_getEvent

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 02-11-2017
!  Based on Helmut Burkhardt's code presented on 05-04-2017
!  "Elastic pp scattering estimates and simulation relevant for burn-off"
!  https://indico.cern.ch/event/625576/
! =================================================================================================
real(kind=fPrec) function scatter_generator_getPPElastic(a, b1, b2, phi, tmin) result(t)

  use crcoall

  implicit none

  real(kind=fPrec), intent(in) :: a, b1, b2, phi, tmin

  ! Temp Variables
  integer          nItt, maxItt
  real(kind=fPrec) g1, g2, g3, gg, prob3, invB1, invB2, rndArr(3)

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
    write(lout,"(a)")      "SCATTER> ERROR in generator PPBEAMELASTIC"
    write(lout,"(a,i0,a)") "SCATTER>       Limit of ",maxItt," misses in generator loop reached."
    call prror(-1)
  end if

end function scatter_generator_getPPElastic

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-26
!  Called from CRCHECK; reads the _CR arrays for scatter from file-
!  Sets readerr=.true. if something goes wrong while reading.
! =================================================================================================
#ifdef CR
subroutine scatter_crcheck_readdata(fileUnit, readErr)

  use crcoall
  use mod_alloc

  implicit none

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: readErr

  integer j

  read(fileUnit, err=10, end=10) scatter_logFilePos_CR, scatter_sumFilePos_CR
  read(fileUnit, err=10, end=10) scatter_seed1_CR, scatter_seed2_CR

  readErr = .false.
  return

10 continue
  write(lout,"(a,i0)") "READERR in scatter_crcheck; fileUnit = ",fileUnit
  write(93,  "(a,i0)") "READERR in scatter_crcheck; fileUnit = ",fileUnit
  readErr = .true.

end subroutine scatter_crcheck_readdata

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-26
!  Called from CRCHECK; resets the position of scatter_log.dat
! =================================================================================================
subroutine scatter_crcheck_positionFiles

  use crcoall
  use file_units
  use mod_units

  implicit none

  logical isOpen, fErr
  integer iError
  integer j
  character(len=1024) aRecord

  if(scatter_logFile == -1) call funit_requestUnit("scatter_log.dat",scatter_logFile)
  inquire(unit=scatter_logFile, opened=isOpen)
  if(isOpen) then
    write(93,"(a)")      "SIXTRACR> ERROR CRCHECK FAILED while repsositioning 'scatter_log.dat'"
    write(93,"(a,i0,a)") "SIXTRACR>       UNIT ",scatter_logFile," already in use!"

    endfile(93,iostat=iError)
    backspace(93,iostat=iError)

    write(lout,"(a)") "SIXTRACR> CRCHECK failure positioning 'scatter_log.dat'"
    call prror(-1)
  end if

  if(scatter_logFilePos_CR /= -1) then
    call units_openUnit(unit=scatter_logFile,fileName="scatter_log.dat",formatted=.true.,mode="rw",err=fErr,status="old")
    if(fErr) goto 10
    scatter_logFilePos = 0
    do j=1, scatter_logFilePos_CR
      read(scatter_logFile,"(a1024)",end=10,err=10,iostat=iError) aRecord
      scatter_logFilePos = scatter_logFilePos + 1
    end do
    endfile(scatter_logFile,iostat=iError)
    close(scatter_logFile)

    call units_openUnit(unit=scatter_logFile,fileName="scatter_log.dat",formatted=.true.,mode="w+",err=fErr,status="old")
    if(fErr) goto 10
    write(97,"(2(a,i0))") "SIXTRACR> CRCHECK sucessfully repositioned 'scatter_log.dat': "//&
      "scatter_logFilePos = ",scatter_logFilePos,", scatter_logFilePos_CR = ",scatter_logFilePos_CR
    endfile(93,iostat=iError)
    backspace(93,iostat=iError)

  else
    write(93,"(a,i0)") "SIXTRACR> CRCHECK did not attempt repositioning "// &
      "of 'scatter_log.dat', scatter_logFilePos_CR = ",scatter_logFilePos_CR
    write(93,"(a)")    "SIXTRACR> If anything has been written to the file, "// &
      "it will be correctly truncated in scatter_initialise."
    endfile(93,iostat=iError)
    backspace(93,iostat=iError)
  end if

  if(scatter_sumFile == -1) call funit_requestUnit("scatter_summary.dat",scatter_sumFile)
  inquire(unit=scatter_sumFile, opened=isOpen)
  if(isOpen) then
    write(93,"(a)")      "SIXTRACR> ERROR CRCHECK FAILED while repsositioning 'scatter_summary.dat'"
    write(93,"(a,i0,a)") "SIXTRACR>       UNIT ",scatter_sumFile," already in use!"

    endfile(93,iostat=iError)
    backspace(93,iostat=iError)

    write(lout,"(a)") "SIXTRACR> CRCHECK failure positioning 'scatter_summary.dat'"
    call prror(-1)
  end if

  if(scatter_sumFilePos_CR /= -1) then
    call units_openUnit(unit=scatter_sumFile,fileName="scatter_summary.dat",formatted=.true.,mode="rw",err=fErr,status="old")
    if(fErr) goto 10
    scatter_sumFilePos = 0
    do j=1, scatter_sumFilePos_CR
      read(scatter_sumFile,"(a1024)",end=10,err=10,iostat=iError) aRecord
      scatter_sumFilePos = scatter_sumFilePos + 1
    end do
    endfile(scatter_sumFile,iostat=iError)
    close(scatter_sumFile)

    call units_openUnit(unit=scatter_sumFile,fileName="scatter_summary.dat",formatted=.true.,mode="w+",err=fErr,status="old")
    if(fErr) goto 10
    write(97,"(2(a,i0))") "SIXTRACR> CRCHECK sucessfully repositioned 'scatter_summary.dat': "//&
      "scatter_sumFilePos = ",scatter_sumFilePos,", scatter_sumFilePos_CR = ",scatter_sumFilePos_CR
    endfile(93,iostat=iError)
    backspace(93,iostat=iError)

  else
    write(93,"(a,i0)") "SIXTRACR> CRCHECK did not attempt repositioning "// &
      "of 'scatter_summary.dat', scatter_sumFilePos_CR = ",scatter_sumFilePos_CR
    write(93,"(a)")    "SIXTRACR> If anything has been written to the file, "// &
      "it will be correctly truncated in scatter_initialise."
    endfile(93,iostat=iError)
    backspace(93,iostat=iError)
  end if

  return

10 continue
  write(93,"(a,i0)")    "SIXTRACR> ERROR reading 'scatter_log.dat' or 'scatter_summary.dat', iostat=",iError
  write(93,"(2(a,i0))") "SIXTRACR> scatter_logFilePos = ",scatter_logFilePos,", scatter_logFilePos_CR = ",scatter_logFilePos_CR
  write(93,"(2(a,i0))") "SIXTRACR> scatter_sumFilePos = ",scatter_sumFilePos,", scatter_sumFilePos_CR = ",scatter_sumFilePos_CR
  endfile(93,iostat=iError)
  backspace(93,iostat=iError)
  write(lout,"(a)")"SIXTRACR> ERROR CRCHECK failure positioning 'scatter_log.dat' or 'scatter_summary.dat'."
  call prror(-1)

end subroutine scatter_crcheck_positionFiles

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-26
!  Called from CRPOINT; write checkpoint data to fort.95/96
! =================================================================================================
subroutine scatter_crpoint(fileUnit, writeErr, iError)

  use crcoall

  implicit none

  integer, intent(in)    :: fileUnit
  logical, intent(out)   :: writeErr
  integer, intent(inout) :: iError

  integer j

  write(fileunit,err=10,iostat=iError) scatter_logFilePos, scatter_sumFilePos
  write(fileunit,err=10,iostat=iError) scatter_seed1, scatter_seed2
  endfile(fileUnit,iostat=iError)
  backspace(fileUnit,iostat=iError)

  return

10 continue
  writeErr = .true.

end subroutine scatter_crpoint
#endif
! End of CR

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-20
! =================================================================================================
subroutine scatter_dumpData

  use crcoall

  implicit none

  integer i

  write(lout,"(a)")       "SCATTER> DEBUG Data Dump"

  write(lout,"(a)")       "Options:"
  write(lout,"(a,l2)")    "scatter_active =", scatter_active
  write(lout,"(a,l2)")    "scatter_debug  =", scatter_debug
  write(lout,"(2(a,i8))") "Seeds          =", scatter_seed1, ";", scatter_seed2

  write(lout,"(a)")       "Arrays:"

  write(lout,"(a,2(i3,a))") "scatter_ELEM: (",scatter_nELEM,",",5,"):"
  do i=1,scatter_nELEM
    write(lout,"(i4,a,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)") i,":",scatter_ELEM(i,:)
  end do

  write(lout,"(a,i3,a)") "scatter_ELEM_scale: (",scatter_nELEM,"):"
  do i=1,scatter_nELEM
    write(lout,"(i4,a,e14.7)") i,":",scatter_ELEM_scale(i)
  end do

  write(lout,"(a,2(i3,a))") "scatter_PROFILE: (", scatter_nPROFILE,",",5,"):"
  do i=1,scatter_nPROFILE
    write(lout,"(i4,a,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)") i,":",scatter_PROFILE(i,:)
  end do

  write(lout,"(a,2(i3,a))") "scatter_GENERATOR: (",scatter_nGENERATOR,",",5,"):"
  do i=1,scatter_nGENERATOR
    write(lout,"(i4,a,5(1x,i3))") i,":",scatter_GENERATOR(i,:)
  end do

  write(lout,"(a,i3,a)") "scatter_iData: (",scatter_niData,"):"
  do i=1,scatter_niData
    write(lout,"(i4,a,i4)") i,":",scatter_iData(i)
  end do

  write(lout,"(a,i3,a)") "scatter_fData: (",scatter_nfData,"):"
  do i=1,scatter_nfData
    write(lout,"(i4,a,e14.7)") i,":",scatter_fData(i)
  end do

  write(lout,"(a,i3,a)") "scatter_cData: (",scatter_ncData,"):"
  do i=1,scatter_ncData
    write(lout,"(i4,a)") i,": '"//chr_trimZero(scatter_cData(i))//"'"
  end do

  write(lout,"(a)") "SCATTER> DEBUG END DUMP"

end subroutine scatter_dumpData

end module scatter
