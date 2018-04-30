#ifdef HDF5

! ================================================================================================ !
!  HDF5 Module
! ~~~~~~~~~~~~~
!  Written by:    Veronica Berglyd Olsen, BE-ABP-HSS, 2018
!  Last Modified: 2018-04-13
!
!  Alternative output format.
!  Controlled by the HDF5 input block.
! ================================================================================================ !
module hdf5_output
  
  use hdf5
  use floatPrecision
  use end_sixtrack
  use crcoall
  use strings
  use mod_alloc
  
  implicit none
  
  ! Common Settings
  logical,      public,  save :: h5_isActive    ! Existence of the HDF5 block
  logical,      public,  save :: h5_debugOn     ! HDF5 debug flag present
  logical,      public,  save :: h5_isReady     ! HDF5 file is open and ready for input
  logical,      private, save :: h5_useDouble   ! Whether to use double precision or not
  logical,      private, save :: h5_doTruncate  ! Whether or not to truncate previous file if it exists
  type(string), private, save :: h5_fileName    ! The HDF5 output file name
  type(string), private, save :: h5_rootPath    ! The root group where the data for this session is stored
  
  ! Input Block Switches
  logical, public, save :: h5_useForCOLL
  logical, public, save :: h5_useForDUMP
  logical, public, save :: h5_useForSCAT
  
  ! Runtime Variables
  logical, private, save :: h5_fileIsOpen  ! True if file is open.
  integer, public,  save :: h5_fileError   ! For errors related to file essentials (critical)
  integer, public,  save :: h5_dataError   ! For errors related to datasets
  
  ! HDF5 File/Group IDs
  integer(HID_T), public,  save :: h5_fileID  ! The internal ID of the file
  integer(HID_T), public,  save :: h5_rootID  ! The internal ID of the root group
  integer(HID_T), public,  save :: h5_collID  ! The internal ID of the collimation group
  integer(HID_T), public,  save :: h5_dumpID  ! The internal ID of the dump group
  integer(HID_T), public,  save :: h5_scatID  ! The internal ID of the scatter group
  
  ! HDF5 Internals
  integer(HID_T), private, save :: h5_plistID ! Dataset transfer property
  
  ! DataSet ID Mappings
  type(string),   allocatable, private, save :: h5_dataSetName(:)
 !type(string),   allocatable, private, save :: h5_dataSetPath(:)
  integer(HID_T), allocatable, private, save :: h5_dataSetMap(:)
  integer,                     private, save :: h5_nDataSets
  
  ! Default Group Names
  character(len=11), parameter :: h5_collGroup = "collimation"
  character(len=4),  parameter :: h5_dumpGroup = "dump"
  character(len=7),  parameter :: h5_scatGroup = "scatter"
  
  integer, parameter :: h5_typeInt  = 1 ! Size is fixed to fortran native integer
  integer, parameter :: h5_typeReal = 2 ! Size is determined by fPrec
  integer, parameter :: h5_typeChar = 3 ! Any size is valid
  
  ! Field Definitions
  type, public :: h5_dataField
    character(len=:), allocatable, public :: name
    integer,                       public :: type
    integer,                       public :: size
    integer(HID_T),                public :: typeH5 = 0
    integer(HID_T),                public :: typeID = 0
  end type h5_dataField
  
  ! Dataset Container
  type, public :: h5_dataSet
    character(len=:),   allocatable, public :: name
    character(len=:),   allocatable, public :: path
    integer(HID_T),                  public :: dataID  = 0
    integer(HID_T),                  public :: spaceID = 0
    integer(HID_T),                  public :: dtypeID = 0
    type(h5_dataField), allocatable, public :: fields(:)
  end type h5_dataSet
  
  interface h5_writeValue
    module procedure h5_writeValue_real32
    module procedure h5_writeValue_real64
    module procedure h5_writeValue_int
  end interface h5_writeValue
  
contains

! ================================================================================================ !
!  Set Initial Values
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-20
! ================================================================================================ !
subroutine h5_comnul
  
  h5_isActive   = .false.
  h5_debugOn    = .false.
  h5_isReady    = .false.
  h5_useDouble  = .true.
  h5_doTruncate = .false.
  h5_fileName   = ""
  h5_rootPath   = ""
  
  h5_useForCOLL = .false.
  h5_useForDUMP = .false.
  h5_useForSCAT = .false.
  
  h5_fileIsOpen = .false.
  h5_fileError  = 0
  h5_dataError  = 0
  
  h5_fileID     = 0
  h5_rootID     = 0
  h5_collID     = 0
  h5_dumpID     = 0
  h5_scatID     = 0
  
  h5_nDataSets  = 0
  
end subroutine h5_comnul

! ================================================================================================ !
!  HDF5 Initialisation
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
! ================================================================================================ !
subroutine h5_initHDF5()
  
  use end_sixtrack
  
  implicit none
  
  call h5open_f(h5_fileError)
  call h5pcreate_f(H5P_DATASET_XFER_F, h5_plistID, h5_fileError)
  call h5pset_preserve_f(h5_plistID, .true., h5_fileError)
  if(h5_fileError == -1) then
    write(lout,"(a)") "HDF5> ERROR Failed to initialise Fortran HDF5."
    call prror(-1)
  end if
  
end subroutine h5_initHDF5

! ================================================================================================ !
!  HDF5 File Initialisation
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
! ================================================================================================ !
subroutine h5_openFile()
  
  use end_sixtrack
  
  implicit none
  
  integer accessFlag
  
  if(.not. h5_isActive) return
  
  if(h5_doTruncate) then
    accessFlag = H5F_ACC_TRUNC_F
    if(h5_debugOn) then
      write(lout,"(a)") "HDF5> DEBUG Truncating previous file if it exists."
    end if
  else
    accessFlag = H5F_ACC_EXCL_F
  end if
  
  call h5fcreate_f(h5_fileName%chr, accessFlag, h5_fileID, h5_fileError)
  if(h5_fileError == -1) then
    write(lout,"(3a)") "HDF5> ERROR Failed to open HDF5 file '",h5_fileName%chr,"'."
    call prror(-1)
  end if
  write(lout,"(3a)") "HDF5> Created/opened HDF5 file '",h5_fileName%chr,"'."
  
  ! If a root group was requested, create it and save the rootID
  ! Otherwise, use the fileID as the rootID
  if(h5_rootPath%chr /= "") then
    call h5gcreate_f(h5_fileID, h5_rootPath%chr, h5_rootID, h5_fileError)
    if(h5_fileError == -1) then
      write(lout,"(3a)") "HDF5> ERROR Failed to create root group '",h5_rootPath%chr,"'."
      call prror(-1)
    end if
    write(lout,"(3a)") "HDF5> Created root group '",h5_rootPath%chr,"'."
  else
    h5_rootID = h5_fileID
  end if
  
  h5_isReady = .true.
  
end subroutine h5_openFile

! ================================================================================================ !
!  HDF5 Finalise
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
! ================================================================================================ !
subroutine h5_closeHDF5()
  
  use end_sixtrack
  
  implicit none
  
  if(.not. h5_isReady) return
  
  call h5fclose_f(h5_fileID, h5_fileError)
  if(h5_fileError == -1) then
    write(lout,"(a)") "HDF5> ERROR Failed to close HDF5 file."
    call prror(-1)
  end if
  
  call h5close_f(h5_fileError)
  if(h5_fileError == -1) then
    write(lout,"(a)") "HDF5> ERROR Failed to close Fortran HDF5."
    call prror(-1)
  end if
  
  write(lout,"(a)") "HDF5> Closed HDF5 file."
  
end subroutine h5_closeHDF5

! ================================================================================================ !
!  Block-Wise Initialisation
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
! ================================================================================================ !
subroutine h5_initForScatter()
  
  call h5gcreate_f(h5_rootID, h5_scatGroup, h5_scatID, h5_fileError)
  if(h5_fileError == -1) then
    write(lout,"(a)") "HDF5> ERROR Failed to create scatter group '"//h5_scatGroup//"'."
    call prror(-1)
  end if
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Group initialised for scatter."
  end if

end subroutine h5_initForScatter

! ================================================================================================ !
!  Create DataSet
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-30
! ================================================================================================ !
subroutine h5_createDataSet(setName, setGroup, setFields, dataSet)
  
  use mod_alloc
  use end_sixtrack
  
  implicit none
  
  ! Routine Variables
  character(len=*) ,               intent(in)    :: setName
  integer(HID_T),                  intent(in)    :: setGroup
  type(h5_dataField), allocatable, intent(inout) :: setFields(:)
  type(h5_dataSet),                intent(out)   :: dataSet
  
  ! Internal Variables
  integer(HID_T),  allocatable :: fieldType(:)
  integer(SIZE_T), allocatable :: fieldSize(:)
  
  integer(HID_T)   spaceID, dtypeID, dataID, tmpID
  integer(HSIZE_T) spaceSize(1)
  integer(HSIZE_T) spaceMaxSize(1)
  integer(SIZE_T)  memOffset, tmpOffset, memSize, tmpSize
  
  integer i, nFields
  
  ! Init Variables
  spaceSize    = (/1000/)
  ! spaceMaxSize = (/H5S_UNLIMITED_F/)
  spaceMaxSize = (/1000/)
  memOffset    = 0
  memSize      = 0
  tmpOffset    = 0
  tmpSize      = 0
  
  ! ***** Routine Code *****
  
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Creating dataset '"//setName//"'"
  end if
  
  ! h5_nDataSets = h5_nDataSets + 1
  
  ! Translate data types
  nFields = size(setFields)
  allocate(fieldType(nFields))
  allocate(fieldSize(nFields))
  do i=1,nFields
    select case (setFields(i)%type)
    case(h5_typeInt)
      fieldType(i) = H5T_NATIVE_INTEGER
      call h5tget_size_f(fieldType(i), fieldSize(i), h5_dataError)
    case(h5_typeReal)
      if(h5_useDouble) then
        fieldType(i) = H5T_NATIVE_DOUBLE
      else
        fieldType(i) = H5T_NATIVE_REAL
      end if
      call h5tget_size_f(fieldType(i), fieldSize(i), h5_dataError)
    case(h5_typeChar)
      tmpSize = setFields(i)%size
      call h5tcopy_f(H5T_NATIVE_CHARACTER, tmpID, h5_dataError)
      call h5tset_size_f(tmpID, tmpSize, h5_dataError)
      call h5tget_size_f(tmpID, fieldSize(i), h5_dataError)
      fieldType(i) = tmpID
    end select
    memSize = memSize + fieldSize(i)
  end do
  
  if(h5_debugOn) then
    write(lout,"(a,i0,a)") "HDF5> DEBUG Data set size is ",memSize," bytes per row."
  end if
  
  ! Create the dataset
  call h5screate_simple_f(1, spaceSize, spaceID, h5_dataError, spaceMaxSize)
  call h5tcreate_f(H5T_COMPOUND_F, memSize, dtypeID, h5_dataError)
  
  do i=1,nFields
    
    ! Create in File
    call h5tinsert_f(dtypeID, setFields(i)%name, memOffset, fieldType(i), h5_dataError)
    if(h5_dataError /= 0) write(lout,"(a)") "ERROR!"
    memOffset = memOffset + fieldSize(i)
    
    ! Create in Memory
    setFields(i)%typeH5 = fieldType(i)
    call h5tcreate_f(H5T_COMPOUND_F, fieldSize(i), setFields(i)%typeID, h5_dataError)
    call h5tinsert_f(setFields(i)%typeID, setFields(i)%name, tmpOffset, setFields(i)%typeH5, h5_dataError)
    
  end do
  
  call h5dcreate_f(setGroup, setName, dtypeID, spaceID, dataID, h5_dataError)
  
  dataSet%name    = setName
  dataSet%path    = ""
  dataSet%dataID  = dataID
  dataSet%spaceID = spaceID
  dataSet%dtypeID = dtypeID
  dataSet%fields  = setFields
  
  ! ! Save the new dataset
  ! call resize(h5_dataSetName, h5_nDataSets, string(setName), "h5_dataSetName")
  ! call resize(h5_dataSetMap,  h5_nDataSets, dataID,          "h5_dataSetMap")
  
  ! Clean up
  deallocate(fieldType)
  deallocate(fieldSize)
  
end subroutine h5_createDataSet

! ================================================================================================ !
!  Close DataSet
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-30
! ================================================================================================ !
subroutine h5_closeDataSet(dataSet)
  
  type(h5_dataSet), intent(inout) :: dataSet
  
  integer i
  
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Closing dataset '"//dataSet%name//"'"
  end if
  
  call h5dclose_f(dataSet%dataID, h5_dataError)
  call h5sclose_f(dataSet%spaceID, h5_dataError)
  call h5tclose_f(dataSet%dtypeID, h5_dataError)
  do i=1,size(dataSet%fields)
    call h5tclose_f(dataSet%fields(i)%typeID, h5_dataError)
  end do
  
end subroutine h5_closeDataSet

! ================================================================================================ !
!  Write Single Value
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-30
! ================================================================================================ !
subroutine h5_writeValue_real32(dataSet, colID, theValue)

  type(h5_dataSet),  intent(inout) :: dataSet
  integer,           intent(in)    :: colID
  real(kind=real32), intent(in)    :: theValue
  
  integer(HSIZE_T) dataSize(1)
  dataSize(1) = 1
  
  call h5dwrite_f(dataSet%dataID, dataSet%fields(colID)%typeID, theValue, dataSize, h5_dataError, xfer_prp=h5_plistID)
  
end subroutine h5_writeValue_real32

subroutine h5_writeValue_real64(dataSet, colID, theValue)

  type(h5_dataSet),  intent(inout) :: dataSet
  integer,           intent(in)    :: colID
  real(kind=real64), intent(in)    :: theValue
  
  integer(HSIZE_T) dataSize(1)
  dataSize(1) = 1
  
  call h5dwrite_f(dataSet%dataID, dataSet%fields(colID)%typeID, theValue, dataSize, h5_dataError, xfer_prp=h5_plistID)
  
end subroutine h5_writeValue_real64

subroutine h5_writeValue_int(dataSet, colID, theValue)
  
  type(h5_dataSet),  intent(inout) :: dataSet
  integer,           intent(in)    :: colID
  integer,           intent(in)    :: theValue
  
  integer(HSIZE_T) dataSize(1)
  dataSize(1) = 1
  
  call h5dwrite_f(dataSet%dataID, dataSet%fields(colID)%typeID, theValue, dataSize, h5_dataError, xfer_prp=h5_plistID)
  
end subroutine h5_writeValue_int

! ================================================================================================ !
!  HDF5 Input File Parsing
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-13
! ================================================================================================ !
subroutine h5_parseInputLine(inLine)
  
  use string_tools
  use end_sixtrack
  
  implicit none
  
  type(string), intent(in)  :: inLine
  
  type(string), allocatable :: lnSplit(:)
  integer nSplit
  
  integer i
  
  ! Split the input line
  call str_split(inLine,lnSplit,nSplit)
  
  if(nSplit == 0) then
    if(h5_debugOn) then
      write (lout,"(a,i3,a)") "HDF5> DEBUG Input line len=",len(inLine),": '"//inLine%strip()//"'."
      write (lout,"(a)")      "HDF5> DEBUG  * No fields found."
    end if
    return
  end if
  
  ! Report if debugging is ON
  if(h5_debugOn) then
    write (lout,"(a,i3,a)")  "HDF5> DEBUG Input line len=",len(inLine),": '"//inLine%strip()//"'."
    write (lout,"(a,i2,a)") ("HDF5> DEBUG  * Field(",i,") = '"//lnSplit(i)//"'",i=1,nSplit)
  end if
  
  select case(lnSplit(1)%chr)
  
  case("DEBUG")
    h5_debugOn = .true.
    write(lout,"(a)") "HDF5> HDF5 block debugging is ON."
  
  case("SINGLE")
    h5_useDouble = .false.
    write(lout,"(a)") "HDF5> HDF5 will use single precision."
  
  case("DOUBLE")
    h5_useDouble = .true.
    write(lout,"(a)") "HDF5> HDF5 will use double precision."
  
  case("FILE")
    if(nSplit < 2 .or. nSplit > 3) then
      write(lout,"(a,i2,a)") "HDF5> ERROR: FILE statement takes 1 or 2 input parameters, ",(nSplit-1)," given."
      write(lout,"(a)")      "HDF5>        Valid input is FILE filename [truncate]"
      call prror(-1)
    end if
    if(nSplit == 3) then
      read(lnSplit(3)%chr,*) h5_doTruncate
    else
      h5_doTruncate = .false.
    end if
    h5_fileName = str_stripQuotes(lnSplit(2))
    write(lout, "(a)") "HDF5> Output file name set to: '"//h5_fileName//"'."
  
  case("ROOT")
    if(nSplit /= 2) then
      write(lout,"(a,i2,a)") "HDF5> ERROR: ROOT statement takes 1 input parameter, ",(nSplit-1)," given."
      call prror(-1)
    end if
    if(str_inStr(lnSplit(2)," ") /= 0) then
      write(lout,"(a)") "HDF5> ERROR: ROOT group name cannot contain a space."
      call prror(-1)
    end if
    if(str_inStr(lnSplit(2),"/") /= 0) then
      write(lout,"(a)") "HDF5> ERROR: ROOT group name cannot contain a slash."
      call prror(-1)
    end if
    h5_rootPath = str_stripQuotes(lnSplit(2))
    write(lout, "(a)") "HDF5> Root group set to: '"//h5_rootPath//"'."
  
  case("ENABLE")
  
    if(nSplit /= 2) then
      write(lout,"(a,i2,a)") "HDF5> ERROR: ENABLE statement takes 1 input parameter, ",(nSplit-1)," given."
      call prror(-1)
    end if
    if(len(lnSplit(2)%chr) < 4) then
      write(lout,"(a,i2,a)") "HDF5> ERROR: ENABLE argument must be at least 4 characters."
      call prror(-1)
    end if
    
    select case(lnSplit(2)%chr(1:4))
    case("COLL")
      h5_useForCOLL = .true.
      write(lout,"(3a)") "HDF5> HDF5 is enabled for COLLIMATION."
    case("DUMP")
      h5_useForDUMP = .true.
      write(lout,"(3a)") "HDF5> HDF5 is enabled for DUMP."
    case("SCAT")
      h5_useForSCAT = .true.
      write(lout,"(a)") "HDF5> HDF5 is enabled for SCATTER."
    case default
      write(lout,"(a)") "HDF5> ERROR: HDF5 output is not available for "//lnSplit(2)%chr(1:4)//" blocks."
      call prror(-1)
    end select
  
  case default
    write(lout,"(a)") "HDF5> ERROR: Unrecognised statement '"//lnSplit(1)//"'."
    call prror(-1)
  
  end select
  
end subroutine h5_parseInputLine

! ================================================================================================ !
end module hdf5_output

#endif
