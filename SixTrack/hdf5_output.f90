#ifdef HDF5

! ================================================================================================ !
!  HDF5 Module
! ~~~~~~~~~~~~~
!  Written by:    Veronica Berglyd Olsen, BE-ABP-HSS, May 2018
!  Last Modified: 2018-05-08
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
  logical,          public,  save :: h5_isActive    ! Existence of the HDF5 block
  logical,          public,  save :: h5_debugOn     ! HDF5 debug flag present
  logical,          public,  save :: h5_isReady     ! HDF5 file is open and ready for input
  logical,          private, save :: h5_useDouble   ! Whether to use double precision or not
  logical,          private, save :: h5_doTruncate  ! Whether or not to truncate previous file if it exists
  type(string),     private, save :: h5_fileName    ! The HDF5 output file name
  type(string),     private, save :: h5_rootPath    ! The root group where the data for this session is stored
  integer,          private, save :: h5_gzipLevel   ! The level of compression used: 0 for none to 9 for maximum
  integer(HSIZE_T), private, save :: h5_defChunk    ! The default size of chunks, used for mainly logging output
  
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
  
  ! Default Group Names
  character(len=11), parameter :: h5_collGroup = "collimation"
  character(len=4),  parameter :: h5_dumpGroup = "dump"
  character(len=7),  parameter :: h5_scatGroup = "scatter"
  
  ! Atomic Data Type Constants
  integer, parameter :: h5_typeInt  = 1 ! Size is fixed to fortran native integer
  integer, parameter :: h5_typeReal = 2 ! Size is determined by fPrec
  integer, parameter :: h5_typeChar = 3 ! Any size is valid, must be specified
  
  ! Field Definitions
  type, public :: h5_dataField
    character(len=:), allocatable, public  :: name
    integer,                       public  :: type   = 0
    integer(HID_T),                public  :: size   = 0
    integer(HID_T),                private :: typeH5 = 0
    integer(HID_T),                private :: typeID = 0
  end type h5_dataField
  
  ! Format Container
  type, public :: h5_dataFmt
    character(len=:),   allocatable, public  :: name
    integer(HSIZE_T),                private :: memSize = 0
    integer(HID_T),                  private :: dtypeID = 0
    type(h5_dataField), allocatable, public  :: fields(:)
  end type h5_dataFmt
  
  ! Dataset Container
  type, public :: h5_dataSet
    character(len=:),   allocatable, public  :: name
    character(len=:),   allocatable, public  :: path
    integer,                         private :: format  = 0
    integer(HSIZE_T),                private :: records = 0
    integer(HID_T),                  private :: groupID = 0
    integer(HID_T),                  private :: dataID  = 0
    integer(HID_T),                  private :: spaceID = 0
    integer(HID_T),                  private :: memID   = 0
    integer(HID_T),                  private :: propID  = 0
  end type h5_dataSet
  
  ! Storage Arrays
  type(h5_dataFmt), allocatable, private, save :: h5_fmtList(:)
  integer,                       private, save :: h5_fmtCount
  type(h5_dataSet), allocatable, private, save :: h5_setList(:)
  integer,                       private, save :: h5_setCount
  
  ! Interface for Data Writing
  interface h5_writeData
    module procedure h5_writeArray_real32
    module procedure h5_writeValue_real32
    module procedure h5_writeArray_real64
    module procedure h5_writeValue_real64
    module procedure h5_writeArray_real128
    module procedure h5_writeValue_real128
    module procedure h5_writeArray_int
    module procedure h5_writeValue_int
    module procedure h5_writeArray_char
    module procedure h5_writeValue_char
  end interface h5_writeData
  
  private :: h5_writeArray_real32
  private :: h5_writeValue_real32
  private :: h5_writeArray_real64
  private :: h5_writeValue_real64
  private :: h5_writeArray_real128
  private :: h5_writeValue_real128
  private :: h5_writeArray_int
  private :: h5_writeValue_int
  private :: h5_writeArray_char
  private :: h5_writeValue_char
  
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
  h5_gzipLevel  = -1
  h5_defChunk   = 10
  
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
  
  h5_fmtCount   = 0
  h5_setCount   = 0
  
end subroutine h5_comnul

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
  
  case("GZIP")
    if(nSplit /= 2) then
      write(lout,"(a,i2,a)") "HDF5> ERROR GZIP level takes 1 input parameter, ",(nSplit-1)," given."
      call prror(-1)
    end if
    read(lnSplit(2)%chr,*) h5_gzipLevel
    if(h5_gzipLevel < -1 .or. h5_gzipLevel > 9) then
      write(lout,"(a,i2)") "HDF5> ERROR Illegal value for GZIP: ",h5_gzipLevel
      write(lout,"(a,i2)") "HDF5> ERROR   Allowed values are -1 for disabled, and 0-9 for none to max compression."
      call prror(-1)
    end if
  
  case("CHUNK")
    if(nSplit /= 2) then
      write(lout,"(a,i2,a)") "HDF5> ERROR CHUNK takes 1 input parameter, ",(nSplit-1)," given."
      call prror(-1)
    end if
    read(lnSplit(2)%chr,*) h5_defChunk
    if(h5_defChunk < 1) then
      write(lout,"(a,i2)") "HDF5> ERROR Illegal value for CHUNK: ",h5_gzipLevel
      write(lout,"(a,i2)") "HDF5> ERROR   Value must be larger than 0."
      call prror(-1)
    end if
  
  case("FILE")
    if(nSplit < 2 .or. nSplit > 3) then
      write(lout,"(a,i2,a)") "HDF5> ERROR FILE statement takes 1 or 2 input parameters, ",(nSplit-1)," given."
      write(lout,"(a)")      "HDF5> ERROR   Valid input is FILE filename [truncate]"
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
      write(lout,"(a,i2,a)") "HDF5> ERROR ROOT statement takes 1 input parameter, ",(nSplit-1)," given."
      call prror(-1)
    end if
    if(str_inStr(lnSplit(2)," ") /= 0) then
      write(lout,"(a)") "HDF5> ERROR ROOT group name cannot contain a space."
      call prror(-1)
    end if
    if(str_inStr(lnSplit(2),"/") /= 0) then
      write(lout,"(a)") "HDF5> ERROR ROOT group name cannot contain a slash."
      call prror(-1)
    end if
    h5_rootPath = str_stripQuotes(lnSplit(2))
    write(lout, "(a)") "HDF5> Root group set to: '"//h5_rootPath//"'."
  
  case("ENABLE")
  
    if(nSplit /= 2) then
      write(lout,"(a,i2,a)") "HDF5> ERROR ENABLE statement takes 1 input parameter, ",(nSplit-1)," given."
      call prror(-1)
    end if
    if(len(lnSplit(2)%chr) < 4) then
      write(lout,"(a,i2,a)") "HDF5> ERROR ENABLE argument must be at least 4 characters."
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
      write(lout,"(a)") "HDF5> ERROR HDF5 output is not available for "//lnSplit(2)%chr(1:4)//" blocks."
      call prror(-1)
    end select
  
  case default
    write(lout,"(a)") "HDF5> ERROR Unrecognised statement '"//lnSplit(1)//"'."
    call prror(-1)
  
  end select
  
end subroutine h5_parseInputLine

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
  if(h5_fileError < 0) then
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
  logical doesExist
  
  if(.not. h5_isActive) return
  
  inquire(file=h5_fileName%chr, exist=doesExist)
  
  if(doesExist) then
    if(h5_doTruncate) then
      call h5fcreate_f(h5_fileName%chr, H5F_ACC_TRUNC_F, h5_fileID, h5_fileError)
      if(h5_fileError < 0) then
        write(lout,"(a)") "HDF5> ERROR Failed to open HDF5 file '"//h5_fileName//"'."
        call prror(-1)
      end if
      write(lout,"(a)") "HDF5> Truncated HDF5 file '"//h5_fileName//"'."
    else
      call h5fopen_f(h5_fileName%chr, H5F_ACC_RDWR_F, h5_fileID, h5_fileError)
      if(h5_fileError < 0) then
        write(lout,"(3a)") "HDF5> ERROR Failed to open HDF5 file '",h5_fileName%chr,"'."
        call prror(-1)
      end if
      write(lout,"(a)") "HDF5> Opened HDF5 file '"//h5_fileName//"'."
    end if
  else
    call h5fcreate_f(h5_fileName%chr, H5F_ACC_EXCL_F, h5_fileID, h5_fileError)
    if(h5_fileError < 0) then
      write(lout,"(a)") "HDF5> ERROR Failed to create HDF5 file '"//h5_fileName//"'."
      call prror(-1)
    end if
    write(lout,"(a)") "HDF5> Created HDF5 file '"//h5_fileName//"'."
  end if
  
  ! If a root group was requested, create it and save the rootID
  ! Otherwise, use the fileID as the rootID
  if(h5_rootPath%chr /= "") then
    call h5gcreate_f(h5_fileID, h5_rootPath%chr, h5_rootID, h5_fileError)
    if(h5_fileError < 0) then
      write(lout,"(a)") "HDF5> ERROR Failed to create root group '"//h5_rootPath//"'."
      write(lout,"(a)") "HDF5> ERROR   If you are writing to an existing file, the root group must be unique."
      call prror(-1)
    end if
    write(lout,"(a)") "HDF5> Created root group '"//h5_rootPath//"'."
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
  
  integer i,j
  
  ! Close all DataType IDs
  do i=1, h5_fmtCount
    do j=1, size(h5_fmtList(i)%fields)
      if(h5_fmtList(i)%fields(j)%typeID  /= 0) call h5tclose_f(h5_fmtList(i)%fields(j)%typeID, h5_dataError)
    end do
    if(h5_fmtList(i)%dtypeID /= 0) call h5tclose_f(h5_fmtList(i)%dtypeID, h5_dataError)
  end do
  
  ! Close all DataSet IDs
  do i=1, h5_setCount
    if(h5_setList(i)%dataID  /= 0) call h5dclose_f(h5_setList(i)%dataID,  h5_dataError)
    if(h5_setList(i)%spaceID /= 0) call h5sclose_f(h5_setList(i)%spaceID, h5_dataError)
    if(h5_setList(i)%memID   /= 0) call h5sclose_f(h5_setList(i)%memID,   h5_dataError)
    if(h5_setList(i)%propID  /= 0) call h5pclose_f(h5_setList(i)%propID,  h5_dataError)
  end do
  
  ! Close groups, if opened
  if(h5_rootID /= 0) call h5gclose_f(h5_rootID, h5_fileError)
  if(h5_collID /= 0) call h5gclose_f(h5_collID, h5_fileError)
  if(h5_dumpID /= 0) call h5gclose_f(h5_dumpID, h5_fileError)
  if(h5_scatID /= 0) call h5gclose_f(h5_scatID, h5_fileError)
  
  ! This closes the file
  call h5fclose_f(h5_fileID, h5_fileError)
  if(h5_fileError < 0) then
    write(lout,"(a)") "HDF5> ERROR Failed to close HDF5 file."
    call prror(-1)
  end if
  
  ! This cleans up everything left over
  call h5close_f(h5_fileError)
  if(h5_fileError < 0) then
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
  if(h5_fileError < 0) then
    write(lout,"(a)") "HDF5> ERROR Failed to create scatter group '"//h5_scatGroup//"'."
    call prror(-1)
  end if
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Group created for SCATTER."
  end if
  
end subroutine h5_initForScatter

subroutine h5_initForDump()
  
  call h5gcreate_f(h5_rootID, h5_dumpGroup, h5_dumpID, h5_fileError)
  if(h5_fileError < 0) then
    write(lout,"(a)") "HDF5> ERROR Failed to create dump group '"//h5_dumpGroup//"'."
    call prror(-1)
  end if
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Group created for DUMP."
  end if
  
end subroutine h5_initForDump

subroutine h5_initForCollimation()
  
  call h5gcreate_f(h5_rootID, h5_collGroup, h5_collID, h5_fileError)
  if(h5_fileError < 0) then
    write(lout,"(a)") "HDF5> ERROR Failed to create collimation group '"//h5_collGroup//"'."
    call prror(-1)
  end if
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Group created for COLLIMATION."
  end if
  
end subroutine h5_initForCollimation

! ================================================================================================ !
!  Create DataType
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
! ================================================================================================ !
subroutine h5_createFormat(formatName, setFields, formatID)
  
  use end_sixtrack
  
  character(len=*),                intent(in)    :: formatName
  type(h5_dataField), allocatable, intent(inout) :: setFields(:)
  integer,                         intent(out)   :: formatID
  
  type(h5_dataFmt), allocatable :: tmpTypes(:)
  integer(HID_T),   allocatable :: fieldType(:)
  integer(HSIZE_T), allocatable :: fieldSize(:)
  
  integer(HID_T)   :: dtypeID, tmpID
  integer(HSIZE_T) :: memSize, memOffset
  integer          :: i, nFields
  
  ! Check inputs
  if(len(formatName) == 0) then
    write(lout,"(a)") "HDF5> ERROR Empty format name given. This is a bug."
    call prror(-1)
  end if
  if(allocated(setFields) .eqv. .false.) then
    write(lout,"(a)") "HDF5> ERROR Fields array not allocated. This is a bug."
    call prror(-1)
  end if
  if(size(setFields) == 0) then
    write(lout,"(a)") "HDF5> ERROR Fields array empty. This is a bug."
    call prror(-1)
  end if
  
  ! First, extend the h5_fmtList array
  if(allocated(h5_fmtList) .eqv. .false.) then
    allocate(h5_fmtList(1))
    h5_fmtCount = 1
  else
    allocate(tmpTypes(h5_fmtCount + 1))
    tmpTypes(1:h5_fmtCount) = h5_fmtList(1:h5_fmtCount)
    h5_fmtCount = h5_fmtCount + 1
    call move_alloc(tmpTypes, h5_fmtList)
  end if
  
  ! Translate DataTypes into HDF5 equivalents
  memSize   = 0
  memOffset = 0
  nFields   = size(setFields)
  allocate(fieldType(nFields))
  allocate(fieldSize(nFields))
  
  do i=1,nFields
    
    select case (setFields(i)%type)
    
    case(h5_typeInt) ! Integer types
      fieldType(i) = H5T_NATIVE_INTEGER
      call h5tget_size_f(fieldType(i), fieldSize(i), h5_dataError)
    
    case(h5_typeReal) ! Real types
      if(h5_useDouble) then
        fieldType(i) = H5T_NATIVE_DOUBLE
      else
        fieldType(i) = H5T_NATIVE_REAL
      end if
      call h5tget_size_f(fieldType(i), fieldSize(i), h5_dataError)
    
    case(h5_typeChar) ! Character types
      call h5tcopy_f(H5T_NATIVE_CHARACTER, tmpID, h5_dataError)
      call h5tset_size_f(tmpID, setFields(i)%size, h5_dataError)
      call h5tget_size_f(tmpID, fieldSize(i), h5_dataError)
      fieldType(i) = tmpID
    
    end select
    
    ! Calculate the total memory allocation needed per record
    memSize = memSize + fieldSize(i)
  end do
  
  if(h5_debugOn) then
    write(lout,"(a,i0)")   "HDF5> DEBUG Creating data format '"//formatName//"' with ID ",h5_fmtCount
    write(lout,"(a,i0,a)") "HDF5> DEBUG   Format requires ",memSize," bytes per record."
  end if
  
  ! Create the compound DataType as laid out in setFields
  call h5tcreate_f(H5T_COMPOUND_F, memSize, dtypeID, h5_dataError)
  
  do i=1,nFields
    
    ! Create in-file map
    call h5tinsert_f(dtypeID, setFields(i)%name, memOffset, fieldType(i), h5_dataError)
    memOffset = memOffset + fieldSize(i)
    
    ! Create in-memory map
    setFields(i)%typeH5 = fieldType(i)
    call h5tcreate_f(H5T_COMPOUND_F, fieldSize(i), setFields(i)%typeID, h5_dataError)
    call h5tinsert_f(setFields(i)%typeID, setFields(i)%name, 0_HSIZE_T, setFields(i)%typeH5, h5_dataError)
    
  end do
  
  if(h5_dataError /= 0) then
    write(lout,"(a)") "HDF5> ERROR Failed to create compund data record '"//formatName//"'"
    call prror(-1)
  end if
  
  ! Save the Resulting HDF5 DataType
  h5_fmtList(h5_fmtCount)%name    = formatName
  h5_fmtList(h5_fmtCount)%memSize = memSize
  h5_fmtList(h5_fmtCount)%dtypeID = dtypeID
  h5_fmtList(h5_fmtCount)%fields  = setFields
  
  formatID = h5_fmtCount
  
  ! Clean up
  deallocate(fieldType)
  deallocate(fieldSize)
  
end subroutine h5_createFormat

! ================================================================================================ !
!  Create DataSet
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-30
! ================================================================================================ !
subroutine h5_createDataSet(setName, groupID, formatID, setID, chunckSize)
  
  use mod_alloc
  use end_sixtrack
  use string_tools
  
  implicit none
  
  ! Routine Variables
  character(len=*),           intent(in)  :: setName
  integer(HID_T),             intent(in)  :: groupID
  integer,                    intent(in)  :: formatID
  integer,                    intent(out) :: setID
  integer,          optional, intent(in)  :: chunckSize
  
  type(h5_dataSet), allocatable :: tmpSets(:)
  character(len=:), allocatable :: cleanName
  
  integer(HSIZE_T) :: spaceSize(1)
  integer(HSIZE_T) :: spaceMaxSize(1)
  integer(HID_T)   :: spaceID, dtypeID, dataID, propID
  
  ! First, extend the h5_setList array
  if(allocated(h5_setList) .eqv. .false.) then
    allocate(h5_setList(1))
    h5_setCount = 1
  else
    allocate(tmpSets(h5_setCount + 1))
    tmpSets(1:h5_setCount) = h5_setList(1:h5_setCount)
    h5_setCount = h5_setCount + 1
    call move_alloc(tmpSets, h5_setList)
  end if
  
  ! Determine the chuncking
  if(present(chunckSize)) then
    spaceSize(1) = int(chunckSize,kind=HSIZE_T)
  else
    spaceSize(1) = h5_defChunk
  end if
  spaceMaxSize(1) = H5S_UNLIMITED_F
  cleanName       = chr_strip(chr_trimZero(setName))
  
  if(h5_debugOn) then
    write(lout,"(a,i0)") "HDF5> DEBUG Creating dataset '"//cleanName//"' with ID ",h5_setCount
  end if
  
  ! Create a 1 by chunckSize DataSpace with a 1 by inf max size
  call h5screate_simple_f(1, spaceSize, spaceID, h5_dataError, spaceMaxSize)
  
  ! Create the chunkcing property
  call h5pcreate_f(H5P_DATASET_CREATE_F, propID, h5_dataError)
  call h5pset_chunk_f(propID, 1, spaceSize, h5_dataError)
  if(h5_dataError /= 0) then
    write(lout,"(a)") "HDF5> ERROR Failed to set chunck size for '"//cleanName//"'"
    call prror(-1)
  end if
  if(h5_gzipLevel > -1) then
    call h5pset_deflate_f(propID, h5_gzipLevel, h5_dataError)
  end if
  
  ! Create the DataSet in-file
  dtypeID = h5_fmtList(formatID)%dtypeID
  call h5dcreate_f(groupID, cleanName, dtypeID, spaceID, dataID, h5_dataError, propID)
  if(h5_dataError /= 0) then
    write(lout,"(a)") "HDF5> ERROR Failed to create dataset '"//cleanName//"'"
    call prror(-1)
  end if
  
  call h5sclose_f(spaceID, h5_dataError)
  call h5dclose_f(dataID, h5_dataError)
  
  h5_setList(h5_setCount)%name    = cleanName
  h5_setList(h5_setCount)%path    = ""
  h5_setList(h5_setCount)%format  = formatID
  h5_setList(h5_setCount)%records = 0
  h5_setList(h5_setCount)%groupID = groupID
  h5_setList(h5_setCount)%dataID  = 0
  h5_setList(h5_setCount)%spaceID = 0
  h5_setList(h5_setCount)%memID   = 0
  h5_setList(h5_setCount)%propID  = propID
  
  setID = h5_setCount
  
end subroutine h5_createDataSet

! ================================================================================================ !
!  BEGIN WRITING TO DATASETS
! ================================================================================================ !

! ================================================================================================ !
!  Prepare to Write
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
!  Creates the necessary memory and file space for writing a appendSize chunck of data.
! ================================================================================================ !
subroutine h5_prepareWrite(setID, appendSize)
  
  use end_sixtrack
  
  integer, intent(in) :: setID
  integer, intent(in) :: appendSize
  
  integer(HID_T)   :: groupID, dataID, propID, spaceID, memID
  
  integer(HSIZE_T) :: oldSize(1)
  integer(HSIZE_T) :: addSize(1)
  integer(HSIZE_T) :: newSize(1)
  
  oldSize(1) = h5_setList(setID)%records
  addSize(1) = int(appendSize,kind=HSIZE_T)
  newSize(1) = oldSize(1) + addSize(1)
  
  groupID = h5_setList(setID)%groupID
  propID  = h5_setList(setID)%propID
  
  call h5dopen_f(groupID, h5_setList(setID)%name, dataID, h5_dataError)
  call h5dextend_f(dataID, newSize, h5_dataError)
  call h5dget_space_f(dataID, spaceID, h5_dataError)
  call h5sselect_hyperslab_f(spaceID, H5S_SELECT_SET_F, oldSize, addSize, h5_dataError) 
  call h5screate_simple_f(1, addSize, memID, h5_dataError)
  
  if(h5_dataError /= 0) then
    write(lout,"(a)") "HDF5> ERROR Failed to extend dataset '"//h5_setList(setID)%name//"'"
    call prror(-1)
  end if
  
  h5_setList(setID)%records = newSize(1)
  h5_setList(setID)%dataID  = dataID
  h5_setList(setID)%spaceID = spaceID
  h5_setList(setID)%memID   = memID
  
end subroutine h5_prepareWrite

! ================================================================================================ !
!  Finalise Writing
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
!  Properly closes the file and memory dataspace.
!  This should always be called after writng a new chunck of data to prevent memory leaks during
!  executuion. Everything is properly cleaned up and closed at the end though with closeHDF5.
! ================================================================================================ !
subroutine h5_finaliseWrite(setID)
  
  integer, intent(in) :: setID
  
  call h5dclose_f(h5_setList(setID)%dataID,  h5_dataError)
  call h5sclose_f(h5_setList(setID)%spaceID, h5_dataError)
  call h5sclose_f(h5_setList(setID)%memID,   h5_dataError)
  
  h5_setList(setID)%dataID  = 0
  h5_setList(setID)%spaceID = 0
  h5_setList(setID)%memID   = 0
  
end subroutine h5_finaliseWrite

! ================================================================================================ !
!  Interfaced Write Routines
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
!  Handles writing of the differrent datatypes
! ================================================================================================ !

! ================================================================================================ !
!  Single Precision Arrays and Single Values
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
!  Converted to double when file double precision is required.
! ================================================================================================ !
subroutine h5_writeArray_real32(setID, colID, arrSize, arrData)
  
  integer,           intent(in)  :: setID
  integer,           intent(in)  :: colID
  integer,           intent(in)  :: arrSize
  real(kind=real32), intent(in)  :: arrData(arrSize)
  
  real(kind=real64), allocatable :: arrSave(:)
  integer(HID_T)                 :: dtypeID
  integer(HSIZE_T)               :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  
  if(h5_useDouble) then
    allocate(arrSave(arrSize), source=real(arrData, kind=real64))
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSave, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSave)
  else
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
  end if
  
end subroutine h5_writeArray_real32

subroutine h5_writeValue_real32(setID, colID, arrSize, valData)
  
  integer,           intent(in)  :: setID
  integer,           intent(in)  :: colID
  integer,           intent(in)  :: arrSize
  real(kind=real32), intent(in)  :: valData
  
  real(kind=real32), allocatable :: arrSaveS(:)
  real(kind=real64), allocatable :: arrSaveD(:)
  integer(HID_T)                 :: dtypeID
  integer(HSIZE_T)               :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  
  if(h5_useDouble) then
    allocate(arrSaveD(arrSize))
    arrSaveD(:) = real(valData, kind=real64)
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSaveD, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(arrSize))
    arrSaveS(:) = real(valData, kind=real32)
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSaveS, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSaveS)
  end if
  
end subroutine h5_writeValue_real32

! ================================================================================================ !
!  Double Precision Arrays and Single Values
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
!  Converted to single when file single precision is required.
! ================================================================================================ !
subroutine h5_writeArray_real64(setID, colID, arrSize, arrData)
  
  integer,           intent(in)  :: setID
  integer,           intent(in)  :: colID
  integer,           intent(in)  :: arrSize
  real(kind=real64), intent(in)  :: arrData(arrSize)
  
  real(kind=real32), allocatable :: arrSave(:)
  integer(HID_T)                 :: dtypeID
  integer(HSIZE_T)               :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  
  if(h5_useDouble) then
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
  else
    allocate(arrSave(arrSize), source=real(arrData, kind=real32))
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSave, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSave)
  end if
  
end subroutine h5_writeArray_real64

subroutine h5_writeValue_real64(setID, colID, arrSize, valData)

  integer,           intent(in)  :: setID
  integer,           intent(in)  :: colID
  integer,           intent(in)  :: arrSize
  real(kind=real64), intent(in)  :: valData
  
  real(kind=real32), allocatable :: arrSaveS(:)
  real(kind=real64), allocatable :: arrSaveD(:)
  integer(HID_T)                 :: dtypeID
  integer(HSIZE_T)               :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  
  if(h5_useDouble) then
    allocate(arrSaveD(arrSize))
    arrSaveD(:) = real(valData, kind=real64)
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSaveD, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(arrSize))
    arrSaveS(:) = real(valData, kind=real32)
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSaveS, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSaveS)
  end if
  
end subroutine h5_writeValue_real64

! ================================================================================================ !
!  Quad Precision Arrays and Single Values
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
!  Converted to double or single depending on file precision.
! ================================================================================================ !
subroutine h5_writeArray_real128(setID, colID, arrSize, arrData)
  
  integer,            intent(in)  :: setID
  integer,            intent(in)  :: colID
  integer,            intent(in)  :: arrSize
  real(kind=real128), intent(in)  :: arrData(arrSize)
  
  real(kind=real32),  allocatable :: arrSaveS(:)
  real(kind=real64),  allocatable :: arrSaveD(:)
  integer(HID_T)                  :: dtypeID
  integer(HSIZE_T)                :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  
  if(h5_useDouble) then
    allocate(arrSaveD(arrSize), source=real(arrData, kind=real64))
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSaveD, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(arrSize), source=real(arrData, kind=real32))
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSaveS, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSaveS)
  end if
  
end subroutine h5_writeArray_real128

subroutine h5_writeValue_real128(setID, colID, arrSize, valData)

  integer,            intent(in)  :: setID
  integer,            intent(in)  :: colID
  integer,            intent(in)  :: arrSize
  real(kind=real128), intent(in)  :: valData
  
  real(kind=real32),  allocatable :: arrSaveS(:)
  real(kind=real64),  allocatable :: arrSaveD(:)
  integer(HID_T)                  :: dtypeID
  integer(HSIZE_T)                :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  
  if(h5_useDouble) then
    allocate(arrSaveD(arrSize))
    arrSaveD(:) = real(valData, kind=real64)
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSaveD, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(arrSize))
    arrSaveS(:) = real(valData, kind=real32)
    call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrSaveS, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
    deallocate(arrSaveS)
  end if
  
end subroutine h5_writeValue_real128

! ================================================================================================ !
!  Integer Arrays and Single Values
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
!  These are assumed to be the internal integer type, handled by HDF5 NATIVE_INT
! ================================================================================================ !
subroutine h5_writeArray_int(setID, colID, arrSize, arrData)
  
  integer,           intent(in)  :: setID
  integer,           intent(in)  :: colID
  integer,           intent(in)  :: arrSize
  integer,           intent(in)  :: arrData(arrSize)
  
  integer(HID_T)                 :: dtypeID
  integer(HSIZE_T)               :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  
  call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
    xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
  
end subroutine h5_writeArray_int

subroutine h5_writeValue_int(setID, colID, arrSize, valData)
  
  integer,           intent(in)  :: setID
  integer,           intent(in)  :: colID
  integer,           intent(in)  :: arrSize
  integer,           intent(in)  :: valData
  
  integer                        :: arrData(arrSize)
  integer(HID_T)                 :: dtypeID
  integer(HSIZE_T)               :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  arrData(:)    = valData
  
  call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
    xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
  
end subroutine h5_writeValue_int

! ================================================================================================ !
!  Character Arrays
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
!  Handled as fixed length character arrays
! ================================================================================================ !
subroutine h5_writeArray_char(setID, colID, arrSize, arrData)
  
  integer,           intent(in)  :: setID
  integer,           intent(in)  :: colID
  integer,           intent(in)  :: arrSize
  character(len=*),  intent(in)  :: arrData(arrSize)
  
  integer(HID_T)                 :: dtypeID
  integer(HSIZE_T)               :: h5_arrSize(1)
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  
  call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
    xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
  
end subroutine h5_writeArray_char

subroutine h5_writeValue_char(setID, colID, arrSize, valData)
  
  integer,           intent(in)  :: setID
  integer,           intent(in)  :: colID
  integer,           intent(in)  :: arrSize
  character(len=*),  intent(in)  :: valData
  
  character(len=:),  allocatable :: arrData(:)
  integer(HID_T)                 :: dtypeID
  integer(HSIZE_T)               :: h5_arrSize(1)
  
  allocate(character(len=len(valData)) :: arrData(arrSize))
  
  h5_arrSize(1) = arrSize
  dtypeID       = h5_fmtList(h5_setList(setID)%format)%fields(colID)%typeID
  arrData(:)    = valData
  
  call h5dwrite_f(h5_setList(setID)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
    xfer_prp=h5_plistID, mem_space_id=h5_setList(setID)%memID, file_space_id=h5_setList(setID)%spaceID)
  
end subroutine h5_writeValue_char

! ================================================================================================ !
!  END WRITING TO DATASETS
! ================================================================================================ !

! ================================================================================================ !
end module hdf5_output

#endif
