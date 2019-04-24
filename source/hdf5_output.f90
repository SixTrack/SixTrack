#ifdef HDF5

! ================================================================================================ !
!  HDF5 Module
! ~~~~~~~~~~~~~
!  Written by:    Veronica Berglyd Olsen, BE-ABP-HSS, May 2018
!  Last Modified: 2018-05-08
!
!  Alternative output format.
!  Controlled by the HDF5 input block.
!
!   Usage:
!  ~~~~~~~~
!   Format:  A format is created from an array of type(h5_dataField). Each entry represents a column
!            with a name and a type, where the type is either h5_typeInt, h5_typeReal or
!            h5_typeChar. In the latter case, a size has to be set as well.
!            The format is created with the h5_createFormat routine, which returns a formatID.
!   DataSet: A dataset is created from a formatID, and returns a setID. Many datasets can use the
!            same formatID, but each dataset is unique.
!   Writing: Before writing, the number of lines one intend to write needs to be declared. This is
!            because a memory map of the data structure needs to be set up. After this, each data
!            type can be written to this map either as full arrays (fastest) or a single value in
!            the case of 1 record being declared. This is done by the h5_prepareWrite and
!            h5_writeData routines. After all the columns have been written, the h5_finaliseWrite
!            routine needs to be called to flush the memory map to file.
!   Buffers: In the case where the size of the dataset is not known before a loop, the option is to
!            either buffer the data locally, or write sets of single records. The latter has a small
!            overhead that may be significant if a large number of lines are written. There are
!            buffering routines in the module tht can handle this at this end, but they seem to have
!            some issues when buffering character arrays, at least for the gfortran 7.3 compiler.
!            The buffers seem to work just fine for integers and reals though.
! ================================================================================================ !
module hdf5_output

  use hdf5
  use floatPrecision
  use crcoall
  use strings
  use mod_alloc

  use, intrinsic :: iso_fortran_env, only : real32, real64, real128

  implicit none

  ! Common Settings
  logical,          public,  save :: h5_isActive   = .false. ! Existence of the HDF5 block
  logical,          public,  save :: h5_debugOn    = .false. ! HDF5 debug flag present
  logical,          public,  save :: h5_isReady    = .false. ! HDF5 file is open and ready for input
  logical,          private, save :: h5_useDouble  = .true.  ! Whether to use double precision or not
  logical,          private, save :: h5_doTruncate = .false. ! Whether or not to truncate previous file if it exists
  integer,          private, save :: h5_simNumber  = 1       ! A simulation number for the current simulation
  integer,          private, save :: h5_gzipLevel  = -1      ! The level of compression used: 0 for none to 9 for maximum
  integer(HSIZE_T), private, save :: h5_defChunk   = 10      ! The default size of chunks, used for mainly logging output
  type(string),     private, save :: h5_fileName             ! The HDF5 output file name
  type(string),     private, save :: h5_rootPath             ! The root group where the data for this session is stored

  ! Input Block Switches
  logical, public,  save :: h5_useForAPER = .false.
  logical, public,  save :: h5_useForCOLL = .false.
  logical, public,  save :: h5_useForDUMP = .false.
  logical, public,  save :: h5_useForSCAT = .false.

  ! Additional Write Flags
  logical, public,  save :: h5_writeOptics  = .false. ! Write the linear optics parameters
  logical, public,  save :: h5_writeTracks2 = .false. ! For backwards compatibility for old tracks2 format

  ! Runtime Variables
  logical, private, save :: h5_fileIsOpen = .false. ! True if file is open.
  integer, public,  save :: h5_fileError  = 0       ! For errors related to file essentials (critical)
  integer, public,  save :: h5_dataError  = 0       ! For errors related to datasets

  ! HDF5 File/Group IDs
  integer(HID_T), public,  save :: h5_fileID = 0 ! The internal ID of the file
  integer(HID_T), public,  save :: h5_rootID = 0 ! The internal ID of the root group
  integer(HID_T), public,  save :: h5_aperID = 0 ! The internal ID of the aperture group
  integer(HID_T), public,  save :: h5_collID = 0 ! The internal ID of the collimation group
  integer(HID_T), public,  save :: h5_dumpID = 0 ! The internal ID of the dump group
  integer(HID_T), public,  save :: h5_scatID = 0 ! The internal ID of the scatter group

  ! HDF5 Internals
  integer(HID_T), private, save :: h5_plistID = 0 ! Dataset transfer property

  ! Default Group Names
  character(len=8),  parameter :: h5_aperGroup = "aperture"
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
  type, private :: h5_dataFmt
    character(len=:),   allocatable, public  :: name
    integer(HSIZE_T),                private :: memSize = 0
    integer(HID_T),                  private :: dtypeID = 0
    type(h5_dataField), allocatable, public  :: fields(:)
  end type h5_dataFmt

  ! Dataset Container
  type, private :: h5_dataSet
    character(len=:),   allocatable, public  :: name
    character(len=:),   allocatable, public  :: path
    integer,                         private :: format  = 0
    integer,                         private :: buffer  = 0
    integer(HSIZE_T),                private :: records = 0
    integer(HID_T),                  private :: groupID = 0
    integer(HID_T),                  private :: dataID  = 0
    integer(HID_T),                  private :: spaceID = 0
    integer(HID_T),                  private :: memID   = 0
    integer(HID_T),                  private :: propID  = 0
  end type h5_dataSet

  ! Data Buffer Container
  type, private :: h5_dataBuf
    character(len=:), allocatable, public  :: name
    integer,                       public  :: nInt  = 0
    integer,                       public  :: nReal = 0
    integer,                       public  :: nChar = 0
    integer,                       public  :: nCols = 0
    integer,                       public  :: nRows = 0
    integer,                       public  :: bSize = 0
    integer,                       public  :: cSize = 0
    integer,          allocatable, public  :: colMap(:,:)
    integer,          allocatable, private :: iData(:,:)
    real(kind=fPrec), allocatable, private :: rData(:,:)
    character(len=:), allocatable, private :: cData(:,:)
  end type h5_dataBuf

  ! Storage Arrays
  type(h5_dataFmt), allocatable, private, save :: h5_fmtList(:)
  integer,                       private, save :: h5_fmtCount = 0
  integer,          parameter,   private       :: h5_fmtOff   = 1000
  type(h5_dataSet), allocatable, private, save :: h5_setList(:)
  integer,                       private, save :: h5_setCount = 0
  integer,          parameter,   private       :: h5_setOff   = 2000
  type(h5_dataBuf), allocatable, private, save :: h5_bufList(:)
  integer,                       private, save :: h5_bufCount = 0
  integer,          parameter,   private       :: h5_bufOff   = 3000

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

  interface h5_writeBuffer
    module procedure h5_writeBuffer_real
    module procedure h5_writeBuffer_int
    module procedure h5_writeBuffer_char
  end interface h5_writeBuffer

  interface h5_writeAttr
    module procedure h5_writeAttr_real32
    module procedure h5_writeAttr_real32_arr
    module procedure h5_writeAttr_real64
    module procedure h5_writeAttr_real64_arr
    module procedure h5_writeAttr_real128
    module procedure h5_writeAttr_real128_arr
    module procedure h5_writeAttr_int
    module procedure h5_writeAttr_int_arr
    module procedure h5_writeAttr_char
    module procedure h5_writeAttr_char_arr
  end interface h5_writeAttr

  interface h5_writeDataSetAttr
    module procedure h5_writeDataSetAttr_real32
    module procedure h5_writeDataSetAttr_real32_arr
    module procedure h5_writeDataSetAttr_real64
    module procedure h5_writeDataSetAttr_real64_arr
    module procedure h5_writeDataSetAttr_real128
    module procedure h5_writeDataSetAttr_real128_arr
    module procedure h5_writeDataSetAttr_int
    module procedure h5_writeDataSetAttr_int_arr
    module procedure h5_writeDataSetAttr_char
    module procedure h5_writeDataSetAttr_char_arr
  end interface h5_writeDataSetAttr

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

  private :: h5_writeBuffer_real
  private :: h5_writeBuffer_int
  private :: h5_writeBuffer_char

  private :: h5_writeAttr_real32
  private :: h5_writeAttr_real32_arr
  private :: h5_writeAttr_real64
  private :: h5_writeAttr_real64_arr
  private :: h5_writeAttr_real128
  private :: h5_writeAttr_real128_arr
  private :: h5_writeAttr_int
  private :: h5_writeAttr_int_arr
  private :: h5_writeAttr_char
  private :: h5_writeAttr_char_arr

  private :: h5_writeDataSetAttr_real32
  private :: h5_writeDataSetAttr_real32_arr
  private :: h5_writeDataSetAttr_real64
  private :: h5_writeDataSetAttr_real64_arr
  private :: h5_writeDataSetAttr_real128
  private :: h5_writeDataSetAttr_real128_arr
  private :: h5_writeDataSetAttr_int
  private :: h5_writeDataSetAttr_int_arr
  private :: h5_writeDataSetAttr_char
  private :: h5_writeDataSetAttr_char_arr

contains

! ================================================================================================ !
!  HDF5 Input File Parsing
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-13
! ================================================================================================ !
subroutine h5_parseInputLine(inLine,iErr)

  use string_tools

  type(string), intent(in)    :: inLine
  logical,      intent(inout) :: iErr

  type(string), allocatable :: lnSplit(:)
  integer i, nSplit
  logical spErr

  ! Split the input line
  call str_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lerr,"(a)") "HDF5> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) then
    if(h5_debugOn) then
      write (lout,"(a,i0,a)") "HDF5> DEBUG Input line len=",len_trim(inLine),": '"//inLine%strip()//"'."
      write (lout,"(a)")      "HDF5> DEBUG  * No fields found."
    end if
    return
  end if

  ! Report if debugging is ON
  if(h5_debugOn) then
    write (lout,"(a,i0,a)")  "HDF5> DEBUG Input line len=",len_trim(inLine),": '"//inLine%strip()//"'."
    write (lout,"(a,i3,a)") ("HDF5> DEBUG  * Field(",i,") = '"//lnSplit(i)//"'",i=1,nSplit)
  end if

  select case(lnSplit(1)%chr)

  case("DEBUG")
    h5_debugOn = .true.
    write(lout,"(a)") "HDF5> HDF5 block debugging is ON."

  case("SIMNO")
    if(nSplit /= 2) then
      write(lerr,"(a,i0,a)") "HDF5> ERROR SIMNO level takes 1 input parameter, ",(nSplit-1)," given."
      iErr = .true.
      return
    end if
    call str_cast(lnSplit(2),h5_simNumber,spErr)

  case("SINGLE")
    h5_useDouble = .false.
    write(lout,"(a)") "HDF5> HDF5 will use single precision."

  case("DOUBLE")
    h5_useDouble = .true.
    write(lout,"(a)") "HDF5> HDF5 will use double precision."

  case("GZIP")
    if(nSplit /= 2) then
      write(lerr,"(a,i0,a)") "HDF5> ERROR GZIP level takes 1 input parameter, ",(nSplit-1)," given."
      iErr = .true.
      return
    end if
    call str_cast(lnSplit(2),h5_gzipLevel,spErr)
    if(h5_gzipLevel < -1 .or. h5_gzipLevel > 9) then
      write(lerr,"(a,i0)") "HDF5> ERROR Illegal value for GZIP: ",h5_gzipLevel
      write(lerr,"(a)")    "HDF5> ERROR   Allowed values are -1 for disabled, and 0-9 for none to max compression."
      iErr = .true.
      return
    end if

  case("CHUNK")
    if(nSplit /= 2) then
      write(lerr,"(a,i0,a)") "HDF5> ERROR CHUNK takes 1 input parameter, ",(nSplit-1)," given."
      iErr = .true.
      return
    end if
    call str_cast(lnSplit(2),h5_defChunk,spErr)
    if(h5_defChunk < 1) then
      write(lerr,"(a,i0)") "HDF5> ERROR Illegal value for CHUNK: ",h5_defChunk
      write(lerr,"(a)")    "HDF5> ERROR   Value must be larger than 0."
      iErr = .true.
      return
    end if

  case("FILE")
    if(nSplit < 2 .or. nSplit > 3) then
      write(lerr,"(a,i0,a)") "HDF5> ERROR FILE statement takes 1 or 2 input parameters, ",(nSplit-1)," given."
      write(lerr,"(a)")      "HDF5> ERROR   Valid input is FILE filename [truncate]"
      iErr = .true.
      return
    end if
    if(nSplit == 3) then
      call str_cast(lnSplit(3),h5_doTruncate,spErr)
    else
      h5_doTruncate = .false.
    end if
    h5_fileName = trim(lnSplit(2))
    write(lout, "(a)") "HDF5> Output file name set to: '"//h5_fileName//"'."

  case("ROOT")
    if(nSplit /= 2) then
      write(lerr,"(a,i0,a)") "HDF5> ERROR ROOT statement takes 1 input parameter, ",(nSplit-1)," given."
      iErr = .true.
      return
    end if
    if(str_inStr(lnSplit(2)," ") /= 0) then
      write(lerr,"(a)") "HDF5> ERROR ROOT group name cannot contain a space."
      iErr = .true.
      return
    end if
    if(str_inStr(lnSplit(2),"/") /= 0) then
      write(lerr,"(a)") "HDF5> ERROR ROOT group name cannot contain a slash."
      iErr = .true.
      return
    end if
    h5_rootPath = trim(lnSplit(2))
    write(lout, "(a)") "HDF5> Root group set to: '"//h5_rootPath//"'."

  case("ENABLE")

    if(nSplit /= 2) then
      write(lerr,"(a,i0,a)") "HDF5> ERROR ENABLE statement takes 1 input parameter, ",(nSplit-1)," given."
      iErr = .true.
      return
    end if
    if(len(lnSplit(2)%chr) < 4) then
      write(lerr,"(a)") "HDF5> ERROR ENABLE argument must be at least 4 characters."
      iErr = .true.
      return
    end if

    select case(lnSplit(2)%chr(1:4))
    case("APER")
      h5_useForAPER = .true.
      write(lout,"(a)") "HDF5> HDF5 is enabled for APERTURE."
    case("COLL")
      h5_useForCOLL = .true.
      write(lout,"(a)") "HDF5> HDF5 is enabled for COLLIMATION."
    case("DUMP")
      h5_useForDUMP = .true.
      write(lout,"(a)") "HDF5> HDF5 is enabled for DUMP."
    case("SCAT")
      h5_useForSCAT = .true.
      write(lout,"(a)") "HDF5> HDF5 is enabled for SCATTER."
    case default
      write(lerr,"(a)") "HDF5> ERROR HDF5 output is not available for "//lnSplit(2)%chr(1:4)//" blocks."
      iErr = .true.
      return
    end select

  case("WRITE")

    if(nSplit /= 2) then
      write(lerr,"(a,i0,a)") "HDF5> ERROR WRITE statement takes 1 input parameter, ",(nSplit-1)," given."
      iErr = .true.
      return
    end if

    select case(lnSplit(2)%chr)
    case("OPTICS")
      h5_writeOptics = .true.
      write(lout,"(a)") "HDF5> Will write linear optics."
    case("TRACKS2")
      h5_writeTracks2 = .true.
      write(lout,"(a)") "HDF5> Will write tracks2."
    case default
      write(lerr,"(a)") "HDF5> ERROR Unrecognised WRITE option '"//lnSplit(2)//"'"
      iErr = .true.
      return
    end select

  case default
    write(lerr,"(a)") "HDF5> ERROR Unrecognised statement '"//lnSplit(1)//"'."
    iErr = .true.
    return

  end select

end subroutine h5_parseInputLine

! ================================================================================================ !
!  HDF5 Initialisation
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
! ================================================================================================ !
subroutine h5_initHDF5

  call h5open_f(h5_fileError)
  call h5pcreate_f(H5P_DATASET_XFER_F, h5_plistID, h5_fileError)
  call h5pset_preserve_f(h5_plistID, .true., h5_fileError)
  if(h5_fileError < 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to initialise Fortran HDF5."
    call prror(-1)
  end if
  write(lout,"(a)") "HDF5> Fortran HDF5 initialised."

end subroutine h5_initHDF5

! ================================================================================================ !
!  HDF5 File Initialisation
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
! ================================================================================================ !
subroutine h5_openFile

  integer accessFlag
  logical doesExist

  if(.not. h5_isActive) then
    write(lerr,"(a)") "HDF5> ERROR HDF5 file open called, but HDF5 is not active. This is a bug."
    call prror(-1)
  end if

  inquire(file=h5_fileName%chr, exist=doesExist)

  if(doesExist) then
    if(h5_doTruncate) then
      call h5fcreate_f(h5_fileName%chr, H5F_ACC_TRUNC_F, h5_fileID, h5_fileError)
      if(h5_fileError < 0) then
        write(lerr,"(a)") "HDF5> ERROR Failed to open HDF5 file '"//h5_fileName//"'."
        call prror(-1)
      end if
      write(lout,"(a)") "HDF5> Truncated HDF5 file '"//h5_fileName//"'."
    else
      call h5fopen_f(h5_fileName%chr, H5F_ACC_RDWR_F, h5_fileID, h5_fileError)
      if(h5_fileError < 0) then
        write(lerr,"(3a)") "HDF5> ERROR Failed to open HDF5 file '",h5_fileName%chr,"'."
        call prror(-1)
      end if
      write(lout,"(a)") "HDF5> Opened HDF5 file '"//h5_fileName//"'."
    end if
  else
    call h5fcreate_f(h5_fileName%chr, H5F_ACC_EXCL_F, h5_fileID, h5_fileError)
    if(h5_fileError < 0) then
      write(lerr,"(a)") "HDF5> ERROR Failed to create HDF5 file '"//h5_fileName//"'."
      call prror(-1)
    end if
    write(lout,"(a)") "HDF5> Created HDF5 file '"//h5_fileName//"'."
  end if

  ! If a root group was requested, create it and save the rootID
  ! Otherwise, use the fileID as the rootID
  if(h5_rootPath%chr /= "") then
    call h5gcreate_f(h5_fileID, h5_rootPath%chr, h5_rootID, h5_fileError)
    if(h5_fileError < 0) then
      write(lerr,"(a)") "HDF5> ERROR Failed to create root group '"//h5_rootPath//"'."
      write(lerr,"(a)") "HDF5> ERROR   If you are writing to an existing file, the root group must be unique."
      call prror(-1)
    end if
    write(lout,"(a)") "HDF5> Created root group '"//h5_rootPath//"'."
  else
    h5_rootID = h5_fileID
  end if

  h5_isReady = .true.

end subroutine h5_openFile

! ================================================================================================ !
!  Write Simulation Info
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-31
! ================================================================================================ !
subroutine h5_writeSimInfo

  use mod_common, only : napx, numl
  use mod_version
  use mod_meta

  character(len=23) timeStamp
  character(len=8)  cDate
  character(len=10) cTime

  ! TimeStamp
  call date_and_time(cDate,cTime)
  timeStamp = cDate(1:4)//"-"//cDate(5:6)//"-"//cDate(7:8)//"T"//cTime(1:2)//":"//cTime(3:4)//":"//cTime(5:10)

  ! Simulation info
  call h5_writeAttr(h5_rootID,"TimeStamp",     timeStamp)
  call h5_writeAttr(h5_rootID,"Particles",     napx*2)
  call h5_writeAttr(h5_rootID,"Turns",         numl)
  call h5_writeAttr(h5_rootID,"SimNumber",     h5_simNumber)

  ! SixTrack Version
  call h5_writeAttr(h5_rootID,"CreatedBy",     "SixTrack "//version)
  call h5_writeAttr(h5_rootID,"GitHash",       git_revision)
  call h5_writeAttr(h5_rootID,"ExecVersion",   version)
  call h5_writeAttr(h5_rootID,"ExecNumVersion",numvers)
  call h5_writeAttr(h5_rootID,"ExecReleased",  moddate)

  ! Write some additional infor to sim meta
  call meta_write("HDF5Active",   h5_isActive)
  call meta_write("HDF5DataFile", h5_fileName%chr)

end subroutine h5_writeSimInfo

! ================================================================================================ !
!  HDF5 Finalise
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
! ================================================================================================ !
subroutine h5_closeHDF5

  integer i,j

  if(.not. h5_isReady) then
    if(h5_debugOn) then
      write(lout,"(a)") "HDF5> DEBUG closeHDF5 called, but no files are open. Nothing to do."
    end if
    return
  end if

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
  if(h5_aperID /= 0) call h5gclose_f(h5_aperID, h5_fileError)
  if(h5_collID /= 0) call h5gclose_f(h5_collID, h5_fileError)
  if(h5_dumpID /= 0) call h5gclose_f(h5_dumpID, h5_fileError)
  if(h5_scatID /= 0) call h5gclose_f(h5_scatID, h5_fileError)
  if(h5_rootID /= h5_fileID) then
    if(h5_rootID /= 0) call h5gclose_f(h5_rootID, h5_fileError)
  end if

  ! This closes the file
  call h5fclose_f(h5_fileID, h5_fileError)
  if(h5_fileError < 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to close HDF5 file."
    call prror(-1)
  end if

  ! This cleans up everything left over
  call h5close_f(h5_fileError)
  if(h5_fileError < 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to close Fortran HDF5."
    call prror(-1)
  end if

  write(lout,"(a)") "HDF5> Closed HDF5 file."

end subroutine h5_closeHDF5

! ================================================================================================ !
!  Block-Wise Initialisation
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
! ================================================================================================ !
subroutine h5_initForAperture

  if(.not. h5_isReady) then
    write(lerr,"(a)") "HDF5> ERROR Block initialisation requested, but no HDF5 file is open. This is a bug."
    call prror(-1)
  end if

  call h5gcreate_f(h5_rootID, h5_aperGroup, h5_aperID, h5_fileError)
  if(h5_fileError < 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to create aperture group '"//h5_aperGroup//"'."
    call prror(-1)
  end if
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Group created for APERTURE."
  end if

end subroutine h5_initForAperture

subroutine h5_initForCollimation

  if(.not. h5_isReady) then
    write(lerr,"(a)") "HDF5> ERROR Block initialisation requested, but no HDF5 file is open. This is a bug."
    call prror(-1)
  end if

  call h5gcreate_f(h5_rootID, h5_collGroup, h5_collID, h5_fileError)
  if(h5_fileError < 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to create collimation group '"//h5_collGroup//"'."
    call prror(-1)
  end if
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Group created for COLLIMATION."
  end if

end subroutine h5_initForCollimation

subroutine h5_initForDump

  if(.not. h5_isReady) then
    write(lerr,"(a)") "HDF5> ERROR Block initialisation requested, but no HDF5 file is open. This is a bug."
    call prror(-1)
  end if

  call h5gcreate_f(h5_rootID, h5_dumpGroup, h5_dumpID, h5_fileError)
  if(h5_fileError < 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to create dump group '"//h5_dumpGroup//"'."
    call prror(-1)
  end if
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Group created for DUMP."
  end if

end subroutine h5_initForDump

subroutine h5_initForScatter

  if(.not. h5_isReady) then
    write(lerr,"(a)") "HDF5> ERROR Block initialisation requested, but no HDF5 file is open. This is a bug."
    call prror(-1)
  end if

  call h5gcreate_f(h5_rootID, h5_scatGroup, h5_scatID, h5_fileError)
  if(h5_fileError < 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to create scatter group '"//h5_scatGroup//"'."
    call prror(-1)
  end if
  if(h5_debugOn) then
    write(lout,"(a)") "HDF5> DEBUG Group created for SCATTER."
  end if

end subroutine h5_initForScatter

! ================================================================================================ !
!  Create DataType
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-07
! ================================================================================================ !
subroutine h5_createFormat(formatName, setFields, fmtID)

  character(len=*),                intent(in)    :: formatName
  type(h5_dataField), allocatable, intent(inout) :: setFields(:)
  integer,                         intent(out)   :: fmtID

  type(h5_dataFmt), allocatable :: tmpTypes(:)
  integer(HID_T),   allocatable :: fieldType(:)
  integer(HSIZE_T), allocatable :: fieldSize(:)

  integer(HID_T)   :: dtypeID, tmpID
  integer(HSIZE_T) :: memSize, memOffset
  integer          :: i, nFields

  if(.not. h5_isActive) then
    write(lerr,"(a)") "HDF5> ERROR HDF5 routine called, but HDF5 is not active. This is a bug."
    call prror(-1)
  end if

  ! Check inputs
  if(len(formatName) == 0) then
    write(lerr,"(a)") "HDF5> ERROR Empty format name given. This is a bug."
    call prror(-1)
  end if
  if(allocated(setFields) .eqv. .false.) then
    write(lerr,"(a)") "HDF5> ERROR Fields array not allocated. This is a bug."
    call prror(-1)
  end if
  if(size(setFields) == 0) then
    write(lerr,"(a)") "HDF5> ERROR Fields array empty. This is a bug."
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
    write(lout,"(a,i0)")   "HDF5> DEBUG Creating data format '"//formatName//"' with ID ",(h5_fmtCount+h5_fmtOff)
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
    write(lerr,"(a)") "HDF5> ERROR Failed to create compund data record '"//formatName//"'"
    call prror(-1)
  end if

  ! Save the Resulting HDF5 DataType
  h5_fmtList(h5_fmtCount)%name    = formatName
  h5_fmtList(h5_fmtCount)%memSize = memSize
  h5_fmtList(h5_fmtCount)%dtypeID = dtypeID
  h5_fmtList(h5_fmtCount)%fields  = setFields

  fmtID = h5_fmtCount + h5_fmtOff

  ! Clean up
  deallocate(fieldType)
  deallocate(fieldSize)

end subroutine h5_createFormat

! ================================================================================================ !
!  Create DataSet
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-30
! ================================================================================================ !
subroutine h5_createDataSet(setName, groupID, fmtID, setID, chunckSize)

  use mod_alloc
  use string_tools

  ! Routine Variables
  character(len=*),           intent(in)  :: setName
  integer(HID_T),             intent(in)  :: groupID
  integer,                    intent(in)  :: fmtID
  integer,                    intent(out) :: setID
  integer,          optional, intent(in)  :: chunckSize

  type(h5_dataSet), allocatable :: tmpSets(:)
  character(len=:), allocatable :: cleanName

  integer(HSIZE_T) :: spaceSize(1)
  integer(HSIZE_T) :: spaceMaxSize(1)
  integer(HID_T)   :: spaceID, dtypeID, dataID, propID

  if(.not. h5_isReady) then
    write(lerr,"(a)") "HDF5> ERROR Dataset creation requested, but no HDF5 file is open. This is a bug."
    call prror(-1)
  end if

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
  cleanName       = chr_strip(setName)

  if(h5_debugOn) then
    write(lout,"(a,i0)") "HDF5> DEBUG Creating dataset '"//cleanName//"' with ID ",(h5_setCount+h5_setOff)
  end if

  ! Create a 1 by chunckSize DataSpace with a 1 by inf max size
  call h5screate_simple_f(1, spaceSize, spaceID, h5_dataError, spaceMaxSize)

  ! Create the chunkcing property
  call h5pcreate_f(H5P_DATASET_CREATE_F, propID, h5_dataError)
  call h5pset_chunk_f(propID, 1, spaceSize, h5_dataError)
  if(h5_dataError /= 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to set chunck size for '"//cleanName//"'"
    call prror(-1)
  end if
  if(h5_gzipLevel > -1) then
    call h5pset_deflate_f(propID, h5_gzipLevel, h5_dataError)
  end if

  ! Create the DataSet in-file
  dtypeID = h5_fmtList(fmtID-h5_fmtOff)%dtypeID
  call h5dcreate_f(groupID, cleanName, dtypeID, spaceID, dataID, h5_dataError, propID)
  if(h5_dataError /= 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to create dataset '"//cleanName//"'"
    call prror(-1)
  end if

  call h5sclose_f(spaceID, h5_dataError)
  call h5dclose_f(dataID, h5_dataError)

  h5_setList(h5_setCount)%name    = cleanName
  h5_setList(h5_setCount)%path    = ""
  h5_setList(h5_setCount)%format  = fmtID
  h5_setList(h5_setCount)%records = 0
  h5_setList(h5_setCount)%groupID = groupID
  h5_setList(h5_setCount)%dataID  = 0
  h5_setList(h5_setCount)%spaceID = 0
  h5_setList(h5_setCount)%memID   = 0
  h5_setList(h5_setCount)%propID  = propID

  setID = h5_setCount + h5_setOff

end subroutine h5_createDataSet

! ================================================================================================ !
!  Create DataSet
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-30
! ================================================================================================ !
subroutine h5_createBuffer(bufName, fmtID, setID, bufSize)

  use mod_alloc

  ! Routine Variables
  character(len=*), intent(in)  :: bufName
  integer,          intent(in)  :: fmtID
  integer,          intent(in)  :: setID
  integer,          intent(in)  :: bufSize

  type(h5_dataBuf), allocatable :: tmpBufs(:)

  ! Temp Variables
  integer, allocatable :: colMap(:,:)
  integer              :: nInt, nReal, nChar, nCols, iCol, cSize

  integer,          allocatable :: iDataBuffer(:,:)
  real(kind=fPrec), allocatable :: rDataBuffer(:,:)
  character(len=:), allocatable :: cDataBuffer(:,:)

  if(.not. h5_isReady) then
    write(lerr,"(a)") "HDF5> ERROR Buffer creation requested, but no HDF5 file is open. This is a bug."
    call prror(-1)
  end if

  ! First, extend the h5_setList array
  if(allocated(h5_bufList) .eqv. .false.) then
    allocate(h5_bufList(1))
    h5_bufCount = 1
  else
    allocate(tmpBufs(h5_bufCount + 1))
    tmpBufs(1:h5_bufCount) = h5_bufList(1:h5_bufCount)
    h5_bufCount = h5_bufCount + 1
    call move_alloc(tmpBufs, h5_bufList)
  end if

  ! Iterate through the fields to set up the buffer
  nCols = size(h5_fmtList(fmtID-h5_fmtOff)%fields)
  allocate(colMap(2,nCols))
  colMap(:,:) = 0
  nInt        = 0
  nReal       = 0
  nChar       = 0
  cSize       = 0

  do iCol=1,nCols
    select case(h5_fmtList(fmtID-h5_fmtOff)%fields(iCol)%type)
    case(h5_typeInt)
      nInt           = nInt + 1
      colMap(1,iCol) = nInt
      colMap(2,iCol) = h5_typeInt
    case(h5_typeReal)
      nReal          = nReal + 1
      colMap(1,iCol) = nReal
      colMap(2,iCol) = h5_typeReal
    case(h5_typeChar)
      nChar          = nChar + 1
      colMap(1,iCol) = nChar
      colMap(2,iCol) = h5_typeChar
      if(h5_fmtList(fmtID-h5_fmtOff)%fields(iCol)%size > cSize) then
        cSize = h5_fmtList(fmtID-h5_fmtOff)%fields(iCol)%size
      end if
    end select
  end do

  h5_bufList(h5_bufCount)%name   = bufName
  h5_bufList(h5_bufCount)%nInt   = nInt
  h5_bufList(h5_bufCount)%nReal  = nReal
  h5_bufList(h5_bufCount)%nChar  = nChar
  h5_bufList(h5_bufCount)%nCols  = nCols
  h5_bufList(h5_bufCount)%nRows  = 1
  h5_bufList(h5_bufCount)%bSize  = bufSize
  h5_bufList(h5_bufCount)%cSize  = cSize
  h5_bufList(h5_bufCount)%colMap = colMap

  if(nInt  > 0) then
    call alloc(iDataBuffer,nInt,bufSize,0,"hdf5_iDataBuffer")
    h5_bufList(h5_bufCount)%iData = iDataBuffer
  end if
  if(nReal > 0) then
    call alloc(rDataBuffer,nReal,bufSize,0.0_fPrec,"hdf5_rDataBuffer")
    h5_bufList(h5_bufCount)%rData = rDataBuffer
  end if
  if(nChar > 0) then
    call alloc(cDataBuffer,cSize,nChar,bufSize,repeat(" ",cSize),"hdf5_cDataBuffer")
    h5_bufList(h5_bufCount)%cData = cDataBuffer
    write(lout,"(a,i0)") "HDF5-ALLOC> Size 1 = ",size(h5_bufList(h5_bufCount)%cData,1)
    write(lout,"(a,i0)") "HDF5-ALLOC> Size 2 = ",size(h5_bufList(h5_bufCount)%cData,2)
    write(lout,"(a,i0)") "HDF5-ALLOC> Size 3 = ",len(h5_bufList(h5_bufCount)%cData(1,1))
  end if

  h5_setList(setID-h5_setOff)%buffer = h5_bufCount+h5_bufOff

  if(h5_debugOn) then
    write(lout,"(a,i0,a)")    "HDF5> DEBUG Created data buffer '"//bufName//"' with ID ",(h5_bufCount+h5_bufOff)," containing:"
    write(lout,"(a,i2,a)")    "HDF5> DEBUG  ",nInt, " integer columns"
    write(lout,"(a,i2,a)")    "HDF5> DEBUG  ",nReal," real columns"
    write(lout,"(a,i2,a,i0)") "HDF5> DEBUG  ",nChar," character columns of length ",cSize
    write(lout,"(a,i0,a)")    "HDF5> DEBUG Buffer length is ",bufSize," rows"
  end if

end subroutine h5_createBuffer

! ================================================================================================ !
!  BEGIN WRITING TO BUFFER
!  Note: The buffering routines seem to work fine for integers and reals, but due to what possibly
!        is compiler related issues (testing on gfortran 7.3), the buffer for character values
!        produces garbage output to the hdf5 file.
! ================================================================================================ !

! ================================================================================================ !
!  Increment Buffer Index
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-16
!  This subroutine also checks if it is necessary to flush the buffer.
!  It must be called after all columns of a buffered dataset has been written to.
! ================================================================================================ !
subroutine h5_incrBuffer(setID)

  integer, intent(in) :: setID

  integer bufID

  bufID = h5_setList(setID-h5_setOff)%buffer
  if(h5_bufList(bufID-h5_bufOff)%nRows >= h5_bufList(bufID-h5_bufOff)%bSize) then
    call h5_flushBuffer(setID)
  else
    h5_bufList(bufID-h5_bufOff)%nRows = h5_bufList(bufID-h5_bufOff)%nRows + 1
  end if

end subroutine h5_incrBuffer

! ================================================================================================ !
!  Flush Buffer
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-16
!  This subroutine flushes the buffer for a given dataset, and resets it.
! ================================================================================================ !
subroutine h5_flushBuffer(setID)

  integer, intent(in) :: setID

  integer bufID, nCols, nRows, iCol, bufCol

  bufID = h5_setList(setID-h5_setOff)%buffer
  nCols = h5_bufList(bufID-h5_bufOff)%nCols
  nRows = h5_bufList(bufID-h5_bufOff)%nRows

  call h5_prepareWrite(setID,nRows)
  do iCol=1,nCols
    bufCol = h5_bufList(bufID-h5_bufOff)%colMap(1,iCol)
    select case(h5_bufList(bufID-h5_bufOff)%colMap(2,iCol))
    case(h5_typeInt)
      call h5_writeData(setID, iCol, nRows, h5_bufList(bufID-h5_bufOff)%iData(bufCol,1:nRows))
      h5_bufList(bufID-h5_bufOff)%iData(bufCol,1:nRows) = 0
    case(h5_typeReal)
      call h5_writeData(setID, iCol, nRows, h5_bufList(bufID-h5_bufOff)%rData(bufCol,1:nRows))
      h5_bufList(bufID-h5_bufOff)%rData(bufCol,1:nRows) = 0.0_fPrec
    case(h5_typeChar)
      call h5_writeData(setID, iCol, nRows, h5_bufList(bufID-h5_bufOff)%cData(bufCol,1:nRows))
      h5_bufList(bufID-h5_bufOff)%cData(bufCol,1:nRows) = repeat(" ",h5_bufList(bufID-h5_bufOff)%cSize)
    end select
  end do
  call h5_finaliseWrite(setID)

  h5_bufList(bufID-h5_bufOff)%nRows = 1

  if(h5_debugOn) then
    write(lout,"(a,i0,a)") "HDF5> DEBUG ",nRows," records were flushed from buffer '"//h5_bufList(bufID-h5_bufOff)%name//"'"
  end if

end subroutine h5_flushBuffer

! ================================================================================================ !
!  Interface Buffer Write Routines
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-16
! ================================================================================================ !
subroutine h5_writeBuffer_real(setID, colID, valData)

  integer,          intent(in) :: setID
  integer,          intent(in) :: colID
  real(kind=fPrec), intent(in) :: valData

  integer bufID, bufCol, nRows

  bufID  = h5_setList(setID-h5_setOff)%buffer
  nRows  = h5_bufList(bufID-h5_bufOff)%nRows
  bufCol = h5_bufList(bufID-h5_bufOff)%colMap(1,colID)

  h5_bufList(bufID-h5_bufOff)%rData(bufCol,nRows) = valData

end subroutine h5_writeBuffer_real

subroutine h5_writeBuffer_int(setID, colID, valData)

  integer, intent(in) :: setID
  integer, intent(in) :: colID
  integer, intent(in) :: valData

  integer bufID, bufCol, nRows

  bufID  = h5_setList(setID-h5_setOff)%buffer
  nRows  = h5_bufList(bufID-h5_bufOff)%nRows
  bufCol = h5_bufList(bufID-h5_bufOff)%colMap(1,colID)

  h5_bufList(bufID-h5_bufOff)%iData(bufCol,nRows) = valData

end subroutine h5_writeBuffer_int

subroutine h5_writeBuffer_char(setID, colID, valData)

  use string_tools

  integer,          intent(in) :: setID
  integer,          intent(in) :: colID
  character(len=*), intent(in) :: valData

  integer bufID, bufCol, nRows, cSize, inSize
  character(len=:), allocatable :: tmpData

  bufID  = h5_setList(setID-h5_setOff)%buffer
  nRows  = h5_bufList(bufID-h5_bufOff)%nRows
  bufCol = h5_bufList(bufID-h5_bufOff)%colMap(1,colID)
  cSize  = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%size
  inSize = len_trim(valData)

  if(inSize > cSize) then
    inSize = cSize
    write(lout,"(a,2(i0,a))") "HDF5> WARNING Trying to write a ",inSize," length string to a ",cSize," length buffer."
  end if

  allocate(character(len=cSize) :: tmpData)
  tmpData = repeat(" ",cSize)
  tmpData(1:inSize) = valData(1:inSize)
  h5_bufList(bufID-h5_bufOff)%cData(bufCol,nRows) = tmpData

  ! write(lout,"(a)")    "HDF5> Writing string: '"//tmpData//"'"
  ! write(lout,"(a)")    "HDF5> Buffer holds:   '"//h5_bufList(bufID-h5_bufOff)%cData(bufCol,nRows)//"'"
  ! write(lout,"(a,i0)") "HDF5>   Size 1: ",size(h5_bufList(bufID-h5_bufOff)%cData,1)
  ! write(lout,"(a,i0)") "HDF5>   Size 2: ",size(h5_bufList(bufID-h5_bufOff)%cData,2)
  ! write(lout,"(a,i0)") "HDF5>   Size 3: ",len(h5_bufList(bufID-h5_bufOff)%cData(1,1))

  deallocate(tmpData)

end subroutine h5_writeBuffer_char

! ================================================================================================ !
!  END WRITING TO BUFFER
! ================================================================================================ !

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

  integer, intent(in) :: setID
  integer, intent(in) :: appendSize

  integer(HID_T)   :: groupID, dataID, propID, spaceID, memID

  integer(HSIZE_T) :: oldSize(1)
  integer(HSIZE_T) :: addSize(1)
  integer(HSIZE_T) :: newSize(1)

  if(.not. h5_isReady) then
    write(lerr,"(a)") "HDF5> ERROR Dataset operation requested, but no HDF5 file is open. This is a bug."
    call prror(-1)
  end if

  oldSize(1) = h5_setList(setID-h5_setOff)%records
  addSize(1) = int(appendSize,kind=HSIZE_T)
  newSize(1) = oldSize(1) + addSize(1)

  groupID = h5_setList(setID-h5_setOff)%groupID
  propID  = h5_setList(setID-h5_setOff)%propID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5dextend_f(dataID, newSize, h5_dataError)
  call h5dget_space_f(dataID, spaceID, h5_dataError)
  call h5sselect_hyperslab_f(spaceID, H5S_SELECT_SET_F, oldSize, addSize, h5_dataError)
  call h5screate_simple_f(1, addSize, memID, h5_dataError)

  if(h5_dataError /= 0) then
    write(lerr,"(a)") "HDF5> ERROR Failed to extend dataset '"//h5_setList(setID-h5_setOff)%name//"'"
    call prror(-1)
  end if

  h5_setList(setID-h5_setOff)%records = newSize(1)
  h5_setList(setID-h5_setOff)%dataID  = dataID
  h5_setList(setID-h5_setOff)%spaceID = spaceID
  h5_setList(setID-h5_setOff)%memID   = memID

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

  if(.not. h5_isReady) then
    write(lerr,"(a)") "HDF5> ERROR Dataset operation requested, but no HDF5 file is open. This is a bug."
    call prror(-1)
  end if

  call h5dclose_f(h5_setList(setID-h5_setOff)%dataID,  h5_dataError)
  call h5sclose_f(h5_setList(setID-h5_setOff)%spaceID, h5_dataError)
  call h5sclose_f(h5_setList(setID-h5_setOff)%memID,   h5_dataError)

  h5_setList(setID-h5_setOff)%dataID  = 0
  h5_setList(setID-h5_setOff)%spaceID = 0
  h5_setList(setID-h5_setOff)%memID   = 0

end subroutine h5_finaliseWrite

! ================================================================================================ !
!  Interface Write Routines
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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID

  if(h5_useDouble) then
    allocate(arrSave(arrSize), source=real(arrData, kind=real64))
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSave, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
    deallocate(arrSave)
  else
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID

  if(h5_useDouble) then
    allocate(arrSaveD(arrSize))
    arrSaveD(:) = real(valData, kind=real64)
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSaveD, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(arrSize))
    arrSaveS(:) = real(valData, kind=real32)
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSaveS, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID

  if(h5_useDouble) then
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
  else
    allocate(arrSave(arrSize), source=real(arrData, kind=real32))
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSave, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID

  if(h5_useDouble) then
    allocate(arrSaveD(arrSize))
    arrSaveD(:) = real(valData, kind=real64)
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSaveD, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(arrSize))
    arrSaveS(:) = real(valData, kind=real32)
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSaveS, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID

  if(h5_useDouble) then
    allocate(arrSaveD(arrSize), source=real(arrData, kind=real64))
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSaveD, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(arrSize), source=real(arrData, kind=real32))
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSaveS, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID

  if(h5_useDouble) then
    allocate(arrSaveD(arrSize))
    arrSaveD(:) = real(valData, kind=real64)
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSaveD, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(arrSize))
    arrSaveS(:) = real(valData, kind=real32)
    call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrSaveS, h5_arrSize, h5_dataError, &
      xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)
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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID

  call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
    xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)

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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID
  arrData(:)    = valData

  call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
    xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)

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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID

  call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
    xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)

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
  dtypeID       = h5_fmtList(h5_setList(setID-h5_setOff)%format-h5_fmtOff)%fields(colID)%typeID
  arrData(:)    = valData

  call h5dwrite_f(h5_setList(setID-h5_setOff)%dataID, dtypeID, arrData, h5_arrSize, h5_dataError, &
    xfer_prp=h5_plistID, mem_space_id=h5_setList(setID-h5_setOff)%memID, file_space_id=h5_setList(setID-h5_setOff)%spaceID)

end subroutine h5_writeValue_char

! ================================================================================================ !
!  END WRITING TO DATASETS
! ================================================================================================ !

! ================================================================================================ !
!  Write Attribute
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-29
! ================================================================================================ !

! Real 32bit Attributes
subroutine h5_writeAttr_real32(attrTarget, attrName, attrValue)

  integer(HID_T),    intent(in)  :: attrTarget
  character(len=*),  intent(in)  :: attrName
  real(kind=real32), intent(in)  :: attrValue
  real(kind=real32), allocatable :: arrSaveS(:)
  real(kind=real64), allocatable :: arrSaveD(:)

  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrDim(1) = 1

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  if(h5_useDouble) then
    allocate(arrSaveD(1))
    arrSaveD(:) = real(attrValue, kind=real64)
    call h5tcopy_f(H5T_NATIVE_DOUBLE, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveD, attrDim, h5_dataError)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(1))
    arrSaveS(:) = real(attrValue, kind=real32)
    call h5tcopy_f(H5T_NATIVE_REAL, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveS, attrDim, h5_dataError)
    deallocate(arrSaveS)
  end if
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_real32

subroutine h5_writeAttr_real32_arr(attrTarget, attrName, attrValue)

  integer(HID_T),    intent(in)  :: attrTarget
  character(len=*),  intent(in)  :: attrName
  real(kind=real32), intent(in)  :: attrValue(:)
  real(kind=real32), allocatable :: arrSaveS(:)
  real(kind=real64), allocatable :: arrSaveD(:)

  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrDim(1) = size(attrValue,1)

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  if(h5_useDouble) then
    allocate(arrSaveD(size(attrValue,1)))
    arrSaveD(:) = real(attrValue, kind=real64)
    call h5tcopy_f(H5T_NATIVE_DOUBLE, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveD, attrDim, h5_dataError)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(size(attrValue,1)))
    arrSaveS(:) = real(attrValue, kind=real32)
    call h5tcopy_f(H5T_NATIVE_REAL, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveS, attrDim, h5_dataError)
    deallocate(arrSaveS)
  end if
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_real32_arr

! Real 64bit Attributes
subroutine h5_writeAttr_real64(attrTarget, attrName, attrValue)

  integer(HID_T),    intent(in)  :: attrTarget
  character(len=*),  intent(in)  :: attrName
  real(kind=real64), intent(in)  :: attrValue
  real(kind=real32), allocatable :: arrSaveS(:)
  real(kind=real64), allocatable :: arrSaveD(:)

  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrDim(1) = 1

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  if(h5_useDouble) then
    allocate(arrSaveD(1))
    arrSaveD(:) = real(attrValue, kind=real64)
    call h5tcopy_f(H5T_NATIVE_DOUBLE, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveD, attrDim, h5_dataError)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(1))
    arrSaveS(:) = real(attrValue, kind=real32)
    call h5tcopy_f(H5T_NATIVE_REAL, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveS, attrDim, h5_dataError)
    deallocate(arrSaveS)
  end if
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_real64

subroutine h5_writeAttr_real64_arr(attrTarget, attrName, attrValue)

  integer(HID_T),    intent(in)  :: attrTarget
  character(len=*),  intent(in)  :: attrName
  real(kind=real64), intent(in)  :: attrValue(:)
  real(kind=real32), allocatable :: arrSaveS(:)
  real(kind=real64), allocatable :: arrSaveD(:)

  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrDim(1) = size(attrValue,1)

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  if(h5_useDouble) then
    allocate(arrSaveD(size(attrValue,1)))
    arrSaveD(:) = real(attrValue, kind=real64)
    call h5tcopy_f(H5T_NATIVE_DOUBLE, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveD, attrDim, h5_dataError)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(size(attrValue,1)))
    arrSaveS(:) = real(attrValue, kind=real32)
    call h5tcopy_f(H5T_NATIVE_REAL, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveS, attrDim, h5_dataError)
    deallocate(arrSaveS)
  end if
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_real64_arr

! Real 128bit Attributes
subroutine h5_writeAttr_real128(attrTarget, attrName, attrValue)

  integer(HID_T),     intent(in)  :: attrTarget
  character(len=*),   intent(in)  :: attrName
  real(kind=real128), intent(in)  :: attrValue
  real(kind=real32),  allocatable :: arrSaveS(:)
  real(kind=real64),  allocatable :: arrSaveD(:)

  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrDim(1) = 1

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  if(h5_useDouble) then
    allocate(arrSaveD(1))
    arrSaveD(:) = real(attrValue, kind=real64)
    call h5tcopy_f(H5T_NATIVE_DOUBLE, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveD, attrDim, h5_dataError)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(1))
    arrSaveS(:) = real(attrValue, kind=real32)
    call h5tcopy_f(H5T_NATIVE_REAL, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveS, attrDim, h5_dataError)
    deallocate(arrSaveS)
  end if
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_real128

subroutine h5_writeAttr_real128_arr(attrTarget, attrName, attrValue)

  integer(HID_T),     intent(in)  :: attrTarget
  character(len=*),   intent(in)  :: attrName
  real(kind=real128), intent(in)  :: attrValue(:)
  real(kind=real32),  allocatable :: arrSaveS(:)
  real(kind=real64),  allocatable :: arrSaveD(:)

  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrDim(1) = size(attrValue,1)

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  if(h5_useDouble) then
    allocate(arrSaveD(size(attrValue,1)))
    arrSaveD(:) = real(attrValue, kind=real64)
    call h5tcopy_f(H5T_NATIVE_DOUBLE, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveD, attrDim, h5_dataError)
    deallocate(arrSaveD)
  else
    allocate(arrSaveS(size(attrValue,1)))
    arrSaveS(:) = real(attrValue, kind=real32)
    call h5tcopy_f(H5T_NATIVE_REAL, typeID, h5_dataError)
    call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
    call h5awrite_f(attrID, typeID, arrSaveS, attrDim, h5_dataError)
    deallocate(arrSaveS)
  end if
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_real128_arr

! Integer Attributes
subroutine h5_writeAttr_int(attrTarget, attrName, attrValue)

  integer(HID_T),   intent(in) :: attrTarget
  character(len=*), intent(in) :: attrName
  integer,          intent(in) :: attrValue

  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrDim(1) = 1

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  call h5tcopy_f(H5T_NATIVE_INTEGER, typeID, h5_dataError)
  call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
  call h5awrite_f(attrID, typeID, attrValue, attrDim, h5_dataError)
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_int

subroutine h5_writeAttr_int_arr(attrTarget, attrName, attrValue)

  integer(HID_T),   intent(in) :: attrTarget
  character(len=*), intent(in) :: attrName
  integer,          intent(in) :: attrValue(:)

  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrDim(1) = size(attrValue,1)

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  call h5tcopy_f(H5T_NATIVE_INTEGER, typeID, h5_dataError)
  call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
  call h5awrite_f(attrID, typeID, attrValue, attrDim, h5_dataError)
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_int_arr

! Character Attributes
subroutine h5_writeAttr_char(attrTarget, attrName, attrValue)

  integer(HID_T),   intent(in) :: attrTarget
  character(len=*), intent(in) :: attrName
  character(len=*), intent(in) :: attrValue

  integer(HSIZE_T) :: attrLen
  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrLen    = len(attrValue,kind=HSIZE_T)
  attrDim(1) = 1

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, h5_dataError)
  call h5tset_size_f(typeID, attrLen, h5_dataError)
  call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
  call h5awrite_f(attrID, typeID, attrValue, attrDim, h5_dataError)
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_char

subroutine h5_writeAttr_char_arr(attrTarget, attrName, attrValue)

  integer(HID_T),   intent(in) :: attrTarget
  character(len=*), intent(in) :: attrName
  character(len=*), intent(in) :: attrValue(:)

  integer(HSIZE_T) :: attrLen
  integer(HSIZE_T) :: attrDim(1)
  integer(HID_T)   :: spaceID, typeID, attrID

  attrLen    = len(attrValue,kind=HSIZE_T)
  attrDim(1) = size(attrValue,1)

  call h5screate_simple_f(1, attrDim, spaceID, h5_dataError)
  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, h5_dataError)
  call h5tset_size_f(typeID, attrLen, h5_dataError)
  call h5acreate_f(attrTarget, attrName, typeID, spaceID, attrID, h5_dataError)
  call h5awrite_f(attrID, typeID, attrValue, attrDim, h5_dataError)
  call h5aclose_f(attrID, h5_dataError)
  call h5sclose_f(spaceID, h5_dataError)

end subroutine h5_writeAttr_char_arr

! ================================================================================================ !
!  Write DataSet Attribute
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-09-25
! ================================================================================================ !


! Real 32bit Attributes
subroutine h5_writeDataSetAttr_real32(setID, attrName, attrValue)

  integer,           intent(in) :: setID
  character(len=*),  intent(in) :: attrName
  real(kind=real32), intent(in) :: attrValue

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_real32(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_real32

subroutine h5_writeDataSetAttr_real32_arr(setID, attrName, attrValue)

  integer,           intent(in) :: setID
  character(len=*),  intent(in) :: attrName
  real(kind=real32), intent(in) :: attrValue(:)

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_real32_arr(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_real32_arr

! Real 64bit Attributes
subroutine h5_writeDataSetAttr_real64(setID, attrName, attrValue)

  integer,           intent(in) :: setID
  character(len=*),  intent(in) :: attrName
  real(kind=real64), intent(in) :: attrValue

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_real64(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_real64

subroutine h5_writeDataSetAttr_real64_arr(setID, attrName, attrValue)

  integer,           intent(in) :: setID
  character(len=*),  intent(in) :: attrName
  real(kind=real64), intent(in) :: attrValue(:)

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_real64_arr(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_real64_arr

! Real 128bit Attributes
subroutine h5_writeDataSetAttr_real128(setID, attrName, attrValue)

  integer,            intent(in) :: setID
  character(len=*),   intent(in) :: attrName
  real(kind=real128), intent(in) :: attrValue

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_real128(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_real128

subroutine h5_writeDataSetAttr_real128_arr(setID, attrName, attrValue)

  integer,            intent(in) :: setID
  character(len=*),   intent(in) :: attrName
  real(kind=real128), intent(in) :: attrValue(:)

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_real128_arr(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_real128_arr

! Integer Attributes
subroutine h5_writeDataSetAttr_int(setID, attrName, attrValue)

  integer,          intent(in) :: setID
  character(len=*), intent(in) :: attrName
  integer,          intent(in) :: attrValue

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_int(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_int

subroutine h5_writeDataSetAttr_int_arr(setID, attrName, attrValue)

  integer,          intent(in) :: setID
  character(len=*), intent(in) :: attrName
  integer,          intent(in) :: attrValue(:)

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_int_arr(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_int_arr

! Character Attributes
subroutine h5_writeDataSetAttr_char(setID, attrName, attrValue)

  integer,          intent(in) :: setID
  character(len=*), intent(in) :: attrName
  character(len=*), intent(in) :: attrValue

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_char(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_char

subroutine h5_writeDataSetAttr_char_arr(setID, attrName, attrValue)

  integer,          intent(in) :: setID
  character(len=*), intent(in) :: attrName
  character(len=*), intent(in) :: attrValue(:)

  integer(HID_T)   :: groupID, dataID

  groupID = h5_setList(setID-h5_setOff)%groupID

  call h5dopen_f(groupID, h5_setList(setID-h5_setOff)%name, dataID, h5_dataError)
  call h5_writeAttr_char_arr(dataID, attrName, attrValue)
  call h5dclose_f(dataID, h5_dataError)

end subroutine h5_writeDataSetAttr_char_arr

! ================================================================================================ !
end module hdf5_output

#endif
