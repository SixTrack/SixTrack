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
  
  implicit none
  
  ! Common Settings
  logical,      public,  save :: h5_isActive    ! Existence of the HDF5 block
  logical,      public,  save :: h5_debugOn     ! HDF5 debug flag present
  logical,      public,  save :: h5_isReady     ! HDF5 file is open and ready for input
  logical,      private, save :: h5_doTruncate  ! Whether or not to truncate previous file if it exists
  type(string), private, save :: h5_fileName    ! The HDF5 output file name
  type(string), private, save :: h5_rootPath    ! The root group where the data for this session is stored
  
  ! Input Block Switches
  logical, public, save :: h5_useForCOLL
  logical, public, save :: h5_useForDUMP
  logical, public, save :: h5_useForSCAT
  
  ! Runtime Variables
  integer, public,  save :: h5_fileError
  logical, private, save :: h5_fileIsOpen
  
  ! HDF5 File/Group IDs
  integer(HID_T), public, save :: h5_fileID ! The internal ID of the file
  integer(HID_T), public, save :: h5_rootID ! The internal ID of the root group
  integer(HID_T), public, save :: h5_collID ! The internal ID of the collimation group
  integer(HID_T), public, save :: h5_dumpID ! The internal ID of the dump group
  integer(HID_T), public, save :: h5_scatID ! The internal ID of the scatter group
  
  ! Default Group Names
  character(len=11), parameter :: h5_collGroup = "collimation"
  character(len=4),  parameter :: h5_dumpGroup = "dump"
  character(len=7),  parameter :: h5_scatGroup = "scatter"
  
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
  h5_doTruncate = .false.
  h5_fileName   = string("")
  h5_rootPath   = string("")
  
  h5_useForCOLL = .false.
  h5_useForDUMP = .false.
  h5_useForSCAT = .false.
  
  h5_fileError  = 0
  h5_fileIsOpen = .false.
  
  h5_fileID     = 0
  h5_rootID     = 0
  h5_collID     = 0
  h5_dumpID     = 0
  h5_scatID     = 0
  
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
      write (lout,"(a,i3,a)") "HDF5> DEBUG Input line len=",len(inLine),": '"//trim(inLine)//"'."
      write (lout,"(a)")      "HDF5> DEBUG  * No fields found."
    end if
    return
  end if
  
  ! Report if debugging is ON
  if(h5_debugOn) then
    write (lout,"(a,i3,a)")  "HDF5> DEBUG Input line len=",len(inLine),": '"//trim(inLine)//"'."
    write (lout,"(a,i2,a)") ("HDF5> DEBUG  * Field(",i,") = '"//lnSplit(i)//"'",i=1,nSplit)
  end if
  
  select case(lnSplit(1)%chr)
  
  case("DEBUG")
    h5_debugOn = .true.
    write(lout,"(a)") "HDF5> HDF5 block debugging is ON."
  
  case("FILE")
    if(nSplit < 2 .or. nSplit > 3) then
      write(lout,"(a,i2,a)") "HDF5> ERROR in FILE. Statement takes 1 or 2 input parameters, ",(nSplit-1)," given."
      write(lout,"(a)")      "HDF5> Valid input is FILE filename [truncate]"
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
      write(lout,"(a,i2,a)") "HDF5> ERROR in ROOT. Statement takes 1 input parameter, ",(nSplit-1)," given."
      call prror(-1)
    end if
    if(str_inStr(lnSplit(2)," ") /= 0) then
      write(lout,"(a)") "HDF5> ERROR in ROOT. Group name cannot contain a space."
      call prror(-1)
    end if
    if(str_inStr(lnSplit(2),"/") /= 0) then
      write(lout,"(a)") "HDF5> ERROR in ROOT. Group name cannot contain a slash."
      call prror(-1)
    end if
    h5_rootPath = str_stripQuotes(lnSplit(2))
    write(lout, "(a)") "HDF5> Root group set to: '"//h5_rootPath//"'."
  
  case("ENABLE")
  
    if(nSplit /= 2) then
      write(lout,"(a,i2,a)") "HDF5> ERROR in ENABLE. Statement takes 1 input parameter, ",(nSplit-1)," given."
      call prror(-1)
    end if
    if(len(lnSplit(2)%chr) < 4) then
      write(lout,"(a,i2,a)") "HDF5> ERROR in ENABLE. Argument must be at least 4 characters."
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
      write(lout,"(a)") "HDF5> ERROR HDF5 is not available for "//lnSplit(2)%chr(1:4)//" blocks."
      call prror(-1)
    end select
  
  case default
    write(lout,"(a)") "HDF5> ERROR Unrecognised statement '"//lnSplit(1)//"'."
    call prror(-1)
  
  end select
  
end subroutine h5_parseInputLine

! ================================================================================================ !
end module hdf5_output

#endif
