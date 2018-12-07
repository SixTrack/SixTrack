! ================================================================================================ !
!  DYNAMIC KICKS
!  A.Mereghetti, for the FLUKA Team
!  K.Sjobak, A. Santamaria, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-28
! ================================================================================================ !
module dynk

  use floatPrecision
  use mathlib_bouncer
  use mod_hions
  use numerical_constants, only : zero, one, two, c1e3
  use parpro, only : nele, mNameLen, str_dSpace
  use mod_alloc
  use string_tools

  implicit none

  ! General-purpose variables
  logical, public,  save :: dynk_enabled      = .false. ! DYNK input bloc issued in the fort.3 file
  logical, public,  save :: dynk_debug        = .false. ! Print debug messages in main output
  logical, public,  save :: dynk_noDynkSets   = .false. ! Disable writing dynksets.dat?
  integer, private, save :: dynk_fileUnit               ! The file unit for dynksets.dat
  integer, private, save :: dynk_fileUnitFUN            ! File unit for parseFUN files

  ! Max Array Sizes
  integer, private, save :: dynk_maxFuncs
  integer, private, save :: dynk_maxSets
  integer, private, save :: dynk_maxiData
  integer, private, save :: dynk_maxfData
  integer, private, save :: dynk_maxcData

  ! Number of used positions in arrays
  integer, private, save :: dynk_nFuncs       = 0
  integer, private, save :: dynk_niData       = 0
  integer, private, save :: dynk_nfData       = 0
  integer, private, save :: dynk_ncData       = 0
  integer, public,  save :: dynk_nSets        = 0
  integer, public,  save :: dynk_nSets_unique = 0

  ! 1 row/FUN, cols are:
  ! (1) = function name in fort.3 (points within dynk_cData),
  ! (2) = indicates function type
  ! (3,4,5) = arguments (often pointing within other arrays {i|f|c}expr_dynk)
  integer,          allocatable, private, save :: dynk_funcs(:,:)

  ! Data for DYNK FUNs
  integer,          allocatable, private, save :: dynk_iData(:)
  real(kind=fPrec), allocatable, private, save :: dynk_fData(:)
  character(len=:), allocatable, private, save :: dynk_cData(:)

  ! 1 row/SET, cols are:
  ! (1) = function index (points within dynk_funcs)
  ! (2) = first turn num. where it is active
  ! (3) =  last turn num. where it is active
  ! (4) = Turn shift - number added to turn before evaluating the FUN
  integer,          allocatable, private, save :: dynk_sets(:,:)

  ! 1 row/SET (same ordering as dynk_sets), cols are:
  ! (1) element name
  ! (2) attribute name
  character(len=:), allocatable, public,  save :: dynk_cSets(:,:)

  ! Similar to dynk_cSets, but only one entry per elem/attr
  character(len=:), allocatable, public,  save :: dynk_cSets_unique(:,:)

  ! Some elements (multipoles) overwrites the general settings info when initialized.
  ! Store this information on the side.
  ! Also used by setvalue and getvalue
  integer,          allocatable, public,  save :: dynk_izuIndex(:)
  real(kind=fPrec), allocatable, public,  save :: dynk_elemData(:,:)

#ifdef CR
  ! Number of records written to dynkfile (dynksets.dat)
  integer, public,  save :: dynk_filePos   = -1
  integer, private, save :: dynk_filePosCR

  ! Data for DYNK FUNs
  integer,          allocatable, private, save :: dynk_iData_cr(:)
  real(kind=fPrec), allocatable, private, save :: dynk_fData_cr(:)
  character(len=:), allocatable, private, save :: dynk_cData_cr(:)

  ! Number of used positions in arrays
  integer, private, save :: dynk_niData_cr
  integer, private, save :: dynk_nfData_cr
  integer, private, save :: dynk_ncData_cr

  ! Store current settings from dynk
  real(kind=fPrec), allocatable, public, save :: dynk_fSets_cr(:)
#endif

contains

subroutine dynk_allocate_arrays
  call alloc(dynk_izuIndex,nele,0,     "dynk_izuIndex")
  call alloc(dynk_elemData,nele,3,zero,"dynk_elemData")
end subroutine dynk_allocate_arrays

subroutine dynk_expand_arrays(nele_new)
  integer, intent(in) :: nele_new
  call alloc(dynk_izuIndex,nele_new,0,     "dynk_izuIndex")
  call alloc(dynk_elemData,nele_new,3,zero,"dynk_elemData")
end subroutine dynk_expand_arrays

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-28
! =================================================================================================
subroutine dynk_allocate

  use crcoall
  use file_units

  ! Setting inital allocations
  ! These values are increased if needed when dynk_checkspace is called
  dynk_maxiData = 500
  dynk_maxfData = 500
  dynk_maxcData = 100
  dynk_maxFuncs = 10
  dynk_maxSets  = 10

  call alloc(dynk_iData,               dynk_maxiData,  0,         "dynk_iData")
  call alloc(dynk_fData,               dynk_maxfData,  zero,      "dynk_fData")
  call alloc(dynk_cData,       mStrLen,dynk_maxcData,  str_dSpace,"dynk_cData")
  call alloc(dynk_funcs,               dynk_maxFuncs,5,0,         "dynk_funcs")

  call alloc(dynk_cSets,       mStrLen,dynk_maxSets,2, str_dSpace,"dynk_cSets")
  call alloc(dynk_cSets_unique,mStrLen,dynk_maxSets,2, str_dSpace,"dynk_cSets_unique")
#ifdef CR
  call alloc(dynk_fSets_cr,            dynk_maxSets,   zero,      "dynk_fSets_cr")
#endif
  call alloc(dynk_sets,                dynk_maxSets,4, 0,         "dynk_sets")

  ! Set file units for I/O files
  call funit_requestUnit("dynksets.dat",    dynk_fileUnit)
  call funit_requestUnit("dynk_parseFUN_IO",dynk_fileUnitFUN)

end subroutine dynk_allocate

! =================================================================================================
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-28
!  Parse input line
! =================================================================================================
subroutine dynk_parseInputLine(inLine,iErr)

  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable :: lnSplit(:)
  integer i, nSplit
  logical spErr

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lout,"(a)") "DYNK> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) then
    if(dynk_debug) then
      write(lout,"(a,i0,a)") "DYNK> DEBUG Input line len = ",len_trim(inLine),": '"//trim(inLine)//"'."
      write(lout,"(a)")      "DYNK> DEBUG  * No fields found."
    end if
    return
  end if

  ! Report if debugging is ON
  if(dynk_debug) then
    write(lout,"(a,i0,a)")  "DYNK> DEBUG Input line len = ",len_trim(inLine),": '"//trim(inLine)//"'."
    write(lout,"(a,i3,a)") ("DYNK> DEBUG  * Field(",i,") = '"//trim(lnSplit(i))//"'",i=1,nSplit)
  end if

  select case(trim(lnSplit(1)))

  case("DEBUG")
    dynk_debug = .true.
    write(lout,"(a)") "DYNK> Debugging is ENABLED"

  case("NOFILE")
    dynk_noDynkSets = .true.
    write(lout,*) "DYNK> Disabled writing dynksets.dat"

  case("FUN")
    call dynk_parseFUN(inLine,iErr)

  case("SET")
    call dynk_parseSET(inLine,iErr)

  case default
    write(lout,"(a)") "DYNK> ERROR Unrecognised statement '"//trim(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

end subroutine dynk_parseInputLine

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-28
!  Parse FUN lines in the fort.3 input file.
! =================================================================================================
subroutine dynk_parseFUN(inLine, iErr)

  use crcoall
  use file_units

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable :: lnSplit(:), lnFile(:)
  character(len=mInputLn)       :: fLine
  integer nSplit, nFile
  logical spErr, cErr

  ! Temp variables
  integer          ii,t,nLines,ioStat
  real(kind=fPrec) x,y,z,u                                ! FILE, FILELIN, FIR/IIR
  real(kind=fPrec) x1,x2,y1,y2,deriv                      ! LINSEG, QUADSEG,
  real(kind=fPrec) tinj,Iinj,Inom,A,D,R,te                ! PELP (input)
  real(kind=fPrec) derivI_te,I_te,bexp,aexp,t1,I1,td,tnom ! PELP (calc)

  logical isFIR ! FIR/IIR
  logical isOpen

#ifdef BOINC
  character(len=256) filename
#endif

  if(dynk_nFuncs+1 > dynk_maxFuncs) then
    dynk_maxFuncs = dynk_maxFuncs + 10
    call alloc(dynk_funcs,dynk_maxFuncs,5,0,"dynk_funcs")
  end if

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lout,"(a)") "DYNK> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(len_trim(lnSplit(2)) > 20) then
    write(lout,"(a,i0)") "DYNK> ERROR Max length of a FUN name is 20, got ",len_trim(lnSplit(2))
    iErr = .true.
    return
  end if

  ! Parse the different functions
  !===============================

  ! System functions: #0-19
  select case(trim(lnSplit(3)))

  case("GET")
    ! GET: Store the value of an element/value

    call dynk_checkargs(nSplit,5,"FUN funname GET elementName attribute")
    call dynk_checkspace(0,1,3)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1
    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData !NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 0           !TYPE (GET)
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData !ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          !ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          !ARG3

    ! Sanity checks, length of BEZ elements
    if(len_trim(lnSplit(4)) > mNameLen) then
      write (lout,"(2(a,i0))") "DYNK> ERROR FUN:GET got an element name with length ",len_trim(lnSplit(4)),&
        ", but max is ", mNameLen
      call prror(-1)
    end if
    if(len_trim(lnSplit(5)) > mStrLen-1) then
      write (lout,"(2(a,i0))") "DYNK> ERROR FUN:GET got an attribute name with length ",len_trim(lnSplit(5)),&
        ", but max is ",mStrLen-1
      call prror(-1)
    end if

    ! Store data
    dynk_cData(dynk_ncData)   = trim(lnSplit(2)) ! NAME
    dynk_cData(dynk_ncData+1) = trim(lnSplit(4)) ! ELEMENT_NAME
    dynk_cData(dynk_ncData+2) = trim(lnSplit(5)) ! ATTRIBUTE_NAME
    dynk_ncData               = dynk_ncData+2
    dynk_fData(dynk_nfData)   = -1.0 ! Initialize a place in the array to store the value

  ! END CASE GET

  case("FILE")
    ! FILE: Load the contents from a file
    ! File format: two ASCII columns of numbers,
    ! first  column = turn number (all turns should be there, starting from 1)
    ! second column = value (as a double)

    call dynk_checkargs(nSplit,4,"FUN funname FILE filename")
    call dynk_checkspace(0,0,2)

    ! Set pointers to start of funs data blocks (dynk_nfData handled when reading data)
    dynk_nFuncs = dynk_nFuncs+1
    dynk_ncData = dynk_ncData+1
    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData   ! NAME     (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 1             ! TYPE     (FILE)
    dynk_funcs(dynk_nFuncs,3) = dynk_ncData+1 ! Filename (in dynk_cData)
    dynk_funcs(dynk_nFuncs,4) = dynk_nfData+1 ! Data     (in dynk_fData)
    dynk_funcs(dynk_nFuncs,5) = -1            ! Below: Length of file

    ! Store data
    dynk_cData(dynk_ncData)   = trim(lnSplit(2)) ! NAME
    dynk_cData(dynk_ncData+1) = trim(lnSplit(4)) ! FILE NAME
    dynk_ncData = dynk_ncData+1

    ! Open the file
    inquire(unit=dynk_fileUnitFUN,opened=isOpen)
    if(isOpen) then
      write(lout,"(a)") "DYNK> ERROR FUN:FILE ould not open file '"//trim(dynk_cData(dynk_ncData))//"'"
      iErr = .true.
      return
    end if

#ifdef BOINC
    call boincrf(dynk_cData(dynk_ncData),filename)
    open(unit=dynk_fileUnitFUN,file=filename,action="read",iostat=ioStat,status="old")
#else
    open(unit=dynk_fileUnitFUN,file=dynk_cData(dynk_ncData),action="read",iostat=ioStat,status="old")
#endif
    if(ioStat /= 0) then
      write(lout,"(a)") "DYNK> ERROR FUN:FILE ould not open file '"//trim(dynk_cData(dynk_ncData))//"'"
      iErr = .true.
      return
    end if

    ! Count number of lines and allocate space
    nLines = 0
    do
      read(dynk_fileUnitFUN,"(a)",iostat=ioStat)
      if(ioStat /= 0) exit
      nLines = nLines + 1
    end do
    rewind(dynk_fileUnitFUN)
    call dynk_checkspace(0,nLines,0)

    ii = 0 ! Number of data lines read
    do
      read(dynk_fileUnitFUN,"(a)",iostat=ioStat) fLine
      if(ioStat /= 0) exit

      call chr_split(fLine,lnFile,nFile,spErr)
      if(spErr) then
        write(lout,"(a)") "DYNK> ERROR FUN:FILE Failed to parse input line."
        iErr = .true.
        return
      end if
      if(nFile /= 2) then
        write(lout,"(a,i0)") "DYNK> ERROR FUN:FILE Expected two values per line, got ",nFile
        iErr = .true.
        return
      end if

      call chr_cast(lnFile(1),t,cErr)
      call chr_cast(lnFile(2),y,cErr)

      ii = ii+1
      if(t /= ii) then
        write(lout,"(a)")       "DYNK> ERROR FUN:FILE Reading file '"//trim(dynk_cData(dynk_ncData))//"'"
        write(lout,"(2(a,i0))") "DYNK>       Missing turn number ",ii,", got turn ",t
        iErr = .true.
        return
      end if

      dynk_nfData = dynk_nfData+1
      dynk_fData(dynk_nfData) = y
    end do
    dynk_funcs(dynk_nFuncs,5) = ii
    close(dynk_fileUnitFUN)

  ! END CASE FILE

  case("FILELIN")
    ! FILELIN: Load the contents from a file, linearly interpolate
    ! File format: two ASCII columns of numbers,
    ! first  column = turn number (as a double)
    ! second column = value (as a double)

    call dynk_checkargs(nSplit,4,"FUN funname FILELIN filename")
    call dynk_checkspace(0,0,2)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_ncData = dynk_ncData+1
    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData   !NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 2             !TYPE (FILELIN)
    dynk_funcs(dynk_nFuncs,3) = dynk_ncData+1 !Filename (in dynk_cData)
    dynk_funcs(dynk_nFuncs,4) = dynk_nfData+1 !Data     (in dynk_fData)
    dynk_funcs(dynk_nFuncs,5) = -1            !Below: Length of file (number of x,y sets)

    ! Store data
    dynk_cData(dynk_ncData)   = trim(lnSplit(2)) ! NAME
    dynk_cData(dynk_ncData+1) = trim(lnSplit(4)) ! FILE NAME
    dynk_ncData = dynk_ncData+1

    ! Open the file
    inquire(unit=dynk_fileUnitFUN, opened=isOpen)
    if(isOpen) then
      write(lout,"(a)") "DYNK> ERROR FUN:FILELIN Could not open file '"//trim(dynk_cData(dynk_ncData))//"'"
      iErr = .true.
      return
    end if

#ifdef BOINC
    call boincrf(dynk_cData(dynk_ncData),filename)
    open(unit=dynk_fileUnitFUN,file=filename,action="read",iostat=ioStat,status="old")
#else
    open(unit=dynk_fileUnitFUN,file=dynk_cData(dynk_ncData),action="read",iostat=ioStat,status="old")
#endif

    if(ioStat /= 0) then
      write(lout,"(a)") "DYNK> ERROR FUN:FILELIN Could not open file '"//trim(dynk_cData(dynk_ncData))//"'"
      iErr = .true.
      return
    end if

    ! Find the size of the file
    ii = 0 ! Number of data lines read
    do
      read(dynk_fileUnitFUN,"(a)",iostat=ioStat) fLine
      if(ioStat /= 0) exit

      call chr_split(fLine,lnFile,nFile,spErr)
      if(spErr) then
        write(lout,"(a)") "DYNK> ERROR FUN:FILELIN Failed to parse input line."
        iErr = .true.
        return
      end if
      if(nFile /= 2) then
        write(lout,"(a,i0)") "DYNK> ERROR FUN:FILELIN Expected two values per line, got ",nFIle
        iErr = .true.
        return
      end if

      call chr_cast(lnFile(1),x,cErr)
      call chr_cast(lnFile(2),y,cErr)

      if(ii > 0 .and. x <= x2) then ! Insane: Decreasing x
        write(lout,"(a)") "DYNK> ERROR FUN:FILELIN Reading file '"//trim(dynk_cData(dynk_ncData))//"'"
        write(lout,"(a)") "DYNK>       x values must be in increasing order"
        iErr = .true.
        return
      end if
      x2 = x
      ii = ii+1
    end do
    t = ii
    rewind(dynk_fileUnitFUN)

    call dynk_checkspace(0,2*t,0)
    ! Read the file
    ii = 0
    do
      read(dynk_fileUnitFUN,"(a)",iostat=ioStat) fLine
      if(ioStat /= 0) then ! EOF
        if(ii /= t) then
          write(lout,"(a)")       "DYNK> ERROR FUN:FILELIN Unexpected when reading file '"//trim(dynk_cData(dynk_ncData))//"'"
          write(lout,"(2(a,i0))") "DYNK>       ii = ",ii,", t = ",t
          iErr = .true.
          return
        end if
        exit
      end if

      call chr_split(fLine,lnFile,nFile,spErr)
      if(spErr) then
        write(lout,"(a)") "DYNK> ERROR FUN:FILELIN Failed to parse input line."
        iErr = .true.
        return
      end if
      if(nFile /= 2) then
        write(lout,"(a,i0)") "DYNK> ERROR FUN:FILELIN Expected two values per line, got ",nFIle
        iErr = .true.
        return
      end if

      call chr_cast(lnFile(1),x,cErr)
      call chr_cast(lnFile(2),y,cErr)

      ii = ii+1

      dynk_fData(dynk_nfData+ii)   = x
      dynk_fData(dynk_nfData+ii+t) = y
    end do

    dynk_nfData = dynk_nfData + 2*t
    dynk_funcs(dynk_nFuncs,5) = t
    close(dynk_fileUnitFUN)

  ! END CASE FILELIN

  case("PIPE")
    ! PIPE: Use a pair of UNIX FIFOs.
    ! Another program is expected to hook onto the other end of the pipe,
    ! and will recieve a message when SixTrack's dynk_computeFUN() is called.
    ! That program should then send a value back (in ASCII), which will be the new setting.

    call dynk_checkargs(nSplit,6,"FUN funname PIPE inPipeName outPipeName ID" )
    call dynk_checkspace(1,0,4)

#ifdef CR
    write(lout,"(a)") "DYNK> ERROR FUN PIPE not supported in CR version."
    iErr = .true.
    return
#endif

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_niData = dynk_niData+1
    dynk_ncData = dynk_ncData+1
    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData   ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 3             ! TYPE (PIPE)
    dynk_funcs(dynk_nFuncs,3) = dynk_niData   ! UnitNR (set below)
    dynk_funcs(dynk_nFuncs,4) = -1            ! Not used
    dynk_funcs(dynk_nFuncs,5) = -1            ! Not used

    ! Sanity checks
    if(len(lnSplit(4)) > mStrLen-1 .or. len(lnSplit(5)) > mStrLen-1 .or. len(lnSplit(6)) > mStrLen-1 ) then
      write(lout,"(a,i0)") "DYNK> ERROR FUN:PIPE got one or more strings which is too long. Max length is ",mStrLen-1
      iErr = .true.
      return
    end if

    ! Store data
    dynk_cData(dynk_ncData)   = trim(lnSplit(2)) ! NAME
    dynk_cData(dynk_ncData+1) = trim(lnSplit(4)) ! inPipe
    dynk_cData(dynk_ncData+2) = trim(lnSplit(5)) ! outPipe
    dynk_cData(dynk_ncData+3) = trim(lnSplit(6)) ! ID
    dynk_ncData = dynk_ncData+3

    call funit_requestUnit(dynk_cData(dynk_ncData+1),dynk_iData(dynk_niData))   ! fileUnit 1
    call funit_requestUnit(dynk_cData(dynk_ncData+2),dynk_iData(dynk_niData+1)) ! fileUnit 2
    dynk_niData = dynk_niData+1

    ! Look if the filenames are used in a different FUN PIPE
    t=0 ! Used to hold the index of the other pipe; t=0 if no older pipe -> open files.
    do ii=1,dynk_nFuncs-1
      if (dynk_funcs(ii,2) .eq. 3) then ! It's a PIPE
        ! Does any of the settings match?
        if(dynk_cData(dynk_funcs(ii,1)+1) == dynk_cData(dynk_ncData-2) .or. & ! InPipe filename
           dynk_cData(dynk_funcs(ii,1)+2) == dynk_cData(dynk_ncData-1)) then  ! OutPipe filename
          ! Does *all* of the settings match?
          if(dynk_cData(dynk_funcs(ii,1)+1) == dynk_cData(dynk_ncData-2) .and. & ! InPipe filename
             dynk_cData(dynk_funcs(ii,1)+2) == dynk_cData(dynk_ncData-1)) then   ! OutPipe filename
            t=ii
            write(lout,"(a)") "DYNK> FUN:PIPE '"//trim(dynk_cData(dynk_funcs(dynk_nFuncs,1)))// &
              "' using same settings as previously defined FUN '"//trim(dynk_cData(dynk_funcs(ii,1)))// &
              "' -> reusing files"
            if(dynk_cData(dynk_funcs(ii,1)+3) == dynk_cData(dynk_ncData)) then ! ID
              write(lout,"(a)") "DYNK> ERROR FUN:PIPE IDs must be different when sharing PIPEs."
              iErr = .true.
              return
            end if
            exit ! Break loop
          else ! Partial match
            write(lout,"(a)") "DYNK> ERROR FUN:PIPE Partial match of inPipe/outPipe between '"// &
              trim(dynk_cData(dynk_funcs(dynk_nFuncs,1)))//"' and '"//trim(dynk_cData(dynk_funcs(ii,1)))//"'"
              iErr = .true.
              return
          end if
        end if
      end if
    end do

    if(t == 0) then ! Must open a new set of files
      ! Open the inPipe
      write(lout,"(a)") "DYNK> Opening input pipe '"//&
        trim(dynk_cData(dynk_ncData-2))//"' for FUN '"//&
        trim(dynk_cData(dynk_ncData-3))//"', ID='"//&
        trim(dynk_cData(dynk_ncData))//"'"

      inquire(unit=dynk_iData(dynk_niData), opened=isOpen)
      if(isOpen) then
        write(lout,"(a)")"DYNK> ERROR FUN:PIPE File '"//trim(dynk_cData(dynk_ncData-2))//"' is already open"
        iErr = .true.
        return
      end if

      ! DYNK PIPE does not support the CR version, so BOINC support (call boincrf()) isn't needed
      open(unit=dynk_iData(dynk_niData),file=dynk_cData(dynk_ncData-2),action="read",iostat=ioStat,status="old")
      if(ioStat /= 0) then
        write(lout,"(a,i0)") "DYNK> ERROR FUN:PIPE Could not open file '"//trim(dynk_cData(dynk_ncData-2))//"' stat = ",ioStat
        iErr = .true.
        return
      end if

      ! Open the outPipe
      write(lout,"(a)") "DYNK> Opening output pipe '"//&
        trim(dynk_cData(dynk_ncData-1))//"' for FUN '"//&
        trim(dynk_cData(dynk_ncData-3))//"', ID='"//&
        trim(dynk_cData(dynk_ncData))//"'"

      inquire(unit=dynk_iData(dynk_niData+1), opened=isOpen)
      if(isOpen) then
        write(lout,"(a)")"DYNK> ERROR FUN:PIPE File '"//trim(dynk_cData(dynk_ncData-1))//"' is already open"
        iErr = .true.
        return
      end if

      ! DYNK PIPE does not support the CR version, so BOINC support (call boincrf()) isn't needed
      open(unit=dynk_iData(dynk_niData+1),file=dynk_cData(dynk_ncData-1),action="write",iostat=ioStat,status="old")
      if(ioStat /= 0) then
        write(lout,"(a)") "DYNK> ERROR FUN:PIPE Could not open file '"//trim(dynk_cData(dynk_ncData-1))//"' stat = ",ioStat
        iErr = .true.
        return
      end if
      write(dynk_iData(dynk_niData+1),"(a)") "DYNKPIPE !******************!" ! Once per file

    end if ! End "if (t == 0)"/must open new files

    write(dynk_iData(dynk_niData+1),"(a)") "INIT ID="//trim(dynk_cData(dynk_ncData))//" for FUN="//trim(dynk_cData(dynk_ncData-3))

  ! END CASE PIPE

  case("RANDG")
    ! RANDG: Gausian random number with mu, sigma, and optional cutoff

    call dynk_checkargs(nSplit,8,"FUN funname RANDG seed1 seed2 mu sigma cut")
    call dynk_checkspace(5,2,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_niData = dynk_niData+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData !NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 6           !TYPE (RANDG)
    dynk_funcs(dynk_nFuncs,3) = dynk_niData !seed1(initial),seed2(initial),mcut,seed1(current),seed2(current) (in dynk_iData)
    dynk_funcs(dynk_nFuncs,4) = dynk_nfData !mu, sigma (in dynk_fData)
    dynk_funcs(dynk_nFuncs,5) = -1          !ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2))               ! NAME
    call chr_cast(lnSplit(4),dynk_iData(dynk_niData),  cErr) ! seed1 (initial)
    call chr_cast(lnSplit(5),dynk_iData(dynk_niData+1),cErr) ! seed2 (initial)
    call chr_cast(lnSplit(6),dynk_fData(dynk_nfData),  cErr) ! mu
    call chr_cast(lnSplit(7),dynk_fData(dynk_nfData+1),cErr) ! sigma
    call chr_cast(lnSplit(8),dynk_iData(dynk_niData+2),cErr) ! mcut

    dynk_iData(dynk_niData+3) = dynk_iData(dynk_niData)      ! seed1 (current)
    dynk_iData(dynk_niData+4) = dynk_iData(dynk_niData+1)    ! seed2 (current)

    dynk_niData = dynk_niData+4
    dynk_nfData = dynk_nfData+1

    if(dynk_iData(dynk_funcs(dynk_nFuncs,3)+2) < 0) then
      ! mcut < 0
      write (lout,"(a)") "DYNK> ERROR FUN:RANDG mcut must be >= 0"
      iErr = .true.
      return
    end if

  ! END CASE RANDG

  case("RANDU")
    ! RANDU: Uniform random number

    call dynk_checkargs(nSplit,5,"FUN funname RANDU seed1 seed2")
    call dynk_checkspace(4,0,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_niData = dynk_niData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData !NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 7           !TYPE (RANDU)
    dynk_funcs(dynk_nFuncs,3) = dynk_niData !seed1(initial), seed2(initial), seed1(current), seed2(current)
    dynk_funcs(dynk_nFuncs,4) = -1          !ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          !ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2))               ! NAME
    call chr_cast(lnSplit(4),dynk_iData(dynk_niData),  cErr) ! seed1 (initial)
    call chr_cast(lnSplit(5),dynk_iData(dynk_niData+1),cErr) ! seed2 (initial)

    dynk_iData(dynk_niData+2) = dynk_iData(dynk_niData)      ! seed1 (current)
    dynk_iData(dynk_niData+3) = dynk_iData(dynk_niData+1)    ! seed2 (current)

    dynk_niData = dynk_niData+3

  ! END CASE RANDU

  case("RANDON")
    ! RANDON: Turn by turn ON for one turn with the probability P, else OFF

    call dynk_checkargs(nSplit,6,"FUN funname RANDON seed1 seed2 P")
    call dynk_checkspace(4,1,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_niData = dynk_niData+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData !NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 8           !TYPE (RANDON)
    dynk_funcs(dynk_nFuncs,3) = dynk_niData !seed1(initial), seed2(initial), seed1(current), seed2(current)
    dynk_funcs(dynk_nFuncs,4) = dynk_nfData !P (in dynk_fData)
    dynk_funcs(dynk_nFuncs,5) = -1          !ARG2 (unused)

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2))               ! NAME
    call chr_cast(lnSplit(4),dynk_iData(dynk_niData),  cErr) ! seed1 (initial)
    call chr_cast(lnSplit(5),dynk_iData(dynk_niData+1),cErr) ! seed2 (initial)
    call chr_cast(lnSplit(6),dynk_fData(dynk_nfData),  cErr) ! P

    dynk_iData(dynk_niData+2) = dynk_iData(dynk_niData)      ! seed1 (current)
    dynk_iData(dynk_niData+3) = dynk_iData(dynk_niData+1)    ! seed2 (current)

    dynk_niData = dynk_niData+3

  ! END CASE RANDON

  case("FIR","IIR")
    ! FIR: Finite Impulse Response filter
    ! y[n] = \sum_{i=0}^N b_i*x[n-i]
    ! where N is the order of the filter, x[] is the results from previous calls to the input
    ! function, and b_i is a set of coefficients. The coefficients are loaded from an ASCII
    ! file, formatted with three columns: the first one being the index 0...N, the second being
    ! the coefficients b_0...b_N, and the third one being the initial values of x[n]..x[n-N].
    ! When running, the values x[n]...x[n-N] are the N last results from calling baseFUN.
    ! Note that this means that at the first call, x[n-0] is pushed into x[n-1] etc., and x[n-N]
    ! is deleted; i.e. the initial x[n-N] is never used.
    !
    ! Format in dynk_fData:
    ! b_0 <- dynk_funcs(<this>,3)
    ! x[n]
    ! x_init[n] (holds the x[n]s from the input file, used to reset the FIR at the first turn)
    ! b_1
    ! x[n-1]
    ! x_init[n-1]
    ! (etc., repeat dynk_funcs(<this>,4)+1 times)
    !
    ! IIR: Infinite Impulse Response filter
    ! y[n] = \sum_{i=0}^N b_i*x[n-i] \sum_{i=1}^M a_i*y[i-n]
    ! where N=M. This is the same as FIR, except that it also uses previous values of it's own
    ! output. The input file is also identical, except adding two extra columns: One for the
    ! coefficients a_0...a_N, and one for the initial values of y[n]...y[n-N]. For both these
    ! columns, the first row (a_0 and y[n]) are ignored. For the first of these columns, the
    ! first value (a_0) is ignored and never used, while y[n-0] is pushed into y[n-1] at the
    ! first evaluation, such that the initial x[n-N] is never used (just like for x[n-N]).
    !
    ! Format in dynk_fData:
    ! b_0 <- dynk_funcs(<this>,3)
    ! x[n]
    ! x_init[n]
    ! a_0  (a_0 is never used)
    ! y[n] (zeroed for computation, used to hold previously returned value)
    ! y_init[n] (holds the y[n]s from the input file, used to reset the FIR at the first turn)
    ! b_1
    ! x[n-1]
    ! x_init[n-1]
    ! a_1
    ! y[n-1]
    ! y_init[n-1]
    ! (etc., repeat dynk_funcs(<this>,4) times)

    call dynk_checkargs(nSplit,6,"FUN funname {FIR|IIR} N filename baseFUN")

    select case(trim(lnSplit(3)))
    case("FIR")
      isFIR = .true.
    case("IIR")
      isFIR = .false.
    end select

    call chr_cast(lnSplit(4),t,cErr) ! N
    if(isFIR) then
      call dynk_checkspace(0,3*(t+1),2)
    else
      call dynk_checkspace(0,6*(t+1),2)
    end if

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers and metadata
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData   ! NAME (in dynk_cData)
    if(isFIR) then
      dynk_funcs(dynk_nFuncs,2) = 10 ! TYPE (FIR)
    else
      dynk_funcs(dynk_nFuncs,2) = 11 ! TYPE (IIR)
    end if
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData+1 ! ARG1 (start of float storage)
    dynk_funcs(dynk_nFuncs,4) = t             ! ARG2 (filter order N)
    dynk_funcs(dynk_nFuncs,5) = &             ! ARG3 (filtered function)
         dynk_findFUNindex(trim(lnSplit(6)), 1)

    dynk_cData(dynk_ncData) = trim(lnSplit(2))             ! NAME

    ! Sanity check
    if(dynk_funcs(dynk_nFuncs,5) == -1) then
      call dynk_dumpdata
      write(lout,"(a)") "DYNK> ERROR FUN:FIR/IIR Requesting function '"//trim(lnSplit(6))//"'. This FUN is unknown."
      iErr = .true.
      return
    end if
    if(len(lnSplit(5)) > mStrLen-1) then
      write(lout,"(2(a,i0))") "DYNK> ERROR FUN:FIR/IIR Got a filename name with length = ",len(lnSplit(5))," > ",mStrLen-1
      iErr = .true.
      return
    end if
    if(dynk_funcs(dynk_nFuncs,4) <= 0) then
      write(lout,"(a)") "DYNK> ERROR FUN:FIR/IIR Got N <= 0, this is not valid."
      iErr = .true.
      return
    end if

    ! More metadata
    dynk_ncData = dynk_ncData+1
    dynk_cData(dynk_ncData) = trim(lnSplit(5)) ! FILE NAME

    ! Read the file
    inquire(unit=dynk_fileUnitFUN, opened=isOpen)
    if(isOpen) then
      write(lout,"(a)") "DYNK> ERROR FUN:FIR/IIR Could not open file '"//trim(dynk_cData(dynk_ncData))//"'"
      iErr = .true.
      return
      end if
#ifdef BOINC
    call boincrf(dynk_cData(dynk_ncData),filename)
    open(unit=dynk_fileUnitFUN,file=filename,action="read",iostat=ioStat,status="old")
#else
    open(unit=dynk_fileUnitFUN,file=dynk_cData(dynk_ncData),action="read",iostat=ioStat,status="old")
#endif
    if(ioStat /= 0) then
      write(lout,"(a,i0)") "DYNK> ERROR FUN:FIR/IIR Could not open file '"//trim(dynk_cData(dynk_ncData))//"', stat = ",ioStat
      iErr = .true.
      return
    end if

    do ii=0, dynk_funcs(dynk_nFuncs,4)
      ! Reading the FIR/IIR file without CRLIBM
      read(dynk_fileUnitFUN,"(a)",iostat=ioStat) fLine
      if(ioStat /= 0) then ! EOF
        write(lout,"(a)") "DYNK> ERROR FUN:FIR/IIR Unexpected when reading file '"//trim(dynk_cData(dynk_ncData))//"'"
        iErr = .true.
        return
      end if

      call chr_split(fLine,lnFile,nFile,spErr)
      if(spErr) then
        write(lout,"(a)") "DYNK> ERROR FUN:FIR/IIR Failed to parse input line."
        iErr = .true.
        return
      end if
      if(isFIR) then
        if(nFile /= 3) then
          write(lout,"(a,i0)") "DYNK> ERROR FUN:FIR/IIR Expected 3 values per line, got ",nFIle
          iErr = .true.
          return
        end if
        call chr_cast(lnFile(1),t,cErr)
        call chr_cast(lnFile(2),x,cErr)
        call chr_cast(lnFile(3),y,cErr)
      else
        if(nFile /= 5) then
          write(lout,"(a,i0)") "DYNK> ERROR FUN:FIR/IIR Expected 5 values per line, got ",nFIle
          iErr = .true.
          return
        end if
        call chr_cast(lnFile(1),t,cErr)
        call chr_cast(lnFile(2),x,cErr)
        call chr_cast(lnFile(3),y,cErr)
        call chr_cast(lnFile(4),z,cErr)
        call chr_cast(lnFile(5),u,cErr)
      end if

      ! More sanity checks
      if(t /= ii) then
        write(lout,"(a)")       "DYNK> ERROR FUN:FIR/IIR Reading file '"//trim(dynk_cData(dynk_ncData))//"'"
        write(lout,"(2(a,i0))") "DYNK> Got line t = ",t,", expected ",ii
        iErr = .true.
        return
      end if

      ! Save data to arrays
      ! Store coefficients (x) and initial/earlier values (y) in interlaced order
      dynk_nfData = dynk_nfData+1
      dynk_fData(dynk_nfData) = x   ! b_i
      dynk_nfData = dynk_nfData+1
      dynk_fData(dynk_nfData) = y   ! x[n-1]
      dynk_nfData = dynk_nfData+1
      dynk_fData(dynk_nfData) = y   ! x_init[n-i] (Not really needed anymore, but fixing allignment is painfull)
      if(.not.isFIR) then
        dynk_nfData = dynk_nfData+1
        dynk_fData(dynk_nfData) = z ! a_i
        dynk_nfData = dynk_nfData+1
        dynk_fData(dynk_nfData) = u ! y[n-i]
        dynk_nfData = dynk_nfData+1
        dynk_fData(dynk_nfData) = u ! y_init[n-i]  (Not really needed anymore, but fixing allignment is painfull)
      end if
    end do
    close(dynk_fileUnitFUN)

  ! END CASES FIR & IIR

  case("ADD","SUB","MUL","DIV","POW") ! Operators: #20-39
    ! Two-argument operators  y = OP(f1, f2)

    call dynk_checkargs(nSplit,5,"FUN funname {ADD|SUB|MUL|DIV|POW} funname1 funname2")
    call dynk_checkspace(0,0,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData
    select case(trim(lnSplit(3)))
    case("ADD")
      dynk_funcs(dynk_nFuncs,2) = 20 ! TYPE (ADD)
    case("SUB")
      dynk_funcs(dynk_nFuncs,2) = 21 ! TYPE (SUB)
    case("MUL")
      dynk_funcs(dynk_nFuncs,2) = 22 ! TYPE (MUL)
    case("DIV")
      dynk_funcs(dynk_nFuncs,2) = 23 ! TYPE (DIV)
    case("POW")
      dynk_funcs(dynk_nFuncs,2) = 24 ! TYPE (POW)
    end select
    dynk_funcs(dynk_nFuncs,3) = dynk_findFUNindex(trim(lnSplit(4)), 1) ! Index to f1
    dynk_funcs(dynk_nFuncs,4) = dynk_findFUNindex(trim(lnSplit(5)), 1) ! Index to f2
    dynk_funcs(dynk_nFuncs,5) = -1  ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME
    ! Sanity check (string lengths are done inside dynk_findFUNindex)
    if(dynk_funcs(dynk_nFuncs,3) == -1 .or. dynk_funcs(dynk_nFuncs,4) == -1) then
      write(lout,"(a)") "DYNK> ERROR TWO ARG OPERATOR wanting functions '"//trim(lnSplit(4))//"' and '"//trim(lnSplit(5))//"', "
      write(lout,"(2(a,i0))") "DYNK>       Calculated indices: ",dynk_funcs(dynk_nFuncs,3)," and ",dynk_funcs(dynk_nFuncs,4)
      write(lout,"(a)") "DYNK>       One or both of these are unknown."
      call dynk_dumpdata
      iErr = .true.
      return
    end if

  ! END CASES ADD, SUB, MUL, DIV & POW

  case("MINUS","SQRT","SIN","COS","LOG","LOG10","EXP","ABS")
    ! One-argument operators  y = OP(f1)

    call dynk_checkargs(nSplit,4,"FUN funname {MINUS|SQRT|SIN|COS|LOG|LOG10|EXP} funname")
    call dynk_checkspace(0,0,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    select case(trim(lnSplit(3)))
    case("MINUS")
      dynk_funcs(dynk_nFuncs,2) = 30 ! TYPE (MINUS)
    case("SQRT")
      dynk_funcs(dynk_nFuncs,2) = 31 ! TYPE (SQRT)
    case("SIN")
      dynk_funcs(dynk_nFuncs,2) = 32 ! TYPE (SIN)
    case("COS")
      dynk_funcs(dynk_nFuncs,2) = 33 ! TYPE (COS)
    case("LOG")
      dynk_funcs(dynk_nFuncs,2) = 34 ! TYPE (LOG)
    case("LOG10")
      dynk_funcs(dynk_nFuncs,2) = 35 ! TYPE (LOG10)
    case("EXP")
      dynk_funcs(dynk_nFuncs,2) = 36 ! TYPE (EXP)
    case("ABS")
      dynk_funcs(dynk_nFuncs,2) = 37 ! TYPE (EXP)
    end select

    ! Index to f1
    dynk_funcs(dynk_nFuncs,3) = dynk_findFUNindex(trim(lnSplit(4)),1)
    dynk_funcs(dynk_nFuncs,5) = -1 ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME
    ! Sanity check (string lengths are done inside dynk_findFUNindex)
    if(dynk_funcs(dynk_nFuncs,3) == -1) then
      write(lout,"(a)")    "DYNK> ERROR SINGLE OPERATOR FUNC wanting function '"//trim(lnSplit(4))//"'"
      write(lout,"(a,i0)") "DYNK>       Calculated index: ",dynk_funcs(dynk_nFuncs,3)
      write(lout,"(a)")    "DYNK>       This function is unknown."
      call dynk_dumpdata
      iErr = .true.
      return
    end if

  ! END CASES MINUS, SQRT, SIN, COS, LOG, LOG10 & EXP

  ! Polynomial & Elliptical functions: # 40-59

  case("CONST")
    ! CONST: Just a constant value

    call dynk_checkargs(nSplit,4,"FUN funname CONST value")
    call dynk_checkspace(0,1,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 40          ! TYPE (CONST)
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData ! ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME
    call chr_cast(lnSplit(4),dynk_fData(dynk_nfData),cErr) ! value

  ! END CASE CONST

  case("TURN")
    ! TURN: Just the current turn number

    call dynk_checkargs(nSplit,3,"FUN funname TURN")
    call dynk_checkspace(0,0,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 41          ! TYPE (TURN)
    dynk_funcs(dynk_nFuncs,3) = -1          ! ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME

  ! END CASE TURN

  case("LIN")
    ! LIN: Linear ramp y = dy/dt*T+b

    call dynk_checkargs(nSplit,5,"FUN funname LIN dy/dt b")
    call dynk_checkspace(0,2,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 42          ! TYPE (LIN)
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData ! ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME

    call chr_cast(lnSplit(4),dynk_fData(dynk_nfData),  cErr) ! dy/dt
    call chr_cast(lnSplit(5),dynk_fData(dynk_nfData+1),cErr) ! b
    dynk_nfData = dynk_nfData + 1

  ! END CASE LIN

  case("LINSEG")
    ! LINSEG: Linear ramp between points (x1,y1) and (x2,y2)

    call dynk_checkargs(nSplit,7,"FUN funname LINSEG x1 x2 y1 y2")
    call dynk_checkspace(0,4,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 43          ! TYPE (LINSEG)
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData ! ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME

    call chr_cast(lnSplit(4),dynk_fData(dynk_nfData),  cErr) ! x1
    call chr_cast(lnSplit(5),dynk_fData(dynk_nfData+1),cErr) ! x2
    call chr_cast(lnSplit(6),dynk_fData(dynk_nfData+2),cErr) ! y1
    call chr_cast(lnSplit(7),dynk_fData(dynk_nfData+3),cErr) ! y2
    dynk_nfData = dynk_nfData + 3

    if(dynk_fData(dynk_nfData-3) == dynk_fData(dynk_nfData-2)) then
      write(lout,"(a)") "DYNK> ERROR FUN:LINSEG x1 and x2 must be different."
      iErr = .true.
      return
    end if

  ! END CASE LINSEG

  case("QUAD")
    ! QUAD: Quadratic ramp y = a*T^2 + b*T + c

    call dynk_checkargs(nSplit,6,"FUN funname QUAD a b c")
    call dynk_checkspace(0,3,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 44          ! TYPE (QUAD)
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData ! ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME

    call chr_cast(lnSplit(4),dynk_fData(dynk_nfData),  cErr) ! a
    call chr_cast(lnSplit(5),dynk_fData(dynk_nfData+1),cErr) ! b
    call chr_cast(lnSplit(6),dynk_fData(dynk_nfData+2),cErr) ! c
    dynk_nfData = dynk_nfData + 2

  ! END CASE QUAD

  case("QUADSEG")
    ! QUADSEG: Quadratic ramp y = a*T^2 + b*T + c,
    ! input as start point (x1,y1), end point (x2,y2), derivative at at x1

    call dynk_checkargs(nSplit,8,"FUN funname QUADSEG x1 x2 y1 y2 deriv")
    call dynk_checkspace(0,8,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 45          ! TYPE (QUADSEG)
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData ! ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME

    call chr_cast(lnSplit(4),x1,   cErr)
    call chr_cast(lnSplit(5),x2,   cErr)
    call chr_cast(lnSplit(6),y1,   cErr)
    call chr_cast(lnSplit(7),y2,   cErr)
    call chr_cast(lnSplit(8),deriv,cErr)

    if(x1 == x2) then
      write(lout,"(a)") "DYNK> ERROR FUN:QUADSEG x1 and x2 must be different."
      iErr = .true.
      return
    end if

    ! Compute a:
    dynk_fData(dynk_nfData)   = deriv/(x1-x2) + (y2-y1)/((x1-x2)**2)
    ! Compute b:
    dynk_fData(dynk_nfData+1) = (y2-y1)/(x2-x1) - (x1+x2)*dynk_fData(dynk_nfData)
    ! Compute c:
    dynk_fData(dynk_nfData+2) = y1 + (- x1**2 * dynk_fData(dynk_nfData) - x1 * dynk_fData(dynk_nfData+1))

    ! Store input data:
    dynk_fData(dynk_nfData+3) = x1
    dynk_fData(dynk_nfData+4) = x2
    dynk_fData(dynk_nfData+5) = y1
    dynk_fData(dynk_nfData+6) = y2
    dynk_fData(dynk_nfData+7) = deriv

    dynk_nfData = dynk_nfData + 7

  ! END CASE QUADSEG

  case("SINF","COSF","COSF_RIPP")
    ! Trancedental functions: #60-79
    ! SINF     : Sin functions y = A*sin(omega*T+phi)
    ! COSF     : Cos functions y = A*cos(omega*T+phi)
    ! COSF_RIPP: Cos functions y = A*cos(2*pi*(T-1)/period+phi)

    call dynk_checkargs(nSplit,6,"FUN funname {SINF|COSF|COSF_RIPP} amplitude {omega|period} phase")
    call dynk_checkspace(0,3,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData
    select case(trim(lnSplit(3)))
    case("SINF")
      dynk_funcs(dynk_nFuncs,2) = 60 ! TYPE (SINF)
    case("COSF")
      dynk_funcs(dynk_nFuncs,2) = 61 ! TYPE (COSF)
    case ("COSF_RIPP")
      dynk_funcs(dynk_nFuncs,2) = 62 ! TYPE (COSF_RIPP)
    end select
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData ! ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME

    call chr_cast(lnSplit(4),dynk_fData(dynk_nfData),  cErr) ! A
    call chr_cast(lnSplit(5),dynk_fData(dynk_nfData+1),cErr) ! omega
    call chr_cast(lnSplit(6),dynk_fData(dynk_nfData+2),cErr) ! phi
    dynk_nfData = dynk_nfData + 2

  ! END CASE SINF, COSF & COSF_RIPP

  case("PELP")
    ! PELP: Parabolic/exponential/linear/parabolic
    ! From "Field Computation for Accelerator Magnets:
    ! Analytical and Numerical Methods for Electromagnetic Design and Optimization"
    ! By Dr.-Ing. Stephan Russenschuck
    ! Appendix C: "Ramping the LHC Dipoles"

    call dynk_checkargs(nSplit,10,"FUN funname PELP tinj Iinj Inom A D R te")
    call dynk_checkspace(0,13,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_nfData = dynk_nfData+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 80          ! TYPE (PELP)
    dynk_funcs(dynk_nFuncs,3) = dynk_nfData ! ARG1
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME

    ! Read and calculate parameters
    call chr_cast(lnSplit(4), tinj,cErr)
    call chr_cast(lnSplit(5), Iinj,cErr)
    call chr_cast(lnSplit(6), Inom,cErr)
    call chr_cast(lnSplit(7), A,   cErr)
    call chr_cast(lnSplit(8), D,   cErr)
    call chr_cast(lnSplit(9), R,   cErr)
    call chr_cast(lnSplit(10),te,  cErr)

    derivI_te = A*(te-tinj)                 ! nostore
    I_te      = (A/2.0)*(te-tinj)**2 + Iinj ! nostore
    bexp      = derivI_te/I_te
    aexp      = exp_mb(-bexp*te)*I_te
    t1        = log_mb(R/(aexp*bexp))/bexp
    I1        = aexp*exp_mb(bexp*t1)
    td        = (Inom-I1)/R + (t1 - R/(2*D))
    tnom      = td + R/D

    if(dynk_debug) then
      write(lout,"(a,f14.6)") "DYNK> DEBUG FUN:PELP"
      write(lout,"(a,f14.6)") "DYNK> DEBUG tinj = ", tinj
      write(lout,"(a,f14.6)") "DYNK> DEBUG Iinj = ", Iinj
      write(lout,"(a,f14.6)") "DYNK> DEBUG Inom = ", Inom
      write(lout,"(a,f14.6)") "DYNK> DEBUG A    = ", A
      write(lout,"(a,f14.6)") "DYNK> DEBUG D    = ", D
      write(lout,"(a,f14.6)") "DYNK> DEBUG R    = ", R
      write(lout,"(a,f14.6)") "DYNK> DEBUG te   = ", te
      write(lout,"(a,f14.6)") "DYNK> DEBUG "
      write(lout,"(a,f14.6)") "DYNK> DEBUG derivI_te = ", derivI_te
      write(lout,"(a,f14.6)") "DYNK> DEBUG I_te      = ", I_te
      write(lout,"(a,f14.6)") "DYNK> DEBUG bexp      = ", bexp
      write(lout,"(a,f14.6)") "DYNK> DEBUG aexp      = ", aexp
      write(lout,"(a,f14.6)") "DYNK> DEBUG t1        = ", t1
      write(lout,"(a,f14.6)") "DYNK> DEBUG I1        = ", I1
      write(lout,"(a,f14.6)") "DYNK> DEBUG td        = ", td
      write(lout,"(a,f14.6)") "DYNK> DEBUG tnom      = ", tnom
    end if

    if(.not. (tinj < te .and. te < t1 .and. t1 < td .and. td < tnom)) then
      write(lout,"(a)") "DYNK> ERROR FUN:PELP Order of times not correct."
      iErr = .true.
      return
    end if

    ! Store: Times
    dynk_fData(dynk_nfData)    = tinj
    dynk_fData(dynk_nfData+ 1) = te
    dynk_fData(dynk_nfData+ 2) = t1
    dynk_fData(dynk_nfData+ 3) = td
    dynk_fData(dynk_nfData+ 4) = tnom
    ! Store: Parameters / section1 (parabola)
    dynk_fData(dynk_nfData+ 5) = Iinj
    dynk_fData(dynk_nfData+ 6) = A
    ! Store: Parameters / section2 (exponential)
    dynk_fData(dynk_nfData+ 7) = aexp
    dynk_fData(dynk_nfData+ 8) = bexp
    ! Store: Parameters / section3 (linear)
    dynk_fData(dynk_nfData+ 9) = I1
    dynk_fData(dynk_nfData+10) = R
    ! Store: Parameters / section4 (parabola)
    dynk_fData(dynk_nfData+11) = D
    dynk_fData(dynk_nfData+12) = Inom

    dynk_nfData = dynk_nfData + 12

  ! END CASE PELP

  case("ONOFF")
    ! ONOFF: On for p1 turns, then off for the rest of the period p2

    call dynk_checkargs(nSplit,5,"FUN funname ONOFF p1 p2")
    call dynk_checkspace(0,0,1)

    ! Set pointers to start of funs data blocks
    dynk_nFuncs = dynk_nFuncs+1
    dynk_ncData = dynk_ncData+1

    ! Store pointers
    dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
    dynk_funcs(dynk_nFuncs,2) = 81          ! TYPE (ONOFF)
    dynk_funcs(dynk_nFuncs,3) = -1          ! ARG1 (p1)
    dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2 (p2)
    dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3 (unused)

    ! Store data
    dynk_cData(dynk_ncData) = trim(lnSplit(2)) ! NAME

    call chr_cast(lnSplit(4),dynk_funcs(dynk_nFuncs,3),cErr) ! p1
    call chr_cast(lnSplit(5),dynk_funcs(dynk_nFuncs,4),cErr) ! p2

    ! Check for bad input
    if(dynk_funcs(dynk_nFuncs,3) < 0 .or. dynk_funcs(dynk_nFuncs,4) <= 1 .or. &
       dynk_funcs(dynk_nFuncs,4) < dynk_funcs(dynk_nFuncs,3)) then
      write(lout,"(a)") "DYNK> ERROR FUN:ONOFF Expected p1 >= 0, p2 > 1, p1 <= p2"
      iErr = .true.
      return
    end if

  ! END CASE ONOFF

  case default
    ! UNKNOWN function
    write(lout,"(a)") "DYNK> ERROR Unknown function in FUN"
    iErr = .true.
    return
  end select

end subroutine dynk_parseFUN

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-07-19
!  Check the arguments
! ================================================================================================ !
subroutine dynk_checkargs(nActual, nExpected, correctSyntax)

  use crcoall
  implicit none

  integer,          intent(in) :: nActual
  integer,          intent(in) :: nExpected
  character(len=*), intent(in) :: correctSyntax

  if(nActual /= nExpected) then
    write(lout,"(2(a,i0))") "DYNK> ERROR Function expected ",nExpected," arguments, got ",nActual
    write(lout,"(a)")       "CYNK>       Correct Syntax: "//correctSyntax
    call prror(-1)
  end if

end subroutine dynk_checkargs

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-28
!  Checks if the expr_dynk arrays are large enough, and allocate more memory if they're not.
!  - If a small increase is requested, just allocate a new chunk of 500 numbers or 100 strings
!  - If a large number of new elements is requested, allocate this in one chunk
! ================================================================================================ !
subroutine dynk_checkspace(iReq,fReq,cReq)

  use crcoall
  use numerical_constants

  implicit none

  integer, intent(in) :: iReq
  integer, intent(in) :: fReq
  integer, intent(in) :: cReq

  integer iNeeded, fNeeded, cNeeded

  iNeeded = dynk_niData + iReq
  fNeeded = dynk_nfData + fReq
  cNeeded = dynk_ncData + cReq

  if(iNeeded > dynk_maxiData) then
    if(iReq < 500) then
      dynk_maxiData = dynk_maxiData + 500
    else
      dynk_maxiData = iNeeded
    end if
    call alloc(dynk_iData,dynk_maxiData,0,"dynk_iData")
  end if
  if(fNeeded > dynk_maxfData) then
    if(fReq < 500) then
      dynk_maxfData = dynk_maxfData + 500
    else
      dynk_maxfData = fNeeded
    end if
    call alloc(dynk_fData,dynk_maxfData,zero,"dynk_fData")
  end if
  if(cNeeded > dynk_maxcData) then
    if(cReq < 200) then
      dynk_maxcData = dynk_maxcData + 200
    else
      dynk_maxcData = cNeeded
    end if
    call alloc(dynk_cData,mStrLen,dynk_maxcData,str_dSpace,"dynk_cData")
  end if

end subroutine dynk_checkspace

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 15-10-2014
!  - Parse SET lines in the fort.3 input file.
! ================================================================================================ !
subroutine dynk_parseSET(inLine, iErr)

  use crcoall
  use mod_alloc

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable :: lnSplit(:)
  integer nSplit, ii
  logical spErr, cErr

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lout,"(a)") "DYNK> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(dynk_nSets+1 > dynk_maxSets) then
    dynk_maxSets = dynk_maxSets + 10
    call alloc(dynk_cSets,       mStrLen,dynk_maxSets,2, str_dSpace,"dynk_cSets")
    call alloc(dynk_cSets_unique,mStrLen,dynk_maxSets,2, str_dSpace,"dynk_cSets_unique")
#ifdef CR
    call alloc(dynk_fSets_cr,            dynk_maxSets,   zero,      "dynk_fSets_cr")
#endif
    call alloc(dynk_sets,                dynk_maxSets,4, 0,         "dynk_sets")
  end if

  if(nSplit /= 7) then
    write(lout,"(a,i0)") "DYNK> ERROR Expected 7 fields on line while parsing SET, got ",nSplit
    iErr = .true.
    return
  end if

  dynk_nSets = dynk_nSets + 1

  ! function_name -> function index
  dynk_sets(dynk_nSets,1) = dynk_findFUNindex(lnSplit(4),1)
  call chr_cast(lnSplit(5),dynk_sets(dynk_nSets,2),cErr) ! startTurn
  call chr_cast(lnSplit(6),dynk_sets(dynk_nSets,3),cErr) ! endTurn
  call chr_cast(lnSplit(7),dynk_sets(dynk_nSets,4),cErr) ! turnShift

  ! Sanity check on string lengths
  if(len_trim(lnSplit(2)) > mNameLen) then
    write(lout,"(a,i0)") "DYNK> ERROR SET got an element name with length ",len_trim(lnSplit(2)),&
      ", but max is ",mNameLen
    iErr = .true.
    return
  end if
  if(len_trim(lnSplit(3)) > mStrLen-1) then
    write(lout,"(a,i0)") "DYNK> ERROR The attribute name '"//trim(lnSplit(3))//"' is too long. Max length is ",mStrLen-1
    iErr = .true.
    return
  end if

  ! OK -- save them!
  dynk_cSets(dynk_nSets,1) = trim(lnSplit(2)) ! element_name
  dynk_cSets(dynk_nSets,2) = trim(lnSplit(3)) ! attribute_name

  ! Sanity check
  if(dynk_sets(dynk_nSets,1) == -1) then
    write(lout,"(a)")    "DYNK> ERROR SET wanting function '"//lnSplit(4)//"'"
    write(lout,"(a,i0)") "DYNK>       Calculated index ",dynk_sets(dynk_nSets,1)
    write(lout,"(a)")    "DYNK>       This function is not known"
    iErr = .true.
    return
  end if

  if(dynk_sets(dynk_nSets,3) /= -1 .and. dynk_sets(dynk_nSets,2) > dynk_sets(dynk_nSets,3)) then
    write(lout,"(a)")    "DYNK> ERROR SET got first turn number > last turn number."
    write(lout,"(a,i0)") "DYNK>       first = ",dynk_sets(dynk_nSets,2)
    write(lout,"(a,i0)") "DYNK>       last  = ",dynk_sets(dynk_nSets,3)
    write(lout,"(a,i0)") "DYNK>       SET #   ",dynk_nSets
    iErr = .true.
    return
  end if

  if(dynk_sets(dynk_nSets,2) <= 0 .or. dynk_sets(dynk_nSets,3) < -1 .or. dynk_sets(dynk_nSets,3) == 0) then
    write(lout,"(a)") "DYNK> ERROR SET got turn number <= 0 (not last = -1 meaning infinity)"
    write(lout,"(a,i0)") "DYNK>       first = ",dynk_sets(dynk_nSets,2)
    write(lout,"(a,i0)") "DYNK>       last  = ",dynk_sets(dynk_nSets,3)
    write(lout,"(a,i0)") "DYNK>       SET #   ",dynk_nSets
    iErr = .true.
    return
  end if

end subroutine dynk_parseSET

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-07-19
!  - Find and return the index in the ifuncs array to the function with name funName,
!    which should be zero-padded.
!  - Return -1 if nothing was found.
!
!  Note: It is expected that the length of funName_input is equal or less than mStrLen,
!        and if it equal, that it is a zero-terminated string.
! ================================================================================================ !
integer function dynk_findFUNindex(funName, startFrom)

  use crcoall
  implicit none

  character(len=*), intent(in) :: funName
  integer,          intent(in) :: startFrom

  integer i

  dynk_findFUNindex = -1

  do i=startFrom,dynk_nFuncs
    if(dynk_cData(dynk_funcs(i,1)) == funName) then
      dynk_findFUNindex = i
      return
    end if
  end do

end function dynk_findFUNindex

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-11
!  - Find and return the index in the sets array to the set which matches element_name and att_name,
!    which should be zero-padded.
!  - Return -1 if nothing was found.
!
!  Note: It is expected that the length of element_name and att_name is exactly mStrLen .
! ================================================================================================ !
integer function dynk_findSETindex(elementName, attName, startFrom)

  implicit none

  character(mStrLen), intent(in) :: elementName
  character(mStrLen), intent(in) :: attName
  integer,            intent(in) :: startFrom

  integer ii

  dynk_findSETindex = -1

  do ii=startfrom, dynk_nSets
    if(dynk_cSets(ii,1) == elementName .and. dynk_cSets(ii,2) == attName) then
      dynk_findSETindex = ii
      exit
    end if
  end do

end function dynk_findSETindex

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-11
!  - Check that DYNK block input in fort.3 was sane
! ================================================================================================ !
subroutine dynk_inputSanityCheck

  use crcoall
  implicit none

  integer ii, jj
  integer biggestTurn ! Used as a replacement for ending turn -1 (infinity)
  logical sane

  sane = .true.

  ! Check that there are no doubly-defined function names
  do ii=1, dynk_nFuncs-1
    jj = dynk_findFUNindex(dynk_cData(dynk_funcs(ii,1)),ii+1)
    if(jj /= -1) then
      sane = .false.
      write(lout,"(2(a,i0))") "DYNK> Insane function ",ii," has the same name as ",jj
    end if
  end do

  ! Check that no SETS work on same elem/att at same time
  biggestTurn = 1
  do ii=1, dynk_nSets
    if(dynk_sets(ii,3) > biggestTurn) then
      biggestTurn = dynk_sets(ii,3)
    end if
  end do
  biggestTurn = biggestTurn+1 ! Make sure it is unique
  if(biggestTurn <= 0) then
    ! In case of integer overflow
    write(lout,"(a)") "DYNK> ERROR Integer overflow in sanity check"
    call prror(-1)
  end if

  ! Do the search!
  do ii=1, dynk_nSets-1
    if(dynk_sets(ii,3) == -1) dynk_sets(ii,3) = biggestTurn

    jj = ii
    do
      ! Only check SETs affecting the same elem/att
      jj = dynk_findSETindex(dynk_cSets(ii,1),dynk_cSets(ii,2),jj+1)

      if(jj == -1) exit ! next outer loop
      if(dynk_sets(jj,3) == -1) dynk_sets(jj,3) = biggestTurn
      if(dynk_sets(jj,2) <= dynk_sets(ii,2) .and. dynk_sets(jj,3) >= dynk_sets(ii,2)) then
        sane = .false.
        write(lout,"(a,i4,a,i8,a,i4,a,i8,a,i4,a,i8,a,i4)") "DYNK> Insane: Lower edge of SET #", jj, &
          " =", dynk_sets(jj,2)," <= lower edge of SET #",ii, &
          " =", dynk_sets(ii,2),"; and also higer edge of SET #",jj, &
          " =", dynk_sets(jj,3)," >= lower edge of SET #", ii
      else if(dynk_sets(jj,3) >= dynk_sets(ii,3) .and. dynk_sets(jj,2) <= dynk_sets(ii,3)) then
        sane = .false.
        write(lout,"(a,i4,a,i8,a,i4,a,i8,a,i4,a,i8,a,i4)") "DYNK> Insane: Upper edge of SET #", jj, &
          " =", dynk_sets(jj,3)," >= upper edge of SET #",ii, &
          " =", dynk_sets(ii,3),"; and also lower edge of SET #",jj, &
          " =", dynk_sets(jj,2)," <= upper edge of SET #", ii
      else if(dynk_sets(jj,2) >= dynk_sets(ii,2) .and. dynk_sets(jj,3) <= dynk_sets(ii,3)) then
        ! (other way round gets caugth by the first "if")
        sane = .false.
        write(lout,"(a,i4,a,i8,a,i8,a,a,i4,a,i8,a,i8,a)") "DYNK> Insane: SET #", jj, &
          " = (", dynk_sets(jj,2),", ", dynk_sets(jj,3), ")", &
          " is inside SET #", ii, " = (",  &
          dynk_sets(ii,2),", ", dynk_sets(ii,3), ")"
      end if
      if(dynk_sets(jj,3) == biggestTurn) dynk_sets(jj,3) = -1
    end do

    if (dynk_sets(ii,3) == biggestTurn) dynk_sets(ii,3) = -1
  end do

  if(.not. sane) then
    write(lout,"(a)") "DYNK> ERROR Input was insane"
    ! call dynk_dumpdata
    call prror(-1)
  else if(sane .and. dynk_debug) then
    write(lout,"(a)") "DYNK> DYNK input was sane"
  end if

end subroutine dynk_inputSanityCheck

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-11
!  - Dump arrays with DYNK FUN and SET data to the std. output for debugging
! ================================================================================================ !
subroutine dynk_dumpdata

  use crcoall
  implicit none

  integer ii

  write(lout,"(a)")      "DYNK> DEBUG OPTIONS:"
  write(lout,"(a,l1)")   "DYNK> DEBUG  * dynk_enabled    = ",dynk_enabled
  write(lout,"(a,l1)")   "DYNK> DEBUG  * dynk_debug      = ",dynk_debug
  write(lout,"(a,l1)")   "DYNK> DEBUG  * dynk_noDynkSets = ",dynk_noDynkSets

  write(lout,"(a)")      "DYNK> DEBUG FUN:"
  write(lout,"(a,i0,a)") "DYNK> DEBUG ifuncs(",dynk_nFuncs,"):"
  do ii=1,dynk_nFuncs
    write(lout,"(a,i5,a,6(1x,i5))") "DYNK> DEBUG  * ",ii,": ",dynk_funcs(ii,:)
  end do
  write(lout,"(a,i0,a)") "DYNK> DEBUG dynk_iData(",dynk_niData,"):"
  do ii=1,dynk_niData
    write(lout,"(a,i5,a,i0)") "DYNK> DEBUG  * ",ii,": ",dynk_iData(ii)
  end do
  write(lout,"(a,i0,a)") "DYNK> DEBUG dynk_fData(",dynk_nfData,"):"
  do ii=1,dynk_nfData
    write(lout,"(a,i5,a,e16.9)") "DYNK> DEBUG  * ",ii,": ",dynk_fData(ii)
  end do
  write(lout,"(a,i0,a)") "DYNK> DEBUG dynk_cData(",dynk_ncData,"):"
  do ii=1,dynk_ncData
    write(lout,"(a,i5,a)") "DYNK> DEBUG  * ",ii,": '"//trim(dynk_cData(ii))//"'"
  end do

  write(lout,"(a)")         "DYNK> DEBUG SET:"
  write(lout,"(a,3(i0,a))") "DYNK> DEBUG sets(",dynk_nSets,",:) csets(",dynk_nSets,",1) csets(",dynk_nSets,",2):"
  do ii=1,dynk_nSets
    write(lout,"(a,i5,a,4(1x,i5),a)") "DYNK> DEBUG  * ",ii,":",dynk_sets(ii,:), &
      " '"//trim(dynk_cSets(ii,1))//"' '"//trim(dynk_cSets(ii,2))//"'"
  end do

end subroutine dynk_dumpdata

! ================================================================================================ !
!  K. Sjobak, BE-ABP-HSS
!  Last modified: 2018-08-12
!  - Save original values for GET functions and sanity check that elements/attributes for SET
!    actually exist.
! ================================================================================================ !
subroutine dynk_pretrack

  use crcoall
  use mod_common
  use mod_commond

  implicit none

  ! Temp variables
  integer ii,jj
  character(mStrLen) element_name_s, att_name_s
  logical found, badelem
  if(dynk_debug) then
    write(lout,"(a)") "DYNK> DEBUG In pretrack"
  end if

  ! Find which elem/attr combos are affected by SET
  dynk_nSets_unique = 0 !Assuming this is only run once
  do ii=1,dynk_nSets
    if (dynk_findSETindex(dynk_cSets(ii,1),dynk_cSets(ii,2), ii+1 ) .eq. -1) then
      ! Last SET which has this attribute, store it
      dynk_nSets_unique = dynk_nSets_unique+1

      dynk_cSets_unique(dynk_nSets_unique,1) = dynk_cSets(ii,1)
      dynk_cSets_unique(dynk_nSets_unique,2) = dynk_cSets(ii,2)

      ! Sanity check: Does the element actually exist?
      element_name_s = trim(dynk_cSets_unique(dynk_nSets_unique,1))
      att_name_s     = trim(dynk_cSets_unique(dynk_nSets_unique,2))
      found          = .false.

      ! Special case: the element name GLOBAL-VARS (not a real element)
      ! can be used to redefine a global variable by some function.
      if(element_name_s == "GLOBAL-VARS") then
        found   = .true.
        badelem = .false.

        if(att_name_s == "E0") then
          if(idp == 0 .or. ition == 0) then ! 4d tracking..
            write(lout,"(a)") "DYNK> ERROR Attribute '"//att_name_s//"' is not valid for 'GLOBAL-VARS' when doing 4d tracking"
            call prror(-1)
          end if
        else
          badelem=.true.
        end if

        if(badelem) then
          write(lout,"(a)") "DYNK> ERROR Attribute '"//att_name_s//"' is not valid for 'GLOBAL-VARS'"
          call prror(-1)
        end if
      end if

      do jj=1,il
        if(bez(jj) == element_name_s) then
          found = .true.

          ! Check that the element type and attribute is supported
          ! Check that the element can be used now
          badelem = .false.
          if(abs(kz(jj)) >= 1 .and. abs(kz(jj)) <= 10) then !thin kicks
            if(att_name_s /= "average_ms") then
              badelem = .true.
            end if
          else if(abs(kz(jj)) == 12) then ! cavity
            if(.not.(att_name_s == "voltage"  .or. att_name_s == "harmonic" .or. att_name_s == "lag_angle")) then
              badelem = .true.
            end if
            if(kp(jj) /= 6) then
              write(lout,"(a)") "DYNK> ERROR Want to modify DISABLED RF cavity named '"//element_name_s//"'"
              write(lout,"(a)") "DYNK>       Please make sure that the voltage and harmonic number in the "//&
                "SINGLE ELEMENTS block is not 0!"
              call prror(-1)
            end if
            if(nvar == 5) then
              write(lout,"(a)") "DYNK> ERROR Want to modify RF cavity named '"//element_name_s//"', but nvars=5 (from DIFF block)."
            end if
          else if(abs(kz(jj)) == 23 .or. abs(kz(jj)) == 26 .or. abs(kz(jj)) == 27 .or. abs(kz(jj)) == 28) then
            if(.not.(att_name_s == "voltage" .or. att_name_s == "frequency" .or. att_name_s == "phase")) then
              badelem = .true.
            end if
          end if

          ! Special case:
          ! Should the error only occur if we actually have a GLOBAL-VARS element?
          if(bez(jj) == "GLOBAL-VARS") then
            write(lout,"(a)") "DYNK> ERROR Element found 'GLOBAL-VARS' is not a valid element name, it is reserved"
            call prror(-1)
          end if

          if(badelem) then
            write(lout,"(a,i0)") "DYNK> ERROR Attribute '"//att_name_s//"' is not valid for element '"//element_name_s//"'"//&
              " which is of type ",kz(jj)
            call prror(-1)
          end if
        end if
      end do

      if(.not. found) then
        write(lout,"(a)") "DYNK> ERROR Element '",element_name_s,"' was not found"
        call prror(-1)
      end if
    end if
  end do

  ! Save original values for GET functions
  do ii=1,dynk_nFuncs
    if(dynk_funcs(ii,2) == 0) then ! GET
      dynk_fData(dynk_funcs(ii,3)) = dynk_getvalue(dynk_cData(dynk_funcs(ii,1)+1),dynk_cData(dynk_funcs(ii,1)+2))
    end if
  end do

  if(dynk_debug) call dynk_dumpdata

end subroutine dynk_pretrack

! ================================================================================================ !
!  A.Mereghetti, for the FLUKA Team
!  K.Sjobak, A. Santamaria, V.K Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-12
!  - Actually apply dynamic kicks
!  - Always in main code
!
!  For each element (group) flagged with SET(R), compute the new value using dynk_computeFUN() at
!  the given (shifted) turn number using the specified FUN function. The values are stored in the
!  element using dynk_setvalue().
!
!  Also writes the file "dynksets.dat", only on the first turn.
! ================================================================================================ !
subroutine dynk_apply(turn)

  use crcoall
  use mod_common
  use mod_commont
  use mod_commonmn
  use string_tools

  implicit none

#ifdef BOINC
  character(len=256) filename
#endif

  ! interface variables
  integer, intent(in) :: turn

  ! temporary variables
  integer ii, jj, shiftedTurn
  logical isOpen
  real(kind=fPrec) getvaldata, newValue

  ! Used for output dynksets.dat
  character(len=mStrLen) whichFUN(dynk_nSets_unique) ! Which function was used to set a given elem/attr?
  integer whichSET(dynk_nSets_unique)                ! Which SET was used for a given elem/attr?

  ! Temp variable for padding the strings for output to dynksets.dat
  character(20) outstring_tmp1,outstring_tmp2,outstring_tmp3

  if(dynk_debug) then
    write (lout,"(a,i0)") "DYNK> DEBUG In apply at turn ",turn
  end if

  ! Initialize variables (every call)
  do jj=1, dynk_nSets_unique
    whichSET(jj) = -1
    whichFUN(jj) = " "
  end do

  ! First-turn initialization, including some parts which are specific for collimat.
  if(turn == 1) then
    ! Open dynksets.dat
#ifdef CR
    ! Could have loaded a CR just before tracking starts;
    ! In this case, the dynksets is already open and positioned,
    ! so don't try to open the file again.
    if(dynk_filePos  == -1) then
#endif
      inquire(unit=dynk_fileUnit, opened=isOpen)
      if(isOpen) then
        write(lout,"(a)") "DYNK> ERROR Could not open file 'dynksets.dat'"
        call prror(-1)
      end if
#ifdef BOINC
      call boincrf("dynksets.dat",filename)
      open(unit=dynk_fileUnit,file=filename,status="replace",action="write")
#else
      open(unit=dynk_fileUnit,file="dynksets.dat",status="replace",action="write")
#endif

      if(dynk_noDynkSets) then
        write(dynk_fileUnit,"(a)") "### DYNK file output was disabled with flag NOFILE in fort.3 ###"
      else
        write(dynk_fileUnit,"(a1,1x,a10,2(1x,a20),1x,a4,1x,a20,a16)") "#",&
          "turn", chr_rPad("element",20),chr_rPad("attribute",20),"idx",chr_rPad("funname",20),"value"
      end if
#ifdef CR
      ! Note: To be able to reposition, each line should be shorter than 255 chars
      dynk_filePos = 1

      ! Flush the unit
      endfile(dynk_fileUnit,iostat=ierro)
      backspace(dynk_fileUnit,iostat=ierro)
    end if ! END if(dynk_filePos == -1)
#endif
  end if ! END "if (turn == 1) then"

  ! Apply the sets
  do ii=1,dynk_nSets
    ! Sanity check already confirms that only a single SET
    ! is active on a given element:attribute on a given turn.

    ! Active in this turn?
    if(turn >= dynk_sets(ii,2) .and. (turn <= dynk_sets(ii,3) .or. dynk_sets(ii,3) == -1)) then

      ! Shifting
      shiftedTurn = turn + dynk_sets(ii,4)

      ! Set the value
      newValue = dynk_computeFUN(dynk_sets(ii,1),shiftedTurn)
      if(dynk_debug) then
        write(lout, "(a,i5,a,i8,a,e16.9)") "DYNK> DEBUG Applying set #", ii, " on '"//trim(dynk_cSets(ii,1))// &
          "':'"//trim(dynk_cSets(ii,2))//"', shiftedTurn = ",shiftedTurn,", value = ",newValue
      end if
      call dynk_setvalue(dynk_cSets(ii,1),dynk_cSets(ii,2),newValue)

      if(dynk_debug) then
        getvaldata = dynk_getvalue(dynk_cSets(ii,1),dynk_cSets(ii,2))
        write(lout, "(a,e16.9)") "DYNK> DEBUG Read back value = ", getvaldata
        if(getvaldata /= newValue) then
          write(lout,"(a)") "DYNK> DEBUG WARNING Read back value differs from set!"
        end if
      end if

      ! For the output file: Which function was used?
      do jj=1, dynk_nSets_unique
        if(dynk_cSets(ii,1) == dynk_cSets_unique(jj,1) .and. dynk_cSets(ii,2) == dynk_cSets_unique(jj,2)) then
          whichSET(jj)=ii
          whichFUN(jj)=dynk_cData(dynk_funcs(dynk_sets(ii,1),1))
        end if
      end do
    end if
  end do

  ! Write output file
  if(.not.dynk_noDynkSets) then
    do jj=1,dynk_nSets_unique
      getvaldata = dynk_getvalue(dynk_cSets_unique(jj,1),dynk_cSets_unique(jj,2))

      if(whichSET(jj) == -1) then
        whichFUN(jj) = "N/A"
      end if

      write(dynk_fileUnit,"(i12,2(1x,a20),1x,i4,1x,a20,e16.9)") turn, &
        chr_rPad(dynk_cSets_unique(jj,1),20),chr_rPad(dynk_cSets_unique(jj,2),20),&
        whichSET(jj),chr_rPad(whichFUN(jj),20),getvaldata
    end do

#ifdef CR
    ! Note: To be able to reposition, each line should be shorter than 255 chars
    dynk_filePos = dynk_filePos+dynk_nSets_unique
#endif
    ! Flush the unit
    endfile(dynk_fileUnit,iostat=ierro)
    backspace(dynk_fileUnit,iostat=ierro)
  end if

end subroutine dynk_apply

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-12
!  - Compute the value of a given DYNK function (funNum) for the given turn
! ================================================================================================ !
recursive real(kind=fPrec) function dynk_computeFUN(funNum, turn) result(retval)

  use crcoall
  use mod_common
  use mod_ranecu
  use numerical_constants, only : pi
  use utils

  implicit none

  integer, intent(in) :: funNum, turn

  ! Temporaries for FILELIN
  integer filelin_start, filelin_xypoints

  character(len=:), allocatable :: lnFile(:)
  character(len=mInputLn)       :: fLine
  integer nFile, ioStat
  logical spErr, cErr

  ! Temporaries for random generator functions
  integer tmpseed1, tmpseed2
  real(kind=fPrec) ranecu_rvec(1)

  ! General temporaries
  integer foff  ! Base offset into fexpr array
  integer ii,jj ! Loop variable

  if(funNum < 1 .or. funNum > dynk_nFuncs) then
    write(lout,"(a,i0)") "DYNK> ERROR computeFUN: funNum = ",funNum
    write(lout,"(a,i0)") "DYNK>       Invalid funNum, dynk_nFuncs = ", dynk_nFuncs
    if(dynk_debug) call dynk_dumpdata
    call prror(-1)
  end if

  select case(dynk_funcs(funNum,2)) ! WHICH FUNCTION TYPE?

  case(0) ! GET
    retval = dynk_fData(dynk_funcs(funNum,3))

  case(1) ! FILE
    if(turn > dynk_funcs(funNum,5)) then
      write(lout,"(2(a,i0))")"DYNK> ERROR computeFUN:FILE funNum = ",funNum,"turn = ",turn
      write(lout,"(a,i0)")   "DYNK>       Turn > length of file = ",dynk_funcs(funNum,5)
      if(dynk_debug) call dynk_dumpdata
      call prror(-1)
    else if(turn < 1) then
      write(lout,"(2(a,i0))")"DYNK> ERROR computeFUN:FILE funNum = ",funNum,"turn = ",turn
      write(lout,"(a)")      "DYNK>       Turn < 1, check your turn-shift!"
      if(dynk_debug) call dynk_dumpdata
      call prror(-1)
    end if
    retval = dynk_fData(dynk_funcs(funNum,4)+turn-1)

  case(2) ! FILELIN
    filelin_start    = dynk_funcs(funNum,4)
    filelin_xypoints = dynk_funcs(funNum,5)
    ! Pass the correct array views/sections to lininterp
    retval = lininterp( real(turn,fPrec), &
             dynk_fData(filelin_start:filelin_start+filelin_xypoints-1), &
             dynk_fData(filelin_start +  filelin_xypoints: &
                        filelin_start +2*filelin_xypoints-1), &
                        filelin_xypoints )

  case(3) ! PIPE
    write(dynk_iData(dynk_funcs(funNum,3)+1),"(a,i7)") &
      "GET ID="//trim(dynk_cData(dynk_funcs(funNum,1)+3))//" TURN=",turn

    read(dynk_iData(dynk_funcs(funNum,3)),"(a)",iostat=ioStat) fLine
    if(ioStat /= 0) then
      write(lout,"(a)") "DYNK> ERROR computeFUN:PIPE Failed to open file."
      call prror(-1)
    end if

    call chr_split(fLine,lnFile,nFile,spErr)
    if(spErr) then
      write(lout,"(a)") "DYNK> ERROR computeFUN:PIPE Failed to parse input line."
      call prror(-1)
    end if
    if(nFile /= 1) then
      write(lout,"(a,i0)") "DYNK> ERROR computeFUN:PIPE Expected one values per line, got ",nFIle
      call prror(-1)
    end if

    call chr_cast(lnFile(1),retval,cErr)

  case(6) ! RANDG
    ! Save old seeds and load our current seeds
    call recuut(tmpseed1,tmpseed2)
    call recuin(dynk_iData(dynk_funcs(funNum,3)+3),dynk_iData(dynk_funcs(funNum,3)+4))
    ! Run generator for 1 value with current mcut
    call ranecu(ranecu_rvec,1,dynk_iData(dynk_funcs(funNum,3)+2))
    ! Save our current seeds and load old seeds
    call recuut(dynk_iData(dynk_funcs(funNum,3)+3),dynk_iData(dynk_funcs(funNum,3)+4))
    call recuin(tmpseed1,tmpseed2)
    ! Change to mu, sigma
    retval = dynk_fData(dynk_funcs(funNum,4)) + dynk_fData(dynk_funcs(funNum,4)+1)*ranecu_rvec(1)

  case(7) ! RANDU
    ! Save old seeds and load our current seeds
    call recuut(tmpseed1,tmpseed2)
    call recuin(dynk_iData(dynk_funcs(funNum,3)+2),dynk_iData(dynk_funcs(funNum,3)+3))
    ! Run generator for 1 value with mcut=-1
    call ranecu( ranecu_rvec, 1, -1 )
    ! Save our current seeds and load old seeds
    call recuut(dynk_iData(dynk_funcs(funNum,3)+2),dynk_iData(dynk_funcs(funNum,3)+3))
    call recuin(tmpseed1,tmpseed2)
    retval = ranecu_rvec(1)

  case (8) ! RANDON
    ! Save old seeds and load our current seeds
    call recuut(tmpseed1,tmpseed2)
    call recuin(dynk_iData(dynk_funcs(funNum,3)+2),dynk_iData(dynk_funcs(funNum,3)+3))
    ! Run generator for 1 value with mcut=-1
    call ranecu( ranecu_rvec, 1, -1 )
    ! Save our current seeds and load old seeds
    call recuut(dynk_iData(dynk_funcs(funNum,3)+2),dynk_iData(dynk_funcs(funNum,3)+3))
    call recuin(tmpseed1,tmpseed2)
    ! routine for switching element (orginially the electron lens) ON or OFF
    ! when random value is less than P, set ON, else OFF
    if (ranecu_rvec(1) .lt. dynk_fData(dynk_funcs(funNum,4))) then
      retval = 1.0
    else
      retval = 0.0
    end if

  case(10) ! FIR
    foff = dynk_funcs(funNum,3)
    ! Shift storage 1 back
    do ii=dynk_funcs(funNum,4)-1,0,-1
      jj = ii*3
      dynk_fData(foff+jj+4) = dynk_fData(foff+jj+1)
    end do
    ! Evaluate the next input function
    dynk_fData(foff+1) = dynk_computeFUN(dynk_funcs(funNum,5),turn)
    ! Compute the filtered value
    retval = 0.0
    do ii=0,dynk_funcs(funNum,4)
      jj = ii*3
      retval = retval + dynk_fData(foff+jj)*dynk_fData(foff+jj+1)
    end do

  case(11) ! IIR
    foff = dynk_funcs(funNum,3)
    ! Shift storage 1 back
    do ii=dynk_funcs(funNum,4)-1,0,-1
      jj = ii*6
      dynk_fData(foff+jj+7) = dynk_fData(foff+jj+1)
      dynk_fData(foff+jj+10) = dynk_fData(foff+jj+4)
    end do
    ! Evaluate the next input function
    dynk_fData(foff+1) = dynk_computeFUN(dynk_funcs(funNum,5),turn)
    dynk_fData(foff+4) = 0.0
    ! Compute the filtered value
    retval = 0.0
    do ii=0,dynk_funcs(funNum,4)
      jj = ii*6
      retval = retval + &
        dynk_fData(foff+jj  ) * dynk_fData(foff+jj+1) + &
        dynk_fData(foff+jj+3) * dynk_fData(foff+jj+4)
    end do
    ! To be shifted at the next evaluation
    dynk_fData(foff+4) = retval

  case(20) ! ADD
    retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) + dynk_computeFUN(dynk_funcs(funNum,4),turn)

  case(21) ! SUB
    retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) - dynk_computeFUN(dynk_funcs(funNum,4),turn)

  case(22) ! MUL
    retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) * dynk_computeFUN(dynk_funcs(funNum,4),turn)

  case(23) ! DIV
    retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) / dynk_computeFUN(dynk_funcs(funNum,4),turn)

  case(24) ! POW
    retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) ** dynk_computeFUN(dynk_funcs(funNum,4),turn)

  case(30) ! MINUS
    retval = (-1)*dynk_computeFUN(dynk_funcs(funNum,3),turn)

  case(31) ! SQRT
    retval = sqrt(dynk_computeFUN(dynk_funcs(funNum,3),turn))

  case(32) ! SIN
    retval = sin_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

  case(33) ! COS
    retval = cos_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

  case(34) ! LOG
    retval = log_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

  case(35) ! LOG10
    retval = log10_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

  case(36) ! EXP
    retval = exp_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

  case(37) ! ABS
    retval = abs(dynk_computeFUN(dynk_funcs(funNum,3),turn))

  case(40) ! CONST
    retval = dynk_fData(dynk_funcs(funNum,3))

  case(41) ! TURN
    retval = turn

  case(42) ! LIN
    retval = turn*dynk_fData(dynk_funcs(funNum,3)) + dynk_fData(dynk_funcs(funNum,3)+1)

  case(43) ! LINSEG
    filelin_start    = dynk_funcs(funNum,3)
    filelin_xypoints = 2
    ! Pass the correct array views/sections to lininterp
    retval = lininterp( real(turn,fPrec), &
             dynk_fData(filelin_start:filelin_start+1), &
             dynk_fData(filelin_start+2:filelin_xypoints+3), &
             filelin_xypoints )

  case(44,45) ! QUAD/QUADSEG
    retval = (turn**2)*dynk_fData(dynk_funcs(funNum,3))   + ( &
              turn*dynk_fData(dynk_funcs(funNum,3)+1) + &
                    dynk_fData(dynk_funcs(funNum,3)+2) )

  case(60) ! SINF
    retval = dynk_fData(dynk_funcs(funNum,3)) &
            * sin_mb( dynk_fData(dynk_funcs(funNum,3)+1) * turn  &
                    + dynk_fData(dynk_funcs(funNum,3)+2) )

  case(61) ! COSF
    retval = dynk_fData(dynk_funcs(funNum,3)) &
            * cos_mb( dynk_fData(dynk_funcs(funNum,3)+1) * turn  &
                    + dynk_fData(dynk_funcs(funNum,3)+2) )

  case(62) ! COSF_RIPP
    retval = dynk_fData(dynk_funcs(funNum,3)) &
            * cos_mb( (two*pi)*real(turn-1,fPrec)/dynk_fData(dynk_funcs(funNum,3)+1) &
                    + dynk_fData(dynk_funcs(funNum,3)+2) )

  case(80) ! PELP
    foff = dynk_funcs(funNum,3)
    if(turn <= dynk_fData(foff)) then ! <= tinj
      ! Constant Iinj
      retval = dynk_fData(foff+5)
    else if(turn <= dynk_fData(foff+1)) then ! <= te
      ! Parabola (accelerate)
      retval = (dynk_fData(foff+6)*(turn-dynk_fData(foff))**2)/2.0 + dynk_fData(foff+5)
    else if(turn <= dynk_fData(foff+2)) then ! <= t1
      ! Exponential
      retval = dynk_fData(foff+7)*exp_mb(dynk_fData(foff+8)*turn)
    else if(turn <= dynk_fData(foff+3)) then ! <= td
      ! Linear (max ramp rate)
      retval = dynk_fData(foff+10)*(turn-dynk_fData(foff+2)) + dynk_fData(foff+9)
    else if(turn <= dynk_fData(foff+4)) then ! <= tnom
      ! Parabola (decelerate)
      retval = -((dynk_fData(foff+11) * (dynk_fData(foff+4)-turn)**2) ) / 2.0 + dynk_fData(foff+12)
    else ! > tnom
      ! Constant Inom
      retval = dynk_fData(foff+12)
    end if

  case (81)                                                         ! ONOFF
      ii=mod(turn-1,dynk_funcs(funNum,4))
      if (ii .lt. dynk_funcs(funNum,3)) then
          retval = 1.0
      else
          retval = 0.0
      end if

  case default
    write(lout,"(2(a,i0))") "DYNK> ERROR computeFUN() funNum = ",funNum," turn = ",turn
    write(lout,"(a,i0)")    "DYNK>       Unknown function type ",dynk_funcs(funNum,2)
    if(dynk_debug) call dynk_dumpdata
    call prror(-1)

  end select

end function dynk_computeFUN

! ================================================================================================ !
!  A. Santamaria, K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-12
!  - Set the value of the element's attribute
! ================================================================================================ !
subroutine dynk_setvalue(element_name, att_name, newValue)

  use scatter, only : scatter_ELEM_scale, scatter_elemPointer
  use crcoall
  use mod_common
  use mod_commont
  use mod_commonmn
  use mod_particles

  use elens
  use parbeam, only : beam_expflag
  implicit none

  character(mStrLen), intent(in) :: element_name, att_name
  real(kind=fPrec),   intent(in) :: newValue

  ! Temp variables
  integer el_type, ii, j, orderMult, im,k, range

  ! Original energies before energy update
  real(kind=fPrec) e0fo, e0o, r0a, r0

  ! For sanity check
  logical ldoubleElement, iErr
  ldoubleElement = .false.
  iErr = .false.

  if(dynk_debug) then
    write(lout,"(a,e16.9)") "DYNK> DEBUG setvalue Element_name = '"//trim(element_name)//"', "//&
      "att_name = '"//trim(att_name)//"', newValue = ", newValue
  end if

  ! Here comes the logic for setting the value of the attribute for all instances of the element.

  ! Special non-physical elements
  if(element_name == "GLOBAL-VARS") then
    if(att_name == "E0" ) then
      ! Modify the reference particle
      call part_updateRefEnergy(newValue)
      ! Modify energy-dependent element parameters
      call eLensThetas
    end if
    ldoubleElement = .true.
  end if

  ! Normal SINGLE ELEMENTs
  do ii=1,il
    ! TODO: Here one could find the right ii in dynk_pretrack,
    ! and then avoid this loop / string-comparison
    if(element_name == bez(ii)) then ! name found
      el_type=kz(ii)      ! type found

      if(ldoubleElement) then ! Sanity check
        write(lout,"(a)") "DYNK> ERROR Two elements with the same BEZ"
        call prror(-1)
      end if
      ldoubleElement = .true.

      select case(abs(el_type))

      case(1,2,3,4,5,6,7,8,9,10)
        ! horizontal bending kick, quadrupole kick, sextupole kick, octupole kick, decapole kick,
        ! dodecapole kick, 14th pole kick, 16th pole kick, 18th pole kick,20th pole kick
        if(att_name == "average_ms") then
          ed(ii) = newValue
        else
          goto 100 ! ERROR
        end if
        call initialize_element(ii, .false.)

      case(11)
        im = irm(ii)
        if(att_name == "scaleall") then
          scalemu(im) = newValue
        else if(att_name(1:1) == "a" .or. att_name(1:1) == "b") then
          if(len_trim(att_name) == 5) then
            range = 3
            call chr_cast(att_name(2:2), orderMult, iErr)
          else if(len_trim(att_name) == 6) then
            call chr_cast(att_name(2:3), orderMult, iErr)
            range = 4
          else
            goto 100
          end if
          if(iErr) goto 100
          if(nmu(ii) < orderMult) then
            nmu(ii) = orderMult
          end if

          r0  = r00(im)
          r0a = one
          do k=2,orderMult
            r0a = r0a*r0
          end do
          if(att_name(1:1) == "a" .and. att_name(range:range+3) == "rms") then
            aka(im,orderMult) = newValue*benkc(im)/r0a
          else if(att_name(1:1) == "b" .and. att_name(range:range+3) == "rms") then
            bka(im,orderMult) = newValue*benkc(im)/r0a
          else if(att_name(1:1) == "a" .and. att_name(range:range+3) == "str") then
            ak0(im,orderMult) = newValue*benkc(im)/r0a
          else if(att_name(1:1) == "b" .and. att_name(range:range+3) == "str") then
            bk0(im,orderMult) = newValue*benkc(im)/r0a
          else
            goto 100 ! ERROR
          end if
        end if
        call initialize_element(ii, .false.)

      case(12)
        if(att_name == "voltage") then ! [MV]
          ed(ii) = newValue
        else if(att_name == "harmonic") then
          ek(ii) = newValue
          el(ii) = dynk_elemData(ii,3) ! Need to reset el before calling initialize_element()
          call initialize_element(ii, .false.)
        else if(att_name == "lag_angle") then ! [deg]
          el(ii) = newValue
          ! Note: el is set to 0 in initialize_element and in daten.
          ! Calling initialize element on a cavity without setting el
          ! will set phasc = 0!
          call initialize_element(ii, .false.)
        else
          goto 100 ! ERROR
        end if

      ! Not yet supported : AC dipole (16)

      case(20)
        if(beam_expflag == 1) then
          if(att_name == "h-sep") then ! [mm]
            parbe(ii,5)  = newValue
          else if(att_name == "v-sep") then ! [mm]
            parbe(ii,6)  = newValue
          else if(att_name == "4dSxx") then ! strong I think
            parbe(ii,1)  = newValue
          else if(att_name == "4dSyy") then
            parbe(ii,3)  = newValue
          else if(att_name == "4dSxy") then
            parbe(ii,13) = newValue
          else if(att_name == "strength") then
            ptnfac(ii)   = newValue
            parbe(ii,4)  = (((-one*crad)*ptnfac(ii))*half)*c1m6
          else if(att_name == "Sxx") then
            parbe(ii,7)   = newValue
          else if(att_name == "Sxxp") then
            parbe(ii,8)   = newValue
          else if(att_name == "Sxpxp") then
            parbe(ii,9) = newValue
          else if(att_name == "Syy") then
            parbe(ii,10) = newValue
          else if(att_name == "Syyp") then
            parbe(ii,11) = newValue
          else if(att_name == "Sypyp") then
            parbe(ii,12) = newValue
          else if(att_name == "Sxy") then
            parbe(ii,13) = newValue
          else if(att_name == "Sxyp") then
            parbe(ii,14) = newValue
          else if(att_name == "Sxpy") then
            parbe(ii,15) = newValue
          else if(att_name == "Sxpyp") then
            parbe(ii,16) = newValue
          else
            goto 100
          end if
          call initialize_element(ii, .false.)
        else
          goto 102
        end if

      case(23,26,27,28)
        ! crab cavity, cc mult. kick order 2,3 and 4
        if(att_name == "voltage") then ! [MV]
          ed(ii) = newValue
        else if(att_name == "frequency") then ! [MHz]
          ek(ii) = newValue
        else if(att_name == "phase") then ! [rad]
          ! Note: el is set to 0 in initialize_element and in daten.
          ! Calling initialize element on a crab without setting el
          ! will set crabph = 0!
          el(ii) = newValue
          call initialize_element(ii, .false.)
        else
          goto 100 ! ERROR
        end if

      case(29) ! Electron lens
        if(att_name == "theta_r2") then ! [mrad]
          elens_theta_r2(ielens(ii)) = newValue
        elseif(att_name == "elens_I") then ! [A]
          elens_I(ielens(ii)) = newValue
          call eLensTheta(ielens(ii))
        elseif(att_name == "elens_Ek") then ! [keV]
          elens_Ek(ielens(ii)) = newValue
          call eLensTheta(ielens(ii))
        else
          goto 100 ! ERROR
        end if

      case(40) ! Scatter
        if(att_name == "scaling") then
          scatter_ELEM_scale(scatter_elemPointer(ii)) = newValue
        else
          goto 100 ! ERROR
        end if

      case default
        write(lout,"(a,i0,a)") "DYNK> ERROR setValu Unsupported element type ",el_type," element name = '"//trim(element_name)//"'"
        call prror(-1)

      end select
    end if
  end do

  ! Sanity check
  if(.not.ldoubleElement) then
    goto 101
  end if

  return

  ! Error handlers
100 continue
  write(lout,"(a,i0)")"DYNK> ERROR setValue Attribute '"//trim(att_name)//"' does not exist for type = ",el_type
  call prror(-1)

101 continue
  write(lout,"(a)") "DYNK> ERROR setValue The element named '"//trim(element_name)//"' was not found."
  call prror(-1)

102 continue
  write(lout,"(a)") "DYNK> ERROR setValue Only Beam-beam expert mode is supported for DYNK"
  call prror(-1)

end subroutine dynk_setvalue

! ================================================================================================ !
!  A. Santamaria, K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-12
!  - Returns the original value currently set by an element.
! ================================================================================================ !
real(kind=fPrec) function dynk_getvalue(element_name, att_name)

  use scatter, only : scatter_ELEM_scale, scatter_elemPointer
  use crcoall
  use mod_common
  use mod_commont
  use mod_commonmn
  use elens
  use parbeam, only : beam_expflag

  implicit none

  character(mStrLen), intent(in) :: element_name, att_name

  integer el_type, ii, orderMult, im, range

  logical ldoubleElement, iErr
  ldoubleElement = .false.  ! For sanity check
  iErr = .false.
  if(dynk_debug) then
    write(lout,"(a)") "DYNK> DEBUG In getValue, element_name = '"//trim(element_name)//"'"//&
      ", att_name = '"//trim(att_name)//"'"
  end if

  ! Special non-physical elements
  if(element_name == "GLOBAL-VARS") then
    if(att_name == "E0" ) then
      ! Return the energy
      dynk_getvalue = e0
    end if
    ldoubleElement = .true.
  end if

  ! Normal SINGLE ELEMENTs
  do ii=1,il
    ! TODO: Here one could find the right ii in dynk_pretrack,
    ! and then avoid this loop / string-comparison
    if(element_name == bez(ii)) then ! name found
      el_type=kz(ii)
      if(ldoubleElement) then
        write (lout,"(a)") "DYNK> ERROR Two elements with the same BEZ"
        call prror(-1)
      end if
      ldoubleElement = .true.

      select case(abs(el_type))

      case(1,2,3,4,5,6,7,8,9,10) ! Nonlinear elements
        if(att_name == "average_ms") then
          dynk_getvalue = ed(ii)
        else
          goto 100 ! ERROR
        end if


      case(11)
        im=irm(ii)
        if(att_name=="scaleall") then
          dynk_getvalue = scalemu(im)
        else if(att_name(1:1) == "a" .or. att_name(1:1) == "b") then
          if(LEN_TRIM(att_name) .eq. 5) then
             range = 3
            call chr_cast(att_name(2:2), orderMult, iErr)
          else if(LEN_TRIM(att_name) .eq. 6) then
            call chr_cast(att_name(2:3), orderMult, iErr)
            range = 4
          else
            goto 100
          endif
          if(iErr) goto 100

          if(att_name(1:1) == "a" .and. att_name(range:range+3) == "rms") then
            dynk_getvalue = aka(im,orderMult)
          else if(att_name(1:1) == "b" .and. att_name(range:range+3) == "rms") then
            dynk_getvalue = bka(im,orderMult)
          else if(att_name(1:1) == "a" .and. att_name(range:range+3) == "str") then
            dynk_getvalue = ak0(im,orderMult)
          else if(att_name(1:1) == "b" .and. att_name(range:range+3) == "str") then
            dynk_getvalue = bk0(im,orderMult)
          else
            goto 100 ! ERROR
          endif
        endif

      case(12) ! Cavities
        if(att_name == "voltage"  ) then ! MV
          dynk_getvalue = ed(ii)
        else if(att_name == "harmonic" ) then ! harmonic number
          dynk_getvalue = ek(ii)
        else if(att_name == "lag_angle") then ! [deg]
          dynk_getvalue = dynk_elemData(ii,3)
        else
          goto 100 ! ERROR
        end if

      ! Not yet supported : AC dipole (16)

      case(20)
        if(beam_expflag.eq.1) then
          if (att_name == "h-sep") then ! [mm]
            dynk_getvalue = parbe(ii,5)
          else if (att_name == "v-sep") then ! [mm]
            dynk_getvalue = parbe(ii,6)
          else if (att_name=="4dSxx") then !strong I think
            dynk_getvalue = parbe(ii,1)
          else if (att_name=="4dSyy") then
            dynk_getvalue = parbe(ii,3)
          else if (att_name=="4dSxy") then
            dynk_getvalue = parbe(ii,13)
          else if (att_name == "strength") then !
            dynk_getvalue = ptnfac(ii)
          else if(att_name == "Sxx") then
            dynk_getvalue = parbe(ii,7)
          else if(att_name == "Sxxp") then
            dynk_getvalue = parbe(ii,8)
          else if(att_name == "Sxpxp") then
            dynk_getvalue = parbe(ii,9)
          else if(att_name == "Syy") then
            dynk_getvalue = parbe(ii,10)
          else if(att_name == "Syyp") then
            dynk_getvalue = parbe(ii,11)
          else if(att_name == "Sypyp") then
            dynk_getvalue = parbe(ii,12)
          else if(att_name == "Sxy") then
            dynk_getvalue = parbe(ii,13)
          else if(att_name == "Sxyp") then
            dynk_getvalue = parbe(ii,14)
          else if(att_name == "Sxpy") then
            dynk_getvalue = parbe(ii,15)
          else if(att_name == "Sxpyp") then
            dynk_getvalue = parbe(ii,16)
          else
            go to 100
          endif
        else
          go to 102
        end if

      case(23,26,27,28) ! crab cavity, cc mult. kick order 2, 3 and 4
        if(att_name == "voltage") then ! [MV]
          dynk_getvalue = ed(ii)
        else if(att_name == "frequency") then ! [MHz]
          dynk_getvalue = ek(ii)
        else if(att_name == "phase") then ! [rad]
          if(abs(el_type) == 23) then
              dynk_getvalue = crabph(ii)
          else if(abs(el_type) == 26) then
              dynk_getvalue = crabph2(ii)
          else if(abs(el_type) == .27) then
              dynk_getvalue = crabph3(ii)
          else if(abs(el_type) == 28) then
              dynk_getvalue = crabph4(ii)
          end if
        else
          goto 100 ! ERROR
        end if

      case(29) ! Electron lens
        if(att_name == "theta_r2") then ! [mrad]
          dynk_getvalue = elens_theta_r2(ielens(ii))
        elseif(att_name == "elens_I") then ! [A]
          dynk_getvalue = elens_I(ielens(ii))
        elseif(att_name == "elens_Ek") then ! [keV]
          dynk_getvalue = elens_Ek(ielens(ii))
        else
          goto 100 ! ERROR
        end if

      case(40) ! Scatter
        if(att_name == "scaling") then
          dynk_getvalue = scatter_ELEM_scale(scatter_elemPointer(ii))
        else
          goto 100 ! ERROR
        end if

      end select
    end if ! bez
  end do

  if(dynk_debug) then
    write(lout,"(a,e16.9)") "DYNK> DEBUG getValue, returning = ",dynk_getvalue
  end if

  return

  ! Error handlers
100 continue
  write(lout,"(a,i0,a)") "DYNK> ERROR getValueUnknown attribute '"//trim(att_name)//"'"//&
    " for type ",el_type," name '"//trim(bez(ii))//"'"
  call prror(-1)

102 continue
  write(lout,"(a)") "DYNK> ERROR  --- Only Beam-beam expert mode is supported for DYNK"
  call prror(-1)

end function dynk_getvalue

! ================================================================================================ !
!  K. Sjobak, BE-ABP-HSS,
!  Last modified: 2018-08-12
!  - Indicates whether a structure element is in use by DYNK
! ================================================================================================ !
logical function dynk_isused(i)

  use crcoall
  use mod_common

  implicit none

  integer, intent(in) :: i
  integer ix,k
  character(mStrLen) element_name

  ! Sanity check
  if(i > iu .or. i <= 0) then
    write(lout,"(a,i0,a)") "DYNK> ERROR isused: i=",i," out of range"
    call prror(-1)
  end if
  ix = ic(i)-nblo
  if(i <= 0) then
    write(lout,"(a,i0,a)") "DYNK> ERROR isused: ix-nblo = ",ix," is a block?"
    call prror(-1)
  end if

  do k=1,dynk_nSets
    element_name = trim(dynk_cSets(k,1))
    if(bez(ix) == element_name) then
      dynk_isused = .true.
      if(dynk_debug) then
        write(lout,"(a)") "DYNK> DEBUG dynk_isused = TRUE, bez='"//bez(ix)//"', element_name_stripped='"//element_name//"'"
      end if
      return
    end if
  end do

  if(dynk_debug) then
    write(lout,"(a)") "DYNK> DEBUG dynk_isused = FALSE, bez='"//bez(ix)//"'"
  end if

  dynk_isused = .false.
  return

end function dynk_isused

! =================================================================================================
!  A. Mereghetti, for the FLUKA Team
!  Last Modified: 2018-04-17
!  Close units for logging dynks
! =================================================================================================
subroutine dynk_closeFiles

  implicit none

  integer i
  logical isOpen

  if(.not. dynk_enabled) return

  do i=1,dynk_nFuncs
    if (dynk_funcs(i,2) == 3) then ! PIPE FUN
      ! InPipe
      inquire(unit=dynk_iData(dynk_funcs(i,3)), opened=isOpen)
      if(isOpen) close(dynk_iData(dynk_funcs(i,3)))

      ! OutPipe
      inquire(unit=dynk_iData(dynk_funcs(i,3)+1), opened=isOpen)
      if(isOpen) then
        write(dynk_iData(dynk_funcs(i,3))+1,"(a)") "CLOSEUNITS"
        close(dynk_iData(dynk_funcs(i,3))+1)
      end if
    end if
  end do

end subroutine dynk_closeFiles

! ================================================================================================ !
!  BEGIN Checkpoint Restart
! ================================================================================================ !
#ifdef CR

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen BE-ABP-HSS
!  Last modified: 2018-05-28
!  - Called from CRCHECK; reads the _cr arrays for scatter from file.
!  - Sets readerr=.true. if something goes wrong while reading.
! ================================================================================================ !
subroutine dynk_crcheck_readdata(fileunit,readerr)

  use crcoall

  implicit none

  integer, intent(in)  :: fileunit
  logical, intent(out) :: readerr

  integer j, stat

  read(fileunit,err=100,end=100) dynk_filePosCR, dynk_niData_cr, dynk_nfData_cr, dynk_ncData_cr

  call alloc(dynk_iData_cr,        dynk_niData_cr, 0,         "dynk_iData_cr")
  call alloc(dynk_fData_cr,        dynk_nfData_cr, zero,      "dynk_fData_cr")
  call alloc(dynk_cData_cr,mStrLen,dynk_ncData_cr, str_dSpace,"dynk_cData_cr")

  read(fileunit,err=100,end=100) &
    (dynk_iData_cr(j),j=1,dynk_niData_cr), (dynk_fData_cr(j),j=1,dynk_nfData_cr), &
    (dynk_cData_cr(j),j=1,dynk_ncData_cr), (dynk_fSets_cr(j),j=1,dynk_maxSets)

  readerr=.false.
  return

100 continue

  write(lout,"(a,i0)") "READERR in scatter_crcheck; fileunit=",fileunit
  write(93,*)          "READERR in scatter_crcheck; fileunit=",fileunit
  readerr=.true.

end subroutine dynk_crcheck_readdata

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen BE-ABP-HSS
!  Last modified: 25-09-2017
!  - Called from CRCHECK; resets the position of dynksets.txt
! ================================================================================================ !
subroutine dynk_crcheck_positionFiles

  use crcoall
  implicit none

  logical isOpen
  integer ierro
#ifdef BOINC
  character(len=256) filename
#endif
  integer j
  character(len=1024) arecord

  inquire(unit=dynk_fileUnit, opened=isOpen)
  if (isOpen) then
    write(93,*) "SIXTRACR CRCHECK FAILED while repositioning 'dynksets.dat'"
    write(93,*) "Could not open file!"
    endfile (93,iostat=ierro)
    backspace (93,iostat=ierro)

    write(lout,"(a)") "SIXTRACR> CRCHECK failure positioning 'dynksets.dat'"
    call prror(-1)
  end if

  if (dynk_filePosCR /= -1) then
#ifdef BOINC
    call boincrf("dynksets.dat",filename)
    open(unit=dynk_fileUnit,file=filename,status="old",action="readwrite",err=110)
#else
    open(unit=dynk_fileUnit,file='dynksets.dat',status="old",action="readwrite",err=110)
#endif
    dynk_filePos = 0     ! Start counting lines at 0, not -1
    do j=1,dynk_filePosCR
      read(dynk_fileUnit,'(a1024)',end=110,err=110,iostat=ierro) arecord
      dynk_filePos=dynk_filePos+1
    end do

    endfile(dynk_fileUnit,iostat=ierro)
    close(dynk_fileUnit)
#ifdef BOINC
    call boincrf("dynksets.dat",filename)
    open(unit=dynk_fileUnit,file=filename,status="old",position='append',action="write")
#else
    open(unit=dynk_fileUnit,file="dynksets.dat",status="old",position='append',action="write")
#endif

    write(93,*) "SIXTRACR CRCHECK sucessfully repositioned 'dynksets.dat', "// &
                "dynk_filePos=",dynk_filePos, "dynk_filePosCR=",dynk_filePosCR
    endfile (93,iostat=ierro)
    backspace (93,iostat=ierro)
  else
    write(93,*) "SIXTRACR CRCHECK did not attempt repositioning "// &
                "of dynksets.dat, dynk_filePosCR=",dynk_filePosCR
    write(93,*) "If anything has been written to the file, "// &
                "it will be correctly truncated in dynk_apply on the first turn."
    endfile (93,iostat=ierro)
    backspace (93,iostat=ierro)
  end if

  return

110 continue
  write(93,*) "SIXTRACR CRCHECK *** ERROR *** reading 'dynksets.dat', iostat=",ierro
  write(93,*) "dynk_filePos=",dynk_filePos," dynk_filePosCR=",dynk_filePosCR
  endfile   (93,iostat=ierro)
  backspace (93,iostat=ierro)
  write(lout,"(a)") "SIXTRACR> CRCHECK failure positioning 'dynksets.dat'"
  call prror(-1)

end subroutine dynk_crcheck_positionFiles

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen BE-ABP-HSS
!  Last modified: 2018-05-28
!  - Called from CRPOINT; write checkpoint data to fort.95/96
! ================================================================================================ !
subroutine dynk_crpoint(fileunit,fileerror,ierro)

  implicit none

  integer, intent(in)    :: fileunit
  logical, intent(inout) :: fileerror
  integer, intent(inout) :: ierro

  integer j

  !Note: dynk_fSets_cr is set in global `crpoint` routine, in order to avoid
  ! that it is filled twice (requiring loop over all dynk_fsets_unique and call to dynk_getvalue)
  write(fileunit,err=100,iostat=ierro) dynk_filePos, dynk_niData, dynk_nfData, dynk_ncData
  write(fileunit,err=100,iostat=ierro) &
      (dynk_iData(j),j=1,dynk_niData), (dynk_fData(j),j=1,dynk_nfData), &
      (dynk_cData(j),j=1,dynk_ncData), (dynk_fSets_cr(j),j=1,dynk_maxSets)
  endfile (fileunit,iostat=ierro)
  backspace (fileunit,iostat=ierro)

  return

100 continue
    fileerror=.true.
end subroutine dynk_crpoint

! ================================================================================================ !
!  K. Sjobak, V.K. Berglyd Olsen BE-ABP-HSS
!  Last modified: 2018-05-28
!  - Called from CRSTART; copies the _cr arrays into the normal arrays used during tracking in
!    order to recreate the state of the SCATTER block at the time of the checkpoint.
! ================================================================================================ !
subroutine dynk_crstart

  implicit none

  integer j

  call dynk_checkspace(dynk_niData_cr-dynk_niData, dynk_nfData_cr-dynk_nfData, dynk_ncData_cr-dynk_ncData)

  dynk_niData = dynk_niData_cr
  dynk_nfData = dynk_nfData_cr
  dynk_ncData = dynk_ncData_cr

  dynk_iData(1:dynk_niData) = dynk_iData_cr(1:dynk_niData)
  dynk_fData(1:dynk_nfData) = dynk_fData_cr(1:dynk_nfData)
  dynk_cData(1:dynk_ncData) = dynk_cData_cr(1:dynk_ncData)

  call dealloc(dynk_iData_cr,        "dynk_iData_cr")
  call dealloc(dynk_fData_cr,        "dynk_fData_cr")
  call dealloc(dynk_cData_cr,mStrLen,"dynk_cData_cr")

  do j=1,dynk_nSets_unique
    ! It is OK to write to lout from here
    call dynk_setvalue(dynk_cSets_unique(j,1),dynk_cSets_unique(j,2),dynk_fSets_cr(j))
  end do

end subroutine dynk_crstart

! ================================================================================================ !
!  END Checkpoint Restart
! ================================================================================================ !
#endif

! ================================================================================================ !
end module dynk
