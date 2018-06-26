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
  use parpro, only : nele, mNameLen
  use mod_alloc
  use string_tools

  implicit none

  ! General-purpose variables
  logical, save :: ldynk            ! dynamic kick requested, i.e. DYNK input bloc issued in the fort.3 file
  logical, save :: ldynkdebug       ! print debug messages in main output
  logical, save :: ldynkfiledisable ! Disable writing dynksets.dat?
  integer, save :: ldynkfileunit    ! The file unit for dynksets.dat
  integer, save :: ldynkFUNfileunit ! File unit for parseFUN files

  ! Store the FUN statements
  integer, save :: dynk_maxFuncs
  integer, save :: dynk_maxiData
  integer, save :: dynk_maxfData
  integer, save :: dynk_maxcData

  ! 1 row/FUN, cols are:
  ! (1) = function name in fort.3 (points within dynk_cData),
  ! (2) = indicates function type
  ! (3,4,5) = arguments (often pointing within other arrays {i|f|c}expr_dynk)
  integer, allocatable, save :: dynk_funcs(:,:) ! (dynk_maxFuncs,5)

  ! Data for DYNK FUNs
  integer,          allocatable, save :: dynk_iData(:)
  real(kind=fPrec), allocatable, save :: dynk_fData(:)
  character(len=:), allocatable, save :: dynk_cData(:)

  ! Number of used positions in arrays
  integer, save :: dynk_nFuncs
  integer, save :: dynk_niData
  integer, save :: dynk_nfData
  integer, save :: dynk_ncData

  ! Store the SET statements
  integer, parameter :: dynk_maxSets = 200

  ! 1 row/SET, cols are:
  ! (1) = function index (points within dynk_funcs)
  ! (2) = first turn num. where it is active
  ! (3) =  last turn num. where it is active
  ! (4) = Turn shift - number added to turn before evaluating the FUN
  integer, save :: dynk_sets(dynk_maxSets, 4)

  ! 1 row/SET (same ordering as dynk_sets), cols are:
  ! (1) element name
  ! (2) attribute name
  character(mStrLen), save :: dynk_cSets(dynk_maxSets,2)

  ! Number of used positions in arrays
  integer, save :: dynk_nSets

  ! Similar to dynk_cSets, but only one entry per elem/attr
  character(mStrLen), save :: dynk_cSets_unique(dynk_maxSets,2)

  ! Store original value from dynk
  real(kind=fPrec), save :: dynk_fSets_orig(dynk_maxSets)

  ! Number of used positions in arrays
  integer, save :: dynk_nSets_unique

  ! Some elements (multipoles) overwrites the general settings info when initialized.
  ! Store this information on the side.
  ! Also used by setvalue and getvalue
  integer, allocatable, save :: dynk_izuIndex(:) !(nele
  real(kind=fPrec), allocatable, save :: dynk_elemdata(:,:) !(nele,3)

#ifdef CR
  ! Block with data/fields needed for checkpoint/restart of DYNK

  ! Number of records written to dynkfile (dynksets.dat)
  integer, save :: dynkfilepos
  integer, save :: dynkfilepos_cr

  ! Data for DYNK FUNs
  integer,          allocatable, save :: dynk_iData_cr(:)
  real(kind=fPrec), allocatable, save :: dynk_fData_cr(:)
  character(len=:), allocatable, save :: dynk_cData_cr(:)

  ! Number of used positions in arrays
  integer, save :: dynk_niData_cr
  integer, save :: dynk_nfData_cr
  integer, save :: dynk_ncData_cr

  ! Store current settings from dynk
  real(kind=fPrec), save :: dynk_fSets_cr(dynk_maxSets)

#endif

contains

subroutine dynk_allocate_arrays
  call alloc(dynk_izuIndex,nele,0,'dynk_izuIndex')
  call alloc(dynk_elemdata,nele,3,zero,'dynk_elemdata')
end subroutine dynk_allocate_arrays

subroutine dynk_expand_arrays(nele_new)
  integer, intent(in) :: nele_new
  call alloc(dynk_izuIndex,nele_new,0,'dynk_izuIndex')
  call alloc(dynk_elemdata,nele_new,3,zero,'dynk_elemdata')
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

  call alloc(dynk_iData,dynk_maxiData,0,"dynk_iData")
  call alloc(dynk_fData,dynk_maxfData,zero,"dynk_fData")
  call alloc(dynk_cData,mStrLen,dynk_maxcData,str_dZeros,"dynk_cData")
  call alloc(dynk_funcs,dynk_maxFuncs,5,0,"dynk_funcs")

  ! Set file units for I/O files
  call funit_requestUnit("dynksets.dat",ldynkfileunit)
  call funit_requestUnit("dynk_parseFUN_IO",ldynkFUNfileunit)

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

  character gFields(str_maxFields)*(mStrLen)
  integer   nFields
  integer   lFields(str_maxFields)
  logical   eFields

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lout,"(a)") "DYNK> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) then
    if(ldynkdebug) then
      write (lout,"(a,i3,a)") "DYNK> DEBUG Input line len=",len(inLine),": '"//chr_strip(inLine)//"'."
      write (lout,"(a)")      "DYNK> DEBUG  * No fields found."
    end if
    return
  end if

  ! Report if debugging is ON
  if(ldynkdebug) then
    write (lout,"(a,i3,a)")  "DYNK> DEBUG Input line len=",len(inLine),": '"//chr_strip(inLine)//"'."
    write (lout,"(a,i3,a)") ("DYNK> DEBUG  * Field(",i,") = '"//chr_trimZero(lnSplit(i))//"'",i=1,nSplit)
  end if

  select case(chr_trimZero(lnSplit(1)))

  case("DEBUG")
    ldynkdebug = .true.
    write(lout,"(a)") "DYNK> DYNK block debugging is ON."

  case("NOFILE")
    ldynkfiledisable = .true.
    write (lout,*) "DYNK> Disabled writing dynksets.dat"

  case("FUN")
    call getfields_split(inLine,gFields,lFields,nFields,eFields)
    call dynk_parseFUN(gFields,lFields,nFields)

  case("SET")
    call getfields_split(inLine,gFields,lFields,nFields,eFields)
    call dynk_parseSET(gFields,lFields,nFields)

  case default
    write(lout,"(a)") "DYNK> ERROR Unrecognised statement '"//chr_trimZero(lnSplit(1))//"'."
    iErr = .true.
    return

  end select

end subroutine dynk_parseInputLine

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-28
!  Parse FUN lines in the fort.3 input file.
! =================================================================================================
subroutine dynk_parseFUN(gFields, lFields, nFields)

  use crcoall
  implicit none

  character, intent(in) :: gFields(str_maxFields)*(mStrLen)
  integer,   intent(in) :: nFields
  integer,   intent(in) :: lFields(str_maxFields)

  ! Temp variables
  integer ii, stat, t, nlines
  real(kind=fPrec) x,y,z,u                                  ! FILE, FILELIN, FIR/IIR
  real(kind=fPrec) x1,x2,y1,y2,deriv                        ! LINSEG, QUADSEG,
  real(kind=fPrec) tinj,Iinj,Inom,A,D,R,te                  ! PELP (input)
  real(kind=fPrec) derivI_te,I_te,bexp,aexp, t1,I1, td,tnom ! PELP (calc)

  logical isFIR ! FIR/IIR
  logical lopen

#ifdef CRLIBM
  integer nchars
  parameter(nchars=160) ! Same as in daten
  character(len=nchars) ch

  character filefields_fields(str_maxFields)*(mStrLen)
  integer filefields_nfields
  integer filefields_lfields(str_maxFields)
  logical filefields_lerr

  real(kind=fPrec) round_near
  integer errno
#endif

#ifdef BOINC
  character(len=256) filename
#endif

  if(dynk_nFuncs+1 > dynk_maxFuncs) then
    dynk_maxFuncs = dynk_maxFuncs + 10
    call alloc(dynk_funcs,dynk_maxFuncs,5,0,"dynk_funcs")
  end if

    if (lFields(2).gt.mStrLen-1 .or. lFields(2) .gt. 20) then
        write(lout,*) "ERROR in DYNK block parsing (fort.3):"
        write(lout,*) "Max length of a FUN name is the smallest of", mStrLen-1, "and", 20, "."
        write(lout,*) "The limitation of 20 comes from the output to dynksets.dat."
        write(lout,*) "Offending FUN: '"//gFields(2)(1:lFields(2))//"'"
        write(lout,*) "length:", lFields(2)
        call prror(51)
    endif

    ! Parse the different functions
    !===============================

    !!! System functions: #0-19 !!!
    select case ( gFields(3)(1:lFields(3)) )

    case ("GET")
        ! GET: Store the value of an element/value

        call dynk_checkargs(nFields,5,"FUN funname GET elementName attribute")
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
        if (lFields(4) .gt. mNameLen .or. lFields(4) .gt. mStrLen-1 ) then
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "FUN GET got an element name with     "
            write (lout,*) "length =", lFields(4), "> ", mNameLen
            write (lout,*) "or > ",mStrLen-1
            write (lout,*) "The name was: '",gFields(4)(1:lFields(4)),"'"
            write (lout,*) "*************************************"
            call prror(51)
        end if
        if (lFields(5) .gt. mStrLen-1) then
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "FUN GET got an attribute name with   "
            write (lout,*) "length =", lFields(5)
            write (lout,*) "> ",mStrLen-1
            write (lout,*) "The name was: '",gFields(5)(1:lFields(5)),"'"
            write (lout,*) "*************************************"
            call prror(51)
        end if

        ! Store data
        ! NAME
        dynk_cData(dynk_ncData  )(1:lFields(2)) = gFields(2)(1:lFields(2))
        ! ELEMENT_NAME
        dynk_cData(dynk_ncData+1)(1:lFields(4)) = gFields(4)(1:lFields(4))
        ! ATTRIBUTE_NAME
        dynk_cData(dynk_ncData+2)(1:lFields(5)) = gFields(5)(1:lFields(5))

        dynk_ncData = dynk_ncData+2

        dynk_fData(dynk_nfData) = -1.0 ! Initialize a place in the array to store the value

    ! END CASE GET

    case ("FILE")
        ! FILE: Load the contents from a file
        ! File format: two ASCII columns of numbers,
        ! first  column = turn number (all turns should be there, starting from 1)
        ! second column = value (as a double)

        call dynk_checkargs(nFields,4,"FUN funname FILE filename")
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

        ! Sanity checks
        if (lFields(4) .gt. mStrLen-1) then
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "FUN FILE got a filename name with   "
            write (lout,*) "length =", lFields(4)
            write (lout,*) "> ",mStrLen-1
            write (lout,*) "The name was: '",gFields(4)(1:lFields(4)),"'"
            write (lout,*) "*************************************"
            call prror(51)
        end if

        ! Store data
        ! NAME
        dynk_cData(dynk_ncData  )(1:lFields(2)) = gFields(2)(1:lFields(2))
        ! FILE NAME
        dynk_cData(dynk_ncData+1)(1:lFields(4)) = gFields(4)(1:lFields(4))

        dynk_ncData = dynk_ncData+1

        ! Open the file
        inquire(unit=ldynkFUNfileunit,opened=lopen)
        if (lopen) then
            write(lout,*)"DYNK> **** ERROR in dynk_parseFUN():FILE ****"
            write(lout,*)"DYNK> Could not open file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
            call prror(-1)
        end if

#ifdef BOINC
        call boincrf(dynk_cData(dynk_ncData),filename)
        open(unit=ldynkFUNfileunit,file=filename,action='read',iostat=stat,status="OLD")
#else
        open(unit=ldynkFUNfileunit,file=dynk_cData(dynk_ncData),action='read',iostat=stat,status="OLD")
#endif
        if (stat .ne. 0) then
            write(lout,*) "DYNK> dynk_parseFUN():FILE"
            write(lout,*) "DYNK> Error opening file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))// "'"
            call prror(51)
        end if

        ! Count number of lines and allocate space
        nlines = 0
        do
            read(ldynkFUNfileunit,*, iostat=stat)
            if (stat .ne. 0) exit
            nlines = nlines + 1
        end do
        rewind(ldynkFUNfileunit)
        call dynk_checkspace(0,nlines,0)

        ii = 0 ! Number of data lines read
        do
#ifndef CRLIBM
            read(ldynkFUNfileunit,*, iostat=stat) t,y
            if (stat .ne. 0) exit !EOF
#else
            read(ldynkFUNfileunit,'(a)', iostat=stat) ch
            if (stat .ne. 0) exit !EOF
            call getfields_split(ch,filefields_fields,filefields_lfields,filefields_nfields,filefields_lerr)
            if (filefields_lerr) then
                write(lout,*) "DYNK> dynk_parseFUN():FILE"
                write(lout,*) "DYNK> Error reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> Error in getfields_split"
                call prror(-1)
            end if

            if (filefields_nfields .ne. 2) then
                write(lout,*) "DYNK> dynk_parseFUN():FILE"
                write(lout,*) "DYNK> Error reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> expected 2 fields, got", filefields_nfields, "ch =",ch
                call prror(-1)
            end if

            read(filefields_fields(1)(1:filefields_lfields(1)),*) t
            y = round_near(errno, filefields_lfields(2)+1, filefields_fields(2))
            if (errno.ne.0) call rounderr(errno,filefields_fields,2,y)
#endif

            ii = ii+1
            if (t .ne. ii) then
                write(lout,*) "DYNK> dynk_parseFUN():FILE"
                write(lout,*) "DYNK> Error reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> Missing turn number", ii,", got turn", t
                call prror(51)
            end if
            ! call dynk_checkspace(0,1,0)

            dynk_nfData = dynk_nfData+1
            dynk_fData(dynk_nfData) = y
        end do
        dynk_funcs(dynk_nFuncs,5) = ii

        close(ldynkFUNfileunit)

    ! END CASE FILE

    case ("FILELIN")
        ! FILELIN: Load the contents from a file, linearly interpolate
        ! File format: two ASCII columns of numbers,
        ! first  column = turn number (as a double)
        ! second column = value (as a double)

        call dynk_checkargs(nFields,4,"FUN funname FILELIN filename")
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

        !Sanity checks
        if (lFields(4) .gt. mStrLen-1) then
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "FUN FILELIN got a filename name with   "
            write (lout,*) "length =", lFields(4)
            write (lout,*) "> ",mStrLen-1
            write (lout,*) "The name was: '",gFields(4)(1:lFields(4)),"'"
            write (lout,*) "*************************************"
            call prror(51)
        end if

        ! Store data
        ! NAME
        dynk_cData(dynk_ncData  )(1:lFields(2)) = gFields(2)(1:lFields(2))
        ! FILE NAME
        dynk_cData(dynk_ncData+1)(1:lFields(4)) = gFields(4)(1:lFields(4))

        dynk_ncData = dynk_ncData+1

        ! Open the file
        inquire(unit=ldynkFUNfileunit, opened=lopen)
        if (lopen) then
            write(lout,*) "DYNK> **** ERROR in dynk_parseFUN():FILELIN ****"
            write(lout,*) "DYNK> Could not open file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
            call prror(-1)
        end if

#ifdef BOINC
        call boincrf(dynk_cData(dynk_ncData),filename)
        open(unit=ldynkFUNfileunit,file=filename,action='read',iostat=stat,status='OLD')
#else
        open(unit=ldynkFUNfileunit,file=dynk_cData(dynk_ncData),action='read',iostat=stat,status='OLD')
#endif

        if (stat .ne. 0) then
            write(lout,*) "DYNK> dynk_parseFUN():FILELIN"
            write(lout,*) "DYNK> Error opening file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
            call prror(51)
        end if

        ! Find the size of the file
        ii = 0 ! Number of data lines read
        do
#ifndef CRLIBM
            read(ldynkFUNfileunit,*, iostat=stat) x,y
            if (stat .ne. 0) exit !EOF
#else
            read(ldynkFUNfileunit,'(a)', iostat=stat) ch
            if (stat .ne. 0) exit !EOF
            call getfields_split(ch,filefields_fields,filefields_lfields,filefields_nfields,filefields_lerr)
            if (filefields_lerr) then
                write(lout,*) "DYNK> dynk_parseFUN():FILELIN"
                write(lout,*) "DYNK> Error reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> Error in getfields_split"
                call prror(-1)
            end if

            if (filefields_nfields  .ne. 2) then
                write(lout,*) "DYNK> dynk_parseFUN():FILELIN"
                write(lout,*) "DYNK> Error reading file '"// &
                              trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> expected 2 fields, got",filefields_nfields,"ch =",ch
                call prror(-1)
            end if

            x = round_near(errno, filefields_lfields(1)+1, filefields_fields(1))
            if (errno.ne.0) call rounderr(errno,filefields_fields,1,x)
            y = round_near(errno, filefields_lfields(2)+1, filefields_fields(2))
            if (errno.ne.0) call rounderr(errno,filefields_fields,2,y)
#endif

            if (ii .gt. 0 .and. x .le. x2) then ! Insane: Decreasing x
                write (lout,*) "DYNK> dynk_parseFUN():FILELIN"
                write (lout,*) "DYNK> Error while reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write (lout,*) "DYNK> x values must be in increasing order"
                call prror(-1)
            end if
            x2 = x
            ii = ii+1
        end do
        t = ii
        rewind(ldynkFUNfileunit)

        call dynk_checkspace(0,2*t,0)
        ! if (dynk_nfData+2*t .gt. maxdata_dynk) then
        !     write (lout,*) "DYNK> dynk_parseFUN():FILELIN"
        !     write (lout,*) "DYNK> Error reading file '"// &
        !                    trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
        !     write (lout,*) "DYNK> Not enough space in dynk_fData, need", 2*t
        !     write (lout,*) "DYNK> Please increase maxdata_dynk"
        !     call prror(51)
        ! end if

        ! Read the file
        ii = 0
        do
#ifndef CRLIBM
            read(ldynkFUNfileunit,*, iostat=stat) x,y
            if (stat .ne. 0) then ! EOF
                if (ii .ne. t) then
                    write (lout,*)"DYNK> dynk_parseFUN():FILELIN"
                    write (lout,*)"DYNK> Unexpected when reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                    write (lout,*)"DYNK> ii=",ii,"t=",t
                    call prror(51)
                end if
                exit
            end if
#else
            read(ldynkFUNfileunit,'(a)', iostat=stat) ch
            if (stat .ne. 0) then !EOF
                if (ii .ne. t) then
                    write (lout,*)"DYNK> dynk_parseFUN():FILELIN"
                    write (lout,*)"DYNK> Unexpected when reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                    write (lout,*)"DYNK> ii=",ii,"t=",t
                    call prror(51)
                end if
                exit
            end if

            call getfields_split(ch,filefields_fields,filefields_lfields,filefields_nfields,filefields_lerr)
            if (filefields_lerr) then
                write(lout,*) "DYNK> dynk_parseFUN():FILELIN"
                write(lout,*) "DYNK> Error reading file '"// &
                              trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> Error in getfields_split"
                call prror(-1)
            end if

            if (filefields_nfields .ne. 2) then
                write(lout,*) "DYNK> dynk_parseFUN():FILELIN"
                write(lout,*) "DYNK> Error reading file '"// &
                              trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> expected 2 fields, got",filefields_nfields,"ch =",ch
                call prror(-1)
            end if

            x = round_near(errno, filefields_lfields(1)+1, filefields_fields(1))
            if (errno.ne.0) call rounderr(errno,filefields_fields,1,x)
            y = round_near(errno, filefields_lfields(2)+1, filefields_fields(2))
            if (errno.ne.0) call rounderr(errno,filefields_fields,2,y)
            ! write(*,*) "DBGDBG: ch=",ch
            ! write(*,*) "DBGDBG: filefields_fields(1)=", filefields_fields(1)
            ! write(*,*) "DBGDBG: filefields_fields(2)=", filefields_fields(2)
#endif
            ! write(*,*) "DBGDBG: x,y = ",x,y

            ! Current line number
            ii = ii+1

            dynk_fData(dynk_nfData + ii    ) = x
            dynk_fData(dynk_nfData + ii + t) = y
        end do

        dynk_nfData = dynk_nfData + 2*t
        dynk_funcs(dynk_nFuncs,5) = t
        close(ldynkFUNfileunit)

    ! END CASE FILELIN

    case ("PIPE")
        ! PIPE: Use a pair of UNIX FIFOs.
        ! Another program is expected to hook onto the other end of the pipe,
        ! and will recieve a message when SixTrack's dynk_computeFUN() is called.
        ! That program should then send a value back (in ASCII), which will be the new setting.

        call dynk_checkargs(nFields,7,"FUN funname PIPE inPipeName outPipeName ID fileUnit" )
        call dynk_checkspace(1,0,4)

#ifdef CR
        write(lout,*) "DYNK FUN PIPE not supported in CR version"
        write(lout,*) "Sorry :("
        call prror(-1)
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
        if (lFields(4) .gt. mStrLen-1 .or. &
            lFields(5) .gt. mStrLen-1 .or. &
            lFields(6) .gt. mStrLen-1 ) then

            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "FUN PIPE got one or more strings which "
            write (lout,*) "was too long (>",mStrLen-1,")"
            write (lout,*) "Strings: '",                                          &
                           gFields(4)(1:lFields(4)),"' and '", &
                           gFields(5)(1:lFields(5)),"' and '", &
                           gFields(6)(1:lFields(6)),"'."
            write (lout,*) "lengths =",                  &
                           lFields(4),", ",    &
                           lFields(5)," and ", &
                           lFields(6)
            write (lout,*) "*************************************"
            call prror(51)
        end if

        ! Store data
        ! NAME
        dynk_cData(dynk_ncData  )(1:lFields(2)) = gFields(2)(1:lFields(2))
        ! inPipe
        dynk_cData(dynk_ncData+1)(1:lFields(4)) = gFields(4)(1:lFields(4))
        ! outPipe
        dynk_cData(dynk_ncData+2)(1:lFields(5)) = gFields(5)(1:lFields(5))
        ! ID
        dynk_cData(dynk_ncData+3)(1:lFields(6)) = gFields(6)(1:lFields(6))

        dynk_ncData = dynk_ncData+3

        ! fileUnit
        read(gFields(7)(1:lFields(7)),*) dynk_iData(dynk_niData)

        ! Look if the fileUnit or filenames are used in a different FUN PIPE
        t=0 !Used to hold the index of the other pipe; t=0 if no older pipe -> open files.
        do ii=1,dynk_nFuncs-1
            if (dynk_funcs(ii,2) .eq. 3) then !It's a PIPE
                !Does any of the settings match?
                if (dynk_iData(dynk_funcs(ii,3))  .eq.dynk_iData(dynk_niData)   .or. & ! Unit number
                    dynk_cData(dynk_funcs(ii,1)+1).eq.dynk_cData(dynk_ncData-2) .or. & ! InPipe filename
                    dynk_cData(dynk_funcs(ii,1)+2).eq.dynk_cData(dynk_ncData-1)) then  ! OutPipe filename
                    ! Does *all* of the settings match?
                    if (dynk_iData(dynk_funcs(ii,3))  .eq.dynk_iData(dynk_niData)   .and. & ! Unit number
                        dynk_cData(dynk_funcs(ii,1)+1).eq.dynk_cData(dynk_ncData-2) .and. & ! InPipe filename
                        dynk_cData(dynk_funcs(ii,1)+2).eq.dynk_cData(dynk_ncData-1)) then   ! OutPipe filename

                        t=ii
                        write(lout,*) "DYNK> PIPE FUN '" // &
                                      trim(stringzerotrim(dynk_cData(dynk_funcs(dynk_nFuncs,1))))// &
                                      "' using same settings as previously defined FUN '"// &
                                      trim(stringzerotrim(dynk_cData(dynk_funcs(ii,1))))// &
                                      "' -> reusing files !"
                        if (dynk_cData(dynk_funcs(ii,1)+3).eq.dynk_cData(dynk_ncData)) then ! ID
                            write(lout,*) "DYNK> ERROR: IDs must be different when sharing PIPEs."
                            call prror(-1)
                        end if
                        exit ! Break loop

                    else ! Partial match
                        write(lout,*) "DYNK> *** Error in dynk_parseFUN():PIPE ***"
                        write(lout,*) "DYNK> Partial match of inPipe/outPipe/unit number"
                        write(lout,*) "DYNK> between PIPE FUN '"// &
                                      trim(stringzerotrim(dynk_cData(dynk_funcs(dynk_nFuncs,1))))// &
                                      "' and '"// &
                                      trim(stringzerotrim(dynk_cData(dynk_funcs(ii,1))))//"'"
                        call prror(-1)
                    end if
                end if
            end if
        end do

        if (t .eq. 0) then ! Must open a new set of files
            ! Open the inPipe
            inquire(unit=dynk_iData(dynk_niData), opened=lopen)
            if (lopen) then
                write(lout,*)"DYNK> **** ERROR in dynk_parseFUN():PIPE ****"
                write(lout,*)"DYNK> unit",dynk_iData(dynk_niData),"for file '"// &
                             trim(stringzerotrim(dynk_cData(dynk_ncData-2)))// &
                             "' was already taken"
                call prror(-1)
            end if

            write(lout,*) "DYNK> Opening input pipe '"// &
                          trim(stringzerotrim(dynk_cData(dynk_ncData-2)))//"' for FUN '"// &
                          trim(stringzerotrim(dynk_cData(dynk_ncData-3)))//"', ID='"// &
                          trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"

            ! DYNK PIPE does not support the CR version, so BOINC support (call boincrf()) isn't needed
            open(unit=dynk_iData(dynk_niData),file=dynk_cData(dynk_ncData-2),action='read',iostat=stat,status="OLD")
            if (stat .ne. 0) then
                write(lout,*) "DYNK> dynk_parseFUN():PIPE"
                write(lout,*) "DYNK> Error opening file '"// &
                              trim(stringzerotrim(dynk_cData(dynk_ncData-2)))//"' stat=",stat
                call prror(51)
            end if

            ! Open the outPipe
            write(lout,*) "DYNK> Opening output pipe '"// &
                          trim(stringzerotrim(dynk_cData(dynk_ncData-1)))//"' for FUN '"// &
                          trim(stringzerotrim(dynk_cData(dynk_ncData-3)))//"', ID='"// &
                          trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
            inquire(unit=dynk_iData(dynk_niData)+1, opened=lopen)
            if (lopen) then
                write(lout,*)"DYNK> **** ERROR in dynk_parseFUN():PIPE ****"
                write(lout,*)"DYNK> unit",dynk_iData(dynk_niData)+1,"for file '"// &
                             trim(stringzerotrim(dynk_cData(dynk_ncData-1)))//"' was already taken"
                call prror(-1)
            end if

            ! DYNK PIPE does not support the CR version, so BOINC support (call boincrf()) isn't needed
            open(unit=dynk_iData(dynk_niData)+1,file=dynk_cData(dynk_ncData-1),action='write',iostat=stat,status="OLD")
            if (stat .ne. 0) then
                write(lout,*) "DYNK> dynk_parseFUN():PIPE"
                write(lout,*) "DYNK> Error opening file '"// &
                              trim(stringzerotrim(dynk_cData(dynk_ncData-1)))//"' stat=",stat
                call prror(51)
            end if
            write(dynk_iData(dynk_niData)+1,'(a)') "DYNKPIPE !******************!" ! Once per file

        end if ! End "if (t.eq.0)"/must open new files

        write(dynk_iData(dynk_niData)+1,'(a)') "INIT ID="// & ! Once per ID
                                               trim(stringzerotrim(dynk_cData(dynk_ncData)))// &
                                               " for FUN="// &
                                               trim(stringzerotrim(dynk_cData(dynk_ncData-3)))

    ! END CASE PIPE

    case ("RANDG")
        ! RANDG: Gausian random number with mu, sigma, and optional cutoff

        call dynk_checkargs(nFields,8,"FUN funname RANDG seed1 seed2 mu sigma cut")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

        read(gFields(4)(1:lFields(4)),*) dynk_iData(dynk_niData)   ! seed1 (initial)
        read(gFields(5)(1:lFields(5)),*) dynk_iData(dynk_niData+1) ! seed2 (initial)
#ifndef CRLIBM
        read(gFields(6)(1:lFields(6)),*) dynk_fData(dynk_nfData)   ! mu
        read(gFields(7)(1:lFields(7)),*) dynk_fData(dynk_nfData+1) ! sigma
#else
        dynk_fData(dynk_nfData) = round_near(errno,lFields(6)+1,gFields(6))   ! mu
        if (errno.ne.0) call rounderr(errno,gFields,6,dynk_fData(dynk_nfData))
        dynk_fData(dynk_nfData+1) = round_near(errno,lFields(7)+1,gFields(7)) ! sigma
        if (errno.ne.0) call rounderr( errno,gFields,7,dynk_fData(dynk_nfData+1))
#endif
        read(gFields(8)(1:lFields(8)),*) dynk_iData(dynk_niData+2) ! mcut

        dynk_iData(dynk_niData+3) = 0 ! seed1 (current)
        dynk_iData(dynk_niData+4) = 0 ! seed2 (current)

        dynk_niData = dynk_niData+4
        dynk_nfData = dynk_nfData+1

        if (dynk_iData(dynk_funcs(dynk_nFuncs,3)+2) .lt. 0) then
            ! mcut < 0
            write (lout,*) "DYNK> dynk_parseFUN():RANDG"
            write (lout,*) "DYNK> ERROR in DYNK block parsing (fort.3)"
            write (lout,*) "DYNK> mcut must be >= 0"
            call prror(51)
        end if

    ! END CASE RANDG

    case ("RANDU")
        ! RANDU: Uniform random number

        call dynk_checkargs(nFields,5,"FUN funname RANDU seed1 seed2")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

        read(gFields(4)(1:lFields(4)),*) dynk_iData(dynk_niData)   ! seed1 (initial)
        read(gFields(5)(1:lFields(5)),*) dynk_iData(dynk_niData+1) ! seed2 (initial)

        dynk_iData(dynk_niData+2) = 0 ! seed1 (current)
        dynk_iData(dynk_niData+3) = 0 ! seed2 (current)

        dynk_niData = dynk_niData+3

    ! END CASE RANDU

    case("RANDON")
        ! RANDON: Turn by turn ON for one turn with the probability P, else OFF
        call dynk_checkargs(nFields,6,"FUN funname RANDON seed1 seed2 P")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

        read(gFields(4)(1:lFields(4)),*) dynk_iData(dynk_niData)   ! seed1 (initial)
        read(gFields(5)(1:lFields(5)),*) dynk_iData(dynk_niData+1) ! seed2 (initial)
        read(gFields(6)(1:lFields(6)),*) dynk_fData(dynk_nfData)   ! P

        dynk_iData(dynk_niData+2) = 0 ! seed1 (current)
        dynk_iData(dynk_niData+3) = 0 ! seed2 (current)

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

        call dynk_checkargs(nFields,6,"FUN funname {FIR|IIR} N filename baseFUN")

        select case(gFields(3)(1:lFields(3)))
        case("FIR")
            isFIR = .true.
        case("IIR")
            isFIR = .false.
        case default
            write (lout,*) "DYNK> dynk_parseFUN():FIR/IIR"
            write (lout,*) "DYNK> non-recognized type in inner switch?"
            write (lout,*) "DYNK> Got: '"//gFields(3)(1:lFields(3))//"'"
            call prror(-1)
        end select

        read(gFields(4)(1:lFields(4)),*) t ! N
        if (isFIR) then
            call dynk_checkspace(0,3*(t+1),2)
        else
            call dynk_checkspace(0,6*(t+1),2)
        end if

        ! Set pointers to start of funs data blocks
        dynk_nFuncs = dynk_nFuncs+1
        dynk_ncData = dynk_ncData+1

        ! Store pointers
        dynk_funcs(dynk_nFuncs,1) = dynk_ncData   ! NAME (in dynk_cData)
        if (isFIR) then
            dynk_funcs(dynk_nFuncs,2) = 10 ! TYPE (FIR)
        else
            dynk_funcs(dynk_nFuncs,2) = 11 ! TYPE (IIR)
        end if
        dynk_funcs(dynk_nFuncs,3) = dynk_nfData+1 ! ARG1 (start of float storage)
        dynk_funcs(dynk_nFuncs,4) = t             ! ARG2 (filter order N)
        dynk_funcs(dynk_nFuncs,5) = &             ! ARG3 (filtered function)
              dynk_findFUNindex( gFields(6) (1:lFields(6)), 1)

        ! Store metadata
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))
        read(gFields(4)(1:lFields(4)),*) dynk_iData(dynk_niData) ! N

        ! Sanity check
        if (dynk_funcs(dynk_nFuncs,5).eq.-1) then
            call dynk_dumpdata
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "FIR/IIR function wanting function '", &
                           gFields(6)(1:lFields(6)), "'"
            write (lout,*) "This FUN is unknown!"
            write (lout,*) "*************************************"
            call prror(51)
        end if
        if (lFields(5) .gt. mStrLen-1) then
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "FUN FIR/IIR got a filename name with "
            write (lout,*) "length =", lFields(5)
            write (lout,*) "> ",mStrLen-1
            write (lout,*) "The name was: '",gFields(5)(1:lFields(5)),"'"
            write (lout,*) "*************************************"
            call prror(51)
        end if
        if (dynk_iData(dynk_niData) .le. 0) then
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "FUN FIR/IIR got N <= 0, this is not valid"
            write (lout,*) "*************************************"
            call prror(51)
        end if

        ! More metadata
        dynk_ncData = dynk_ncData+1
        ! FILE NAME
        dynk_cData(dynk_ncData)(1:lFields(5)) = gFields(5)(1:lFields(5))

        ! Read the file
        inquire(unit=ldynkFUNfileunit, opened=lopen)
        if (lopen) then
            write(lout,*) "DYNK> **** ERROR in dynk_parseFUN():FIR/IIR ****"
            write(lout,*) "DYNK> Could not open file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
            call prror(-1)
        end if
#ifdef BOINC
        call boincrf(dynk_cData(dynk_ncData),filename)
        open(unit=ldynkFUNfileunit,file=filename,action='read',iostat=stat,status="OLD")
#else
        open(unit=ldynkFUNfileunit,file=dynk_cData(dynk_ncData),action='read',iostat=stat,status="OLD")
#endif
        if (stat .ne. 0) then
            write(lout,*) "DYNK> dynk_parseFUN():FIR/IIR"
            write(lout,*) "DYNK> Error opening file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
            call prror(51)
        end if

        do ii=0, dynk_funcs(dynk_nFuncs,4)

            ! Reading the FIR/IIR file without CRLIBM
#ifndef CRLIBM
            if (isFIR) then
                read(ldynkFUNfileunit,*,iostat=stat) t, x, y
            else
                read(ldynkFUNfileunit,*,iostat=stat) t, x, y, z, u
            end if
            if (stat.ne.0) then
                write(lout,*) "DYNK> dynk_parseFUN():FIR/IIR"
                write(lout,*) "DYNK> Error reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> File ended unexpectedly at ii =",ii
                call prror(-1)
            end if
#else
            ! Reading the FIR/IIR file with CRLIBM
            read(ldynkFUNfileunit,'(a)', iostat=stat) ch
            if (stat.ne.0) then
                write(lout,*) "DYNK> dynk_parseFUN():FIR/IIR"
                write(lout,*) "DYNK> Error reading file '"//trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> File ended unexpectedly at ii =",ii
                call prror(-1)
            end if

            call getfields_split(ch,filefields_fields,filefields_lfields,filefields_nfields,filefields_lerr)

            ! Sanity checks
            if (filefields_lerr) then
                write(lout,*) "DYNK> dynk_parseFUN():FIR/IIR"
                write(lout,*) "DYNK> Error reading file '", &
                              trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> Error in getfields_split()"
                call prror(-1)
            end if
            if ( (      isFIR .and.filefields_nfields .ne. 3) .or. &
                 ((.not.isFIR).and.filefields_nfields .ne. 5)     ) then
                write(lout,*) "DYNK> dynk_parseFUN():FIR/IIR"
                write(lout,*) "DYNK> Error reading file '"// &
                              trim(stringzerotrim(dynk_cData(dynk_ncData)))//"', line =", ii
                write(lout,*) "DYNK> Expected 3[5] fields (idx, fac, init, selfFac, selfInit), ", &
                              "got ",filefields_nfields
                call prror(-1)
            endif

            ! Read the data into t,x,y(,z,u):
            read(filefields_fields(1)(1:filefields_lfields(1)),*) t

            x = round_near(errno,filefields_lfields(2)+1,filefields_fields(2))
            if (errno.ne.0) call rounderr(errno,filefields_fields,2,x)

            y = round_near(errno,filefields_lfields(3)+1,filefields_fields(3))
            if (errno.ne.0) call rounderr(errno,filefields_fields,3,y)

            if (.not.isFIR) then
                z = round_near(errno,filefields_lfields(4)+1,filefields_fields(4))
                if (errno.ne.0) call rounderr(errno,filefields_fields,4,z)
                u = round_near(errno,filefields_lfields(5)+1,filefields_fields(5))
                if (errno.ne.0) call rounderr(errno,filefields_fields,5,u)
            end if
#endif
            ! More sanity checks
            if (t .ne. ii) then
                write(lout,*) "DYNK> dynk_parseFUN():FIR/IIR"
                write(lout,*) "DYNK> Error reading file '"// &
                              trim(stringzerotrim(dynk_cData(dynk_ncData)))//"'"
                write(lout,*) "DYNK> Got line t =",t, ", expected ", ii
                call prror(-1)
            end if
            ! Save data to arrays
            ! Store coefficients (x) and initial/earlier values (y) in interlaced order
            dynk_nfData = dynk_nfData+1
            dynk_fData(dynk_nfData) = x      ! b_i
            dynk_nfData = dynk_nfData+1
            dynk_fData(dynk_nfData) = 0.0    ! x[n-1], will be initialized in dynk_apply()
            dynk_nfData = dynk_nfData+1
            dynk_fData(dynk_nfData) = y      ! x_init[n-i]
            if (.not.isFIR) then
                dynk_nfData = dynk_nfData+1
                dynk_fData(dynk_nfData) = z   ! a_i
                dynk_nfData = dynk_nfData+1
                dynk_fData(dynk_nfData) = 0.0 ! y[n-i], will be initialized in dynk_apply()
                dynk_nfData = dynk_nfData+1
                dynk_fData(dynk_nfData) = u   ! y_init[n-i]
            end if
        end do
        close(ldynkFUNfileunit)

    ! END CASES FIR & IIR

    case("ADD","SUB","MUL","DIV","POW") ! Operators: #20-39
        ! Two-argument operators  y = OP(f1, f2)

        call dynk_checkargs(nFields,5,"FUN funname {ADD|SUB|MUL|DIV|POW} funname1 funname2")
        call dynk_checkspace(0,0,1)

        ! Set pointers to start of funs data blocks
        dynk_nFuncs = dynk_nFuncs+1
        dynk_ncData = dynk_ncData+1

        ! Store pointers
        ! NAME (in dynk_cData)
        dynk_funcs(dynk_nFuncs,1) = dynk_ncData
        select case (gFields(3)(1:lFields(3)))
        case ("ADD")
            dynk_funcs(dynk_nFuncs,2) = 20 ! TYPE (ADD)
        case ("SUB")
            dynk_funcs(dynk_nFuncs,2) = 21 ! TYPE (SUB)
        case ("MUL")
            dynk_funcs(dynk_nFuncs,2) = 22 ! TYPE (MUL)
        case ("DIV")
            dynk_funcs(dynk_nFuncs,2) = 23 ! TYPE (DIV)
        case ("POW")
            dynk_funcs(dynk_nFuncs,2) = 24 ! TYPE (POW)
        case default
            write (lout,*) "DYNK> dynk_parseFUN() : 2-arg function"
            write (lout,*) "DYNK> non-recognized type in inner switch"
            write (lout,*) "DYNK> Got: '"//gFields(3)(1:lFields(3))//"'"
            call prror(51)
        end select
        dynk_funcs(dynk_nFuncs,3) =  &
              dynk_findFUNindex( gFields(4)(1:lFields(4)), 1) ! Index to f1
        dynk_funcs(dynk_nFuncs,4) =  &
              dynk_findFUNindex( gFields(5)(1:lFields(5)), 1) ! Index to f2
        dynk_funcs(dynk_nFuncs,5) = -1  ! ARG3

        ! Store data
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))
        ! Sanity check (string lengths are done inside dynk_findFUNindex)
        if (dynk_funcs(dynk_nFuncs,3) .eq. -1 .or. dynk_funcs(dynk_nFuncs,4) .eq. -1) then
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "TWO ARG OPERATOR wanting functions '", &
                           gFields(4)(1:lFields(4)), "' and '",  &
                           gFields(5)(1:lFields(5)), "'"
            write (lout,*) "Calculated indices:", &
                           dynk_funcs(dynk_nFuncs,3), dynk_funcs(dynk_nFuncs,4)
            write (lout,*) "One or both of these are not known (-1)."
            write (lout,*) "*************************************"
            call dynk_dumpdata
            call prror(51)
        end if

    ! END CASES ADD, SUB, MUL, DIV & POW

    case ("MINUS","SQRT","SIN","COS","LOG","LOG10","EXP")
        ! One-argument operators  y = OP(f1)

        call dynk_checkargs(nFields,4,"FUN funname {MINUS|SQRT|SIN|COS|LOG|LOG10|EXP} funname")
        call dynk_checkspace(0,0,1)

        ! Set pointers to start of funs data blocks
        dynk_nFuncs = dynk_nFuncs+1
        dynk_ncData = dynk_ncData+1

        ! Store pointers
        dynk_funcs(dynk_nFuncs,1) = dynk_ncData ! NAME (in dynk_cData)
        select case (gFields(3)(1:lFields(3)))
        case ("MINUS")
            dynk_funcs(dynk_nFuncs,2) = 30 ! TYPE (MINUS)
        case ("SQRT")
            dynk_funcs(dynk_nFuncs,2) = 31 ! TYPE (SQRT)
        case ("SIN")
            dynk_funcs(dynk_nFuncs,2) = 32 ! TYPE (SIN)
        case ("COS")
            dynk_funcs(dynk_nFuncs,2) = 33 ! TYPE (COS)
        case ("LOG")
            dynk_funcs(dynk_nFuncs,2) = 34 ! TYPE (LOG)
        case ("LOG10")
            dynk_funcs(dynk_nFuncs,2) = 35 ! TYPE (LOG10)
        case ("EXP")
            dynk_funcs(dynk_nFuncs,2) = 36 ! TYPE (EXP)
        case default
            write (lout,*) "DYNK> dynk_parseFUN() : 1-arg function"
            write (lout,*) "DYNK> non-recognized type in inner switch?"
            write (lout,*) "DYNK> Got: '"//gFields(3)(1:lFields(3))//"'"
            call prror(51)
        end select

        ! Index to f1
        dynk_funcs(dynk_nFuncs,3) = dynk_findFUNindex(gFields(4)(1:lFields(4)),1)
        dynk_funcs(dynk_nFuncs,5) = -1 ! ARG3

        ! Store data
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))
        ! Sanity check (string lengths are done inside dynk_findFUNindex)
        if (dynk_funcs(dynk_nFuncs,3) .eq. -1) then
            write (lout,*) "*************************************"
            write (lout,*) "ERROR in DYNK block parsing (fort.3):"
            write (lout,*) "SINGLE OPERATOR FUNC wanting function '", &
                           gFields(4)(1:lFields(4)), "'"
            write (lout,*) "Calculated index:",dynk_funcs(dynk_nFuncs,3)
            write (lout,*) "One or both of these are not known (-1)."
            write (lout,*) "*************************************"
            call dynk_dumpdata
            call prror(51)
        end if

    ! END CASES MINUS, SQRT, SIN, COS, LOG, LOG10 & EXP

    ! Polynomial & Elliptical functions: # 40-59

    case("CONST")
        ! CONST: Just a constant value

        call dynk_checkargs(nFields,4,"FUN funname CONST value")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

#ifndef CRLIBM
        read(gFields(4)(1:lFields(4)),*) dynk_fData(dynk_nfData) ! value
#else
        dynk_fData(dynk_nfData) = round_near(errno,lFields(4)+1,gFields(4)) ! value
        if (errno.ne.0) call rounderr( errno,gFields,4,dynk_fData(dynk_nfData))
#endif

    ! END CASE CONST

    case ("TURN")
        ! TURN: Just the current turn number

        call dynk_checkargs(nFields,3,"FUN funname TURN")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

    ! END CASE TURN

    case ("LIN")
        ! LIN: Linear ramp y = dy/dt*T+b

        call dynk_checkargs(nFields,5,"FUN funname LIN dy/dt b")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

#ifndef CRLIBM
        read(gFields(4)(1:lFields(4)),*) dynk_fData(dynk_nfData)   ! dy/dt
        read(gFields(5)(1:lFields(5)),*) dynk_fData(dynk_nfData+1) ! b
#else
        dynk_fData(dynk_nfData)   = round_near(errno,lFields(4)+1,gFields(4))  ! dy/dt
        if (errno.ne.0) call rounderr( errno,gFields,4,dynk_fData(dynk_nfData))
        dynk_fData(dynk_nfData+1) = round_near(errno,lFields(5)+1, gFields(5)) ! b
        if (errno.ne.0) call rounderr( errno,gFields,5,dynk_fData(dynk_nfData+1))
#endif
        dynk_nfData = dynk_nfData + 1

    ! END CASE LIN

    case ("LINSEG")
        ! LINSEG: Linear ramp between points (x1,y1) and (x2,y2)

        call dynk_checkargs(nFields,7,"FUN funname LINSEG x1 x2 y1 y2")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))
#ifndef CRLIBM
        read(gFields(4)(1:lFields(4)),*) dynk_fData(dynk_nfData)   ! x1
        read(gFields(5)(1:lFields(5)),*) dynk_fData(dynk_nfData+1) ! x2
        read(gFields(6)(1:lFields(6)),*) dynk_fData(dynk_nfData+2) ! y1
        read(gFields(7)(1:lFields(7)),*) dynk_fData(dynk_nfData+3) ! y2
#else
        dynk_fData(dynk_nfData)   = round_near(errno,lFields(4)+1,gFields(4)) ! x1
        if (errno.ne.0) call rounderr( errno,gFields,4,dynk_fData(dynk_nfData))
        dynk_fData(dynk_nfData+1) = round_near(errno,lFields(5)+1,gFields(5)) ! x2
        if (errno.ne.0) call rounderr( errno,gFields,5,dynk_fData(dynk_nfData+1))
        dynk_fData(dynk_nfData+2) = round_near(errno,lFields(6)+1,gFields(6)) ! y1
        if (errno.ne.0) call rounderr( errno,gFields,6,dynk_fData(dynk_nfData+2))
        dynk_fData(dynk_nfData+3) = round_near(errno,lFields(7)+1,gFields(7)) ! y2
        if (errno.ne.0) call rounderr( errno,gFields,7,dynk_fData(dynk_nfData+3))
#endif
        dynk_nfData = dynk_nfData + 3

        if (dynk_fData(dynk_nfData-3).eq.dynk_fData(dynk_nfData-2)) then
            write (lout,*) "ERROR in DYNK block parsing (fort.3)"
            write (lout,*) "LINSEG: x1 and x2 must be different."
            call prror(51)
        end if

    ! END CASE LINSEG

    case ("QUAD")
        ! QUAD: Quadratic ramp y = a*T^2 + b*T + c

        call dynk_checkargs(nFields,6,"FUN funname QUAD a b c")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

#ifndef CRLIBM
        read(gFields(4)(1:lFields(4)),*) dynk_fData(dynk_nfData)   ! a
        read(gFields(5)(1:lFields(5)),*) dynk_fData(dynk_nfData+1) ! b
        read(gFields(6)(1:lFields(6)),*) dynk_fData(dynk_nfData+2) ! c
#else
        dynk_fData(dynk_nfData)   = round_near(errno,lFields(4)+1,gFields(4)) ! a
        if (errno.ne.0) call rounderr(errno,gFields,4,dynk_fData(dynk_nfData))
        dynk_fData(dynk_nfData+1) = round_near(errno,lFields(5)+1,gFields(5)) ! b
        if (errno.ne.0) call rounderr(errno,gFields,5,dynk_fData(dynk_nfData+1))
        dynk_fData(dynk_nfData+2) = round_near(errno,lFields(6)+1,gFields(6)) ! c
        if (errno.ne.0) call rounderr( errno,gFields,6,dynk_fData(dynk_nfData+2))
#endif
        dynk_nfData = dynk_nfData + 2

    ! END CASE QUAD

    case ("QUADSEG")
        ! QUADSEG: Quadratic ramp y = a*T^2 + b*T + c,
        ! input as start point (x1,y1), end point (x2,y2), derivative at at x1

        call dynk_checkargs(nFields,8,"FUN funname QUADSEG x1 x2 y1 y2 deriv")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))
#ifndef CRLIBM
        read(gFields(4)(1:lFields(4)),*) x1
        read(gFields(5)(1:lFields(5)),*) x2
        read(gFields(6)(1:lFields(6)),*) y1
        read(gFields(7)(1:lFields(7)),*) y2
        read(gFields(8)(1:lFields(8)),*) deriv
#else
        x1 = round_near(errno,lFields(4)+1,gFields(4))    ! x1
        if (errno.ne.0) call rounderr(errno,gFields,4,x1)
        x2 = round_near(errno,lFields(5)+1,gFields(5))    ! x2
        if (errno.ne.0) call rounderr(errno,gFields,5,x2)
        y1 = round_near(errno,lFields(6)+1,gFields(6))    ! y1
        if (errno.ne.0) call rounderr(errno,gFields,6,y1)
        y2 = round_near(errno,lFields(7)+1,gFields(7))    ! y2
        if (errno.ne.0) call rounderr(errno,gFields,7,y2)
        deriv = round_near(errno,lFields(8)+1,gFields(8)) ! deriv
        if (errno.ne.0) call rounderr(errno,gFields,8,deriv)
#endif
        if (x1 .eq. x2) then
            write (lout,*) "ERROR in DYNK block parsing (fort.3)"
            write (lout,*) "QUADSEG: x1 and x2 must be different."
            call prror(51)
        end if

        ! Compute a:
        dynk_fData(dynk_nfData) = deriv/(x1-x2) + (y2-y1)/((x1-x2)**2)
        ! Compute b:
        dynk_fData(dynk_nfData+1) = (y2-y1)/(x2-x1) - (x1+x2)*dynk_fData(dynk_nfData)
        ! Compute c:
        dynk_fData(dynk_nfData+2) = y1 + (            &
                - x1**2 * dynk_fData(dynk_nfData)     &
                - x1    * dynk_fData(dynk_nfData+1) )

        ! Store input data:
        dynk_fData(dynk_nfData+3) = x1
        dynk_fData(dynk_nfData+4) = x2
        dynk_fData(dynk_nfData+5) = y1
        dynk_fData(dynk_nfData+6) = y2
        dynk_fData(dynk_nfData+7) = deriv

        dynk_nfData = dynk_nfData + 7

    ! END CASE QUADSEG

    case ("SINF","COSF","COSF_RIPP")
        ! Trancedental functions: #60-79
        ! SINF     : Sin functions y = A*sin(omega*T+phi)
        ! COSF     : Cos functions y = A*cos(omega*T+phi)
        ! COSF_RIPP: Cos functions y = A*cos(2*pi*(T-1)/period+phi)

        call dynk_checkargs(nFields,6,"FUN funname {SINF|COSF|COSF_RIPP} amplitude {omega|period} phase")
        call dynk_checkspace(0,3,1)

        ! Set pointers to start of funs data blocks
        dynk_nFuncs = dynk_nFuncs+1
        dynk_nfData = dynk_nfData+1
        dynk_ncData = dynk_ncData+1

        ! Store pointers
        ! NAME
        dynk_funcs(dynk_nFuncs,1) = dynk_ncData
        select case (gFields(3)(1:lFields(3)))
        case("SINF")
            dynk_funcs(dynk_nFuncs,2) = 60 ! TYPE (SINF)
        case("COSF")
            dynk_funcs(dynk_nFuncs,2) = 61 ! TYPE (COSF)
        case ("COSF_RIPP")
            dynk_funcs(dynk_nFuncs,2) = 62 ! TYPE (COSF_RIPP)
        case default
            write (lout,*) "DYNK> dynk_parseFUN() : SINF/COSF"
            write (lout,*) "DYNK> non-recognized type in inner switch"
            write (lout,*) "DYNK> Got: '"//gFields(3)(1:lFields(3))//"'"
            call prror(51)
        end select
        dynk_funcs(dynk_nFuncs,3) = dynk_nfData ! ARG1
        dynk_funcs(dynk_nFuncs,4) = -1          ! ARG2
        dynk_funcs(dynk_nFuncs,5) = -1          ! ARG3

        ! Store data
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

#ifndef CRLIBM
        read(gFields(4)(1:lFields(4)),*) dynk_fData(dynk_nfData)   ! A
        read(gFields(5)(1:lFields(5)),*) dynk_fData(dynk_nfData+1) ! omega
        read(gFields(6)(1:lFields(6)),*) dynk_fData(dynk_nfData+2) ! phi
#else
        dynk_fData(dynk_nfData)   = round_near(errno,lFields(4)+1,gFields(4))  ! A
        if (errno.ne.0) call rounderr(errno,gFields,4,dynk_fData(dynk_nfData))
        dynk_fData(dynk_nfData+1) = round_near(errno,lFields(5)+1,gFields(5)) ! omega
        if (errno.ne.0) call rounderr(errno,gFields,5,dynk_fData(dynk_nfData+1))
        dynk_fData(dynk_nfData+2) = round_near(errno,lFields(6)+1,gFields(6)) ! phi
        if (errno.ne.0) call rounderr(errno,gFields,6,dynk_fData(dynk_nfData+2))
#endif
        dynk_nfData = dynk_nfData + 2

    ! END CASE SINF, COSF & COSF_RIPP

    case ("PELP")
        ! PELP: Parabolic/exponential/linear/parabolic
        ! From "Field Computation for Accelerator Magnets:
        ! Analytical and Numerical Methods for Electromagnetic Design and Optimization"
        ! By Dr.-Ing. Stephan Russenschuck
        ! Appendix C: "Ramping the LHC Dipoles"

        call dynk_checkargs(nFields,10,"FUN funname PELP tinj Iinj Inom A D R te")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

        ! Read and calculate parameters
#ifndef CRLIBM
        read(gFields(4) (1:lFields( 4)),*) tinj
        read(gFields(5) (1:lFields( 5)),*) Iinj
        read(gFields(6) (1:lFields( 6)),*) Inom
        read(gFields(7) (1:lFields( 7)),*) A
        read(gFields(8) (1:lFields( 8)),*) D
        read(gFields(9) (1:lFields( 9)),*) R
        read(gFields(10)(1:lFields(10)),*) te
#else
        tinj = round_near(errno,lFields(4)+1,gFields(4))
        if (errno.ne.0) call rounderr(errno,gFields,4,tinj)
        Iinj = round_near(errno,lFields(5)+1,gFields(5))
        if (errno.ne.0) call rounderr(errno,gFields,5,Iinj)
        Inom = round_near(errno,lFields(6)+1,gFields(6))
        if (errno.ne.0) call rounderr(errno,gFields,6,Inom)
        A    = round_near(errno,lFields(7)+1,gFields(7))
        if (errno.ne.0) call rounderr(errno,gFields,7,A)
        D    = round_near(errno,lFields(8)+1,gFields(8))
        if (errno.ne.0) call rounderr(errno,gFields,8,D)
        R    = round_near(errno,lFields(9)+1,gFields(9))
        if (errno.ne.0) call rounderr(errno,gFields,9,R)
        te   = round_near(errno,lFields(10)+1,gFields(10))
        if (errno.ne.0) call rounderr(errno,gFields,10,te)
#endif
        derivI_te = A*(te-tinj)                 ! nostore
        I_te      = (A/2.0)*(te-tinj)**2 + Iinj ! nostore
        bexp      = derivI_te/I_te
        aexp      = exp_mb(-bexp*te)*I_te
        t1        = log_mb(R/(aexp*bexp))/bexp
        I1        = aexp*exp_mb(bexp*t1)
        td        = (Inom-I1)/R + (t1 - R/(2*D))
        tnom      = td + R/D

        if (ldynkdebug) then
            write (lout,*) "DYNKDEBUG> *** PELP SETTINGS: ***"
            write (lout,*) "DYNKDEBUG> tinj =", tinj
            write (lout,*) "DYNKDEBUG> Iinj =", Iinj
            write (lout,*) "DYNKDEBUG> Inom =", Inom
            write (lout,*) "DYNKDEBUG> A    =", A
            write (lout,*) "DYNKDEBUG> D    =", D
            write (lout,*) "DYNKDEBUG> R    =", R
            write (lout,*) "DYNKDEBUG> te   =", te
            write (lout,*) "DYNKDEBUG> "
            write (lout,*) "DYNKDEBUG> derivI_te =", derivI_te
            write (lout,*) "DYNKDEBUG> I_te      =", I_te
            write (lout,*) "DYNKDEBUG> bexp      =", bexp
            write (lout,*) "DYNKDEBUG> aexp      =", aexp
            write (lout,*) "DYNKDEBUG> t1        =", t1
            write (lout,*) "DYNKDEBUG> I1        =", I1
            write (lout,*) "DYNKDEBUG> td        =", td
            write (lout,*) "DYNKDEBUG> tnom      =", tnom
            write (lout,*) "DYNKDEBUG> **********************"
        end if

        if (.not. (tinj .lt. te .and. te .lt. t1 .and. t1 .lt. td .and. td .lt. tnom)) then
            WRITE(lout,*) "DYNK> ********************************"
            WRITE(lout,*) "DYNK> ERROR***************************"
            write(lout,*) "DYNK> PELP: Order of times not correct"
            WRITE(lout,*) "DYNK> ********************************"
            call prror(51)
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

        call dynk_checkargs(nFields,5,"FUN funname ONOFF p1 p2")
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
        ! NAME
        dynk_cData(dynk_ncData)(1:lFields(2)) = gFields(2)(1:lFields(2))

        read(gFields(4)(1:lFields(4)),*) dynk_funcs(dynk_nFuncs,3) ! p1
        read(gFields(5)(1:lFields(5)),*) dynk_funcs(dynk_nFuncs,4) ! p2

        ! Check for bad input
        if (dynk_funcs(dynk_nFuncs,3) .lt. 0 .or.                     &    ! p1 <  1 ?
            dynk_funcs(dynk_nFuncs,4) .le. 1 .or.                     &    ! p2 <= 1 ?
            dynk_funcs(dynk_nFuncs,4) .lt. dynk_funcs(dynk_nFuncs,3)) then ! p2 < p1 ?

            write(lout,*) "DYNK> Error in ONOFF: Expected p1 >= 0, p2 > 1, p1 <= p2"
            call prror(-1)
        end if

    ! END CASE ONOFF

    case default
        ! UNKNOWN function
        write (lout,*) "*************************************"
        write (lout,*) "ERROR in DYNK block parsing (fort.3):"
        write (lout,*) "Unkown function to dynk_parseFUN()   "
        write (lout,*) "Got fields:"
        do ii=1,nFields
            write (lout,*) "Field(",ii,") ='",gFields(ii)(1:lFields(ii)),"'"
        end do
        write (lout,*) "*************************************"

        call dynk_dumpdata
        call prror(51)
    end select

end subroutine dynk_parseFUN

! ================================================================================================ !

subroutine dynk_checkargs(nfields,nfields_expected,funsyntax)

    use crcoall
    implicit none

    integer nfields, nfields_expected
    character(*) funsyntax
    intent(in) nfields, nfields_expected, funsyntax

    if (nfields .ne. nfields_expected) then
        write (lout,*) "ERROR in DYNK block parsing (fort.3)"
        write (lout,*) "The function expected",nfields_expected,"arguments, got",nfields
        write (lout,*) "Expected syntax:"
        write (lout,*) funsyntax(:)
        call prror(51)
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
    call alloc(dynk_fData,dynk_maxfData,0.0_fprec,"dynk_fData")
  end if
  if(cNeeded > dynk_maxcData) then
    if(cReq < 200) then
      dynk_maxcData = dynk_maxcData + 200
    else
      dynk_maxcData = cNeeded
    end if
    call alloc(dynk_cData,mStrLen,dynk_maxcData,str_dZeros,"dynk_cData")
  end if

end subroutine dynk_checkspace

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 15-10-2014
!  - Parse SET lines in the fort.3 input file.
! ================================================================================================ !
subroutine dynk_parseSET(gFields,lFields,nFields)

    use crcoall
    implicit none

    character, intent(in) :: gFields(str_maxFields)*(mStrLen)
    integer,   intent(in) :: nFields
    integer,   intent(in) :: lFields(str_maxFields)

    integer ii

    if (dynk_nSets+1 .gt. dynk_maxSets) then
        write (lout,*) "ERROR in DYNK block parsing (fort.3):"
        write (lout,*) "Maximum number of SET exceeded, please increase parameter dynk_maxSets."
        write (lout,*) "Current value of dynk_maxSets:", dynk_maxSets
        call prror(51)
    end if

    if (nFields .ne. 7) then
        write (lout,*) "ERROR in DYNK block parsing (fort.3):"
        write (lout,*) "Expected 7 fields on line while parsing SET."
        write (lout,*) "Correct syntax:"
        write (lout,*) "SET element_name attribute_name function_name startTurn endTurn turnShift"
        write (lout,*) "got field:"

        do ii=1,nFields
            write (lout,*) "Field(",ii,") ='",gFields(ii)(1:lFields(ii)),"'"
        end do

        call prror(51)
    end if

    dynk_nSets = dynk_nSets + 1

    ! function_name -> function index
    dynk_sets(dynk_nSets,1) = dynk_findFUNindex(gFields(4)(1:lFields(4)),1)
    read(gFields(5)(1:lFields(5)),*) dynk_sets(dynk_nSets,2) ! startTurn
    read(gFields(6)(1:lFields(6)),*) dynk_sets(dynk_nSets,3) ! endTurn
    read(gFields(7)(1:lFields(7)),*) dynk_sets(dynk_nSets,4) ! turnShift

    ! Sanity check on string lengths
    if (lFields(2).gt.mNameLen .or. lFields(2).gt.mStrLen-1) then
        write (lout,*) "*************************************"
        write (lout,*) "ERROR in DYNK block parsing (fort.3):"
        write (lout,*) "SET got an element name with length =",lFields(2),"> ", mNameLen, " or > mStrLen-1."
        write (lout,*) "The name was: '",gFields(2)(1:lFields(2)),"'"
        write (lout,*) "*************************************"
        call prror(51)
    end if

    if (lFields(3).gt.mStrLen-1) then
        write(lout,*) "ERROR in DYNK block parsing (fort.3) (SET):"
        write(lout,*) "The attribute name '"//gFields(2)(1:lFields(2))//"'"
        write(lout,*) "is too long! Max length is",mStrLen-1
        call prror(51)
    end if

    ! OK -- save them!
    dynk_cSets(dynk_nSets,1)(1:lFields(2)) = gFields(2)(1:lFields(2)) ! element_name
    dynk_cSets(dynk_nSets,2)(1:lFields(3)) = gFields(3)(1:lFields(3)) ! attribute_name

    ! Sanity check
    if (dynk_sets(dynk_nSets,1).eq.-1) then
        write (lout,*) "*************************************"
        write (lout,*) "ERROR in DYNK block parsing (fort.3):"
        write (lout,*) "SET wanting function '",gFields(4)(1:lFields(4)),"'"
        write (lout,*) "Calculated index:",dynk_sets(dynk_nSets,1)
        write (lout,*) "This function is not known."
        write (lout,*) "*************************************"
        call prror(51)
    end if

    if ((dynk_sets(dynk_nSets,3) .ne. -1) .and. (dynk_sets(dynk_nSets,2) .gt. dynk_sets(dynk_nSets,3))) then
        write (lout,*) "*************************************"
        write (lout,*) "ERROR in DYNK block parsing (fort.3):"
        write (lout,*) "SET got first turn num > last turn num"
        write (lout,*) "first=",dynk_sets(dynk_nSets,2)
        write (lout,*) "last =",dynk_sets(dynk_nSets,3)
        write (lout,*) "SET #", dynk_nSets
        write (lout,*) "*************************************"
        call prror(51)
    end if

    if ((dynk_sets(dynk_nSets,2) .le. 0 ) .or. &
        (dynk_sets(dynk_nSets,3) .lt. -1) .or. &
        (dynk_sets(dynk_nSets,3) .eq. 0 )) then

        write (lout,*) "*************************************"
        write (lout,*) "ERROR in DYNK block parsing (fort.3):"
        write (lout,*) "SET got turn number <= 0 "
        write (lout,*) "(not last = -1 meaning infinity)"
        write (lout,*) "first=",dynk_sets(dynk_nSets,2)
        write (lout,*) "last =",dynk_sets(dynk_nSets,3)
        write (lout,*) "SET #", dynk_nSets
        write (lout,*) "*************************************"
        call prror(51)
    end if

end subroutine dynk_parseSET

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 14-07-2015
!  - Find and return the index in the ifuncs array to the function with name funName,
!    which should be zero-padded.
!  - Return -1 if nothing was found.
!
!  Note: It is expected that the length of funName_input is equal or less than mStrLen,
!        and if it equal, that it is a zero-terminated string.
! ================================================================================================ !
integer function dynk_findFUNindex(funName_input, startfrom)

    use crcoall
    implicit none


    character(*) funName_input
    character(mStrLen) funName
    integer startfrom
    intent(in) funName_input, startfrom

    integer ii

    ! write(*,*)"DBGDBG input: '"//funName_input//"'",len(funName_input)

    if (len(funName_input).gt.mStrLen) then
        write (lout,*) "ERROR in dynk_findFUNindex"
        write (lout,*) "len(funName_input) = ",len(funName_input), &
                       ".gt. mStrLen-1 = ", mStrLen-1
        call prror(-1)
    end if

    ! If the length is exactly mStrLen, it should be zero-terminated.
    if (( len(funName_input).eq.mStrLen ) .and. &
        ( funName_input(len(funName_input):len(funName_input)) .ne.char(0)) ) then

        write (lout,*) "ERROR in dynk_findFUNindex"
        write (lout,*) "Expected funName_input[-1]=NULL"
        call prror(-1)
    end if

    do ii=1,len(funName_input)
        ! write(*,*) "DBGDBG a:", ii
        funName(ii:ii) = funName_input(ii:ii)
    end do

    funName(1:len(funName_input)) = funName_input
    do ii=len(funName_input)+1,mStrLen
        ! write(*,*) "DBGDBG b:", ii
        funName(ii:ii) = char(0)
    end do
    ! write(*,*) "DBGDBG c:", funName, len(funName)

    dynk_findFUNindex = -1

    do ii=startfrom, dynk_nFuncs
        if (dynk_cData(dynk_funcs(ii,1)).eq.funName) then
            dynk_findFUNindex = ii
            exit
        end if
    end do

end function dynk_findFUNindex

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 23-10-2014
!  - Find and return the index in the sets array to the set which matches element_name and att_name,
!    which should be zero-padded.
!  - Return -1 if nothing was found.
!
!  Note: It is expected that the length of element_name and att_name is exactly mStrLen .
! ================================================================================================ !
integer function dynk_findSETindex(element_name, att_name, startfrom)

    implicit none

    character(mStrLen) element_name, att_name
    integer startfrom
    intent(in) element_name, att_name, startfrom

    integer ii

    dynk_findSETindex = -1

    do ii=startfrom, dynk_nSets
        if ( dynk_cSets(ii,1) .eq. element_name .and. dynk_cSets(ii,2) .eq. att_name ) then
            dynk_findSETindex = ii
            exit
        end if
    end do

end function dynk_findSETindex

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 14-10-2014
!  - Check that DYNK block input in fort.3 was sane
! ================================================================================================ !
subroutine dynk_inputsanitycheck

    use crcoall
    implicit none

    integer ii, jj
    integer biggestTurn ! Used as a replacement for ending turn -1 (infinity)
    logical sane
    sane = .true.

    ! Check that there are no doubly-defined function names
    do ii=1, dynk_nFuncs-1
        jj = dynk_findFUNindex(dynk_cData(dynk_funcs(ii,1)),ii+1)
        if ( jj.ne. -1) then
            sane = .false.
            write (lout,*) "DYNK> Insane: function ",ii,"has the same name as",jj
        end if
    end do

    ! Check that no SETS work on same elem/att at same time
    biggestTurn = 1
    do ii=1, dynk_nSets
        if (dynk_sets(ii,3) .gt. biggestTurn) then
            biggestTurn = dynk_sets(ii,3)
        end if
    end do
    biggestTurn = biggestTurn+1 ! Make sure it is unique
    if (biggestTurn .le. 0) then
        ! In case of integer overflow
        write(lout,*) "FATAL ERROR: Integer overflow in dynk_inputsanitycheck!"
        call prror(-1)
    end if

    ! Do the search!
    do ii=1, dynk_nSets-1
        if (dynk_sets(ii,3).eq.-1) dynk_sets(ii,3) = biggestTurn
        ! write(*,*) "DBG: ii=",ii,dynk_cSets(ii,1)," ", dynk_cSets(ii,2)
        ! write(*,*)"DBG:", dynk_sets(ii,2),dynk_sets(ii,3)

        jj = ii
        do while (.true.)
            ! Only check SETs affecting the same elem/att
            jj = dynk_findSETindex(dynk_cSets(ii,1),dynk_cSets(ii,2),jj+1)
            ! write(*,*)" DBG: jj=",jj,dynk_cSets(jj,1)," ", dynk_cSets(jj,2)

            if (jj .eq. -1) exit ! next outer loop
            if (dynk_sets(jj,3).eq.-1) dynk_sets(jj,3) = biggestTurn

            ! write(*,*)" DBG:", dynk_sets(jj,2),dynk_sets(jj,3)
            if (dynk_sets(jj,2) .le. dynk_sets(ii,2) .and. &
                dynk_sets(jj,3) .ge. dynk_sets(ii,2)) then
                sane = .false.
                write (lout,"(A,I4,A,I8,A,I4,A,I8,A,I4,A,I8,A,I4)") &
                    " DYNK> Insane: Lower edge of SET #", jj, &
                    " =", dynk_sets(jj,2)," <= lower edge of SET #",ii, &
                    " =", dynk_sets(ii,2),"; and also higer edge of SET #",jj, &
                    " =", dynk_sets(jj,3)," >= lower edge of SET #", ii
            else if (dynk_sets(jj,3) .ge. dynk_sets(ii,3) .and. &
                     dynk_sets(jj,2) .le. dynk_sets(ii,3)) then
                sane = .false.
                write(lout, "(A,I4,A,I8,A,I4,A,I8,A,I4,A,I8,A,I4)") &
                    " DYNK> Insane: Upper edge of SET #", jj, &
                    " =", dynk_sets(jj,3)," >= upper edge of SET #",ii, &
                    " =", dynk_sets(ii,3),"; and also lower edge of SET #",jj, &
                    " =", dynk_sets(jj,2)," <= upper edge of SET #", ii
            else if (dynk_sets(jj,2) .ge. dynk_sets(ii,2) .and. &
                     dynk_sets(jj,3) .le. dynk_sets(ii,3)) then
                    ! (other way round gets caugth by the first "if")
                sane = .false.
                write(lout, "(A,I4,A,I8,A,I8,A,A,I4,A,I8,A,I8,A)") &
                    " DYNK> Insane: SET #", jj, &
                    " = (", dynk_sets(jj,2),", ", dynk_sets(jj,3), ")", &
                    " is inside SET #", ii, " = (",  &
                    dynk_sets(ii,2),", ", dynk_sets(ii,3), ")"
            end if
            if (dynk_sets(jj,3).eq.biggestTurn) dynk_sets(jj,3) = -1

        end do

        if (dynk_sets(ii,3).eq.biggestTurn) dynk_sets(ii,3) = -1

    end do

    if (.not. sane) then
        write (lout,*) "****************************************"
        write (lout,*) "******** DYNK input was insane *********"
        write (lout,*) "****************************************"
        call dynk_dumpdata
        call prror(-11)
    else if (sane .and. ldynkdebug) then
        write (lout,*) "DYNK> DYNK input was sane"
    end if

end subroutine dynk_inputsanitycheck

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 14-10-2014
!  - Dump arrays with DYNK FUN and SET data to the std. output for debugging
! ================================================================================================ !
subroutine dynk_dumpdata

    use crcoall
    implicit none

    integer ii
    write (lout,*) "**************** DYNK parser knows: ****************"
    write (lout,*) "OPTIONS:"
    write (lout,*) " ldynk            =", ldynk
    write (lout,*) " ldynkdebug       =", ldynkdebug
    write (lout,*) " ldynkfiledisable =", ldynkfiledisable

    write (lout,*) "FUN:"
    write (lout,*) "ifuncs: (",dynk_nFuncs,")"
    do ii=1,dynk_nFuncs
        write (lout,*) ii, ":", dynk_funcs(ii,:)
    end do
    write (lout,*) "dynk_iData: (",dynk_niData,")"
    do ii=1,dynk_niData
        write (lout,*) ii, ":", dynk_iData(ii)
    end do
    write (lout,*) "dynk_fData: (",dynk_nfData,")"
    do ii=1,dynk_nfData
        write (lout, '(1x,I8,1x,A,1x,E16.9)') ii, ":", dynk_fData(ii)
    end do
    write (lout,*) "dynk_cData: (",dynk_ncData,")"
    do ii=1,dynk_ncData
        write(lout,*) ii, ":", "'"//trim(stringzerotrim(dynk_cData(ii)))//"'"
    end do

    write (lout,*) "SET:"
    write (lout,*) "sets(,:) csets(,1) csets(,2): (",dynk_nSets,")"
    do ii=1,dynk_nSets
        write (lout,*) ii, ":", dynk_sets(ii,:), &
                       "'"//trim(stringzerotrim(dynk_cSets(ii,1)))// &
                       "' ", "'"//trim(stringzerotrim(dynk_cSets(ii,2)))//"'"
    end do
    write (lout,*) "dynk_cSets_unique: (",dynk_nSets_unique,")"
    do ii=1,dynk_nSets_unique
        write(lout, '(1x,I8,1x,A,1x,E16.9)') ii, ": '"// &
            trim(stringzerotrim(dynk_cSets_unique(ii,1)))//"' '"// &
            trim(stringzerotrim(dynk_cSets_unique(ii,2)))//"' = ", &
            dynk_fSets_orig(ii)
    end do

    write (lout,*) "*************************************************"

end subroutine dynk_dumpdata

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 21-10-2014
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
    integer ix
    if (ldynkdebug) then
        write(lout,*) "DYNKDEBUG> In dynk_pretrack()"
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
            element_name_s = trim(stringzerotrim(dynk_cSets_unique(dynk_nSets_unique,1)))
            att_name_s     = trim(stringzerotrim(dynk_cSets_unique(dynk_nSets_unique,2)))
            found          = .false.

            ! Special case: the element name GLOBAL-VARS (not a real element)
            ! can be used to redefine a global variable by some function.
            if (element_name_s .eq. "GLOBAL-VARS") then
                found   = .true.
                badelem = .false.

                if (att_name_s .eq. "E0") then
                    if (idp.eq.0 .or. ition.eq.0) then ! 4d tracking..
                        write(lout,*) "DYNK> Insane - attribute '", &
                                      att_name_s, "' is not valid for 'GLOBAL-VARS' ", &
                                      "when doing 4d tracking"
                        call prror(-1)
                    end if
                else
                    badelem=.true.
                end if

                if (badelem) then
                    write(lout,*) "DYNK> Insane - attribute '", &
                                  att_name_s, "' is not valid for 'GLOBAL-VARS'"
                    call prror(-1)
                end if
            end if

            do jj=1,il
                if ( bez(jj).eq. element_name_s) then

                    found = .true.

                    ! Check that the element type and attribute is supported
                    ! Check that the element can be used now
                    badelem = .false.
                    if (abs(kz(jj)).ge.1 .and. abs(kz(jj)).le.10) then !thin kicks
                        if (att_name_s .ne. "average_ms") then
                            badelem = .true.
                        end if
                    else if (abs(kz(jj)).eq.12) then ! cavity
                        if (.not. (att_name_s.eq."voltage"  .or. &
                                   att_name_s.eq."harmonic" .or. &
                                   att_name_s.eq."lag_angle")) then
                            badelem = .true.
                        end if
                        if (kp(jj).ne.6) then
                            write(lout,*) "DYNK> Insane - want to modify ", &
                                          "DISABLED RF cavity named '",element_name_s, &
                                          ". Please make sure that the voltage and ", &
                                          "harmonic number in the SINGLE ELEMENTS ", &
                                          "block is not 0!"
                            call prror(-1)
                        end if
                        if (nvar .eq. 5) then
                            write(lout,*) "DYNK> Insane - want to modify ", &
                                          "RF cavity named '", element_name_s, "', ", &
                                          "but nvars=5 (from DIFF block)."
                        end if
                    else if (abs(kz(jj)).eq.23 .or. & ! crab
                             abs(kz(jj)).eq.26 .or. & ! cc multipole,  order 2
                             abs(kz(jj)).eq.27 .or. & ! cc multipole,  order 3
                             abs(kz(jj)).eq.28) then  ! cc muiltipole, order 4
                        if (.not. (att_name_s.eq."voltage"   .or. &
                                   att_name_s.eq."frequency" .or. &
                                   att_name_s.eq."phase"     )) then
                            badelem = .true.
                        end if
                    end if

                    ! Special case:
                    ! Should the error only occur if we actually have a GLOBAL-VARS element?
                    if (bez(jj) .eq. "GLOBAL-VARS") then
                        write(lout,*) "DYNK> Insane - element found '", &
                                      "GLOBAL-VARS' is not a valid element name, ", &
                                      "it is reserved"
                        call prror(-1)
                    end if

                    if (badelem) then
                        write(lout,*) "DYNK> Insane - attribute '", &
                                      att_name_s, "' is not valid for element '", &
                                      element_name_s, "' which is of type",kz(jj)
                        call prror(-1)
                    end if
                end if
            end do

            if (.not. found) then
                write (lout,*) "DYNK> Insane: Element '",element_name_s,"' was not found"
                call prror(-1)
            end if

            ! Store original value of data point
            dynk_fSets_orig(dynk_nSets_unique) = dynk_getvalue(dynk_cSets(ii,1),dynk_cSets(ii,2))
        end if
    end do

    ! Save original values for GET functions
    do ii=1,dynk_nFuncs
        if (dynk_funcs(ii,2) .eq. 0) then ! GET
            dynk_fData(dynk_funcs(ii,3)) = dynk_getvalue( dynk_cData(dynk_funcs(ii,1)+1), &
                                           dynk_cData(dynk_funcs(ii,1)+2))
        end if
    end do

    if (ldynkdebug) call dynk_dumpdata

end subroutine dynk_pretrack

! ================================================================================================ !
!  A.Mereghetti, for the FLUKA Team
!  K.Sjobak & A. Santamaria, BE-ABP/HSS
!  Last modified: 30-10-2014
!  - Actually apply dynamic kicks
!  - Always in main code
!
!  For each element (group) flagged with SET(R), compute the new value using dynk_computeFUN() at
!  the given (shifted) turn number using the specified FUN function. The values are stored in the
!  element using dynk_setvalue().
!
!  Also resets the values at the beginning of each pass through the turn loop (for COLLIMATION).
!
!  Also writes the file "dynksets.dat", only on the first turn.
! ================================================================================================ !
subroutine dynk_apply(turn)

    use crcoall
    use mod_common
    use mod_commont
    use mod_commonmn

#ifdef COLLIMAT
    use collimation
#endif

    implicit none

#ifdef BOINC
    character(len=256) filename
#endif

    ! interface variables
    integer turn  ! current turn number
    intent(in) turn

    ! temporary variables
    integer ii, jj, shiftedTurn
    logical lopen
    real(kind=fPrec) getvaldata, newValue

    character(mStrLen) whichFUN(dynk_maxSets) ! Which function was used to set a given elem/attr?
    integer whichSET(dynk_maxSets)                   ! Which SET was used for a given elem/attr?

    ! Temp variable for padding the strings for output to dynksets.dat
    character(20) outstring_tmp1,outstring_tmp2,outstring_tmp3

    integer, parameter :: samplenumber = 1

    if ( ldynkdebug ) then
      write (lout,*) 'DYNKDEBUG> In dynk_apply(), turn = ',turn
    end if

    ! Initialize variables (every call)
    do jj=1, dynk_nSets_unique
        whichSET(jj) = -1
        do ii=1,mStrLen
            whichFUN(jj)(ii:ii) = char(0)
        end do
    end do

    ! First-turn initialization, including some parts which are specific for collimat.
    if (turn .eq. 1) then
        ! Reset RNGs and filters
        do ii=1, dynk_nFuncs
            if (dynk_funcs(ii,2) .eq. 6) then ! RANDG
                if (ldynkdebug) then
                    write (lout,*) "DYNKDEBUG> Resetting RANDG for FUN named '", &
                                   trim(stringzerotrim(dynk_cData(dynk_funcs(ii,1)))),"'"
                end if
                dynk_iData(dynk_funcs(ii,3)+3) = dynk_iData(dynk_funcs(ii,3))
                dynk_iData(dynk_funcs(ii,3)+4) = dynk_iData(dynk_funcs(ii,3)+1)
            else if (dynk_funcs(ii,2) .eq. 7) then ! RANDU
                if (ldynkdebug) then
                    write (lout,*) "DYNKDEBUG> Resetting RANDU for FUN named '", &
                                   trim(stringzerotrim(dynk_cData(dynk_funcs(ii,1)))),"'"
                end if
                dynk_iData(dynk_funcs(ii,3)+2) = dynk_iData(dynk_funcs(ii,3))
                dynk_iData(dynk_funcs(ii,3)+3) = dynk_iData(dynk_funcs(ii,3)+1)
            else if (dynk_funcs(ii,2) .eq. 8) then ! RANDON
                if (ldynkdebug) then
                    write (lout,*) "DYNKDEBUG> Resetting RANDON for FUN named '", &
                                   trim(stringzerotrim(dynk_cData(dynk_funcs(ii,1)))),"'"
                end if
                dynk_iData(dynk_funcs(ii,3)+2) = dynk_iData(dynk_funcs(ii,3))
                dynk_iData(dynk_funcs(ii,3)+3) = dynk_iData(dynk_funcs(ii,3)+1)
            else if (dynk_funcs(ii,2) .eq. 10) then ! FIR
                if (ldynkdebug) then
                    write (lout,*) "DYNKDEBUG> Resetting FIR named '", &
                                   trim(stringzerotrim(dynk_cData(dynk_funcs(ii,1)))),"'"
                end if
                do jj=0, dynk_funcs(ii,4)
                    dynk_fData(dynk_funcs(ii,3)+jj*3+1) = dynk_fData(dynk_funcs(ii,3)+jj*3+2)
                end do
            else if (dynk_funcs(ii,2) .eq. 11) then ! IIR
                if (ldynkdebug) then
                    write (lout,*) "DYNKDEBUG> Resetting IIR named '", &
                                   trim(stringzerotrim(dynk_cData(dynk_funcs(ii,1)))),"'"
                end if
                do jj=0, dynk_funcs(ii,4)
                    dynk_fData(dynk_funcs(ii,3)+jj*6+1) = dynk_fData(dynk_funcs(ii,3)+jj*6+2)
                    dynk_fData(dynk_funcs(ii,3)+jj*6+4) = dynk_fData(dynk_funcs(ii,3)+jj*6+5)
                end do
            end if

        end do ! END "do ii=1, dynk_nFuncs"

        ! Open dynksets.dat
#ifdef COLLIMAT
        if (samplenumber.eq.1) then
#endif
#ifdef CR
        ! Could have loaded a CR just before tracking starts;
        ! In this case, the dynksets is already open and positioned,
        ! so don't try to open the file again.
        if (dynkfilepos .eq.-1) then
#endif
            inquire(unit=ldynkfileunit, opened=lopen)
            if (lopen) then
                write(lout,*) "DYNK> **** ERROR in dynk_apply() ****"
                write(lout,*) "DYNK> Could not open file 'dynksets.dat'"
                call prror(-1)
            end if
#ifdef BOINC
            call boincrf("dynksets.dat",filename)
            open(unit=ldynkfileunit,file=filename,status="replace",action="write")
#else
            open(unit=ldynkfileunit,file="dynksets.dat",status="replace",action="write")
#endif

            if (ldynkfiledisable) then
                write(ldynkfileunit,*) "### DYNK file output was disabled with flag NOFILE in fort.3 ###"
            else
                write(ldynkfileunit,*) "# turn element attribute SETidx funname value"
            end if
#ifdef CR
            ! Note: To be able to reposition, each line should be shorter than 255 chars
            dynkfilepos = 1

            ! Flush the unit
            endfile(ldynkfileunit,iostat=ierro)
            backspace(ldynkfileunit,iostat=ierro)
          end if ! END if(dynkfilepos.eq.-1)
#endif
#ifdef COLLIMAT
        end if !END if(samplenumber.eq.1)

        ! Reset values to original settings in turn 1
        if (samplenumber.gt.1) then
            if (ldynkdebug) then
                write (lout,*) "DYNKDEBUG> New collimat sample, ", &
                               "samplenumber = ",samplenumber,"resetting the SET'ed values."
            end if
            do ii=1, dynk_nSets_unique
                newValue = dynk_fSets_orig(ii)
                if (ldynkdebug) then
                    write (lout,*) "DYNKDEBUG> Resetting: '", &
                                   trim(stringzerotrim(dynk_cSets_unique(ii,1))), &
                                   "':'",trim(stringzerotrim(dynk_cSets_unique(ii,2))), &
                                   "', newValue=", newValue
                end if

                call dynk_setvalue(dynk_cSets_unique(ii,1),dynk_cSets_unique(ii,2),newValue)
            end do
        end if !END "if (samplenumber.gt.1) then"
#endif
    end if ! END "if (turn .eq. 1) then"

    ! Apply the sets
    do ii=1,dynk_nSets
        ! Sanity check already confirms that only a single SET
        ! is active on a given element:attribute on a given turn.

        ! Active in this turn?
        if (turn .ge. dynk_sets(ii,2) .and. (turn .le. dynk_sets(ii,3) .or. dynk_sets(ii,3) .eq. -1)) then

            ! Shifting
            shiftedTurn = turn + dynk_sets(ii,4)

            ! Set the value
            newValue = dynk_computeFUN(dynk_sets(ii,1),shiftedTurn)
            if (ldynkdebug) then
                write (lout, '(1x,A,I5,A,I8,A,E16.9)') "DYNKDEBUG> Applying set #", ii, " on '"// &
                                trim(stringzerotrim(dynk_cSets(ii,1)))// &
                                "':'"// trim(stringzerotrim(dynk_cSets(ii,2)))// &
                                "', shiftedTurn=",shiftedTurn,", value=",newValue
            end if
            call dynk_setvalue(dynk_cSets(ii,1),dynk_cSets(ii,2),newValue)

            if (ldynkdebug) then
                getvaldata = dynk_getvalue(dynk_cSets(ii,1),dynk_cSets(ii,2))
                write (lout, '(1x,A,E16.9)') "DYNKDEBUG> Read back value = ", getvaldata

                if (getvaldata .ne. newValue) then
                    write(lout,*) "DYNKDEBUG> WARNING Read back value differs from set!"
                end if
            end if

            ! For the output file: Which function was used?
            do jj=1, dynk_nSets_unique
                if (dynk_cSets(ii,1) .eq. dynk_cSets_unique(jj,1) .and. &
                    dynk_cSets(ii,2) .eq. dynk_cSets_unique(jj,2)) then
                    whichSET(jj)=ii
                    whichFUN(jj)=dynk_cData(dynk_funcs(dynk_sets(ii,1),1))
                end if
            end do
        end if
    end do

    ! Write output file
    if (.not.ldynkfiledisable) then
#ifdef COLLIMAT
      if (samplenumber.eq.1) then
#endif
        do jj=1,dynk_nSets_unique
            getvaldata = dynk_getvalue(dynk_cSets_unique(jj,1),dynk_cSets_unique(jj,2))

            if (whichSET(jj) .eq. -1) then
                whichFUN(jj) = "N/A"
            end if

            !For compatibility with old output, the string output to dynksets.dat should be left-adjusted within each column.
            !Previously, the dynk_cSets_unique etc. strings could maximally be 20 long each.
            !Note that the length of each string is limited by the max length of element names (16), attribute names, and FUN names.
            write(outstring_tmp1,'(A20)') stringzerotrim(dynk_cSets_unique(jj,1))
            outstring_tmp1(len(outstring_tmp1)+1:) = ' ' ! Pad with trailing blanks
            write(outstring_tmp2,'(A20)') stringzerotrim(dynk_cSets_unique(jj,2))
            outstring_tmp2(len(outstring_tmp2)+1:) = ' '
            write(outstring_tmp3,'(A20)') stringzerotrim(whichFUN(jj))
            outstring_tmp3(len(outstring_tmp3)+1:) = ' '

            write(ldynkfileunit,'(I12,1x,A20,1x,A20,1x,I4,1x,A20,E16.9)') &
                 turn,outstring_tmp1,outstring_tmp2,whichSET(jj),outstring_tmp3,getvaldata
        end do

#ifdef CR
        ! Note: To be able to reposition, each line should be shorter than 255 chars
        dynkfilepos = dynkfilepos+dynk_nSets_unique
#endif
        ! Flush the unit
        endfile(ldynkfileunit,iostat=ierro)
        backspace(ldynkfileunit,iostat=ierro)
#ifdef COLLIMAT
      end if
#endif
    end if

end subroutine dynk_apply

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 17-10-2014
!  - Compute the value of a given DYNK function (funNum) for the given turn
! ================================================================================================ !
recursive real(kind=fPrec) function dynk_computeFUN(funNum, turn) result(retval)

    use crcoall
    use mod_common
    use mod_ranecu

    implicit none

    character gFields(str_maxFields)*(mStrLen)
    integer   nFields
    integer   lFields(str_maxFields)
    logical   getfields_lerr

    integer funNum, turn
    intent (in) funNum, turn

    ! Functions to call
#ifdef CRLIBM
    real(kind=fPrec) round_near
#endif

    ! Temporaries for FILELIN
    integer filelin_start, filelin_xypoints

    ! Temporaries for random generator functions
    integer tmpseed1, tmpseed2
    real(kind=fPrec) ranecu_rvec(1)

    ! General temporaries
    integer foff  ! Base offset into fexpr array
    integer ii,jj ! Loop variable

#ifdef CRLIBM
    ! String handling tempraries for PIPE, preformatting for round_near
    integer errno !for round_near
    integer nchars
    parameter(nchars=160)
    character(len=nchars) ch
#endif

! Usefull constants (pi and two)

    if (funNum .lt. 1 .or. funNum .gt. dynk_nFuncs) then
        write(lout,*) "DYNK> **** ERROR in dynk_computeFUN() ****"
        write(lout,*) "DYNK> funNum =", funNum
        write(lout,*) "DYNK> Invalid funNum, dynk_nFuncs=", dynk_nFuncs
        call dynk_dumpdata
        call prror(-1)
    end if

    select case ( dynk_funcs(funNum,2) )                              ! WHICH FUNCTION TYPE?

    case (0)                                                          ! GET
        retval = dynk_fData(dynk_funcs(funNum,3))

    case (1)                                                          ! FILE
        if (turn .gt. dynk_funcs(funNum,5) ) then
            write(lout,*)"DYNK> ****ERROR in dynk_computeFUN():FILE****"
            write(lout,*)"DYNK> funNum =", funNum, "turn=", turn
            write(lout,*)"DYNK> Turn > length of file = ", dynk_funcs(funNum,5)
            call dynk_dumpdata
            call prror(-1)
        else if (turn .lt. 1) then
            write(lout,*)"DYNK> ****ERROR in dynk_computeFUN():FILE****"
            write(lout,*)"DYNK> funNum =", funNum, "turn=", turn
            write(lout,*)"DYNK> Turn < 1, check your turn-shift!"
            call dynk_dumpdata
            call prror(-1)
        end if

        retval = dynk_fData(dynk_funcs(funNum,4)+turn-1)

    case(2)                                                           ! FILELIN
        filelin_start    = dynk_funcs(funNum,4)
        filelin_xypoints = dynk_funcs(funNum,5)
        ! Pass the correct array views/sections to dynk_lininterp
        retval = dynk_lininterp( real(turn,fPrec), &
                 dynk_fData(filelin_start:filelin_start+filelin_xypoints-1), &
                 dynk_fData(filelin_start +  filelin_xypoints: &
                            filelin_start +2*filelin_xypoints-1), &
                            filelin_xypoints )

    case(3)                                                           ! PIPE
        write(dynk_iData(dynk_funcs(funNum,3))+1,"(a,i7)") &
            "GET ID="//trim(stringzerotrim(dynk_cData(dynk_funcs(funNum,1)+3)))//" TURN=",turn
#ifndef CRLIBM
        read(dynk_iData(dynk_funcs(funNum,3)),*) retval
#else
        read(dynk_iData(dynk_funcs(funNum,3)),"(a)") ch
        call getfields_split(ch,gFields,lFields,nFields,getfields_lerr)
        if ( getfields_lerr ) then
            write(lout,*)"DYNK> ****ERROR in dynk_computeFUN():PIPE****"
            write(lout,*)"DYNK> getfields_lerr=", getfields_lerr
            call prror(-1)
        end if
        if (nFields .ne. 1) then
            write(lout,*)"DYNK> ****ERROR in dynk_computeFUN():PIPE****"
            write(lout,*)"DYNK> nFields=", nFields
            write(lout,*)"DYNK> Expected a single number."
            call prror(-1)
        end if
        retval = round_near(errno,lFields(1)+1,gFields(1))
        if (errno.ne.0) call rounderr(errno,gFields,1,retval)
#endif

    case (6)                                                          ! RANDG
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

    case (7)                                                          ! RANDU
        ! Save old seeds and load our current seeds
        call recuut(tmpseed1,tmpseed2)
        call recuin(dynk_iData(dynk_funcs(funNum,3)+2),dynk_iData(dynk_funcs(funNum,3)+3))
        ! Run generator for 1 value with mcut=-1
        call ranecu( ranecu_rvec, 1, -1 )
        ! Save our current seeds and load old seeds
        call recuut(dynk_iData(dynk_funcs(funNum,3)+2),dynk_iData(dynk_funcs(funNum,3)+3))
        call recuin(tmpseed1,tmpseed2)
        retval = ranecu_rvec(1)

    case (8)                                                         ! RANDON
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

    case(10)                                                          ! FIR
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

    case(11)                                                          ! IIR
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

    case (20)                                                         ! ADD
        retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) &
               + dynk_computeFUN(dynk_funcs(funNum,4),turn)

    case (21)                                                         ! SUB
        retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) &
               - dynk_computeFUN(dynk_funcs(funNum,4),turn)

    case (22)                                                         ! MUL
        retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) &
               * dynk_computeFUN(dynk_funcs(funNum,4),turn)

    case (23)                                                         ! DIV
        retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) &
               / dynk_computeFUN(dynk_funcs(funNum,4),turn)

    case (24)                                                         ! POW
        retval = dynk_computeFUN(dynk_funcs(funNum,3),turn) &
              ** dynk_computeFUN(dynk_funcs(funNum,4),turn)

    case (30)                                                         ! MINUS
        retval = (-1)*dynk_computeFUN(dynk_funcs(funNum,3),turn)

    case (31)                                                         ! SQRT
        retval = sqrt(dynk_computeFUN(dynk_funcs(funNum,3),turn))

    case (32)                                                         ! SIN
        retval = sin_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

    case (33)                                                         ! COS
        retval = cos_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

    case (34)                                                         ! LOG
        retval = log_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

    case (35)                                                         ! LOG10
        retval = log10_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

    case (36)                                                         ! EXP
        retval = exp_mb(dynk_computeFUN(dynk_funcs(funNum,3),turn))

    case (40)                                                         ! CONST
        retval = dynk_fData(dynk_funcs(funNum,3))

    case (41)                                                         ! TURN
        retval = turn

    case (42)                                                         ! LIN
        retval = turn*dynk_fData(dynk_funcs(funNum,3)) +  &
                      dynk_fData(dynk_funcs(funNum,3)+1)

    case (43)                                                         ! LINSEG
        filelin_start    = dynk_funcs(funNum,3)
        filelin_xypoints = 2
        ! Pass the correct array views/sections to dynk_lininterp
        retval = dynk_lininterp( real(turn,fPrec), &
                 dynk_fData(filelin_start:filelin_start+1), &
                 dynk_fData(filelin_start+2:filelin_xypoints+3), &
                 filelin_xypoints )

    case (44,45)                                                      ! QUAD/QUADSEG
        retval = (turn**2)*dynk_fData(dynk_funcs(funNum,3))   + ( &
                      turn*dynk_fData(dynk_funcs(funNum,3)+1) + &
                           dynk_fData(dynk_funcs(funNum,3)+2) )

      case (60)                                                         ! SINF
        retval = dynk_fData(dynk_funcs(funNum,3)) &
               * SIN_MB( dynk_fData(dynk_funcs(funNum,3)+1) * turn  &
                       + dynk_fData(dynk_funcs(funNum,3)+2) )


    case (61)                                                         ! COSF
        retval = dynk_fData(dynk_funcs(funNum,3)) &
               * COS_MB( dynk_fData(dynk_funcs(funNum,3)+1) * turn  &
                       + dynk_fData(dynk_funcs(funNum,3)+2) )

    case (62)                                                         ! COSF_RIPP
        retval = dynk_fData(dynk_funcs(funNum,3)) &
               * COS_MB( (two*pi)*real(turn-1,fPrec)/dynk_fData(dynk_funcs(funNum,3)+1) &
                       + dynk_fData(dynk_funcs(funNum,3)+2) )

    case (80)                                                         ! PELP
        foff = dynk_funcs(funNum,3)
        if (turn .le. dynk_fData(foff)) then ! <= tinj
            ! Constant Iinj
            retval = dynk_fData(foff+5)
        else if (turn .le. dynk_fData(foff+1)) then ! <= te
            ! Parabola (accelerate)
            retval = ( dynk_fData(foff+6)*(turn-dynk_fData(foff))**2 ) / 2.0 + dynk_fData(foff+5)
        else if (turn .le. dynk_fData(foff+2)) then ! <= t1
            ! Exponential
            retval = dynk_fData(foff+7)*exp_mb( dynk_fData(foff+8)*turn ) !!!!!! WTF ????? EXP, should be EXP_MB !!!!!!!!!!!!!
        else if (turn .le. dynk_fData(foff+3)) then ! <= td
            ! Linear (max ramp rate)
            retval = dynk_fData(foff+10)*(turn-dynk_fData(foff+2)) + dynk_fData(foff+9)
        else if (turn .le. dynk_fData(foff+4)) then ! <= tnom
            ! Parabola (decelerate)
            retval =  - ( (dynk_fData(foff+11) * &
                          (dynk_fData(foff+4)-turn)**2) ) / 2.0 &
                         + dynk_fData(foff+12)
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
        write(lout,*) "DYNK> **** ERROR in dynk_computeFUN(): ****"
        write(lout,*) "DYNK> funNum =", funNum, "turn=", turn
        write(lout,*) "DYNK> Unknown function type ", dynk_funcs(funNum,2)
        call dynk_dumpdata
        call prror(-1)
    end select

end function dynk_computeFUN

! ================================================================================================ !
!  A.Santamaria & K.Sjobak, BE-ABP/HSS
!  Last modified: 31-10-2014
!  - Set the value of the element's attribute
! ================================================================================================ !
subroutine dynk_setvalue(element_name, att_name, newValue)

    use scatter, only : scatter_ELEM_scale, scatter_elemPointer
    use crcoall
    use mod_common
    use mod_commont
    use mod_commonmn
    use elens
    implicit none

    character(mStrLen) element_name, att_name
    real(kind=fPrec) newValue
    intent (in) element_name, att_name, newValue

    ! Temp variables
    integer el_type, ii, j
    character(mStrLen) element_name_stripped
    character(mStrLen) att_name_stripped

    !Original energies before energy update
    real(kind=fPrec) e0fo, e0o

    ! For sanity check
    logical ldoubleElement
    ldoubleElement = .false.

    element_name_stripped = trim(stringzerotrim(element_name))
    att_name_stripped = trim(stringzerotrim(att_name))

    if (ldynkdebug) then
        write (lout, '(1x,A,E16.9)') "DYNKDEBUG> In dynk_setvalue(), element_name = '"// &
                                     trim(element_name_stripped)//"', att_name = '"// &
                                     trim(att_name_stripped)//"', newValue =", newValue
    end if

    ! Here comes the logic for setting the value of the attribute for all instances of the element.

    ! Special non-physical elements
    if (element_name_stripped .eq. "GLOBAL-VARS") then
        if (att_name_stripped .eq. "E0" ) then
            ! Modify the reference particle

            e0o  = e0
            e0fo = e0f

            e0     = newValue
            e0f    = sqrt(e0**2 - nucm0**2)
            gammar = nucm0/e0

            ! Modify the Energy
            do j = 1, napx
                dpsv(j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
                dpsv1(j) = (dpsv(j)*c1e3)/(one + dpsv(j))
                dpd(j) = one + dpsv(j)
                dpsq(j) = sqrt(dpd(j))
                oidpsv(j) = one/(one + dpsv(j))
                moidpsv(j) = mtc(j)/(one + dpsv(j))
                rvv(j) = (ejv(j)*e0f)/(e0*ejfv(j))

                !Also update sigmv with the new beta0 = e0f/e0
                sigmv(j)=((e0f*e0o)/(e0fo*e0))*sigmv(j)
              end do
              if(ithick.eq.1) call synuthck
        end if
        ldoubleElement = .true.
    end if

    ! Normal SINGLE ELEMENTs
    do ii=1,il
        ! TODO: Here one could find the right ii in dynk_pretrack,
        ! and then avoid this loop / string-comparison
        if (element_name_stripped.eq.bez(ii)) then ! name found
            el_type=kz(ii)      ! type found

            if (ldoubleElement) then ! Sanity check
                write(lout,*) "DYNK> ERROR: two elements with the same BEZ?"
                call prror(-1)
            end if
            ldoubleElement = .true.

            if ((abs(el_type).eq.1).or. & ! horizontal bending kick
                (abs(el_type).eq.2).or. & ! quadrupole kick
                (abs(el_type).eq.3).or. & ! sextupole kick
                (abs(el_type).eq.4).or. & ! octupole kick
                (abs(el_type).eq.5).or. & ! decapole kick
                (abs(el_type).eq.6).or. & ! dodecapole kick
                (abs(el_type).eq.7).or. & ! 14th pole kick
                (abs(el_type).eq.8).or. & ! 16th pole kick
                (abs(el_type).eq.9).or. & ! 18th pole kick
                (abs(el_type).eq.10)) then ! 20th pole kick

                if (att_name_stripped.eq."average_ms") then !
                    ed(ii) = newValue
                else
                    goto 100 ! ERROR
                end if
                call initialize_element(ii, .false.)

            ! Not yet supported
            ! else if (abs(el_type).eq.11) then ! MULTIPOLES
            !     if (att_name_stripped.eq."bending_str") then
            !         ed(ii) = newValue
            !     else
            !         goto 100 ! ERROR
            !     end if
            !     call initialize_element(ii, .false.)

            ! Cavities
            else if (abs(el_type).eq.12) then
                if (att_name_stripped.eq."voltage") then ! [MV]
                    ed(ii) = newValue
                else if (att_name_stripped.eq."harmonic") then !
                    ek(ii) = newValue
                    el(ii) = dynk_elemdata(ii,3) ! Need to reset el before calling initialize_element()
                    call initialize_element(ii, .false.)
                else if (att_name_stripped.eq."lag_angle") then ! [deg]
                    el(ii) = newValue
                    ! Note: el is set to 0 in initialize_element and in daten.
                    ! Calling initialize element on a cavity without setting el
                    ! will set phasc = 0!
                    call initialize_element(ii, .false.)
                else
                    goto 100 ! ERROR
                end if

            ! Not yet supported
            ! AC dipole
            ! else if (abs(el_type).eq.16) then
            !     if (att_name_stripped.eq."amplitude") then ! [T.m]
            !         ed(ii) = dynk_computeFUN(funNum,turn)
            !     else if (att_name_stripped.eq."frequency") then ! [2pi]
            !         ek(ii) = dynk_computeFUN(funNum,turn)
            !     else if (att_name_stripped.eq."phase") then ! [rad]
            !         el(ii) = dynk_computeFUN(funNum,turn)
            !     else
            !         goto 100 ! ERROR
            !     end if

            ! Not yet supported
            ! beam-beam separation
            ! else if (abs(el_type).eq.20) then
            !     if (att_name_stripped.eq."horizontal") then ! [mm]
            !         ed(ii) = dynk_computeFUN(funNum,turn)
            !     else if (att_name_stripped.eq."vertical") then ! [mm]
            !         ek(ii) = dynk_computeFUN(funNum,turn)
            !     else if (att_name_stripped.eq."strength") then ! [m]
            !         el(ii) = dynk_computeFUN(funNum,turn)
            !     else
            !         goto 100 ! ERROR
            !     end if

            else if ((abs(el_type).eq.23).or. &  ! crab cavity
                     (abs(el_type).eq.26).or. &  ! cc mult. kick order 2
                     (abs(el_type).eq.27).or. &  ! cc mult. kick order 3
                     (abs(el_type).eq.28)) then  ! cc mult. kick order 4
                if (att_name_stripped.eq."voltage") then ! [MV]
                    ed(ii) = newValue
                else if (att_name_stripped.eq."frequency") then ! [MHz]
                    ek(ii) = newValue
                else if (att_name_stripped.eq."phase") then ! [rad]
                    ! Note: el is set to 0 in initialize_element and in daten.
                    ! Calling initialize element on a crab without setting el
                    ! will set crabph = 0!
                    el(ii) = newValue
                    call initialize_element(ii, .false.)
                else
                    goto 100 ! ERROR
                end if

            ! Electron lens
            else if (el_type.eq.29) then
                if (att_name_stripped.eq."theta_r2") then ! [mrad]
                    elens_theta_r2(ii) = newValue
                else
                    goto 100 ! ERROR
                end if

            ! Scatter
            else if (el_type.eq.40) then
                if(att_name_stripped.eq."scaling") then
                    scatter_ELEM_scale(scatter_elemPointer(ii)) = newValue
                else
                    goto 100 ! ERROR
                end if

            else
                WRITE (lout,*) "DYNK> *** ERROR in dynk_setvalue() ***"
                write (lout,*) "DYNK> Unsupported element type", el_type
                write (lout,*) "DYNK> element name = '",element_name_stripped,"'"
                call prror(-1)
            end if
        end if
    end do

    ! Sanity check
    if (.not.ldoubleElement) then
        goto 101
    end if

    return

    ! Error handlers
100 continue
        write (lout,*)"DYNK> *** ERROR in dynk_setvalue() ***"
        write (lout,*)"DYNK> Attribute'", att_name_stripped, &
                      "' does not exist for type =", el_type
        call prror(-1)

101 continue
        write (lout,*)"DYNK> *** ERROR in dynk_setvalue() ***"
        write (lout,*)"DYNK> The element named '",element_name_stripped, &
                             "' was not found."
        call prror(-1)

end subroutine dynk_setvalue

! ================================================================================================ !
!  A.Santamaria & K. Sjobak, BE-ABP/HSS
!  Last modified: 2101-2015
!  - Returns the original value currently set by an element.
!
!  Note: Expects that arguments element_name and att_name are zero-terminated strings of
!        length mStrLen!
! ================================================================================================ !
real(kind=fPrec) function dynk_getvalue(element_name, att_name)

    use scatter, only : scatter_ELEM_scale, scatter_elemPointer
    use crcoall
    use mod_common
    use mod_commont
    use mod_commonmn
    use elens
    implicit none

    character(mStrLen) element_name, att_name
    intent(in) element_name, att_name

    integer el_type, ii
    character(mStrLen) element_name_s, att_name_s

    logical ldoubleElement
    ldoubleElement = .false.  ! For sanity check

    element_name_s = trim(stringzerotrim(element_name))
    att_name_s = trim(stringzerotrim(att_name))

    if (ldynkdebug) then
        write(lout,*) "DYNKDEBUG> In dynk_getvalue(), element_name = '"// &
                      trim(element_name_s)//"', att_name = '"//trim(att_name_s)//"'"
    end if

    ! Special non-physical elements
    if (element_name_s .eq. "GLOBAL-VARS") then
        if (att_name_s .eq. "E0" ) then
            ! Return the energy
            dynk_getvalue = e0
        end if
        ldoubleElement = .true.
    end if

    ! Normal SINGLE ELEMENTs
    do ii=1,il
        ! TODO: Here one could find the right ii in dynk_pretrack,
        ! and then avoid this loop / string-comparison
        if (element_name_s.eq.bez(ii)) then ! name found
            el_type=kz(ii)
            if (ldoubleElement) then
               write (lout,*) "DYNK> ERROR: two elements with the same BEZ"
               call prror(-1)
            end if
            ldoubleElement = .true.

            ! Nonlinear elements
            if ((abs(el_type).eq.1).or. &
                (abs(el_type).eq.2).or. &
                (abs(el_type).eq.3).or. &
                (abs(el_type).eq.4).or. &
                (abs(el_type).eq.5).or. &
                (abs(el_type).eq.6).or. &
                (abs(el_type).eq.7).or. &
                (abs(el_type).eq.8).or. &
                (abs(el_type).eq.9).or. &
                (abs(el_type).eq.10)) then
                if (att_name_s.eq."average_ms") then
                    dynk_getvalue = ed(ii)
                else
                    goto 100 ! ERROR
                end if

            ! Multipoles (Not yet supported)
            ! else if (abs(el_type).eq.11) then
            !     if (att_name_s.eq."bending_str") then
            !         dynk_getvalue = dynk_elemdata(ii,2)
            !     elseif (att_name_s.eq."radius") then
            !         dynk_getvalue = dynk_elemdata(ii,3)
            !     else
            !         goto 100 ! ERROR
            !     end if

            ! Cavities
            else if (abs(el_type).eq.12) then
                if (att_name_s.eq."voltage"  ) then ! MV
                    dynk_getvalue = ed(ii)
                else if (att_name_s.eq."harmonic" ) then ! harmonic number
                    dynk_getvalue = ek(ii)
                else if (att_name_s.eq."lag_angle") then ! [deg]
                    dynk_getvalue = dynk_elemdata(ii,3)
                else
                    goto 100 ! ERROR
                end if

            ! Not yet supported
            ! AC dipole
            ! else if (abs(el_type).eq.16) then
            !     if (att_name_s.eq."amplitude") then ! [T.m]
            !         nretdata = nretdata+1
            !         retdata(nretdata) = ed(ii)
            !     else if (att_name_s.eq."frequency") then !  [2pi]
            !         nretdata = nretdata+1
            !         retdata(nretdata) = ek(ii)
            !     else if (att_name_s.eq."phase") then !  [rad]
            !         nretdata = nretdata+1
            !         retdata(nretdata) = el(ii)
            !     else
            !         goto 100 ! ERROR
            !     end if

            ! Not yet supported
            ! beam-beam separation
            ! else if (abs(el_type).eq.20) then
            !     if (att_name_s.eq."horizontal") then ! [mm]
            !         nretdata = nretdata+1
            !         retdata(nretdata) = ed(ii)
            !     else if (att_name_s.eq."vertical") then ! [mm]
            !         nretdata = nretdata+1
            !         retdata(nretdata) = ek(ii)
            !     else if (att_name_s.eq."strength") then ! [m]
            !         nretdata = nretdata+1
            !         retdata(nretdata) = el(ii)
            !     else
            !         goto 100 ! ERROR
            !     end if

            else if ((abs(el_type).eq.23).or. & ! crab cavity
                     (abs(el_type).eq.26).or. & ! cc mult. kick order 2
                     (abs(el_type).eq.27).or. & ! cc mult. kick order 3
                     (abs(el_type).eq.28)) then ! cc mult. kick order 4
                if (att_name_s.eq."voltage") then ! [MV]
                    dynk_getvalue = ed(ii)
                else if (att_name_s.eq."frequency") then ! [MHz]
                    dynk_getvalue = ek(ii)
                else if (att_name_s.eq."phase") then ! [rad]
                    if (abs(el_type).eq.23) then
                        dynk_getvalue = crabph(ii)
                    else if (abs(el_type).eq.26) then
                        dynk_getvalue = crabph2(ii)
                    else if (abs(el_type).eq.27) then
                        dynk_getvalue = crabph3(ii)
                    else if (abs(el_type).eq.28) then
                        dynk_getvalue = crabph4(ii)
                    end if
                else
                    goto 100 ! ERROR
                end if

            ! Electron lens
            else if (el_type.eq.29) then
                if(att_name_s.eq."theta_r2") then ! [mrad]
                    dynk_getvalue = elens_theta_r2(ii)
                else
                    goto 100 ! ERROR
                end if

            ! Scatter
            else if (el_type.eq.40) then
                if(att_name_s.eq."scaling") then
                    dynk_getvalue = scatter_ELEM_scale(scatter_elemPointer(ii))
                else
                    goto 100 ! ERROR
                end if

            end if ! el_type
        end if ! bez
    end do

    if (ldynkdebug) then
        write(lout,*) "DYNKDEBUG> In dynk_getvalue(), returning =", dynk_getvalue
    end if

    return

    ! Error handlers
100 continue
        write(lout,*) "DYNK> *** ERROR in dynk_getvalue() ***"
        write(lout,*) "DYNK> Unknown attribute '", trim(att_name_s),"'", &
                      " for type",el_type," name '", trim(bez(ii)), "'"

        call prror(-1)

end function dynk_getvalue

! ================================================================================================ !
!  A.Mereghetti, for the FLUKA Team and K.Sjobak for BE-ABP/HSS
!  Last modified: 29-10-2014
!
!  - Define a linear function with a set of x,y-coordinates xvals, yvals
!  - Return this function evaluated at the point x.
!  - The length of the arrays xvals and yvals should be given in datalen.
!
!  - xvals should be in increasing order, if not then program is aborted.
!  - If x < min(xvals) or x>max(xvals), program is aborted.
!  - If datalen <= 0, program is aborted.
! ================================================================================================ !
real(kind=fPrec) function dynk_lininterp(x,xvals,yvals,datalen)

      use crcoall
    implicit none


    real(kind=fPrec) x, xvals(*),yvals(*)
    integer datalen
    intent(in) x,xvals,yvals,datalen

    integer ii
    real(kind=fPrec) dydx, y0

    ! Sanity checks
    if (datalen .le. 0) then
        write(lout,*) "DYNK> **** ERROR in dynk_lininterp() ****"
        write(lout,*) "DYNK> datalen was 0!"
        call prror(-1)
    end if
    if (x .lt. xvals(1) .or. x .gt. xvals(datalen)) then
        write(lout,*) "DYNK> **** ERROR in dynk_lininterp() ****"
        write(lout,*) "x =",x, "outside range", xvals(1),xvals(datalen)
        call prror(-1)
    end if

    ! Find the right indexes i1 and i2
    ! Special case: first value at first point
    if (x .eq. xvals(1)) then
        dynk_lininterp = yvals(1)
        return
    end if

    do ii=1, datalen-1
        if (xvals(ii) .ge. xvals(ii+1)) then
            write (lout,*) "DYNK> **** ERROR in dynk_lininterp() ****"
            write (lout,*) "DYNK> xvals should be in increasing order"
            write (lout,*) "DYNK> xvals =", xvals(:datalen)
            call prror(-1)
        end if

        if (x .le. xvals(ii+1)) then
            ! We're in the right interval
            dydx = (yvals(ii+1)-yvals(ii)) / (xvals(ii+1)-xvals(ii))
            y0   = yvals(ii) - dydx*xvals(ii)
            dynk_lininterp = dydx*x + y0
            return
        end if
    end do

    ! We didn't return yet: Something wrong
    write (lout,*) "DYNK> ****ERROR in dynk_lininterp() ****"
    write (lout,*) "DYNK> Reached the end of the function"
    write (lout,*) "DYNK> This should not happen, please contact developers"
    call prror(-1)

end function dynk_lininterp

! ================================================================================================ !
!  K. Sjobak, ABP-HSS,
!  Last modified: 23-01-2015
!  - Indicates whether a structure element is in use by DYNK
! ================================================================================================ !
logical function dynk_isused(i)

    use crcoall
    use mod_common
    implicit none

    integer, intent(in) :: i
    integer ix,k
    character(mStrLen) element_name_stripped

    ! Sanity check
    if (i .gt. iu .or. i .le. 0) then
        write (lout,*) "Error in dynk_isused(): i=",i,"out of range"
        call prror(-1)
    end if
    ix = ic(i)-nblo
    if (i .le. 0) then
        write (lout,*) "Error in dynk_isused(): ix-nblo=",ix,"is a block?"
        call prror(-1)
    end if

    do k=1,dynk_nSets
        element_name_stripped = trim(stringzerotrim(dynk_cSets(k,1)))
        if (bez(ix) .eq. element_name_stripped) then
            dynk_isused = .true.
            if (ldynkdebug) then
                write(lout,*) "DYNKDEBUG> dynk_isused = TRUE, bez='"//bez(ix)// &
                              "', element_name_stripped='"//element_name_stripped//"'"
            end if
            return
        end if
    end do

    if (ldynkdebug) then
        write(lout,*) "DYNKDEBUG> dynk_isused = FALSE, bez='"//bez(ix)//"'"
    end if

    dynk_isused = .false.
    return

end function dynk_isused

! =================================================================================================
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 18-11-2017
!  Initialize the memory for DYNK
! =================================================================================================
subroutine dynk_comnul

  use numerical_constants, only : zero

  implicit none

  integer i,j

  ldynk = .false.
  ldynkdebug = .false.
  ldynkfiledisable = .false.

  dynk_nFuncs = 0
  dynk_niData = 0
  dynk_nfData = 0
  dynk_ncData = 0

  do i=1,dynk_maxFuncs
    dynk_funcs(i,1)= 0 !FUN name ( index in dynk_cData; 0 is invalid )
    dynk_funcs(i,2)=-1 !FUN type (-1 is invalid)
    dynk_funcs(i,3)= 0
    dynk_funcs(i,4)= 0
    dynk_funcs(i,5)= 0
  end do

  dynk_nSets = 0

  do i=1, dynk_maxSets
    dynk_sets(i,1) = 0 !FUN idx ( index in dynk_funcs; 0 is invalid )
    dynk_sets(i,2) = 0
    dynk_sets(i,3) = 0
    dynk_sets(i,4) = 0

    do j=1, mStrLen
      dynk_cSets(i,1)(j:j) = char(0)
      dynk_cSets(i,2)(j:j) = char(0)
      dynk_cSets_unique(i,1)(j:j) = char(0)
      dynk_cSets_unique(i,2)(j:j) = char(0)
    end do
    dynk_fSets_orig(i) = zero
  end do

  do i=1,nele
    dynk_izuIndex(i) = 0
    dynk_elemdata(i,1) = 0
    dynk_elemdata(i,2) = 0
    dynk_elemdata(i,3) = 0
  end do
#ifdef CR
  dynkfilepos = -1 ! This line counter becomes >= 0 once the file is opened.
#endif
end subroutine dynk_comnul

! =================================================================================================
!  A. Mereghetti, for the FLUKA Team
!  Last Modified: 2018-04-17
!  Close units for logging dynks
! =================================================================================================
subroutine dynk_closeFiles

  implicit none

  integer i
  logical lopen

  if(.not. ldynk) return

  do i=1,dynk_nFuncs
    if (dynk_funcs(i,2) == 3) then ! PIPE FUN
      ! InPipe
      inquire(unit=dynk_iData(dynk_funcs(i,3)), opened=lopen)
      if(lopen) close(dynk_iData(dynk_funcs(i,3)))

      ! OutPipe
      inquire(unit=dynk_iData(dynk_funcs(i,3)+1), opened=lopen)
      if(lopen) then
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

  read(fileunit,err=100,end=100) dynkfilepos_cr, dynk_niData_cr, dynk_nfData_cr, dynk_ncData_cr

  call alloc(dynk_iData_cr,dynk_niData_cr,0,"dynk_iData_cr")
  call alloc(dynk_fData_cr,dynk_nfData_cr,zero,"dynk_fData_cr")
  call alloc(dynk_cData_cr,mStrLen,dynk_ncData_cr,str_dZeros,"dynk_cData_cr")

  read(fileunit,err=100,end=100) &
    (dynk_iData_cr(j),j=1,dynk_niData_cr), (dynk_fData_cr(j),j=1,dynk_nfData_cr), &
    (dynk_cData_cr(j),j=1,dynk_ncData_cr), (dynk_fSets_cr(j),j=1,dynk_maxSets)

  readerr=.false.
  return

100 continue

  write(lout,*) "READERR in scatter_crcheck; fileunit=",fileunit
  write(93,*)   "READERR in scatter_crcheck; fileunit=",fileunit
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


    logical lopen
    integer ierro
#ifdef BOINC
    character(len=256) filename
#endif
    integer j
    character(len=1024) arecord

    inquire(unit=ldynkfileunit, opened=lopen)
    if (lopen) then
        write(93,*) "SIXTRACR CRCHECK FAILED while repositioning 'dynksets.dat'"
        write(93,*) "Could not open file!"
        endfile (93,iostat=ierro)
        backspace (93,iostat=ierro)

        write(lout,*) "SIXTRACR CRCHECK failure positioning 'dynksets.dat'"
        call prror(-1)
    end if

    if (dynkfilepos_cr .ne. -1) then
#ifdef BOINC
        call boincrf("dynksets.dat",filename)
        open(unit=ldynkfileunit,file=filename,status="old",action="readwrite",err=110)
#else
        open(unit=ldynkfileunit,file='dynksets.dat',status="old",action="readwrite",err=110)
#endif
        dynkfilepos = 0     ! Start counting lines at 0, not -1
        do j=1,dynkfilepos_cr
            read(ldynkfileunit,'(a1024)',end=110,err=110,iostat=ierro) arecord
            dynkfilepos=dynkfilepos+1
        end do

        endfile(ldynkfileunit,iostat=ierro)
        close(ldynkfileunit)
#ifdef BOINC
        call boincrf("dynksets.dat",filename)
        open(unit=ldynkfileunit,file=filename,status="old",position='append',action="write")
#else
        open(unit=ldynkfileunit,file="dynksets.dat",status="old",position='append',action="write")
#endif

        write(93,*) "SIXTRACR CRCHECK sucessfully repositioned 'dynksets.dat', "// &
                    "dynkfilepos=",dynkfilepos, "dynkfilepos_cr=",dynkfilepos_cr
        endfile (93,iostat=ierro)
        backspace (93,iostat=ierro)
    else
        write(93,*) "SIXTRACR CRCHECK did not attempt repositioning "// &
                    "of dynksets.dat, dynkfilepos_cr=",dynkfilepos_cr
        write(93,*) "If anything has been written to the file, "// &
                    "it will be correctly truncated in dynk_apply on the first turn."
        endfile (93,iostat=ierro)
        backspace (93,iostat=ierro)
    end if

    return

110 continue
    write(93,*) "SIXTRACR CRCHECK *** ERROR *** reading 'dynksets.dat', iostat=",ierro
    write(93,*) "dynkfilepos=",dynkfilepos," dynkfilepos_cr=",dynkfilepos_cr
    endfile   (93,iostat=ierro)
    backspace (93,iostat=ierro)
    write(lout,*) "SIXTRACR CRCHECK failure positioning 'dynksets.dat'"
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
  logical, intent(out)   :: fileerror
  integer, intent(inout) :: ierro

  integer j

  write(fileunit,err=100,iostat=ierro) dynkfilepos, dynk_niData, dynk_nfData, dynk_ncData
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

  call dealloc(dynk_iData_cr,"dynk_iData_cr")
  call dealloc(dynk_fData_cr,"dynk_fData_cr")
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
