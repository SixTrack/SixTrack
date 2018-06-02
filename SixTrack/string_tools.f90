! ================================================================================================ !
!  SixTrack String Tools
! ~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
!
!  Old method: getfields_split, stringzerotrim
!  A. Mereghetti, for the FLUKA Team
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-19
!
!  Note: Careful with adding use <module> to this module as it easily creates circular dependencies
! ================================================================================================ !
module string_tools
  
  use strings
  
  implicit none
  
  ! "Standard" string and name length, +1 for \0
  integer, parameter :: str_maxName   = 48
  integer, parameter :: str_maxLen    = 161
  integer, parameter :: str_maxFields = 15
  
  ! Dummy empty strings
  character(len=str_maxLen),  parameter :: str_dSpace  = repeat(" ",str_maxLen)
  character(len=str_maxLen),  parameter :: str_dZeros  = repeat(char(0),str_maxLen)
  character(len=str_maxName), parameter :: str_nmSpace = repeat(" ",str_maxName)
  character(len=str_maxName), parameter :: str_nmZeros = repeat(char(0),str_maxName)
  
  public str_strip, chr_strip, chr_trimZero
  public str_stripQuotes, chr_stripQuotes
  public str_sub, chr_expandBrackets
  public chr_padZero, chr_padSpace
  public str_inStr, chr_inStr
  
  interface str_cast
    module procedure str_toReal
    module procedure str_toInt
    module procedure str_toLog
  end interface str_cast
  
  interface chr_cast
    module procedure chr_toReal
    module procedure chr_toInt
    module procedure chr_toLog
  end interface chr_cast
  
  !
  ! Old stuff added for backwards compatibility
  !
  
  integer, parameter :: getfields_n_max_fields = str_maxFields ! Max number of returned fields
  integer, parameter :: getfields_l_max_string = str_maxLen    ! Max string length
  integer, parameter :: stringzerotrim_maxlen  = str_maxLen    ! Max string length
  
  public stringzerotrim
  
contains

! ================================================================================================ !
!  String Split Routine
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-18
!  Splits a string into an array of strings by one or more space,
!    but not within a pair of single our double quotes.
!  Note: Copying data from one string to another string, by slicing the char array, does not work
!        with the intel compiler unless -assume realloc-lhs is enabled.
! ================================================================================================ !
subroutine str_split(toSplit, sArray, nArray, hasErr, fixLen, nIndent)
  
  use crcoall
  
  type(string),               intent(in)  :: toSplit
  type(string), allocatable,  intent(out) :: sArray(:)
  integer,                    intent(out) :: nArray
  logical,                    intent(out) :: hasErr
  integer,      optional,     intent(in)  :: fixLen
  integer,      optional,     intent(out) :: nIndent
  
  character(len=:), allocatable :: mskSplit, tVal
  integer mLen, nVal, nInd
  integer i, iVal, iSt
  logical hErr
  
  call chr_scanString(toSplit%chr, mskSplit, mLen, nVal, nInd, hErr)
  if(present(fixLen))  mLen    = fixLen
  if(present(nIndent)) nIndent = nInd
  nArray = nVal
  hasErr = hErr
  
  if(nVal == 0) return
  
  if(allocated(sArray)) deallocate(sArray)
  allocate(sArray(nVal))
  
  iVal = 0
  iSt  = 0
  tVal = ""
  do i=1,len(mskSplit)
    if(mskSplit(i:i) == "X") then
      iSt = iSt + 1
      if(iSt > mLen) then
        write(lout,"(2(a,i0))") "SPLIT> ERROR Split element ",iVal," is longer than the buffer of ",mLen,"."
        hasErr = .true.
        exit
      end if
      tVal = tVal//toSplit%chr(i:i)
    else
      if(iSt > 0) then
        iSt  = 0
        iVal = iVal + 1
        if(iVal > nVal) exit
        sArray(iVal) = tVal
        tVal = ""
      end if
    end if
  end do
  
end subroutine str_split

subroutine chr_split(toSplit, sArray, nArray, hasErr, fixLen, nIndent)
  
  use crcoall
  
  character(len=*),              intent(in)  :: toSplit
  character(len=:), allocatable, intent(out) :: sArray(:)
  integer,                       intent(out) :: nArray
  logical,                       intent(out) :: hasErr
  integer,          optional,    intent(in)  :: fixLen
  integer,          optional,    intent(out) :: nIndent
  
  character(len=:), allocatable :: mskSplit
  integer mLen, nVal, nInd
  integer i, iVal, iSt
  logical hErr
  
  call chr_scanString(toSplit, mskSplit, mLen, nVal, nInd, hErr)
  if(present(fixLen))  mLen    = fixLen
  if(present(nIndent)) nIndent = nInd
  nArray = nVal
  hasErr = hErr
  
  if(nVal == 0) return
  
  if(allocated(sArray)) deallocate(sArray)
  allocate(character(len=mLen) :: sArray(nVal))
  do i=1,nVal
    sArray(i) = repeat(" ",mLen)
  end do
  
  iVal = 1
  iSt  = 0
  do i=1,len(mskSplit)
    if(mskSplit(i:i) == "X") then
      iSt = iSt + 1
      if(iSt > mLen) then
        write(lout,"(2(a,i0))") "SPLIT> ERROR Split element ",iVal," is longer than the buffer of ",mLen,"."
        hasErr = .true.
        exit
      end if
      sArray(iVal)(iSt:iSt) = toSplit(i:i)
    else
      if(iSt > 0) then
        iSt  = 0
        iVal = iVal + 1
        if(iVal > nVal) exit
      end if
    end if
  end do
  
end subroutine chr_split

! ================================================================================================ !
!  Scan a String and Mark it for Splitting
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-06-01
! ================================================================================================ !
subroutine chr_scanString(theString, theMask, maxLen, nValues, nIndent, hasErr)
  
  use crcoall
  
  character(len=*),              intent(in)  :: theString
  character(len=:), allocatable, intent(out) :: theMask
  integer,                       intent(out) :: maxLen
  integer,                       intent(out) :: nValues
  integer,                       intent(out) :: nIndent
  logical,                       intent(out) :: hasErr
  
  character ch, pCh
  integer   i, nCh, qSt, nIn, vSt, vLn, nVl
  logical   isI
  
  nCh = len(theString)
  allocate(character(len=nCh) :: theMask)
  theMask = repeat(" ",nCh)
  hasErr  = .false.
  maxLen  = 0
  
  nIn = 0         ! Counter for indents
  vSt = 0         ! 0 = no value, 1 = in value, 2 = comment char
  isI = .true.    ! If we are in the beginning of the line (no values yet)
  qSt = 0         ! 0 = not in quote, 1 = in single quote, 2 = in double quote
  
  do i=1,nCh
    ch  = theString(i:i)
    vSt = 0                               ! Default to treat everything as not a value
    if(ch == " " .and. isI) nIn = nIn + 1 ! Count indents, but only spaces
    if(ch /= " ") isI = .false.           ! Stop counting indents
    if(ch == char(0)) ch = " "            ! Treat null as space
    if(ch == char(9)) ch = " "            ! Treat tab as space
    if(ch == "'" .and. qSt == 0) qSt = 1  ! Entering single quoted region
    if(ch == '"' .and. qSt == 0) qSt = 3  ! Entering double quoted region
    if(ch /= " " .and. qSt == 0) vSt = 1  ! This is a value if it is not in quotes
    if(ch /= "'" .and. qSt == 2) vSt = 1  ! This is a value in single quotes
    if(ch /= '"' .and. qSt == 4) vSt = 1  ! This is a value in double quotes
    if(ch == "'" .and. qSt == 2) qSt = 0  ! Exiting single quoted region
    if(ch == '"' .and. qSt == 4) qSt = 0  ! Exiting double quoted region
    if(ch == "!" .and. qSt == 0) vSt = 2  ! Comment character encountered
    if(qSt == 1) qSt = 2                  ! Flag the newly entered quoted region for saving values next time
    if(qSt == 3) qSt = 4                  ! Flag the newly entered quoted region for saving values next time
    if(vSt == 1) theMask(i:i) = "X"       ! Mark character as a value
    if(vSt == 2) exit                     ! We've reached a comment character, exit
    if(ichar(ch) < 32) then               ! This is a control character, we don't want those
      write(lout,"(2(a,i0))") "SPLIT> ERROR Control character char(",ichar(ch),") encountered at position ",i
      hasErr = .true.
      return
    end if
  end do
  nIndent = nIn
  
  ! Report un-closed quotes
  if(qSt > 0) then
    write(lout,"(a,i0,a)") "SPLIT> ERROR Reached end of line with quotes still open."
    hasErr = .true.
    return
  end if
  
  ! Sum everything up
  vLn = 0
  nVl = 0
  pCh = " "
  do i=1,nCh
    ch = theMask(i:i)
    if(ch == "X") then                   ! This index is part of a value
      vLn    = vLn + 1                   ! Increment the value length
      maxLen = max(maxLen, vLn)          ! Update max value length
      if(pCh == " ") nVl = nVl + 1       ! If previous char was space, count the "edge" as a new value
    else                                 ! This index is not part of a value
      vLn = 0                            ! Reset the value length counter
    end if
    pCh = ch                             ! Record the current char for next round
  end do
  nValues = nVl
  
end subroutine chr_scanString

! ================================================================================================ !
!  Expand Brackets
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-18
!  Will return a string with the content of brackets N(something something) repeated N times.
!  This routine is NOT recursive.
! ================================================================================================ !
function chr_expandBrackets(theString) result(theResult)
  
  implicit none
  
  character(len=*), intent(in)  :: theString
  character(len=:), allocatable :: theResult, theBuffer
  
  integer, allocatable :: bPos(:,:)
  integer ch, nLB, nRB, nB, tSP, lSP, lLB, iSet, nSet, iPos, iMult
  logical iErr
  
  ! Count the brackets
  nLB = 0
  nRB = 0
  do ch=1,len(theString)
    if(theString(ch:ch) == "(") nLB = nLB + 1
    if(theString(ch:ch) == ")") nRB = nRB + 1
  end do
  
  ! If there are none, then just return
  if(nLB == 0 .or. nLB /= nRB) then
    theResult = theString
    return
  end if
  
  ! Otherwise, get all the positions for slicing
  allocate(bPos(3,nLB))
  theBuffer = " "//theString//" "
  bPos(:,:) = 0
  lSP  = 0
  lLB  = 0
  iSet = 0
  do ch=1,len(theBuffer)
    if(theBuffer(ch:ch) == " ") tSP = ch
    if(theBuffer(ch:ch) == "(") then
      lSP = tSP
      lLB = ch
    end if
    if(theBuffer(ch:ch) == ")") then
      if(lSP > 0 .and. lLB > 0 .and. lLB-lSP > 1 .and. ch-lSP > 2) then
        iSet = iSet + 1
        bPos(1,iSet) = lSP
        bPos(2,iSet) = lLB
        bPos(3,iSet) = ch
        lSP = 0
        lLB = 0
      end if
    end if
  end do
  nSet = iSet
  
  ! Then combine all the pieces
  iPos      = 1
  theResult = ""
  do iSet=1,nSet
    theResult = theResult//theBuffer(iPos:bPos(1,iSet))
    call chr_toInt(theBuffer(bPos(1,iSet)+1:bPos(2,iSet)-1),iMult,iErr)
    theResult = theResult//repeat(theBuffer(bPos(2,iSet)+1:bPos(3,iSet)-1)//" ",iMult)
    iPos      = bPos(3,iSet)+1
  end do
  theResult = trim(adjustl(theResult//theBuffer(iPos:)))
  
  deallocate(bPos)
  
end function chr_expandBrackets

! ================================================================================================ !
!  SubString Routine
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-14
!  Returns a substring of a string
!  This is a safer way of extracting a substring than slicing string%chr.
!  Writing a substring of a string to another string produces garbage with the intel compiler.
! ================================================================================================ !
function str_sub(theString, iA, iB) result(retString)
  
  type(string), intent(in)      :: theString
  type(string)                  :: retString
  character(len=:), allocatable :: tmpString
  
  integer iA, iB, oldLen, newLen
  
  oldLen = len(theString%chr)
  newLen = iB - iA + 1
  
  if(newLen <= 0) then
    retString = string("")
    return
  end if
  if(iA < 1)      iA = 1
  if(iB > oldLen) iB = oldLen
  
  allocate(character(newLen) :: tmpString)
  tmpString(1:newLen) = theString%chr(iA:iB)
  retString           = string(tmpString)
  deallocate(tmpString)
  
end function str_sub

! ================================================================================================ !
!  Strip String Routine
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-14
!  Trims leading and trailing white spaces from a string
!  str_trim is for strings in and out
!  chr_trim is for char arrays in and out
! ================================================================================================ !
function str_strip(theString) result(retString)
  type(string), intent(in) :: theString
  type(string)             :: retString
  retString = trim(adjustl(theString))
end function str_strip

function chr_strip(theString) result(retString)
  character(len=*), intent(in)  :: theString
  character(len=:), allocatable :: retString
  retString = trim(adjustl(theString))
end function chr_strip

! ================================================================================================ !
!  Pad String with Zeros or Spaces
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-18
! ================================================================================================ !
function chr_padZero(theString, theSize) result(retString)
  character(len=*), intent(in)  :: theString
  integer,          intent(in)  :: theSize
  character(len=:), allocatable :: retString
  integer                       :: inSize
  inSize = len(theString)
  if(inSize > 0 .and. inSize < theSize) then
    retString = theString(1:inSize)//repeat(char(0),theSize-inSize)
  else
    retString = theString
  end if
end function chr_padZero

function chr_padSpace(theString, theSize) result(retString)
  character(len=*), intent(in)  :: theString
  integer,          intent(in)  :: theSize
  character(len=:), allocatable :: retString
  integer                       :: inSize
  inSize = len(theString)
  if(inSize > 0 .and. inSize < theSize) then
    retString = theString(1:inSize)//repeat(" ",theSize-inSize)
  else
    retString = theString
  end if
end function chr_padSpace

! ================================================================================================ !
!  Trim Zero String Routine
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-14
!  Removes zero chars from any length string
! ================================================================================================ !
function chr_trimZero(theString) result(retString)
  
  character(len=*), intent(in)  :: theString
  character(len=:), allocatable :: retString
  
  integer ch
  
  allocate(character(len=len(theString)) :: retString)
  do ch=1, len(theString)
    if(theString(ch:ch) == char(0)) then
      retString(ch:ch) = " "
    else
      retString(ch:ch) = theString(ch:ch)
    end if
  end do
  retString = trim(retString)
  
end function chr_trimZero

! ================================================================================================ !
!  Count the Occurence of a Single Character in a String
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-16
! ================================================================================================ !
function str_inStr(theString, theNeedle) result(nNeedle)
  
  type(string),     intent(in)  :: theString
  character(len=1), intent(in)  :: theNeedle
  integer                       :: nNeedle
  
  integer ch
  
  nNeedle = 0
  do ch=1, len(theString%chr)
    if(theString%chr(ch:ch) == theNeedle) then
      nNeedle = nNeedle + 1
    end if
  end do
  
end function str_inStr

function chr_inStr(theString, theNeedle) result(nNeedle)
  
  character(len=*), intent(in)  :: theString
  character(len=1), intent(in)  :: theNeedle
  integer                       :: nNeedle
  
  integer ch
  
  nNeedle = 0
  do ch=1, len(theString)
    if(theString(ch:ch) == theNeedle) then
      nNeedle = nNeedle + 1
    end if
  end do
  
end function chr_inStr

! ================================================================================================ !
!  Strip Quotes Routine
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-14
!  Removes a matching pair of double quotes at the beginning and end of a string
! ================================================================================================ !
function str_stripQuotes(theString) result(retString)
  
  type(string), intent(in) :: theString
  type(string)             :: tmpString
  type(string)             :: retString
  
  integer strLen
  
  tmpString = str_strip(theString)
  strLen    = len(tmpString)
  
  if(strLen < 2) then
    retString = tmpString
    return
  end if
  
  if(tmpString%chr(1:1) == '"' .and. tmpString%chr(strLen:strLen) == '"') then
    retString = str_sub(tmpString,2,strLen-1)
    return
  else
    retString = tmpString
    return
  end if
  
end function str_stripQuotes

function chr_stripQuotes(theString) result(retString)
  
  character(len=*), intent(in)  :: theString
  character(len=:), allocatable :: tmpString
  character(len=:), allocatable :: retString
  
  integer strLen
  
  tmpString = chr_strip(theString)
  strLen    = len(tmpString)
  
  if(strLen < 2) then
    retString = tmpString
    return
  end if
  
  if(tmpString(1:1) == '"' .and. tmpString(strLen:strLen) == '"') then
    retString = tmpString(2:strLen-1)
    return
  else
    retString = tmpString
    return
  end if
  
end function chr_stripQuotes

! ================================================================================================ !
!  Rounding Routines
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-20
!  A wrapper for round_near for strings and character arrays
! ================================================================================================ !
subroutine str_toReal(theString, theValue, rErr)
  
  use floatPrecision
  use crcoall
  
  implicit none
  
  type(string), intent(in)        :: theString
  real(kind=fPrec), intent(out)   :: theValue
  logical,          intent(inout) :: rErr
  
  character(len=:), allocatable   :: tmpString
  integer                         :: readErr
  
#if !defined(FIO) && defined(CRLIBM)
  real(kind=fPrec)                :: round_near
  character(len=:), allocatable   :: cString
  integer                         :: cLen, cErr
  
  tmpString = chr_trimZero(theString%chr)
  cLen      = len(tmpString) + 1
  cString   = tmpString//char(0)
  theValue  = round_near(cErr,cLen,cString)
  
  if(cErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input with CRLIBM"
    write (lout,"(a)")    "TYPECAST> Overfow/Underflow in round_near"
    write (lout,"(a,i2)") "TYPECAST> Error value: ",cErr
    rErr = .true.
  end if
#endif
  
#if defined(FIO) && defined(CRLIBM)
  call enable_xp()
  tmpString = chr_trimZero(theString%chr)
  read(tmpString,*,round="nearest",iostat=readErr) theValue
  call disable_xp()
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input with FIO overriding CRLIBM"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to real"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
#endif
  
#if defined(FIO) && !defined(CRLIBM)
  tmpString = chr_trimZero(theString%chr)
  read(tmpString,*,round="nearest",iostat=readErr) theValue
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input with FIO"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to real"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
#endif
  
#if !defined(FIO) && !defined(CRLIBM)
  tmpString = chr_trimZero(theString%chr)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to real"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
#endif
  
end subroutine str_toReal

subroutine chr_toReal(theString, theValue, rErr)
  
  use floatPrecision
  use crcoall
  
  implicit none
  
  character(len=*), intent(in)    :: theString
  real(kind=fPrec), intent(out)   :: theValue
  logical,          intent(inout) :: rErr
  
  character(len=:), allocatable   :: tmpString
  integer                         :: readErr
  
#if !defined(FIO) && defined(CRLIBM)
  real(kind=fPrec)                :: round_near
  character(len=:), allocatable   :: cString
  integer                         :: cLen, cErr
  
  tmpString = chr_trimZero(theString)
  cLen      = len(tmpString) + 1
  cString   = tmpString//char(0)
  theValue  = round_near(cErr,cLen,cString)
  
  if(cErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input with CRLIBM"
    write (lout,"(a)")    "TYPECAST> Overfow/Underflow in round_near"
    write (lout,"(a,i2)") "TYPECAST> Error value: ",cErr
    rErr = .true.
  end if
#endif
  
#if defined(FIO) && defined(CRLIBM)
  call enable_xp()
  tmpString = chr_trimZero(theString)
  read(tmpString,*,round="nearest",iostat=readErr) theValue
  call disable_xp()
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input with FIO overriding CRLIBM"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to real"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
#endif
  
#if defined(FIO) && !defined(CRLIBM)
  tmpString = chr_trimZero(theString)
  read(tmpString,*,round="nearest",iostat=readErr) theValue
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input with FIO"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to real"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
#endif
  
#if !defined(FIO) && !defined(CRLIBM)
  tmpString = chr_trimZero(theString)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to real"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
#endif
  
end subroutine chr_toReal

! ================================================================================================ !
!  String to Integer
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-18
! ================================================================================================ !

subroutine str_toInt(theString, theValue, rErr)
  
  use crcoall
  
  type(string),     intent(in)    :: theString
  integer,          intent(out)   :: theValue
  logical,          intent(inout) :: rErr
  
  character(len=:), allocatable   :: tmpString
  integer                         :: readErr
  
  tmpString = chr_trimZero(theString%chr)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to integer"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
  
end subroutine str_toInt

subroutine chr_toInt(theString, theValue, rErr)
  
  use crcoall
  
  character(len=*), intent(in)    :: theString
  integer,          intent(out)   :: theValue
  logical,          intent(inout) :: rErr
  
  character(len=:), allocatable   :: tmpString
  integer                         :: readErr
  
  tmpString = chr_trimZero(theString)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to integer"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
  
end subroutine chr_toInt

! ================================================================================================ !
!  String to Logical
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-18
! ================================================================================================ !

subroutine str_toLog(theString, theValue, rErr)
  
  use crcoall
  
  type(string),     intent(in)    :: theString
  logical,          intent(out)   :: theValue
  logical,          intent(inout) :: rErr
  
  character(len=:), allocatable   :: tmpString
  integer                         :: readErr
  
  tmpString = chr_trimZero(theString%chr)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to logical"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
  
end subroutine str_toLog

subroutine chr_toLog(theString, theValue, rErr)
  
  use crcoall
  
  character(len=*), intent(in)    :: theString
  logical,          intent(out)   :: theValue
  logical,          intent(inout) :: rErr
  
  character(len=:), allocatable   :: tmpString
  integer                         :: readErr
  
  tmpString = chr_trimZero(theString)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write (lout,"(a)")    "TYPECAST> ERROR Data Input"
    write (lout,"(a)")    "TYPECAST> Failed to cast '"//tmpString//"' to logical"
    write (lout,"(a,i2)") "TYPECAST> iostat value: ",readErr
    rErr = .true.
  end if
  
end subroutine chr_toLog

! ================================================================================================ !
!  HERE FOLLOWS THE OLD ROUTINES
! ================================================================================================ !

! ================================================================================================ !
!  A.Mereghetti, for the FLUKA Team
!  K.Sjobak and A.Santamaria, BE-ABP-HSS
!  Last modified: 2018-05-14
!
!  Parse a line and split it into its fields fields are returned as 0-terminated and padded string.
!
!  input:
!    inLine: usually line read in from fort.2 or fort.3. Values must be separated by spaces
!  Output:
!    Array of values with 
!      gFields(i): (char) value of field
!      lFields(i): (int) length of field
!      nFields:    (int) number of fields
!      errFields:  (logical)
! ================================================================================================ !
subroutine getfields_split(inLine, gFields, lFields, nFields, errFields)
  
  use crcoall
  
  implicit none
  
  character inLine*(str_maxLen-1)               ! nchars in daten is 160
  character gFields(str_maxFields)*(str_maxLen) ! Array of fields
  integer   nFields                             ! Number of identified fields
  integer   lFields(str_maxFields)              ! Length of each what:
  logical   errFields                           ! An error flag
  
  intent(in)  inLine
  intent(out) gFields, lFields, nFields, errFields
  
  ! Runtime variables
  integer ii, jj
  logical lchar
  integer lenstr, istart
  
  ! Initialise output variables
  errFields = .false.
  nFields   = 0
  lenstr    = 0
  istart    = 0
  
  do ii=1,str_maxFields
    do jj=1,str_maxLen
      gFields(ii)(jj:jj) = char(0) ! ZERO terminate/pad
    end do
    lFields(ii) = 0
  end do
  
  ! Parse the line
  lchar = .false.
  do ii=1, str_maxLen-1 ! For \0 termination
    if(inLine(ii:ii) == " ") then
      ! Blank char
      if(lchar) then
        ! End of a string: record it
        lFields(nFields) = lenstr
        gFields(nFields)(1:lenstr) = inLine(istart:istart+lenstr)
        lchar = .false.
      end if
    else
      ! Non-blank char
      if(.not.lchar) then
        ! A new what starts
        nFields = nFields + 1
        if(nFields > str_maxFields) then
          write (lout,"(a)") "ERROR: Too many fields in line:"
          write (lout,"(a)") inLine
          write (lout,"(a)") "please increase str_maxFields"
          errFields = .true.
          exit ! Break do
        end if
        istart = ii
        lchar  = .true.
        lenstr = 0
      end if
      lenstr = lenstr + 1
    end if
  end do
  
end subroutine getfields_split

! ================================================================================================ !
!  K. Sjobak, BE-ABP/HSS
!  Last modified: 30-10-2014
!  Replace "\0" with ' ' in strings.
!  Usefull before output, else "write (*,*)" will actually write all the \0s
!  Warning: Do not add any write(*,*) inside this function:
!  if this function is called by a write(*,*) and then does a write, the program may deadlock!
! ================================================================================================ !
function stringzerotrim(instring) result(retval)
  
  implicit none
  
  character(len=stringzerotrim_maxlen), intent(in) :: instring
  character(len=stringzerotrim_maxlen) :: retval
  
  integer ii
  
  do ii=1,stringzerotrim_maxlen
    if (instring(ii:ii) /= char(0)) then
      retval(ii:ii) = instring(ii:ii)
    else 
      retval(ii:ii) = " "
    end if
  end do
  retval = trim(retval)
  
end function stringzerotrim

end module string_tools
