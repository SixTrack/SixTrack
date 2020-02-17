! ================================================================================================ !
!  SixTrack String Tools
! ~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-16
!
!  Note: Careful with adding use <module> to this module as it easily creates circular dependencies
! ================================================================================================ !
module string_tools

  use strings
  use parpro, only : mNameLen, mFileName, mPathName
  use, intrinsic :: iso_fortran_env, only : int16, int32, int64, real32, real64, real128

  implicit none

  public str_strip, chr_strip, chr_trimZero
  public str_stripQuotes, chr_stripQuotes
  public str_sub, chr_expandBrackets
  public chr_lPad, chr_rPad, chr_lPadCut, chr_rPadCut
  public chr_toUpper, chr_toLower, chr_isNumeric
  public str_inStr, chr_inStr

  interface str_cast
    module procedure str_toReal32
    module procedure str_toReal64
    module procedure str_toReal128
    module procedure str_toInt16
    module procedure str_toInt32
    module procedure str_toInt64
    module procedure str_toLog
  end interface str_cast

  interface chr_cast
    module procedure chr_toReal32
    module procedure chr_toReal64
    module procedure chr_toReal128
    module procedure chr_toInt16
    module procedure chr_toInt32
    module procedure chr_toInt64
    module procedure chr_toLog
  end interface chr_cast

  private :: str_toReal32
  private :: str_toReal64
  private :: str_toReal128
  private :: str_toInt16
  private :: str_toInt32
  private :: str_toInt64
  private :: str_toLog

  private :: chr_toReal32
  private :: chr_toReal64
  private :: chr_toReal128
  private :: chr_toInt16
  private :: chr_toInt32
  private :: chr_toInt64
  private :: chr_toLog

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
        write(lerr,"(2(a,i0))") "SPLIT> ERROR Split element ",iVal," is longer than the buffer of ",mLen,"."
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

  character(len=*),              intent(in)    :: toSplit
  character(len=:), allocatable, intent(inout) :: sArray(:)
  integer,                       intent(out)   :: nArray
  logical,                       intent(out)   :: hasErr
  integer,          optional,    intent(in)    :: fixLen
  integer,          optional,    intent(out)   :: nIndent

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
        write(lerr,"(2(a,i0))") "SPLIT> ERROR Split element ",iVal," is longer than the buffer of ",mLen,"."
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
      write(lerr,"(2(a,i0))") "SPLIT> ERROR Control character char(",ichar(ch),") encountered at position ",i
      hasErr = .true.
      return
    end if
  end do
  nIndent = nIn

  ! Report un-closed quotes
  if(qSt > 0) then
    write(lerr,"(a,i0,a)") "SPLIT> ERROR Reached end of line with quotes still open."
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
  integer ch, nLB, nRB, tSP, lSP, lLB, iSet, nSet, iPos, iMult
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
  tSP  = 0
  lLB  = 0
  iSet = 0
  do ch=1,len(theBuffer)
    if(theBuffer(ch:ch) == " ") then
      tSP = ch
    elseif(theBuffer(ch:ch) == "(") then
      lSP = tSP
      lLB = ch
    elseif(theBuffer(ch:ch) == ")") then
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
    call chr_cast(theBuffer(bPos(1,iSet)+1:bPos(2,iSet)-1),iMult,iErr)
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
!  Pad String With Spaces
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-27
!   * chr_lPad    :: Pad a string to the left to specified size. Does not trunctae if longer.
!   * chr_rPad    :: Pad a string to the right to specified size. Does not trunctae if longer.
!   * chr_lPadCut :: Pad a string to the left to specified size. Does trunctae if longer.
!   * chr_rPadCut :: Pad a string to the right to specified size. Does trunctae if longer.
! ================================================================================================ !
function chr_lPad(theString, theSize) result(retString)
  character(len=*), intent(in)  :: theString
  integer,          intent(in)  :: theSize
  character(len=:), allocatable :: retString
  integer                       :: inSize
  inSize = len(theString)
  if(inSize > 0 .and. inSize < theSize) then
    retString = repeat(" ",theSize-inSize)//theString(1:inSize)
  else
    retString = theString
  end if
end function chr_lPad

function chr_rPad(theString, theSize) result(retString)
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
end function chr_rPad

function chr_lPadCut(theString, theSize) result(retString)
  character(len=*), intent(in) :: theString
  integer,          intent(in) :: theSize
  character(len=theSize) :: retString
  integer                :: inSize
  inSize = len(theString)
  if(inSize > theSize) then
    retString = theString(1:theSize)
  elseif(inSize > 0) then
    retString = repeat(" ",theSize-inSize)//theString(1:inSize)
  else
    retString = repeat(" ",theSize)
  end if
end function chr_lPadCut

function chr_rPadCut(theString, theSize) result(retString)
  character(len=*), intent(in) :: theString
  integer,          intent(in) :: theSize
  character(len=theSize) :: retString
  integer                :: inSize
  inSize = len(theString)
  if(inSize > theSize) then
    retString = theString(1:theSize)
  elseif(inSize > 0) then
    retString = theString(1:inSize)//repeat(" ",theSize-inSize)
  else
    retString = repeat(" ",theSize)
  end if
end function chr_rPadCut

! ================================================================================================ !
!  Convert to Lower/Upper Case
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-04
! ================================================================================================ !
function chr_toUpper(theString) result(retString)
  character(len=*), intent(in)  :: theString
  character(len=:), allocatable :: retString
  character          :: ch
  integer, parameter :: ulOffset = ichar("A") - ichar("a")
  integer            :: i
  allocate(character(len(theString)) :: retString)
  do i = 1,len(theString)
    ch = theString(i:i)
    if(ch >= "a" .and. ch <= "z") ch = char(ichar(ch)+ulOffset)
    retString(i:i) = ch
  end do
end function chr_toUpper

function chr_toLower(theString) result(retString)
  character(len=*), intent(in)  :: theString
  character(len=:), allocatable :: retString
  character          :: ch
  integer, parameter :: ulOffset = ichar("A") - ichar("a")
  integer            :: i
  allocate(character(len(theString)) :: retString)
  do i = 1,len(theString)
    ch = theString(i:i)
    if(ch >= "A" .and. ch <= "Z") ch = char(ichar(ch)-ulOffset)
    retString(i:i) = ch
  end do
end function chr_toLower

! ================================================================================================ !
!  Check if String is a Number
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-04-05
! ================================================================================================ !
logical function chr_isNumeric(theString)
  character(len=*), intent(in) :: theString
  integer ioStat
  real dumDum
  read(theString,"(e15.6)",iostat=ioStat) dumDum
  if(ioStat == 0) then
    chr_isNumeric = .true.
  else
    chr_isNumeric = .false.
  end if
end function chr_isNumeric

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
subroutine str_toReal32(theString, theValue, rErr)

  use floatPrecision
  use crcoall

  implicit none

  type(string),      intent(in)    :: theString
  real(kind=real32), intent(out)   :: theValue
  logical,           intent(inout) :: rErr

  character(len=:),  allocatable   :: tmpString
  real(kind=real64)                :: tmpValue
  integer                          :: readErr

#if !defined(FIO) && defined(CRLIBM)
  real(kind=fPrec)                 :: round_near
  character(len=:),  allocatable   :: cString
  integer                          :: cLen, cErr

  tmpString = trim(theString%chr)
  cLen      = len(tmpString) + 1
  cString   = tmpString//char(0)
  tmpValue  = round_near(cErr,cLen,cString)

  if(cErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR CRLIBM Failed to cast '"//tmpString//"' to real32 width error ",cErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && defined(CRLIBM)
  call enable_xp()
  tmpString = trim(theString%chr)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  call disable_xp()
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real32 width error ",readErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString%chr)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real32 width error ",readErr
    rErr = .true.
  end if
#endif

#if !defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString%chr)
  read(tmpString,*,iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to real32 width error ",readErr
    rErr = .true.
  end if
#endif

  theValue = real(tmpValue,kind=real32)

end subroutine str_toReal32

subroutine chr_toReal32(theString, theValue, rErr)

  use floatPrecision
  use crcoall

  implicit none

  character(len=*),  intent(in)    :: theString
  real(kind=real32), intent(out)   :: theValue
  logical,           intent(inout) :: rErr

  character(len=:),  allocatable   :: tmpString
  real(kind=real64)                :: tmpValue
  integer                          :: readErr

#if !defined(FIO) && defined(CRLIBM)
  real(kind=fPrec)                 :: round_near
  character(len=:),  allocatable   :: cString
  integer                          :: cLen, cErr

  tmpString = trim(theString)
  cLen      = len(tmpString) + 1
  cString   = tmpString//char(0)
  tmpValue  = round_near(cErr,cLen,cString)

  if(cErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR CRLIBM Failed to cast '"//tmpString//"' to real32 width error ",cErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && defined(CRLIBM)
  call enable_xp()
  tmpString = trim(theString)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  call disable_xp()
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real32 width error ",readErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real32 width error ",readErr
    rErr = .true.
  end if
#endif

#if !defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString)
  read(tmpString,*,iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to real32 width error ",readErr
    rErr = .true.
  end if
#endif

  theValue = real(tmpValue,kind=real32)

end subroutine chr_toReal32

subroutine str_toReal64(theString, theValue, rErr)

  use floatPrecision
  use crcoall

  implicit none

  type(string),      intent(in)    :: theString
  real(kind=real64), intent(out)   :: theValue
  logical,           intent(inout) :: rErr

  character(len=:),  allocatable   :: tmpString
  real(kind=real64)                :: tmpValue
  integer                          :: readErr

#if !defined(FIO) && defined(CRLIBM)
  real(kind=fPrec)                 :: round_near
  character(len=:),  allocatable   :: cString
  integer                          :: cLen, cErr

  tmpString = trim(theString%chr)
  cLen      = len(tmpString) + 1
  cString   = tmpString//char(0)
  tmpValue  = round_near(cErr,cLen,cString)

  if(cErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR CRLIBM Failed to cast '"//tmpString//"' to real64 width error ",cErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && defined(CRLIBM)
  call enable_xp()
  tmpString = trim(theString%chr)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  call disable_xp()
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real64 width error ",readErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString%chr)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real64 width error ",readErr
    rErr = .true.
  end if
#endif

#if !defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString%chr)
  read(tmpString,*,iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to real64 width error ",readErr
    rErr = .true.
  end if
#endif

  theValue = real(tmpValue,kind=real64)

end subroutine str_toReal64

subroutine chr_toReal64(theString, theValue, rErr)

  use floatPrecision
  use crcoall

  implicit none

  character(len=*),  intent(in)    :: theString
  real(kind=real64), intent(out)   :: theValue
  logical,           intent(inout) :: rErr

  character(len=:),  allocatable   :: tmpString
  real(kind=real64)                :: tmpValue
  integer                          :: readErr

#if !defined(FIO) && defined(CRLIBM)
  real(kind=fPrec)                 :: round_near
  character(len=:),  allocatable   :: cString
  integer                          :: cLen, cErr

  tmpString = trim(theString)
  cLen      = len(tmpString) + 1
  cString   = tmpString//char(0)
  tmpValue  = round_near(cErr,cLen,cString)

  if(cErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR CRLIBM Failed to cast '"//tmpString//"' to real64 width error ",cErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && defined(CRLIBM)
  call enable_xp()
  tmpString = trim(theString)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  call disable_xp()
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real64 width error ",readErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real64 width error ",readErr
    rErr = .true.
  end if
#endif

#if !defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString)
  read(tmpString,*,iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to real64 width error = ",readErr
    rErr = .true.
  end if
#endif

  theValue = real(tmpValue,kind=real64)

end subroutine chr_toReal64

subroutine str_toReal128(theString, theValue, rErr)

  use floatPrecision
  use crcoall

  implicit none

  type(string),       intent(in)    :: theString
  real(kind=real128), intent(out)   :: theValue
  logical,            intent(inout) :: rErr

  character(len=:),   allocatable   :: tmpString
  real(kind=real64)                 :: tmpValue
  integer                           :: readErr

#if !defined(FIO) && defined(CRLIBM)
  real(kind=fPrec)                  :: round_near
  character(len=:),   allocatable   :: cString
  integer                           :: cLen, cErr

  tmpString = trim(theString%chr)
  cLen      = len(tmpString) + 1
  cString   = tmpString//char(0)
  tmpValue  = round_near(cErr,cLen,cString)

  if(cErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR CRLIBM Failed to cast '"//tmpString//"' to real128 width error ",cErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && defined(CRLIBM)
  call enable_xp()
  tmpString = trim(theString%chr)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  call disable_xp()
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real128 width error ",readErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString%chr)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real128 width error ",readErr
    rErr = .true.
  end if
#endif

#if !defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString%chr)
  read(tmpString,*,iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to real128 width error ",readErr
    rErr = .true.
  end if
#endif

  theValue = real(tmpValue,kind=real128)

end subroutine str_toReal128

subroutine chr_toReal128(theString, theValue, rErr)

  use floatPrecision
  use crcoall

  implicit none

  character(len=*),   intent(in)    :: theString
  real(kind=real128), intent(out)   :: theValue
  logical,            intent(inout) :: rErr

  character(len=:),   allocatable   :: tmpString
  real(kind=real64)                 :: tmpValue
  integer                           :: readErr

#if !defined(FIO) && defined(CRLIBM)
  real(kind=fPrec)                  :: round_near
  character(len=:),   allocatable   :: cString
  integer                           :: cLen, cErr

  tmpString = trim(theString)
  cLen      = len(tmpString) + 1
  cString   = tmpString//char(0)
  tmpValue  = round_near(cErr,cLen,cString)

  if(cErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR CRLIBM Failed to cast '"//tmpString//"' to real128 width error ",cErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && defined(CRLIBM)
  call enable_xp()
  tmpString = trim(theString)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  call disable_xp()
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real128 width error ",readErr
    rErr = .true.
  end if
#endif

#if defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString)
  read(tmpString,*,round="nearest",iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR FIO Failed to cast '"//tmpString//"' to real128 width error ",readErr
    rErr = .true.
  end if
#endif

#if !defined(FIO) && !defined(CRLIBM)
  tmpString = trim(theString)
  read(tmpString,*,iostat=readErr) tmpValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to real128 width error ",readErr
    rErr = .true.
  end if
#endif

  theValue = real(tmpValue,kind=real128)

end subroutine chr_toReal128

! ================================================================================================ !
!  String to Integer
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-18
! ================================================================================================ !

subroutine str_toInt16(theString, theValue, rErr)

  use crcoall

  type(string),        intent(in)    :: theString
  integer(kind=int16), intent(out)   :: theValue
  logical,             intent(inout) :: rErr

  character(len=:),    allocatable   :: tmpString
  integer                            :: readErr

  tmpString = trim(theString%chr)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to int16 width error ",readErr
    rErr = .true.
  end if

end subroutine str_toInt16

subroutine chr_toInt16(theString, theValue, rErr)

  use crcoall

  character(len=*),    intent(in)    :: theString
  integer(kind=int16), intent(out)   :: theValue
  logical,             intent(inout) :: rErr

  character(len=:),    allocatable   :: tmpString
  integer                            :: readErr

  tmpString = trim(theString)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to int16 width error ",readErr
    rErr = .true.
  end if

end subroutine chr_toInt16

subroutine str_toInt32(theString, theValue, rErr)

  use crcoall

  type(string),        intent(in)    :: theString
  integer(kind=int32), intent(out)   :: theValue
  logical,             intent(inout) :: rErr

  character(len=:),    allocatable   :: tmpString
  integer                            :: readErr

  tmpString = trim(theString%chr)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to int32 width error ",readErr
    rErr = .true.
  end if

end subroutine str_toInt32

subroutine chr_toInt32(theString, theValue, rErr)

  use crcoall

  character(len=*),    intent(in)    :: theString
  integer(kind=int32), intent(out)   :: theValue
  logical,             intent(inout) :: rErr

  character(len=:),    allocatable   :: tmpString
  integer                            :: readErr

  tmpString = trim(theString)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to int32 width error ",readErr
    rErr = .true.
  end if

end subroutine chr_toInt32

subroutine str_toInt64(theString, theValue, rErr)

  use crcoall

  type(string),        intent(in)    :: theString
  integer(kind=int64), intent(out)   :: theValue
  logical,             intent(inout) :: rErr

  character(len=:),    allocatable   :: tmpString
  integer                            :: readErr

  tmpString = trim(theString%chr)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to int64 width error ",readErr
    rErr = .true.
  end if

end subroutine str_toInt64

subroutine chr_toInt64(theString, theValue, rErr)

  use crcoall

  character(len=*),    intent(in)    :: theString
  integer(kind=int64), intent(out)   :: theValue
  logical,             intent(inout) :: rErr

  character(len=:),    allocatable   :: tmpString
  integer                            :: readErr

  tmpString = trim(theString)
  read(tmpString,*,iostat=readErr) theValue
  if(readErr /= 0) then
    write(lerr,"(a,i0)") "TYPECAST> ERROR Failed to cast '"//tmpString//"' to int64 width error ",readErr
    rErr = .true.
  end if

end subroutine chr_toInt64

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

  tmpString = chr_toLower(trim(theString%chr))
  select case(tmpString)
  case("true",".true.","on","1")
    theValue = .true.
  case("false",".false.","off","0")
    theValue = .false.
  case default
    write(lerr,"(a)") "TYPECAST> ERROR Failed to cast '"//trim(theString)//"' to logical"
    rErr = .true.
  end select

end subroutine str_toLog

subroutine chr_toLog(theString, theValue, rErr)

  use crcoall

  character(len=*), intent(in)    :: theString
  logical,          intent(out)   :: theValue
  logical,          intent(inout) :: rErr

  character(len=:), allocatable   :: tmpString

  tmpString = chr_toLower(trim(theString))
  select case(tmpString)
  case("true",".true.","on","1")
    theValue = .true.
  case("false",".false.","off","0")
    theValue = .false.
  case default
    write(lerr,"(a)") "TYPECAST> ERROR Failed to cast '"//trim(theString)//"' to logical"
    rErr = .true.
  end select

end subroutine chr_toLog

! ================================================================================================ !
!  Real to String
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-07-13
! ================================================================================================ !
subroutine chr_fromReal(theValue, theString, nPrec, ePrec, rErr)

  use crcoall
  use floatPrecision

  real(kind=fPrec),             intent(in)    :: theValue
  integer,                      intent(in)    :: nPrec, ePrec
  character(len=nPrec+ePrec+4), intent(out)   :: theString
  logical,                      intent(inout) :: rErr

  character(len=11) dFmt
#ifdef CRLIBM
  character :: dStr(nPrec)
  integer   :: dtoaf, fLen, nStr, dPoint, iSign, i
#endif
  write(dFmt,"(a3,i2.2,a1,i2.2,a1,i1,a1)") "(es",(nPrec+ePrec+4),".",(nPrec-1),"e",ePrec,")"

#ifdef CRLIBM
  fLen = nPrec
  nStr = dtoaf(theValue, 2, nPrec, dPoint, iSign, dStr(1), 1)

  theString = " "
  if(dPoint == 9999) then
    ! This is infinity or nan, so just use default output
    write(theString,dFmt) theValue
  else
    if(iSign /= 0) theString(1:1) = "-"
    theString(2:3) = dStr(1)//"."
    do i=2,nStr ! Get the numbers returned from dtoaf
      theString(i+2:i+2) = dStr(i)
    end do
    do i=nStr+3,nPrec+2 ! Pad the rest with 0
      theString(i:i) = "0"
    end do
    if(dPoint < 1) then
      theString(nPrec+3:nPrec+4) = "E-"
    else
      theString(nPrec+3:nPrec+4) = "E+"
    end if
    select case(ePrec)
    case(2)
      write(theString(nPrec+5:nPrec+6),"(i2.2)") abs(dPoint-1)
    case(3)
      write(theString(nPrec+5:nPrec+7),"(i3.3)") abs(dPoint-1)
    case default
      write(lerr,"(a)") "DTOAF> ERROR Exponent must be either 2 or 3. This is a bug."
      rErr = .true.
      return
    end select
  end if

  ! write(lout,"(a)")           "DTOAF> TESTING:"
  ! write(lout,"(a,i0)")        "DTOAF>  * nStr      = ",nStr
  ! write(lout,"(a,i0)")        "DTOAF>  * dPoint    = ",dPoint
  ! write(lout,"(a,i0)")        "DTOAF>  * iSign     = ",iSign
  ! write(lout,"(a,es24.16e3)") "DTOAF>  * theValue  = ",theValue
  ! write(lout,"(a)")           "DTOAF>  * theString = "//theString

#else
  write(theString,dFmt) theValue
#endif

end subroutine chr_fromReal

! ================================================================================================ !
!  HERE FOLLOWS THE OLD ROUTINES
! ================================================================================================ !

end module string_tools
