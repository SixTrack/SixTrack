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
  
  ! "Standard" string length, +1 for \0
  integer, parameter :: str_maxLen    = 161
  integer, parameter :: str_maxFields = 15
  
  ! Dummy empty strings
  character(len=str_maxLen), parameter :: str_dSpace = repeat(" ",str_maxLen)
  character(len=str_maxLen), parameter :: str_dZeros = repeat(char(0),str_maxLen)
  
  public str_strip, chr_strip, chr_trimZero
  public str_stripQuotes, chr_stripQuotes
  public str_sub
  public chr_padZero
  public str_inStr, chr_inStr
  public str_toReal, chr_toReal
  public str_toInt, chr_toInt
  
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
subroutine str_split(toSplit, returnArray, nArray)
  
  implicit none
  
  type(string),              intent(in)  :: toSplit
  type(string), allocatable, intent(out) :: returnArray(:)
  integer,                   intent(out) :: nArray
  
  integer ch, newBit
  logical sngQ, dblQ
  
  if(allocated(returnArray)) deallocate(returnArray)
  
  newBit = 0
  nArray = 0
  sngQ   = .false.
  dblQ   = .false.
  do ch=1, len(toSplit%chr)
    if(toSplit%chr(ch:ch) == "'") sngQ = .not. sngQ
    if(toSplit%chr(ch:ch) == '"') dblQ = .not. dblQ
    if((toSplit%chr(ch:ch) == " " .or. toSplit%chr(ch:ch) == char(9)) .and. .not. sngQ .and. .not. dblQ) then
      if(newBit == 0) then
        cycle
      else
        call str_arrAppend(returnArray, str_sub(toSplit, newBit, ch-1))
        newBit = 0
        nArray = nArray + 1
      end if
    else
      if(newBit == 0) newBit = ch
    end if
  end do
  
end subroutine str_split

subroutine chr_split(toSplit, returnArray, nArray)
  
  implicit none
  
  character(len=*),              intent(in)  :: toSplit
  character(len=:), allocatable, intent(out) :: returnArray(:)
  integer,                       intent(out) :: nArray
  
  integer ch, newBit
  logical sngQ, dblQ
  
  if(allocated(returnArray)) deallocate(returnArray)
  
  newBit = 0
  nArray = 0
  sngQ   = .false.
  dblQ   = .false.
  do ch=1, len(toSplit)
    if(toSplit(ch:ch) == "'") sngQ = .not. sngQ
    if(toSplit(ch:ch) == '"') dblQ = .not. dblQ
    if((toSplit(ch:ch) == " " .or. toSplit(ch:ch) == char(9)) .and. .not. sngQ .and. .not. dblQ) then
      if(newBit == 0) then
        cycle
      else
        call chr_arrAppend(returnArray, toSplit(newBit:ch-1))
        newBit = 0
        nArray = nArray + 1
      end if
    else
      if(newBit == 0) newBit = ch
    end if
  end do
  
end subroutine chr_split

! ================================================================================================ !
!  Safe Append to Array
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-18
!  Appends a string to a string array.
! ================================================================================================ !
subroutine str_arrAppend(theArray, theString)
  
  use crcoall
  
  implicit none
  
  type(string), allocatable, intent(inout) :: theArray(:)
  type(string),              intent(in)    :: theString
  
  type(string), allocatable :: tmpArray(:)
  integer                   :: allocErr
  integer                   :: arrSize, arrElem
  
  if(allocated(theArray)) then
    
    arrSize = size(theArray,1)
    allocate(tmpArray(arrSize + 1), stat=allocErr)
    if(allocErr /= 0) then
      write(lout,"(a)") "STRING_TOOLS> ERROR Appending of string array failed."
      stop 1
    end if
    
    do arrElem=1, arrSize
      tmpArray(arrElem) = theArray(arrElem)
    end do
    tmpArray(arrSize + 1) = theString
    
    call move_alloc(tmpArray,theArray)
    
  else
    
    allocate(theArray(1), stat=allocErr)
    if(allocErr /= 0) then
      write(lout,"(a)") "STRING_TOOLS> ERROR Allocation of string array failed."
      stop 1
    end if
    theArray(1) = theString
    
  end if
  
end subroutine str_arrAppend

subroutine chr_arrAppend(theArray, theString)
  
  use crcoall
  
  implicit none
  
  character(len=:), allocatable, intent(inout) :: theArray(:)
  character(len=*),              intent(in)    :: theString
  
  character(len=:), allocatable :: tmpArray(:)
  integer :: allocErr
  integer :: arrSize, arrElem
  integer :: inLen, maxLen, elemLen
  
  inLen = len(theString)
  
  if(allocated(theArray)) then
    
    maxLen = inLen
    do arrElem=1, len(theArray)
      elemLen = len(theArray(arrElem))
      if(elemLen > maxLen) maxLen = elemLen
    end do
    
    arrSize = size(theArray,1)
    allocate(character(len=maxLen) :: tmpArray(arrSize + 1), stat=allocErr)
    if(allocErr /= 0) then
      write(lout,"(a)") "STRING_TOOLS> ERROR Appending of string array failed."
      stop 1
    end if
    
    do arrElem=1, arrSize
      elemLen = len(theArray(arrElem))
      tmpArray(arrElem)            = repeat(char(0),maxLen)
      tmpArray(arrElem)(1:elemLen) = theArray(arrElem)(1:elemLen)
    end do
    tmpArray(arrSize + 1)          = repeat(char(0),maxLen)
    tmpArray(arrSize + 1)(1:inLen) = theString(1:inLen)
    
    call move_alloc(tmpArray,theArray)
    
  else
    
    allocate(character(len=inLen) :: theArray(1), stat=allocErr)
    if(allocErr /= 0) then
      write(lout,"(a)") "STRING_TOOLS> ERROR Allocation of character array failed."
      stop 1
    end if
    theArray(1) = theString
    
  end if
  
end subroutine chr_arrAppend

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
!  Pad String with Zeros
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

! ================================================================================================ !
!  Trim Zero String Routine
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-14
!  Cuts the string at first char(0)
! ================================================================================================ !
function chr_trimZero(theString) result(retString)
  
  character(len=*), intent(in)  :: theString
  character(len=:), allocatable :: retString
  
  integer ch, cut
  do ch=1, len(theString)
    if(theString(ch:ch) == char(0)) then
      cut = ch-1
      exit
    end if
  end do
  
  if(cut > 0 .and. cut < len(theString)) then
    retString = theString(1:cut)
  else
    retString = theString
  end if
  
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
!  Removes a matching pair of single or double quotes at the beginning or end of a string
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
  elseif(tmpString%chr(1:1) == "'" .and. tmpString%chr(strLen:strLen) == "'") then
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
  elseif(tmpString(1:1) == "'" .and. tmpString(strLen:strLen) == "'") then
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
function str_toReal(theString) result(theValue)
  
  use floatPrecision
#ifdef CRLIBM
  use crcoall
#endif
  
  implicit none
  
  type(string), intent(in) :: theString
  real(kind=fPrec)         :: theValue
  
#ifdef CRLIBM
  real(kind=fPrec)              :: round_near
  character(len=:), allocatable :: cString
  integer                       :: cLen, cErr
  
  cLen     = len(theString) + 1
  cString  = theString%chr//char(0)
  theValue = round_near(cErr,cLen,cString)
  
  if(cErr /= 0) then
    write (lout,"(a)")    "ERROR Data Input Error"
    write (lout,"(a)")    "Overfow/Underflow in string_tools->str_toReal"
    write (lout,"(a,i2)") "Errno: ",cErr
    stop 1
  end if
#else
#ifdef FIO
#ifdef CRLIBM
  call enable_xp()
#endif
  read(theString%chr,*,round="nearest") theValue
#ifdef CRLIBM
  call disable_xp()
#endif
#else
  read(theString%chr,*) theValue
#endif
#endif
  
end function str_toReal

function chr_toReal(theString) result(theValue)
  
  use floatPrecision
#ifdef CRLIBM
  use crcoall
#endif
  
  implicit none
  
  character(len=*), intent(in) :: theString
  real(kind=fPrec)             :: theValue
  
#ifdef CRLIBM
  real(kind=fPrec)              :: round_near
  character(len=:), allocatable :: cString
  integer                       :: cLen, cErr
  
  cLen     = len(theString) + 1
  cString  = theString//char(0)
  theValue = round_near(cErr,cLen,cString)
  
  if(cErr /= 0) then
    write (lout,"(a)")    "++++++++++++++++++++++++"
    write (lout,"(a)")    "+    ERROR DETECTED    +"
    write (lout,"(a)")    "++++++++++++++++++++++++"
    write (lout,"(a)")    "Data Input Error"
    write (lout,"(a)")    "Overfow/Underflow in string_tools->chr_toReal"
    write (lout,"(a,i2)") "Errno: ",cErr
    stop -1
  end if
#else
#ifdef FIO
#ifdef CRLIBM
  call enable_xp()
#endif
  read(theString,*,round="nearest") theValue
#ifdef CRLIBM
  call disable_xp()
#endif
#else
  read(theString,*) theValue
#endif
#endif
  
end function chr_toReal

! ================================================================================================ !
!  String to Integer
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-18
! ================================================================================================ !

function str_toInt(theString) result(theValue)
  type(string), intent(in) :: theString
  integer                  :: theValue
  read(theString%chr,*) theValue
end function str_toInt

function chr_toInt(theString) result(theValue)
  character(len=*), intent(in) :: theString
  integer                      :: theValue
  read(theString,*) theValue
end function chr_toInt

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
