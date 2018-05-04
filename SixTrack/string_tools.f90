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
  public str_inStr, chr_inStr
  public str_toReal, chr_toReal
  
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
!  Last modified: 2018-04-14
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
  logical sngQuote, dblQuote
  
  newBit   = 0
  nArray   = 0
  sngQuote = .false.
  dblQuote = .false.
  do ch=1, len(toSplit%chr)
    if(toSplit%chr(ch:ch) == "'") sngQuote = .not. sngQuote
    if(toSplit%chr(ch:ch) == '"') dblQuote = .not. dblQuote
    if(toSplit%chr(ch:ch) == " " .and. .not. sngQuote .and. .not. dblQuote) then
      if(newBit == 0) then
        cycle
      else
        call str_arrAppend(returnArray, str_sub(toSplit, newBit, ch-1))
        ! call str_arrAppend(returnArray, string(toSplit%chr(newBit:ch-1)))
        newBit = 0
        nArray = nArray + 1
      end if
    else
      if(newBit == 0) newBit = ch
    end if
  end do
  
end subroutine str_split

! ================================================================================================ !
!  Safe Append to Array
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-04-14
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
      write(lout,"(a)") "ERROR Allocation of string array failed."
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
      write(lout,"(a)") "ERROR Allocation of string array failed."
      stop 1
    end if
    theArray(1) = theString
    
  end if
  
end subroutine str_arrAppend

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
!  Rounting Routines
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
  
  cLen     = len(theString%chr) + 1
  cString  = theString%chr//char(0)
  theValue = round_near(cErr,cLen,cString)
  
  if(cErr /= 0) then
    write (lout,"(a)")    "ERROR Data Input Error"
    write (lout,"(a)")    "Overfow/Underflow in string_tools->str_toReal"
    write (lout,"(a,i2)") "Errno: ",cErr
    stop 1
  end if
#else
  read(theString%chr,*) theValue
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
  read(theString,*) theValue
#endif
  
end function chr_toReal

! ================================================================================================ !
!  HERE FOLLOWS THE OLD ROUTINES
! ================================================================================================ !

! ================================================================================================ !
!  A.Mereghetti, for the FLUKA Team
!  K.Sjobak and A.Santamaria, BE-ABP-HSS
!  Last modified: 2018-04-13
!
!  Parse a line and split it into its fields fields are returned as 0-terminated and padded string.
!
!  input:
!    tmpline: usually line read in from fort.2 or fort.3. Values must be separated by spaces
!  Output:
!    Array of values with 
!      getfields_fields(i):  (char) value of field
!      getfields_lfields(i): (int) length of field
!      getfields_nfields:    (int) number of fields
!      getfields_lerr:       (logical)
! ================================================================================================ !
subroutine getfields_split(tmpline, getfields_fields, getfields_lfields, getfields_nfields, getfields_lerr)
  
  use crcoall
  
  implicit none
  
  character tmpline*(getfields_l_max_string-1)                                ! nchars in daten is 160
  character getfields_fields(getfields_n_max_fields)*(getfields_l_max_string) ! Array of fields
  integer   getfields_nfields                                                 ! Number of identified fields
  integer   getfields_lfields(getfields_n_max_fields)                         ! Length of each what:
  logical   getfields_lerr                                                    ! An error flag
  
  intent(in) tmpline
  intent(out) getfields_fields, getfields_lfields, getfields_nfields, getfields_lerr
  
  ! Runtime variables
  integer ii, jj
  logical lchar
  integer lenstr, istart
  
  ! Initialise output variables
  getfields_lerr    = .false.
  getfields_nfields = 0
  
  do ii=1,getfields_n_max_fields
    do jj=1,getfields_l_max_string
      getfields_fields(ii)(jj:jj) = char(0) ! ZERO terminate/pad
    end do
    getfields_lfields(ii) = 0
  end do
  
  ! Parse the line
  lchar = .false.
  do ii=1, getfields_l_max_string-1 ! For \0 termination
    if(tmpline(ii:ii) == " ") then
      ! Blank char
      if(lchar) then
        ! End of a string: record it
        getfields_lfields(getfields_nfields) = lenstr
        getfields_fields(getfields_nfields)(1:lenstr) = tmpline(istart:istart+lenstr)
        lchar = .false.
      end if
    else
      ! Non-blank char
      if(.not.lchar) then
        ! A new what starts
        getfields_nfields = getfields_nfields +1
        if(getfields_nfields > getfields_n_max_fields) then
          write (lout,*) "ERROR: Too many fields in line:"
          write (lout,*) tmpline
          write (lout,*) "please increase getfields_n_max_fields"
          getfields_lerr = .true.
          exit ! Break do
        end if
        istart = ii
        lchar  = .true.
        lenstr = 0
      end if
      lenstr = lenstr+1
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
