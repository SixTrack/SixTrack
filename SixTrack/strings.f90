! ================================================================================================ !
!  SixTrack String Module
! ~~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-25
!
!  Usage
!    Declaration:   type(string) aString
!    Get:           aString%get()                   => Returns chracater array
!    Set:           aString%set("foo")              => Returns string
!    Assignment:    aString = "foo"                 => Returns string
!    Addition:      aString = aString + "bar"       => Returns combined string "foobar"
!    Concatenation: write(*,*) "This is "//aString  => On-the-fly conversion to character array
!    Comparison:    aString == "foobar"             => Compares a string/char to a string/char
!    Comparison:    aString /= "barfoo"             => Compares a string/char to a string/char
! ================================================================================================ !
module strings
  
  implicit none
  
  type, public :: string
    
    character(len=:), allocatable, public :: chr
    
  contains
    
    procedure, public,  pass(this)  :: get     => getStr
    procedure, public,  pass(this)  :: set     => setStr
    procedure, public,  pass(this)  :: len     => lenStr
    procedure, public,  pass(this)  :: trim    => trimStr
    procedure, public,  pass(this)  :: adjustl => adjLStr
    procedure, public,  pass(this)  :: adjustr => adjRStr
    
    procedure, private, pass(this)  :: getStr
    procedure, private, pass(this)  :: setStr
    procedure, private, pass(this)  :: lenStr
    procedure, private, pass(this)  :: trimStr
    procedure, private, pass(this)  :: adjLStr
    procedure, private, pass(this)  :: adjRStr
    
    procedure, private, pass(left)  :: assignStrStr
    procedure, private, pass(left)  :: assignStrChr
    procedure, private, pass(right) :: assignChrStr
    
    procedure, private, pass(left)  :: appendStrStr
    procedure, private, pass(left)  :: appendStrChr
    
    procedure, private, pass(left)  :: concatStrStr
    procedure, private, pass(left)  :: concatStrChr
    procedure, private, pass(right) :: concatChrStr
    
    procedure, private, pass(left)  :: compStrStr
    procedure, private, pass(left)  :: compStrChr
    procedure, private, pass(right) :: compChrStr
    
    procedure, private, pass(left)  :: compNStrStr
    procedure, private, pass(left)  :: compNStrChr
    procedure, private, pass(right) :: compNChrStr
    
  end type string
  
  interface string
    module procedure constructStr
  end interface string
  
  interface len
    module procedure lenStr
  end interface len
  
  interface trim
    module procedure trimStr
  end interface trim
  
  interface adjustl
    module procedure adjLStr
  end interface adjustl
  
  interface adjustr
    module procedure adjRStr
  end interface adjustr
  
  interface assignment(=)
    module procedure assignStrStr
    module procedure assignStrChr
    module procedure assignChrStr
  end interface
  
  interface operator(+)
    module procedure appendStrStr
    module procedure appendStrChr
  end interface
  
  interface operator(//)
    module procedure concatStrStr
    module procedure concatStrChr
    module procedure concatChrStr
  end interface
  
  interface operator(==)
    module procedure compStrStr
    module procedure compStrChr
    module procedure compChrStr
  end interface
  
  interface operator(/=)
    module procedure compNStrStr
    module procedure compNStrChr
    module procedure compNChrStr
  end interface
  
contains
  
  type(string) function constructStr()
    constructStr%chr = ""
  end function constructStr
  
  ! ================================================================ !
  !  Set and Get
  ! ================================================================ !
  pure function getStr(this) result(retValue)
    class(string),    intent(in)  :: this
    character(len=:), allocatable :: retValue
    retValue = this%chr
  end function getStr
  
  pure subroutine setStr(this, setValue)
    class(string),    intent(inout) :: this
    character(len=*), intent(in)    :: setValue
    if(allocated(this%chr)) this%chr = setValue
  end subroutine setStr
  
  ! ================================================================ !
  !  Standard String Information and Operations
  ! ================================================================ !
  elemental function lenStr(this) result(retValue)
    class(string), intent(in) :: this
    integer                   :: retValue
    if(allocated(this%chr)) then
      retValue = len(this%chr)
    else
      retValue = 0
    end if
  end function lenStr
  
  pure function trimStr(this) result(retValue)
    class(string),    intent(in)  :: this
    character(len=:), allocatable :: retValue
    if(allocated(this%chr)) then
      retValue = trim(this%chr)
    end if
  end function trimStr
  
  pure function adjLStr(this) result(retValue)
    class(string),    intent(in)  :: this
    character(len=:), allocatable :: retValue
    if(allocated(this%chr)) then
      retValue = adjustl(this%chr)
    end if
  end function adjLStr
  
  pure function adjRStr(this) result(retValue)
    class(string),    intent(in)  :: this
    character(len=:), allocatable :: retValue
    if(allocated(this%chr)) then
      retValue = adjustr(this%chr)
    end if
  end function adjRStr
  
  ! ================================================================ !
  !  String Assignment
  ! ================================================================ !
  elemental subroutine assignStrStr(left, right)
    class(string), intent(inout) :: left
    class(string), intent(in)    :: right
    left%chr = right%chr
  end subroutine assignStrStr
  
  elemental subroutine assignStrChr(left, right)
    class(string),    intent(inout) :: left
    character(len=*), intent(in)    :: right
    left%chr = right
  end subroutine assignStrChr
  
  elemental subroutine assignChrStr(left, right)
    character(len=*), intent(inout) :: left
    class(string),    intent(in)    :: right
    left = right%chr
  end subroutine assignChrStr
  
  ! ================================================================ !
  !  Append Strings with +, Returning String
  ! ================================================================ !
  type(string) elemental function appendStrStr(left, right)
    class(string),    intent(in)  :: left
    class(string),    intent(in)  :: right
    character(len=:), allocatable :: tmpChr
    integer nOld, nNew
    nOld =        len(left%chr)
    nNew = nOld + len(right%chr)
    allocate(character(nNew) :: tmpChr)
    tmpChr(1:nOld)      = left%chr
    tmpChr(nOld+1:nNew) = right%chr
    appendStrStr%chr    = tmpChr
    deallocate(tmpChr)
  end function appendStrStr
  
  type(string) elemental function appendStrChr(left, right)
    class(string),    intent(in)  :: left
    character(len=*), intent(in)  :: right
    character(len=:), allocatable :: tmpChr
    integer nOld, nNew
    nOld =        len(left%chr)
    nNew = nOld + len(right)
    allocate(character(nNew) :: tmpChr)
    tmpChr(1:nOld)      = left%chr
    tmpChr(nOld+1:nNew) = right
    appendStrChr%chr    = tmpChr
    deallocate(tmpChr)
  end function appendStrChr
  
  ! ================================================================ !
  !  Concat Strings with //, Returning Character
  ! ================================================================ !
  pure function concatStrStr(left, right) result(retValue)
    class(string),    intent(in)  :: left
    class(string),    intent(in)  :: right
    character(len=:), allocatable :: retValue
    integer nOld, nNew
    nOld =        len(left%chr)
    nNew = nOld + len(right%chr)
    allocate(character(nNew) :: retValue)
    retValue(1:nOld)      = left%chr
    retValue(nOld+1:nNew) = right%chr
  end function concatStrStr
  
  pure function concatChrStr(left, right) result(retValue)
    character(len=*), intent(in)  :: left
    class(string),    intent(in)  :: right
    character(len=:), allocatable :: retValue
    integer nOld, nNew
    nOld =        len(left)
    nNew = nOld + len(right%chr)
    allocate(character(nNew) :: retValue)
    retValue(1:nOld)      = left
    retValue(nOld+1:nNew) = right%chr
  end function concatChrStr
  
  pure function concatStrChr(left, right) result(retValue)
    class(string),    intent(in)  :: left
    character(len=*), intent(in)  :: right
    character(len=:), allocatable :: retValue
    integer nOld, nNew
    nOld =        len(left%chr)
    nNew = nOld + len(right)
    allocate(character(nNew) :: retValue)
    retValue(1:nOld)      = left%chr
    retValue(nOld+1:nNew) = right
  end function concatStrChr
  
  ! ================================================================ !
  !  Compare Strings, Equal
  ! ================================================================ !
  elemental function compStrStr(left, right) result(retValue)
    class(string), intent(in) :: left
    class(string), intent(in) :: right
    logical retValue
    retValue = left%chr == right%chr
  end function compStrStr
  
  elemental function compStrChr(left, right) result(retValue)
    class(string),    intent(in) :: left
    character(len=*), intent(in) :: right
    logical retValue
    retValue = left%chr == right
  end function compStrChr
  
  elemental function compChrStr(left, right) result(retValue)
    character(len=*), intent(in) :: left
    class(string),    intent(in) :: right
    logical retValue
    retValue = left == right%chr
  end function compChrStr
  
  ! ================================================================ !
  !  Compare Strings, Not Equal
  ! ================================================================ !
  elemental function compNStrStr(left, right) result(retValue)
    class(string), intent(in) :: left
    class(string), intent(in) :: right
    logical retValue
    retValue = left%chr /= right%chr
  end function compNStrStr
  
  elemental function compNStrChr(left, right) result(retValue)
    class(string),    intent(in) :: left
    character(len=*), intent(in) :: right
    logical retValue
    retValue = left%chr /= right
  end function compNStrChr
  
  elemental function compNChrStr(left, right) result(retValue)
    character(len=*), intent(in) :: left
    class(string),    intent(in) :: right
    logical retValue
    retValue = left /= right%chr
  end function compNChrStr
  
end module strings
