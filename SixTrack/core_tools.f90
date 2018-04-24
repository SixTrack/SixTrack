! =================================================================================================
!  STANDARD OUTPUT MODULE
!  Last modified: 2018-03-22
!  For CR version, this is the "buffer file" fort.92;
!  Otherwise write directly to "*" aka iso_fortran_env::output_unit (usually unit 6)
! =================================================================================================
module crcoall
  
  implicit none
  
  integer lout
  save lout
  
end module crcoall

! ================================================================================================ !
!  SixTrack String Type
! ~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-24
!  Based on: https://github.com/bceverly/fortranString
! ================================================================================================ !
module strings
  
  implicit none
  
  type, public :: string
    character(len=:), allocatable, public :: chr
  contains
    procedure, public, pass(this) :: get => getString
    procedure, public, pass(this) :: set => setString
  end type string
  
  interface string
    module procedure stringConstructor
  end interface string
  
  interface assignment(=)
    module procedure assignStrStr
    module procedure assignStrChr
    module procedure assignChrStr
  end interface
  
  interface operator(+)
    module procedure appendStrStr
    module procedure appendStrChr
  end interface
  
  interface operator(==)
    module procedure compStrStr
    module procedure compChrStr
    module procedure compStrChr
  end interface
  
contains
  
  type(string) function stringConstructor()
    stringConstructor%chr = ""
  end function stringConstructor
  
  pure function getString(this) result(retValue)
    class(string),    intent(in)  :: this
    character(len=:), allocatable :: retValue
    retValue = this%chr
  end function getString
  
  pure subroutine setString(this, setValue)
    class(string),    intent(inout) :: this
    character(len=*), intent(in)    :: setValue
    this%chr = setValue
  end subroutine setString
  
  ! ================================================================ !
  !  String Assignment
  ! ================================================================ !
  subroutine assignStrStr(left, right)
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
  !  Append Strings
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
  !  Compare Strings
  ! ================================================================ !
  logical elemental function compStrStr(left, right)
    class(string), intent(in) :: left
    class(string), intent(in) :: right
    compStrStr = left%chr == right%chr
  end function compStrStr
  
  logical elemental function compChrStr(left, right)
    character(len=*), intent(in) :: left
    class(string),    intent(in) :: right
    compChrStr = left == right%chr
  end function compChrStr
  
  logical elemental function compStrChr(left, right)
    class(string),    intent(in) :: left
    character(len=*), intent(in) :: right
    compStrChr = left%chr == right
  end function compStrChr
  
end module strings
