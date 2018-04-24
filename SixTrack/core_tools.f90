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
    procedure, public,  pass(this) :: get => getString
    procedure, public,  pass(this) :: set => setString
    generic,   public              :: assignment(=) => assignStrStr, assignStrChr
    generic,   public              :: operator(+)   => appendStrStr, appendStrChr
    procedure, private, pass(left) :: assignStrStr, assignStrChr
    procedure, private, pass(left) :: appendStrStr, appendStrChr
  end type string
  
  interface string
    module procedure stringConstructor
  end interface string
  
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
  
  type(string) elemental function appendStrStr(left, right)
    class(string),    intent(in)  :: left
    class(string),    intent(in)  :: right
    character(len=:), allocatable :: tmpChr
    integer nL, nR
    nL = len(left%chr)
    nR = len(right%chr)
    allocate(character(nL+nR) :: tmpChr)
    tmpChr(1:nL)     = left%chr
    tmpChr(nL+1:nR)  = right%chr
    appendStrStr%chr = tmpChr
    deallocate(tmpChr)
  end function appendStrStr
  
  type(string) elemental function appendStrChr(left, right)
    class(string),    intent(in)  :: left
    character(len=*), intent(in)  :: right
    character(len=:), allocatable :: tmpChr
    integer nL, nR
    nL = len(left%chr)
    nR = len(right)
    allocate(character(nL+nR) :: tmpChr)
    tmpChr(1:nL)     = left%chr
    tmpChr(nL+1:nR)  = right
    appendStrChr%chr = tmpChr
    deallocate(tmpChr)
  end function appendStrChr
  
end module strings
