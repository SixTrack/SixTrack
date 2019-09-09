! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  A collection of common search routines for SixTrack
!  Created: 2019-06-05
!  Updated: 2019-06-07
! ================================================================================================ !
module mod_find

  use floatPrecision

  implicit none

  private

  public :: find_singElemFromName
  public :: find_struElemFromName
  public :: find_struElemFromPos

contains

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Returns the single element index of a given element name. Is -1 if not found.
!  Created: 2019-06-05
!  Updated: 2019-06-05
! ================================================================================================ !
integer function find_singElemFromName(elemName)

  use mod_common

  character(len=*), intent(in) :: elemName
  integer i

  find_singElemFromName = -1
  do i=1,il
    if(bez(i) == elemName) then
      find_singElemFromName = i
      exit
    end if
  end do

end function find_singElemFromName

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Returns the structure element index of a given element name. Is -1 if not found.
!  Requires a start and end index as these are not necessarily the same during parsing as after
!  reshuffling of the lattice.
!  Created: 2019-06-05
!  Updated: 2019-06-05
! ================================================================================================ !
integer function find_struElemFromName(elemName, iStart, iEnd)

  use mod_common

  character(len=*), intent(in) :: elemName
  integer,          intent(in) :: iStart
  integer,          intent(in) :: iEnd

  integer i

  find_struElemFromName = -1

  ! If we don't have a multicolumn struct block, the name is not unique
  if(strumcol .eqv. .false.) return

  do i=iStart,iEnd
    if(bezs(i) == elemName) then
      find_struElemFromName = i
      exit
    end if
  end do

end function find_struElemFromName

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Returns the structure element index of a given element position. Is -1 if not found.
!  Created: 2019-06-05
!  Updated: 2019-06-05
! ================================================================================================ !
integer function find_struElemFromPos(sPos, selectMethod)

  use mod_common, only : dcum, iu

  real(kind=fPrec), intent(in) :: sPos
  integer,          intent(in) :: selectMethod

  integer i, iA, iB, iC, iD, iN

  find_struElemFromPos = -1

  iA = 0
  iB = iu+1
  iC = 0
  iD = 0
  iN = 2*nint(log(real(iu))) + 5 ! Give some extra room, but 2*log(n) should be enough

  ! Binary search
  do i=1,iN
    iD = iB - iA
    if(iD == 1) exit
    iC = iA + iD/2
    if(sPos > dcum(iC)) then
      iA = iC
    else
      iB = iC
    end if
  end do

  if(selectMethod == -1) then    ! Return lowest index
    find_struElemFromPos = iA
  elseif(selectMethod == 1) then ! Return highest index
    find_struElemFromPos = iB
  elseif(selectMethod == 0) then ! Return nearest index
    if(sPos - dcum(iA) < dcum(iB) - sPos) then
      find_struElemFromPos = iA
    else
      find_struElemFromPos = iB ! If they are equal, we get the last element
    end if
  end if

end function find_struElemFromPos

end module mod_find
