! ================================================================================================ !
! A.Mereghetti (CERN, 2018-03-01)
! A general module, collecting some utility functions
! ================================================================================================ !
module mod_utils

  implicit none

contains

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
real(kind=fPrec) function lininterp(x,xvals,yvals,datalen)

  use crcoall
  use numerical_constants, only : zero
  use floatPrecision

  integer,          intent(in) :: datalen
  real(kind=fPrec), intent(in) :: x, xvals(1:datalen), yvals(1:datalen)

  integer ii
  real(kind=fPrec) dydx, y0

  lininterp = zero ! -Wmaybe-uninitialized

  ! Sanity checks
  if(datalen <= 0) then
    write(lerr,"(a)") "UTILS> ERROR lininterp: datalen was 0!"
    call prror
  end if
  if(x < xvals(1) .or. x > xvals(datalen)) then
    write(lerr,"(3(a,e16.9))") "UTILS> ERROR lininterp: x = ",x," outside range",xvals(1)," - ",xvals(datalen)
    call prror
  end if

  ! Find the right indexes i1 and i2
  ! Special case: first value at first point
  if(x == xvals(1)) then
    lininterp = yvals(1)
    return
  end if

  do ii=1, datalen-1
    if(xvals(ii) >= xvals(ii+1)) then
      write(lerr,"(a)") "UTILS> ERROR lininterp: xvals should be in increasing order"
      call prror
    end if

    if(x <= xvals(ii+1)) then
      ! We're in the right interval
      dydx = (yvals(ii+1)-yvals(ii)) / (xvals(ii+1)-xvals(ii))
      y0   = yvals(ii) - dydx*xvals(ii)
      lininterp = dydx*x + y0
      return
    end if
  end do

  ! We didn't return yet: Something wrong
  write(lerr,"(a)") "UTILS> ERROR lininterp: Reached the end of the function. This should not happen, please contact developers"
  call prror

end function lininterp

end module mod_utils
