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

! ================================================================================================ !
!  K. Koelbig, CERN
!  Modification of wwerf, real complex error function, taken from mad8
!  Input:
!    xx, yy    (real)    argument to cerf.
!  Output:
!    wx, wy    (real)    function result.
! ================================================================================================ !
subroutine errf(xx, yy, wx, wy)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer

  real(kind=fPrec), intent(in)  :: xx
  real(kind=fPrec), intent(in)  :: yy
  real(kind=fPrec), intent(out) :: wx
  real(kind=fPrec), intent(out) :: wy

  real(kind=fPrec), parameter :: cc   = 1.12837916709551_fPrec ! FIXME: Should use two/pisqrt
  real(kind=fPrec), parameter :: xlim = 5.33_fPrec
  real(kind=fPrec), parameter :: ylim = 4.29_fPrec

  real(kind=fPrec) h,q,rx(33),ry(33),saux,sx,sy,tn,tx,ty,x,xh,xl,y,yh
  integer n,nc,nu

  x = abs(xx)
  y = abs(yy)
  if(y < ylim .and. x < xlim) then
    q  = (one-y/ylim)*sqrt(one-(x/xlim)**2)
    h  = one/(3.2_fPrec*q)
    nc = 7+int(23.0_fPrec*q)
    xl = exp_mb((1-nc)*log_mb(h))
    xh = y+half/h
    yh = x
    nu = 10+int(21.0_fPrec*q)
    rx(nu+1) = zero
    ry(nu+1) = zero
    do n=nu,1,-1
      tx = xh+real(n,fPrec)*rx(n+1)
      ty = yh-real(n,fPrec)*ry(n+1)
      tn = tx**2+ty**2
      rx(n) = (half*tx)/tn
      ry(n) = (half*ty)/tn
    end do
    sx = zero
    sy = zero
    do n=nc,1,-1
      saux = sx+xl
      sx   = rx(n)*saux-ry(n)*sy
      sy   = rx(n)*sy+ry(n)*saux
      xl   = h*xl
    end do
    wx = cc*sx
    wy = cc*sy
  else
    xh = y
    yh = x
    rx(1) = zero
    ry(1) = zero
    do n=9,1,-1
      tx    = xh+real(n,fPrec)*rx(1)
      ty    = yh-real(n,fPrec)*ry(1)
      tn    = tx**2+ty**2
      rx(1) = (half*tx)/tn
      ry(1) = (half*ty)/tn
    end do
    wx = cc*rx(1)
    wy = cc*ry(1)
  end if

  if(yy < zero) then
    wx = (two*exp_mb(y**2-x**2))*cos_mb((two*x)*y)-wx
    wy = ((-one*two)*exp_mb(y**2-x**2))*sin_mb((two*x)*y)-wy
    if(xx > zero) wy = -one*wy
  else
    if(xx < zero) wy = -one*wy
  end if

end subroutine errf

end module mod_utils
