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
!  A.Mereghetti, CERN, BE-ABP-HSS
!  Last modified: 15-04-2019
!  Make some sanity checks on arrays to be interpolated
! ================================================================================================ !
logical function checkArray(xvals,datalen,lprint)

  use crcoall
  use floatPrecision

  integer,           intent(in) :: datalen
  real(kind=fPrec),  intent(in) :: xvals(1:datalen)
  logical, optional, intent(in) :: lprint

  integer ii
  logical llprint

  llprint = .true.
  if(present(lprint)) llprint = lprint
  checkArray=.true.

  if(datalen <= 0) then
    checkArray=.false.
    if (llprint) write(lerr,"(a)") "UTILS> ERROR checkArray: datalen was <=0!"
  end if
  do ii=1, datalen-1
    if(xvals(ii) >= xvals(ii+1)) then
      checkArray=.false.
      if (llprint) then
        write(lerr,"(a)") "UTILS> ERROR checkArray: xvals should be in increasing order"
        write(lerr,"(a,i0)") "UTILS>       it is not the case at ii=",ii
      end if
    end if
  end do
end function checkArray

! ================================================================================================ !
!  A.Mereghetti, CERN, BE-ABP-HSS
!  Last modified: 24-04-2019
!  Given an array xvals(1:datalen), and given a value x, it returns the index in the xvals array
!      such that x is between xvals(huntBin) and xvals(huntBin+1). xvals(1:datalen) must be
!      monotonic in increasing/decreasing order - not checked!
!  -1 is returned to indicate that x is out of range. jlo in input is taken as initial guess
!  subroutine from Numerical Recipes in Fortran 77
! ================================================================================================ !
integer function huntBin(x,xvals,datalen,jlo)
  
  use floatPrecision
  
  integer,           intent(in) :: datalen
  real(kind=fPrec),  intent(in) :: x, xvals(1:datalen)
  integer, optional, intent(in) :: jlo

  integer inc, jhi, jm, klo
  logical ascnd

  huntBin=-1
  klo = -1
  if(present(jlo)) klo = jlo
  ascnd=xvals(datalen).ge.xvals(1) ! True if ascending order of table, false otherwise.
  
  if(klo <= 0 .or. datalen<klo)then
     ! Input guess not useful. Go immediately to bisection.
     klo=0
     jhi=datalen+1
     goto 3
  endif
  
  ! use guessed position
  
  ! set the hunting increment.
  inc=1
  if(x.ge.xvals(klo).eqv.ascnd)then
     ! hunt up
1    jhi=klo+inc
     if(jhi.gt.datalen)then
        ! Done hunting, since off end of table.
        jhi=datalen+1
     else if(x.ge.xvals(jhi).eqv.ascnd) then
        ! Not done hunting, so double the increment and try again.
        klo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     ! hunt down
     jhi=klo
2    klo=jhi-inc
     if(klo.lt.1)then
        ! Done hunting, since off end of table.
        klo=0
     else if(x.lt.xvals(klo).eqv.ascnd)then
        ! Not done hunting, so double the increment and try again.
        jhi=klo
        inc=inc+inc
        goto 2
     endif
  endif

  ! begin the final bisection phase:
3 if(jhi-klo.eq.1)then
     if(x.eq.xvals(datalen)) klo=datalen-1
     if(x.eq.xvals(1)) klo=1
     huntBin=klo
     return
  endif
  jm=(jhi+klo)/2
  if(x.ge.xvals(jm).eqv.ascnd)then
     klo=jm
  else
     jhi=jm
  endif
  goto 3
end function huntBin

! ================================================================================================ !
!  A.Mereghetti, CERN, BE-ABP-HSS
!  Last modified: 24-04-2019
!  Given arrays xa and ya, each of length n, and given a value x, this function returns a
!        value y, and an error estimate dy. If P(x) is the polynomial of degree N−1 such that
!        P(xa_i)=ya_i, i=1,...,n , then the returned value y = P ( x ).
!  This function implements Neville's method
!  function from Numerical Recipes in Fortran 77
! ================================================================================================ !
real(kind=fPrec) function nevilleMethod(xa,ya,nn,x,dy)
  
  use crcoall
  use numerical_constants, only : zero
  use floatPrecision
  
  integer,           intent(in)  :: nn
  real(kind=fPrec),  intent(in)  :: x, xa(1:nn), ya(1:nn)
  real(kind=fPrec),  intent(out) :: dy

  real(kind=fPrec) den, dif, dift, ho, hp, ww, Cs(1:nn), Ds(1:nn)
  integer ii, mm, ns

  ! find the index ns of the closest table entry,
  !    and initialize the tableau of Cs and Ds.
  ns=1
  dif=abs(x-xa(1))
  do ii=1,nn
    dift=abs(x-xa(ii))
    if (dift.lt.dif) then
      ns=ii
      dif=dift
    endif
    Cs(ii)=ya(ii)
    Ds(ii)=ya(ii)
  enddo
  
  nevilleMethod=ya(ns) ! This is the initial approximation to y.
  ns=ns-1
  do mm=1,nn-1      ! For each column of the tableau,
    do ii=1,nn-mm   ! we loop over the current c’s and d’s and update them.
      ho=xa(ii)-x
      hp=xa(ii+mm)-x
      ww=Cs(ii+1)-Ds(ii)
      den=ho-hp
      if (den.eq.zero) then
        ! This error can occur only if two input xa’s are (to within roundoff) identical.
        write(lerr,"(A)") "UTILS> ERROR nevilleMethod den.eq.zero"
        call prror
      end if
      den=ww/den
      ! update Cs and Ds
      Ds(ii)=hp*den
      Cs(ii)=ho*den
    enddo
    ! After each column in the tableau is completed, we decide
    !   which correction, C or D, we want to add to our accu-
    !   mulating value of y, i.e., which path to take through
    !   the tableau—forking up or down. We do this in such a
    !   way as to take the most “straight line” route through the
    !   tableau to its apex, updating ns accordingly to keep track
    !   of where we are. This route keeps the partial approxima-
    !   tions centered (insofar as possible) on the target x. The
    !   last dy added is thus the error indication.
    if (2*ns.lt.nn-mm)then
      dy=Cs(ns+1)
    else
      dy=Ds(ns)
      ns=ns-1
    endif
    nevilleMethod=nevilleMethod+dy
  enddo

end function nevilleMethod

! ================================================================================================ !
!  A.Mereghetti, CERN, BE-ABP-HSS
!  Last modified: 15-04-2019
!  Interpolate (xvals,yvals) with a polynom through mpoints data, and find the value y at x
! ================================================================================================ !
real(kind=fPrec) function polinterp(x,xvals,yvals,datalen,mpoints,jguess)

  use crcoall, only : lout
  use numerical_constants, only : zero
  use floatPrecision

  integer,          intent(in)    :: datalen, mpoints
  integer,          intent(inout) :: jguess
  real(kind=fPrec), intent(in)    :: x, xvals(1:datalen), yvals(1:datalen)

  integer ii, jj, jMin, jMax
  real(kind=fPrec) dy

  polinterp = zero ! -Wmaybe-uninitialized
  
  ! get bin
  jj=huntBin(x,xvals,datalen,jguess)
  ! get bin extremes (stay within range)
  jMin=min(max(jj-(mpoints-1)/2,1),datalen+1-mpoints)
  jMax=min(jMin+mpoints-1,datalen)

  ! actually interpolate
  polinterp = nevilleMethod(xvals(jMin:jMax),yvals(jMin:jMax),mpoints,x,dy)

#ifdef DEBUG
  write(lout,'(/a,1pe16.9)')   'UTILS> interpolating arrays for x= ',x
  write(lout,'(a,i0)')         'UTILS>    arrays of length: ',datalen
  write(lout,'(a,i0)')         'UTILS>    initial guess: jguess=',jguess
  write(lout,'(a,i0)')         'UTILS>    index found at:  jj=',jj
  write(lout,'(2(a,1pe16.9))') 'UTILS>    i.e. between x= ',xvals(jj),' and x= ',xvals(jj+1)
  write(lout,'(3(a,i0))')      'UTILS>    using ',mpoints,' points: jMin= ',jMin,' and jMax= ',jMax
  write(lout,'(2(a,1pe16.9))') 'UTILS>    i.e. between x= ',xvals(jMin),'and x= ',xvals(jMax)
  write(lout,'(a,1pe16.9)')    'UTILS> calculated value: ',polinterp
  write(lout,'(a)')            'UTILS> arrays:'
  do ii=1,datalen
    write(lout,'(i0,2(1X,1pe16.9))') ii, xvals(ii), yvals(ii)
  enddo
#endif
  
  ! update guess for next call
  jguess=jj

end function polinterp

! ================================================================================================ !
!  A.Mereghetti, CERN, BE-ABP-HSS
!  Last modified: 15-04-2019
!  Given arrays xa(1:nn) and ya(1:nn) of length nn containing a tabulated function ya_i=f(xa_i), 
!        this routine returns an array of coefficients cof(1:nn), also of length nn, such that
!    ya_i = sum_j cof_j xa^{j−1}_i
!  from Numerical Recipes in Fortran 77
! ================================================================================================ !
subroutine polcof(xa,ya,nn,cof)
  
  use numerical_constants, only : zero, c1e38
  use floatPrecision

  integer,          intent(in)    :: nn
  real(kind=fPrec), intent(in)    :: xa(1:nn), ya(1:nn)
  real(kind=fPrec), intent(out)   :: cof(1:nn)

  integer ii, jj, kk
  real(kind=fPrec) dy, xmin, x(1:nn), y(1:nn)
  
  x(1:nn)=xa(1:nn)
  y(1:nn)=ya(1:nn)
  do jj=1,nn
    cof(jj)=nevilleMethod(x,y,nn+1-jj,zero,dy)
    xmin=c1e38
    kk=0
    do ii=1,nn+1-jj
      ! Find the remaining x_i of smallest absolute value,
      if (abs(x(ii)).lt.xmin)then
        xmin=abs(x(ii))
        kk=ii
      endif
      if(x(ii).ne.0.) y(ii)=(y(ii)-cof(jj))/x(ii) ! (meanwhile reducing all the terms)
    enddo
    ! and eliminate it.
    y(kk:nn-jj)=y(kk+1:nn+1-jj)
    x(kk:nn-jj)=x(kk+1:nn+1-jj)
  enddo

end subroutine polcof

! ================================================================================================ !
!  A.Mereghetti, CERN, BE-ABP-HSS
!  Last modified: 16-04-2019
!  Integrate (xvals,yvals) with polynoms through mpoints data, and find the value y at x
! ================================================================================================ !
real(kind=fPrec) function polintegrate(xvals,yvals,datalen,mpoints,order,cumul,rmin,rmax)

  use crcoall, only : lout, lerr
  use numerical_constants, only : zero, one, two, four, pi
  use floatPrecision

  integer,                    intent(in)    :: datalen, mpoints, order
  real(kind=fPrec),           intent(in)    :: xvals(1:datalen), yvals(1:datalen)
  real(kind=fPrec),           intent(out)   :: cumul(1:datalen)
  real(kind=fPrec), optional, intent(in)    :: rmin, rmax

  integer jj, jMin, jMax, kk, kMin, kMax
  real(kind=fPrec) xmin, xmax, fmax, fmin, tmin, tmax, tmpInt
  real(kind=fPrec) cof(1:mpoints)

  polintegrate = zero ! -Wmaybe-uninitialized
  cumul(1:datalen)=zero
  
  ! get bins
  jMin=1
  if(present(rmin)) jMin=huntBin(rmin,xvals,datalen,-1)
  jMax=datalen
  if(present(rmax)) jMax=huntBin(rmax,xvals,datalen,-1)

#ifdef DEBUG
  write(lout,"(/a,5(1X,i0))") "UTILS> DEBUG polintegrate start: ",datalen,mpoints,order,jMin,jMax
#endif

  do jj=jMin+1,jMax ! loop over data points
     kMin=min(max(jj-1-(mpoints-1)/2,1),datalen+1-mpoints)
     kMax=min(kMin+mpoints-1,datalen)
     call polcof(xvals(kMin:kMax),yvals(kMin:kMax),mpoints,cof)
     xmin=xvals(jj-1)
     if (present(rmin).and.jj==jMin+1) xmin=rmin
     xmax=xvals(jj)
     if (present(rmax).and.jj==jMax) xmax=rmax
     fmax=zero
     fmin=zero
     tmin=one
     tmax=one
     do kk=1,order-1
        tmin=tmin*xmin
        tmax=tmax*xmax
     end do
     tmpInt=zero
#ifdef DEBUG
     write(lout,"(/a,i4,6(1X,1pe16.9))") "UTILS> DEBUG polintegrate jj-step: ", &
            jj, tmax, tmin, fmax, fmin, xmax, xmin
#endif
     do kk=1,kMax-kMin+1
        tmin=tmin*xmin
        tmax=tmax*xmax
        fmin=fmin+(tmin*cof(kk))/real(kk+order-1,fPrec)
        fmax=fmax+(tmax*cof(kk))/real(kk+order-1,fPrec)
#ifdef DEBUG
        write(lout,"(a,i4,5(1X,1pe16.9))") "UTILS> DEBUG polintegrate kk-step: ", &
               kk, tmax, tmin, fmax, fmin, cof(kk)
#endif
     end do
     select case(order)
     case(1)
        tmpInt=fmax-fmin
     case(2)
        tmpInt=two*(pi*(fmax-fmin))
     case(3)
        tmpInt=four*(pi*(fmax-fmin))
     case default
        write(lerr,"(a,i0)") "UTILS> ERROR polintegrate: order not recognised:",order
        write(lerr,"(a)")    "UTILS>       recognised values: 1D (1), 2D (2), 3D (3)"
        call prror
        return
     end select
     polintegrate=polintegrate+tmpInt
     cumul(jj)=polintegrate
#ifdef DEBUG
     write(lout,"(a,i4,6(1X,1pe16.9))") "UTILS> DEBUG polintegrate jj-step: ", &
            jj, tmax, tmin, fmax, fmin, tmpInt, cumul(jj)
#endif
  end do

  ! fill outer range
  cumul(jMax+1:datalen)=cumul(jMax)

#ifdef DEBUG
  write(lout,"(a,1pe16.9)") "UTILS> DEBUG polintegrate end: ",polintegrate
#endif
    
end function polintegrate
  
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

subroutine wwerf(x, y, wr, wi)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  real(kind=fPrec), intent(in)  :: x
  real(kind=fPrec), intent(in)  :: y
  real(kind=fPrec), intent(out) :: wr
  real(kind=fPrec), intent(out) :: wi

  integer n
  real(kind=fPrec) rr(37),ri(37),sr0,sr,si,tr,ti,vi,vr,xa,xl,ya,zhi,zhr

  real(kind=fPrec), parameter :: c1 = 74.0_fPrec/c1e1
  real(kind=fPrec), parameter :: c2 = 83.0_fPrec/c1e1
  real(kind=fPrec), parameter :: c3 = c1e1/32.0_fPrec
  real(kind=fPrec), parameter :: c4 = 16.0_fPrec/c1e1
  real(kind=fPrec), parameter :: c  = two/pisqrt
  real(kind=fPrec), parameter :: p  = 46768052394588893.3825_fPrec

  xa = abs(x)
  ya = abs(y)
  if(ya < c1 .and. xa < c2) then
    zhr = ya+c4
    zhi = xa
    rr(37) = zero
    ri(37) = zero
    do n=36,1,-1
      tr = zhr+real(n,fPrec)*rr(n+1)
      ti = zhi-real(n,fPrec)*ri(n+1)
      rr(n) = (half*tr)/(tr**2+ti**2)
      ri(n) = (half*ti)/(tr**2+ti**2)
    end do
    xl = p
    sr = zero
    si = zero
    do n=33,1,-1
      xl  = c3*xl
      sr0 = rr(n)*(sr+xl)-ri(n)*si
      si  = rr(n)*si+ri(n)*(sr+xl)
      sr  = sr0
    end do
    vr = c*sr
    vi = c*si
  else
    zhr = ya
    zhi = xa
    rr(1) = zero
    ri(1) = zero
    do n=9,1,-1
      tr = zhr+real(n,fPrec)*rr(1)
      ti = zhi-real(n,fPrec)*ri(1)
      rr(1) = (half*tr)/(tr**2+ti**2)
      ri(1) = (half*ti)/(tr**2+ti**2)
    end do
    vr = c*rr(1)
    vi = c*ri(1)
  end if
  if(ya == zero) then
    vr = exp_mb(-one*xa**2)
  end if
  if(y < zero) then
    vr = (two*exp_mb(ya**2-xa**2))*cos_mb((two*xa)*ya)-vr
    vi = (-two*exp_mb(ya**2-xa**2))*sin_mb((two*xa)*ya)-vi
    if(x > zero) vi = -one*vi
  else
    if(x < zero) vi = -one*vi
  end if
  wr = vr
  wi = vi

end subroutine wwerf

! ================================================================================================ !
!  subroutine wzsubv
!
!  This subroutine sets u=real(w(z)) and v=imag(w(z)), where z=x+i*y and
!  where w(z) is the complex error function defined by formula 7.1.3 in
!  "Handbook of Mathematical functions [eds. M.Abramowitz & I.A.Stegun,
!  Washington, 1966].  The absolute error of the computed value is less
!  than 1E-8.
!
!  *** Note.  Subroutine WZSET must have been called before this sub-
!  routine can be used.
!
!  For (x,y) inside the rectangle with opposite corners (xcut,0) and
!  (0,ycut), where xcut and ycut have been set by WZSET, an interpo-
!  lation formula is used.  For (x,y) outside this rectangle, a two-
!  term rational approximation is used.
!
!  (G.A.Erskine, 29.09.1997)
!
!  Vectorised for up to 64 argument values by E.McIntosh, 30.10.1997.
!  Much impoved using short vector buffers Eric 1st May, 2014.
!
!  Third-order divided-difference interpolation over the corners of a
!  square [e.g. formula (2.5.1) in "Introduction to Numerical Analysis"
!  (F.B.Hildebrand New York, 1957), but with complex nodes and
!  function values].
!
!  In the interpolation formula the corners of the grid square contain-
!  ing (x,y) are numbered (0,0)=3, (h,0)=4, (h,h)=1, (0,h)=2.
!  Identifiers d, dd and ddd denote divided-differences of orders 1, 2
!  and 3 respectively, and a preceding 't' indicates twice the value.
!
!
!  Two-term rational approximation to w(z) [Footnote to Table 7.9
!  in "Handbook of Mathematical Functions (eds. M.Abramowitz &
!  I.A.Stegun, Washington, 1966), but with additional digits in
!  the constants]:
!              u+i*v = i*z*( a1/(z**2-b1) + a2/(z**2-b2) ).
!  Maximum absolute error:
!        <1.E-6  for  x>=4.9  or  y>=4.4
!        <1.E-7  for  x>=6.1  or  y>=5.7
!        <1.E-8  for  x>=7.8  or  y>=7.5
!
! ================================================================================================ !
subroutine wzsubv(n,vx,vy,vu,vv)

  use parpro, only : npart
  use parbeam
  use floatPrecision
  use numerical_constants

  integer,          intent(in)  :: n
  real(kind=fPrec), intent(in)  :: vx(*)
  real(kind=fPrec), intent(in)  :: vy(*)
  real(kind=fPrec), intent(out) :: vu(*)
  real(kind=fPrec), intent(out) :: vv(*)

  integer i,j,k,vmu,vnu
  integer in,out,ins(npart),outs(npart)
  real(kind=fPrec) vd12i,vd12r,vd23i,vd23r,vd34i,vd34r,vp,vq,vqsq,vr,vsimag,vsreal,vt,vtdd13i,      &
    vtdd13r,vtdd24i,vtdd24r,vtdddi,vtdddr,vti,vtr,vusum,vusum3,vvsum,vvsum3,vw1i,vw1r,vw2i,vw2r,    &
    vw3i,vw3r,vw4i,vw4r,vxh,vxhrel,vyh,vyhrel,xx,yy

  real(kind=fprec), parameter :: a1 = 0.5124242248_fPrec
  real(kind=fprec), parameter :: a2 = 0.0517653588_fPrec
  real(kind=fprec), parameter :: b1 = 0.2752551286_fPrec
  real(kind=fprec), parameter :: b2 = 2.7247448714_fPrec
  real(kind=fPrec), parameter :: xm = 1e16_fPrec

  in  = 0
  out = 0
  do i=1,n
    if (vx(i).ge.xcut.or.vy(i).ge.ycut) then
      out=out+1
      outs(out)=i
      if(out == npart) then
        ! Everything outside the rectangle so approximate
        do j=1,out
          xx = vx(outs(j))
          yy = vy(outs(j))
          if(xx >= xm) xx = xm
          if(yy >= xm) yy = xm
          vp = xx**2-yy**2
          vq = (two*xx)*yy
          vqsq = vq**2
          ! First term.
          vt = vp-b1
          vr = a1/(vt**2+vqsq)
          vsreal = vr*vt
          vsimag = -vr*vq
          ! Second term
          vt = vp-b2
          vr = a2/(vt**2+vqsq)
          vsreal = vsreal+vr*vt
          vsimag = vsimag-vr*vq
          ! Multiply by i*z.
          vu(outs(j)) = -(yy*vsreal+xx*vsimag)
          vv(outs(j)) = xx*vsreal-yy*vsimag
        end do
        out = 0
      end if
    else
      in = in+1
      ins(in) = i
      if(in == npart) then
      ! Everything inside the square, so interpolate
        do j=1,in
          vxh = hrecip*vx(ins(j))
          vyh = hrecip*vy(ins(j))
          vmu = int(vxh)
          vnu = int(vyh)
          !  Compute divided differences.
          k = 2 + vmu + vnu*kstep
          vw4r = wtreal(k)
          vw4i = wtimag(k)
          k = k - 1
          vw3r = wtreal(k)
          vw3i = wtimag(k)
          vd34r = vw4r - vw3r
          vd34i = vw4i - vw3i
          k = k + kstep
          vw2r = wtreal(k)
          vw2i = wtimag(k)
          vd23r = vw2i - vw3i
          vd23i = vw3r - vw2r
          vtr = vd23r - vd34r
          vti = vd23i - vd34i
          vtdd24r = vti - vtr
          vtdd24i = -one* ( vtr + vti )
          k = k + 1
          vw1r = wtreal(k)
          vw1i = wtimag(k)
          vd12r = vw1r - vw2r
          vd12i = vw1i - vw2i
          vtr = vd12r - vd23r
          vti = vd12i - vd23i
          vtdd13r = vtr + vti
          vtdd13i = vti - vtr
          vtdddr = vtdd13i - vtdd24i
          vtdddi = vtdd24r - vtdd13r
          ! Evaluate polynomial.
          vxhrel = vxh - real(vmu,fPrec)
          vyhrel = vyh - real(vnu,fPrec)
          vusum3=half*(vtdd13r+(vxhrel*vtdddr-vyhrel*vtdddi))
          vvsum3=half*(vtdd13i+(vxhrel*vtdddi+vyhrel*vtdddr))
          vyhrel = vyhrel - one
          vusum=vd12r+(vxhrel*vusum3-vyhrel*vvsum3)
          vvsum=vd12i+(vxhrel*vvsum3+vyhrel*vusum3)
          vxhrel = vxhrel - one
          vu(ins(j))=vw1r+(vxhrel*vusum-vyhrel*vvsum)
          vv(ins(j))=vw1i+(vxhrel*vvsum+vyhrel*vusum)
        end do
        in = 0
      end if
    end if
  end do

  ! Everything outside the rectangle so approximate
  do j=1,out
    xx = vx(outs(j))
    yy = vy(outs(j))
    if(xx >= xm) xx = xm
    if(yy >= xm) yy = xm
    vp = xx**2-yy**2
    vq = (two*xx)*yy
    vqsq = vq**2
    ! First term
    vt = vp-b1
    vr = a1/(vt**2+vqsq)
    vsreal = vr*vt
    vsimag = -vr*vq
    ! Second term
    vt = vp-b2
    vr = a2/(vt**2+vqsq)
    vsreal = vsreal+vr*vt
    vsimag = vsimag-vr*vq
    ! Multiply by i*z.
    vu(outs(j)) = -(yy*vsreal+xx*vsimag)
    vv(outs(j)) = xx*vsreal-yy*vsimag
  end do

  ! Everything inside the square, so interpolate
  do j=1,in
    vxh = hrecip*vx(ins(j))
    vyh = hrecip*vy(ins(j))
    vmu = int(vxh)
    vnu = int(vyh)
    ! Compute divided differences
    k = 2 + vmu + vnu*kstep
    vw4r = wtreal(k)
    vw4i = wtimag(k)
    k = k - 1
    vw3r = wtreal(k)
    vw3i = wtimag(k)
    vd34r = vw4r - vw3r
    vd34i = vw4i - vw3i
    k = k + kstep
    vw2r = wtreal(k)
    vw2i = wtimag(k)
    vd23r = vw2i - vw3i
    vd23i = vw3r - vw2r
    vtr = vd23r - vd34r
    vti = vd23i - vd34i
    vtdd24r = vti - vtr
    vtdd24i = -one* ( vtr + vti )
    k = k + 1
    vw1r = wtreal(k)
    vw1i = wtimag(k)
    vd12r = vw1r - vw2r
    vd12i = vw1i - vw2i
    vtr = vd12r - vd23r
    vti = vd12i - vd23i
    vtdd13r = vtr + vti
    vtdd13i = vti - vtr
    vtdddr = vtdd13i - vtdd24i
    vtdddi = vtdd24r - vtdd13r
    ! Evaluate polynomial
    vxhrel = vxh - real(vmu,fPrec)
    vyhrel = vyh - real(vnu,fPrec)
    vusum3=half*(vtdd13r+(vxhrel*vtdddr-vyhrel*vtdddi))
    vvsum3=half*(vtdd13i+(vxhrel*vtdddi+vyhrel*vtdddr))
    vyhrel = vyhrel - one
    vusum=vd12r+(vxhrel*vusum3-vyhrel*vvsum3)
    vvsum=vd12i+(vxhrel*vvsum3+vyhrel*vusum3)
    vxhrel = vxhrel - one
    vu(ins(j))=vw1r+(vxhrel*vusum-vyhrel*vvsum)
    vv(ins(j))=vw1i+(vxhrel*vvsum+vyhrel*vusum)
  end do

end subroutine wzsubv

! ================================================================================================ !
!  subroutine wzsub
!
!  This subroutine sets u=real(w(z)) and v=imag(w(z)), where z=x+i*y and
!  where w(z) is the complex error function defined by formula 7.1.3 in
!  "Handbook of Mathematical functions [eds. M.Abramowitz & I.A.Stegun,
!  Washington, 1966].  The absolute error of the computed value is less
!  than 1E-8.
!
!  *** Note.  Subroutine WZSET must have been called before this sub-
!  routine can be used.
!
!  For (x,y) inside the rectangle with opposite corners (xcut,0) and
!  (0,ycut), where xcut and ycut have been set by WZSET, an interpo-
!  lation formula is used.  For (x,y) outside this rectangle, a two-
!  term rational approximation is used.
!
!  (G.A.Erskine, 29.09.1997)
!
!
!  Third-order divided-difference interpolation over the corners of a
!  square [e.g. formula (2.5.1) in "Introduction to Numerical Analysis"
!  (F.B.Hildebrand New York, 1957), but with complex nodes and
!  function values].
!
!  In the interpolation formula the corners of the grid square contain-
!  ing (x,y) are numbered (0,0)=3, (h,0)=4, (h,h)=1, (0,h)=2.
!  Identifiers d, dd and ddd denote divided-differences of orders 1, 2
!  and 3 respectively, and a preceding 't' indicates twice the value.
!
! ================================================================================================ !
subroutine wzsub(x,y,u,v)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use parpro
  use parbeam

  real(kind=fPrec), intent(in)  :: x
  real(kind=fPrec), intent(in)  :: y
  real(kind=fPrec), intent(out) :: u
  real(kind=fPrec), intent(out) :: v

  integer k,mu,nu
  real(kind=fPrec) d12i,d12r,d23i,d23r,d34i,d34r,p,q,qsq,r,simag,sreal,t,tdd13i,tdd13r,tdd24i,      &
    tdd24r,tdddi,tdddr,ti,tr,usum,usum3,vsum,vsum3,w1i,w1r,w2i,w2r,w3i,w3r,w4i,w4r,xh,xhrel,yh,yhrel

  real(kind=fPrec), parameter :: a1 = 0.5124242248_fPrec
  real(kind=fPrec), parameter :: a2 = 0.0517653588_fPrec
  real(kind=fPrec), parameter :: b1 = 0.2752551286_fPrec
  real(kind=fPrec), parameter :: b2 = 2.7247448714_fPrec

  if(x >= xcut .or. y >= ycut) goto 10

  xh = hrecip*x
  yh = hrecip*y
  mu = int(xh)
  nu = int(yh)
  ! Compute divided differences.
  k = 2 + mu + nu*kstep
  w4r = wtreal(k)
  w4i = wtimag(k)
  k = k - 1
  w3r = wtreal(k)
  w3i = wtimag(k)
  d34r = w4r - w3r
  d34i = w4i - w3i
  k = k + kstep
  w2r = wtreal(k)
  w2i = wtimag(k)
  d23r = w2i - w3i
  d23i = w3r - w2r
  tr = d23r - d34r
  ti = d23i - d34i
  tdd24r = ti - tr
  tdd24i = -one* ( tr + ti )
  k = k + 1
  w1r = wtreal(k)
  w1i = wtimag(k)
  d12r = w1r - w2r
  d12i = w1i - w2i
  tr = d12r - d23r
  ti = d12i - d23i
  tdd13r = tr + ti
  tdd13i = ti - tr
  tdddr = tdd13i - tdd24i
  tdddi = tdd24r - tdd13r
  ! Evaluate polynomial.
  xhrel = xh - real(mu,fPrec)
  yhrel = yh - real(nu,fPrec)
  usum3 = half*( tdd13r + ( xhrel*tdddr - yhrel*tdddi ) )
  vsum3 = half*( tdd13i + ( xhrel*tdddi + yhrel*tdddr ) )
  yhrel = yhrel - one
  usum = d12r + ( xhrel*usum3 - yhrel*vsum3 )
  vsum = d12i + ( xhrel*vsum3 + yhrel*usum3 )
  xhrel = xhrel - one
  u = w1r + ( xhrel*usum - yhrel*vsum )
  v = w1i + ( xhrel*vsum + yhrel*usum )
  return

10 continue
  !  Two-term rational approximation to w(z) [Footnote to Table 7.9
  !  in "Handbook of Mathematical Functions (eds. M.Abramowitz &
  !  I.A.Stegun, Washington, 1966), but with additional digits in
  !  the constants]:
  !              u+i*v = i*z*( a1/(z**2-b1) + a2/(z**2-b2) ).
  !  Maximum absolute error:
  !        <1.E-6  for  x>=4.9  or  y>=4.4
  !        <1.E-7  for  x>=6.1  or  y>=5.7
  !        <1.E-8  for  x>=7.8  or  y>=7.5
  !
  p = x**2-y**2
  q = (2.d0*x)*y
  qsq = q**2
  ! First term
  t = p-b1
  r = a1/(t**2+qsq)
  sreal = r*t
  simag = (-one*r)*q
  ! Second term
  t = p-b2
  r = a2/(t**2+qsq)
  sreal = sreal+r*t
  simag = simag-r*q
  ! Multiply by i*z.
  u = -one*(y*sreal+x*simag)
  v = x*sreal-y*simag

end subroutine wzsub

!  *********************************************************************
!
!  This subroutine must be called before subroutine WZSUB can be used to
!  compute values of the complex error function w(z).
!
!  Parameters xcut and ycut specify the opposite corners (xcut,0) and
!  (0,ycut) of the rectangle inside which interpolation is to be used
!  by subroutine WZSUB.
!
!  Parameter h is the side of the squares of the interpolation grid.
!
!  Parameters nx and ny must be set to the nearest integers to xcut/h
!  and ycut/h respectively (or to larger values).
!
!  Calls wwerf new version of (CERN library) WWERF (C335)
!
!  (G.A.Erskine, 29.09.1995)
!
!  *********************************************************************
subroutine wzset

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use parbeam

  integer i,j,k
  real(kind=fPrec) wi,wr,x,y

  hrecip = one/h
  kstep = nx+2
  k = 0
  do j=0,ny+1
    do i=0,nx+1
      k = k+1
      x = real(i,fPrec)*h
      y = real(j,fPrec)*h
      call wwerf(x,y,wr,wi)
      wtreal(k) = wr
      wtimag(k) = wi
    end do
  end do

end subroutine wzset

end module mod_utils
