module utils

  ! A.Mereghetti (CERN, 2018-03-01)
  ! a general module, collecting some utility functions
  use floatPrecision

  logical :: ldebug=.true.

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
    implicit none

    integer,          intent(in) :: datalen
    real(kind=fPrec), intent(in) :: x, xvals(1:datalen), yvals(1:datalen)

    integer ii
    real(kind=fPrec) dydx, y0

    lininterp = zero ! -Wmaybe-uninitialized

    ! Sanity checks
    if(datalen <= 0) then
      write(lerr,"(a)") "UTILS> ERROR lininterp: datalen was 0!"
      call prror(-1)
    end if
    if(x < xvals(1) .or. x > xvals(datalen)) then
      write(lerr,"(3(a,e16.9))") "UTILS> ERROR lininterp: x = ",x," outside range",xvals(1)," - ",xvals(datalen)
      call prror(-1)
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
        call prror(-1)
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
    call prror(-1)

  end function lininterp

  ! ================================================================================================ !
  !  A.Mereghetti, CERN, BE-ABP-HSS
  !  Last modified: 15-04-2019
  !  Make some sanity checks on arrays to be interpolated
  ! ================================================================================================ !
  logical function checkArray(xvals,datalen,lprint)

    use crcoall
    implicit none

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
      if (llprint) write(lerr,"(a)") "UTILS> ERROR checkArray: datalen was 0!"
    end if
    do ii=1, datalen-1
      if(xvals(ii) >= xvals(ii+1)) then
        checkArray=.false.
        if (llprint) write(lerr,"(a)") "UTILS> ERROR checkArray: xvals should be in increasing order"
      end if
    end do
  end function checkArray

  ! ================================================================================================ !
  !  A.Mereghetti, CERN, BE-ABP-HSS
  !  Last modified: 15-04-2019
  !  Given an array xvals(1:datalen), and given a value x, it returns a value such that x is
  !      between xvals(huntBin) and xvals(huntBin+1). xvals(1:datalen) must be monotonic in
  !      increasing/decreasing order - not checked!
  !  -1 is returned to indicate that x is out of range. jlo in input is taken as initial guess
  !  from Numerical Recipes in Fortran 77
  ! ================================================================================================ !
  integer function huntBin(x,xvals,datalen,jlo)
    
    implicit none
    
    integer,           intent(in) :: datalen
    real(kind=fPrec),  intent(in) :: x, xvals(1:datalen)
    integer, optional, intent(in) :: jlo

    integer inc, jhi, jm, klo
    logical ascnd

    huntBin=-1
    klo = -1
    if(present(jlo)) klo = jlo
    ascnd=xvals(datalen).ge.xvals(1) ! True if ascending order of table, false otherwise.
    
    if(klo.le.0.or.klo.gt.datalen)then
       ! Input guess not useful. Go immediately to bisection.
       klo=0
       jhi=datalen+1
       goto 3
    endif
    
    ! set the hunting increment.
    inc=1
    if(x.ge.xvals(klo).eqv.ascnd)then
       ! hunt up
1      jhi=klo+inc
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
2      klo=jhi-inc
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
3   if(jhi-klo.eq.1)then
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
  !  Last modified: 15-04-2019
  !  Given arrays xa and ya, each of length n, and given a value x, this function returns a
  !        value y, and an error estimate dy. If P(x) is the polynomial of degree N−1 such that
  !        P(xa_i)=ya_i, i=1,...,n , then the returned value y = P ( x ).
  !  This function implements Neville's method
  !  from Numerical Recipes in Fortran 77
  ! ================================================================================================ !
  real(kind=fPrec) function nevilleMethod(xa,ya,nn,x,dy)
    
    use crcoall
    use numerical_constants, only : zero
    implicit none
    
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

    use numerical_constants, only : zero
    implicit none

    integer,          intent(in)    :: datalen, mpoints
    integer,          intent(inout) :: jguess
    real(kind=fPrec), intent(in)    :: x, xvals(1:datalen), yvals(1:datalen)

    integer jj, jMin, jMax
    real(kind=fPrec) dy

    polinterp = zero ! -Wmaybe-uninitialized
    
    ! get bin
    jj=huntBin(x,xvals,datalen,jguess)
    ! get bin extremes (stay within range)
    jMin=min(max(jj-(mpoints-1)/2,1),datalen+1-mpoints)
    jMax=min(jMin+mpoints-1,datalen)

    ! actually interpolate
    polinterp = nevilleMethod(xvals(jMin:jMax),yvals(jMin:jMax),mpoints,x,dy)

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
    implicit none

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

    use crcoall
    use numerical_constants, only : zero, one, two, four, pi
    implicit none

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

    if (ldebug) write(lout,"(/a,5(1X,i0))") "UTILS> DEBUG polintegrate start: ",datalen,mpoints,order,jMin,jMax

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
       if (ldebug) write(lout,"(/a,i4,6(1X,1pe16.9))") "UTILS> DEBUG polintegrate jj-step: ", &
            jj, tmax, tmin, fmax, fmin, xmax, xmin
       do kk=1,kMax-kMin+1
          tmin=tmin*xmin
          tmax=tmax*xmax
          fmin=fmin+(tmin*cof(kk))/real(kk+order-1,fPrec)
          fmax=fmax+(tmax*cof(kk))/real(kk+order-1,fPrec)
          if (ldebug) write(lout,"(a,i4,5(1X,1pe16.9))") "UTILS> DEBUG polintegrate kk-step: ", &
               kk, tmax, tmin, fmax, fmin, cof(kk)
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
       if (ldebug) write(lout,"(a,i4,6(1X,1pe16.9))") "UTILS> DEBUG polintegrate jj-step: ", &
            jj, tmax, tmin, fmax, fmin, tmpInt, cumul(jj)
    end do

    ! fill outer range
    cumul(jMax+1:datalen)=cumul(jMax)

    if (ldebug) write(lout,"(a,1pe16.9)") "UTILS> DEBUG polintegrate end: ",polintegrate
    
  end function polintegrate
  
end module utils
