module utils

  ! A.Mereghetti (CERN, 2018-03-01)
  ! a general module, collecting some utility functions
  
contains
  
  ! ========================================================================== !
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
  ! ========================================================================== !
  real(kind=fPrec) function lininterp(x,xvals,yvals,datalen)
    
    use floatPrecision
    use numerical_constants
    use crcoall
    implicit none


    real(kind=fPrec) x, xvals(*),yvals(*)
    integer datalen
    intent(in) x,xvals,yvals,datalen
    
    integer ii
    real(kind=fPrec) dydx, y0
    
    ! Sanity checks
    if (datalen .le. 0) then
       write(lout,*) "**** ERROR in lininterp() ****"
       write(lout,*) "datalen was 0!"
       call prror(-1)
    end if
    if (x .lt. xvals(1) .or. x .gt. xvals(datalen)) then
       write(lout,*) "**** ERROR in lininterp() ****"
       write(lout,*) "x =",x, "outside range", xvals(1),xvals(datalen)
       call prror(-1)
    end if

    ! Find the right indexes i1 and i2
    ! Special case: first value at first point
    if (x .eq. xvals(1)) then
       lininterp = yvals(1)
       return
    end if
    
    do ii=1, datalen-1
       if (xvals(ii) .ge. xvals(ii+1)) then
          write (lout,*) "**** ERROR in lininterp() ****"
          write (lout,*) "xvals should be in increasing order"
          write (lout,*) "xvals =", xvals(:datalen)
          call prror(-1)
       end if
           
       if (x .le. xvals(ii+1)) then
          ! We're in the right interval
          dydx = (yvals(ii+1)-yvals(ii)) / (xvals(ii+1)-xvals(ii))
          y0   = yvals(ii) - dydx*xvals(ii)
          lininterp = dydx*x + y0
          return
       end if
    end do
        
    ! We didn't return yet: Something wrong
    write (lout,*) "****ERROR in lininterp() ****"
    write (lout,*) "Reached the end of the function"
    write (lout,*) "This should not happen, please contact developers"
    call prror(-1)
    
  end function lininterp

end module utils
