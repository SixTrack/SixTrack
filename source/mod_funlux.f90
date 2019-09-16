module mod_funlux

  use floatPrecision
  use numerical_constants, only : zero, two

  implicit none

  real(kind=fPrec), private, save :: tftot

contains

!cccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine funlxp (func,xfcum,x2low,x2high)
!         F. JAMES,   Sept, 1994
!
!         Prepares the user function FUNC for FUNLUX
!         Inspired by and mostly copied from FUNPRE and FUNRAN
!         except that
!    1. FUNLUX uses RANLUX underneath,
!    2. FUNLXP expands the first and last bins to cater for
!              functions with long tails on left and/or right,
!    3. FUNLXP calls FUNPCT to do the actual finding of percentiles.
!    4. both FUNLXP and FUNPCT use RADAPT for Gaussian integration.
!
      use crcoall
      implicit none
      external func
      integer ifunc,ierr
      real(kind=fPrec) x2high,x2low,xfcum,rteps,xhigh,xlow,xrange,uncert,x2,tftot1,x3,tftot2,func
      dimension xfcum(200)
      parameter (rteps=0.0002)
      save ifunc
      data ifunc/0/
      ifunc = ifunc + 1
!         FIND RANGE WHERE FUNCTION IS NON-ZERO.
      call funlz(func,x2low,x2high,xlow,xhigh)
      xrange = xhigh-xlow
      if(xrange .le. 0.)  then
        write(lout,'(A,2G15.5)') ' FUNLXP finds function range .LE.0',xlow,xhigh
        go to 900
      endif
      call radapt(func,xlow,xhigh,1,rteps,zero,tftot ,uncert)
!      WRITE(6,1003) IFUNC,XLOW,XHIGH,TFTOT
! 1003 format(' FUNLXP: integral of USER FUNCTION', i3,' from ',e12.5,' to ',e12.5,' is ',e14.6)
!
!      WRITE (6,'(A,A)') ' FUNLXP preparing ',
!     + 'first the whole range, then left tail, then right tail.'
      call funpct(func,ifunc,xlow,xhigh,xfcum,1,99,tftot,ierr)
      if (ierr .gt. 0)  go to 900
      x2 = xfcum(3)
      call radapt(func,xlow,x2,1,rteps,zero,tftot1 ,uncert)
      call funpct(func,ifunc,xlow,x2 ,xfcum,101,49,tftot1,ierr)
      if (ierr .gt. 0)  go to 900
      x3 = xfcum(98)
      call radapt(func,x3,xhigh,1,rteps,zero,tftot2 ,uncert)
      call funpct(func,ifunc,x3,xhigh,xfcum,151,49,tftot2,ierr)
      if (ierr .gt. 0)  go to 900
!      WRITE(6,1001) IFUNC,XLOW,XHIGH
! 1001 format(' FUNLXP has prepared USER FUNCTION', i3, ' between',g12.3,' and',g12.3,' for FUNLUX')

      return
  900 continue
      write(lout,*) ' Fatal error in FUNLXP. FUNLUX will not work.'
end subroutine funlxp

subroutine funpct(func,ifunc,xlow,xhigh,xfcum,nlo,nbins,tftot,ierr)
!        Array XFCUM is filled from NLO to NLO+NBINS, which makes
!        the number of values NBINS+1, or the number of bins NBINS
      use crcoall
      implicit none
      external func
      integer ierr,nbins,nlo,ifunc,nz,ibin,maxz,iz,nitmax,ihome
      real(kind=fPrec) tftot,xhigh,xlow,func,xfcum,rteps,tpctil,tz,tzmax,x,f,tcum,  &
     &x1,f1,dxmax,fmin,fminz,xincr,tincr,xbest,dtbest,tpart,x2,precis,  &
     &refx,uncert,tpart2,dtpar2,dtabs,aberr
      dimension xfcum(*)
      parameter (rteps=0.005, nz=10, maxz=20, nitmax=6,precis=1e-6)
!      DOUBLE PRECISION TPCTIL, TZ, TCUM, XINCR, DTABS,
!     &  TINCR, TZMAX, XBEST, DTBEST, DTPAR2
!
      ierr = 0
      if (tftot .le. 0.) go to 900
      tpctil = tftot/real(nbins)
      tz = tpctil/real(nz)
      tzmax = tz * 2.
      xfcum(nlo) = xlow
      xfcum(nlo+nbins) = xhigh
      x = xlow
      f = func(x)
      if (f .lt. 0.) go to 900
!         Loop over percentile bins
      do 600 ibin = nlo, nlo+nbins-2
      tcum = 0.
      x1 = x
      f1 = f
      dxmax = (xhigh -x) / nz
      fmin = tz/dxmax
      fminz = fmin
!         Loop over trapezoids within a supposed percentil
      do 500 iz= 1, maxz
      xincr = tz/max(f1,fmin,fminz)
  350 x = x1 + xincr
      f = func(x)
      if (f .lt. 0.) go to 900
      tincr = ((x-x1) * 0.5) * (f+f1)
      if (tincr .lt. tzmax) go to 370
      xincr = xincr * 0.5
      go to 350
  370 continue
      tcum = tcum + tincr
      if (tcum .ge. tpctil*0.99) go to 520
      fminz = (tz*f)/ (tpctil-tcum)
      f1 = f
      x1 = x
  500 continue
      write(lout,*) ' FUNLUX:  WARNING. FUNPCT fails trapezoid.'
!         END OF TRAPEZOID LOOP
!         Adjust interval using Gaussian integration with
!             Newton corrections since F is the derivative
  520 continue
      x1 = xfcum(ibin)
      xbest = x
      dtbest = tpctil
      tpart = tpctil
!         Allow for maximum NITMAX more iterations on RADAPT
      do 550 ihome= 1, nitmax
  535 xincr = (tpctil-tpart) / max(f,fmin)
      x = xbest + xincr
      x2 = x
        if (ihome .gt. 1 .and. x2 .eq. xbest) then
        write(lout,'(A,G12.3)') ' FUNLUX: WARNING from FUNPCT: insufficient precision at X=',x
        go to 580
        endif
      refx = abs(x)+precis
      call radapt(func,x1,x2,1,rteps,zero,tpart2,uncert)
      dtpar2 = tpart2-tpctil
      dtabs = abs(dtpar2)
      if(abs(xincr)/refx .lt. precis) goto 545
      if(dtabs .lt. dtbest) goto 545
      xincr = xincr * 0.5
      goto 535
  545 dtbest = dtabs
      xbest = x
      tpart = tpart2
      f = func(x)
      if(f .lt. 0.) goto 900
      if(dtabs .lt. rteps*tpctil) goto 580
  550 continue
      write(lout,'(A,I4)') ' FUNLUX: WARNING from FUNPCT: cannot converge, bin',ibin

  580 continue
      xincr = (tpctil-tpart) / max(f,fmin)
      x = xbest + xincr
      xfcum(ibin+1) = x
      f = func(x)
      if(f .lt. 0.) goto 900
  600 continue
!         END OF LOOP OVER BINS
      x1 = xfcum((nlo+nbins)-1)
      x2 = xhigh
      call radapt(func,x1,x2,1,rteps,zero,tpart ,uncert)
      aberr = abs(tpart-tpctil)/tftot
!      WRITE(6,1001) IFUNC,XLOW,XHIGH
      if(aberr .gt. rteps)  write(lout,1002) aberr
      return
  900 write(lout,1000) x,f
      ierr = 1
      return
 1000 format(/' FUNLUX fatal error in FUNPCT: function negative:'/      &
&,' at X=',e15.6,', F=',e15.6/)
! 1001 FORMAT(' FUNPCT has prepared USER FUNCTION',I3,
!     + ' between',G12.3,' and',G12.3,' for FUNLUX.')
 1002 format(' WARNING: Relative error in cumulative distribution may be as big as',f10.7)

end subroutine funpct

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine funlux(array,xran,len)
!         Generation of LEN random numbers in any given distribution,
!         by 4-point interpolation in the inverse cumulative distr.
!         which was previously generated by FUNLXP
  use mod_ranlux
  implicit none
      integer len,ibuf,j,j1
      real(kind=fPrec) array,xran,gap,gapinv,tleft,bright,gaps,gapins,x,p,a,b
      dimension array(200)
      dimension xran(len)
!  Bin width for main sequence, and its inverse
      parameter (gap= 1./99.,  gapinv=99.)
!  Top of left tail, bottom of right tail (each tail replaces 2 bins)
      parameter (tleft= 2./99.,bright=97./99.)
!  Bin width for minor sequences (tails), and its inverse
      parameter (gaps=tleft/49.,  gapins=1./gaps)
!
!   The array ARRAY is assumed to have the following structure:
!        ARRAY(1-100) contains the 99 bins of the inverse cumulative
!                     distribution of the entire function.
!        ARRAY(101-150) contains the 49-bin blowup of main bins
!                       1 and 2 (left tail of distribution)
!        ARRAY(151-200) contains the 49-bin blowup of main bins
!                       98 and 99 (right tail of distribution)
!
      x = zero ! -Wmaybe-uninitialized
      call ranlux(xran,len)
!      call ranecu(xran,len,-1)

      do 500 ibuf= 1, len
      x = xran(ibuf)
      j = int(  x    *gapinv) + 1
      if (j .lt. 3)  then
         j1 = int( x *gapins)
             j = j1 + 101
             j = max(j,102)
             j = min(j,148)
         p = (   x -gaps*real(j1-1)) * gapins
         a = (p+1.0) * array(j+2) - (p-2.0)*array(j-1)
         b = (p-1.0) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.0))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5
      else if (j .gt. 97)  then
         j1 = int((x-bright)*gapins)
             j = j1 + 151
             j = max(j,152)
             j = min(j,198)
         p = ((x -bright) -gaps*(j1-1)) * gapins
         a = (p+1.0) * array(j+2) - (p-2.0)*array(j-1)
         b = (p-1.0) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.0))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5
      else
!      J = MAX(J,2)
!      J = MIN(J,98)
         p = (   x -gap*real(j-1)) * gapinv
         a = (p+1.) * array(j+2) - (p-2.)*array(j-1)
         b = (p-1.) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5
      endif
  500 continue
      tftot = x
      return
end subroutine funlux

subroutine funlz(func,x2low,x2high,xlow,xhigh)
! FIND RANGE WHERE FUNC IS NON-ZERO.
! WRITTEN 1980, F. JAMES
! MODIFIED, NOV. 1985, TO FIX BUG AND GENERALIZE
! TO FIND SIMPLY-CONNECTED NON-ZERO REGION (XLOW,XHIGH)
! ANYWHERE WITHIN THE GIVEN REGION (X2LOW,H2HIGH).
!    WHERE 'ANYWHERE' MEANS EITHER AT THE LOWER OR UPPER
!    EDGE OF THE GIVEN REGION, OR, IF IN THE MIDDLE,
!    COVERING AT LEAST 1% OF THE GIVEN REGION.
! OTHERWISE IT IS NOT GUARANTEED TO FIND THE NON-ZERO REGION.
! IF FUNCTION EVERYWHERE ZERO, FUNLZ SETS XLOW=XHIGH=0.
      use crcoall
      implicit none
      external func
      integer logn,nslice,i,k
      real(kind=fPrec) xhigh,xlow,x2high,x2low,func,xmid,xh,xl,xnew
      xlow = x2low
      xhigh = x2high
!         FIND OUT IF FUNCTION IS ZERO AT ONE END OR BOTH
      xmid = xlow
      if (func(xlow) .gt. 0.) go to 120
      xmid = xhigh
      if (func(xhigh) .gt. 0.)  go to 50
!         FUNCTION IS ZERO AT BOTH ENDS,
!         LOOK FOR PLACE WHERE IT IS NON-ZERO.
      do 30 logn= 1, 7
      nslice = 2**logn
      do 20 i= 1, nslice, 2
      xmid = xlow + (real(i) * (xhigh-xlow)) / real(nslice)
      if (func(xmid) .gt. 0.)  go to 50
   20 continue
   30 continue
!         FALLING THROUGH LOOP MEANS CANNOT FIND NON-ZERO VALUE
      write(lout,554)
      write(lout,555) xlow, xhigh
      xlow = 0.
      xhigh = 0.
      go to 220
!
   50 continue
!         DELETE 'LEADING' ZERO RANGE
      xh = xmid
      xl = xlow
      do 70 k= 1, 20
      xnew = 0.5*(xh+xl)
      if (func(xnew) .eq. 0.) go to 68
      xh = xnew
      go to 70
   68 xl = xnew
   70 continue
      xlow = xl
      write(lout,555) x2low,xlow
  120 continue
      if (func(xhigh) .gt. 0.) go to 220
!         DELETE 'TRAILING' RANGE OF ZEROES
      xl = xmid
      xh = xhigh
      do 170 k= 1, 20
      xnew = 0.5*(xh+xl)
      if (func(xnew) .eq. 0.) go to 168
      xl = xnew
      go to 170
  168 xh = xnew
  170 continue
      xhigh = xh
      write(lout,555) xhigh, x2high
!
  220 continue
      return
  554 format('0CANNOT FIND NON-ZERO FUNCTION VALUE')
  555 format(' FUNCTION IS ZERO FROM X=',e12.5,' TO ',e12.5)
end subroutine funlz

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine radapt(f,a,b,nseg,reltol,abstol,res,err)

! RES = Estimated Integral of F from A to B,
! ERR = Estimated absolute error on RES.
! NSEG  specifies how the adaptation is to be done:
!    =0   means use previous binning,
!    =1   means fully automatic, adapt until tolerance attained.
!    =n>1 means first split interval into n equal segments,
!         then adapt as necessary to attain tolerance.
! The specified tolerances are:
!        relative: RELTOL ;  absolute: ABSTOL.
!    It stop s when one OR the other is satisfied, or number of
!    segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
!    can be set to zero, in which case only the other is used.

      implicit none

      external f
      integer nseg,ndim,nter,nsegd,i,iter,ibig
      real(kind=fPrec) err,res,abstol,reltol,b,a,xlo,xhi,tval,ters,te,root,xhib,bin,xlob,bige,hf,xnew,r1,f
      real(kind=fPrec) tvals,terss

      parameter (ndim=100)
      parameter (r1 = 1., hf = r1/2.)

      dimension xlo(ndim),xhi(ndim),tval(ndim),ters(ndim)
      save xlo,xhi,tval,ters,nter
      data nter /0/

      if(nseg .le. 0)  then
       if(nter .eq. 0) then
        nsegd=1
        go to 2
       endif
       tvals=zero
       terss=zero
       do 1 i = 1,nter
       call rgs56p(f,xlo(i),xhi(i),tval(i),te)
       ters(i)=te**2
       tvals=tvals+tval(i)
       terss=terss+ters(i)
    1  continue
       root= sqrt(two*terss)
       go to 9
      endif
      nsegd=min(nseg,ndim)
    2 xhib=a
      bin=(b-a)/real(nsegd,fPrec)
      do 3 i = 1,nsegd
      xlo(i)=xhib
      xlob=xlo(i)
      xhi(i)=xhib+bin
      if(i .eq. nsegd) xhi(i)=b
      xhib=xhi(i)
      call rgs56p(f,xlob,xhib,tval(i),te)
      ters(i)=te**2
    3 continue
      nter=nsegd
      do 4 iter = 1,ndim
      tvals=tval(1)
      terss=ters(1)
      do 5 i = 2,nter
      tvals=tvals+tval(i)
      terss=terss+ters(i)
    5 continue
      root=sqrt(two*terss)

      if(root .le. abstol .or. root .le. reltol*abs(tvals)) then
        goto 9
      end if

      if(nter .eq. ndim) go to 9
      bige=ters(1)
      ibig=1
      do 6 i = 2,nter
      if(ters(i) .gt. bige) then
       bige=ters(i)
       ibig=i
      endif
    6 continue
      nter=nter+1
      xhi(nter)=xhi(ibig)
      xnew=hf*(xlo(ibig)+xhi(ibig))
      xhi(ibig)=xnew
      xlo(nter)=xnew
      call rgs56p(f,xlo(ibig),xhi(ibig),tval(ibig),te)
      ters(ibig)=te**2
      call rgs56p(f,xlo(nter),xhi(nter),tval(nter),te)
      ters(nter)=te**2
    4 continue
    9 res=tvals
      err=root
      return
end subroutine radapt

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine rgs56p(f,a,b,res,err)

  implicit none

  integer i
  real(kind=fPrec) err,res,b,a,f,w6,x6,w5,x5,rang,r1,hf
  real(kind=fPrec) e5,e6

  parameter (r1 = 1., hf = r1/2.)
  dimension x5(5),w5(5),x6(6),w6(6)

  data (x5(i),w5(i),i=1,5)                                          &
 &/4.6910077030668004e-02, 1.1846344252809454e-01,                  &
 &2.3076534494715846e-01, 2.3931433524968324e-01,                   &
 &5.0000000000000000e-01, 2.8444444444444444e-01,                   &
 &7.6923465505284154e-01, 2.3931433524968324e-01,                   &
 &9.5308992296933200e-01, 1.1846344252809454e-01/

  data (x6(i),w6(i),i=1,6)                                          &
 &/3.3765242898423989e-02, 8.5662246189585178e-02,                  &
 &1.6939530676686775e-01, 1.8038078652406930e-01,                   &
 &3.8069040695840155e-01, 2.3395696728634552e-01,                   &
 &6.1930959304159845e-01, 2.3395696728634552e-01,                   &
 &8.3060469323313225e-01, 1.8038078652406930e-01,                   &
 &9.6623475710157601e-01, 8.5662246189585178e-02/

  rang=b-a
  e5=zero
  e6=zero
  do i = 1,5
    e5=e5+dble(w5(i)*f(a+rang*x5(i)))
    e6=e6+dble(w6(i)*f(a+rang*x6(i)))
  end do

  e6=e6+dble(w6(6)*f(a+rang*x6(6)))
  res=real((dble(hf)*(e6+e5))*dble(rang))
  err=real(abs((e6-e5)*dble(rang)))
  return
end subroutine rgs56p

end module mod_funlux
