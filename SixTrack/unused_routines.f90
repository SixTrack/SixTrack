! The subroutines in this file were moved from sixtrack.s90, but are never called from any current code.

      subroutine rinv(n,a,idim,ir,ifail)
!-----------------------------------------------------------------------
!
!     ******************************************************************
!
!     REPLACES A BY ITS INVERSE.
!
!     (PARAMETERS AS FOR REQINV.)
!
!     CALLS ... RFACT, RFINV, F010PR, ABEND.
!
!     ******************************************************************
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      implicit none
      integer idim,ifail,ir,jfail,k,kprnt,n
      real(kind=fPrec) t1,t2,t3,a,det,temp,s,                           &
     &c11,c12,c13,c21,c22,c23,c31,c32,c33

      character(len=6) name
      dimension ir(n),a(idim,n)
      data name/'RINV'/,kprnt/0/
      save
!-----------------------------------------------------------------------
!
!  TEST FOR PARAMETER ERRORS.
!
      if((n.lt.1).or.(n.gt.idim)) goto 7
!
!  TEST FOR N.LE.3.
!
      if(n.gt.3) goto 6
      ifail=0
      if(n.lt.3) goto 4
!
!  N=3 CASE.
!
!     COMPUTE COFACTORS.
      c11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
      c12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
      c13=a(2,1)*a(3,2)-a(2,2)*a(3,1)
      c21=a(3,2)*a(1,3)-a(3,3)*a(1,2)
      c22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
      c23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
      c31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
      c32=a(1,3)*a(2,1)-a(1,1)*a(2,3)
      c33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      t1=abs(a(1,1))
      t2=abs(a(2,1))
      t3=abs(a(3,1))
!
!     (SET TEMP=PIVOT AND DET=PIVOT*DET.)
      if(t1.ge.t2) goto 1
         if(t3.ge.t2) goto 2
!        (PIVOT IS A21)
            temp=a(2,1)
            det=c13*c32-c12*c33
            goto 3
    1 if(t3.ge.t1) goto 2
!     (PIVOT IS A11)
         temp=a(1,1)
         det=c22*c33-c23*c32
         goto 3
!     (PIVOT IS A31)
    2    temp=a(3,1)
         det=c23*c12-c22*c13
!
!     SET ELEMENTS OF INVERSE IN A.
    3 if(det.eq.zero) goto 8
      s=temp/det
      a(1,1)=s*c11
      a(1,2)=s*c21
      a(1,3)=s*c31
      a(2,1)=s*c12
      a(2,2)=s*c22
      a(2,3)=s*c32
      a(3,1)=s*c13
      a(3,2)=s*c23
      a(3,3)=s*c33
      return
!
    4 if(n.lt.2) goto 5
!
!  N=2 CASE BY CRAMERS RULE.
!
      det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if(det.eq.zero) goto 8
      s=one/det                                                          !hr07
      c11   =s*a(2,2)
      a(1,2)=(-one*s)*a(1,2)                                             !hr07
      a(2,1)=(-one*s)*a(2,1)                                             !hr07
      a(2,2)=s*a(1,1)
      a(1,1)=c11
      return
!
!  N=1 CASE.
!
    5 if(a(1,1).eq.zero) goto 8
      a(1,1)=one/a(1,1)                                                  !hr07
      return
!
!  N.GT.3 CASES.  FACTORIZE MATRIX AND INVERT.
!
    6 call rfact(n,a,idim,ir,ifail,det,jfail)
      if(ifail.ne.0) return
      call rfinv(n,a,idim,ir)
      return
!
!  ERROR EXITS.
!
    7 ifail=+1
      call f010pr(name,n,idim,k,kprnt)
      return
!
    8 ifail=-1
      return
!
      end

      subroutine rfeqn(n,a,idim,ir,k,b)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants

      implicit none

      integer i,idim,ij,im1,ir,j,k,l,m,n,nm1,nmi,nmjp1,nxch
      real(kind=fPrec) a,b,te,x,y
      real(kind=fPrec) s21,s22
      character(len=6) hname
      dimension ir(*),a(idim,*),b(idim,*)
      data      hname               /  ' RFEQN'  /
      save
!-----------------------------------------------------------------------
      if(idim .ge. n  .and.  n .gt. 0  .and.  k .gt. 0)  goto 210
      call tmprnt(hname,n,idim,k)
      return
 210  nxch  =  ir(n)
      if(nxch .eq. 0)  goto 220
      do 212    m  =  1, nxch
         ij  =  ir(m)
         i   =  ij / 4096
         j   =  mod(ij,4096)
         do 211   l  =  1, k
            te      =  b(i,l)
            b(i,l)  =  b(j,l)
            b(j,l)  =  te
 211        continue
 212     continue
 220  do 221    l  =  1, k
         b(1,l)  =  a(1,1)*b(1,l)
 221     continue
      if(n .eq. 1)  goto 299
      do 243    l  =  1, k
         do 232   i  =  2, n
            im1  =  i-1
            s21  =  -one*b(i,l)                                          !hr07
            do 231   j  =  1, im1
               s21  =  dotf(a(i,j),b(j,l),s21)
 231           continue
            b(i,l)  = (-one*a(i,i))*s21                                  !hr07
 232        continue
         nm1  =  n-1
         do 242   i  =  1, nm1
            nmi  =  n-i
            s22  =  -one*b(nmi,l)                                        !hr07
            do 241   j  =  1, i
               nmjp1  =  n - j+1
               s22    =  dotf(a(nmi,nmjp1),b(nmjp1,l),s22)
 241           continue
            b(nmi,l) = -one*s22
 242        continue
 243     continue
 299  continue
      return

      contains

        real(kind=fPrec) function dotf(x,y,s21)
          real(kind=fPrec), intent(in) :: x, y, s21
          dotf = x * y + s21
        end function dotf

      end subroutine rfeqn

      subroutine dfeqn(n,a,idim,ir,k,b)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none

      integer i,idim,ij,im1,ir,j,k,l,m,n,nm1,nmi,nmjp1,nxch
      real(kind=fPrec) a,b,x,y,te
      real(kind=fPrec) s21,s22
      character(len=6) hname
      dimension ir(*),a(idim,*),b(idim,*)
      data      hname               /  ' DFEQN'  /
      save
!-----------------------------------------------------------------------
      if(idim .ge. n  .and.  n .gt. 0  .and.  k .gt. 0)  goto 210
      call tmprnt(hname,n,idim,k)
      return
 210  nxch  =  ir(n)
      if(nxch .eq. 0)  goto 220
      do 212    m  =  1, nxch
         ij  =  ir(m)
         i   =  ij / 4096
         j   =  mod(ij,4096)
         do 211   l  =  1, k
            te      =  b(i,l)
            b(i,l)  =  b(j,l)
            b(j,l)  =  te
 211        continue
 212     continue
 220  do 221    l  =  1, k
         b(1,l)  =  a(1,1)*b(1,l)
 221     continue
      if(n .eq. 1)  goto 299
      do 243    l  =  1, k
         do 232   i  =  2, n
            im1  =  i-1
            s21  =  - b(i,l)
            do 231   j  =  1, im1
               s21  =  dotf(a(i,j),b(j,l),s21)
 231           continue
            b(i,l)  =  (-one*a(i,i))*s21                                 !hr07
 232        continue
         nm1  =  n-1
         do 242   i  =  1, nm1
            nmi  =  n-i
            s22  =  -one*b(nmi,l)                                        !hr07
            do 241   j  =  1, i
               nmjp1  =  n - j+1
               s22    =  dotf(a(nmi,nmjp1),b(nmjp1,l),s22)
 241           continue
            b(nmi,l)  =  -one*s22                                        !hr07
 242        continue
 243     continue
 299  continue
      return

      contains

        real(kind=fPrec) function dotf(x,y,s21)
          real(kind=fPrec), intent(in) :: x, y, s21
          dotf = x * y + s21
        end function dotf

      end subroutine dfeqn
