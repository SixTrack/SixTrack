

subroutine lieinit(no1,nv1,nd1,ndpt1,iref1,nis)
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use mod_lie_dab, only : ndim,ndim2,ntt,nreso,nd,nd2,no,nv,ifilt,idpr,iref,itu,lienot,iflow,   &
        jtune,ndc,ndc2,ndpt,ndt,nplane,ista,idsta,mx,nres,epsplane,xplane,sta,dsta,angle,ps,    &
        rads,xintex,iscrda

      implicit none
      integer i,iref1,nd1,ndc1,ndpt1,nis,no1,nv1
      real(kind=fPrec) ang,ra,st
!! Lieinit initializes AD Package and Lielib
      dimension st(ndim),ang(ndim),ra(ndim)
      call daexter

      do i=1,ndim
        nplane(i)=2*i-1
        ang(i)=zero
        ra(i)=zero
        st(i)=one
      end do

      no=no1
      nv=nv1
      nd=nd1
      nd2=2*nd1

      do i=1,100
        iscrda(i)=0
      end do

      call daini(no,nv,0)

      if(nis.gt.0)call etallnom(iscrda,nis,'$$IS      ')
       if(ndpt1.eq.0) then
         ndpt=0
         ndt=0
         ndc1=0
       else
         ndpt=ndpt1
         ndc1=1
            if(ndpt.eq.nd2) then
              ndt=nd2-1
            else
              ndt=nd2
              if(ndpt.ne.nd2-1) then
                write(lout,*) ' LETHAL ERROR IN LIEINIT'
                call prror
              endif
            endif
       endif

      ndc=ndc1
      ndc2=2*ndc1
      iref=0
      call initpert(st,ang,ra)
      iref=iref1
        if(iref1.eq.0) then
         itu=0
        else
         itu=1
        endif
      if(iref1.eq.0) iref=-1

      if(idpr.eq.1)write(lout,*) ' NO = ',no,' IN DA-CALCULATIONS '

      do i=0,20
        xintex(i)=zero
      enddo

      xintex(          0)=       1.000000000000000_fPrec                     !hr11
      xintex(          1)=  5.000000000000000e-001_fPrec                       !hr11
      xintex(          2)=  8.333333333333334e-002_fPrec                       !hr11
      xintex(          3)=  0.000000000000000e+000_fPrec                       !hr11
      xintex(          4)= -1.388888888888898e-003_fPrec                       !hr11
      xintex(          5)=  0.000000000000000e+000_fPrec                       !hr11
      xintex(          6)=  3.306878306878064e-005_fPrec                       !hr11
      xintex(          7)= zero
      xintex(          8)= -8.267195767165669e-007_fPrec                       !hr11
      xintex(          9)=  zero
      xintex(         10)=  4.592886537931051e-008_fPrec                       !hr11
      return
end subroutine lieinit

subroutine flowpara(ifl,jtu)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : iflow,jtune
      implicit none
      integer ifl,jtu
      iflow=ifl
      jtune=jtu
      return
end subroutine flowpara

subroutine pertpeek(st,ang,ra)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,ndim,angle,radn,sta
      implicit none
      integer i
      real(kind=fPrec) ang,ra,st
      dimension st(ndim),ang(ndim),ra(ndim)
      do i=1,nd
        st(i)=sta(i)
        ang(i)=angle(i)
        ra(i)=radn(i)
      end do
      return
end subroutine pertpeek

subroutine inputres(mx1,nres1)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : ndim,nreso,mx,nres
      implicit none
      integer i,j
      integer mx1(ndim,nreso),nres1

      nres=nres1
      do i=1,nreso
        do j=1,ndim
          mx(j,i)=0
        end do
      end do

      do i=1,nres
        do j=1,ndim
          mx(j,i)=mx1(j,i)
        end do
      end do

      return
end subroutine inputres

subroutine respoke(mres,nre,ire)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,iref,ndc,ndc2,ndpt,ndt,ndim,nreso,idsta,ista,angle,dsta,sta,mx,nres
      implicit none
      integer i,ire,j,nre
      real(kind=fPrec) ang,ra,st
      integer mres(ndim,nreso)
      dimension ang(ndim),ra(ndim),st(ndim)
      iref=ire
      nres=nre
      do j=1,nreso
        do i=1,nd
          mx(i,j)=mres(i,j)
        end do
      end do

      call initpert(st,ang,ra)
      return
end subroutine respoke

subroutine liepeek(iia,icoast)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,no,nv,ndc,ndc2,ndpt,ndt
      implicit none
      integer iia(*),icoast(*)

      iia(1)=no
      iia(2)=nv
      iia(3)=nd
      iia(4)=nd2
      icoast(1)=ndc
      icoast(2)=ndc2
      icoast(3)=ndt
      icoast(4)=ndpt

      return
end subroutine liepeek

subroutine lienot(not)
      use floatPrecision
      use mathlib_bouncer
      implicit none
      integer no,not

      call danot(not)
      no=not

      return
end subroutine lienot

subroutine etallnom(x,n,nom)
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      implicit none
      integer i,n,nd2
! CREATES A AD-VARIABLE WHICH CAN BE DESTROYED BY DADAL
! allocates vector of n polynomials and give it the name NOM=A10
      integer x(*),i1(4),i2(4)
      character(len=10) nom

      do  i=1,abs(n)                                                    !hr11
        x(i)=0
      end do

      call daallno(x,abs(n),nom)                                         !hr11

      if(n.lt.0) then
        call liepeek(i1,i2)
        nd2=i1(4)
        do i=nd2+1,-n
          call davar(x(i),zero,i)
        end do
      endif

      return
end subroutine etallnom

subroutine etall(x,n)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
      integer i,n,nd2
! allocates vector of n polynomials
      integer x(*),i1(4),i2(4)

      do i=1,iabs(n)
        x(i)=0
      end do

      call daallno(x,iabs(n),'ETALL     ')
         if(n.lt.0) then
          call liepeek(i1,i2)
          nd2=i1(4)
           do i=nd2+1,-n
             call davar(x(i),zero,i)
           end do
         endif
      return
end subroutine etall

subroutine etall1(x)
      use floatPrecision
      use mathlib_bouncer
      implicit none
      integer x
      x=0
      call daallno(x,1,'ETALL     ')
      return
end subroutine etall1

subroutine dadal1(x)
      use floatPrecision
      use mathlib_bouncer
      implicit none
      integer x
      call dadal(x,1)
      return
end subroutine dadal1

subroutine etppulnv(x,xi,xff)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i
      integer  x(*)
      real(kind=fPrec) xi(*),xff(*),xf(ntt),xii(ntt)

      do i=1,nv
        xii(i)=xi(i)
      end do

      do i=nv+1,ntt
        xii(i)=zero
      end do

      call ppush(x,nv,xii,xf)

      do i=1,nv
        xff(i)=xf(i)
      end do

      return
end subroutine etppulnv

subroutine etmtree(y,x)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i,ie,iv,nt
! ROUTINES USING THE MAP IN AD-FORM
      dimension ie(ntt),iv(ntt)
      integer  x(*),y(*)

      nt=nv-nd2
      if(nt.gt.0) then
        call etallnom(ie,nt,'IE        ')
        do i=nd2+1,nv
          call davar(ie(i-nd2),zero,i)
        end do

        do i=nd2+1,nv
          iv(i)=ie(i-nd2)
        end do
      endif

      do i=1,nd2
        iv(i)=y(i)
      end do

      call mtree(iv,nv,x,nv)

      if(nt.gt.0) then
        call dadal(ie,nt)
      endif

      return
end subroutine etmtree

subroutine etppush(x,xi)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i
      integer  x(*)
      real(kind=fPrec) xi(*),xf(ntt),xii(ntt)

      do i=1,nd2
        xii(i)=xi(i)
      end do

      call ppush(x,nv,xii,xf)

      do i=1,nd2
        xi(i)=xf(i)
      end do

      return
end subroutine etppush

subroutine etppush2(x,xi,xff)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i
      integer  x(*)
      real(kind=fPrec) xi(*),xff(*),xf(ntt),xii(ntt)

      do i=1,nd2
        xii(i)=xi(i)
      end do

      call ppush(x,nv,xii,xf)

      do i=1,nd2
      xff(i)=xf(i)
      enddo

      return
end subroutine etppush2

subroutine ppushlnv(x,xi,xff,nd1)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i,nd1
      integer  x(*)
      real(kind=fPrec) xi(*),xff(*),xf(ntt),xii(ntt)

      do i=1,nd1
        xii(i)=xi(i)
      end do

      do i=nd1+1,ntt
        xii(i)=zero
      end do

      call ppush(x,nv,xii,xf)

      do i=1,nd1
        xff(i)=xf(i)
      end do

      return
end subroutine ppushlnv

subroutine etcct(x,y,z)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i,ie,iv,nt
!  Z=XoY
      dimension ie(ntt),iv(ntt)
      integer  x(*),y(*),z(*)

      nt=nv-nd2
      if(nt.gt.0) then
        call etallnom(ie,nt,'IE        ')
        do i=nd2+1,nv
          call davar(ie(i-nd2),zero,i)
        end do

        do i=nd2+1,nv
          iv(i)=ie(i-nd2)
        end do
      endif

      do i=1,nd2
       iv(i)=y(i)
      end do

      call dacct(x,nd2,iv,nv,z,nd2)

      if(nt.gt.0) then
        call dadal(ie,nt)
      endif

      return
end subroutine etcct

subroutine trx(h,rh,y)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i,ie,iv,nt
!  :RH: = Y :H: Y^-1 =  :HoY:
      dimension ie(ntt),iv(ntt)
      integer h,rh
      integer y(*)

      nt=nv-nd2
      if(nt.gt.0) then
        call etallnom(ie,nt,'IE        ')

        do i=nd2+1,nv
          call davar(ie(i-nd2),zero,i)
        end do

        do i=nd2+1,nv
          iv(i)=ie(i-nd2)
        end do
      endif

      do i=1,nd2
        iv(i)=y(i)
      end do

      call dacct(h,1,iv,nv,rh,1)

      if(nt.gt.0) then
        call dadal(ie,nt)
      endif

      return
end subroutine trx

subroutine trxflo(h,rh,y)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndim2,ntt
      implicit none
      integer j,k
!  *RH* = Y *H* Y^-1  CHANGE OF A VECTOR FLOW OPERATOR
      integer h(*),rh(*),y(*)
      integer yi(ndim2),ht(ndim2),b1(1),b2(1)


      call etallnom(yi,nd2  ,'YI        ')
      call etallnom(ht,nd2  ,'HT        ')
      call etallnom(b1(1),1,'B1        ')
      call etallnom(b2(1),1,'B2        ')

      call etinv(y,yi)
!----- HT= H o Y
      call etcct(h,y,ht)
!----
      call daclrd(rh)

      do j=1,nd2
        do k=1,nd2
          call dader(k,yi(j),b1(1))
          call trx(b1(1),b2(1),y)
          call damul(b2,ht(k),b1(1))
          call daadd(b1(1),rh(j),b2(1))
          call dacop(b2(1),rh(j))
        end do
      end do

      call dadal(b2(1),1)
      call dadal(b1(1),1)
      call dadal(ht,nd2)
      call dadal(yi,nd2)
      return
end subroutine trxflo

subroutine simil(a,x,ai,y)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
!  Y= AoXoAI
      integer  x(*),y(*),a(*),ai(*)

      integer w(ndim2),v(ndim2)
!
      call etallnom(w,nd2  ,'W         ')
      call etallnom(v,nd2  ,'V         ')

      call etcct(a,x,w)
      call etcct(w,ai,v)

      call dacopd(v,y)

      call dadal(v,nd2)
      call dadal(w,nd2)
      return
end subroutine simil

subroutine etini(x)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,no,nv,ndim2
      implicit none
      integer i
!  X=IDENTITY
      integer x(*)
!*DAEXT(NO,NV) X(NDIM2)
      do i=1,nd2
        call davar(x(i),zero,i)
      end do

      return
end subroutine etini

subroutine etinv(x,y)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i,ie1,ie2,iv1,iv2,nt
! Y=X^-1
      dimension ie1(ntt),ie2(ntt),iv1(ntt),iv2(ntt)
      integer x(*),y(*)

      nt=nv-nd2
      if(nt.gt.0) then
       do i=1,nt
         ie1(i)=0
         ie2(i)=0
       end do

       call etallnom(ie1,nt,'IE1       ')
       call etallnom(ie2,nt,'IE2       ')

       do i=nd2+1,nv
         call davar(ie1(i-nd2),zero,i)
       end do

       do i=nd2+1,nv
         iv1(i)=ie1(i-nd2)
         iv2(i)=ie2(i-nd2)
       end do
      endif

      do i=1,nd2
        iv1(i)=x(i)
        iv2(i)=y(i)
      end do

      call dainv(iv1,nv,iv2,nv)

      if(nt.gt.0) then
        call dadal(ie2,nt)
        call dadal(ie1,nt)
      endif

      return
end subroutine etinv

subroutine etpin(x,y,jj)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i,ie1,ie2,iv1,iv2,jj,nt
!  Y=PARTIAL INVERSION OF X SEE BERZ'S PACKAGE
      dimension ie1(ntt),ie2(ntt),iv1(ntt),iv2(ntt),jj(*)

      integer x(*),y(*)

      nt=nv-nd2
      if(nt.gt.0) then
       do i=1,nt
         ie1(i)=0
         ie2(i)=0
       end do

       call etallnom(ie1,nt,'IE1       ')
       call etallnom(ie2,nt,'IE2       ')

       do i=nd2+1,nv
         call davar(ie1(i-nd2),zero,i)
       end do

       do i=nd2+1,nv
         iv1(i)=ie1(i-nd2)
         iv2(i)=ie2(i-nd2)
       end do
      endif

      do i=1,nd2
        iv1(i)=x(i)
        iv2(i)=y(i)
      end do

      call dapin(iv1,nv,iv2,nv,jj)

      if(nt.gt.0) then
        call dadal(ie2,nt)
        call dadal(ie1,nt)
      endif

      return
end subroutine etpin

subroutine dapek0(v,x,jj)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : ntt
      implicit none
      integer i,jj
      real(kind=fPrec) x
!- MORE EXTENSIONS OF BASIC BERZ'S PACKAGE
      integer v(*),jd(ntt)
      dimension x(*)

      do i=1,ntt
        jd(i)=0
      end do

      do i=1,jj
        call dapek(v(i),jd,x(i))
      end do

      return
end subroutine dapek0

subroutine dapok0(v,x,jj)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : ntt
      implicit none
      integer i,jj
      real(kind=fPrec) x
      integer v(*),jd(ntt)
      dimension x(*)

      do i=1,ntt
        jd(i)=0
      end do

      do i=1,jj
        call dapok(v(i),jd,x(i))
      end do

      return
end subroutine dapok0

subroutine dapokzer(v,jj)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : ntt
      implicit none
      integer i,jj
      integer v(*),jd(ntt)

      do i=1,ntt
        jd(i)=0
      end do

      do i=1,jj
        call dapok(v(i),jd,zero)
      end do

      return
end subroutine dapokzer

subroutine davar0(v,x,jj)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : ntt
      implicit none
      integer i,jj
      real(kind=fPrec) x
      integer v(*)
      dimension x(*)

      do i=1,jj
        call davar(v(i),x(i),i)
      end do

      return
end subroutine davar0

subroutine comcfu(b,f1,f2,c)
      use floatPrecision
      use mathlib_bouncer
      implicit none
      real(kind=fPrec) f1,f2
      external f1,f2
! Complex dacfu
      integer b(*),c(*),t(4)
      call etall(t,4)

      call dacfu(b(1),f1,t(1))
      call dacfu(b(1),f2,t(2))
      call dacfu(b(2),f1,t(3))
      call dacfu(b(2),f2,t(4))

      call dasub(t(1),t(4),c(1))
      call daadd(t(2),t(3),c(2))
      call dadal(t,4)
      return
end subroutine comcfu

subroutine take(h,m,ht)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,no,nv,ntt
      implicit none
      integer i,m
      real(kind=fPrec) r
!  HT= H_M  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
      integer h,ht,j(ntt)

      integer b1(1),b2(1),b3(1)

      call etallnom(b1(1),1,'B1        ')
      call etallnom(b2(1),1,'B2        ')
      call etallnom(b3(1),1,'B3        ')

      if(no.ge.2) then
       if(m.eq.0) then
         do i=1,ntt
           j(i)=0
         enddo
         call dapek(h,j,r)
         call dacon(ht,r)
       else
         call danot(m)
         call dacop(h,b1(1))
         call danot(m-1)
         call dacop(b1(1),b2(1))
         call danot(no)
         call dasub(b1(1),b2(1),b3(1))
         call dacop(b3(1),ht)
       endif
      else
        do i=1,ntt
          j(i)=0
        enddo

        if(m.eq.0) then
          call dapek(h,j,r)
          call dacon(ht,r)
        elseif(m.eq.1)  then
          do i=1,nv
            j(i)=1
            call dapek(h,j,r)
            call dapok(b3(1),j,r)
            j(i)=0
          enddo
          call dacop(b3(1),ht)
        else
          call daclr(ht)
        endif
      endif

      call dadal(b3(1),1)
      call dadal(b2(1),1)
      call dadal(b1(1),1)
      return
end subroutine take

subroutine taked(h,m,ht)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndim2,ntt
      implicit none
      integer i,m
!  \VEC{HT}= \VEC{H_M}  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
      integer h(*),ht(*),j(ntt)

      integer b1(1),b2(1),x(ndim2)

      call etallnom(b1(1),1,'B1        ')
      call etallnom(b2(1),1,'B2        ')
      call etallnom(x,nd2  ,'X         ')

      do i=1,ntt
        j(i)=0
      enddo

      do i=1,nd2
        call take(h(i),m,ht(i))
      enddo

      call dadal(x,nd2)
      call dadal(b2(1),1)
      call dadal(b1(1),1)
      return
end subroutine taked

subroutine daclrd(h)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd2
      implicit none
      integer i
! clear a map : a vector of nd2 polynomials
      integer h(*)
      do i=1,nd2
        call daclr(h(i))
      end do
      return
end subroutine daclrd

subroutine dacopd(h,ht)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd2
      implicit none
      integer i
!    H goes into HT  (nd2 array)
      integer h(*),ht(*)
      do i=1,nd2
        call dacop(h(i),ht(i))
      end do
      return
end subroutine dacopd

subroutine dacmud(h,sca,ht)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd2
      implicit none
      integer i
      real(kind=fPrec) sca
      integer h(*),ht(*)
      do i=1,nd2
        call dacmu(h(i),sca,ht(i))
      end do
      return
end subroutine dacmud

subroutine dalind(h,rh,ht,rt,hr)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer i
      real(kind=fPrec) rh,rt
      integer h(*),ht(*),hr(*)

      integer b(ndim2)
!
      call etallnom(b,nd2  ,'B         ')

      do i=1,nd2
        call dalin(h(i),rh,ht(i),rt,b(i))
      end do
      call dacopd(b,hr)
      call dadal(b,nd2)
      return
end subroutine dalind

subroutine daread(h,nd1,mfile,xipo)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : ntt
      implicit none
      integer i,mfile,nd1
      real(kind=fPrec) rx,xipo
!  read a map
      integer h(*),j(ntt)

      do i=1,ntt
        j(i)=0
      end do

      do i=1,nd1
        call darea(h(i),mfile)
        call dapek(h(i),j,rx)
        rx=rx*xipo
        call dapok(h(i),j,rx)
      end do
      return
end subroutine daread

subroutine daprid(h,n1,n2,mfile)
      use floatPrecision
      implicit none
      integer i,mfile,n1,n2
!  print a map
      integer  h(*)
      if(mfile.le.0) return

      do i=n1,n2
        call dapri(h(i),mfile)
      end do

      return
end subroutine daprid

subroutine prresflo(h,eps,mfile)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ifilt,ndim2
      implicit none
      integer i,mfile
      real(kind=fPrec) deps,eps,filtres
!  print a map   in resonance basis for human consumption (useless)
      integer  h(*) ,b(ndim2) ,c(ndim2)
      external filtres
      call etall(b,nd2)
      call etall(c,nd2)
      call dacopd(h,c)
      do i=1,nd2
        ifilt=(-1)**i
        call  dacfu(c(i),filtres,h(i))
      enddo

      deps=-one
      call daeps(deps)
      call daeps(eps)

      call dacopd(c,h)
      call daeps(deps)
      call  dadal(c,nd2)
      call  dadal(b,nd2)
      return
end subroutine prresflo

real(kind=fPrec) function filtres(j)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ifilt,ndc,ndc2,ndpt,ndt
      implicit none
      integer i,ic
      integer j(*)
      filtres=one
      ic=0
      do i=1,(nd2-ndc2)
        ic=ic+j(i)*(-1)**(i+1)
      enddo

      ic=ic+ifilt

      if(ic.lt.0) filtres=zero

      if(ic.eq.0.and.ifilt.eq.1) then
        filtres=zero
      endif

      return
end function filtres

subroutine daflo(h,x,y)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2
      implicit none
      integer i
! LIE EXPONENT ROUTINES WITH FLOW OPERATORS

!     \VEC{H}.GRAD X =Y
      integer  h(*),x,y
      integer b1(1),b2(1),b3(1)
!
      call etallnom(b1(1),1,'B1        ')
      call etallnom(b2(1),1,'B2        ')
      call etallnom(b3(1),1,'B3        ')

      call daclr(b1(1))
      call daclr(b2(1))

      do i=1,nd2
        call dader(i,x,b2(1))
        call damul(b2(1),h(i),b3(1))
        call daadd(b3(1),b1(1),b2(1))
        call dacop(b2(1),b1(1))
      end do

      call dacop(b1(1),y)
      call dadal(b3(1),1)
      call dadal(b2(1),1)
      call dadal(b1(1),1)
      return
end subroutine daflo

subroutine daflod(h,x,y)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer i
      integer  h(*),x(*),y(*)
      integer b1(ndim2),b2(ndim2)

      call etall(b1,nd2)
      call etall(b2,nd2)

      call dacopd(h,b1)
      call dacopd(x,b2)

      do i=1,nd2
        call daflo(b1,b2(i),y(i))
      end do

      call dadal(b1,nd2)
      call dadal(b2,nd2)
      return
end subroutine daflod

subroutine intd(v,h,sca)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer i
      real(kind=fPrec) dlie,sca
! IF SCA=-1.D0
!     \VEC{V}.GRAD   = J GRAD H . GRAD = :H:

! IF SCA=1.D0
!     \VEC{V}.GRAD  = GRAD H . GRAD
      external dlie
      integer v(*),h

      integer b1(1),b2(1),b3(1),b4(1),x(ndim2)

      call etallnom(b1(1),1,'B1        ')
      call etallnom(b2(1),1,'B2        ')
      call etallnom(b3(1),1,'B3        ')
      call etallnom(b4(1),1,'B4        ')
      call etallnom(x,nd2  ,'X         ')

      call daclr(b4(1))
      call daclr(h)
      call etini(x)

      do i=1,nd
        call dacfu(v(2*i-1),dlie,b3(1))
        call dacfu(v(2*i),dlie,b1(1))
        call damul(b1(1),x(2*i-1),b2(1))
        call damul(b3(1),x(2*i),b1(1))
        call dalin(b2(1),one,b1(1),sca,b3(1))
        call daadd(b3(1),b4(1),b2(1))
        call dacop(b2(1),b4(1))
      end do

      call dacop(b4(1),h)
      call dadal(x,nd2)
      call dadal(b4(1),1)
      call dadal(b3(1),1)
      call dadal(b2(1),1)
      call dadal(b1(1),1)
      return
end subroutine intd

subroutine difd(h1,v,sca)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2
      implicit none
      integer i
      real(kind=fPrec) sca
! INVERSE OF INTD ROUTINE
      integer  v(*),h1
      integer b1(1),h(1)
      call etall(b1(1),1)
      call etall(h(1),1)
      call dacop(h1,h(1))

      do i=1,nd
        call dader(2*i-1,h(1),v(2*i))
        call dader(2*i,h(1),b1(1))
        call dacmu(b1(1),sca,v(2*i-1))
      end do

      call dadal(h(1),1)
      call dadal(b1(1),1)
      return
end subroutine difd

subroutine expflo(h,x,y,eps,nrmax)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : idpr
      implicit none
      integer i,nrmax
      real(kind=fPrec) coe,eps,r,rbefore
! DOES EXP( \VEC{H} ) X = Y
      integer h(*),x,y
      integer b1(1),b2(1),b3(1),b4(1)
      logical more

      call etallnom(b1(1),1,'B1        ')
      call etallnom(b2(1),1,'B2        ')
      call etallnom(b3(1),1,'B3        ')
      call etallnom(b4(1),1,'B4        ')

      call dacop(x,b4(1))
      call dacop(x,b1(1))
      more=.true.
      rbefore=1.0e30_fPrec

      do i=1,nrmax
      coe=one/real(i,fPrec)
      call dacmu(b1(1),coe,b2(1))
      call daflo(h,b2(1),b1(1))
      call daadd(b4(1),b1(1),b3(1))
      call daabs(b1(1),r)
         if(more) then
          if(r.gt.eps) then
       rbefore=r
       goto 100
          else
       rbefore=r
           more=.false.
          endif
         else
        if(r.ge.rbefore) then
        call dacop(b3(1),y)
        call dadal(b4(1),1)
        call dadal(b3(1),1)
        call dadal(b2(1),1)
        call dadal(b1(1),1)
        return
        endif
       rbefore=r
       endif
100   continue
      call dacop(b3(1),b4(1))
      end do

      if(idpr.ge.0) then
      write(lout,*) ' NORM ',eps,' NEVER REACHED IN EXPFLO '
      write(lout,*) 'NEW IDPR '
      read(5,*) idpr
      endif
      call dacop(b3(1),y)
      call dadal(b4(1),1)
      call dadal(b3(1),1)
      call dadal(b2(1),1)
      call dadal(b1(1),1)
      return
end subroutine expflo

subroutine expflod(h,x,w,eps,nrmax)
      use floatPrecision
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer j,nrmax
      real(kind=fPrec) eps
! DOES EXP( \VEC{H} ) \VEC{X} = \VEC{Y}
      integer x(*),w(*),h(*)
      integer b0(1),v(ndim2)
!
!
      call etallnom(b0(1),1,'B0        ')
      call etallnom(v,nd2  ,'V         ')

      call dacopd(x,v)

      do j=1,nd2
        call expflo(h,v(j),b0(1),eps,nrmax)
        call dacop(b0(1),v(j))
      end do

      call dacopd(v,w)
      call dadal(v,nd2)
      call dadal(b0(1),1)
      return
end subroutine expflod

subroutine facflo(h,x,w,nrmin,nrmax,sca,ifac)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer i,ifac,nmax,nrmax,nrmin
      real(kind=fPrec) eps,sca
! IFAC=1
! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX ) X= Y
! IFAC=-1
! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) X= Y
      integer x,w,h(*)
      integer bm(ndim2),b0(ndim2),v(1)
!
      call etallnom(bm,nd2  ,'BM        ')
      call etallnom(b0,nd2  ,'B0        ')
      call etallnom(v(1),1  ,'V         ')

      call dacop(x,v(1))

      eps=-one
      call daeps(eps)
      nmax=100
!
! IFAC =1 ---> V = EXP(:SCA*H(NRMAX):)...EXP(:SCA*H(NRMIN):)X
      if(ifac.eq.1) then
         do 1 i=nrmax,nrmin,-1
         call taked(h,i,b0)
         call dacmud(b0,sca,bm)

         call expflo(bm,v(1),b0(1),eps,nmax)
         call dacop(b0(1),v(1))
 1       continue
      else
! IFAC =-1 ---> V = EXP(:SCA*H(NRMIN):)...EXP(:SCA*H(NRMAX):)X
      do 11 i=nrmin,nrmax
      call taked(h,i,b0)
      call dacmud(b0,sca,bm)

         call expflo(bm,v(1),b0(1),eps,nmax)
         call dacop(b0(1),v(1))
 11   continue
      endif
      call dacop(v(1),w)
      call dadal(v(1),1)
      call dadal(b0,nd2)
      call dadal(bm,nd2)
      return
end subroutine facflo

subroutine facflod(h,x,w,nrmin,nrmax,sca,ifac)
      use floatPrecision
      use mod_lie_dab, only : nd,nd2
      implicit none
      integer i,ifac,nrmax,nrmin
      real(kind=fPrec) sca
! IFAC=1
! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX )  \VEC{X}= \VEC{Y}
! IFAC=-1
! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) \VEC{X}= \VEC{Y}
      integer x(*),w(*),h(*)

      do 1 i=1,nd2
      call facflo(h,x(i),w(i),nrmin,nrmax,sca,ifac)
 1    continue

      return
end subroutine facflod

subroutine fexpo(h,x,w,nrmin,nrmax,sca,ifac)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer ifac,nrma,nrmax,nrmi,nrmin
      real(kind=fPrec) sca
!   WRAPPED ROUTINES FOR THE OPERATOR  \VEC{H}=:H:
! WRAPPING FACFLOD
      integer x(*),w(*),h

      integer v(ndim2)

      nrmi=nrmin-1
      nrma=nrmax-1
      call etall(v,nd2)
      call difd(h,v,-one)
      call facflod(v,x,w,nrmi,nrma,sca,ifac)

      call dadal(v,nd2)

      return
end subroutine fexpo

subroutine etcom(x,y,h)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer i,j
! ETCOM TAKES THE BRACKET OF TWO VECTOR FIELDS.
      integer h(*),x(*),y(*),t1(1),t2(1),t3(ndim2)

      call etall(t1(1),1)
      call etall(t2(1),1)
      call etall(t3,nd2)

      do 2 j=1,nd2
      do 1 i=1,nd2

      call dader(i,x(j),t1(1))
      call dader(i,y(j),t2(1))
      call damul(x(i),t2(1),t2(1))
      call damul(y(i),t1(1),t1(1))
      call dalin(t2(1),one,t1(1),-one,t1(1))
      call daadd(t1(1),t3(j),t3(j))

 1    continue
 2    continue

      call dacopd(t3,h)

      call dadal(t1(1),1)
      call dadal(t2(1),1)
      call dadal(t3,nd2)
      return
end subroutine etcom

subroutine etpoi(x,y,h)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2
      implicit none
      integer i
! ETPOI TAKES THE POISSON BRACKET OF TWO FUNCTIONS
      integer h,x,y,t1(1),t2(1),t3(1)

      call etall(t1(1),1)
      call etall(t2(1),1)
      call etall(t3(1),1)

      do 1 i=1,nd

      call dader(2*i-1,x,t1(1))
      call dader(2*i,y,t2(1))
      call damul(t1(1),t2(1),t1(1))

      call dalin(t1(1),one,t3(1),one,t3(1))
      call dader(2*i-1,y,t1(1))
      call dader(2*i,x,t2(1))
      call damul(t1(1),t2(1),t1(1))

      call dalin(t1(1),-one,t3(1),one,t3(1))

 1    continue

      call dacop(t3(1),h)

      call dadal(t1(1),1)
      call dadal(t2(1),1)
      call dadal(t3(1),1)
      return
end subroutine etpoi

subroutine exp1d(h,x,y,eps,non)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,idpr,ndim2
      implicit none
      integer non
      real(kind=fPrec) eps
! WRAPPING EXPFLO
      integer h,x,y

      integer v(ndim2)

      call etall(v,nd2)
      call difd(h,v,-one)
      call expflo(v,x,y,eps,non)

      call dadal(v,nd2)

      return
end subroutine exp1d

subroutine expnd2(h,x,w,eps,nrmax)
      use floatPrecision
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer j,nrmax
      real(kind=fPrec) eps
! WRAPPING EXPFLOD USING EXP1D
      integer x(*),w(*),h

      integer b0(1),v(ndim2)

      call etallnom(b0(1),1,'B0        ')
      call etallnom(v,nd2  ,'V         ')

      call dacopd(x,v)

      do j=1,nd2
        call exp1d(h,v(j),b0(1),eps,nrmax)
        call dacop(b0(1),v(j))
      end do

      call dacopd(v,w)
      call dadal(v,nd2)
      call dadal(b0(1),1)

      return
end subroutine expnd2

subroutine flofacg(xy,h,epsone)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nd,nd2,idpr,ndim2,ntt,xintex
      implicit none
      integer i,k,kk,nrmax
      real(kind=fPrec) eps,epsone,r,xn,xnbefore,xnorm,xnorm1,xx
! GENERAL ONE EXPONENT FACTORIZATION
      logical more
      integer xy(*),x(ndim2),h(*)
      integer v(ndim2),w(ndim2),t(ndim2), z(ndim2)
      integer jj(ntt)
      jj(1)=1
!
      call etallnom(v,nd2  ,'V         ')
      call etallnom(w,nd2  ,'W         ')
      call etallnom(t,nd2  ,'T         ')
      call etallnom(x,nd2  ,'Z         ')
      call etallnom(z,nd2  ,'Z         ')

      call etini(v)
      call daclrd(w)
      xnorm1=zero
      do i=1,nd2
        call daabs(xy(i),r)
      xnorm1=xnorm1+r
       enddo
      xnbefore=1.0e36_fPrec
      more=.false.
      eps=c1m9
      nrmax=1000
      xn=c1e4
      do 333 k=1,nrmax
      call dacmud(h,-one,t)
      call expflod(t,xy,x,eps,nrmax)
      call dalind(x,one,v,-one,t)
! call daprid(t,1,1,20)
       if(xn.lt.epsone) then
            if(idpr.ge.0) write(lout,*) "xn quadratic",xn
        call daflod(t,t,w)
        call dalind(t,one,w,-half,t)
        call dacopd(t,z)
        call dacopd(t,w)
!  second order in W
        call etcom(h,w,x)
           call etcom(x,w,x)
!  END OF  order in W

              do kk=1,10
        call etcom(h,w,w)
        call dalind(z,one,w,xintex(kk),z)
              enddo
        call dacopd(z,t)
      xx=one/12.0_fPrec
      call dalind(x,xx,h,one,h)
       endif

      call dalind(t,one,h,one,h)
      xnorm=zero
      do i=1,nd2
        call daabs(t(i),r)
      xnorm=xnorm+r
       enddo
      xn=xnorm/xnorm1
      if(xn.ge.epsone.and.(idpr.ge.0)) write(lout,*)" xn linear ",xn
      if(xn.lt.eps.or.more) then
      more=.true.
      if(xn.ge.xnbefore) goto 1000
      xnbefore=xn
      endif
 333  continue
1000  write(lout,*) " iteration " , k
      call dadal(x,nd2)
      call dadal(w,nd2)
      call dadal(v,nd2)
      call dadal(t,nd2)
      call dadal(z,nd2)
      return
end subroutine flofacg

subroutine flofac(xy,x,h)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,no,ndim2
      implicit none
      integer k
! GENERAL DRAGT-FINN FACTORIZATION
      integer xy(*),x(*),h(*)
      integer v(ndim2),w(ndim2)
!
      call etallnom(v,nd2  ,'V         ')
      call etallnom(w,nd2  ,'W         ')

      call dacopd(xy,x)
      call dacopd(x,v)
      call daclrd(w)
      call danot(1)
      call etinv(v,w)
      call danot(no)
      call etcct(x,w,v)
      call danot(1)
      call dacopd(xy,x)
      call danot(no)
      call dacopd(v,w)
      call daclrd(h)
      do 333 k=2,no
      call taked(w,k,v)
      call dalind(v,one,h,one,h)
      call facflod(h,w,v,k,k,-one,-1)
      call dacopd(v,w)
 333  continue
      call dadal(w,nd2)
      call dadal(v,nd2)
      return
end subroutine flofac

subroutine liefact(xy,x,h)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
! SYMPLECTIC DRAGT-FINN FACTORIZATION WRAPPING FLOFAC
      integer xy(*),x(*),h

      integer v(ndim2)

      call etall(v,nd2)

      call flofac(xy,x,v)
      call intd(v,h,-one)
!
      call dadal(v,nd2)

      return
end subroutine liefact

subroutine mapnorm(x,ft,a2,a1,xy,h,nord)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer isi,nord
!--NORMALIZATION ROUTINES OF LIELIB
!- WRAPPING MAPNORMF
      integer x(*),a1(*),a2(*),ft,xy(*),h,hf(ndim2),ftf(ndim2)

      call etall(ftf,nd2)
      call etall(hf,nd2)
      isi=0
      call mapnormf(x,ftf,a2,a1,xy,hf,nord,isi)
      call intd(hf,h,-one)
      call intd(ftf,ft,-one)
      call dadal(ftf,nd2)
      call dadal(hf,nd2)

      return
end subroutine mapnorm

subroutine gettura(psq,radsq)
      use floatPrecision
      use mod_lie_dab, only : nd,nd2,ndim,ps,rads
      implicit none
      integer ik
      real(kind=fPrec) psq(ndim),radsq(ndim)

      do ik=1,nd
      psq(ik)=ps(ik)
      radsq(ik)=rads(ik)
      enddo

      return
end subroutine gettura

subroutine setidpr(idprint,nplan)
      use floatPrecision
      use mod_lie_dab, only : nd,nd2,idpr,nplane,epsplane,xplane,ndim
      implicit none
      integer idprint,ik,nplan
      dimension nplan(ndim)

      do ik=1,nd
      nplane(ik)=nplan(ik)
      enddo
      idpr=idprint

      return
end subroutine setidpr

subroutine idprset(idprint)
      use floatPrecision
      use mod_lie_dab, only : idpr
      implicit none
      integer idprint
      idpr=idprint
      return
end subroutine idprset

subroutine mapnormf(x,ft,a2,a1,xy,h,nord,isi)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nd,nd2,idpr,itu,iflow,jtune,ndc,ndc2,ndpt,ndt,ndim,ndim2,ps,rads
      implicit none
      integer ij,isi,nord
      real(kind=fPrec) angle,p,radn,st,x2pi,x2pii
      dimension angle(ndim),st(ndim),p(ndim),radn(ndim)
      integer x(*),a1(*),a2(*),ft(*),xy(*),h(*)
      integer a1i(ndim2),a2i(ndim2)
!
      call etallnom(a1i,nd2  ,'A1I       ')
      call etallnom(a2i,nd2  ,'A2I       ')
!     frank/etienne
      do itu=1,ndim
        angle(itu)=zero
        p(itu)=zero
        st(itu)=zero
        radn(itu)=zero
        ps(itu)=zero
        rads(itu)=zero
      enddo
      jtune=isi
      x2pii=(one/atan_mb(one))/eight                                    !hr11
      x2pi=atan_mb(one)*eight
      call dacopd(x,xy)
! goto fix point in the parameters + pt to order nord>=1
      call gofix(xy,a1,a1i,nord)
      call simil(a1i,xy,a1,xy)
! linear part
      call midbflo(xy,a2,a2i,angle,radn,st)
      do ij=1,nd-ndc
        p(ij)=angle(ij)*(st(ij)*(x2pii-one)+one)
      enddo
      if(ndc.eq.1) p(nd)=angle(nd)
      if(idpr.ge.0) then
        write(lout,*) 'tune    ',(p(ij),ij=1,nd)
        write(lout,*) 'damping ', (radn(ij),ij=1,nd)
      endif
      do ij=1,nd       !  -ndc    frank
        ps(ij)=p(ij)
        rads(ij)=radn(ij)
      enddo
      call initpert(st,angle,radn)
      call simil(a2i,xy,a2,xy)
      call dacopd(xy,a2i)
      call orderflo(h,ft,xy,angle,radn)
      do ij=1,nd-ndc
        p(ij)=angle(ij)
        if(angle(ij).gt.x2pi/two.and.st(ij).gt.zero.and.itu.eq.1)then
          p(ij)=angle(ij)-x2pi
          write(lout,*) ij,' TH TUNE MODIFIED IN H2 TO ',p(ij)/x2pi
        endif
      enddo
      call h2pluflo(h,p,radn)
!      CALL TAKED(A2I,1,XY)
      call taked(a2i,1,a1i)
      call etcct(xy,a1i,xy)

      call dadal(a2i,nd2)
      call dadal(a1i,nd2)
      return
end subroutine mapnormf

subroutine gofix(xy,a1,a1i,nord)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,no,ndc,ndc2,ndpt,ndt,ndim2
      implicit none
      integer i,nord
      real(kind=fPrec) xic
! GETTING TO THE FIXED POINT AND CHANGING TIME APPROPRIATELY IN THE
! COASTING BEAM CASE

!****************************************************************
! X = A1 XY A1I WHERE X IS TO THE FIXED POINT TO ORDER NORD
! for ndpt not zero, works in all cases. (coasting beam: eigenvalue
!1 in Jordan form)
!****************************************************************
      integer xy(*),a1(*),a1i(*)

      integer x(ndim2),w(ndim2),v(ndim2),rel(ndim2)
!
      call etallnom(x,nd2  ,  'X         ')
      call etallnom(w,nd2  ,  'W         ')
      call etallnom(v,nd2  ,  'V         ')
      call etallnom(rel,nd2  ,'REL       ')

! COMPUTATION OF A1 AND A1I USING DAINV
      call etini(rel)

      call danot(nord)

      call etini(v)

      do i=1,nd2-ndc2
        call dacop(xy(i),x(i))
        call dalin(x(i),one,rel(i),-one,v(i))
      enddo
      call etinv(v,w)
      call daclrd(x)
      if(ndc.eq.1) then
        call davar(x(ndpt),zero,ndpt)
      endif
      call etcct(w,x,v)
      if(ndc.eq.1) then
        call daclr(v(nd2))
        call daclr(v(nd2-ndc))
      endif
      call dalind(rel,one,v,one,a1)
      call dalind(rel,one,v,-one,a1i)

      if(ndpt.ne.0) then

!  CORRECTIONS
        call daclrd(w)
        call daclrd(v)
        call daclrd(x)

        do i=1,nd2-ndc2
          call dalin(a1(i),one,rel(i),-one,w(i))
        enddo

!      COMPUTE Deta/Ddelta
        call dacopd(w,a1)

        do i=1,nd2-ndc2
          call dader(ndpt,w(i),w(i))
        enddo
!      COMPUTE J*Deta/dDELTA

        do i=1,nd-ndc
          call dacmu(w(2*i),one,v(2*i-1) )
          call dacmu(w(2*i-1),-one,v(2*i) )
        enddo

        xic=real((-1)**(ndt),fPrec)                                            !hr11

        do i=1,nd2-ndc2
          call damul(v(i),rel(i),x(1))
          call daadd(x(1),w(ndt),w(ndt))
          call dacop(a1(i),w(i))
        enddo
        call dacmu(w(ndt),xic,w(ndt))

        call expflod(w,rel,a1,c1m7,10000)
! END OF  CORRECTIONS

        call etinv(a1,a1i)

      endif

      call danot(no)

      call dadal(rel,nd2)
      call dadal(v,nd2)
      call dadal(w,nd2)
      call dadal(x,nd2)
      return
end subroutine gofix

real(kind=fPrec) function transver(j)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt
      implicit none
      integer i,ic
! USED IN A DACFU CALL OF GOFIX
!      INTEGER J(NTT)
      integer j(*)

      transver=one
      ic=0
      do i=1,nd2-ndc2
        ic=ic+j(i)
      enddo
      if(ic.ne.1) transver=zero
      return
end function transver

subroutine orderflo(h,ft,x,ang,ra)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nd,nd2,no,idpr,ndim,ndim2
      implicit none
      integer k
      real(kind=fPrec) ang,ra
!-   NONLINEAR NORMALIZATION PIECE OF MAPNORMF
      dimension ang(ndim),ra(ndim)
      integer x(*),h(*),ft(*)
      integer w(ndim2),v(ndim2),rel(ndim2)
      integer roi(ndim2)
      integer b1(ndim2),b5(ndim2),b6(ndim2),b9(ndim2)
!
      call etallnom(w,nd2  ,'W         ')
      call etallnom(v,nd2  ,'V         ')
      call etallnom(rel,nd2  ,'REL       ')
      call etallnom(roi,nd2  ,'ROI       ')
      call etallnom(b1,nd2  ,'B1        ')
      call etallnom(b5,nd2  ,'B5        ')
      call etallnom(b6,nd2  ,'B6        ')
      call etallnom(b9,nd2  ,'B9        ')
      call rotiflo(roi,ang,ra)
      call etini(rel)
      call daclrd(h)
      call daclrd(ft)
      call etcct(x,roi,x)
      do 33 k=2,no
! IF K>2 V = H(K)^-1 X(K)
      call facflod(h,x,v,2,k-1,-one,-1)
! EXTRACTING K TH DEGREE OF V ----> W
      call taked(v,k,w)
! W = EXP(B5) + ...
       call dacopd(w,b5)
!      CALL INTD(W,B5,-1.D0)
! B5 ON EXIT IS THE NEW CONTRIBUTION TO H
! B6 IS THE NEW CONTRIBUTION TO FT
       call nuanaflo(b5,b6)
      call dalind(b5,one,h,one,b1)
      call dacopd(b1,h)
! EXP(B9) = EXP( : ROTI B6 :)
      call trxflo(b6,b9,roi)

! V = EXP(-B6) REL
      call facflod(b6,rel,v,k,k,-one,1)
! W = V o X
      call etcct(v,x,w)
      if(idpr.ge.0) then
        write(lout,*) ' ORDERFLO K= ', k
      endif
! X = EXP(B9) W
      call facflod(b9,w,x,k,k,one,1)
! B6 IS THE NEW CONTRIBUTION TO FT
      call dalind(b6,one,ft,one,b1)
      call dacopd(b1,ft)
 33   continue
      call dadal(b9,nd2)
      call dadal(b6,nd2)
      call dadal(b5,nd2)
      call dadal(b1,nd2)
      call dadal(roi,nd2)
      call dadal(rel,nd2)
      call dadal(v,nd2)
      call dadal(w,nd2)
      return
end subroutine orderflo

subroutine nuanaflo(h,ft)
      use floatPrecision
      use mod_lie_dab, only : nd,nd2,iflow,jtune,ndim2
      implicit none
      integer i
      real(kind=fPrec) dfilt,filt,xgam,xgbm
! RESONANCE DENOMINATOR OPERATOR (1-R^-1)^-1
      external xgam,xgbm,dfilt,filt
      integer h(*),ft(*),br(ndim2),bi(ndim2),c(ndim2),ci(ndim2)
      integer t1(2),t2(2)

      call etall(br,nd2)
      call etall(bi,nd2)
      call etall(c,nd2)
      call etall(ci,nd2)

      call ctorflo(h,br,bi)

! FILTERING RESONANCES AND TUNE SHIFTS
! ASSUMING REALITY I.E. B(2*I-1)=CMPCJG(B(2*I))

      do 2 i=1,nd2
      iflow=i
       call dacfu(br(i),filt,c(i))
       call dacfu(bi(i),filt,ci(i))
 2    continue
      call rtocflo(c,ci,h)

      do 22 i=1,nd2

       iflow=i
       call dacfu(br(i),dfilt,br(i))
       call dacfu(bi(i),dfilt,bi(i))
 22    continue
!  NOW WE MUST REORDER C AND CI TO SEPARATE THE REAL AND IMAGINARY PART
! THIS IS NOT NECESSARY WITH :H: OPERATORS

      do 3 i=1,nd2
      t1(1)=br(i)
      t1(2)=bi(i)
      t2(1)=c(i)
      t2(2)=ci(i)
      iflow=i
      call comcfu(t1,xgam,xgbm,t2)
 3    continue

      call rtocflo(c,ci,ft)

      call dadal(br,nd2)
      call dadal(bi,nd2)
      call dadal(c,nd2)
      call dadal(ci,nd2)

      return
end subroutine nuanaflo

real(kind=fPrec) function xgam(j)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,iflow,jtune,ndc,ndc2,ndpt,ndt,ndim,ndim2,idsta,ista,angle,dsta,radn,sta
      implicit none
      integer i,ic,ij,ik
      real(kind=fPrec) ad,ans,as,ex,exh
! XGAM AND XGBM ARE THE EIGENVALUES OF THE OPERATOR NEWANAFLO
!      INTEGER J(NTT),JJ(NDIM),JP(NDIM)
      integer j(*),jj(ndim),jp(ndim)
      xgam=zero
      ad=zero
      as=zero
      ic=0
      do i=1,nd-ndc
        ik=2*i-1
        ij=2*i
        jp(i)=j(ik)+j(ij)
        jj(i)=j(ik)-j(ij)
        if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
          jp(i)=jp(i)-1
        endif
        ic=ic+abs(jj(i))                                                 !hr11
      enddo

      do i=1,nd-ndc
        ad=((dsta(i)*real(jj(i),fPrec))*angle(i)-real(jp(i),fPrec)*radn(i))+ad
        as=(sta(i)*real(jj(i),fPrec))*angle(i)+as                              !hr11
      enddo

      exh=exp_mb(ad/two)
      ex=exh**2
      ans=(four*ex)*(sinh_mb(ad/two)**2+sin_mb(as/two)**2)             !hr11
      xgam=(two*(ex*sin_mb(as/two)**2-exh*sinh_mb(ad/two)))/ans       !hr11

      return
end function xgam

real(kind=fPrec) function xgbm(j)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,iflow,jtune,ndc,ndc2,ndpt,ndt,ndim,ndim2,idsta,ista,angle,dsta,radn,sta
      implicit none
      integer i,ic,ij,ik
      real(kind=fPrec) ad,ans,as,ex,exh
!      INTEGER J(NTT),JJ(NDIM),JP(NDIM)
      integer j(*),jj(ndim),jp(ndim)
      xgbm=zero
      ad=zero
      as=zero
      ic=0
      do i=1,nd-ndc
        ik=2*i-1
        ij=2*i
        jp(i)=j(ik)+j(ij)
        jj(i)=j(ik)-j(ij)
        if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
          jp(i)=jp(i)-1
        endif
        ic=ic+iabs(jj(i))
      enddo

      do i=1,nd-ndc
        ad=((dsta(i)*real(jj(i),fPrec))*angle(i)-real(jp(i),fPrec)*radn(i))+ad
        as=(sta(i)*real(jj(i),fPrec))*angle(i)+as                              !hr11
      enddo

      exh=exp_mb(ad/two)
      ex=exh**2
      ans=(four*ex)*(sinh_mb(ad/two)**2+sin_mb(as/two)**2)             !hr11
      xgbm=(sin_mb(as)*ex)/ans                                           !hr11

      return
end function xgbm

real(kind=fPrec) function filt(j)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,iflow,jtune,ndc,ndc2,ndpt,ndt,ndim,idsta,ista,angle,dsta,sta,mx,nres
      implicit none
      integer i,ic,ic1,ic2,ij,ik,ji
!  PROJECTION FUNCTIONS ON THE KERNEL ANMD RANGE OF (1-R^-1)
!-  THE KERNEL OF (1-R^-1)
!      INTEGER J(NTT),JJ(NDIM)
      integer j(*),jj(ndim)

      filt=one

      ic=0
      do i=1,nd-ndc
        ik=2*i-1
        ij=2*i
        jj(i)=j(ik)-j(ij)
        if(ik.eq.iflow.or.ij.eq.iflow) then
          jj(i)=jj(i)+(-1)**iflow
        endif
        ic=ic+iabs(jj(i))
      enddo

      if(ic.eq.0.and.jtune.eq.0) return

      do i=1,nres
        ic1=1
        ic2=1
        do ji=1,nd-ndc
          if(mx(ji,i).ne.jj(ji)) ic1=0
          if(mx(ji,i).ne.-jj(ji)) ic2=0
          if(ic1.eq.0.and.ic2.eq.0) goto 3
        enddo
        return
 3      continue
      enddo

      filt=zero
      return
end function filt

real(kind=fPrec) function dfilt(j)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt,mx,nres
      implicit none
      real(kind=fPrec) fil,filt
!-  THE RANGE OF (1-R^-1)^1
!- CALLS FILT AND EXCHANGES 1 INTO 0 AND 0 INTO 1.
      external filt
!      INTEGER J(NTT)
      integer j(*)

      fil=filt(j)
      if(fil.gt.half) then
       dfilt=zero
      else
       dfilt=one
      endif
      return
end function dfilt

subroutine dhdjflo(h,t)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt,ndim2
      implicit none
      integer i
      real(kind=fPrec) coe,x2pi
! CONVENIENT TUNE SHIFT FINDED FOR SYMPLECTIC CASE (NU,DL)(H)=T
      integer h(*),t(*)

      integer b1(ndim2),b2(ndim2),bb1(1),bb2(1)
!
      call etall(b1,nd2)
      call etall(b2,nd2)
      call etall(bb1(1),1)
      call etall(bb2(1),1)

      x2pi=atan_mb(one)*eight
      call ctorflo(h,b1,b2)
      coe=one/x2pi

      do i=1,nd-ndc
        call datra(2*i,b2(2*i),bb1(1))
        call dacmu(bb1(1),coe,t(i+nd))
        call dacop(t(i+nd),bb1(1))
        call daclr(bb2(1))
        call rtoc(bb1(1),bb2(1),bb1(1))
        call dacop(bb1(1),t(i))
      enddo

      if(ndpt.ne.0) then
        call dacop(h(ndt),t(nd))
        call dacop(b1(ndt),t(nd2))
      endif

      call dadal(bb2(1),1)
      call dadal(bb1(1),1)
      call dadal(b2,nd2)
      call dadal(b1,nd2)
      return
end subroutine dhdjflo

subroutine dhdj(h,t)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt
      implicit none
      integer i
      real(kind=fPrec) coe,x2pi
      integer h,t(*)

      integer b1(1),b2(1),bb1(1),bb2(1)
!
      call etallnom(b1(1),1,'B1        ')
      call etallnom(b2(1),1,'B2        ')
      call etallnom(bb1(1),1,'BB1       ')
      call etallnom(bb2(1),1,'BB2       ')

      x2pi=atan_mb(one)*eight
      call ctor(h,b1(1),b2(1))
      coe=-two/x2pi
      do i=1,nd-ndc
        call dader(2*i-1,b1(1),b2(1))
        call datra(2*i,b2(1),bb2(1))
        call dacmu(bb2(1),coe,t(i+nd))
        call dacop(t(i+nd),bb2(1))
        call daclr(b2(1))
        call rtoc(bb2(1),b2(1),bb1(1))
        call dacop(bb1(1),t(i))
      enddo

      if(ndpt.eq.nd2) then
        call dader(ndpt,h,t(nd))
        call dader(ndpt,b1(1),t(nd2))
        call dacmu(t(nd),-one,t(nd))
        call dacmu(t(nd2),-one,t(nd2))
      endif
      if(ndt.eq.nd2) then
        call dader(ndpt,h,t(nd))
        call dader(ndpt,b1,t(nd2))
      endif
      call dadal(bb2(1),1)
      call dadal(bb1(1),1)
      call dadal(b2(1),1)
      call dadal(b1(1),1)
      return
end subroutine dhdj

subroutine h2pluflo(h,ang,ra)
      use floatPrecision
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,no,nv,ndc,ndc2,ndpt,ndt,ndim,ntt,angle,dsta,sta
      implicit none
      integer i,j
      real(kind=fPrec) ang,r1,r2,ra,st
! POKES IN \VEC{H}  ANGLES AND DAMPING COEFFFICIENTS
      dimension ang(ndim),st(ndim),ra(ndim),j(ntt)
      integer h(*)
!*DAEXT(NO,NV) H

      do i=1,nd
        st(i)=two*sta(i)-one
      enddo

      do i=1,ntt
        j(i)=0
      enddo

      do i=1,nd-ndc
        j(2*i-1)=1
        r1=-one*ang(i)                                                   !hr11
!-----
        call dapok(h(2*i),j,r1)

        r2=ra(i)
        call dapok(h(2*i-1),j,r2)
        j(2*i-1)=0

        j(2*i)=1
        r1=ang(i)*st(i)
        call dapok(h(2*i-1),j,r1)
        call dapok(h(2*i),j,r2)
        j(2*i)=0

      enddo

      if(ndpt.eq.nd2-1) then
        j(ndpt)=1
        call dapok(h(ndt),j,ang(nd))
      elseif(ndpt.eq.nd2) then
        j(ndpt)=1
        call dapok(h(ndt),j,-one*ang(nd))                                !hr11
      endif
      return
end subroutine h2pluflo

subroutine rotflo(ro,ang,ra)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt,ndim,ntt,idsta,ista
      implicit none
      integer i
      real(kind=fPrec) ang,ch,co,ra,sh,si,sim,xx
! CREATES R AND R^-1 USING THE EXISTING ANGLES AND DAMPING
! COULD BE REPLACED BY A CALL H2PLUFLO FOLLOWED BY EXPFLOD
! CREATES R
      dimension co(ndim),si(ndim),ang(ndim),ra(ndim)
      integer j(ntt)
      integer ro(*)
      call daclrd(ro)
      do i=1,nd-ndc
        xx=exp_mb(ra(i))
        if(ista(i).eq.0) then
          call hyper(ang(i),ch,sh)
          co(i)=ch*xx
          si(i)=(-one*sh)*xx                                             !hr11
        else
          co(i)=cos_mb(ang(i))*xx
          si(i)=sin_mb(ang(i))*xx
        endif
      enddo
      do i=1,nd-ndc
        if(ista(i).eq.0)then
          sim=si(i)
        else
          sim=-one*si(i)                                                 !hr11
        endif
        j(2*i-1)=1
        call dapok(ro(2*i-1),j,co(i))
        call dapok(ro(2*i),j,sim)
        j(2*i-1)=0
        j(2*i)=1
        call dapok(ro(2*i),j,co(i))
        call dapok(ro(2*i-1),j,si(i))
        j(2*i)=0
      enddo

      if(ndc.eq.1) then
        j(ndt)=1
        call dapok(ro(ndt),j,one)
        call dapok(ro(ndpt),j,zero)
        j(ndt)=0
        j(ndpt)=1
        call dapok(ro(ndt),j,ang(nd))
        call dapok(ro(ndpt),j,one)
        j(ndpt)=0
      endif

      return
end subroutine rotflo

subroutine rotiflo(roi,ang,ra)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt,ndim,ntt,idsta,ista
      implicit none
      integer i
      real(kind=fPrec) ang,ch,co,ra,sh,si,sim,simv,xx
! CREATES  R^-1
      dimension co(ndim),si(ndim),ang(ndim),ra(ndim)
      integer j(ntt)
      integer roi(*)

      do i=1,10
        j(i)=0
      enddo

      call daclrd(roi)
      do i=1,nd-ndc
        xx=exp_mb(-one*ra(i))                                            !hr11
        if(ista(i).eq.0) then
          call hyper(ang(i),ch,sh)
          co(i)=ch*xx
          si(i)=(-one*sh)*xx
        else
          co(i)=cos_mb(ang(i))*xx
          si(i)=sin_mb(ang(i))*xx
        endif
      enddo
      do i=1,nd-ndc
        if(ista(i).eq.0)then
          sim=si(i)
        else
          sim=-one*si(i)                                                 !hr11
        endif
        j(2*i-1)=1
        call dapok(roi(2*i-1),j,co(i))
        simv=-one*sim                                                    !hr11
        call dapok(roi(2*i),j,simv)
        j(2*i-1)=0
        j(2*i)=1
        simv=-one*si(i)                                                  !hr11
        call dapok(roi(2*i),j,co(i))
        call dapok(roi(2*i-1),j,simv)
        j(2*i)=0
      enddo

      if(ndc.eq.1) then
        j(ndt)=1
        call dapok(roi(ndt),j,one)
        call dapok(roi(ndpt),j,zero)
        j(ndt)=0
        j(ndpt)=1
        call dapok(roi(ndt),j,-ang(nd))
        call dapok(roi(ndpt),j,one)
        j(ndpt)=0
      endif

      return
end subroutine rotiflo

subroutine hyper(a,ch,sh)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
      real(kind=fPrec) a,ch,sh,x,xi
!   USED IN ROTIFLO AND ROTFLO
      x=exp_mb(a)
      xi=one/x
      ch=(x+xi)/two
      sh=(x-xi)/two
      return
end subroutine hyper

subroutine ctor(c1,r2,i2)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
! CHANGES OF BASIS
!   C1------> R2+I R1
      integer c1,r2,i2
      integer b1(1),b2(1),x(ndim2)
!
!
      call etallnom(b1(1),1,'B1        ')
      call etallnom(b2(1),1,'B2        ')
      call etallnom(x,nd2  ,'X         ')

      call ctoi(c1,b1(1))
      call etcjg(x)
      call trx(b1(1),b2(1),x)
      call dalin(b1(1),half,b2(1),half,r2)
      call dalin(b1(1),half,b2(1),-half,i2)
      call dadal(x,nd2)
      call dadal(b2(1),1)
      call dadal(b1(1),1)
      return
end subroutine ctor

subroutine rtoc(r1,i1,c2)
      use floatPrecision
      use mathlib_bouncer
      implicit none
!  INVERSE OF CTOR
      integer c2,r1,i1

      integer b1(1)
!
      call etallnom(b1(1),1,'B1        ')

      call daadd(r1,i1,b1(1))
      call itoc(b1(1),c2)
      call dadal(b1(1),1)
      return
end subroutine rtoc

subroutine ctorflo(c,dr,di)
      use floatPrecision
      use mathlib_bouncer
      implicit none
! FLOW CTOR
      integer dr(*),di(*),c(*)

      call ctord(c,dr,di)
      call resvec(dr,di,dr,di)

      return
end subroutine ctorflo

subroutine rtocflo(dr,di,c)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
! FLOW RTOC
      integer dr(*),di(*),c(*),er(ndim2),ei(ndim2)

      call etall(er,nd2)
      call etall(ei,nd2)

      call reelflo(dr,di,er,ei)
      call rtocd(er,ei,c)

      call dadal(er,nd2)
      call dadal(ei,nd2)

      return
end subroutine rtocflo

subroutine ctord(c,cr,ci)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd2
      implicit none
      integer i
! ROUTINES USED IN THE INTERMEDIATE STEPS OF CTORFLO AND RTOCFLO
! SAME AS CTOR  OVER ARRAYS CONTAINING ND2 COMPONENTS
! ROUTINE USEFUL IN INTERMEDIATE FLOW CHANGE OF BASIS
      integer c(*),ci(*),cr(*)
      do 1 i=1,nd2
      call ctor(c(i),cr(i),ci(i))
 1    continue
      return
end subroutine ctord

subroutine rtocd(cr,ci,c)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd2
      implicit none
      integer i
!  INVERSE OF CTORD
      integer c(*),ci(*),cr(*)
      do 1 i=1,nd2
      call rtoc(cr(i),ci(i),c(i))
 1    continue
      return
end subroutine rtocd

subroutine resvec(cr,ci,dr,di)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt,idsta,ista,angle,dsta,sta
      implicit none
      integer i
! DOES THE SPINOR PART IN CTORFLO
      integer dr(*),di(*),ci(*),cr(*),tr(2),ti(2)

      call etall(tr,2)
      call etall(ti,2)

      do i=1,nd-ndc
        if(ista(i).eq.1) then
          call dasub(cr(2*i-1),ci(2*i),tr(1))
          call daadd(ci(2*i-1),cr(2*i),ti(1))
          call daadd(cr(2*i-1),ci(2*i),tr(2))
          call dasub(ci(2*i-1),cr(2*i),ti(2))
          call dacop(tr(1),dr(2*i-1))
          call dacop(tr(2),dr(2*i))
          call dacop(ti(1),di(2*i-1))
          call dacop(ti(2),di(2*i))
        else
          call daadd(cr(2*i-1),cr(2*i),tr(1))
          call daadd(ci(2*i-1),ci(2*i),ti(1))
          call dasub(cr(2*i-1),cr(2*i),tr(2))
          call dasub(ci(2*i-1),ci(2*i),ti(2))
          call dacop(tr(1),dr(2*i-1))
          call dacop(tr(2),dr(2*i))
          call dacop(ti(1),di(2*i-1))
          call dacop(ti(2),di(2*i))
        endif
      enddo

      do i=nd2-ndc2+1,nd2
        call dacop(cr(i),dr(i))
        call dacop(ci(i),di(i))
      enddo

      call dadal(tr,2)
      call dadal(ti,2)
      return
end subroutine resvec

subroutine reelflo(c,ci,f,fi)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt,ndim2,idsta,ista,angle,dsta,sta
      implicit none
      integer i
! DOES THE SPINOR PART IN RTOCFLO
      integer c(*),ci(*),f(*),fi(*),e(ndim2),ei(ndim2)

      call etall(e,nd2)
      call etall(ei,nd2)

      do i=1,nd-ndc
        call dalin(c(2*i-1),half,c(2*i),half,e(2*i-1))
        call dalin(ci(2*i-1),half,ci(2*i),half,ei(2*i-1))
        if(ista(i).eq.1) then
          call dalin(ci(2*i-1),half,ci(2*i),-half,e(2*i))
          call dalin(c(2*i-1),-half,c(2*i),half,ei(2*i))
        else
          call dalin(ci(2*i-1),half,ci(2*i),-half,ei(2*i))
          call dalin(c(2*i-1),half,c(2*i),-half,e(2*i))
        endif
      enddo

      do i=nd2-ndc2+1,nd2
        call dacop(c(i),e(i))
        call dacop(ci(i),ei(i))
      enddo

      call dacopd(e,f)
      call dacopd(ei,fi)

      call dadal(e,nd2)
      call dadal(ei,nd2)
      return
end subroutine reelflo

subroutine compcjg(cr,ci,dr,di)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd2,ndim2
      implicit none
! TAKES THE COMPLEX CONJUGATE IN RESONANCE BASIS OF A POLYNOMIAL
      integer dr,di,ci,cr,x(ndim2)

      call etall(x,nd2)

      call etcjg(x)
      call trx(cr,dr,x)
      call trx(ci,di,x)
      call dacmu(di,-one,di)

      call dadal(x,nd2)
      return
end subroutine compcjg

subroutine midbflo(c,a2,a2i,q,a,st)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,no,nv,ndc,ndc2,ndpt,ndt,ndim,ndim2,ntt
      implicit none
      integer i,j
      real(kind=fPrec) a,ch,cm,cr,q,r,sa,sai,shm,st,x2pi
! LINEAR EXACT NORMALIZATION USING EIGENVALUE PACKAGE OF NERI
      integer jx(ntt)
      dimension cr(ndim2,ndim2),st(ndim),q(ndim),a(ndim)
      dimension sa(ndim2,ndim2),sai(ndim2,ndim2),cm(ndim2,ndim2)
      integer c(*),a2(*),a2i(*)
!*DAEXT(NO,NV) C(NDIM2),A2(NDIM2),A2I(NDIM2)
      x2pi=atan_mb(one)*eight

      do i=1,ntt
        jx(i)=0
      enddo

!     frank/etienne
      do i=1,ndim
        st(i)=zero
        q(i)=zero
        a(i)=zero
      enddo
!     frank/etienne
      do i=1,ndim2
!     frank/etienne
        do j=1,ndim2
          sai(i,j)=zero
          sa(i,j)=zero
          cm(i,j)=zero
          cr(i,j)=zero
        enddo
      enddo

      do i=1,nd2
        do j=1,nd2
          jx(j)=1
          call  dapek(c(i),jx,r)
          jx(j)=0
          cm(i,j)=r
        enddo
      enddo

      call mapflol(sa,sai,cr,cm,st)
      do i=1,nd-ndc
        if(st(i)+c1m3.gt.one) then                                   !hr11
          a(i)=sqrt(cr(2*i-1,2*i-1)**2+cr(2*i-1,2*i)**2)                 !hr11
          q(i)=acos_mb(cr(2*i-1,2*i-1)/a(i))
          a(i)=log_mb(a(i))
          if(cr(2*i-1,2*i).lt.zero) q(i)=x2pi-q(i)
        else
          a(i)=sqrt(cr(2*i-1,2*i-1)**2-cr(2*i-1,2*i)**2)                 !hr11
          ch=cr(2*i-1,2*i-1)/a(i)
          shm=cr(2*i-1,2*i)/a(i)
!       CH=CH+DSQRT(CH**2-1.D0)
!       q(i)=DLOG(CH)
          q(i)=-one*log_mb(ch+shm)                                       !hr11
!       IF(cr(2*i-1,2*i).gt.0.d0) Q(I)=-Q(I)
          a(i)=log_mb(a(i))
        endif
      enddo

      if(ndc.eq.0) then
        if(st(3)+c1m3.gt.one.and.nd.eq.3.and.q(nd).gt.half)             &!hr11
     &q(3)=q(3)-x2pi
      else
        q(nd)=cr(ndt,ndpt)
      endif

      call daclrd(a2)
      call daclrd(a2i)

      do i=1,nd2
        do j=1,nd2
          jx(j)=1
          r=sa(i,j)
          if(r.ne.zero)call  dapok(a2(i),jx,r)
          jx(j)=1
          r=sai(i,j)
          if(r.ne.zero)call  dapok(a2i(i),jx,r)
          jx(j)=0
        enddo
      enddo

      return
end subroutine midbflo

subroutine mapflol(sa,sai,cr,cm,st)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,no,idpr,ndc,ndc2,ndpt,ndt,nplane,epsplane,xplane,ndim,ndim2
      use crcoall
      use mod_units
      use mod_meta
      implicit none

      integer i,ier,iunst,j,l,n,n1
      real(kind=fPrec) ap,ax,cm,cr,p,rd,rd1,ri,rr,s1,sa,sai,st,vi,vr,w,x,x2pi,xd,xj,xsu,xx
!---- FROM TRACKING CODE
! ---------------------
      dimension cr(ndim2,ndim2),xj(ndim2,ndim2),n(ndim),x(ndim)
      dimension rr(ndim2),ri(ndim2),sa(ndim2,ndim2),xx(ndim),sai(ndim2,ndim2),cm(ndim2,ndim2),w(ndim2,ndim2),st(ndim)
      dimension vr(ndim2,ndim2),vi(ndim2,ndim2),s1(ndim2,ndim2),p(ndim2)
      logical lopen
      x2pi=atan_mb(one)*eight
      n1=0
!     frank/etienne
      do i=1,ndim2
        do j=1,ndim2
          cr(j,i)=cm(i,j)
          xj(i,j)=zero
          s1(i,j)=zero
        enddo
      enddo
!     frank/etienne
      do i=1,ndim
        n(i)=0
        xj(2*i-1,2*i)=one
        xj(2*i,2*i-1)=-one
      enddo
!     frank/etienne
      do i=1,ndim2
        do j=1,ndim2
          sai(i,j)=zero
          w(i,j)=cm(i,j)
        enddo
      enddo
      if(ndc.eq.1) then
        s1(nd2-ndc,nd2-ndc)=one
        s1(nd2,nd2)=one
        sai(nd2-ndc,nd2-ndc)=one
        sai(nd2,nd2)=one
      endif
      call mulnd2(xj,w)
      call mulnd2(cr,w)
      if(idpr.ge.0.or.idpr.eq.-102) then
        write(lout,"(a)") "LIELIB> Check of the symplectic condition on the linear part:"
        xsu=zero
        do i=1,nd2
          write(lout,"(3x,6(1x,1pe17.10))") (w(i,j), j=1, nd2)
          do j=1,nd2
            xsu=xsu+abs(w(i,j))
          enddo
        enddo
        ! Report
        meta_sympCheck = (xsu-nd2)/xsu
        write(lout,"(a,es13.6,a)") "LIELIB> Deviation for symplecticity = ",c1e2*meta_sympCheck," %"
      endif
      call eig6(cr,rr,ri,vr,vi)
      if(idpr.ge.0) then
        write(lout,*) '   '
        write(lout,*) '       Index         Real Part  ','       ArcSin(Imaginary Part)/2/pi'
        write(lout,*) '   '
        do i=1,nd-ndc
          rd1=sqrt(rr(2*i-1)**2+ri(2*i-1)**2)
          rd=sqrt(rr(2*i)**2+ri(2*i)**2)
          write(lout,*) 2*i-1,rr(2*i-1),asin_mb(ri(2*i-1)/rd1)/x2pi
          write(lout,*) 2*i,rr(2*i),asin_mb(ri(2*i)/rd)/x2pi
          write(lout,*) ' alphas ', log_mb(sqrt(rd*rd1))
        enddo
        write(lout,*) ' select ',nd-ndc,' eigenplanes (odd integers <0 real axis)'
        read(5,*) (n(i),i=1,nd-ndc)
      elseif(idpr.eq.-100) then
        do i=1,nd-ndc
          n(i)=nplane(i)
        enddo
      elseif(idpr.eq.-101.or.idpr.eq.-102) then
        do i=1,nd-ndc
          if(ri(2*i).ne.zero) then
            n(i)=2*i-1
          else
            n(i)=-2*i+1
          endif
        enddo
      else
        do i=1,nd-ndc
          n(i)=2*i-1
        enddo
      endif
      iunst=0
      do i=1,nd-ndc                  ! frank ndc  kept
        if(n(i).lt.0) then
          n(i)=-n(i)
          st(i)=zero
          iunst=1
        else
          st(i)=one
        endif
        x(i)=zero
        xx(i)=one
        do j=1,nd-ndc
          x(i)=vr(2*j-1,n(i))*vi(2*j,n(i))-vr(2*j,n(i))*vi(2*j-1,n(i))+x(i)
        enddo
      enddo

      do i=1,nd-ndc
        if(x(i).lt.zero) xx(i)=-one
        x(i)=sqrt(abs(x(i)))
      enddo
      do i=1,nd2-ndc2
        do j=1,nd-ndc
          if(st(j)+c1m3.gt.one) then
            sai(2*j-1,i)=vr(i,n(j))*xx(j)/x(j)
            sai(2*j,i)=vi(i,n(j))/x(j)
          else
            ax=vr(i,n(j))*xx(j)/x(j)
            ap=vi(i,n(j))/x(j)
            sai(2*j-1,i)=(ax+ap)/sqrt(two)
            sai(2*j,i)=(ap-ax)/sqrt(two)
          endif
        enddo
      enddo
      if(idpr.eq.-101.or.idpr.eq.-102) then
        call movearou(sai)
      endif
! adjust sa such that sa(1,2)=0 and sa(3,4)=0. (courant-snyder-edwards-teng
! phase advances)
      if(iunst.ne.1) then
        do i=1,nd-ndc
          p(i)=atan_mb(-sai(2*i-1,2*i)/sai(2*i,2*i))
          s1(2*i-1,2*i-1)=cos_mb(p(i))
          s1(2*i,2*i)=cos_mb(p(i))
          s1(2*i-1,2*i)=sin_mb(p(i))
          s1(2*i,2*i-1)=-sin_mb(p(i))
        enddo
        call mulnd2(s1,sai)
! adjust sa to have sa(1,1)>0 and sa(3,3)>0 rotate by pi if necessary.
        do i=1,nd-ndc
          xd=one
          if(sai(2*i-1,2*i-1).lt.zero) xd=-one
          s1(2*i-1,2*i-1)=xd
          s1(2*i-1,2*i)=zero
          s1(2*i,2*i-1)=zero
          s1(2*i,2*i)=xd
        enddo
        call mulnd2(s1,sai)
! sa is now uniquely and unambigeously determined.
      endif
      do i=1,nd2
        do l=1,nd2
          sa(i,l)=sai(i,l)
        enddo
      enddo
      call matinv(sai,sa,nd2,6,ier)

      call mulnd2(sai,cm)
      do i=1,nd2
        do j=1,nd2
          cr(i,j)=sa(i,j)
        enddo
      enddo

      call mulnd2(cm,cr)

      return
end subroutine mapflol

subroutine mulnd2(rt,r)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer i,ia,j
      real(kind=fPrec) r,rt,rtt
      dimension rt(ndim2,ndim2),r(ndim2,ndim2),rtt(ndim2,ndim2)
      do i=1,nd2
        do j=1,nd2
          rtt(i,j)=zero
        end do
      end do

      do i=1,nd2
        do j=1,nd2
          do ia=1,nd2
            rtt(i,ia)=rt(i,j)*r(j,ia)+rtt(i,ia)
          end do
        end do
      end do

      do i=1,nd2
        do j=1,nd2
          r(i,j)=rtt(i,j)
        end do
      end do

      return
end subroutine mulnd2

subroutine movearou(rt)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nd,nd2,no,idpr,ndim2
      implicit none
      integer i,ic,j
      real(kind=fPrec) rt,rto,s,xr,xrold,xy,xyz,xz,xzy,yz
      dimension rt(ndim2,ndim2),rto(ndim2,ndim2)
      dimension xy(ndim2,ndim2),xz(ndim2,ndim2),yz(ndim2,ndim2)
      dimension xyz(ndim2,ndim2),xzy(ndim2,ndim2)
      dimension s(ndim2,ndim2)

      do i=1,nd2
        do j=1,nd2
          s(i,j)=zero
          s(i,i)=one
          xy(i,j)=zero
          xz(i,j)=zero
          yz(i,j)=zero
          xyz(i,j)=zero
          xzy(i,j)=zero
        end do
      end do

      xy(1,3)=one
      xy(3,1)=one
      xy(2,4)=one
      xy(4,2)=one
      xy(5,5)=one
      xy(6,6)=one

      xz(1,5)=one
      xz(5,1)=one
      xz(2,6)=one
      xz(6,2)=one
      xz(3,3)=one
      xz(4,4)=one

      yz(3,5)=one
      yz(5,3)=one
      yz(4,6)=one
      yz(6,4)=one
      yz(1,1)=one
      yz(2,2)=one

      xyz(1,3)=one
      xyz(3,5)=one
      xyz(5,1)=one
      xyz(2,4)=one
      xyz(4,6)=one
      xyz(6,2)=one

      xzy(1,5)=one
      xzy(5,3)=one
      xzy(3,1)=one
      xzy(2,6)=one
      xzy(6,4)=one
      xzy(4,2)=one

      ic=0
      xrold=1000000000.0_fPrec
      call movemul(rt,s,rto,xr)
      if(xr.lt.xrold) then
        xrold=xr
      endif

      if(nd.ge.2) then
        call movemul(rt,xy,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=1
        endif
      endif

      if(nd.eq.3) then
        call movemul(rt,xz,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=2
        endif
        call movemul(rt,yz,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=3
        endif
        call movemul(rt,xyz,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=4
        endif
        call movemul(rt,xzy,rto,xr)
        if(xr.lt.xrold) then
          xrold=xr
          ic=5
        endif
      endif

      if(ic.eq.0) then
        call movemul(rt,s,rto,xr)
        if(idpr.gt.-101) write(lout,*) " no exchanged"
      elseif(ic.eq.1) then
        call movemul(rt,xy,rto,xr)
        if(idpr.gt.-101) write(lout,*) " x-y exchanged"
      elseif(ic.eq.2) then
        call movemul(rt,xz,rto,xr)
        if(idpr.gt.-101) write(lout,*) " x-z exchanged"
      elseif(ic.eq.3) then
        call movemul(rt,yz,rto,xr)
        if(idpr.gt.-101) write(lout,*) " y-z exchanged"
      elseif(ic.eq.4) then
        call movemul(rt,xyz,rto,xr)
        if(idpr.gt.-101) write(lout,*) " x-y-z permuted"
      elseif(ic.eq.5) then
        call movemul(rt,xzy,rto,xr)
        if(idpr.gt.-101) write(lout,*) " x-z-y permuted"
      endif

      do i=1,nd2
        do j=1,nd2
          rt(i,j)=rto(i,j)
        enddo
      enddo

      return
end subroutine movearou

subroutine movemul(rt,xy,rto,xr)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer i,j,k
      real(kind=fPrec) rt,rto,xr,xy
      dimension rt(ndim2,ndim2)
      dimension xy(ndim2,ndim2),rto(ndim2,ndim2)

      do i=1,nd2
        do j=1,nd2
          rto(i,j)=zero
        end do
      end do

      do i=1,nd2
        do j=1,nd2
          do k=1,nd2
            rto(i,k)=xy(i,j)*rt(j,k)+rto(i,k)
          end do
        end do
      end do

      xr=zero
      do i=1,nd2
        do j=1,nd2
          xr=xr+abs(rto(i,j))
        end do
      end do

      do i=1,nd
        xr=xr-abs(rto(2*i-1,2*i-1))
        xr=xr-abs(rto(2*i-1,2*i))
        xr=xr-abs(rto(2*i,2*i))
        xr=xr-abs(rto(2*i,2*i-1))
      enddo

      return
end subroutine movemul

subroutine initpert(st,ang,ra)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nd,nd2,iref,ndc,ndc2,ndpt,ndt,ndim,ndim2,nreso,idsta,ista,angle,dsta,radn,sta,mx,nres
      implicit none
      integer i,nn
      real(kind=fPrec) ang,ra,st
!   X-RATED
!- SETS UP ALL THE COMMON BLOCKS RELEVANT TO NORMAL FORM AND THE BASIS
!- CHANGES INSIDE  MAPNORMF
      dimension st(ndim),ang(ndim),ra(ndim)
      if(iref.gt.0) then
      write(lout,*) iref
      read(iref,*) nres
      if(nres.ge.nreso) then
       write(lout,*) ' NRESO IN LIELIB TOO SMALL '
       write(lout,'(a)') "ERROR 999 in initpert"
       call prror
      endif
      elseif(iref.eq.0) then
      nres=0
      endif
      if(nres.ne.0) write(lout,*)' warning resonances left in the map'
      if(iref.gt.0) then
      do i=1,nres
        read(iref,*) (mx(nn,i),nn=1,nd-ndc)
      enddo
      endif
      do i=nres+1,nreso
        do nn=1,ndim
          mx(nn,i)=0
        enddo
      enddo
!      frank/Etienne
      do i=1,ndim
        angle(i)=zero
        radn(i)=zero
        sta(i)=zero
        dsta(i)=one-sta(i)
        ista(i)=0
        idsta(i)=0
      enddo
      do i=1,nd        !  frank          -ndc
        angle(i)=ang(i)
        radn(i)=ra(i)
        sta(i)=st(i)
        dsta(i)=one-sta(i)
      enddo
      do i=1,nd
        ista(i)=int(sta(i)+c1m2)
        idsta(i)=int(dsta(i)+c1m2)
      enddo
      return
end subroutine initpert

real(kind=fPrec) function dlie(j)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2
      implicit none
      integer i
!      INTEGER J(NTT)
      integer j(*)
      dlie=zero

      do i=1,nd
        dlie=real(j(2*i-1)+j(2*i),fPrec)+dlie
      end do

      dlie=dlie+one
      dlie=one/dlie
      return
end function dlie

real(kind=fPrec) function rext(j)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt,idsta,ista
      implicit none
      integer i,lie,mo
      integer j(*)
      lie=0
      do i=1,nd-ndc
        lie=ista(i)*j(2*i)+lie
      enddo
      mo=mod(lie,4)+1
      rext = zero ! -Wmaybe-uninitialized

      select case (mo)
      case (1)
         rext = one
      case(2)
         rext = -one
      case(3)
         rext = -one
      case(4)
         rext = one
      end select
      return
end function rext

subroutine cpart(h,ch)
      use floatPrecision
      use mathlib_bouncer
      implicit none
      real(kind=fPrec) rext
      external rext
      integer h,ch
      call dacfu(h,rext,ch)
      return
end subroutine cpart

subroutine ctoi(f1,f2)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd2,ndim2
      implicit none
      integer f1,f2
      integer b1(1),x(ndim2)
!
!
      call etallnom(b1(1),1,'B1        ')
      call etallnom(x,nd2  ,'X         ')

      call cpart(f1,b1(1))
      call etctr(x)
      call trx(b1(1),f2,x)
      call dadal(x,nd2)
      call dadal(b1(1),1)
      return
end subroutine ctoi

subroutine itoc(f1,f2)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd2,ndim2
      implicit none
      integer f1,f2
      integer b1(1),x(ndim2)
!
      call etallnom(b1(1),1,'B1        ')
      call etallnom(x,nd2  ,'X         ')

      call etrtc(x)
      call trx(f1,b1(1),x)
      call cpart(b1(1),f2)
      call dadal(x,nd2)
      call dadal(b1(1),1)
      return
end subroutine itoc

subroutine etrtc(x)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndc,ndim2
      implicit none
      integer i
      integer x(*)

      integer rel(ndim2)

      call etallnom(rel,nd2  ,'REL       ')

      call etini(rel)
      call etini(x)
      do i=1,nd-ndc
        call daadd(rel(2*i-1),rel(2*i),x(2*i-1))
        call dasub(rel(2*i-1),rel(2*i),x(2*i))
      enddo
      call dadal(rel,nd2)
      return
end subroutine etrtc

subroutine etctr(x)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndim2
      implicit none
      integer i
      integer x(*)
      integer rel(ndim2)
!
!
      call etallnom(rel,nd2  ,'REL       ')

      call etini(rel)
      call etini(x)
      do i=1,nd-ndc
        call dalin(rel(2*i-1),half,rel(2*i),half,x(2*i-1))
        call dalin(rel(2*i-1),half,rel(2*i),-half,x(2*i))
      enddo
      call dadal(rel,nd2)
      return
end subroutine etctr

subroutine etcjg(x)
      use floatPrecision
      use mathlib_bouncer
      use mod_lie_dab, only : nd,nd2,ndc,ndim2,idsta,ista
      implicit none
      integer i
      integer x(*)

      integer rel(ndim2)
!
!
      call etallnom(rel,nd2  ,'REL       ')

      call etini(rel)
      call etini(x)
      do i=1,nd-ndc
        if(ista(i).eq.1) then
          call dacop(rel(2*i-1),x(2*i))
          call dacop(rel(2*i),x(2*i-1))
        else
          call dacop(rel(2*i-1),x(2*i-1))
          call dacop(rel(2*i),x(2*i))
        endif
      enddo
      call dadal(rel,nd2)
      return
end subroutine etcjg

subroutine eig6(fm,reval,aieval,revec,aievec)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt,ndim2
      implicit none
      integer jet
!**************************************************************************

!  Diagonalization routines of NERI

!ccccccccccccccccc
!
!  this routine finds the eigenvalues and eigenvectors
!  of the full matrix fm.
!  the eigenvectors are normalized so that the real and
!  imaginary part of vectors 1, 3, and 5 have +1 antisymmetric
!  product:
!      revec1 J aivec1 = 1 ; revec3 J aivec3 = 1 ;
!      revec5 J aivec5 = 1.
!  the eigenvectors 2 ,4, and 6 have the opposite normalization.
!  written by F. Neri, Feb 26 1986.
!
      integer nn
      integer ilo,ihi,mdim,info
      real(kind=fPrec) reval(ndim2),aieval(ndim2),                      &
     &revec(ndim2,ndim2),aievec(ndim2,ndim2)
      real(kind=fPrec) fm(ndim2,ndim2),aa(ndim2,ndim2)
      integer i,i1
      real(kind=fPrec) ort(ndim2),vv(ndim2,ndim2)
!  copy matrix to temporary storage (the matrix aa is destroyed)
      do i=1,nd2-ndc2
        do i1=1,nd2-ndc2
          aa(i1,i) = fm(i1,i)
        enddo
      enddo
      ilo = 1
      ihi = nd2-ndc2
      mdim = ndim2
      nn = nd2-ndc2
!  compute eigenvalues and eigenvectors using double
!  precision Eispack routines:
      call ety(mdim,nn,ilo,ihi,aa,ort)
      call etyt(mdim,nn,ilo,ihi,aa,ort,vv)
      call ety2(mdim,nn,ilo,ihi,aa,reval,aieval,vv,info)
      if ( info .ne. 0 ) then
        write(lout,*) '  ERROR IN EIG6'
        return
      endif
!      call neigv(vv,pbkt)
      do i=1,nd-ndc
        do jet=1,nd2-ndc2
          revec(jet,2*i-1)=vv(jet,2*i-1)
          revec(jet,2*i)=vv(jet,2*i-1)
          aievec(jet,2*i-1)=vv(jet,2*i)
          aievec(jet,2*i)=-vv(jet,2*i)
        enddo
      enddo
      do i=1,nd2-ndc2
        if(abs(reval(i)**2+aieval(i)**2 -one).gt.c1m10) then
          write(lout,*) ' EIG6: Eigenvalues off the unit circle!'
        endif
      enddo
      return
end subroutine eig6

subroutine ety(nm,n,low,igh,a,ort)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      real(kind=fPrec) a(nm,n),ort(igh)
      real(kind=fPrec) f,g,h,scale
!
!     this subroutine is a translation of the algol procedure orthes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     orthogonal similarity transformations.
!
!     on input-
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement,
!
!        n is the order of the matrix,
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n,
!
!        a contains the input matrix.
!
!     on output-
!
!        a contains the hessenberg matrix.  information about
!          the orthogonal transformations used in the reduction
!          is stored in the remaining triangle under the
!          hessenberg matrix,
!
!        ort contains further information about the transformations.
!          only elements low through igh are used.
!
!     fortran routine by b. s. garbow
!     modified by filippo neri.
!
!

      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) goto 200

      do 180 m = kp1, la
         h = zero
         ort(m) = zero
         scale = zero
!     ********** scale column (algol tol then not needed) **********
         do i = m, igh
           scale = scale + abs(a(i,m-1))
         end do

         if (scale .eq. zero) goto 180
         mp = m + igh
!     ********** for i=igh step -1 until m do -- **********
         do ii = m, igh
            i = mp - ii
            ort(i) = a(i,m-1) / scale
            h = h + ort(i) * ort(i)
         end do

         g = -sign(sqrt(h),ort(m))
         h = h - ort(m) * g
         ort(m) = ort(m) - g
!     ********** form (i-(u*ut)/h) * a **********
         do j = m, n
            f = zero
!     ********** for i=igh step -1 until m do -- **********
            do ii = m, igh
               i = mp - ii
               f = f + ort(i) * a(i,j)
            end do

            f = f / h

            do i = m, igh
              a(i,j) = a(i,j) - f * ort(i)
            end do
         end do
!     ********** form (i-(u*ut)/h)*a*(i-(u*ut)/h) **********
         do i = 1, igh
            f = zero
!     ********** for j=igh step -1 until m do -- **********
            do jj = m, igh
              j = mp - jj
              f = f + ort(j) * a(i,j)
            end do

            f = f / h

            do j = m, igh
              a(i,j) = a(i,j) - f * ort(j)
            end do
         end do

         ort(m) = scale * ort(m)
         a(m,m-1) = scale * g
  180 continue

  200 return
!     ********** last card of ety **********
end subroutine ety

subroutine etyt(nm,n,low,igh,a,ort,z)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      real(kind=fPrec) a(nm,igh),ort(igh),z(nm,n)
      real(kind=fPrec) g
!
!     this subroutine is a translation of the algol procedure ortrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the orthogonal similarity
!     transformations used in the reduction of a real general
!     matrix to upper hessenberg form by  ety.
!
!     on input-
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement,
!
!        n is the order of the matrix,
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n,
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  orthes
!          in its strict lower triangle,
!
!          ort contains further information about the trans-
!          formations used in the reduction by  ety.
!          only elements low through igh are used.
!
!     on output-
!
!        z contains the transformation matrix produced in the
!          reduction by  ety,
!
!        ort has been altered.
!
!     fortran routine by b. s. garbow.
!     modified by f. neri.


!     ********** initialize z to identity matrix **********
      do i = 1, n

         do j = 1, n
           z(i,j) = zero
         end do
         z(i,i) = one
      end do

      kl = igh - low - 1
      if (kl .lt. 1) goto 200
!     ********** for mp=igh-1 step -1 until low+1 do -- **********
      do 140 mm = 1, kl
         mp = igh - mm
         if (a(mp,mp-1) .eq. zero) goto 140
         mp1 = mp + 1

         do i = mp1, igh
           ort(i) = a(i,mp-1)
         end do

         do j = mp, igh
            g = zero

            do i = mp, igh
              g = g + ort(i) * z(i,j)
            end do
!     ********** divisor below is negative of h formed in orthes.
!                double division avoids possible underflow **********
            g = (g / ort(mp)) / a(mp,mp-1)

            do i = mp, igh
              z(i,j) = z(i,j) + g * ort(i)
            end do
         end do
  140 continue

  200 return
!     ********** last card of etyt **********
end subroutine etyt

subroutine ety2(nm,n,low,igh,h,wr,wi,z,ierr)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,igh,its,low,mp2,enm2,ierr
      real(kind=fPrec) h(nm,n),wr(n),wi(n),z(nm,n)
      real(kind=fPrec) p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,machep
      logical notlas
      real(kind=fPrec) z3r,z3i

      m = 0    ! -Wmaybe-uninitialized
      p = zero ! -Wmaybe-uninitialized
      r = zero ! -Wmaybe-uninitialized
      s = zero ! -Wmaybe-uninitialized
!
!
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to hessenberg form
!     and to accumulate the similarity transformations.
!
!     on input-
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement,
!
!        n is the order of the matrix,
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n,
!
!        h contains the upper hessenberg matrix,
!
!        z contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the hessenberg matrix are desired, z must contain the
!          identity matrix.
!
!     on output-
!
!        h has been destroyed,
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n,
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of z
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of z contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found,
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 200 iterations.
!
!     arithmetic is real(kind=fPrec). complex division
!     is simulated by routin etdiv.
!
!     fortran routine by b. s. garbow.
!     modified by f. neri.
!
!
!     ********** machep is a machine dependent parameter specifying
!                the relative precision of floating point arithmetic.
!
!                **********
      machep = 1.0e-17_fPrec
!     machep = r1mach(4)

      ierr = 0
      norm = zero
      k = 1
!     ********** store roots isolated by balanc
!                and compute matrix norm **********
      do 50 i = 1, n

         do j = k, n
           norm = norm + abs(h(i,j))
         end do

         k = i
         if (i .ge. low .and. i .le. igh) goto 50
         wr(i) = h(i,i)
         wi(i) = zero
   50 continue

      en = igh
      t = zero
!     ********** search for next eigenvalues **********
   60 if (en .lt. low) goto 340
      its = 0
      na = en - 1
      enm2 = na - 1
!     ********** look for single small sub-diagonal element
!                for l=en step -1 until low do -- **********
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) goto 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. zero) s = norm
         if (abs(h(l,l-1)) .le. machep * s) goto 100
   80 continue
!     ********** form shift **********
  100 x = h(en,en)
      if (l .eq. en) goto 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) goto 280
      if (its .eq. 200) goto 1000
      if (its .ne. 10 .and. its .ne. 20) goto 130
!     ********** form exceptional shift **********
      t = t + x
!
      do i = low, en
        h(i,i) = h(i,i) - x
      end do

      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75_fPrec * s
      y = x
      w = -0.4375_fPrec * s * s
  130 its = its + 1
!     ********** look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- **********
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) goto 150
         if (abs(h(m,m-1)) * (abs(q) + abs(r)) .le. machep * abs(p)     &
     &* (abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))) goto 150
  140 continue
!
  150 mp2 = m + 2
!
      do 160 i = mp2, en
         h(i,i-2) = zero
         if (i .eq. mp2) goto 160
         h(i,i-3) = zero
  160 continue
!     ********** double qr step involving rows l to en and
!                columns m to en **********
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) goto 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = zero
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. zero) goto 260
         p = p / x
         q = q / x
         r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) goto 180
         h(k,k-1) = -s * x
         goto 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
!     ********** row modification **********
         do 210 j = k, n
            p = h(k,j) + q * h(k+1,j)
            if (.not. notlas) goto 200
            p = p + r * h(k+2,j)
            h(k+2,j) = h(k+2,j) - p * zz
  200       h(k+1,j) = h(k+1,j) - p * y
            h(k,j) = h(k,j) - p * x
  210    continue
!
         j = min0(en,k+3)
!     ********** column modification **********
         do 230 i = 1, j
            p = x * h(i,k) + y * h(i,k+1)
            if (.not. notlas) goto 220
            p = p + zz * h(i,k+2)
            h(i,k+2) = h(i,k+2) - p * r
  220       h(i,k+1) = h(i,k+1) - p * q
            h(i,k) = h(i,k) - p
  230    continue
!     ********** accumulate transformations **********
         do 250 i = low, igh
            p = x * z(i,k) + y * z(i,k+1)
            if (.not. notlas) goto 240
            p = p + zz * z(i,k+2)
            z(i,k+2) = z(i,k+2) - p * r
  240       z(i,k+1) = z(i,k+1) - p * q
            z(i,k) = z(i,k) - p
  250    continue
!
  260 continue
!
      goto 70
!     ********** one root found **********
  270 h(en,en) = x + t
      wr(en) = h(en,en)
      wi(en) = zero
      en = na
      goto 60
!     ********** two roots found **********
  280 p = (y - x) / two
      q = p * p + w
      zz = sqrt(abs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. zero) goto 320
!     ********** real pair **********
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. zero) wr(en) = x - w / zz
      wi(na) = zero
      wi(en) = zero
      x = h(en,na)
      s = abs(x) + abs(zz)
      p = x / s
      q = zz / s
      r = sqrt(p*p+q*q)
      p = p / r
      q = q / r
!     ********** row modification **********
      do 290 j = na, n
         zz = h(na,j)
         h(na,j) = q * zz + p * h(en,j)
         h(en,j) = q * h(en,j) - p * zz
  290 continue
!     ********** column modification **********
      do 300 i = 1, en
         zz = h(i,na)
         h(i,na) = q * zz + p * h(i,en)
         h(i,en) = q * h(i,en) - p * zz
  300 continue
!     ********** accumulate transformations **********
      do 310 i = low, igh
         zz = z(i,na)
         z(i,na) = q * zz + p * z(i,en)
         z(i,en) = q * z(i,en) - p * zz
  310 continue
!
      goto 330
!     ********** complex pair **********
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      goto 60
!     ********** all roots found.  backsubstitute to find
!                vectors of upper triangular form **********
  340 if (norm .eq. zero) goto 1001
!     ********** for en=n step -1 until 1 do -- **********
      do 800 nn = 1, n
         en = n + 1 - nn
         p = wr(en)
         q = wi(en)
         na = en - 1
         if (q.lt.0) goto 710
         if (q.eq.0) goto 600
         if (q.gt.0) goto 800
!     ********** real vector **********
  600    m = en
         h(en,en) = one
         if (na .eq. 0) goto 800
!     ********** for i=en-1 step -1 until 1 do -- **********
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = h(i,en)
            if (m .gt. na) goto 620
!
            do j = m, na
              r = r + h(i,j) * h(j,en)
            end do

  620       if (wi(i) .ge. zero) goto 630
            zz = w
            s = r
            goto 700
  630       m = i
            if (wi(i) .ne. zero) goto 640
            t = w
            if (w .eq. zero) t = machep * norm
            h(i,en) = -r / t
            goto 700
!     ********** solve real equations **********
  640       x = h(i,i+1)
            y = h(i+1,i)
            q = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i)
            t = (x * s - zz * r) / q
            h(i,en) = t
            if (abs(x) .le. abs(zz)) goto 650
            h(i+1,en) = (-r - w * t) / x
            goto 700
  650       h(i+1,en) = (-s - y * t) / zz
  700    continue
!     ********** end real vector **********
         goto 800
!     ********** complex vector **********
  710    m = na
!     ********** last vector component chosen imaginary so that
!                eigenvector matrix is triangular **********
         if (abs(h(en,na)) .le. abs(h(na,en))) goto 720
         h(na,na) = q / h(en,na)
         h(na,en) = -(h(en,en) - p) / h(en,na)
         goto 730
! 720    z3 = cmplx(0.0,-h(na,en)) / cmplx(h(na,na)-p,q)
!        h(na,na) = real(z3)
!        h(na,en) = aimag(z3)
  720    call etdiv(z3r,z3i,zero,-h(na,en),h(na,na)-p,q)
         h(na,na) = z3r
         h(na,en) = z3i
  730    h(en,na) = zero
         h(en,en) = one
         enm2 = na - 1
         if (enm2 .eq. 0) goto 800
!     ********** for i=en-2 step -1 until 1 do -- **********
         do 790 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = zero
            sa = h(i,en)
!
            do 760 j = m, na
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
!
            if (wi(i) .ge. zero) goto 770
            zz = w
            r = ra
            s = sa
            goto 790
  770       m = i
            if (wi(i) .ne. zero) goto 780
!           z3 = cmplx(-ra,-sa) / cmplx(w,q)
!           h(i,na) = real(z3)
!           h(i,en) = aimag(z3)
            call etdiv(z3r,z3i,-ra,-sa,w,q)
            h(i,na) = z3r
            h(i,en) = z3i
            goto 790
!     ********** solve complex equations **********
  780       x = h(i,i+1)
            y = h(i+1,i)
            vr = (wr(i) - p) * (wr(i) - p) + wi(i) * wi(i) - q * q
            vi = (wr(i) - p) * two * q
            if (vr .eq. zero .and. vi .eq. zero) vr = machep * norm       &
     &* (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz))
!           z3 = cmplx(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra) / cmplx(vr,vi)
!           h(i,na) = real(z3)
!           h(i,en) = aimag(z3)
            call etdiv(z3r,z3i,x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi)
            h(i,na) = z3r
            h(i,en) = z3i
            if (abs(x) .le. abs(zz) + abs(q)) goto 785
            h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
            h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
            goto 790
! 785       z3 = cmplx(-r-y*h(i,na),-s-y*h(i,en)) / cmplx(zz,q)
!           h(i+1,na) = real(z3)
!           h(i+1,en) = aimag(z3)
  785       call etdiv(z3r,z3i,-r-y*h(i,na),-s-y*h(i,en),zz,q)
            h(i+1,na) = z3r
            h(i+1,en) = z3i
  790    continue
!     ********** end complex vector **********
  800 continue
!     ********** end back substitution.
!                vectors of isolated roots **********
      do 840 i = 1, n
         if (i .ge. low .and. i .le. igh) goto 840
!
         do j = i, n
           z(i,j) = h(i,j)
         end do
  840 continue
!     ********** multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- **********
      do jj = low, n
         j = n + low - jj
         m = min0(j,igh)
!
         do i = low, igh
            zz = zero
!
            do k = low, m
              zz = zz + z(i,k) * h(k,j)
            end do

            z(i,j) = zz
        end do
      end do

      goto 1001
!     ********** set error -- no convergence to an
!                eigenvalue after 200 iterations **********
 1000 ierr = en
 1001 return
!     ********** last card of ety2 **********
end subroutine ety2

subroutine etdiv(a,b,c,d,e,f)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
!   computes the complex division
!     a + ib = (c + id)/(e + if)
!  very slow, but tries to be as accurate as
!  possible by changing the order of the
!  operations, so to avoid under(over)flow
!  problems.
!  Written by F. Neri Feb. 12 1986
!
      real(kind=fPrec) a,b,c,d,e,f
      real(kind=fPrec) s,t
      real(kind=fPrec) cc,dd,ee,ff
      real(kind=fPrec) temp
      integer flip
      flip = 0
      cc = c
      dd = d
      ee = e
      ff = f
      if( abs(f).ge.abs(e) ) then
        ee = f
        ff = e
        cc = d
        dd = c
        flip = 1
      endif
      s = one/ee
      t = one/(ee+ ff*(ff*s))
      if ( abs(ff) .ge. abs(s) ) then
        temp = ff
        ff = s
        s = temp
      endif
      if( abs(dd) .ge. abs(s) ) then
        a = t*(cc + s*(dd*ff))
      else if ( abs(dd) .ge. abs(ff) ) then
        a = t*(cc + dd*(s*ff))
      else
        a = t*(cc + ff*(s*dd))
      endif
      if ( abs(cc) .ge. abs(s)) then
        b = t*(dd - s*(cc*ff))
      else if ( abs(cc) .ge. abs(ff)) then
        b = t*(dd - cc*(s*ff))
      else
        b = t*(dd - ff*(s*cc))
      endif
      if (flip.ne.0 ) then
        b = -b
      endif
      return
end subroutine etdiv

subroutine sympl3(m)
!**********************************************************
!
!    SYMPL3
!
!
!   On return ,the matrix m(*,*), supposed to be almost
!   symplectic on entry is made exactly symplectic by
!   using a non iterative, constructive method.
!
!**********************************************************
!
!  Written by F. Neri  Feb 7 1986
!
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
      integer n
      parameter ( n = 3 )
      integer kp,kq,lp,lq,jp,jq,i
      real(kind=fPrec) m(2*n,2*n)
      real(kind=fPrec) qq,pq,qp,pp
!
      do 100 kp=2,2*n,2
        kq = kp-1
        do 200 lp=2,kp-2,2
          lq = lp-1
          qq = zero
          pq = zero
          qp = zero
          pp = zero
          do 300 jp=2,2*n,2
            jq = jp-1
            qq = qq + m(lq,jq)*m(kq,jp) - m(lq,jp)*m(kq,jq)
            pq = pq + m(lp,jq)*m(kq,jp) - m(lp,jp)*m(kq,jq)
            qp = qp + m(lq,jq)*m(kp,jp) - m(lq,jp)*m(kp,jq)
            pp = pp + m(lp,jq)*m(kp,jp) - m(lp,jp)*m(kp,jq)
  300     continue
          do 400 i=1,2*n
            m(kq,i) = m(kq,i) - qq*m(lp,i) + pq*m(lq,i)
            m(kp,i) = m(kp,i) - qp*m(lp,i) + pp*m(lq,i)
  400     continue
  200   continue
        qp = zero
        do 500 jp=2,2*n,2
          jq = jp-1
          qp = qp + m(kq,jq)*m(kp,jp) - m(kq,jp)*m(kp,jq)
  500   continue
        do 600 i=1,2*n
          m(kp,i) = m(kp,i)/qp
  600   continue
!
!  Maybe the following is a better idea ( uses sqrt and is slower )
!       sign = 1.d0
!       if ( qp.lt.0.d0 ) sign = -1.d0
!  OR, BETTER:
!       if ( qp.le.0.d0 ) then complain
!       qp = abs(qp)
!       qp = dsqrt(qp)
!       do 600 i=1,2*n
!         m(kq,i) = m(kq,i)/qp
!         m(kp,i) = sign*m(kp,i)/qp
! 600   continue
  100 continue
      return
end subroutine sympl3

subroutine averaged(f,a,flag,fave)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,no,idpr,ndim2
      implicit none
      integer isi,nord
      real(kind=fPrec) avepol
!      TAKES THE AVERAGE OF A FUNCTION F
!  FLAG TRUE A=ONE TURN MAP
!       FALSE A=A_SCRIPT
!
      integer f,fave,a(*)
      integer cosi(1),sine(1)
      logical flag
      external avepol

      integer a1(ndim2),a2(ndim2),xy(ndim2),hf(ndim2),ftf(ndim2)

      if(.not.flag) then
      call etall(cosi(1),1)
      call etall(sine(1),1)
      call trx(f,f,a)
      call ctor(f,cosi(1),sine(1))
      call dacfu(cosi(1),avepol,fave)
      call dadal(cosi(1),1)
      call dadal(sine(1),1)
      else

      call etall(cosi(1),1)
      call etall(sine(1),1)
      call etall(ftf,nd2)
      call etall(hf,nd2)
      call etall(a2,nd2)
      call etall(a1,nd2)
      call etall(xy,nd2)

      isi=0
      nord=1
      call mapnormf(a,ftf,a2,a1,xy,hf,nord,isi)
      nord=no
      call etcct(a1,a2,xy)
      call facflod(ftf,xy,a1,2,nord,one,-1)
      call trx(f,f,a1)
      call ctor(f,cosi(1),sine(1))
      call dacfu(cosi(1),avepol,fave)

      call dadal(cosi(1),1)
      call dadal(sine(1),1)
      call dadal(ftf,nd2)
      call dadal(hf,nd2)
      call dadal(a2,nd2)
      call dadal(a1,nd2)
      call dadal(xy,nd2)

      endif

      return
end subroutine averaged

real(kind=fPrec) function avepol(j)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt
      implicit none
      integer i
!      INTEGER J(NTT)
      integer j(*)
      avepol=one
      do i=1,(nd-ndc)
        if(j(2*i).ne.j(2*i-1)) then
          avepol=zero
          return
        endif
      enddo

      return
end function avepol

subroutine couplean(map1,tune,map2,oneturn)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use mod_lie_dab, only : nd,nd2,no,nv,ndc,ndc2,ndpt,ndt,ndim,ndim2
      implicit none
      integer i,no1,nord
      real(kind=fPrec) crazy,tpi
!  map1 ascript a1 not there
!  tune 2 or 3 tunes

!   map2 ascript with a couple parameter in nv
!  oneturn map created with tunes and map2
      integer map1(*),oneturn(*),map2(*),ftf(1),hf(1)
      integer xy(ndim2),m1(ndim2),m2(ndim2),a2(ndim2),a1(ndim2)
      integer cs(1),h(1)

      real(kind=fPrec) killnonl,planar,psq(ndim),radsq(ndim)
      real(kind=fPrec) tune(ndim)
      external killnonl,planar

      call etall(ftf(1),1)
      call etall(hf(1),1)
      call etall(a1,nd2)
      call etall(a2,nd2)
      call etall(m1,nd2)
      call etall(m2,nd2)
      call etall(xy,nd2)
      call etall(cs(1),1)
      call etall(h(1),1)

!     map1 is an a-script, the last nv entry should be empty
!  this a-script should around the fixed point to all orders
!     one order is lost because I use PB-field

      tpi=atan_mb(one)*eight
      do i=1,nd2
        call dacfu(map1(i),killnonl,m1(i))
      enddo

      call etini(xy)
      call daclr(cs(1))

      do i=1,nd-ndc
        call dasqr(xy(2*i-1),a2(2*i-1))
        call dasqr(xy(2*i),a2(2*i))
        call daadd(a2(2*i-1),a2(2*i),ftf(1))
        crazy=-tune(i)*tpi/two
        call dacmu(ftf(1),crazy,ftf(1))
        call daadd(ftf(1),cs(1),cs(1))
      enddo

      call etinv(m1,m2)
      call trx(cs(1),h(1),m2)

      call dacfu(h(1),planar,cs(1))
      call dasub(h(1),cs(1),h(1))
      call davar(a2(1),zero,nv)

      call damul(a2(1),h(1),h(1))
      call daadd(cs(1),h(1),h(1))
      call expnd2(h(1),xy,xy,c1m9,1000)

      call dacopd(xy,oneturn)

      nord=1
      call mapnorm(xy,ftf(1),a2,a1,m2,hf(1),nord)

      call gettura(psq,radsq)
      write(lout,*) (psq(i),i=1,nd)

      call etini(xy)
      no1=no
      call fexpo(ftf(1),xy,xy,3,no1,one,-1)
      call etcct(a2,xy,map2)
      call etcct(a1,map2,map2)

      call dadal(ftf(1),1)
      call dadal(hf(1),1)
      call dadal(a1,nd2)
      call dadal(a2,nd2)
      call dadal(m1,nd2)
      call dadal(m2,nd2)
      call dadal(xy,nd2)
      call dadal(cs(1),1)
      call dadal(h(1),1)

      return
end subroutine couplean

real(kind=fPrec) function planar(j)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndc,ndc2,ndpt,ndt
      implicit none
      integer i
!      INTEGER J(NTT)
      integer j(*)
      planar=zero
      do i=1,(nd-ndc)
        if(j(2*i).eq.j(2*i-1)) then
          planar=one
          return
        endif
        if(j(2*i).eq.2) then
          planar=one
          return
        endif
        if(j(2*i-1).eq.2) then
          planar=one
          return
        endif
      enddo

      return
end function planar

real(kind=fPrec) function killnonl(j)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ndc,ndc2,ndpt,ndt
      implicit none
      integer i,ic
!      INTEGER J(NTT)
      integer j(*)

      killnonl=one

      ic=0
      do i=1,nd2-ndc2
        ic=ic+j(i)
      enddo
      if(ic.gt.1) killnonl=zero
      if(j(nv).ne.0) killnonl=zero

      return
end function killnonl

subroutine fexpo1(h,x,w,nrmin,nrmax,sca,ifac)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,ndim2
      implicit none
      integer ifac,nrma,nrmax,nrmi,nrmin
      real(kind=fPrec) sca
      integer x,w,h

      integer v(ndim2)

      nrmi=nrmin-1
      nrma=nrmax-1
      call etall(v,nd2)
      call difd(h,v,-one)
      call facflo(v,x,w,nrmi,nrma,sca,ifac)
      call dadal(v,nd2)

      return
end subroutine fexpo1

subroutine etcctpar(x,ix,xj,z)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use mod_lie_dab, only : nd,nd2,nv,ntt
      implicit none
      integer i,ie,ix
      real(kind=fPrec) xj
      dimension xj(*)
      dimension ie(ntt)
      integer  x(*),z(*)

      call etallnom(ie,nv,'IE        ')
      do i=1,nd2
      call davar(ie(i),zero,i)
      enddo
      do  i=nd2+1,nv
      call dacon(ie(i),xj(i-nd2))
      enddo

      call dacct(x,ix,ie,nv,z,ix)

      call dadal(ie,nv)
      return
end subroutine etcctpar
