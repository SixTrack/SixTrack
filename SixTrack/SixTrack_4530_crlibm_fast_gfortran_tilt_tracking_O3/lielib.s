+cd crlibco
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
+cd crcoall
      integer lout
      common /crflags/lout
+cd choice
      integer nplane
      double precision epsplane,xplane
      common /choice/ xplane(ndim),epsplane,nplane(ndim)
+cd coast
      integer ndc,ndc2,ndpt,ndt
      common /coast/ndc,ndc2,ndt,ndpt
+cd dano
      integer lienot
      common /dano/lienot
+cd dascr
      integer idao,is,iscrri
      double precision rs
      common/dascr/is(100),rs(100),iscrri(100),idao
+cd filtr
      integer ifilt
      common /filtr/ ifilt
+cd ii
      integer nd,nd2,no,nv
      common /ii/no,nv,nd,nd2
+cd integratedex
      double precision xintex
      common /integratedex/ xintex(0:20)
+cd istable
      integer idsta,ista
      common /istable/ista(ndim),idsta(ndim)
+cd printing
      integer idpr
      common /printing/ idpr
+cd resfile
      integer iref
      common /resfile/iref
+cd reson
      integer mx,nres
      common /reson/mx(ndim,nreso),nres
+cd stable
      double precision angle,dsta,rad,sta
      common /stable/sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
+cd tunedef
      integer itu
      common /tunedef/itu
+cd tunerad
      double precision ps,rads
      common /tunerad/ ps(ndim),rads(ndim)
+cd vecflow
      integer iflow,jtune
      common /vecflow/ iflow,jtune
+dk lieinit
      subroutine lieinit(no1,nv1,nd1,ndpt1,iref1,nis)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,iref1,nd1,ndc1,ndim,ndpt1,nis,no1,nv1
      double precision ang,ra,st
!! Lieinit initializes AD Package and Lielib
      parameter (ndim=3)
      dimension st(ndim),ang(ndim),ra(ndim)
+ca ii
+ca coast
+ca resfile
+ca tunedef
+ca dano
+ca printing
+ca integratedex
+ca dascr
+ca choice
!+CA DASCR
      call daexter
      do 1 i=1,ndim
      nplane(i)=2*i-1
      ang(i)=0.d0
      ra(i)=0.d0
 1    st(i)=1.d0
      no=no1
      nv=nv1
      nd=nd1
      nd2=2*nd1
      do 2 i=1,100
 2    is(i)=0
      call daini(no,nv,0)
      if(nis.gt.0)call etallnom(is,nis,'$$IS      ')
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
+if cr
                write(lout,*) ' LETHAL ERROR IN LIEINIT'
+ei
+if .not.cr
                write(*,*) ' LETHAL ERROR IN LIEINIT'
+ei
+if cr
      call abend('                                                  ')
+ei
+if .not.cr
                stop 
+ei
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

+if cr
      if(idpr.eq.1)write(lout,*) ' NO = ',no,' IN DA-CALCULATIONS '
+ei
+if .not.cr
      if(idpr.eq.1)write(*,*) ' NO = ',no,' IN DA-CALCULATIONS '
+ei

      do i=0,20
      xintex(i)=0.d0
      enddo
!hr11 xintex(          0)=       1.000000000000000
      xintex(          0)=       1.000000000000000d0                     !hr11
!hr11 xintex(          1)=  5.000000000000000e-001
      xintex(          1)=  5.000000000000000d-001                       !hr11
!hr11 xintex(          2)=  8.333333333333334e-002
      xintex(          2)=  8.333333333333334d-002                       !hr11
!hr11 xintex(          3)=  0.000000000000000e+000
      xintex(          3)=  0.000000000000000d+000                       !hr11
!hr11 xintex(          4)= -1.388888888888898e-003
      xintex(          4)= -1.388888888888898d-003                       !hr11
!hr11 xintex(          5)=  0.000000000000000e+000
      xintex(          5)=  0.000000000000000d+000                       !hr11
!hr11 xintex(          6)=  3.306878306878064e-005
      xintex(          6)=  3.306878306878064d-005                       !hr11
      xintex(          7)= 0.d0
!hr11 xintex(          8)= -8.267195767165669e-007
      xintex(          8)= -8.267195767165669d-007                       !hr11
      xintex(          9)=  0.d0
!hr11 xintex(         10)=  4.592886537931051e-008
      xintex(         10)=  4.592886537931051d-008                       !hr11
      return
      end
+dk flowpara
      subroutine flowpara(ifl,jtu)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
+ca vecflow
      integer ifl,jtu
      iflow=ifl
      jtune=jtu
      return
      end
+dk pertpeek
      subroutine pertpeek(st,ang,ra)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,nreso,ntt
      double precision ang,ra,st
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      parameter (nreso=20)
      dimension st(ndim),ang(ndim),ra(ndim)
+ca stable
+ca istable
+ca ii
+ca coast
+ca reson
+ca resfile
      do 1 i=1,nd
      st(i)=sta(i)
      ang(i)=angle(i)
      ra(i)=rad(i)
 1    continue
      return
      end
+dk inputres
      subroutine inputres(mx1,nres1)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,j,ndim,ndim2,nreso,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      parameter (nreso=20)
      integer mx1(ndim,nreso),nres1
+ca reson

      nres=nres1
      do 1 i=1,nreso
      do 1 j=1,ndim
      mx(j,i)=0
 1    continue

      do 2 i=1,nres
      do 2 j=1,ndim
      mx(j,i)=mx1(j,i)
 2    continue
      return
      end
+dk respoke
      subroutine respoke(mres,nre,ire)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ire,j,ndim,ndim2,nre,nreso,ntt
      double precision ang,ra,st
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      parameter (nreso=20)
      integer mres(ndim,nreso)
+ca stable
+ca istable
+ca ii
+ca coast
+ca reson
+ca resfile
      dimension ang(ndim),ra(ndim),st(ndim)
      iref=ire
      nres=nre
      do 2 j=1,nreso
      do 1 i=1,nd
      mx(i,j)=mres(i,j)
 1    continue
 2    continue
      call initpert(st,ang,ra)
      return
      end
+dk liepeek
      subroutine liepeek(iia,icoast)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ndim2,nreso,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      parameter (nreso=20)
+ca ii
+ca coast
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
      end
+dk lienot
      subroutine lienot(not)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer no,not

      call danot(not)
      no=not

      return
      end
+dk etallnom
      subroutine etallnom(x,n,nom)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,n,nd2
! CREATES A AD-VARIABLE WHICH CAN BE DESTROYED BY DADAL
! allocates vector of n polynomials and give it the name NOM=A10
      integer x(*),i1(4),i2(4)
      character*10 nom
!hr11 do 1 i=1,iabs(n)
      do 1 i=1,abs(n)                                                    !hr11
 1    x(i)=0
!hr11 call daallno(x,iabs(n),nom)
      call daallno(x,abs(n),nom)                                         !hr11
         if(n.lt.0) then
          call liepeek(i1,i2)
          nd2=i1(4)
           do 2 i=nd2+1,-n
 2         call davar(x(i),0.d0,i)
         endif
      return
      end
+dk etall
      subroutine etall(x,n)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,n,nd2
! allocates vector of n polynomials
      integer x(*),i1(4),i2(4)
      do 1 i=1,iabs(n)
 1    x(i)=0
      call daallno(x,iabs(n),'ETALL     ')
         if(n.lt.0) then
          call liepeek(i1,i2)
          nd2=i1(4)
           do 2 i=nd2+1,-n
 2         call davar(x(i),0.d0,i)
         endif
      return
      end
+dk etall1
      subroutine etall1(x)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer x
!frs
      x=0
!frs
      call daallno(x,1,'ETALL     ')
      return
      end
+dk dadal1
      subroutine dadal1(x)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer x
      call dadal(x,1)
      return
      end
+dk etppulnv
      subroutine etppulnv(x,xi,xff)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer  x(*)
      double precision xi(*),xff(*),xf(ntt),xii(ntt)

      do i=1,nv
      xii(i)=xi(i)
      enddo
      do i=nv+1,ntt
      xii(i)=0.d0
      enddo

      call ppush(x,nv,xii,xf)

      do i=1,nv
      xff(i)=xf(i)
      enddo

      return
      end
+dk etmtree
      subroutine etmtree(y,x)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ie,iv,ndim,ndim2,nt,ntt
! ROUTINES USING THE MAP IN AD-FORM
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie(ntt),iv(ntt)
+ca ii
      integer  x(*),y(*)

      nt=nv-nd2
      if(nt.gt.0) then
      call etallnom(ie,nt,'IE        ')
      do 1 i=nd2+1,nv
      call davar(ie(i-nd2),0.d0,i)
 1    continue
      do 3 i=nd2+1,nv
 3    iv(i)=ie(i-nd2)
      endif
      do 2 i=1,nd2
 2    iv(i)=y(i)
      call mtree(iv,nv,x,nv)
      if(nt.gt.0) then
      call dadal(ie,nt)
      endif
      return
      end
+dk etppush
      subroutine etppush(x,xi)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer  x(*)
      double precision xi(*),xf(ntt),xii(ntt)

      do i=1,nd2
      xii(i)=xi(i)
      enddo

      call ppush(x,nv,xii,xf)

      do i=1,nd2
      xi(i)=xf(i)
      enddo

      return
      end
+dk etppush2
      subroutine etppush2(x,xi,xff)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer  x(*)
      double precision xi(*),xff(*),xf(ntt),xii(ntt)

      do i=1,nd2
      xii(i)=xi(i)
      enddo

      call ppush(x,nv,xii,xf)

      do i=1,nd2
      xff(i)=xf(i)
      enddo

      return
      end
+dk ppushlnv
      subroutine ppushlnv(x,xi,xff,nd1)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,nd1,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer  x(*)
      double precision xi(*),xff(*),xf(ntt),xii(ntt)

      do i=1,nd1
      xii(i)=xi(i)
      enddo
      do i=nd1+1,ntt
      xii(i)=0.d0
      enddo

      call ppush(x,nv,xii,xf)

      do i=1,nd1
      xff(i)=xf(i)
      enddo

      return
      end
+dk etcct
      subroutine etcct(x,y,z)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ie,iv,ndim,ndim2,nt,ntt
!  Z=XoY
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie(ntt),iv(ntt)
+ca ii
      integer  x(*),y(*),z(*)

      nt=nv-nd2
      if(nt.gt.0) then
      call etallnom(ie,nt,'IE        ')
      do 1 i=nd2+1,nv
      call davar(ie(i-nd2),0.d0,i)
 1    continue
      do 3 i=nd2+1,nv
 3    iv(i)=ie(i-nd2)
      endif
      do 2 i=1,nd2
 2    iv(i)=y(i)
      call dacct(x,nd2,iv,nv,z,nd2)
      if(nt.gt.0) then
      call dadal(ie,nt)
      endif
      return
      end
+dk trx
      subroutine trx(h,rh,y)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ie,iv,ndim,ndim2,nt,ntt
!  :RH: = Y :H: Y^-1 =  :HoY:
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie(ntt),iv(ntt)
+ca ii
      integer h,rh
      integer y(*)
!
      nt=nv-nd2
      if(nt.gt.0) then
      call etallnom(ie,nt,'IE        ')
      do 1 i=nd2+1,nv
      call davar(ie(i-nd2),0.d0,i)
 1    continue
      do 3 i=nd2+1,nv
 3    iv(i)=ie(i-nd2)
      endif
      do 2 i=1,nd2
 2    iv(i)=y(i)
      call dacct(h,1,iv,nv,rh,1)
      if(nt.gt.0) then
      call dadal(ie,nt)
      endif
      return
      end
+dk trxflo
      subroutine trxflo(h,rh,y)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer j,k,ndim,ndim2,ntt
!  *RH* = Y *H* Y^-1  CHANGE OF A VECTOR FLOW OPERATOR
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer h(*),rh(*),y(*)
      integer yi(ndim2),ht(ndim2),b1,b2
!
!
      call etallnom(yi,nd2  ,'YI        ')
      call etallnom(ht,nd2  ,'HT        ')
      call etallnom(b1,1,'B1        ')
      call etallnom(b2,1,'B2        ')

      call etinv(y,yi)
!----- HT= H o Y
      call etcct(h,y,ht)
!----
      call daclrd(rh)
      do 1 j=1,nd2
      do 2 k=1,nd2
      call dader(k,yi(j),b1)
      call trx(b1,b2,y)
        call damul(b2,ht(k),b1)
        call daadd(b1,rh(j),b2)
        call dacop(b2,rh(j))
 2    continue
 1    continue

      call dadal(b2,1)
      call dadal(b1,1)
      call dadal(ht,nd2)
      call dadal(yi,nd2)
      return
      end
+dk simil
      subroutine simil(a,x,ai,y)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ndim2,ntt
!  Y= AoXoAI
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
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
      end
+dk etini
      subroutine etini(x)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
!  X=IDENTITY
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer x(*)
!*DAEXT(NO,NV) X(NDIM2)
      do 1 i=1,nd2
      call davar(x(i),0.d0,i)
 1    continue
      return
      end
+dk etinv
      subroutine etinv(x,y)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ie1,ie2,iv1,iv2,ndim,ndim2,nt,ntt
! Y=X^-1
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie1(ntt),ie2(ntt),iv1(ntt),iv2(ntt)
+ca ii

      integer x(*),y(*)

      nt=nv-nd2
      if(nt.gt.0) then
       do 100 i=1,nt
       ie1(i)=0
 100   ie2(i)=0
       call etallnom(ie1,nt,'IE1       ')
       call etallnom(ie2,nt,'IE2       ')
       do 1 i=nd2+1,nv
       call davar(ie1(i-nd2),0.d0,i)
 1     continue
       do 3 i=nd2+1,nv
       iv1(i)=ie1(i-nd2)
 3     iv2(i)=ie2(i-nd2)
      endif
      do 2 i=1,nd2
      iv1(i)=x(i)
 2    iv2(i)=y(i)

      call dainv(iv1,nv,iv2,nv)
      if(nt.gt.0) then
      call dadal(ie2,nt)
      call dadal(ie1,nt)
      endif
      return
      end
+dk etpin
      subroutine etpin(x,y,jj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ie1,ie2,iv1,iv2,jj,ndim,ndim2,nt,ntt
!  Y=PARTIAL INVERSION OF X SEE BERZ'S PACKAGE
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension ie1(ntt),ie2(ntt),iv1(ntt),iv2(ntt),jj(*)
+ca ii

      integer x(*),y(*)

      nt=nv-nd2
      if(nt.gt.0) then
       do 100 i=1,nt
       ie1(i)=0
 100   ie2(i)=0
       call etallnom(ie1,nt,'IE1       ')
       call etallnom(ie2,nt,'IE2       ')
       do 1 i=nd2+1,nv
       call davar(ie1(i-nd2),0.d0,i)
 1     continue
       do 3 i=nd2+1,nv
       iv1(i)=ie1(i-nd2)
 3     iv2(i)=ie2(i-nd2)
      endif
      do 2 i=1,nd2
      iv1(i)=x(i)
 2    iv2(i)=y(i)

      call dapin(iv1,nv,iv2,nv,jj)
      if(nt.gt.0) then
      call dadal(ie2,nt)
      call dadal(ie1,nt)
      endif
      return
      end
+dk dapek0
      subroutine dapek0(v,x,jj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,jj,ndim2,ntt
      double precision x
!- MORE EXTENSIONS OF BASIC BERZ'S PACKAGE
      parameter (ndim2=6)
      parameter (ntt=40)
      integer v(*),jd(ntt)
      dimension x(*)
      do 2 i=1,ntt
 2    jd(i)=0
      do 1 i=1,jj
 1    call dapek(v(i),jd,x(i))
      return
      end
+dk dapok0
      subroutine dapok0(v,x,jj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,jj,ndim2,ntt
      double precision x
      parameter (ndim2=6)
      parameter (ntt=40)
      integer v(*),jd(ntt)
      dimension x(*)
      do 2 i=1,ntt
 2    jd(i)=0
      do 1 i=1,jj
 1    call dapok(v(i),jd,x(i))
      return
      end
+dk dapokzer
      subroutine dapokzer(v,jj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,jj,ndim2,ntt
      parameter (ndim2=6)
      parameter (ntt=40)
      integer v(*),jd(ntt)
      do 2 i=1,ntt
 2    jd(i)=0
      do 1 i=1,jj
 1    call dapok(v(i),jd,0.d0)
      return
      end
+dk davar0
      subroutine davar0(v,x,jj)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,jj,ndim2,ntt
      double precision x
      parameter (ndim2=6)
      parameter (ntt=40)
      integer v(*)
      dimension x(*)
      do 1 i=1,jj
 1    call davar(v(i),x(i),i)
      return
      end
+dk comcfu
      subroutine comcfu(b,f1,f2,c)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      double precision f1,f2
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
      end
+dk take
      subroutine take(h,m,ht)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,m,ndim,ntt
      double precision r
!  HT= H_M  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
      parameter (ndim=3)
      parameter (ntt=40)
+ca ii
      integer h,ht,j(ntt)

      integer b1,b2,b3
!
!
      call etallnom(b1,1,'B1        ')
      call etallnom(b2,1,'B2        ')
      call etallnom(b3,1,'B3        ')

      if(no.ge.2) then
       if(m.eq.0) then
       do i=1,ntt
       j(i)=0
       enddo
        call dapek(h,j,r)
        call dacon(ht,r)
       else
      call danot(m)
      call dacop(h,b1)
      call danot(m-1)
      call dacop(b1,b2)
      call danot(no)
      call dasub(b1,b2,b3)
      call dacop(b3,ht)
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
        call dapok(b3,j,r)
       j(i)=0
       enddo
       call dacop(b3,ht)
       else
       call daclr(ht)
       endif
      endif

      call dadal(b3,1)
      call dadal(b2,1)
      call dadal(b1,1)
      return
      end
+dk taked
      subroutine taked(h,m,ht)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,m,ndim2,ntt
!  \VEC{HT}= \VEC{H_M}  (TAKES M^th DEGREE PIECE ALL VARIABLES INCLUDED)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer h(*),ht(*),j(ntt)

      integer b1,b2,x(ndim2)
!
      call etallnom(b1,1,'B1        ')
      call etallnom(b2,1,'B2        ')
      call etallnom(x,nd2  ,'X         ')


      do i=1,ntt
      j(i)=0
      enddo

       do   i=1,nd2
      call take(h(i),m,ht(i))
      enddo
      call dadal(x,nd2)
      call dadal(b2,1)
      call dadal(b1,1)
      return
      end
+dk daclrd
      subroutine daclrd(h)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim2,ntt
! clear a map : a vector of nd2 polynomials
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer h(*)
      do 1 i=1,nd2
 1    call daclr(h(i))
      return
      end
+dk dacopd
      subroutine dacopd(h,ht)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim2,ntt
!    H goes into HT  (nd2 array)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer h(*),ht(*)
      do 1 i=1,nd2
 1    call dacop(h(i),ht(i))
      return
      end
+dk dacmud
      subroutine dacmud(h,sca,ht)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim2,ntt
      double precision sca
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer h(*),ht(*)
      do 1 i=1,nd2
 1    call dacmu(h(i),sca,ht(i))
      return
      end
+dk dalind
      subroutine dalind(h,rh,ht,rt,hr)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim2
      double precision rh,rt
      parameter (ndim2=6)
+ca ii
      integer h(*),ht(*),hr(*)

      integer b(ndim2)
!
      call etallnom(b,nd2  ,'B         ')

      do 1 i=1,nd2
 1    call dalin(h(i),rh,ht(i),rt,b(i))
      call dacopd(b,hr)
      call dadal(b,nd2)
      return
      end
+dk daread
      subroutine daread(h,nd1,mfile,xipo)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,mfile,nd1,ndim2,ntt
      double precision rx,xipo
!  read a map
      parameter (ndim2=6)
      parameter (ntt=40)
      integer h(*),j(ntt)
      do 11 i=1,ntt
 11   j(i)=0
      do 1 i=1,nd1
      call darea(h(i),mfile)
      call dapek(h(i),j,rx)
      rx=rx*xipo
      call dapok(h(i),j,rx)
 1    continue
      return
      end
+dk daprid
      subroutine daprid(h,n1,n2,mfile)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,mfile,n1,n2,ndim2,ntt
!  print a map
      parameter (ndim2=6)
      parameter (ntt=40)
      integer  h(*)
      if(mfile.le.0) return
      do 1 i=n1,n2
 1    call dapri(h(i),mfile)
      return
      end
+dk prresflo
      subroutine prresflo(h,eps,mfile)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,mfile,ndim2,ntt
      double precision deps,eps,filtres
!  print a map   in resonance basis for human consumption (useless)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer  h(*) ,b(ndim2) ,c(ndim2)
+ca ii
+ca filtr
      external filtres
      call etall(b,nd2)
      call etall(c,nd2)
      call dacopd(h,c)
      do i=1,nd2
      ifilt=(-1)**i
      call  dacfu(c(i),filtres,h(i))
      enddo

      deps=-1.d0
      call daeps(deps)
      call daeps(eps)

      call dacopd(c,h)
      call daeps(deps)
      call  dadal(c,nd2)
      call  dadal(b,nd2)
      return
      end
+dk filtres
      double precision function filtres(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      INTEGER J(NTT)
      integer j(*)
+ca ii
+ca coast
+ca filtr
      filtres=1.d0
      ic=0
      do i=1,(nd2-ndc2)
        ic=ic+j(i)*(-1)**(i+1)
      enddo
      ic=ic+ifilt
      if(ic.lt.0) filtres=0.d0
      if(ic.eq.0.and.ifilt.eq.1) then
        filtres=0.0d0
      endif
      return
      end
+dk daflo
      subroutine daflo(h,x,y)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
! LIE EXPONENT ROUTINES WITH FLOW OPERATORS

!     \VEC{H}.GRAD X =Y
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer  h(*),x,y
      integer b1,b2,b3
!
      call etallnom(b1,1,'B1        ')
      call etallnom(b2,1,'B2        ')
      call etallnom(b3,1,'B3        ')

      call daclr(b1)
      call daclr(b2)
      do 1 i=1,nd2
      call dader(i,x,b2)
      call damul(b2,h(i),b3)
      call daadd(b3,b1,b2)
      call dacop(b2,b1)
 1    continue
      call dacop(b1,y)
      call dadal(b3,1)
      call dadal(b2,1)
      call dadal(b1,1)
      return
      end
+dk daflod
      subroutine daflod(h,x,y)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer  h(*),x(*),y(*)
      integer b1(ndim2),b2(ndim2)
!
      call etall(b1,nd2)
      call etall(b2,nd2)

      call dacopd(h,b1)
      call dacopd(x,b2)

      do 1 i=1,nd2
      call daflo(b1,b2(i),y(i))
 1    continue

      call dadal(b1,nd2)
      call dadal(b2,nd2)
      return
      end
+dk intd
      subroutine intd(v,h,sca)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      double precision dlie,sca
! IF SCA=-1.D0
!     \VEC{V}.GRAD   = J GRAD H . GRAD = :H:

! IF SCA=1.D0
!     \VEC{V}.GRAD  = GRAD H . GRAD
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      external dlie
      integer v(*),h

      integer b1,b2,b3,b4,x(ndim2)
!
!
      call etallnom(b1,1,'B1        ')
      call etallnom(b2,1,'B2        ')
      call etallnom(b3,1,'B3        ')
      call etallnom(b4,1,'B4        ')
      call etallnom(x,nd2  ,'X         ')

      call daclr(b4)
      call daclr(h)
      call etini(x)
      do 1 i=1,nd
      call dacfu(v(2*i-1),dlie,b3)
      call dacfu(v(2*i),dlie,b1)
      call damul(b1,x(2*i-1),b2)
      call damul(b3,x(2*i),b1)
      call dalin(b2,1.d0,b1,sca,b3)
      call daadd(b3,b4,b2)
      call dacop(b2,b4)
 1    continue
      call dacop(b4,h)
      call dadal(x,nd2)
      call dadal(b4,1)
      call dadal(b3,1)
      call dadal(b2,1)
      call dadal(b1,1)
      return
      end
+dk difd
      subroutine difd(h1,v,sca)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      double precision sca
! INVERSE OF INTD ROUTINE
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer  v(*),h1
      integer b1,h
      call etall(b1,1)
      call etall(h,1)
      call dacop(h1,h)
      do 1 i=1,nd
      call dader(2*i-1,h,v(2*i))
      call dader(2*i,h,b1)
      call   dacmu(b1,sca,v(2*i-1))
 1    continue
      call dadal(h,1)
      call dadal(b1,1)
      return
      end
+dk expflo
      subroutine expflo(h,x,y,eps,nrmax)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,nrmax,ntt
      double precision coe,eps,r,rbefore
! DOES EXP( \VEC{H} ) X = Y
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca printing
      integer h(*),x,y
      integer b1,b2,b3,b4
      logical more
!
!
      call etallnom(b1,1,'B1        ')
      call etallnom(b2,1,'B2        ')
      call etallnom(b3,1,'B3        ')
      call etallnom(b4,1,'B4        ')

      call dacop(x,b4)
      call dacop(x,b1)
      more=.true.
      rbefore=1.d30
      do 1 i=1,nrmax
      coe=1.d0/dble(i)
      call dacmu(b1,coe,b2)
      call daflo(h,b2,b1)
      call daadd(b4,b1,b3)
      call daabs(b1,r)
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
        call dacop(b3,y)
        call dadal(b4,1)
        call dadal(b3,1)
        call dadal(b2,1)
        call dadal(b1,1)
        return
        endif
       rbefore=r
       endif
100   continue
      call dacop(b3,b4)
 1    continue
      if(idpr.ge.0) then
+if cr
      write(lout,*) ' NORM ',eps,' NEVER REACHED IN EXPFLO '
+ei
+if .not.cr
      write(*,*) ' NORM ',eps,' NEVER REACHED IN EXPFLO '
+ei
+if cr
      write(lout,*) 'NEW IDPR '
+ei
+if .not.cr
      write(*,*) 'NEW IDPR '
+ei
      read(5,*) idpr
      endif
      call dacop(b3,y)
      call dadal(b4,1)
      call dadal(b3,1)
      call dadal(b2,1)
      call dadal(b1,1)
      return
      end
+dk expflod
      subroutine expflod(h,x,w,eps,nrmax)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer j,ndim,ndim2,nrmax,ntt
      double precision eps
! DOES EXP( \VEC{H} ) \VEC{X} = \VEC{Y}
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer x(*),w(*),h(*)
      integer b0,v(ndim2)
!
!
      call etallnom(b0,1,'B0        ')
      call etallnom(v,nd2  ,'V         ')

      call dacopd(x,v)
      do 333 j=1,nd2
      call expflo(h,v(j),b0,eps,nrmax)
 333  call dacop(b0,v(j))
      call dacopd(v,w)
      call dadal(v,nd2)
      call dadal(b0,1)
      return
      end
+dk facflo
      subroutine facflo(h,x,w,nrmin,nrmax,sca,ifac)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ifac,ndim,ndim2,nmax,nrmax,nrmin,ntt
      double precision eps,sca
! IFAC=1
! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX ) X= Y
! IFAC=-1
! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) X= Y
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer x,w,h(*)
      integer bm(ndim2),b0(ndim2),v
!
      call etallnom(bm,nd2  ,'BM        ')
      call etallnom(b0,nd2  ,'B0        ')
      call etallnom(v,1  ,'V         ')

      call dacop(x,v)

      eps=-1.d0
      call daeps(eps)
      nmax=100
!
! IFAC =1 ---> V = EXP(:SCA*H(NRMAX):)...EXP(:SCA*H(NRMIN):)X
      if(ifac.eq.1) then
         do 1 i=nrmax,nrmin,-1
         call taked(h,i,b0)
         call dacmud(b0,sca,bm)

         call expflo(bm,v,b0(1),eps,nmax)
         call dacop(b0(1),v)
 1       continue
      else
! IFAC =-1 ---> V = EXP(:SCA*H(NRMIN):)...EXP(:SCA*H(NRMAX):)X
      do 11 i=nrmin,nrmax
      call taked(h,i,b0)
      call dacmud(b0,sca,bm)

         call expflo(bm,v,b0(1),eps,nmax)
         call dacop(b0(1),v)
 11   continue
      endif
      call dacop(v,w)
      call dadal(v,1)
      call dadal(b0,nd2)
      call dadal(bm,nd2)
      return
      end
+dk facflod
      subroutine facflod(h,x,w,nrmin,nrmax,sca,ifac)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ifac,ndim,ndim2,nrmax,nrmin,ntt
      double precision sca
! IFAC=1
! DOES EXP(SCA \VEC{H}_MRMIN ) ... EXP(SCA \VEC{H}_NRMAX )  \VEC{X}= \VEC{Y}
! IFAC=-1
! DOES EXP(SCA \VEC{H}_NRMAX ) ... EXP(SCA \VEC{H}_MRMIN ) \VEC{X}= \VEC{Y}
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer x(*),w(*),h(*)

      do 1 i=1,nd2
      call facflo(h,x(i),w(i),nrmin,nrmax,sca,ifac)
 1    continue

      return
      end
+dk fexpo
      subroutine fexpo(h,x,w,nrmin,nrmax,sca,ifac)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ifac,ndim,ndim2,nrma,nrmax,nrmi,nrmin,ntt
      double precision sca
!   WRAPPED ROUTINES FOR THE OPERATOR  \VEC{H}=:H:
! WRAPPING FACFLOD
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer x(*),w(*),h

      integer v(ndim2)

      nrmi=nrmin-1
      nrma=nrmax-1
      call etall(v,nd2)
      call difd(h,v,-1.d0)
      call facflod(v,x,w,nrmi,nrma,sca,ifac)

      call dadal(v,nd2)

      return
      end
+dk etcom
      subroutine etcom(x,y,h)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,j,ndim,ndim2,ntt
! ETCOM TAKES THE BRACKET OF TWO VECTOR FIELDS.
      parameter (ndim2=6)
      parameter (ndim=3)
      parameter (ntt=40)
+ca ii
      integer h(*),x(*),y(*),t1,t2,t3(ndim2)

      call etall(t1,1)
      call etall(t2,1)
      call etall(t3,nd2)

      do 2 j=1,nd2
      do 1 i=1,nd2

      call dader(i,x(j),t1)
      call dader(i,y(j),t2)
      call damul(x(i),t2,t2)
      call damul(y(i),t1,t1)
      call dalin(t2,1.d0,t1,-1.d0,t1)
      call daadd(t1,t3(j),t3(j))

 1    continue
 2    continue

      call dacopd(t3,h)

      call dadal(t1,1)
      call dadal(t2,1)
      call dadal(t3,nd2)
      return
      end
+dk etpoi
      subroutine etpoi(x,y,h)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
! ETPOI TAKES THE POISSON BRACKET OF TWO FUNCTIONS
      parameter (ndim2=6)
      parameter (ndim=3)
      parameter (ntt=40)
+ca ii
      integer h,x,y,t1,t2,t3

      call etall(t1,1)
      call etall(t2,1)
      call etall(t3,1)

      do 1 i=1,nd

      call dader(2*i-1,x,t1)
      call dader(2*i,y,t2)
      call damul(t1,t2,t1)

      call dalin(t1,1.d0,t3,1.d0,t3)
      call dader(2*i-1,y,t1)
      call dader(2*i,x,t2)
      call damul(t1,t2,t1)

      call dalin(t1,-1.d0,t3,1.d0,t3)

 1    continue

      call dacop(t3,h)

      call dadal(t1,1)
      call dadal(t2,1)
      call dadal(t3,1)
      return
      end
+dk exp1d
      subroutine exp1d(h,x,y,eps,non)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ndim2,non,ntt
      double precision eps
! WRAPPING EXPFLO
      parameter (ndim2=6)
      parameter (ndim=3)
      parameter (ntt=40)
+ca ii
+ca printing
      integer h,x,y

      integer v(ndim2)

      call etall(v,nd2)
      call difd(h,v,-1.d0)
      call expflo(v,x,y,eps,non)

      call dadal(v,nd2)

      return
      end
+dk expnd2
      subroutine expnd2(h,x,w,eps,nrmax)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer j,ndim,ndim2,nrmax,ntt
      double precision eps
! WRAPPING EXPFLOD USING EXP1D
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer x(*),w(*),h

      integer b0,v(ndim2)
!
!
      call etallnom(b0,1,'B0        ')
      call etallnom(v,nd2  ,'V         ')

        call dacopd(x,v)
      do 333 j=1,nd2
      call exp1d(h,v(j),b0,eps,nrmax)
 333  call dacop(b0,v(j))
      call dacopd(v,w)
      call dadal(v,nd2)
      call dadal(b0,1)
      return
      end
+dk flofacg
      subroutine flofacg(xy,h,epsone)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,k,kk,ndim,ndim2,nrmax,ntt
      double precision eps,epsone,r,xn,xnbefore,xnorm,xnorm1,xx
! GENERAL ONE EXPONENT FACTORIZATION
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca printing
      logical more
+ca ii
+ca integratedex
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
      xnorm1=0.d0
      do i=1,nd2
        call daabs(xy(i),r)
      xnorm1=xnorm1+r
       enddo
      xnbefore=1.d36
      more=.false.
      eps=1.e-9
      nrmax=1000
      xn=10000.d0
      do 333 k=1,nrmax
      call dacmud(h,-1.d0,t)
      call expflod(t,xy,x,eps,nrmax)
      call dalind(x,1.d0,v,-1.d0,t)
! write(20,*) "$$$$$$$$$$$$$$",k,"$$$$$$$$$$$$$$$$$$$$"
! call daprid(t,1,1,20)
       if(xn.lt.epsone) then
+if cr
            if(idpr.ge.0) write(lout,*) "xn quadratic",xn
+ei
+if .not.cr
            if(idpr.ge.0) write(*,*) "xn quadratic",xn
+ei
        call daflod(t,t,w)
        call dalind(t,1.d0,w,-0.5d0,t)
        call dacopd(t,z)
        call dacopd(t,w)
!  second order in W
        call etcom(h,w,x)
           call etcom(x,w,x)
!  END OF  order in W

              do kk=1,10
        call etcom(h,w,w)
        call dalind(z,1.d0,w,xintex(kk),z)
              enddo
        call dacopd(z,t)
      xx=1.d0/12.d0
      call dalind(x,xx,h,1.d0,h)
       endif

      call dalind(t,1.d0,h,1.d0,h)
      xnorm=0.d0
      do i=1,nd2
        call daabs(t(i),r)
      xnorm=xnorm+r
       enddo
      xn=xnorm/xnorm1
+if cr
      if(xn.ge.epsone.and.(idpr.ge.0)) write(lout,*)" xn linear ",xn
+ei
+if .not.cr
      if(xn.ge.epsone.and.(idpr.ge.0)) write(*,*)" xn linear ",xn
+ei
      if(xn.lt.eps.or.more) then
      more=.true.
      if(xn.ge.xnbefore) goto 1000
      xnbefore=xn
      endif
 333  continue
+if cr
1000  write(lout,*) " iteration " , k
+ei
+if .not.cr
1000  write(*,*) " iteration " , k
+ei
      call dadal(x,nd2)
      call dadal(w,nd2)
      call dadal(v,nd2)
      call dadal(t,nd2)
      call dadal(z,nd2)
      return
      end
+dk flofac
      subroutine flofac(xy,x,h)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer k,ndim,ndim2,ntt
! GENERAL DRAGT-FINN FACTORIZATION
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
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
      call dalind(v,1.d0,h,1.d0,h)
      call facflod(h,w,v,k,k,-1.d0,-1)
      call dacopd(v,w)
 333  continue
      call dadal(w,nd2)
      call dadal(v,nd2)
      return
      end
+dk liefact
      subroutine liefact(xy,x,h)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ndim2,ntt
! SYMPLECTIC DRAGT-FINN FACTORIZATION WRAPPING FLOFAC
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer xy(*),x(*),h

      integer v(ndim2)

      call etall(v,nd2)

      call flofac(xy,x,v)
      call intd(v,h,-1.d0)
!
      call dadal(v,nd2)

      return
      end
+dk mapnorm
      subroutine mapnorm(x,ft,a2,a1,xy,h,nord)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer isi,ndim,ndim2,nord,ntt
!--NORMALIZATION ROUTINES OF LIELIB
!- WRAPPING MAPNORMF
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer x(*),a1(*),a2(*),ft,xy(*),h,hf(ndim2),ftf(ndim2)

      call etall(ftf,nd2)
      call etall(hf,nd2)
      isi=0
      call mapnormf(x,ftf,a2,a1,xy,hf,nord,isi)
      call intd(hf,h,-1.d0)
      call intd(ftf,ft,-1.d0)
      call dadal(ftf,nd2)
      call dadal(hf,nd2)

      return
      end
+dk gettura
      subroutine gettura(psq,radsq)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ik,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      double precision psq(ndim),radsq(ndim)
+ca ii
+ca tunerad

      do ik=1,nd
      psq(ik)=ps(ik)
      radsq(ik)=rads(ik)
      enddo

      return
      end
+dk setidpr
      subroutine setidpr(idprint,nplan)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idprint,ik,ndim,ndim2,nplan
      parameter (ndim=3)
      parameter (ndim2=6)
       dimension nplan(ndim)
+ca ii
+ca printing
+ca choice

      do ik=1,nd
      nplane(ik)=nplan(ik)
      enddo
      idpr=idprint

      return
      end
+dk idprset
      subroutine idprset(idprint)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer idprint,ndim,ndim2
      parameter (ndim=3)
      parameter (ndim2=6)
+ca ii
+ca printing
+ca choice

      idpr=idprint

      return
      end
+dk mapnormf
      subroutine mapnormf(x,ft,a2,a1,xy,h,nord,isi)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ij,isi,ndim,ndim2,nord,ntt
      double precision angle,p,rad,st,x2pi,x2pii
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension angle(ndim),st(ndim),p(ndim),rad(ndim)
+ca ii
+ca coast
      integer x(*),a1(*),a2(*),ft(*),xy(*),h(*)
+ca tunedef
+ca printing
+ca vecflow
+ca tunerad
      integer a1i(ndim2),a2i(ndim2)
!
      call etallnom(a1i,nd2  ,'A1I       ')
      call etallnom(a2i,nd2  ,'A2I       ')
!     frank/etienne
      do itu=1,ndim
        angle(itu)=0.d0
        p(itu)=0.d0
        st(itu)=0.d0
        rad(itu)=0.d0
        ps(itu)=0.d0
        rads(itu)=0.d0
      enddo
      jtune=isi
+if crlibm
!hr11 x2pii=1.d0/atan_rn(1.d0)/8.d0
      x2pii=(1.d0/atan_rn(1.d0))/8.d0                                    !hr11
+ei
+if .not.crlibm
!hr11 x2pii=1.d0/atan(1.d0)/8.d0
      x2pii=(1.d0/atan(1.d0))/8.d0                                       !hr11
+ei
+if crlibm
      x2pi=atan_rn(1.d0)*8.d0
+ei
+if .not.crlibm
      x2pi=atan(1.d0)*8.d0
+ei
      call dacopd(x,xy)
! goto fix point in the parameters + pt to order nord>=1
      call gofix(xy,a1,a1i,nord)
      call simil(a1i,xy,a1,xy)
! linear part
      call midbflo(xy,a2,a2i,angle,rad,st)
      do ij=1,nd-ndc
        p(ij)=angle(ij)*(st(ij)*(x2pii-1.d0)+1.d0)
      enddo
      if(ndc.eq.1) p(nd)=angle(nd)
      if(idpr.ge.0) then
+if cr
        write(lout,*) 'tune    ',(p(ij),ij=1,nd)
+ei
+if .not.cr
        write(*,*) 'tune    ',(p(ij),ij=1,nd)
+ei
+if cr
        write(lout,*) 'damping ', (rad(ij),ij=1,nd)
+ei
+if .not.cr
        write(*,*) 'damping ', (rad(ij),ij=1,nd)
+ei
      endif
      do ij=1,nd       !  -ndc    frank
        ps(ij)=p(ij)
        rads(ij)=rad(ij)
      enddo
      call initpert(st,angle,rad)
      call simil(a2i,xy,a2,xy)
      call dacopd(xy,a2i)
!        write(6,*) 'Entering orderflo'
      call orderflo(h,ft,xy,angle,rad)
      do ij=1,nd-ndc
        p(ij)=angle(ij)
        if(angle(ij).gt.x2pi/2.d0.and.st(ij).gt.0.d0.and.itu.eq.1)then
          p(ij)=angle(ij)-x2pi
+if cr
          write(lout,*) ij,' TH TUNE MODIFIED IN H2 TO ',p(ij)/x2pi
+ei
+if .not.cr
          write(*,*) ij,' TH TUNE MODIFIED IN H2 TO ',p(ij)/x2pi
+ei
        endif
      enddo
      call h2pluflo(h,p,rad)
!      CALL TAKED(A2I,1,XY)
      call taked(a2i,1,a1i)
      call etcct(xy,a1i,xy)
      
      call dadal(a2i,nd2)
      call dadal(a1i,nd2)
      return
      end
+dk gofix
      subroutine gofix(xy,a1,a1i,nord)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,nord,ntt
      double precision xic
! GETTING TO THE FIXED POINT AND CHANGING TIME APPROPRIATELY IN THE
! COASTING BEAM CASE

!****************************************************************
! X = A1 XY A1I WHERE X IS TO THE FIXED POINT TO ORDER NORD
! for ndpt not zero, works in all cases. (coasting beam: eigenvalue
!1 in Jordan form)
!****************************************************************
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
+ca coast
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
        call dalin(x(i),1.d0,rel(i),-1.d0,v(i))
      enddo
      call etinv(v,w)
      call daclrd(x)
      if(ndc.eq.1) then
        call davar(x(ndpt),0.d0,ndpt)
      endif
      call etcct(w,x,v)
      if(ndc.eq.1) then
        call daclr(v(nd2))
        call daclr(v(nd2-ndc))
      endif
      call dalind(rel,1.d0,v,1.d0,a1)
      call dalind(rel,1.d0,v,-1.d0,a1i)

      if(ndpt.ne.0) then

!  CORRECTIONS
        call daclrd(w)
        call daclrd(v)
        call daclrd(x)

        do i=1,nd2-ndc2
          call dalin(a1(i),1.d0,rel(i),-1.d0,w(i))
        enddo

!      COMPUTE Deta/Ddelta
        call dacopd(w,a1)

        do i=1,nd2-ndc2
          call dader(ndpt,w(i),w(i))
        enddo
!      COMPUTE J*Deta/dDELTA

        do i=1,nd-ndc
          call dacmu(w(2*i),1.d0,v(2*i-1) )
          call dacmu(w(2*i-1),-1.d0,v(2*i) )
        enddo

!hr11   xic=(-1)**(ndt)
        xic=dble((-1)**(ndt))                                            !hr11

        do i=1,nd2-ndc2
          call damul(v(i),rel(i),x(1))
          call daadd(x(1),w(ndt),w(ndt))
          call dacop(a1(i),w(i))
        enddo
        call dacmu(w(ndt),xic,w(ndt))

        call expflod(w,rel,a1,1.d-7,10000)
! END OF  CORRECTIONS

        call etinv(a1,a1i)

      endif

      call danot(no)

      call dadal(rel,nd2)
      call dadal(v,nd2)
      call dadal(w,nd2)
      call dadal(x,nd2)
      return
      end
+dk transver
      double precision function transver(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,ndim
! USED IN A DACFU CALL OF GOFIX
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      INTEGER J(NTT)
      integer j(*)
+ca ii
+ca coast

      transver=1.d0
      ic=0
      do i=1,nd2-ndc2
        ic=ic+j(i)
      enddo
      if(ic.ne.1) transver=0.d0
      return
      end
+dk orderflo
      subroutine orderflo(h,ft,x,ang,ra)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer k,ndim,ndim2,ntt
      double precision ang,ra
!-   NONLINEAR NORMALIZATION PIECE OF MAPNORMF
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca printing
+ca ii
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
      call facflod(h,x,v,2,k-1,-1.d0,-1)
! EXTRACTING K TH DEGREE OF V ----> W
      call taked(v,k,w)
!  write(16,*) "$$$$$$$$  K  $$$$$$$$$$", k
! W = EXP(B5) + ...
       call dacopd(w,b5)
!      CALL INTD(W,B5,-1.D0)
! B5 ON EXIT IS THE NEW CONTRIBUTION TO H
! B6 IS THE NEW CONTRIBUTION TO FT
       call nuanaflo(b5,b6)
      call dalind(b5,1.d0,h,1.d0,b1)
      call dacopd(b1,h)
! EXP(B9) = EXP( : ROTI B6 :)
      call trxflo(b6,b9,roi)

! V = EXP(-B6) REL
      call facflod(b6,rel,v,k,k,-1.d0,1)
! W = V o X
      call etcct(v,x,w)
      if(idpr.ge.0) then
+if cr
        write(lout,*) ' ORDERFLO K= ', k
+ei
+if .not.cr
        write(*,*) ' ORDERFLO K= ', k
+ei
      endif
! X = EXP(B9) W
      call facflod(b9,w,x,k,k,1.d0,1)
! B6 IS THE NEW CONTRIBUTION TO FT
      call dalind(b6,1.d0,ft,1.d0,b1)
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
      end
+dk nuanaflo
      subroutine nuanaflo(h,ft)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      double precision dfilt,filt,xgam,xgbm
! RESONANCE DENOMINATOR OPERATOR (1-R^-1)^-1
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
+ca vecflow
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
      end
+dk xgam
      double precision function xgam(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,ij,ik,ndim,ndim2
      double precision ad,ans,as,ex,exh
! XGAM AND XGBM ARE THE EIGENVALUES OF THE OPERATOR NEWANAFLO
      parameter (ndim=3)
      parameter (ndim2=6)
!      PARAMETER (NTT=40)
+ca vecflow
+ca stable
+ca istable
+ca ii
+ca coast
!      INTEGER J(NTT),JJ(NDIM),JP(NDIM)
      integer j(*),jj(ndim),jp(ndim)
      xgam=0.d0
      ad=0.d0
      as=0.d0
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
!hr11   ic=ic+iabs(jj(i))
        ic=ic+abs(jj(i))                                                 !hr11
      enddo

      do i=1,nd-ndc
!hr11   ad=dsta(i)*dble(jj(i))*angle(i)-dble(jp(i))*rad(i)+ad
        ad=((dsta(i)*dble(jj(i)))*angle(i)-dble(jp(i))*rad(i))+ad        !hr11
!hr11   as=sta(i)*dble(jj(i))*angle(i)+as
        as=(sta(i)*dble(jj(i)))*angle(i)+as                              !hr11
      enddo

+if crlibm
      exh=exp_rn(ad/2.d0)
+ei
+if .not.crlibm
      exh=exp(ad/2.d0)
+ei
      ex=exh**2
+if crlibm
!hr11 ans=4.d0*ex*(sinh_rn(ad/2.d0)**2+sin_rn(as/2.d0)**2)
      ans=(4.d0*ex)*(sinh_rn(ad/2.d0)**2+sin_rn(as/2.d0)**2)             !hr11
+ei
+if .not.crlibm
!hr11 ans=4.d0*ex*(dsinh(ad/2.d0)**2+sin(as/2.d0)**2)
      ans=(4.d0*ex)*(sinh(ad/2.d0)**2+sin(as/2.d0)**2)                   !hr11
+ei
+if crlibm
!hr11 xgam=2.d0*(-exh*sinh_rn(ad/2.d0)+ex*sin_rn(as/2.d0)**2)/ans
      xgam=(2.d0*(ex*sin_rn(as/2.d0)**2-exh*sinh_rn(ad/2.d0)))/ans       !hr11
+ei
+if .not.crlibm
!hr11 xgam=2.d0*(-exh*dsinh(ad/2.d0)+ex*sin(as/2.d0)**2)/ans
      xgam=(2.d0*(ex*sin(as/2.d0)**2-exh*sinh(ad/2.d0)))/ans            !hr11
+ei

      return
      end
+dk xgbm
      double precision function xgbm(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,ij,ik,ndim,ndim2
      double precision ad,ans,as,ex,exh
      parameter (ndim=3)
      parameter (ndim2=6)
!      PARAMETER (NTT=40)
+ca vecflow
+ca stable
+ca istable
+ca ii
+ca coast
!      INTEGER J(NTT),JJ(NDIM),JP(NDIM)
      integer j(*),jj(ndim),jp(ndim)
      xgbm=0.d0
      ad=0.d0
      as=0.d0
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
!hr11   ad=dsta(i)*dble(jj(i))*angle(i)-dble(jp(i))*rad(i)+ad
        ad=((dsta(i)*dble(jj(i)))*angle(i)-dble(jp(i))*rad(i))+ad        !hr11
!hr11   as=sta(i)*dble(jj(i))*angle(i)+as
        as=(sta(i)*dble(jj(i)))*angle(i)+as                              !hr11
      enddo

+if crlibm
      exh=exp_rn(ad/2.d0)
+ei
+if .not.crlibm
      exh=exp(ad/2.d0)
+ei
      ex=exh**2
+if crlibm
!hr11 ans=4.d0*ex*(sinh_rn(ad/2.d0)**2+sin_rn(as/2.d0)**2)
      ans=(4.d0*ex)*(sinh_rn(ad/2.d0)**2+sin_rn(as/2.d0)**2)             !hr11
+ei
+if .not.crlibm
!hr11 ans=4.d0*ex*(dsinh(ad/2.d0)**2+sin(as/2.d0)**2)
      ans=(4.d0*ex)*(sinh(ad/2.d0)**2+sin(as/2.d0)**2)                  !hr11
+ei
+if crlibm
!hr11 xgbm=sin_rn(as)*ex/ans
      xgbm=(sin_rn(as)*ex)/ans                                           !hr11
+ei
+if .not.crlibm
!hr11 xgbm=sin(as)*ex/ans
      xgbm=(sin(as)*ex)/ans                                              !hr11
+ei

      return
      end
+dk filt
      double precision function filt(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,ic1,ic2,ij,ik,ji,ndim,ndim2,nreso
!  PROJECTION FUNCTIONS ON THE KERNEL ANMD RANGE OF (1-R^-1)
!-  THE KERNEL OF (1-R^-1)
      parameter (ndim=3)
      parameter (ndim2=6)
!      PARAMETER (NTT=40)
      parameter (nreso=20)
+ca stable
+ca istable
+ca ii
+ca coast
+ca reson
+ca vecflow
!      INTEGER J(NTT),JJ(NDIM)
      integer j(*),jj(ndim)

      filt=1.d0

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

      filt=0.d0
      return
       end
+dk dfilt
      double precision function dfilt(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ndim2,nreso
      double precision fil,filt
!-  THE RANGE OF (1-R^-1)^1
!- CALLS FILT AND EXCHANGES 1 INTO 0 AND 0 INTO 1.
      parameter (ndim=3)
      parameter (ndim2=6)
!      PARAMETER (NTT=40)
      parameter (nreso=20)
+ca ii
+ca coast
+ca reson
      external filt
!      INTEGER J(NTT)
      integer j(*)

      fil=filt(j)
      if(fil.gt.0.5d0) then
       dfilt=0.d0
      else
       dfilt=1.d0
      endif
      return
       end
+dk dhdjflo
      subroutine dhdjflo(h,t)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      double precision coe,x2pi
! CONVENIENT TUNE SHIFT FINDED FOR SYMPLECTIC CASE (NU,DL)(H)=T
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
+ca coast
      integer h(*),t(*)

      integer b1(ndim2),b2(ndim2),bb1,bb2
!
      call etall(b1,nd2)
      call etall(b2,nd2)
      call etall(bb1,1)
      call etall(bb2,1)

+if crlibm
      x2pi=atan_rn(1.d0)*8.d0
+ei
+if .not.crlibm
      x2pi=atan(1.d0)*8.d0
+ei
      call ctorflo(h,b1,b2)
      coe=1.d0/x2pi

      do i=1,nd-ndc
        call datra(2*i,b2(2*i),bb1)
        call dacmu(bb1,coe,t(i+nd))
        call dacop(t(i+nd),bb1)
        call daclr(bb2)
        call rtoc(bb1,bb2,bb1)
        call dacop(bb1,t(i))
      enddo

      if(ndpt.ne.0) then
        call dacop(h(ndt),t(nd))
        call dacop(b1(ndt),t(nd2))
      endif

      call dadal(bb2,1)
      call dadal(bb1,1)
      call dadal(b2,nd2)
      call dadal(b1,nd2)
      return
      end
+dk dhdj
      subroutine dhdj(h,t)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      double precision coe,x2pi
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
+ca coast
      integer h,t(*)

      integer b1,b2,bb1,bb2
!
      call etallnom(b1,1,'B1        ')
      call etallnom(b2,1,'B2        ')
      call etallnom(bb1,1,'BB1       ')
      call etallnom(bb2,1,'BB2       ')

+if crlibm
      x2pi=atan_rn(1.d0)*8.d0
+ei
+if .not.crlibm
      x2pi=atan(1.d0)*8.d0
+ei
      call ctor(h,b1,b2)
      coe=-2.d0/x2pi
      do i=1,nd-ndc
        call dader(2*i-1,b1,b2)
        call datra(2*i,b2,bb2)
        call dacmu(bb2,coe,t(i+nd))
        call dacop(t(i+nd),bb2)
        call daclr(b2)
        call rtoc(bb2,b2,bb1)
        call dacop(bb1,t(i))
      enddo

      if(ndpt.eq.nd2) then
        call dader(ndpt,h,t(nd))
        call dader(ndpt,b1,t(nd2))
        call dacmu(t(nd),-1.d0,t(nd))
        call dacmu(t(nd2),-1.d0,t(nd2))
      endif
      if(ndt.eq.nd2) then
        call dader(ndpt,h,t(nd))
        call dader(ndpt,b1,t(nd2))
      endif
      call dadal(bb2,1)
      call dadal(bb1,1)
      call dadal(b2,1)
      call dadal(b1,1)
      return
      end
+dk h2pluflo
      subroutine h2pluflo(h,ang,ra)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,j,ndim,ndim2,ntt
      double precision ang,r1,r2,ra,st
! POKES IN \VEC{H}  ANGLES AND DAMPING COEFFFICIENTS
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca stable
      dimension ang(ndim),st(ndim),ra(ndim),j(ntt)
+ca ii
+ca coast
      integer h(*)
!*DAEXT(NO,NV) H

      do i=1,nd
        st(i)=2.d0*sta(i)-1.d0
      enddo

      do i=1,ntt
        j(i)=0
      enddo

      do i=1,nd-ndc
        j(2*i-1)=1
!hr11   r1=-ang(i)
        r1=-1d0*ang(i)                                                   !hr11
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
!hr11   call dapok(h(ndt),j,-ang(nd))
        call dapok(h(ndt),j,-1d0*ang(nd))                                !hr11
      endif
      return
      end
+dk rotflo
      subroutine rotflo(ro,ang,ra)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      double precision ang,ch,co,ra,sh,si,sim,xx
! CREATES R AND R^-1 USING THE EXISTING ANGLES AND DAMPING
! COULD BE REPLACED BY A CALL H2PLUFLO FOLLOWED BY EXPFLOD
! CREATES R
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca istable
      dimension co(ndim),si(ndim),ang(ndim),ra(ndim)
      integer j(ntt)
+ca ii
+ca coast
      integer ro(*)
      call daclrd(ro)
      do i=1,nd-ndc
+if crlibm
        xx=exp_rn(ra(i))
+ei
+if .not.crlibm
        xx=exp(ra(i))
+ei
        if(ista(i).eq.0) then
          call hyper(ang(i),ch,sh)
          co(i)=ch*xx
!hr11     si(i)=-sh*xx
          si(i)=(-1d0*sh)*xx                                             !hr11
        else
+if crlibm
          co(i)=cos_rn(ang(i))*xx
+ei
+if .not.crlibm
          co(i)=cos(ang(i))*xx
+ei
+if crlibm
          si(i)=sin_rn(ang(i))*xx
+ei
+if .not.crlibm
          si(i)=sin(ang(i))*xx
+ei
        endif
      enddo
      do i=1,nd-ndc
        if(ista(i).eq.0)then
          sim=si(i)
        else
!hr11     sim=-si(i)
          sim=-1d0*si(i)                                                 !hr11
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
        call dapok(ro(ndt),j,1.d0)
        call dapok(ro(ndpt),j,0.d0)
        j(ndt)=0
        j(ndpt)=1
        call dapok(ro(ndt),j,ang(nd))
        call dapok(ro(ndpt),j,1.d0)
        j(ndpt)=0
      endif

      return
      end
+dk rotiflo
      subroutine rotiflo(roi,ang,ra)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      double precision ang,ch,co,ra,sh,si,sim,simv,xx
! CREATES  R^-1
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca istable
      dimension co(ndim),si(ndim),ang(ndim),ra(ndim)
      integer j(ntt)
+ca ii
+ca coast
      integer roi(*)

      do i=1,10
        j(i)=0
      enddo

      call daclrd(roi)
      do i=1,nd-ndc
+if crlibm
!hr11   xx=exp_rn(-ra(i))
        xx=exp_rn(-1d0*ra(i))                                            !hr11
+ei
+if .not.crlibm
!hr11   xx=exp(-ra(i))
        xx=exp(-1d0*ra(i))                                               !hr11
+ei
        if(ista(i).eq.0) then
          call hyper(ang(i),ch,sh)
          co(i)=ch*xx
!hr11     si(i)=-sh*xx
          si(i)=(-1d0*sh)*xx
        else
+if crlibm
          co(i)=cos_rn(ang(i))*xx
+ei
+if .not.crlibm
          co(i)=cos(ang(i))*xx
+ei
+if crlibm
          si(i)=sin_rn(ang(i))*xx
+ei
+if .not.crlibm
          si(i)=sin(ang(i))*xx
+ei
        endif
      enddo
      do i=1,nd-ndc
        if(ista(i).eq.0)then
          sim=si(i)
        else
!hr11     sim=-si(i)
          sim=-1d0*si(i)                                                 !hr11
        endif
        j(2*i-1)=1
        call dapok(roi(2*i-1),j,co(i))
!hr11   simv=-sim
        simv=-1d0*sim                                                    !hr11
        call dapok(roi(2*i),j,simv)
        j(2*i-1)=0
        j(2*i)=1
!hr11   simv=-si(i)
        simv=-1d0*si(i)                                                  !hr11
        call dapok(roi(2*i),j,co(i))
        call dapok(roi(2*i-1),j,simv)
        j(2*i)=0
      enddo

      if(ndc.eq.1) then
        j(ndt)=1
        call dapok(roi(ndt),j,1.d0)
        call dapok(roi(ndpt),j,0.d0)
        j(ndt)=0
        j(ndpt)=1
        call dapok(roi(ndt),j,-ang(nd))
        call dapok(roi(ndpt),j,1.d0)
        j(ndpt)=0
      endif

      return
      end
+dk hyper
      subroutine hyper(a,ch,sh)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      double precision a,ch,sh,x,xi
!   USED IN ROTIFLO AND ROTFLO
+if crlibm
      x=exp_rn(a)
+ei
+if .not.crlibm
      x=exp(a)
+ei
      xi=1.d0/x
      ch=(x+xi)/2.d0
      sh=(x-xi)/2.d0
      return
      end
+dk ctor
      subroutine ctor(c1,r2,i2)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim2,ntt
! CHANGES OF BASIS
!   C1------> R2+I R1
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer c1,r2,i2
      integer b1,b2,x(ndim2)
!
!
      call etallnom(b1,1,'B1        ')
      call etallnom(b2,1,'B2        ')
      call etallnom(x,nd2  ,'X         ')

      call ctoi(c1,b1)
      call etcjg(x)
      call trx(b1,b2,x)
      call dalin(b1,.5d0,b2,.5d0,r2)
      call dalin(b1,.5d0,b2,-.5d0,i2)
      call dadal(x,nd2)
      call dadal(b2,1)
      call dadal(b1,1)
      return
      end
+dk rtoc
      subroutine rtoc(r1,i1,c2)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim2,ntt
!  INVERSE OF CTOR
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer c2,r1,i1

      integer b1
!
      call etallnom(b1,1,'B1        ')

      call daadd(r1,i1,b1)
      call itoc(b1,c2)
      call dadal(b1,1)
      return
      end
+dk ctorflo
      subroutine ctorflo(c,dr,di)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ndim2,ntt
! FLOW CTOR
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      integer dr(*),di(*),c(*)

      call ctord(c,dr,di)
      call resvec(dr,di,dr,di)

      return
      end
+dk rtocflo
      subroutine rtocflo(dr,di,c)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ndim2,ntt
! FLOW RTOC
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer dr(*),di(*),c(*),er(ndim2),ei(ndim2)

      call etall(er,nd2)
      call etall(ei,nd2)

      call reelflo(dr,di,er,ei)
      call rtocd(er,ei,c)

      call dadal(er,nd2)
      call dadal(ei,nd2)

      return
      end
+dk ctord
      subroutine ctord(c,cr,ci)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
! ROUTINES USED IN THE INTERMEDIATE STEPS OF CTORFLO AND RTOCFLO
! SAME AS CTOR  OVER ARRAYS CONTAINING ND2 COMPONENTS
! ROUTINE USEFUL IN INTERMEDIATE FLOW CHANGE OF BASIS
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer c(*),ci(*),cr(*)
      do 1 i=1,nd2
      call ctor(c(i),cr(i),ci(i))
 1    continue
      return
      end
+dk rtocd
      subroutine rtocd(cr,ci,c)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
!  INVERSE OF CTORD
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer c(*),ci(*),cr(*)
      do 1 i=1,nd2
      call rtoc(cr(i),ci(i),c(i))
 1    continue
      return
      end
+dk resvec
      subroutine resvec(cr,ci,dr,di)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
! DOES THE SPINOR PART IN CTORFLO
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca stable
+ca istable
+ca coast
+ca ii
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
      end
+dk reelflo
      subroutine reelflo(c,ci,f,fi)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
! DOES THE SPINOR PART IN RTOCFLO
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca stable
+ca istable
+ca coast
+ca ii
      integer c(*),ci(*),f(*),fi(*),e(ndim2),ei(ndim2)

      call etall(e,nd2)
      call etall(ei,nd2)

      do i=1,nd-ndc
        call dalin(c(2*i-1),0.5d0,c(2*i),0.5d0,e(2*i-1))
        call dalin(ci(2*i-1),0.5d0,ci(2*i),0.5d0,ei(2*i-1))
        if(ista(i).eq.1) then
          call dalin(ci(2*i-1),0.5d0,ci(2*i),-0.5d0,e(2*i))
          call dalin(c(2*i-1),-0.5d0,c(2*i),0.5d0,ei(2*i))
        else
          call dalin(ci(2*i-1),0.5d0,ci(2*i),-0.5d0,ei(2*i))
          call dalin(c(2*i-1),0.5d0,c(2*i),-0.5d0,e(2*i))
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
      end
+dk compcjg
      subroutine compcjg(cr,ci,dr,di)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ndim2,ntt
! TAKES THE COMPLEX CONJUGATE IN RESONANCE BASIS OF A POLYNOMIAL
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer dr,di,ci,cr,x(ndim2)

      call etall(x,nd2)

      call etcjg(x)
      call trx(cr,dr,x)
      call trx(ci,di,x)
      call dacmu(di,-1.d0,di)

      call dadal(x,nd2)
      return
      end
+dk midbflo
      subroutine midbflo(c,a2,a2i,q,a,st)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,j,ndim,ndim2,ntt
      double precision a,ch,cm,cr,q,r,sa,sai,shm,                       &
     &st,x2pi
! LINEAR EXACT NORMALIZATION USING EIGENVALUE PACKAGE OF NERI
      parameter (ntt=40)
      parameter (ndim2=6)
      parameter (ndim=3)
      integer jx(ntt)
+ca ii
+ca coast
      dimension cr(ndim2,ndim2),st(ndim),q(ndim),a(ndim)
      dimension sa(ndim2,ndim2),sai(ndim2,ndim2),cm(ndim2,ndim2)
      integer c(*),a2(*),a2i(*)
!*DAEXT(NO,NV) C(NDIM2),A2(NDIM2),A2I(NDIM2)
+if crlibm
      x2pi=atan_rn(1.d0)*8.d0
+ei
+if .not.crlibm
      x2pi=atan(1.d0)*8.d0
+ei

      do i=1,ntt
        jx(i)=0
      enddo

!     frank/etienne
      do i=1,ndim
        st(i)=0d0
        q(i)=0d0
        a(i)=0d0
      enddo
!     frank/etienne
      do i=1,ndim2
!     frank/etienne
        do j=1,ndim2
          sai(i,j)=0.d0
          sa(i,j)=0.d0
          cm(i,j)=0.d0
          cr(i,j)=0.d0
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
!hr11   if(st(i)+0.001.gt.1.d0) then
        if(st(i)+0.001d0.gt.1.d0) then                                   !hr11
!hr11     a(i)=dsqrt(cr(2*i-1,2*i-1)**2+cr(2*i-1,2*i)**2)
          a(i)=sqrt(cr(2*i-1,2*i-1)**2+cr(2*i-1,2*i)**2)                 !hr11
+if crlibm
          q(i)=acos_rn(cr(2*i-1,2*i-1)/a(i))
+ei
+if .not.crlibm
          q(i)=acos(cr(2*i-1,2*i-1)/a(i))
+ei
+if crlibm
          a(i)=log_rn(a(i))
+ei
+if .not.crlibm
          a(i)=log(a(i))
+ei
          if(cr(2*i-1,2*i).lt.0.d0) q(i)=x2pi-q(i)
        else
!hr11     a(i)=dsqrt(cr(2*i-1,2*i-1)**2-cr(2*i-1,2*i)**2)
          a(i)=sqrt(cr(2*i-1,2*i-1)**2-cr(2*i-1,2*i)**2)                 !hr11
          ch=cr(2*i-1,2*i-1)/a(i)
          shm=cr(2*i-1,2*i)/a(i)
!       CH=CH+DSQRT(CH**2-1.D0)
!       q(i)=DLOG(CH)
+if crlibm
!hr11     q(i)=-log_rn(ch+shm)
          q(i)=-1d0*log_rn(ch+shm)                                       !hr11
+ei
+if .not.crlibm
!hr11     q(i)=-log(ch+shm)
          q(i)=-1d0*log(ch+shm)                                          !hr11
+ei
!       IF(cr(2*i-1,2*i).gt.0.d0) Q(I)=-Q(I)
+if crlibm
          a(i)=log_rn(a(i))
+ei
+if .not.crlibm
          a(i)=log(a(i))
+ei
        endif
      enddo

      if(ndc.eq.0) then
!hr11   if(st(3)+0.001.gt.1.d0.and.nd.eq.3.and.q(nd).gt.0.5d0)          &
        if(st(3)+0.001d0.gt.1.d0.and.nd.eq.3.and.q(nd).gt.0.5d0)        &!hr11
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
          if(r.ne.0.d0)call  dapok(a2(i),jx,r)
          jx(j)=1
          r=sai(i,j)
          if(r.ne.0.d0)call  dapok(a2i(i),jx,r)
          jx(j)=0
        enddo
      enddo

      return
      end
+dk mapflol
      subroutine mapflol(sa,sai,cr,cm,st)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ier,iunst,j,l,n,n1,ndim,ndim2
      double precision ap,ax,cm,cr,                                     &
     &p,rd,rd1,ri,rr,s1,sa,sai,st,vi,vr,w,x,x2pi,xd,xj,xsu,xx
      parameter (ndim2=6)
      parameter (ndim=3)
+ca ii
+ca coast
!---- FROM TRACKING CODE
+ca printing
+ca choice
! ---------------------
      dimension cr(ndim2,ndim2),xj(ndim2,ndim2),n(ndim),x(ndim)
      dimension rr(ndim2),ri(ndim2),sa(ndim2,ndim2),xx(ndim)            &
     &,sai(ndim2,ndim2),cm(ndim2,ndim2),w(ndim2,ndim2),st(ndim)
      dimension vr(ndim2,ndim2),vi(ndim2,ndim2),s1(ndim2,ndim2),p(ndim2)
+if crlibm
      x2pi=atan_rn(1.d0)*8.d0
+ei
+if .not.crlibm
      x2pi=atan(1.d0)*8.d0
+ei
      n1=0
!     frank/etienne
      do i=1,ndim2
        do j=1,ndim2
          cr(j,i)=cm(i,j)
          xj(i,j)=0.d0
          s1(i,j)=0.d0
        enddo
      enddo
!     frank/etienne
      do i=1,ndim
        n(i)=0
        xj(2*i-1,2*i)=1.d0
        xj(2*i,2*i-1)=-1.d0
      enddo
!     frank/etienne
      do i=1,ndim2
        do j=1,ndim2
          sai(i,j)=0.d0
          w(i,j)=cm(i,j)
        enddo
      enddo
      if(ndc.eq.1) then
        s1(nd2-ndc,nd2-ndc)=1.d0
        s1(nd2,nd2)=1.d0
        sai(nd2-ndc,nd2-ndc)=1.d0
        sai(nd2,nd2)=1.d0
      endif
      call mulnd2(xj,w)
      call mulnd2(cr,w)
      if(idpr.ge.0.or.idpr.eq.-102) then
+if cr
        write(lout,*)                                                   &
     &'Check of the symplectic condition on the linear part'
+ei
+if .not.cr
        write(*,*)'Check of the symplectic condition on the linear part'
+ei
        xsu=0.d0
        do i=1,nd2
+if cr
          write(lout,'(6(2x,g23.16))') ( w(i,j), j = 1, nd2 )
+ei
+if .not.cr
          write(*,'(6(2x,g23.16))') ( w(i,j), j = 1, nd2 )
+ei
          do j=1,nd2
            xsu=xsu+abs(w(i,j))
          enddo
        enddo
+if cr
        write(lout,*)'deviation for symplecticity ',                    &
     &100.d0*(xsu-nd2)/xsu,
+ei
+if .not.cr
        write(*,*)'deviation for symplecticity ',100.d0*(xsu-nd2)/xsu,  &
+ei
     &' %'
+if debug
!       call warr('symplcdev',100.d0*(xsu-nd2)/xsu,0,0,0,0)
+ei
      endif
      call eig6(cr,rr,ri,vr,vi)
      if(idpr.ge.0) then
+if cr
        write(lout,*) '   '
+ei
+if .not.cr
        write(*,*) '   '
+ei
+if cr
        write(lout,*) '       Index         Real Part  ',
+ei
+if .not.cr
        write(*,*) '       Index         Real Part  ',                  &
+ei
     &'       ArcSin(Imaginary Part)/2/pi'
+if cr
        write(lout,*) '   '
+ei
+if .not.cr
        write(*,*) '   '
+ei
        do i=1,nd-ndc
          rd1=dsqrt(rr(2*i-1)**2+ri(2*i-1)**2)
          rd=dsqrt(rr(2*i)**2+ri(2*i)**2)
+if crlibm
+if cr
          write(lout,*) 2*i-1,rr(2*i-1),asin_rn(ri(2*i-1)/rd1)/x2pi
+ei
+if .not.cr
          write(*,*) 2*i-1,rr(2*i-1),asin_rn(ri(2*i-1)/rd1)/x2pi
+ei
+ei
+if .not.crlibm
+if cr
          write(lout,*) 2*i-1,rr(2*i-1),asin(ri(2*i-1)/rd1)/x2pi
+ei
+if .not.cr
          write(*,*) 2*i-1,rr(2*i-1),asin(ri(2*i-1)/rd1)/x2pi
+ei
+ei
+if crlibm
+if cr
          write(lout,*) 2*i,rr(2*i),asin_rn(ri(2*i)/rd)/x2pi
+ei
+if .not.cr
          write(*,*) 2*i,rr(2*i),asin_rn(ri(2*i)/rd)/x2pi
+ei
+ei
+if .not.crlibm
+if cr
          write(lout,*) 2*i,rr(2*i),asin(ri(2*i)/rd)/x2pi
+ei
+if .not.cr
          write(*,*) 2*i,rr(2*i),asin(ri(2*i)/rd)/x2pi
+ei
+ei
+if crlibm
+if cr
          write(lout,*) ' alphas ', log_rn(dsqrt(rd*rd1))
+ei
+if .not.cr
          write(*,*) ' alphas ', log_rn(dsqrt(rd*rd1))
+ei
+ei
+if .not.crlibm
+if cr
          write(lout,*) ' alphas ', log(dsqrt(rd*rd1))
+ei
+if .not.cr
          write(*,*) ' alphas ', log(dsqrt(rd*rd1))
+ei
+ei
        enddo
+if cr
        write(lout,*)
+ei
+if .not.cr
        write(*,*)                                                      &
+ei
     &' select ',nd-ndc,                                                &
     &' eigenplanes (odd integers <0 real axis)'
        read(5,*) (n(i),i=1,nd-ndc)
      elseif(idpr.eq.-100) then
        do i=1,nd-ndc
          n(i)=nplane(i)
        enddo
      elseif(idpr.eq.-101.or.idpr.eq.-102) then
        do i=1,nd-ndc
          if(ri(2*i).ne.0.d0) then
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
          st(i)=0.d0
          iunst=1
        else
          st(i)=1.d0
        endif
        x(i)=0.d0
        xx(i)=1.d0
        do j=1,nd-ndc
          x(i)=vr(2*j-1,n(i))*vi(2*j,n(i))-vr(2*j,n(i))*vi(2*j-1,n(i))+ &
     &x(i)
        enddo
      enddo

      do i=1,nd-ndc
        if(x(i).lt.0.d0) xx(i)=-1.d0
        x(i)=dsqrt(abs(x(i)))
      enddo
      do i=1,nd2-ndc2
        do j=1,nd-ndc
          if(st(j)+0.001.gt.1.d0) then
            sai(2*j-1,i)=vr(i,n(j))*xx(j)/x(j)
            sai(2*j,i)=vi(i,n(j))/x(j)
          else
            ax=vr(i,n(j))*xx(j)/x(j)
            ap=vi(i,n(j))/x(j)
            sai(2*j-1,i)=(ax+ap)/dsqrt(2.d0)
            sai(2*j,i)=(ap-ax)/dsqrt(2.d0)
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
+if crlibm
          p(i)=atan_rn(-sai(2*i-1,2*i)/sai(2*i,2*i))
+ei
+if .not.crlibm
          p(i)=atan(-sai(2*i-1,2*i)/sai(2*i,2*i))
+ei
+if crlibm
          s1(2*i-1,2*i-1)=cos_rn(p(i))
+ei
+if .not.crlibm
          s1(2*i-1,2*i-1)=cos(p(i))
+ei
+if crlibm
          s1(2*i,2*i)=cos_rn(p(i))
+ei
+if .not.crlibm
          s1(2*i,2*i)=cos(p(i))
+ei
+if crlibm
          s1(2*i-1,2*i)=sin_rn(p(i))
+ei
+if .not.crlibm
          s1(2*i-1,2*i)=sin(p(i))
+ei
+if crlibm
          s1(2*i,2*i-1)=-sin_rn(p(i))
+ei
+if .not.crlibm
          s1(2*i,2*i-1)=-sin(p(i))
+ei
        enddo
        call mulnd2(s1,sai)
! adjust sa to have sa(1,1)>0 and sa(3,3)>0 rotate by pi if necessary.
        do i=1,nd-ndc
          xd=1.d0
          if(sai(2*i-1,2*i-1).lt.0.d0) xd=-1.d0
          s1(2*i-1,2*i-1)=xd
          s1(2*i-1,2*i)=0.d0
          s1(2*i,2*i-1)=0.d0
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
      end
+dk mulnd2
      subroutine mulnd2(rt,r)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ia,j,ndim,ndim2
      double precision r,rt,rtt
      parameter (ndim2=6)
      parameter (ndim=3)
+ca ii
      dimension rt(ndim2,ndim2),r(ndim2,ndim2),rtt(ndim2,ndim2)
      do 11 i=1,nd2
        do 12 j=1,nd2
 12     rtt(i,j)=0.d0
 11   continue
      do 1 i=1,nd2
        do 2 j=1,nd2
          do 3 ia=1,nd2
            rtt(i,ia)=rt(i,j)*r(j,ia)+rtt(i,ia)
 3        continue
 2      continue
 1    continue

      do 444 i=1,nd2
        do 555 j=1,nd2
          r(i,j)=rtt(i,j)
 555    continue
 444  continue
      return
      end
+dk movearou
      subroutine movearou(rt)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,j,ndim,ndim2
      double precision rt,rto,s,xr,xrold,xy,xyz,xz,xzy,yz
      parameter (ndim2=6)
      parameter (ndim=3)
+ca ii
+ca printing
      dimension rt(ndim2,ndim2),rto(ndim2,ndim2)
      dimension xy(ndim2,ndim2),xz(ndim2,ndim2),yz(ndim2,ndim2)
      dimension xyz(ndim2,ndim2),xzy(ndim2,ndim2)
      dimension s(ndim2,ndim2)
      do 11 i=1,nd2
        do 12 j=1,nd2
          s(i,j)=0.d0
          s(i,i)=1.d0
          xy(i,j)=0.d0
          xz(i,j)=0.d0
          yz(i,j)=0.d0
          xyz(i,j)=0.d0
          xzy(i,j)=0.d0
 12     continue
 11   continue

      xy(1,3)=1.d0
      xy(3,1)=1.d0
      xy(2,4)=1.d0
      xy(4,2)=1.d0
      xy(5,5)=1.d0
      xy(6,6)=1.d0

      xz(1,5)=1.d0
      xz(5,1)=1.d0
      xz(2,6)=1.d0
      xz(6,2)=1.d0
      xz(3,3)=1.d0
      xz(4,4)=1.d0

      yz(3,5)=1.d0
      yz(5,3)=1.d0
      yz(4,6)=1.d0
      yz(6,4)=1.d0
      yz(1,1)=1.d0
      yz(2,2)=1.d0

      xyz(1,3)=1.d0
      xyz(3,5)=1.d0
      xyz(5,1)=1.d0
      xyz(2,4)=1.d0
      xyz(4,6)=1.d0
      xyz(6,2)=1.d0

      xzy(1,5)=1.d0
      xzy(5,3)=1.d0
      xzy(3,1)=1.d0
      xzy(2,6)=1.d0
      xzy(6,4)=1.d0
      xzy(4,2)=1.d0

      ic=0
      xrold=1000000000.d0
      call movemul(rt,s,rto,xr)
! write(6,*) xr,xrold
!  do i=1,6
!       write(6,'(6(1x,1pe12.5))') (RTO(i,j),j=1,6)
!  enddo
!  PAUSE
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
+if cr
        if(idpr.gt.-101) write(lout,*) " no exchanged"
+ei
+if .not.cr
        if(idpr.gt.-101) write(*,*) " no exchanged"
+ei
      elseif(ic.eq.1) then
        call movemul(rt,xy,rto,xr)
+if cr
        if(idpr.gt.-101) write(lout,*) " x-y exchanged"
+ei
+if .not.cr
        if(idpr.gt.-101) write(*,*) " x-y exchanged"
+ei
      elseif(ic.eq.2) then
        call movemul(rt,xz,rto,xr)
+if cr
        if(idpr.gt.-101) write(lout,*) " x-z exchanged"
+ei
+if .not.cr
        if(idpr.gt.-101) write(*,*) " x-z exchanged"
+ei
      elseif(ic.eq.3) then
        call movemul(rt,yz,rto,xr)
+if cr
        if(idpr.gt.-101) write(lout,*) " y-z exchanged"
+ei
+if .not.cr
        if(idpr.gt.-101) write(*,*) " y-z exchanged"
+ei
      elseif(ic.eq.4) then
        call movemul(rt,xyz,rto,xr)
+if cr
        if(idpr.gt.-101) write(lout,*) " x-y-z permuted"
+ei
+if .not.cr
        if(idpr.gt.-101) write(*,*) " x-y-z permuted"
+ei
      elseif(ic.eq.5) then
        call movemul(rt,xzy,rto,xr)
+if cr
        if(idpr.gt.-101) write(lout,*) " x-z-y permuted"
+ei
+if .not.cr
        if(idpr.gt.-101) write(*,*) " x-z-y permuted"
+ei
      endif

      do i=1,nd2
        do j=1,nd2
          rt(i,j)=rto(i,j)
        enddo
      enddo

      return
      end
+dk movemul
      subroutine movemul(rt,xy,rto,xr)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,j,k,ndim,ndim2
      double precision rt,rto,xr,xy
      parameter (ndim2=6)
      parameter (ndim=3)
+ca ii
      dimension rt(ndim2,ndim2)
      dimension xy(ndim2,ndim2),rto(ndim2,ndim2)

      do 11 i=1,nd2
        do 12 j=1,nd2
 12       rto(i,j)=0.d0
 11   continue

      do  i=1,nd2
        do  j=1,nd2
          do  k=1,nd2
            rto(i,k)=xy(i,j)*rt(j,k)+rto(i,k)
          enddo
        enddo
      enddo

      xr=0.d0
      do i=1,nd2
        do j=1,nd2
          xr=xr+abs(rto(i,j))
        enddo
      enddo
      do i=1,nd
        xr=xr-abs(rto(2*i-1,2*i-1))
        xr=xr-abs(rto(2*i-1,2*i))
        xr=xr-abs(rto(2*i,2*i))
        xr=xr-abs(rto(2*i,2*i-1))
      enddo
      return
      end
+dk initpert
      subroutine initpert(st,ang,ra)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,nn,nreso
      double precision ang,ra,st
!   X-RATED
!- SETS UP ALL THE COMMON BLOCKS RELEVANT TO NORMAL FORM AND THE BASIS
!- CHANGES INSIDE  MAPNORMF
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (nreso=20)
      dimension st(ndim),ang(ndim),ra(ndim)
+ca stable
+ca istable
+ca ii
+ca coast
+ca reson
+ca resfile

      if(iref.gt.0) then
+if cr
      write(lout,*) iref
+ei
+if .not.cr
      write(*,*) iref
+ei
      read(iref,*) nres
      if(nres.ge.nreso) then
+if cr
       write(lout,*) ' NRESO IN LIELIB TOO SMALL '
+ei
+if .not.cr
       write(*,*) ' NRESO IN LIELIB TOO SMALL '
+ei
+if cr
      call abend('999                                               ')
+ei
+if .not.cr
       stop 999
+ei
      endif
      elseif(iref.eq.0) then
      nres=0
      endif
+if cr
      if(nres.ne.0) write(lout,*)' warning resonances left in the map'
+ei
+if .not.cr
      if(nres.ne.0) write(*,*)' warning resonances left in the map'
+ei
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
        angle(i)=0.d0
        rad(i)=0.d0
        sta(i)=0.d0
        dsta(i)=1.d0-sta(i)
        ista(i)=0
        idsta(i)=0
      enddo
      do i=1,nd        !  frank          -ndc
        angle(i)=ang(i)
        rad(i)=ra(i)
        sta(i)=st(i)
        dsta(i)=1.d0-sta(i)
      enddo
      do i=1,nd
        ista(i)=idint(sta(i)+.01)
        idsta(i)=idint(dsta(i)+.01)
      enddo
      return
      end
+dk dlie
      double precision function dlie(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      INTEGER J(NTT)
      integer j(*)
+ca ii
      dlie=0.d0
      do 1 i=1,nd
 1    dlie=dble(j(2*i-1)+j(2*i))+dlie
      dlie=dlie+1.d0
      dlie=1.d0/dlie
      return
      end
+dk rext
      double precision function rext(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,lie,mo,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
+ca ii
+ca coast
+ca istable
      integer j(*)
      lie=0
      do i=1,nd-ndc
        lie=ista(i)*j(2*i)+lie
      enddo
      mo=mod(lie,4)+1
      goto(11,12,13,14),mo
 11   rext = 1.d0
      return
 12   rext = -1.d0
      return
 13   rext = -1.d0
      return
 14   rext = 1.d0
      return
      end
+dk cpart
      subroutine cpart(h,ch)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim,ntt
      double precision rext
      parameter (ndim=3)
      parameter (ntt=40)
      external rext
+ca ii
      integer h,ch
      call dacfu(h,rext,ch)
      return
      end
+dk ctoi
      subroutine ctoi(f1,f2)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim2,ntt
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer f1,f2
      integer b1,x(ndim2)
!
!
      call etallnom(b1,1,'B1        ')
      call etallnom(x,nd2  ,'X         ')

      call cpart(f1,b1)
      call etctr(x)
      call trx(b1,f2,x)
      call dadal(x,nd2)
      call dadal(b1,1)
      return
      end
+dk itoc
      subroutine itoc(f1,f2)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ndim2,ntt
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer f1,f2
      integer b1,x(ndim2)
!
      call etallnom(b1,1,'B1        ')
      call etallnom(x,nd2  ,'X         ')

      call etrtc(x)
      call trx(f1,b1,x)
      call cpart(b1,f2)
      call dadal(x,nd2)
      call dadal(b1,1)
      return
      end
+dk etrtc
      subroutine etrtc(x)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
+ca coast
      integer x(*)

      integer rel(ndim2)
!
!
      call etallnom(rel,nd2  ,'REL       ')

      call etini(rel)
      call etini(x)
      do i=1,nd-ndc
        call daadd(rel(2*i-1),rel(2*i),x(2*i-1))
        call dasub(rel(2*i-1),rel(2*i),x(2*i))
      enddo
      call dadal(rel,nd2)
      return
      end
+dk etctr
      subroutine etctr(x)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
+ca coast
      integer x(*)
      integer rel(ndim2)
!
!
      call etallnom(rel,nd2  ,'REL       ')

      call etini(rel)
      call etini(x)
      do i=1,nd-ndc
        call dalin(rel(2*i-1),.5d0,rel(2*i),.5d0,x(2*i-1))
        call dalin(rel(2*i-1),.5d0,rel(2*i),-.5d0,x(2*i))
      enddo
      call dadal(rel,nd2)
      return
      end
+dk etcjg
      subroutine etcjg(x)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,ntt
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca istable
+ca ii
+ca coast
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
      end
+dk eig6
      subroutine eig6(fm,reval,aieval,revec,aievec)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer jet,ndim2
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
      parameter (ndim2=6)
      integer nn
      integer ilo,ihi,mdim,info
+ca ii
+ca coast
      double precision reval(ndim2),aieval(ndim2),                      &
     &revec(ndim2,ndim2),aievec(ndim2,ndim2)
      double precision fm(ndim2,ndim2),aa(ndim2,ndim2)
      integer i,i1
      double precision ort(ndim2),vv(ndim2,ndim2)
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
+if cr
        write(lout,*) '  ERROR IN EIG6'
+ei
+if .not.cr
        write(*,*) '  ERROR IN EIG6'
+ei
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
        if(abs(reval(i)**2+aieval(i)**2 -1.d0).gt.1.d-10) then
+if cr
          write(lout,*) ' EIG6: Eigenvalues off the unit circle!'
+ei
+if .not.cr
          write(*,*) ' EIG6: Eigenvalues off the unit circle!'
+ei
        endif
      enddo
      return
      end
+dk ety
      subroutine ety(nm,n,low,igh,a,ort)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,j,m,n,ii,jj,la,mp,nm,igh,kp1,low
      double precision a(nm,n),ort(igh)
      double precision f,g,h,scale
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
!
      do 180 m = kp1, la
         h = 0.0
         ort(m) = 0.0
         scale = 0.0
!     ********** scale column (algol tol then not needed) **********
         do 90 i = m, igh
   90    scale = scale + abs(a(i,m-1))
!
         if (scale .eq. 0.0) goto 180
         mp = m + igh
!     ********** for i=igh step -1 until m do -- **********
         do 100 ii = m, igh
            i = mp - ii
            ort(i) = a(i,m-1) / scale
            h = h + ort(i) * ort(i)
  100    continue
!
         g = -dsign(dsqrt(h),ort(m))
         h = h - ort(m) * g
         ort(m) = ort(m) - g
!     ********** form (i-(u*ut)/h) * a **********
         do 130 j = m, n
            f = 0.0
!     ********** for i=igh step -1 until m do -- **********
            do 110 ii = m, igh
               i = mp - ii
               f = f + ort(i) * a(i,j)
  110       continue
!
            f = f / h
!
            do 120 i = m, igh
  120       a(i,j) = a(i,j) - f * ort(i)
!
  130    continue
!     ********** form (i-(u*ut)/h)*a*(i-(u*ut)/h) **********
         do 160 i = 1, igh
            f = 0.0
!     ********** for j=igh step -1 until m do -- **********
            do 140 jj = m, igh
               j = mp - jj
               f = f + ort(j) * a(i,j)
  140       continue
!
            f = f / h
!
            do 150 j = m, igh
  150       a(i,j) = a(i,j) - f * ort(j)
!
  160    continue
!
         ort(m) = scale * ort(m)
         a(m,m-1) = scale * g
  180 continue
!
  200 return
!     ********** last card of ety **********
      end
+dk etyt
      subroutine etyt(nm,n,low,igh,a,ort,z)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,j,n,kl,mm,mp,nm,igh,low,mp1
      double precision a(nm,igh),ort(igh),z(nm,n)
      double precision g
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
!
!
!     ********** initialize z to identity matrix **********
      do 80 i = 1, n
!
         do 60 j = 1, n
   60    z(i,j) = 0.0
!
         z(i,i) = 1.0
   80 continue
!
      kl = igh - low - 1
      if (kl .lt. 1) goto 200
!     ********** for mp=igh-1 step -1 until low+1 do -- **********
      do 140 mm = 1, kl
         mp = igh - mm
         if (a(mp,mp-1) .eq. 0.0) goto 140
         mp1 = mp + 1
!
         do 100 i = mp1, igh
  100    ort(i) = a(i,mp-1)
!
         do 130 j = mp, igh
            g = 0.0
!
            do 110 i = mp, igh
  110       g = g + ort(i) * z(i,j)
!     ********** divisor below is negative of h formed in orthes.
!                double division avoids possible underflow **********
            g = (g / ort(mp)) / a(mp,mp-1)
!
            do 120 i = mp, igh
  120       z(i,j) = z(i,j) + g * ort(i)
!
  130    continue
!
  140 continue
!
  200 return
!     ********** last card of etyt **********
      end
+dk ety2
      subroutine ety2(nm,n,low,igh,h,wr,wi,z,ierr)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,                      &
     &igh,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n),z(nm,n)
      double precision p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,machep
      logical notlas
      double precision z3r,z3i
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
!     arithmetic is double precision. complex division
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
      machep = 1.d-17
!     machep = r1mach(4)
!
      ierr = 0
      norm = 0.0
      k = 1
!     ********** store roots isolated by balanc
!                and compute matrix norm **********
      do 50 i = 1, n
!
         do 40 j = k, n
   40    norm = norm + abs(h(i,j))
!
         k = i
         if (i .ge. low .and. i .le. igh) goto 50
         wr(i) = h(i,i)
         wi(i) = 0.0
   50 continue
!
      en = igh
      t = 0.0
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
         if (s .eq. 0.0) s = norm
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
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
!
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75 * s
      y = x
      w = -0.4375 * s * s
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
         h(i,i-2) = 0.0
         if (i .eq. mp2) goto 160
         h(i,i-3) = 0.0
  160 continue
!     ********** double qr step involving rows l to en and
!                columns m to en **********
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) goto 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.0
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.0) goto 260
         p = p / x
         q = q / x
         r = r / x
  170    s = dsign(dsqrt(p*p+q*q+r*r),p)
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
      wi(en) = 0.0
      en = na
      goto 60
!     ********** two roots found **********
  280 p = (y - x) / 2.0
      q = p * p + w
      zz = dsqrt(abs(q))
      h(en,en) = x + t
      x = h(en,en)
      h(na,na) = y + t
      if (q .lt. 0.0) goto 320
!     ********** real pair **********
      zz = p + dsign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.0) wr(en) = x - w / zz
      wi(na) = 0.0
      wi(en) = 0.0
      x = h(en,na)
      s = abs(x) + abs(zz)
      p = x / s
      q = zz / s
      r = dsqrt(p*p+q*q)
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
  340 if (norm .eq. 0.0) goto 1001
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
         h(en,en) = 1.0
         if (na .eq. 0) goto 800
!     ********** for i=en-1 step -1 until 1 do -- **********
         do 700 ii = 1, na
            i = en - ii
            w = h(i,i) - p
            r = h(i,en)
            if (m .gt. na) goto 620
!
            do 610 j = m, na
  610       r = r + h(i,j) * h(j,en)
!
  620       if (wi(i) .ge. 0.0) goto 630
            zz = w
            s = r
            goto 700
  630       m = i
            if (wi(i) .ne. 0.0) goto 640
            t = w
            if (w .eq. 0.0) t = machep * norm
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
  720    call etdiv(z3r,z3i,0.d0,-h(na,en),h(na,na)-p,q)
         h(na,na) = z3r
         h(na,en) = z3i
  730    h(en,na) = 0.0
         h(en,en) = 1.0
         enm2 = na - 1
         if (enm2 .eq. 0) goto 800
!     ********** for i=en-2 step -1 until 1 do -- **********
         do 790 ii = 1, enm2
            i = na - ii
            w = h(i,i) - p
            ra = 0.0
            sa = h(i,en)
!
            do 760 j = m, na
               ra = ra + h(i,j) * h(j,na)
               sa = sa + h(i,j) * h(j,en)
  760       continue
!
            if (wi(i) .ge. 0.0) goto 770
            zz = w
            r = ra
            s = sa
            goto 790
  770       m = i
            if (wi(i) .ne. 0.0) goto 780
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
            vi = (wr(i) - p) * 2.0 * q
            if (vr .eq. 0.0 .and. vi .eq. 0.0) vr = machep * norm       &
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
         do 820 j = i, n
  820    z(i,j) = h(i,j)
!
  840 continue
!     ********** multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- **********
      do 880 jj = low, n
         j = n + low - jj
         m = min0(j,igh)
!
         do 880 i = low, igh
            zz = 0.0
!
            do 860 k = low, m
  860       zz = zz + z(i,k) * h(k,j)
!
            z(i,j) = zz
  880 continue
!
      goto 1001
!     ********** set error -- no convergence to an
!                eigenvalue after 200 iterations **********
 1000 ierr = en
 1001 return
!     ********** last card of ety2 **********
      end
+dk etdiv
      subroutine etdiv(a,b,c,d,e,f)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
!   computes the complex division
!     a + ib = (c + id)/(e + if)
!  very slow, but tries to be as accurate as
!  possible by changing the order of the
!  operations, so to avoid under(over)flow
!  problems.
!  Written by F. Neri Feb. 12 1986
!
      double precision a,b,c,d,e,f
      double precision s,t
      double precision cc,dd,ee,ff
      double precision temp
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
      s = 1.d0/ee
      t = 1.d0/(ee+ ff*(ff*s))
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
      end
+dk sympl3
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
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer n
      parameter ( n = 3 )
      integer kp,kq,lp,lq,jp,jq,i
      double precision m(2*n,2*n)
      double precision qq,pq,qp,pp
!
      do 100 kp=2,2*n,2
        kq = kp-1
        do 200 lp=2,kp-2,2
          lq = lp-1
          qq = 0.d0
          pq = 0.d0
          qp = 0.d0
          pp = 0.d0
          do 300 jp=2,2*n,2
            jq = jp-1
            qq = qq + m(lq,jq)*m(kq,jp) - m(lq,jp)*m(kq,jq)
            pq = pq + m(lp,jq)*m(kq,jp) - m(lp,jp)*m(kq,jq)
            qp = qp + m(lq,jq)*m(kp,jp) - m(lq,jp)*m(kp,jq)
            pp = pp + m(lp,jq)*m(kp,jp) - m(lp,jp)*m(kp,jq)
  300     continue
!         write(6,*) qq,pq,qp,pp
          do 400 i=1,2*n
            m(kq,i) = m(kq,i) - qq*m(lp,i) + pq*m(lq,i)
            m(kp,i) = m(kp,i) - qp*m(lp,i) + pp*m(lq,i)
  400     continue
  200   continue
        qp = 0.d0
        do 500 jp=2,2*n,2
          jq = jp-1
          qp = qp + m(kq,jq)*m(kp,jp) - m(kq,jp)*m(kp,jq)
  500   continue
!       write(6,*) qp
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
      end
+dk averaged
      subroutine averaged(f,a,flag,fave)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer isi,ndim,ndim2,nord,ntt
      double precision avepol
!      TAKES THE AVERAGE OF A FUNCTION F
!  FLAG TRUE A=ONE TURN MAP
!       FALSE A=A_SCRIPT
!
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca printing
+ca ii
      integer f,fave,a(*)
      integer cosi,sine
      logical flag
      external avepol

      integer a1(ndim2),a2(ndim2),xy(ndim2),hf(ndim2),ftf(ndim2)

      if(.not.flag) then
      call etall(cosi,1)
      call etall(sine,1)
      call trx(f,f,a)
      call ctor(f,cosi,sine)
      call dacfu(cosi,avepol,fave)
      call dadal(cosi,1)
      call dadal(sine,1)
      else

      call etall(cosi,1)
      call etall(sine,1)
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
      call facflod(ftf,xy,a1,2,nord,1.d0,-1)
      call trx(f,f,a1)
      call ctor(f,cosi,sine)
      call dacfu(cosi,avepol,fave)

      call dadal(cosi,1)
      call dadal(sine,1)
      call dadal(ftf,nd2)
      call dadal(hf,nd2)
      call dadal(a2,nd2)
      call dadal(a1,nd2)
      call dadal(xy,nd2)

      endif

      return
      end
+dk avepol
      double precision function avepol(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      INTEGER J(NTT)
      integer j(*)
+ca ii
+ca coast
      avepol=1.d0
      do i=1,(nd-ndc)
        if(j(2*i).ne.j(2*i-1)) then
          avepol=0.d0
          return
        endif
      enddo

      return
      end
+dk couplean
      subroutine couplean(map1,tune,map2,oneturn)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim,ndim2,no1,nord,ntt
      double precision crazy,tpi
!  map1 ascript a1 not there
!  tune 2 or 3 tunes

!   map2 ascript with a couple parameter in nv
!  oneturn map created with tunes and map2

      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca coast
+ca ii
      integer map1(*),oneturn(*),map2(*),ftf,hf
      integer xy(ndim2),m1(ndim2),m2(ndim2),a2(ndim2),a1(ndim2)
      integer cs,h

      double precision killnonl,planar,psq(ndim),radsq(ndim)
      double precision tune(ndim)
      external killnonl,planar

      call etall(ftf,1)
      call etall(hf,1)
      call etall(a1,nd2)
      call etall(a2,nd2)
      call etall(m1,nd2)
      call etall(m2,nd2)
      call etall(xy,nd2)
      call etall(cs,1)
      call etall(h,1)

!     map1 is an a-script, the last nv entry should be empty
!  this a-script should around the fixed point to all orders
!     one order is lost because I use PB-field

+if crlibm
      tpi=atan_rn(1.d0)*8.d0
+ei
+if .not.crlibm
      tpi=atan(1.d0)*8.d0
+ei
      do i=1,nd2
        call dacfu(map1(i),killnonl,m1(i))
      enddo

      call etini(xy)
      call daclr(cs)

      do i=1,nd-ndc
        call dasqr(xy(2*i-1),a2(2*i-1))
        call dasqr(xy(2*i),a2(2*i))
        call daadd(a2(2*i-1),a2(2*i),ftf)
        crazy=-tune(i)*tpi/2.d0
        call dacmu(ftf,crazy,ftf)
        call daadd(ftf,cs,cs)
      enddo

      call etinv(m1,m2)
      call trx(cs,h,m2)

      call dacfu(h,planar,cs)
      call dasub(h,cs,h)
      call davar(a2(1),0.d0,nv)

      call damul(a2(1),h,h)
      call daadd(cs,h,h)
      call expnd2(h,xy,xy,1.d-9,1000)

      call dacopd(xy,oneturn)

      nord=1
      call mapnorm(xy,ftf,a2,a1,m2,hf,nord)

      call gettura(psq,radsq)
+if cr
      write(lout,*) (psq(i),i=1,nd)
+ei
+if .not.cr
      write(*,*) (psq(i),i=1,nd)
+ei

      call etini(xy)
      no1=no
      call fexpo(ftf,xy,xy,3,no1,1.d0,-1)
      call etcct(a2,xy,map2)
      call etcct(a1,map2,map2)

      call dadal(ftf,1)
      call dadal(hf,1)
      call dadal(a1,nd2)
      call dadal(a2,nd2)
      call dadal(m1,nd2)
      call dadal(m2,nd2)
      call dadal(xy,nd2)
      call dadal(cs,1)
      call dadal(h,1)

      return
      end
+dk planar
      double precision function planar(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      INTEGER J(NTT)
      integer j(*)
+ca ii
+ca coast
      planar=0.d0
      do i=1,(nd-ndc)
        if(j(2*i).eq.j(2*i-1)) then
          planar=1.d0
          return
        endif
        if(j(2*i).eq.2) then
          planar=1.d0
          return
        endif
        if(j(2*i-1).eq.2) then
          planar=1.d0
          return
        endif
      enddo

      return
      end
+dk killnonl
      double precision function killnonl(j)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ic,ndim
      parameter (ndim=3)
!      PARAMETER (NTT=40)
!      INTEGER J(NTT)
      integer j(*)
+ca ii
+ca coast

      killnonl=1.d0

      ic=0
      do i=1,nd2-ndc2
        ic=ic+j(i)
      enddo
      if(ic.gt.1) killnonl=0.d0
      if(j(nv).ne.0) killnonl=0.d0

      return
      end
+dk fexpo1
      subroutine fexpo1(h,x,w,nrmin,nrmax,sca,ifac)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer ifac,ndim,ndim2,nrma,nrmax,nrmi,nrmin,ntt
      double precision sca
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
+ca ii
      integer x,w,h

      integer v(ndim2)

      nrmi=nrmin-1
      nrma=nrmax-1
      call etall(v,nd2)
      call difd(h,v,-1.d0)
      call facflo(v,x,w,nrmi,nrma,sca,ifac)
      call dadal(v,nd2)

      return
      end
+dk etcctpar
      subroutine etcctpar(x,ix,xj,z)
      implicit none
+if cr
+ca crcoall
+ei
+if crlibm
+ca crlibco
+ei
      integer i,ie,ix,ndim,ndim2,ntt
      double precision xj
      parameter (ndim=3)
      parameter (ndim2=6)
      parameter (ntt=40)
      dimension xj(*)
      dimension ie(ntt)
+ca ii
      integer  x(*),z(*)

      call etallnom(ie,nv,'IE        ')
      do i=1,nd2
      call davar(ie(i),0.d0,i)
      enddo
      do  i=nd2+1,nv
      call dacon(ie(i),xj(i-nd2))
      enddo

      call dacct(x,ix,ie,nv,z,ix)

      call dadal(ie,nv)
      return
      end
