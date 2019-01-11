subroutine daliesix

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_commond
  use mod_time
  use mod_lie_dab, only : mld_allocArrays

  implicit none

  integer i,mf1,mf2,mf3,mf4,mf5,mfile,nd2,ndim,ndpt,nis,no,nv,damap,a1,a1i,a2,a2i,f,fc,fs,rot,xy,h,hc,hs,h4,df,bb1,bb2,haux
  real tlim,time0,time1,time
  real(kind=fPrec) angle,coe,rad,x2pi
  dimension damap(6),a1(6),a1i(6),a2(6),a2i(6)
  dimension rot(6),xy(6),df(6)
  dimension angle(3),rad(3)

  save

  tlim=1e7
  call time_timerStart
  time0=0.
  call time_timerCheck(time0)

  ! Initialization
  x2pi=atan_mb(one)*eight
  coe=(-one*two)/x2pi
  no=nordf
  if(nord1.gt.no) nord1=no
  ndim=nvar2/2
  if(nvarf/2.lt.ndim) ndim=nvarf/2
  if(ndim == 0) then
    write(lout,"(a)") "DALIESIX> ERROR Number of normal form variables have to be: 2, 4, 5, 6 + parameters."
    call prror(-1)
  end if
  nv=nvarf
  nd2=2*ndim
  ndpt=idptr
  nis=0
  mfile=18
  mf1=21
  mf2=22
  mf3=23
  mf4=24
  mf5=25
  call idprset(-102)
  call mld_allocArrays(.true.)
  call lieinit(no,nv,ndim,ndpt,0,nis)
  call etall(damap,nd2)
  call etall(xy,nd2)
  call etall(bb1,1)
  call etall(bb2,1)
  call etall(haux,1)
  call etallnom(a1,nd2,'A1        ')
  call etallnom(a1i,nd2,'A1I       ')
  call etallnom(a2,nd2,'A2        ')
  call etallnom(a2i,nd2,'A2I       ')
  call etallnom(f,1,'F         ')
  call etallnom(fc,1,'FC        ')
  call etallnom(fs,1,'FS        ')
  call etallnom(rot,nd2,'ROT       ')
  call etallnom(h,1,'H         ')
  call etallnom(h4,1,'H4        ')
  call etallnom(hc,1,'HC        ')
  call etallnom(hs,1,'HS        ')
  call etallnom(df,nd2,'DF        ')
  rewind mfile
  rewind 111
  rewind mf1
  rewind mf2
  rewind mf3
  rewind mf4
  rewind mf5
  call daread(damap,nd2,mfile,zero)

  ! Normal Form Analysis
  call mapnorm(damap,f,a2,a1,xy,h,nord1)
  call dainv(a1,nv,a1i,nv)
  call dainv(a2,nv,a2i,nv)
  call ctor(f,fc,fs)
  call gettura(angle,rad)
  call taked(xy,1,rot)
  call take(h,2,haux)
  call dasub(h,haux,h4)
  call ctor(h,hc,hs)
  call dhdj(h,df)

  ! Printing
  call daprid(a1i,1,nd2,mf1)
  call daprid(a2i,1,nd2,mf1)
  call dapri(f,mf1)
  call daprid(rot,1,nd2,mf1)
  call dapri(h4,mf1)
  call dapri(f,mf1)
  call daprid(a2,1,nd2,mf1)
  call daprid(a1,1,nd2,mf1)
  call dapri(h,mf2)
  call dapri(hc,mf2)
  call dapri(hs,mf2)
  call dapri(fc,mf2)
  call dapri(fs,mf2)
  call daprimax(hc,mf3)
  call daprid(df,ndim+1,nd2,mf4)
  call daprid(df,1,ndim,mf5)
  write(lout,10060)
  if(imod1.eq.0) then
    write(lout,10020) nordf
  else
    write(lout,10010) nordf
  endif
  write(lout,10025) nord1
  if(imod2.eq.0) then
    write(lout,10040) nvarf
  else
    write(lout,10030) nvarf
  endif
  write(lout,10050)
  angle(3)=angle(3)*(-one)
  write(lout,*) (angle(i),i=1,ndim)

  ! Clean-Up
  call dadal(damap,nd2)
  call dadal(a1,nd2)
  call dadal(a1i,nd2)
  call dadal(a2,nd2)
  call dadal(a2i,nd2)
  call dadal(f,1)
  call dadal(fc,1)
  call dadal(fs,1)
  call dadal(rot,nd2)
  call dadal(xy,nd2)
  call dadal(h,1)
  call dadal(h4,1)
  call dadal(hc,1)
  call dadal(hs,1)
  call dadal(df,nd2)
  call dadal(bb1,1)
  call dadal(bb2,1)
  call dadal(haux,1)
  time1=0.
  call time_timerCheck(time1)
  time = time1-time0
  write(lout,10000) no,time

  return

10000 format(/10x,'DA-CALCULATION OF ORDER : ',i7,' TAKES ', f12.3,' SECONDS OF CPU TIME'//131('-')//)
10010 format(t10,'THE ORDER IS GREATER THAN THE ONE SPECIFIED IN THE'/t10,'DIFFERENTIAL ALGEBRA BLOCK.'// t10,'NEW ORDER---> ',i3)
10020 format(t10,'ORDER FOR THE NORMAL FORM CALCULATIONS---> ',i3/)
10025 format(/t10,'CLOSED ORBIT ORDER OF THE NORMAL FORM ---> ',i3/)
10030 format(t10,'THE NUMBER OF VARIABLES IS GREATER THAN THE ONE SPECIFIED IN THE'/ &
             t10,'DIFFERENTIAL ALGEBRA BLOCK.'// t10,'NEW NUMBER OF VARIABLES---> ',i3/)
10040 format(t10,'NUMBER OF VARIABLES FOR THE NORMAL FORM CALCULATIONS---> ',i3//)
10050 format(t10,'LINEAR TUNES USED IN THE NORMAL FORM CALCULATIONS:'//)
10060 format(//131('-')//t10,20('O')/t10,2('O'),16x,2('O')/t10,'OO  NORMAL FORMS  OO', /&
             t10,2('O'),16x,2('O')/t10,20('O')//130('-')//)

end subroutine daliesix

!-----------------------------------------------------------------------
!  CALCULATION OF THE 4-DIMENSIONAL CLOSED ORBIT INCLUDING DELTA
!-----------------------------------------------------------------------
subroutine mydaini(ncase,nnord,nnvar,nndim,nnvar2,nnord1)

  use floatPrecision
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_commond
  use mod_common,  only : iqmodc,ichromc,ilinc
  use mod_lie_dab, only : iscrda,mld_allocArrays

  implicit none

  integer idummy,ncase,ndimfo,ndpt,nis,nndim,nnord,nnord1,nnvar,nnvar2,nord1o,nordo,nvar2o,nvaro
  real(kind=fPrec) am
  dimension am(6,6),idummy(6)
  save

  if(nndim < 2 .or. nndim > 3) then
    write(lout,"(a)") "DAINI> ERROR DA corrections implemented for 4D and 6D only."
    call prror(-1)
  end if

  nordo=nord
  nvaro=nvar
  ndimfo=ndimf
  nvar2o=nvar2
  nord1o=nord1

  nord=nnord
  nvar=nnvar
  ndimf=nndim
  nvar2=nnvar2
  nord1=nnord1

  ndpt=0
  nis=0

  call daeps(preda)
  call idprset(-102)
  call mld_allocArrays(.true.)
  call lieinit(nord,nvar,ndimf,ndpt,0,nis)
  write(lout,10000) nord,nvar,nndim
  call daall(iscrda,100,'$$IS      ',nord,nvar)

  ! closed orbit
  if(ncase.eq.1) call clorda(2*ndimf,idummy,am)

  ! tune variation
  if(ncase.eq.2) call umlauda
  rewind 18
  rewind 111

  ! main map calculation
  if(ncase.eq.3) call runda

  ! %*6 map calculation
  if(ncase.eq.4) call runcav
  iqmodc=0
  ichromc=0
  ilinc=0
  call dadal(iscrda,100)

  nord=nordo
  nvar=nvaro
  nvar2=nvar2o
  ndimf=ndimfo
  nord1=nord1o

  return
10000 format(/131('-')/10x,'DA INITIALIZATION: ORDER = ',i2,', # of VARIABLES = ',i2,', DIMENSION = ',i2/)
end subroutine mydaini

!-----------------------------------------------------------------------
!               DIFFERENTIAL ALGEBRA FOR CAVITY
!                          AUGUST 1994
!-----------------------------------------------------------------------
subroutine runcav

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use crcoall
  use parpro
  use mod_time
  use mod_common
  use mod_commonmn, only : e0f
  use mod_commons
  use mod_commont, only : comt_daStart,comt_daEnd
  use mod_commond
  use mod_hions
  use mod_lie_dab, only : idao,rscrri,iscrda

  implicit none

  integer idaa
  real(kind=fPrec) betr0,dare,sigmdac
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  save
!-----------------------------------------------------------------------
#include "include/daini.f90"
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
  call comt_daStart
  betr0=sqrt(one-(nucm0/e0)**2)
  write(lout,*) ' REENTERING MAP '
  call davar(x(1),zero,1)
  call davar(yp(1),zero,2)
  call davar(x(2),zero,3)
  call davar(yp(2),zero,4)
  call davar(deltas,zero,5)
  call davar(dpda1,zero,6)
  call darea(x(1),18)
  call darea(yp(1),18)
  call darea(x(2),18)
  call darea(yp(2),18)
  call darea(deltas,18)
  call darea(dpda1,18)
  rewind 18
!Eric
    rewind 111
  if(ition.ne.0) then
  e0f=sqrt(e0**2-nucm0**2)                                             !hr08
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  EJF1=E0F*(ONE+DPDA)/(NUCM0/NUCMDA) ;
!FOX  EJ1=SQRT(EJF1*EJF1+NUCMDA*NUCMDA) ;
!FOX  EJF0=EJF1 ;
    if(abs(dppoff).gt.pieni) then
      sigmdac=sigmoff(iicav)
!FOX  DELTAS=DELTAS-SIGMDAC ;

    endif
    call synoda
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=YP(1)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=YP(2)*MTCDA/(ONE+DPDA) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;

  endif
  if(nvar2.ge.1) call dapri(x(1),18)
  if(nvar2.ge.2) call dapri(yp(1),18)
  if(nvar2.ge.3) call dapri(x(2),18)
  if(nvar2.ge.4) call dapri(yp(2),18)
  if(nvar2.eq.5) call dapri(dpda1,18)
  if(nvar2.eq.6) call dapri(deltas,18)
  if(nvar2.eq.6) call dapri(dpda1,18)

  write(lout,*) 'END COORDINATES'
  write(lout,*) dare(x(1)),dare(y(1))
  write(lout,*) dare(x(2)),dare(y(2))
  write(lout,*) dare(sigmda),dare(dpda)

  write(12,'(E22.15)') dare(x(1))
  write(12,'(E22.15)') dare(y(1))
  write(12,'(E22.15)') dare(x(2))
  write(12,'(E22.15)') dare(y(2))
  write(12,'(E22.15)') dare(sigmda)
  write(12,'(E22.15)') dare(dpda)

  write(lout,10010)
!-----------------------------------------------------------------------
!     DADAL AUTOMATIC INCLUSION
  time2=0.
  call time_timerCheck(time2)
!     time=time2-time1
  write(lout,10020) time1-time0
  write(lout,10030) nord,time2-time1
  call comt_daEnd
  return
10000 format(/t10,'TRACKING ENDED ABNORMALLY'&
             /t10,'PARTICLE NO. ',i7,' LOST IN REVOLUTION ',i8,' AT ELEMENT ',i4&
             /t10,'HORIZ:  AMPLITUDE = ',ES23.16,'   APERTURE = ',f15.3&
             /t10,'VERT:   AMPLITUDE = ',ES23.16,'   APERTURE = ',f15.3&
             /t10,'ELEMENT - LIST NUMBER ',i4,' TYP NUMBER ',i4,' NAME ',a16/)
10010 format(//t10,30('*')/t10,'**** ONE TURN COMPLETED ****'/ t10,30('*')/)
10020 format(/10x,'The Preparating Calculations took',f12.3,' second(s) of Computing Time')
10030 format(/10x,'DA-Calculation of Order : ',i7,' took ', f12.3,' second(s) of CPU Time'//131('-')//)
end subroutine runcav

!-----------------------------------------------------------------------
!           UMSCHR    DIFFERENTIAL ALGEBRA 5 -> 6
!                          AUGUST 1994
!-----------------------------------------------------------------------
subroutine umschr(iu1,iu2)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall

  implicit none

  integer i,ii,iio,io,ioo,iplus,iu1,iu2,j,jj,nno
  real(kind=fPrec) c,c1
  character(len=80) aaa
  character(len=18) a18
  character(len=58) a58
  dimension jj(100)
  save

  do i=1,100
    jj(i)=0
  enddo
  do j=1,6
    iplus=0
    read(iu1,'(A80)') aaa
    write(iu2,'(A80)') aaa
    read(iu1,'(A18,I4,A58)') a18,nno,a58
    write(iu2,'(A18,I4,A58)') a18,nno,a58
    read(iu1,'(A80)') aaa
    write(iu2,'(A80)') aaa
    read(iu1,'(A80)') aaa
    write(iu2,'(A80)') aaa
    read(iu1,'(A80)') aaa
    write(iu2,'(A80)') aaa

    if(nno.eq.1) then
      do i=1,5
        read(iu1,'(6X,2X,G21.14,I5)') c,ii
        write(iu2,'(6X,2X,G21.14,I5)') c,ii
      end do
      if (j.eq.5) then
        write(iu2,'(6X,2X,G21.14,I5)') one,5
      else
        write(iu2,'(6X,2X,G21.14,I5)') zero,5
      endif
      read(iu1,'(6X,2X,G21.14,I5)') c,ii
      write(iu2,'(6X,2X,G21.14,I5)') c,ii+1
      read(iu1,'(6X,2X,G21.14,I5)') c,ii
      if(ii.ne.0) write(lout,*) ' ERROR IN UMSCHR'
    else
20     read(iu1,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') ii,c,io,(jj(i),i=1,5)
      if(ii.eq.0) then
        if(j.eq.5.and.ioo.lt.2) then
          write(iu2,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') iio+1,one,1, 0,0,0,0,1,0
          write(iu2,*) one
        endif
        goto 30
      endif
      ioo=io
      iio=ii
      read(iu1,*) c1
      if(j.eq.5.and.(io.eq.2.or.jj(5).eq.1).and.iplus.eq.0) then
        iplus=1
        write(iu2,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') ii,one,io,0,0,0,0,1,0
        write(iu2,*) one
      endif
      write(iu2,'(I6,2X,G21.14,I5,4X,18(2I2,1X))') ii+iplus,c,io,(jj(i),i=1,4),0,jj(5)
      write(iu2,*) c1
      goto 20
    endif
30   write(iu2,*)
  end do
  rewind iu1
  rewind iu2
  return

end subroutine umschr

! ================================================================================================ !
!  CENTRAL LOOP FOR NORMAL FORWARD-TRACKING
!  SPECIALLY PREPARED FOR NEW D.A.
!  5 --> 6  AND  ASD6 / ALD6
! ================================================================================================ !
subroutine runda

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use physical_constants
  use crcoall
  use parpro
  use mod_common
  use mod_commonmn, only : e0f,numx
  use mod_commons
  use mod_commont, only : xxtr,yytr,comt_daStart,comt_daEnd
  use mod_commond
  use mod_commond2
  use mod_hions
  use mod_lie_dab, only : idao,iscrri,rscrri,iscrda
  use mod_units
  use mod_time
  use mod_fluc,    only : fluc_errAlign,fluc_writeFort4

  implicit none

  integer i,ich,i11,i480,icav,ien,ifam,iflag,iflag1,iflag2,ii,ip,ipch,irrtr,iverg,ix,j,jb,jj,jmel,  &
    jx,k,kk,kkk,kpz,kzz,n,ncyo,nmz,nsta,nsto,idaa
  real(kind=fPrec) beamoff1,beamoff2, beamoff3, beamoff4,beamoff5,beamoff6,benkcc,betr0,cbxb,       &
    cbzb,cik,crk,crxb,crzb,dare,dpdav,dpdav2,dummy,fake,ox,oxp,oz,ozp,r0,r000,r0a,r2b,r2bf,rb,rbf,  &
    rho2b,rkb,rkbf,scikveb,scrkveb,sigmdac,startco,tkb,xbb,xrb,xs,zbb,zfeld1,zfeld2,zrb,zs,crabfreq,&
    crabpht,crabpht2,crabpht3,crabpht4
  logical fErr
  character(len=300) ch
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  dimension zfeld1(100),zfeld2(100)
  dimension iverg(mcor)
  dimension fake(2,20),dpdav2(6),jj(100)
  save
!-----------------------------------------------------------------------
#include "include/daini.f90"
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
  call comt_daStart
  if(mout2.eq.1) then
    call f_open(unit=99,file="fort.99",formatted=.true.,mode="w",err=fErr,recl=303)
  end if
  do i=1,100
    jj(i)=0
  end do

  if(mout2.eq.1) write(7,*) e0,nucm0
  e0f=sqrt(e0**2-nucm0**2)
  betr0=sqrt(one-(nucm0/e0)**2)
  do i=1,mcor
    iverg(i)=0
  end do
  do i=1,20
    fake(1,i)=zero
    fake(2,i)=zero
  end do
  time1=0.
  call time_timerCheck(time1)
  if(niu(1).gt.1) then
    do i=1,2
      ii=2*i
      xxtr(1,i)=clon(ii-1)
      yytr(1,i)=clon(ii)
    enddo
    sigm(1)=clon(5)
    dps(1)=clon(6)
  endif
  ox=xxtr(1,1)
  oxp=yytr(1,1)
  oz=xxtr(1,2)
  ozp=yytr(1,2)
  if(nvar2.ge.1) call davar(x(1),zero,1)
!FOX  X(1)=X(1)+OX ;
  if(nvar2.ge.2) call davar(yp(1),zero,2)
!FOX  YP(1)=YP(1)+OXP*(ONE+DPS(1)) ;
  if(nvar2.ge.3) call davar(x(2),zero,3)
!FOX  X(2)=X(2)+OZ ;
  if(nvar2.ge.4) call davar(yp(2),zero,4)
!FOX  YP(2)=YP(2)+OZP*(ONE+DPS(1)) ;
  if(nvar2.lt.5) then
!FOX  DPDA1=DPS(1)*C1E3 ;
  endif
  if(nvar2.eq.5) then
    call davar(dpda1,zero,5)
!FOX  DPDA1=DPDA1+DPS(1)*C1E3 ;
  endif
  if(nvar2.eq.6) then
    call davar(deltas,zero,5)
    call davar(dpda1,zero,6)
!FOX  DPDA1=DPDA1+DPS(1)*C1E3 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  EJF1=E0F*(ONE+DPDA)/(NUCM0/NUCMDA) ;
!FOX  EJ1=SQRT(EJF1*EJF1+NUCM0*NUCM0) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DELTAS=DELTAS+SIGM(1)/RV ;
  else
!FOX  DELTAS=SIGM(1) ;
  endif
!FOX  CORROLD(1)=X(1) ;
!FOX  CORROLD(2)=YP(1) ;
!FOX  CORROLD(3)=X(2) ;
!FOX  CORROLD(4)=YP(2) ;
!FOX  CORROLD(5)=DELTAS ;
!FOX  CORROLD(6)=DPDA1 ;
  do kkk=1,6
    dpdav=dare(corrold(kkk))
!FOX  CORROLD(KKK)=CORROLD(KKK)-DPDAV ;
  end do
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=YP(1)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=YP(2)*MTCDA/(ONE+DPDA) ;
!FOX  SIGMDA=DELTAS*RV ;
  write(lout,*) ' ENTERING MAP '
  write(lout,*) 'INITIAL COORDINATES'
  write(lout,*) dare(x(1)),dare(y(1))
  write(lout,*) dare(x(2)),dare(y(2))
  write(lout,*) dare(sigmda),dare(dpda)
  if(ncor.gt.0) then
    do i=1,ncor
      do ii=1,iu
        if(ipar(i).eq.(ic(ii)-nblo).and.iverg(i).eq.0) then
          iverg(i)=1
          call davar(smida(i),ed(ic(ii)-nblo),nvar2+i)
        endif
      end do
    end do
  endif
  nsta=niu(1)
  nsto=niu(2)
  if(niu(2).lt.niu(1)) nsto=nsto+iu
  do 490 n=1,numl
    numx=n-1
!FOX  EJF1=E0F*(ONE+DPDA)/(NUCM0/NUCMDA) ;
!FOX  EJ1=SQRT(EJF1*EJF1+NUCMDA*NUCMDA) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
    ncyo=ncy
    if(ncy.eq.0) ncy=1
    if(ithick.eq.1) call envada
    ncy=ncyo
    iflag=0
    iflag1=0
    iflag2=0
    icav=0
    do 480 i480=nsta,nsto
      if(i480.gt.iu) then
        i=i480-iu
      else
        i=i480
      endif
      if(mout2.eq.1.and.i480.eq.nsta.and.n.eq.1) call fluc_writeFort4
      if(iflag.eq.1) then
!FOX  EJF1=E0F*(ONE+DPDA)/(NUCM0/NUCMDA) ;
!FOX  EJ1=SQRT(EJF1*EJF1+NUCMDA*NUCMDA) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DELTAS=SIGMDA/RV ;
        if(ithick.eq.1) then
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
          if(icav.eq.0) then
!FOX  CORRNEW(1)=X(1) ;
!FOX  CORRNEW(2)=YP(1) ;
!FOX  CORRNEW(3)=X(2) ;
!FOX  CORRNEW(4)=YP(2) ;
!FOX  CORRNEW(5)=DELTAS ;
!FOX  CORRNEW(6)=DPDA1 ;
            do 24 kkk=1,6
              dpdav=dare(corrnew(kkk))
!FOX  CORRNEW(KKK)=CORRNEW(KKK)-DPDAV ;
24           continue
          else
!FOX  CORRAU2(1)=X(1) ;
!FOX  CORRAU2(2)=YP(1) ;
!FOX  CORRAU2(3)=X(2) ;
!FOX  CORRAU2(4)=YP(2) ;
!FOX  CORRAU2(5)=DELTAS ;
!FOX  CORRAU2(6)=DPDA1 ;
            do 25 kkk=1,6
!FOX  CORRAU1(KKK)=CORRNEW(KKK) ;
              dpdav=dare(corrau2(kkk))
!FOX  CORRAU2(KKK)=CORRAU2(KKK)-DPDAV ;
25           continue
            if(ncor.gt.0) then
              do kkk=1,ncor
                kk=6+kkk
!FOX  CORRAU2(KK)=SMIDA(KKK) ;
                dpdav=dare(smida(kkk))
!FOX  CORRAU1(KK)=SMIDA(KKK)-DPDAV ;
              enddo
            endif
            call dacct(corrau2,nvar,corrau1,nvar,corrnew,nvar)
          endif
          dpdav=dare(x(1))
!FOX  X(1)=CORROLD(1)+DPDAV ;
          dpdav=dare(yp(1))
!FOX  YP(1)=CORROLD(2)+DPDAV ;
          dpdav=dare(x(2))
!FOX  X(2)=CORROLD(3)+DPDAV ;
          dpdav=dare(yp(2))
!FOX  YP(2)=CORROLD(4)+DPDAV ;
          dpdav=dare(deltas)
!FOX  DELTAS=CORROLD(5)+DPDAV ;
          dpdav=dare(dpda1)
!FOX  DPDA1=CORROLD(6)+DPDAV ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=YP(1)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=YP(2)*MTCDA/(ONE+DPDA) ;
!FOX  EJF1=E0F*(ONE+DPDA)/(NUCM0/NUCMDA) ;
!FOX  EJ1=SQRT(EJF1*EJF1+NUCMDA*NUCMDA) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  SIGMDA=DELTAS*RV ;
          icav=icav+1
          call envada
        endif
        iflag=0
      endif
      ix=ic(i)
      if(mout2.eq.1) then
        if(i480.eq.nsta.and.n.eq.1) then
          write(ch,*) 'START ',dare(x(1)),dare(y(1)),dare(x(2)),dare(y(2)),dare(sigmda),dare(dpda)
          do ich=300,1,-1
            if(ch(ich:ich).ne.' ') goto 700
          enddo
700          write(99,'(a)') ch(:ich)
        else
          if(ic(i-1).le.nblo) then
            write(ch,*) bez(mtyp(ic(i-1),mel(ic(i-1)))),dare(x(1)),dare(y(1)),dare(x(2)),dare(y(2)),dare(sigmda),dare(dpda)
            do ich=300,1,-1
              if(ch(ich:ich).ne.' ') goto 701
            enddo
701            write(99,'(a)') ch(:ich)
          else
            write(ch,*) bez(ic(i-1)-nblo),dare(x(1)),dare(y(1)),dare(x(2)),dare(y(2)),dare(sigmda),dare(dpda)
            do ich=300,1,-1
              if(ch(ich:ich).ne.' ') goto 702
            enddo
702            write(99,'(a)') ch(:ich)
          endif
        endif
      endif
      if(ix.gt.nblo) goto 70
      if(mout2.eq.1.and.n.eq.1) then
        jmel=mel(ix)
        do jb=1,jmel
          jx=mtyp(ix,jb)
          if(el(jx).eq.zero) then
            write(7,*) '0'
            write(7,*) bez(jx)
          else
            if(ithick.eq.1) then
              ifam=0
              do ip=1,il
                if(kz(ip).eq.kz(jx).and.el(ip).ne.zero) then
                  if((ed(jx).ne.zero.and.ed(ip).ne.zero).or.(ek(jx).ne.zero.and.ek(ip).ne.zero)) then
                    ifam=ifam+1
                    if(bez(ip).eq.bez(jx)) goto 35
                  endif
                endif
              enddo
35               continue
              write(7,*) '1'
              write(7,*) bez(jx)
              if(kz(jx).eq.1.and.abs(ed(jx)).le.pieni) then
                write(7,*) '0 ',ifam
                write(7,*) zero,ek(jx),el(jx)
              elseif(kz(jx).eq.3.and.abs(ed(jx)).le.pieni) then
                write(7,*) '0 ',ifam
                write(7,*) zero,ek(jx),el(jx)
              elseif(kz(jx).eq.4.and.abs(ed(jx)).le.pieni) then
                write(7,*) '0 ',ifam
                write(7,*) zero,ek(jx),el(jx)
              elseif(kz(jx).eq.5.and.abs(ed(jx)).le.pieni) then
                write(7,*) '0 ',ifam
                write(7,*) zero,ek(jx),el(jx)
              elseif(kz(jx).eq.8.and.abs(ed(jx)).le.pieni) then
                write(7,*) '0 ',ifam
                write(7,*) zero,ek(jx),el(jx)
              elseif(kz(jx).eq.2.and.abs(ek(jx)).le.pieni) then
                write(7,*) '0 ',ifam
                write(7,*) ed(jx),zero,el(jx)
              elseif((kz(jx).eq.6.or.kz(jx).eq.7).and.abs(ed(jx)).le.pieni.and.abs(ek(jx)).le.pieni) then
                write(7,*) '0 ',ifam
                write(7,*) zero,zero,el(jx)
              elseif((kz(jx).eq.6.or.kz(jx).eq.7).and.(abs(ed(jx)).le.pieni.and.abs(ek(jx)).gt.pieni)) then
                write(7,*) '2 ',ifam
                write(7,*) zero,ek(jx),el(jx)
              elseif((kz(jx).eq.6.or.kz(jx).eq.7).and.(abs(ed(jx)).gt.pieni.and.abs(ek(jx)).le.pieni)) then
                write(7,*) '3 ',ifam
                write(7,*) ed(jx),zero,el(jx)
              else
                write(7,*) kz(jx),ifam
                write(7,*) ed(jx),ek(jx),el(jx)
              endif
            else
              write(7,*) '2'
              write(7,*) bez(jx)
              write(7,*) el(jx)
            endif
          endif
        enddo
      endif
      if(ix <= 0) then
        write(lout,"(a)") "RUNDA> ERROR Inverted linear blocks not allowed."
        call prror(-1)
      endif
#include "include/dalin1.f90"
#include "include/dalin2.f90"
#include "include/dalin3.f90"
#include "include/dalin4.f90"
#include "include/dalin5.f90"
#include "include/dalin6.f90"
!FOX  SIGMDA=SIGMDA+
#include "include/sqrtfox.f90"
          endif
        enddo
      endif
      goto 480
70     ix=ix-nblo
      if(abs(dare(x(1))).gt.aint(aper(1)).or.abs(dare(x(2))).gt.aint(aper(2))) then
        write(lout,10000) j,numx,i,dare(x(1)),aper(1),dare(x(2)),aper(2),ix, kz(ix),bez(ix)
        goto 520
      endif
      kpz=abs(kp(ix))
      if (kpz.ge.0.and.kpz.le.5) goto 110
      if (kpz.eq.6) goto 90
      goto 480
90     if(nvar2.le.4.or.(nvar2.eq.5.and.nsix.ne.2)) goto 480
      ixcav=ix
      iicav=i
      if(nsix.eq.2) then
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
!FOX  DELTAS=SIGMDA*RV ;
        call dapri(x(1),19)
        call dapri(yp(1),19)
        call dapri(x(2),19)
        call dapri(yp(2),19)
        call dapri(deltas,19)
        call dapri(dpda1,19)
        if(ncor.gt.0) then
          write(lout,*) ' WARNING: in the 5*6 mode no extra parameters allowed'
        endif
        rewind 19
!     DADAL AUTOMATIC INCLUSION
        return
      endif
      if(ition.ne.0) then
!FOX  EJF0=EJF1 ;
        if(abs(dppoff).gt.pieni) then
          sigmdac=sigmoff(i)
!FOX  SIGMDA=SIGMDA-SIGMDAC ;
        endif
        call synoda
        if(mout2.eq.1.and.n.eq.1) then
          write(7,*) '5'
          if(kz(ix).eq.12) then
            write(7,*) bez(ix)
            write(ch,*) ed(ix),hsyc(ix),itionc(ix),phasc(ix)
            do ich=300,1,-1
              if(ch(ich:ich).ne.' ') goto 703
            enddo
703            write(7,'(a)') ch(:ich)
          else
            write(7,*) 'CAV'
            write(ch,*) hsy(1),hsy(3),ition,phas
            do ich=300,1,-1
              if(ch(ich:ich).ne.' ') goto 704
            enddo
704            write(7,'(a)') ch(:ich)
          endif
        endif
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;

      endif
      if(nvar2.eq.6.and.nsix.ne.2) then
        iflag=1
        iflag1=1
        iflag2=1
      endif
      goto 480
110     kzz=kz(ix)
      if(mout2.eq.1.and.n.eq.1) then
        if(kzz.eq.0) then
          write(7,*) '0'
          write(7,*) bez(ix)
        elseif((kzz.ge.1.and.kzz.le.10).or.(kzz.le.-1.and.kzz.ge.-10)) then
          write(7,*) '3'
          write(7,*) bez(ix)
#ifdef TILT
          write(7,*) xsi(i),zsi(i),atan2_mb(tilts(i),tiltc(i))
#endif
#ifndef TILT
          write(7,*) xsi(i),zsi(i),zero
#endif
          write(7,*) kzz,smi(i)
        elseif(kzz.eq.11) then
          nmz=nmu(ix)
          write(7,*) '4'
          write(7,*) bez(ix)
#ifdef TILT
          write(7,*) xsi(i),zsi(i),atan2_mb(tilts(i),tiltc(i))
#endif
#ifndef TILT
          write(7,*) xsi(i),zsi(i),zero
#endif
          if(abs(dki(ix,1)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
              write(7,*) nmz,' 1',' 1',' 1'
              write(7,*) dki(ix,1),dki(ix,3)
            else
              write(7,*) nmz,' 1',' 1',' 0'
              write(7,*) dki(ix,1),dki(ix,3)
            endif
          elseif(abs(dki(ix,2)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
              write(7,*) nmz,' 1',' 0',' 1'
              write(7,*) dki(ix,2),dki(ix,3)
            else
              write(7,*) nmz,' 1',' 0',' 0'
              write(7,*) dki(ix,2),dki(ix,3)
            endif
          else
            write(7,*) nmz,' 0',' 0',' 0'
          endif
          if(nmz.ge.1) then
            do ip=1,nmz
              write(7,*) bbiv(ip,i),aaiv(ip,i)
            enddo
          endif
        endif
      endif
      if(kzz.eq.15) then
!FOX  XX(1)=X(1) ;
!FOX  XX(2)=X(2) ;
!FOX  YY(1)=Y(1) ;
!FOX  YY(2)=Y(2) ;
      call wireda(ix,i)
!FOX  Y(1)=YY(1) ;
!FOX  Y(2)=YY(2) ;
          goto 480
      endif

      if(kzz.eq.20.and.parbe(ix,2).eq.zero) then
        if(nbeam.ge.1) then
          if(sigman(1,imbb(i)).eq.sigman(2,imbb(i))) then
            if(ibeco.eq.1) then
              if(ibbc.eq.0) then
                crk=parbe(ix,5)
                cik=parbe(ix,6)
              else
                crk=parbe(ix,5)*bbcu(imbb(i),11) + parbe(ix,6)*bbcu(imbb(i),12)
                cik=parbe(ix,6)*bbcu(imbb(i),11) - parbe(ix,5)*bbcu(imbb(i),12)
              endif
              rho2b=crk**2+cik**2
              if(rho2b.gt.pieni) then
                if(abs(sigman(1,imbb(i))).lt.pieni) goto 9088
                tkb=rho2b/((two*sigman(1,imbb(i)))*sigman(1,imbb(i)))
                beamoff4=(((crad*ptnfac(ix))*crk)/rho2b)*(one-exp_mb(-one*tkb))
                beamoff5=(((crad*ptnfac(ix))*cik)/rho2b)*(one-exp_mb(-one*tkb))
              endif
            endif
#include "include/beamcof.f90"
!FOX  RHO2BF=CRKVEBF*CRKVEBF+CIKVEBF*CIKVEBF ;
            if(abs(dare(rho2bf)).gt.pieni) then
              if(abs(sigman(1,imbb(i))).lt.pieni) goto 9088
!FOX  TKBF=RHO2BF/(TWO*SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I))) ;
              if(ibbc.eq.0) then
!FOX   Y(1)=Y(1)+(CRAD*CRKVEBF/RHO2BF*
!FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)*MTCDA/(ONE+DPDA) ;
!FOX   Y(2)=Y(2)+(CRAD*CIKVEBF/RHO2BF*
!FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)*MTCDA/(ONE+DPDA) ;
              else
!FOX   CCCC=(CRAD*CRKVEBF/RHO2BF*
!FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)*BBCU(IMBB(I),11)-
!FOX   (CRAD*CIKVEBF/RHO2BF*
!FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)*BBCU(IMBB(I),12) ;
!FOX   Y(1)=Y(1)+CCCC*MTCDA/(ONE+DPDA) ;
!FOX   CCCC=(CRAD*CRKVEBF/RHO2BF*
!FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF4)*BBCU(IMBB(I),12)+
!FOX   (CRAD*CIKVEBF/RHO2BF*
!FOX   PTNFAC(IX)*(ONE-EXP(-TKBF))-BEAMOFF5)*BBCU(IMBB(I),11) ;
!FOX   Y(2)=Y(2)+CCCC*MTCDA/(ONE+DPDA) ;
              endif
            endif
          else if(sigman(1,imbb(i)).gt.sigman(2,imbb(i))) then
            if(ibeco.eq.1) then
              if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
              r2b=two*(sigman(1,imbb(i))**2-sigman(2,imbb(i))**2) !hr08
              rb=sqrt(r2b)
              rkb=((crad*ptnfac(ix))*pisqrt)/rb
              if(ibbc.eq.0) then
                crk=parbe(ix,5)
                cik=parbe(ix,6)
              else
                crk=parbe(ix,5)*bbcu(imbb(i),11) + parbe(ix,6)*bbcu(imbb(i),12)
                cik=parbe(ix,6)*bbcu(imbb(i),11) - parbe(ix,5)*bbcu(imbb(i),12)
              endif
              xrb=abs(crk)/rb
              zrb=abs(cik)/rb
              call errf(xrb,zrb,crxb,crzb)
              if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
              tkb=(crk**2/sigman(1,imbb(i))**2+cik**2/sigman(2,imbb(i))**2)*half
              xbb=(sigman(2,imbb(i))/sigman(1,imbb(i)))*xrb
              zbb=(sigman(1,imbb(i))/sigman(2,imbb(i)))*zrb
              call errf(xbb,zbb,cbxb,cbzb)
              beamoff4=(rkb*(crzb-exp_mb(-one*tkb)*cbzb))*sign(one,crk)
              beamoff5=(rkb*(crxb-exp_mb(-one*tkb)*cbxb))*sign(one,cik)
            endif
            if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
            r2bf=two*(sigman(1,imbb(i))**2-sigman(2,imbb(i))**2) !hr08
            rbf=sqrt(r2bf)
            rkbf=((crad*ptnfac(ix))*pisqrt)/rbf
#include "include/beamcof.f90"
!FOX  XRBF=CRKVEBF/RBF ;
            if(dare(xrbf).lt.zero) then
!FOX  XRBF=-XRBF ;
            endif
!FOX  ZRBF=CIKVEBF/RBF ;
            if(dare(zrbf).lt.zero) then
!FOX  ZRBF=-ZRBF ;
            endif
            call errff(xrbf,zrbf,crxbf,crzbf)
            if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
!FOX  TKBF=(CRKVEBF*CRKVEBF/(SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I)))+
!FOX  CIKVEBF*CIKVEBF/(SIGMAN(2,IMBB(I))*SIGMAN(2,IMBB(I))))*HALF ;
!FOX  XBBF=SIGMAN(2,IMBB(I))/SIGMAN(1,IMBB(I))*XRBF ;
!FOX  ZBBF=SIGMAN(1,IMBB(I))/SIGMAN(2,IMBB(I))*ZRBF ;
            call errff(xbbf,zbbf,cbxbf,cbzbf)
            scrkveb=sign(one,dare(crkvebf))
            scikveb=sign(one,dare(cikvebf))
            if(ibbc.eq.0) then
!FOX  Y(1)=Y(1)+(RKBF*(CRZBF-EXP(-TKBF)*
!FOX  CBZBF)*SCRKVEB-BEAMOFF4)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=Y(2)+(RKBF*(CRXBF-EXP(-TKBF)*
!FOX  CBXBF)*SCIKVEB-BEAMOFF5)*MTCDA/(ONE+DPDA) ;
            else
!FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*
!FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),11)-
!FOX  (RKBF*(CRXBF-EXP(-TKBF)*
!FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),12) ;
!FOX   Y(1)=Y(1)+CCCC*MTCDA/(ONE+DPDA) ;
!FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*
!FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),12)+
!FOX  (RKBF*(CRXBF-EXP(-TKBF)*
!FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),11) ;
!FOX   Y(2)=Y(2)+CCCC*MTCDA/(ONE+DPDA) ;
            endif
          else if(sigman(1,imbb(i)).lt.sigman(2,imbb(i))) then
            if(ibeco.eq.1) then
              if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
              r2b=two*(sigman(2,imbb(i))**2-sigman(1,imbb(i))**2)
              rb=sqrt(r2b)
              rkb=((crad*ptnfac(ix))*pisqrt)/rb
              if(ibbc.eq.0) then
                crk=parbe(ix,5)
                cik=parbe(ix,6)
              else
                crk=parbe(ix,5)*bbcu(imbb(i),11) + parbe(ix,6)*bbcu(imbb(i),12)
                cik=parbe(ix,6)*bbcu(imbb(i),11) - parbe(ix,5)*bbcu(imbb(i),12)
              endif
              xrb=abs(crk)/rb
              zrb=abs(cik)/rb
              call errf(zrb,xrb,crzb,crxb)
              if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
              tkb=(crk**2/sigman(1,imbb(i))**2+cik**2/sigman(2,imbb(i))**2)*half
              xbb=(sigman(2,imbb(i))/sigman(1,imbb(i)))*xrb
              zbb=(sigman(1,imbb(i))/sigman(2,imbb(i)))*zrb
              call errf(zbb,xbb,cbzb,cbxb)
              beamoff4=(rkb*(crzb-exp_mb(-one*tkb)*cbzb))*sign(one,crk)
              beamoff5=(rkb*(crxb-exp_mb(-one*tkb)*cbxb))*sign(one,cik)
            endif
            if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
            r2bf=two*(sigman(2,imbb(i))**2-sigman(1,imbb(i))**2)
            rbf=sqrt(r2bf)
            rkbf=((crad*ptnfac(ix))*pisqrt)/rbf
#include "include/beamcof.f90"
!FOX  XRBF=CRKVEBF/RBF ;
            if(dare(xrbf).lt.zero) then
!FOX  XRBF=-XRBF ;
            endif
!FOX  ZRBF=CIKVEBF/RBF ;
            if(dare(zrbf).lt.zero) then
!FOX  ZRBF=-ZRBF ;
            endif
            call errff(zrbf,xrbf,crzbf,crxbf)
            if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
!FOX  TKBF=(CRKVEBF*CRKVEBF/(SIGMAN(1,IMBB(I))*SIGMAN(1,IMBB(I)))+
!FOX  CIKVEBF*CIKVEBF/(SIGMAN(2,IMBB(I))*SIGMAN(2,IMBB(I))))*HALF ;
!FOX  XBBF=SIGMAN(2,IMBB(I))/SIGMAN(1,IMBB(I))*XRBF ;
!FOX  ZBBF=SIGMAN(1,IMBB(I))/SIGMAN(2,IMBB(I))*ZRBF ;
            call errff(zbbf,xbbf,cbzbf,cbxbf)
            scrkveb=sign(one,dare(crkvebf))
            scikveb=sign(one,dare(cikvebf))
            if(ibbc.eq.0) then
!FOX  Y(1)=Y(1)+(RKBF*(CRZBF-EXP(-TKBF)*
!FOX  CBZBF)*SCRKVEB-BEAMOFF4)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=Y(2)+(RKBF*(CRXBF-EXP(-TKBF)*
!FOX  CBXBF)*SCIKVEB-BEAMOFF5)*MTCDA/(ONE+DPDA) ;
            else
!FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*
!FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),11)-
!FOX  (RKBF*(CRXBF-EXP(-TKBF)*
!FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),12) ;
!FOX   Y(1)=Y(1)+CCCC*MTCDA/(ONE+DPDA) ;
!FOX  CCCC=(RKBF*(CRZBF-EXP(-TKBF)*
!FOX  CBZBF)*SCRKVEB-BEAMOFF4)*BBCU(IMBB(I),12)+
!FOX  (RKBF*(CRXBF-EXP(-TKBF)*
!FOX  CBXBF)*SCIKVEB-BEAMOFF5)*BBCU(IMBB(I),11) ;
!FOX   Y(2)=Y(2)+CCCC*MTCDA/(ONE+DPDA) ;
            endif
          endif
          goto 480
        endif
        goto 480
      endif
      if(kzz.eq.20.and.parbe(ix,2).gt.zero) then
#include "include/beam6dfi.f90"
        goto 480
      endif
      if(kzz.eq.22) then
        irrtr=imtr(ix)
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
!FOX  PUSIG=((EJ1-E0)/E0F)*C1E3*(E0/E0F) ;
!FOX  TEMPI(1) = X(1) ;
!FOX  TEMPI(2) = YP(1) ;
!FOX  TEMPI(3) = X(2) ;
!FOX  TEMPI(4) = YP(2) ;
!FOX  TEMPI(5) = SIGMDA ;
!FOX  TEMPI(6) = PUSIG ;
!FOX  X(1)=COTR(IRRTR,1) +
!FOX  RRTR(IRRTR,1,1)*TEMPI(1)+RRTR(IRRTR,1,2)*TEMPI(2)+
!FOX  RRTR(IRRTR,1,3)*TEMPI(3)+RRTR(IRRTR,1,4)*TEMPI(4)+
!FOX  RRTR(IRRTR,1,5)*TEMPI(5)+RRTR(IRRTR,1,6)*TEMPI(6) ;
!FOX  YP(1)=COTR(IRRTR,2) +
!FOX  RRTR(IRRTR,2,1)*TEMPI(1)+RRTR(IRRTR,2,2)*TEMPI(2)+
!FOX  RRTR(IRRTR,2,3)*TEMPI(3)+RRTR(IRRTR,2,4)*TEMPI(4)+
!FOX  RRTR(IRRTR,2,5)*TEMPI(5)+RRTR(IRRTR,2,6)*TEMPI(6) ;
!FOX  X(2)=COTR(IRRTR,3) +
!FOX  RRTR(IRRTR,3,1)*TEMPI(1)+RRTR(IRRTR,3,2)*TEMPI(2)+
!FOX  RRTR(IRRTR,3,3)*TEMPI(3)+RRTR(IRRTR,3,4)*TEMPI(4)+
!FOX  RRTR(IRRTR,3,5)*TEMPI(5)+RRTR(IRRTR,3,6)*TEMPI(6) ;
!FOX  YP(2)=COTR(IRRTR,4) +
!FOX  RRTR(IRRTR,4,1)*TEMPI(1)+RRTR(IRRTR,4,2)*TEMPI(2)+
!FOX  RRTR(IRRTR,4,3)*TEMPI(3)+RRTR(IRRTR,4,4)*TEMPI(4)+
!FOX  RRTR(IRRTR,4,5)*TEMPI(5)+RRTR(IRRTR,4,6)*TEMPI(6) ;
!FOX  SIGMDA=COTR(IRRTR,5)+
!FOX  RRTR(IRRTR,5,1)*TEMPI(1)+RRTR(IRRTR,5,2)*TEMPI(2)+
!FOX  RRTR(IRRTR,5,3)*TEMPI(3)+RRTR(IRRTR,5,4)*TEMPI(4)+
!FOX  RRTR(IRRTR,5,5)*TEMPI(5)+RRTR(IRRTR,5,6)*TEMPI(6) ;
!FOX  PUSIG=COTR(IRRTR,6)+
!FOX  RRTR(IRRTR,6,1)*TEMPI(1)+RRTR(IRRTR,6,2)*TEMPI(2)+
!FOX  RRTR(IRRTR,6,3)*TEMPI(3)+RRTR(IRRTR,6,4)*TEMPI(4)+
!FOX  RRTR(IRRTR,6,5)*TEMPI(5)+RRTR(IRRTR,6,6)*TEMPI(6) ;
!FOX  EJ1 = E0F*PUSIG/(C1E3*(E0/E0F))+E0 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=YP(1)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=YP(2)*MTCDA/(ONE+DPDA) ;
        end if
        if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22.or.kzz.eq.15) then
        if(bez(ix).eq.'DAMAP') then
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
!FOX  DELTAS=SIGMDA/RV ;
          if(icav.eq.0.or.ithick.ne.1) then
            if(nvar2.ge.1) call dapri(x(1),17)
            if(nvar2.ge.2) call dapri(yp(1),17)
            if(nvar2.ge.3) call dapri(x(2),17)
            if(nvar2.ge.4) call dapri(yp(2),17)
            if(nvar2.eq.5.and.nsix.ne.2) call dapri(dpda1,17)
            if(nvar2.eq.6.or.nsix.eq.2) call dapri(deltas,17)
            if(nvar2.eq.6.or.nsix.eq.2) call dapri(dpda1,17)
          else
!FOX  CORRAU1(1)=X(1) ;
!FOX  CORRAU1(2)=YP(1) ;
!FOX  CORRAU1(3)=X(2) ;
!FOX  CORRAU1(4)=YP(2) ;
!FOX  CORRAU1(5)=DELTAS ;
!FOX  CORRAU1(6)=DPDA1 ;
            do kkk=1,6
              dpdav2(kkk)=dare(corrau1(kkk))
!FOX  CORRAU1(KKK)=CORRAU1(KKK)-DPDAV2(KKK) ;
            end do
            if(ncor.gt.0) then
              do kkk=1,ncor
                kk=6+kkk
!FOX  CORRAU1(KK)=SMIDA(KKK) ;
                dpdav=dare(smida(kkk))
!FOX  CORRNEW(KK)=SMIDA(KKK)-DPDAV ;
              enddo
            endif
            call dacct(corrau1,nvar,corrnew,nvar,corrau2,nvar)
            do kkk=1,6
!FOX  CORRAU2(KKK)=CORRAU2(KKK)+DPDAV2(KKK) ;
            end do
            call dapri(corrau2(1),17)
            call dapri(corrau2(2),17)
            call dapri(corrau2(3),17)
            call dapri(corrau2(4),17)
            call dapri(corrau2(5),17)
            call dapri(corrau2(6),17)
          endif
          if(ncor.gt.0) then
            do i11=1,ncor
              call dapri(smida(i11),17)
            end do
          endif
        endif
        goto 480
      endif

    if(kzz.eq.23) then
!FOX  CRABAMP=ED(IX)*ZZ0 ;

        crabfreq=ek(ix)*c1e3
        crabpht=crabph(ix)

!FOX  Y(1)=Y(1) - CRABAMP*C1E3/E0F*
!FOX  SIN((SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT))*MTCDA/(ONE+DPDA) ;

!FOX  EJ1=EJ1-CRABAMP
!FOX  *CRABFREQ*2D0*PI/(CLIGHT)*X(1)*
!FOX  COS((SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT)) ;

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;

!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;

      goto 440
  endif
    if(kzz.eq.-23) then
!FOX  CRABAMP=ED(IX)*ZZ0 ;
        crabfreq=ek(ix)*c1e3
        crabpht=crabph(ix)
!FOX  Y(2)=Y(2) - CRABAMP*C1E3/E0F*
!FOX  SIN((SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT))*MTCDA/(ONE+DPDA) ;

!FOX  EJ1=EJ1-CRABAMP
!FOX  *CRABFREQ*2D0*PI/(CLIGHT)*X(2)*
!FOX  COS((SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT)) ;

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;

!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
      goto 440
  endif

! JBG RF CC Multipoles
    if(kzz.eq.26) then
        ! JBG bypass this element if 4D/5D case
        if(iclo6.eq.0) then
            goto 440
        endif
      xs=xsi(i) ! JBG change of variables for misal calculations
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP2=ED(IX)*ZZ0 ;

    crabfreq=ek(ix)*c1e3 !JBG Input in MHz changed to kHz
    crabpht2=crabph2(ix)
!FOX  Y(1)=Y(1) + (CRABAMP2*CRKVE)*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI
!FOX  + CRABPHT2)*MTCDA/(ONE+DPDA);
!FOX  Y(2)=Y(2) - (CRABAMP2*CIKVE)*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI
!FOX  + CRABPHT2)*MTCDA/(ONE+DPDA);
!FOX  EJ1=EJ1 - (0.5D0)*(CRABAMP2)*(CRKVE*CRKVE-
!FOX  CIKVE*CIKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*E0F*C1M3*
!FOX  SIN(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI+CRABPHT2) ;

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;

!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
      goto 440
  endif
      if(kzz.eq.-26) then
        ! JBG bypass this element if 4D/5D case
        if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
            goto 440
        endif
      xs=xsi(i) ! JBG change of variables for misal calculations
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP2=ED(IX)*ZZ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht2=crabph2(ix)
!FOX  Y(2)=Y(2) + (CRABAMP2*CRKVE)*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI + CRABPHT2)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=Y(1) + (CRABAMP2*CIKVE)*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI + CRABPHT2)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  EJ1=EJ1 - (0.5D0)*(CRABAMP2)*(CIKVE*CRKVE)
!FOX  *(((CRABFREQ*2D0)*PI)/CLIGHT)*E0F*C1M3*
!FOX  SIN(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI+CRABPHT2) ;

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
      endif
      if(kzz.eq.27) then
        ! JBG bypass this element if 4D/5D case
        if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
            goto 440
        endif
      xs=xsi(i)
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP3=ED(IX)*ZZ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht3=crabph3(ix)
!FOX  Y(1)=Y(1) + 2D0*(1D0/2D0)*CRABAMP3*((CRKVE*CRKVE)-
!FOX  (CIKVE*CIKVE))*C1M3*MTCDA/(ONE+DPDA)*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI + CRABPHT3);
!FOX  Y(2)=Y(2) - 2D0*CRABAMP3*(CRKVE*CIKVE)*C1M3*MTCDA/(ONE+DPDA)*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI + CRABPHT3);
!FOX  EJ1=EJ1 - 2D0*(1/6D0)*(CRABAMP3)*(CRKVE*CRKVE*CRKVE-
!FOX  3*CIKVE*CIKVE*CRKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*
!FOX  C1M6*E0F*
!FOX  SIN(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI+CRABPHT3);

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
      goto 440
      endif
      if(kzz.eq.-27) then
        ! JBG bypass this element if 4D/5D case
        if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
            goto 440
        endif
      xs=xsi(i)
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP3=ED(IX)*ZZ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht3=crabph3(ix)
!FOX  Y(2)=Y(2) - CRABAMP3*((CIKVE*CIKVE)-
!FOX  (CRKVE*CRKVE))*C1M3*MTCDA/(ONE+DPDA)*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI + CRABPHT3);
!FOX  Y(1)=Y(1) + 2D0*CRABAMP3*(CRKVE*CIKVE)*C1M3*MTCDA/(ONE+DPDA)*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI + CRABPHT3);
!FOX  EJ1=EJ1 + 2D0*(1D0/6D0)*(CRABAMP3)*(CIKVE*CIKVE*CIKVE-
!FOX  3*CIKVE*CRKVE*CRKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*
!FOX  C1M6*E0F*
!FOX  SIN(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI+CRABPHT3);

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
      endif
      if(kzz.eq.28) then
        ! JBG bypass this element if 4D/5D case
        if(iclo6.eq.0) then
!                write(*,*)'Bypassing RF mult 4D or 5D case'
            goto 440
        endif
      xs=xsi(i)
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP4=ED(IX)*ZZ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht4=crabph4(ix)
!FOX  Y(1)=Y(1) + (CRABAMP4)*
!FOX  (CRKVE*CRKVE*CRKVE-(3D0*CRKVE*CIKVE*CIKVE))*C1M6*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI + CRABPHT4)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=Y(2) - (CRABAMP4)*
!FOX  (3D0*CIKVE*CRKVE*CRKVE-CIKVE*CIKVE*CIKVE)*C1M6*
!FOX  COS(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI + CRABPHT4)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  EJ1=EJ1-6D0*(1D0/24D0)*(CRABAMP4)*(CRKVE*CRKVE*CRKVE*CRKVE-
!FOX  6*CRKVE*CRKVE*CIKVE*CIKVE+CIKVE*CIKVE*CIKVE*CIKVE)*
!FOX  C1M9*(((CRABFREQ*2D0)*PI)/CLIGHT)*E0F*
!FOX  SIN(SIGMDA/(CLIGHT*(E0F/E0))*CRABFREQ*2D0*PI+CRABPHT4) ;


!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;

      endif
      ipch=0
      if(ncor.gt.0) then
        do i11=1,ncor
          if(ipar(i11).eq.ix) ipch=i11
        end do
      end if
      if(ipch.ne.0) then
!FOX  EKK=(SMIDA(IPCH)+SMIZF(I))*MTCDA/(ONE+DPDA) ;
      else
!FOX  EKK=SMI(I)*MTCDA/(ONE+DPDA) ;
      end if
      xs=xsi(i)
      zs=zsi(i)
      if(mout2 == 1 .and. n == 1 .and. icextal(i) /= 0) then
        write(27,"(a16,2x,1p,2d14.6,d17.9)") bez(ix),&
          fluc_errAlign(1,icextal(i)),fluc_errAlign(2,icextal(i)),fluc_errAlign(3,icextal(i))
      end if
#include "include/alignf.f90"
      select case (kzz)
      case (1)  ! HORIZONTAL DIPOLE
!FOX  EKK=EKK*C1E3 ;
#include "include/kickf01h.f90"
      case (2)  ! NORMAL QUADRUPOLE
#include "include/kickfxxh.f90"
      case (3)  ! NORMAL SEXTUPOLE
!FOX  EKK=EKK*C1M3 ;
#include "include/kickfho.f90"
#include "include/kickfxxh.f90"
      case (4)  ! NORMAL OCTUPOLE
!FOX  EKK=EKK*C1M6 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxh.f90"
      case (5)  ! NORMAL DECAPOLE
!FOX  EKK=EKK*C1M9 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxh.f90"
      case (6)  ! NORMAL DODECAPOL
!FOX  EKK=EKK*C1M12 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxh.f90"
      case (7)  ! NORMAL 14-POL
!FOX  EKK=EKK*C1M15 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxh.f90"
      case (8)  ! NORMAL 16-POL
!FOX  EKK=EKK*C1M18 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxh.f90"
      case (9)  ! NORMAL 18-POL
!FOX  EKK=EKK*C1M21 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxh.f90"
      case (10) ! NORMAL 20-POL
!FOX  EKK=EKK*C1M24 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxh.f90"
      case (11)
        r0  = ek(ix)
        nmz = nmu(ix)
#include "include/multf01.f90"
        if (abs(r0).le.pieni.or.nmz.eq.0) goto 480
        if(mout2.eq.1.and.n.eq.1) then
          benkcc = ed(ix)*benkc(irm(ix))
          r0a    = one
          r000   = r0*r00(irm(ix))
          do j=1,mmul
            fake(1,j)=(bbiv(j,i)*r0a)/benkcc
            fake(2,j)=(aaiv(j,i)*r0a)/benkcc
            r0a=r0a*r000
          end do

          write(9,'(a16)') bez(ix)
          write(9,'(1p,3d23.15)') (fake(1,j), j=1,3)
          write(9,'(1p,3d23.15)') (fake(1,j), j=4,6)
          write(9,'(1p,3d23.15)') (fake(1,j), j=7,9)
          write(9,'(1p,3d23.15)') (fake(1,j), j=10,12)
          write(9,'(1p,3d23.15)') (fake(1,j), j=13,15)
          write(9,'(1p,3d23.15)') (fake(1,j), j=16,18)
          write(9,'(1p,2d23.15)') (fake(1,j), j=19,20)
          write(9,'(1p,3d23.15)') (fake(2,j), j=1,3)
          write(9,'(1p,3d23.15)') (fake(2,j), j=4,6)
          write(9,'(1p,3d23.15)') (fake(2,j), j=7,9)
          write(9,'(1p,3d23.15)') (fake(2,j), j=10,12)
          write(9,'(1p,3d23.15)') (fake(2,j), j=13,15)
          write(9,'(1p,3d23.15)') (fake(2,j), j=16,18)
          write(9,'(1p,2d23.15)') (fake(2,j), j=19,20)

          do j=1,20
            fake(1,j)=zero
            fake(2,j)=zero
          end do
        end if
        if(nmz.ge.2) then
#include "include/multf02.f90"
          do k=3,nmz
#include "include/multf03.f90"
          end do
#include "include/multf04.f90"
        else
#include "include/multf05.f90"
        end if
      case (12,13,14,15,16,17,18,19,20,21,22,23)
        goto 480
      case (24) ! DIPEDGE ELEMENT
#include "include/kickfdpe.f90"
      case (25) ! Solenoid
#include "include/kickfso1.f90"

      !-----------------
      !--SKEW ELEMENTS--
      !-----------------
      case (-1)  ! VERTICAL DIPOLE
!FOX  EKK=EKK*C1E3 ;
#include "include/kickf01v.f90"
      case (-2)  ! SKEW QUADRUPOLE
#include "include/kickfxxv.f90"
      case (-3)  ! SKEW SEXTUPOLE
!FOX  EKK=EKK*C1M3 ;
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
      case (-4)  ! SKEW OCTUPOLE
!FOX  EKK=EKK*C1M6 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
      case (-5)  ! SKEW DECAPOLE
!FOX  EKK=EKK*C1M9 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
      case (-6)  ! SKEW DODECAPOL
!FOX  EKK=EKK*C1M12 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
      case (-7)  ! SKEW 14-POL
!FOX  EKK=EKK*C1M15 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
      case (-8)  ! SKEW 16-POL
!FOX  EKK=EKK*C1M18 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
      case (-9)  ! SKEW 18-POL
!FOX  EKK=EKK*C1M21 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
      case (-10) ! SKEW 20-POL
!FOX  EKK=EKK*C1M24 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
      end select
! ----------------------------------------------------------------------
! This is a repeat of the case for SKEW 14-POL because it is called by
! a number of gotos in the DA version
440   continue
!FOX  EKK=EKK*C1M15 ;
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfho.f90"
#include "include/kickfxxv.f90"
! ----------------------------------------------------------------------
480   continue
    if(mout2.eq.1) then
      if(ic(iu).le.nblo) then
        write(ch,*) bez(mtyp(ic(iu),mel(ic(iu)))),dare(x(1)),       &
  &dare(y(1)),dare(x(2)),dare(y(2)),dare(sigmda),dare(dpda)
        do ich=300,1,-1
          if(ch(ich:ich).ne.' ') goto 705
        enddo
705        write(99,'(a)') ch(:ich)
      else
        write(ch,*) bez(ic(iu)-nblo),dare(x(1)),                    &
  &dare(y(1)),dare(x(2)),dare(y(2)),dare(sigmda),dare(dpda)
        do ich=300,1,-1
          if(ch(ich:ich).ne.' ') goto 706
        enddo
706        write(99,'(a)') ch(:ich)
      endif
    endif
490 continue
500 continue
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
!FOX  DELTAS=SIGMDA/RV ;
  if(nsix.eq.2) nsix=0
  if(icav.eq.0.or.ithick.ne.1) then
    if(nvar2.ge.1) call dapri(x(1),18)
    if(nvar2.ge.2) call dapri(yp(1),18)
    if(nvar2.ge.3) call dapri(x(2),18)
    if(nvar2.ge.4) call dapri(yp(2),18)
    if(nvar2.eq.5) call dapri(dpda1,18)
    if(nvar2.eq.6) call dapri(deltas,18)
    if(nvar2.eq.6) call dapri(dpda1,18)
  else
!FOX  CORRAU1(1)=X(1) ;
!FOX  CORRAU1(2)=YP(1) ;
!FOX  CORRAU1(3)=X(2) ;
!FOX  CORRAU1(4)=YP(2) ;
!FOX  CORRAU1(5)=DELTAS ;
!FOX  CORRAU1(6)=DPDA1 ;
    do 505 kkk=1,6
      dpdav2(kkk)=dare(corrau1(kkk))
!FOX  CORRAU1(KKK)=CORRAU1(KKK)-DPDAV2(KKK) ;
505   continue
    if(ncor.gt.0) then
      do kkk=1,ncor
        kk=6+kkk
!FOX  CORRAU1(KK)=SMIDA(KKK) ;
        dpdav=dare(smida(kkk))
!FOX  CORRNEW(KK)=SMIDA(KKK)-DPDAV ;
      enddo
    endif
    call dacct(corrau1,nvar,corrnew,nvar,corrau2,nvar)
    do 506 kkk=1,6
!FOX  CORRAU2(KKK)=CORRAU2(KKK)+DPDAV2(KKK) ;
506   continue
!FOX  CORRAU1(2)=CORRAU2(2)/(ONE+CORRAU2(6)) ;
!FOX  CORRAU1(4)=CORRAU2(4)/(ONE+CORRAU2(6)) ;
!FOX  X(1)=CORRAU2(1) ;
!FOX  Y(1)=CORRAU1(2) ;
!FOX  X(2)=CORRAU2(3) ;
!FOX  Y(2)=CORRAU1(4) ;
!FOX  SIGMDA=CORRAU2(5)*RV ;
!FOX  DPDA1=CORRAU2(6) ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
    call dapri(corrau2(1),18)
    call dapri(corrau2(2),18)
    call dapri(corrau2(3),18)
    call dapri(corrau2(4),18)
    call dapri(corrau2(5),18)
    call dapri(corrau2(6),18)
  endif
  if(ncor.gt.0) then
    do i11=1,ncor
      call dapri(smida(i11),18)
    end do
  endif
  write(lout,*) 'END COORDINATES'
  write(lout,*) dare(x(1)),dare(y(1))
  write(lout,*) dare(x(2)),dare(y(2))
  write(lout,*) dare(sigmda),dare(dpda)

  write(12,'(E22.15)') dare(x(1))
  write(12,'(E22.15)') dare(y(1))
  write(12,'(E22.15)') dare(x(2))
  write(12,'(E22.15)') dare(y(2))
  write(12,'(E22.15)') dare(sigmda)
  write(12,'(E22.15)') dare(dpda)

  write(lout,10010)

520 continue
!     DADAL AUTOMATIC INCLUSION
  time2=0.
  call time_timerCheck(time2)
!     time=time2-time1
  write(lout,10020) time1-time0
  write(lout,10030) nord,time2-time1
!-----------------------------------------------------------------------
call comt_daEnd
  return

9088 continue
  write(lout,"(a)") "RUNDA> ERROR Either normalized emittances or the resulting sigma values equal to zero for beam-beam/"
  call prror(-1)
  return

10000 format(/t10,'TRACKING ENDED ABNORMALLY'/t10, 'PARTICLE NO. ',     &
  &i7,' LOST IN REVOLUTION ',i8,' AT ELEMENT ',i4/ t10,              &
  &'HORIZ:  AMPLITUDE = ',ES23.16,'   APERTURE = ',f15.3/ t10,       &
  &'VERT:   AMPLITUDE = ',ES23.16,'   APERTURE = ',f15.3/ t10,       &
  &'ELEMENT - LIST NUMBER ',i4,' TYP NUMBER ',i4,' NAME ',a16/)
10010 format(//t10,30('*')/t10,'**** ONE TURN COMPLETED ****'/ t10,30(  &
  &'*')/)
10020 format(/10x,'The Preparating Calculations took',f12.3,' second(s)'&
  &,' of Computing Time')
10030 format(/10x,'DA-Calculation of Order : ',i7,' took ', f12.3,      &
  &' second(s) of CPU Time'//131('-')//)
end subroutine runda

!-----------------------------------------------------------------------
!  CALCULATION OF INITIAL COORDINATES
!-----------------------------------------------------------------------
subroutine anfb(tas)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  implicit none
  integer i,ii,jj,l,ll
  real(kind=fPrec) bet0s1,bet0x2,bet0z2,chi,co,dchi,dpsic,dsign,si,tas,tas56,x1,x11,x13,x2
  dimension tas(6,6),x1(6),x2(6)
  save
!-----------------------------------------------------------------------
  write(lout,10030)
  if(itra.eq.0) goto 60
  tas56=tas(5,6)*c1m3
  bet0x2=tas(1,3)**2+tas(1,4)**2                                     !hr08
  bet0z2=tas(3,1)**2+tas(3,2)**2                                     !hr08
  bet0s1=tas(5,5)**2+tas56**2                                        !hr08
  dsign=one
  if(tas(3,3).lt.-one*pieni) rat=-one*rat                            !hr08
  if(rat.lt.-one*pieni) dsign=-one*one
  x11=amp(1)/(sqrt(bet0(1))+sqrt(abs(rat)*bet0x2))
  x13=(x11*dsign)*sqrt(abs(rat))                                     !hr08
  amp(2)=(dsign*real(1-iver,fPrec))*(abs(x11)*sqrt(bet0z2)+abs(x13)*sqrt(bet0(2)))
  x1(5)=zero
  if(iclo6.eq.1.or.iclo6.eq.2) then
    x1(6)=(dp1-clop6(3))*sqrt(bet0s1)
  else
    x1(6)=dp1*sqrt(bet0s1)
  endif
  chi=chi0*rad
  dchi=chid*rad
  do 50 i=1,itra
    si=sin_mb(chi)
    co=cos_mb(chi)
    x1(1)=x11*co
    x1(2)=x11*si
    x1(3)=x13*co
    x1(4)=x13*si
    do 20 ii=1,6
      x2(ii)=zero
      do 10 jj=1,6
        x2(ii)=x2(ii)+tas(ii,jj)*x1(jj)
10     continue
20   continue
    if(iclo6.eq.1.or.iclo6.eq.2) then
      x2(2)=x2(2)/((one+x2(6))+clop6(3))                             !hr08
      x2(4)=x2(4)/((one+x2(6))+clop6(3))                             !hr08
    endif
    if(abs(bet0s1).le.pieni) x2(6)=dp1
    if(iver.eq.1) then
      x2(3)=zero
      x2(4)=zero
    endif
    do 30 l=1,2
      ll=(l-1)*2
      x(i,l)=x2(1+ll)+exz(i,1+ll)
      y(i,l)=x2(2+ll)+exz(i,2+ll)
30   continue
    sigm(i)=x2(5)+exz(i,5)
    dps(i)=x2(6)
    dpsic=dps(i)+clop6(3)
    if(idp.eq.1.and.abs(ition).eq.1.and.iclo6.eq.0) then
      do 40 l=1,2
        x(i,l)=x(i,l)+di0(l)*dpsic
        y(i,l)=y(i,l)+dip0(l)*dpsic
40     continue
    endif
    chi=chi+dchi
50 continue
  write(lout,10000) itra,amp,chi0,chid
  write(lout,10010) x(1,1), y(1,1), x(1,2), y(1,2), sigm(1), dps(1), x(2,1), y(2,1), x(2,2), y(2,2), sigm(2), dps(2)
  return
60 itra=2
  do 80 i=1,itra
    sigm(i)=exz(i,5)
    dps(i)=exz(i,6)
    do 70 l=1,2
      ll=(l-1)*2
      x(i,l)=exz(i,1+ll)
      y(i,l)=exz(i,2+ll)
70   continue
80 continue
  write(lout,10020)
  write(lout,10010) x(1,1), y(1,1), x(1,2), y(1,2), sigm(1), dps(1), x(2,1), y(2,1), x(2,2), y(2,2), sigm(2), dps(2)
!-----------------------------------------------------------------------
  return
10000 format(t5,'---- ENTRY ANFB ----/ITRA/',i3,' /AMP/ ',f8.3,2x,f8.3, &
  &' /CHI0,CHID/  ',f6.1,2x,f6.1)
10010 format(/5x,'---- TWIN-TRAJECTORIES NO CL.ORBIT ADDED'/ 5x,'/X1  /'&
  &,f47.33/5x,'/XP1 /',f47.33/ 5x,'/Y1  /',f47.33/5x,'/YP1 /',f47.33/&
  &5x,'/SIG1/',f47.33/5x,'/DP1 /',f47.33/ 5x,'/X2  /',f47.33/5x,     &
  &'/XP2 /',f47.33/ 5x,'/Y2  /',f47.33/5x,'/YP2 /',f47.33/ 5x,       &
  &'/SIG2/',f47.33/5x,'/DP2 /',f47.33/)
10020 format(t5,'---- ENTRY ANFB ----/COORDINATE-INPUT')
10030 format(//131('-')//t10,27('O')/t10,2('O'),23x,2('O')/t10,         &
  &'OO  INITIAL COORDINATES  OO'/ t10,2('O'),23x,2('O')/t10,27('O')  &
  &//131('-')//)
end subroutine anfb

! FixMe: Dummy routine for dynk when compiled as SixDA
subroutine synuthck
  return
end subroutine synuthck
