!-----------------------------------------------------------------------
!  CENTRAL LOOP FOR 6-DIMENSIONAL CLOSED ORBIT
!-----------------------------------------------------------------------
subroutine umlauda

  use floatPrecision
  use physical_constants
  use numerical_constants
  use mathlib_bouncer
  use dump, only : dumpclo, dumptas, dumptasinv, ldump
  use crcoall
  use string_tools
  use mod_units
  use parpro
  use parbeam, only : beam_expflag,beam_expfile_open
  use mod_common
  use mod_commons
  use mod_common_track, only : xxtr,yytr,crois,tasm,comt_daStart,comt_daEnd
  use mod_common_da
  use mod_commond2
  use wire
  use mod_lie_dab, only : idao,iscrri,rscrri,iscrda

  implicit none

  integer i,ibb,iii,i2,i3,icav,icoonly,ien,iflag,iflag1,iflag2,ii,ii2,ip,ipch,irrtr,ivar,ivar1,     &
    iwrite,ix,j,j1,jb,jmel,jx,k,kkk,kpz,kzz,mfile,nd2,nmz,idaa,angno,f,h
  real(kind=fPrec) beamoff1,beamoff2,beamoff3,beamoff4,beamoff5,beamoff6,betr0,cbxb,cbzb,           &
    coefh1,cik,coefh2,coefv1,coefv2,crk,crxb,crzb,cx,dare,det1,dpdav,dps1,dps11,dummy,ed1,ed2,ox,   &
    oxp,oxp1,oz,ozp,ozp1,r0,r2b,r2bf,rb,rbf,rho2b,rkb,rkbf,scikveb,scrkveb,sfac1,sfac2,sfac2s,sfac3,&
    sfac4,sfac5,sigm1,sigmdac,startco,sx,tas,tkb,tl,x2pi,xbb,xrb,xs,zbb,zrb,zs,crabfreq,crabpht,    &
    crabpht2,crabpht3,crabpht4,temp_angle,tan_t,sin_t,cos_t
  integer damap(6),damapi(6),damap1(6),aa2(6),aa2r(6),a1(6),a1r(6),xy(6),df(6),jj(100),i4(10,2)
  real(kind=fPrec) zfeld1(100),zfeld2(100),dpdav2(6),rrad(3),rdd(6,6),dicu(20),angnoe(3),angp(2,6), &
    phi(3),dphi(3),b1(3),b2(3),b3(3),al1(3),al2(3),al3(3),g1(3),g2(3),g3(3),d(3),dp(3),c(3),cp(3),  &
    au(6,6),aui(2)
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  character(len=mNameLen) typ
  integer expertUnit

! For treatment and/or conversion of BEAM parameters in/to the new format
  character(len=25) tmpStr(6)
  logical rErr
  integer wire_num_aux ! auxiliary variable to count number of wires
  save
!-----------------------------------------------------------------------
#include "include/daini.f90"
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
  call comt_daStart
  nd2=ndimf*2
  call etall(damap,6)
  call etall(damapi,6)
  call etall(damap1,6)
  call etall(angno,1)
  call etall(f,1)
  call etall(aa2,6)
  call etall(aa2r,6)
  call etall(a1,6)
  call etall(a1r,6)
  call etall(xy,6)
  call etall(h,1)
  call etall(df,6)
  if(iqmodc.eq.1) call danot(2)
  if(iqmodc.eq.3) call danot(1)
  if(ichromc.eq.1) call danot(3)
  icoonly=0
  if(iqmodc.eq.2.or.iqmodc.eq.4.or.ichromc.eq.2) icoonly=1
  do j=1,2
    angnoe(j)=zero
    do i=1,6
      angp(j,i)=zero
    enddo
  enddo
  do i=1,100
    jj(i)=0
  enddo
  x2pi=atan_mb(one)*eight
  i4(1,1)=1
  i4(1,2)=1
  i4(2,1)=3
  i4(2,2)=3
  i4(3,1)=1
  i4(3,2)=3
  i4(4,1)=1
  i4(4,2)=2
  i4(5,1)=1
  i4(5,2)=4
  i4(6,1)=2
  i4(6,2)=2
  i4(7,1)=2
  i4(7,2)=3
  i4(8,1)=2
  i4(8,2)=4
  i4(9,1)=3
  i4(9,2)=4
  i4(10,1)=4
  i4(10,2)=4
  e0f=sqrt(e0**2-nucm0**2)                                             !hr05
  betr0=sqrt(one-(nucm0/e0)**2)
  ox=xxtr(1,1)
  oxp=yytr(1,1)
  oz=xxtr(1,2)
  ozp=yytr(1,2)
  sigm1=sigm(1)
  dps1=dps(1)
  nucmda=nucm0 ! this is doen like this in case we want to change so we can run Ions designed for different optics.
  !mtc     (j) = (nzz(j)*nucm0)/(zz0*nucm(j)) Again we keep this variable in case we want to change the implementation later.
  mtcda = one


  if(iqmodc.eq.1) then
    if(el(iq(1)).le.pieni) then
      ed1=ed(iq(1))
    else
      ed1=ek(iq(1))
    endif
    if(el(iq(2)).le.pieni) then
      ed2=ed(iq(2))
    else
      ed2=ek(iq(2))
    endif
  endif
  if(ichromc.eq.1) then
    ed1=ed(crois(1))
    ed2=ed(crois(2))
  endif
  call davar(x(1),ox,1)
  oxp1=oxp*(one+dps1)
  call davar(yp(1),oxp1,2)
  ivar=2
  if(nvar2.ge.4) then
    call davar(x(2),oz,3)
    ozp1=ozp*(one+dps1)
    call davar(yp(2),ozp1,4)
    ivar=4
  else
!FOX  X(2)=OZ ;
!FOX  YP(2)=OZP*(ONE+DPS1) ;
  endif
  dps11=dps1*c1e3
  if(nvar2.eq.3) then
    call davar(dpda1,dps11,3)
    ivar=ivar+1
  elseif(nvar2.eq.5) then
    call davar(dpda1,dps11,5)
    ivar=ivar+1
  elseif(nvar2.eq.6) then
    call davar(deltas,zero,5)
    call davar(dpda1,zero,6)

!FOX  DPDA1=DPDA1+DPS1*C1E3 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  EJF1=E0F*(ONE+DPDA)/(NUCM0/NUCMDA) ;
!FOX  EJ1=SQRT(EJF1*EJF1+NUCM0*NUCM0) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DELTAS=DELTAS+SIGM1/RV ;
!FOX  SIGMDA=DELTAS*RV ;

    ivar=ivar+2
  else
!FOX  DELTAS=SIGM1 ;
!FOX  DPDA1=DPS1*C1E3 ;
  endif
  ivar1=ivar
  if(iqmodc.eq.1.or.ichromc.eq.1) then
    call davar(smida(1),ed1,ivar+1)
    call davar(smida(2),ed2,ivar+2)
    ivar=ivar+2
  endif
!--Normal Form Analysis for calculation of linear lattice functions
  if(ilinc.eq.1.or.ilinc.eq.2) then
    mfile=18
!Eric
    rewind mfile
    rewind 26
!ERIC HERE
    call daread(damap,nvar,mfile,one)
    call mapnorm(damap,f,aa2,a1,xy,h,nord1)
    do j=1,nvar
      call dacop(damap(j),damap1(j))
      dummy=dare(damap1(j))
      call dacsu(damap1(j),dummy,damap1(j))
    enddo
    if(ndimf.eq.3) then
      call damul(damap1(5),damap1(5),angno)
      call averaged(angno,damap1,.true.,angno, rv)
      jj(5)=1
      jj(6)=1
      call dapek(angno,jj,emitz)
      jj(5)=0
      jj(6)=0
      if(abs(emitz).le.pieni) then
        emitz=zero
      else
        emitz=((sigz**2/emitz)*half)*c1e6                            !hr05
      endif
    endif
    jj(5)=1
    do j=1,nd2
      call dapek(a1(j),jj,dicu(j))
    enddo
    jj(5)=0
  endif
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  CORROLD(1)=X(1) ;
!FOX  CORROLD(2)=YP(1) ;
!FOX  CORROLD(3)=X(2) ;
!FOX  CORROLD(4)=YP(2) ;
!FOX  CORROLD(5)=DELTAS ;
!FOX  CORROLD(6)=DPDA1 ;
        do 5 kkk=1,6
          dpdav=dare(corrold(kkk))
!FOX  CORROLD(KKK)=CORROLD(KKK)-DPDAV ;
5       continue
!FOX  Y(1)=YP(1)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=YP(2)*MTCDA/(ONE+DPDA) ;

  iflag=0
  iflag1=0
  iflag2=0
!FOX  EJF1=E0F*(ONE+DPDA)/(NUCM0/NUCMDA) ;
!FOX  EJ1=SQRT(EJF1*EJF1+NUCMDA*NUCMDA) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  SIGMDA=DELTAS*RV ;
  if(ithick.eq.1) call envada
  icav=0
  typ='START'
  phi(1)=zero
  phi(2)=zero
  phi(3)=zero
  ibb=0
  wire_num_aux=0
!     start loop over single elements
  do 430 i=1,iu
    if(iqmodc.eq.2.or.iqmodc.eq.4) then
      if(i.eq.niu(1)) then
        do ii=1,2
          ii2=2*ii
          clon(ii2-1)=dare(x(ii))
          clon(ii2)=dare(y(ii))
        enddo
        clon(5)=dare(sigmda)
        clon(6)=dare(dpda)
      endif
    endif
    if(ilinc.eq.1.and.i.eq.1) then
      write(lout,10000) nd2
      if(iprint.eq.1) write(lout,10130)
      write(lout,10010)
      write(lout,10020)
      write(lout,10010)
      tl=zero
#include "include/umlalid.f90"
    endif
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
24         continue
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
25         continue
          if(ivar.gt.ivar1) then
!FOX  CORRAU2(7)=SMIDA(1) ;
!FOX  CORRAU2(8)=SMIDA(2) ;
            dpdav=dare(smida(1))
!FOX  CORRAU1(7)=SMIDA(1)-DPDAV ;
            dpdav=dare(smida(2))
!FOX  CORRAU1(8)=SMIDA(2)-DPDAV ;
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
    if(ix.gt.nblo) goto 50
    if(ix <= 0) then
      write(lerr,"(a)") "UMLAUDA> ERROR Inverted linear blocks not allowed."
      call prror
    endif
#include "include/dalin1.f90"
          ipch=0
          if(iqmodc.eq.1.and.kz(jx).eq.2) then
            if(jx.eq.iq(1).or.iratioe(jx).eq.iq(1)) then
              ipch=1
            else if(jx.eq.iq(2).or.iratioe(jx).eq.iq(2)) then
              ipch=2
            endif
          endif
          if(ipch.ne.0) then
            call envquad(jx,ipch)
!FOX  PUX=X(1) ;
!FOX  PUZ=Y(1) ;
!FOX  X(1)=ALDAQ(1,1)*PUX+ALDAQ(1,2)*PUZ+ALDAQ(1,5)*IDZ(1) ;
!FOX  Y(1)=ALDAQ(1,3)*PUX+ALDAQ(1,4)*PUZ+ALDAQ(1,6)*IDZ(1) ;
!FOX  PUX=X(2) ;
!FOX  PUZ=Y(2) ;
!FOX  X(2)=ALDAQ(2,1)*PUX+ALDAQ(2,2)*PUZ+ALDAQ(2,5)*IDZ(2) ;
!FOX  Y(2)=ALDAQ(2,3)*PUX+ALDAQ(2,4)*PUZ+ALDAQ(2,6)*IDZ(2) ;
          else
#include "include/dalin2.f90"
          endif
#include "include/dalin3.f90"
      if(ilinc.eq.1) then
        do jb=1,jmel
          jx=mtyp(ix,jb)
          typ=bez(jx)
          tl=tl+el(jx)
#include "include/umlalid.f90"
          if(i.eq.nt) goto 470
        enddo
      endif
#include "include/dalin4.f90"
          ipch=0
          if(iqmodc.eq.1.and.kz(jx).eq.2) then
            if(jx.eq.iq(1).or.iratioe(jx).eq.iq(1)) then
              ipch=1
            else if(jx.eq.iq(2).or.iratioe(jx).eq.iq(2)) then
              ipch=2
            endif
          endif
          if(ipch.ne.0) then
            call envquad(jx,ipch)
!FOX  PUX=X(1) ;
!FOX  PUZ=Y(1) ;
!FOX  SIGMDA=SIGMDA+ASDAQ(1,1)+ASDAQ(1,2)*PUX+
!FOX  ASDAQ(1,3)*PUZ+ASDAQ(1,4)*PUX*PUZ+ASDAQ(1,5)*PUX*PUX+
!FOX  ASDAQ(1,6)*PUZ*PUZ ;
!FOX  X(1)=ALDAQ(1,1)*PUX+ALDAQ(1,2)*PUZ+ALDAQ(1,5)*IDZ(1) ;
!FOX  Y(1)=ALDAQ(1,3)*PUX+ALDAQ(1,4)*PUZ+ALDAQ(1,6)*IDZ(1) ;
!FOX  PUX=X(2) ;
!FOX  PUZ=Y(2) ;
!FOX  SIGMDA=SIGMDA+ASDAQ(2,1)+ASDAQ(2,2)*PUX+
!FOX  ASDAQ(2,3)*PUZ+ASDAQ(2,4)*PUX*PUZ+ASDAQ(2,5)*PUX*PUX+
!FOX  ASDAQ(2,6)*PUZ*PUZ ;
!FOX  X(2)=ALDAQ(2,1)*PUX+ALDAQ(2,2)*PUZ+ALDAQ(2,5)*IDZ(2) ;
!FOX  Y(2)=ALDAQ(2,3)*PUX+ALDAQ(2,4)*PUZ+ALDAQ(2,6)*IDZ(2) ;
          else
#include "include/dalin5.f90"
          endif
        else
          if(iexact) then
!-----------------------------------------------------------------------
!  EXACT DRIFT
!-----------------------------------------------------------------------
!FOX  X(1)=X(1)*C1M3 ;
!FOX  X(2)=X(2)*C1M3 ;
!FOX  Y(1)=Y(1)*C1M3 ;
!FOX  Y(2)=Y(2)*C1M3 ;
!FOX  SIGMDA=SIGMDA*C1M3 ;
!FOX  PZ=SQRT(ONE-Y(1)*Y(1)-Y(2)*Y(2)) ;
!FOX  X(1)=X(1)+EL(JX)*(Y(1)/PZ) ;
!FOX  X(2)=X(2)+EL(JX)*(Y(2)/PZ) ;
!FOX  SIGMDA=SIGMDA+(ONE-(RV/PZ))*EL(JX) ;
!FOX  X(1)=X(1)*C1E3 ;
!FOX  X(2)=X(2)*C1E3 ;
!FOX  Y(1)=Y(1)*C1E3 ;
!FOX  Y(2)=Y(2)*C1E3 ;
!FOX  SIGMDA=SIGMDA*C1E3 ;
!-----------------------------------------------------------------------
          else
! Regular drift
#include "include/dalin6.f90"
!FOX  SIGMDA=SIGMDA+
#include "include/sqrtfox.f90"
          endif
        endif
        if(ilinc.eq.1) then
          typ=bez(jx)
          tl=tl+el(jx)
#include "include/umlalid.f90"
          if(i.eq.nt) goto 470
        endif
      enddo
    endif
    goto 430
50   ix=ix-nblo
    if(abs(dare(x(1))) > aint(aper(1)) .or. abs(dare(x(2))) > aint(aper(2))) then
      write(lout,10120)j,i,dare(x(1)),aper(1),dare(x(2)),aper(2),ix,kz(ix),bez(ix)
      write(lerr,"(a)") "UMLAUDA> ERROR Unstable closed orbit in DA calculation."
      call prror
    end if
    kpz=abs(kp(ix))
    if(kpz.ge.0 .and. kpz.lt.6) goto 80
    if(kpz.eq.6) goto 70
    goto 430
70   continue
    if(ition.ne.0) then
!FOX  EJF0=EJF1 ;
      ixcav=ix
      if(abs(dppoff).gt.pieni) then
        sigmdac=dare(sigmda)
        sigmoff(i)=sigmdac
!FOX  SIGMDA=SIGMDA-SIGMDAC ;
      endif
      call synoda

!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;

      if(nvar2.eq.6.and.nsix.ne.2) then
        iflag=1
        iflag1=1
        iflag2=1
      endif
    endif
    goto 440
80   kzz=kz(ix)
    if(kzz.eq.15) then
! the same as in umlalid1
      wire_num_aux = wire_num_aux+1
      if(wire_num_aux.gt.wire_max) then
        write(lerr,"(2(a,i0))") "UMLAUDA> ERROR Maximum number of wires exceeded. Max is ",wire_max,", got ",wire_num_aux
        call prror
      endif
      wire_num(i) = wire_num_aux
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
!FOX  DELTAS=SIGMDA/RV ;
!FOX  DPDA1=DPDA*C1E3 ;
      call dacop(x(1),damap(1))
      call dacop(yp(1),damap(2))
      call dacop(x(2),damap(3))
      call dacop(yp(2),damap(4))
      do j=1,2
        ii=2*j
        call dapek(damap(ii-1),jj,c(j))
        call dapek(damap(ii),jj,cp(j))
      enddo
      call dacsu(damap(1),c(1),damap(1))
      call dacsu(damap(2),cp(1),damap(2))
      call dacsu(damap(3),c(2),damap(3))
      call dacsu(damap(4),cp(2),damap(4))
      if(ndimf.eq.3) then
        call dacop(deltas,damap(5))
        call dacop(dpda1,damap(6))
        call dapek(damap(5),jj,c(3))
        call dapek(damap(6),jj,cp(3))
        call dacsu(damap(5),c(3),damap(5))
        call dacsu(damap(6),cp(3),damap(6))
        if(iflag2.eq.1.and.ithick.eq.1) then
          call dacct(damap,nvar,corrnew,nvar,damap,nvar)
        endif
      endif
      call dainv(damap,nvar,damapi,nvar)
      call dacct(damap,nvar,aa2,nvar,aa2r,nvar)
      call dacct(damap,nvar,damap1,nvar,damap,nvar)
      call dacct(damap,nvar,damapi,nvar,damap,nvar)
! the same as in umlalid1

!FOX  XX(1)=X(1) ;
!FOX  XX(2)=X(2) ;
!FOX  YY(1)=Y(1) ;
!FOX  YY(2)=Y(2) ;
  wire_clo(1,wire_num(i))=dare(x(1))
  wire_clo(2,wire_num(i))=dare(x(2))
  wire_clo(4,wire_num(i))=dare(y(1))*(one+dare(DPDA))
  wire_clo(5,wire_num(i))=dare(y(2))*(one+dare(DPDA))
  if(ndimf.eq.3) then
      wire_clo(3,wire_num(i))=dare(SIGMDA)
      wire_clo(6,wire_num(i))=dare(DPDA)
  endif

  call wireda(ix,i)


!FOX  Y(1)=YY(1) ;
!FOX  Y(2)=YY(2) ;
      goto 440
    endif
    if(ilinc.eq.2.and.kzz.eq.20) then
      if(nbeam.ge.1) then
        ibb=ibb+1
        if(ibb > nbb) then
          write(lerr,"(a,i0)") "UMLAUDA> ERROR Maximum element number for beam-beam with coupling exceeded: nbb = ",nbb
          call prror
        end if
        imbb(i)=ibb
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
!FOX  DPDA1=DPDA*C1E3 ;
!FOX  DELTAS=SIGMDA/RV ;
        call dacop(x(1),damap(1))
        call dacop(yp(1),damap(2))
        call dacop(x(2),damap(3))
        call dacop(yp(2),damap(4))
        do j=1,2
          ii=2*j
          call dapek(damap(ii-1),jj,c(j))
          call dapek(damap(ii),jj,cp(j))
        enddo
        call dacsu(damap(1),c(1),damap(1))
        call dacsu(damap(2),cp(1),damap(2))
        call dacsu(damap(3),c(2),damap(3))
        call dacsu(damap(4),cp(2),damap(4))
        if(ndimf.eq.3) then
          call dacop(deltas,damap(5))
          call dacop(dpda1,damap(6))
          call dapek(damap(5),jj,c(3))
          call dapek(damap(6),jj,cp(3))
          call dacsu(damap(5),c(3),damap(5))
          call dacsu(damap(6),cp(3),damap(6))
          if(iflag2.eq.1.and.ithick.eq.1) then
            call dacct(damap,nvar,corrnew,nvar,damap,nvar)
          endif
        endif
        call dainv(damap,nvar,damapi,nvar)
        call dacct(damap,nvar,aa2,nvar,aa2r,nvar)
        call dacct(damap,nvar,damap1,nvar,damap,nvar)
        call dacct(damap,nvar,damapi,nvar,damap,nvar)
        do ii=1,3
          call damul(damap(i4(ii,1)),damap(i4(ii,2)),angno)
          call averaged(angno,aa2r,.false.,angno,rv)
          do j=1,ndimf
            j1=2*j
            jj(j1-1)=1
            jj(j1)=1
            call dapek(angno,jj,angnoe(j))
            jj(j1-1)=0
            jj(j1)=0
          enddo
          if(ndimf.eq.3) then
            bbcu(ibb,ii)=two*((emitx*angnoe(1)+emity*angnoe(2))+emitz*angnoe(3))
          else
            bbcu(ibb,ii)=two*(emitx*angnoe(1)+emity*angnoe(2))
          endif
        enddo
        if (beam_expflag .eq. 0) then !Old-style input
          if(parbe(ix,2).gt.zero) then
            do ii=4,10
              call damul(damap(i4(ii,1)),damap(i4(ii,2)),angno)
              call averaged(angno,aa2r,.false.,angno, rv)
              do j=1,ndimf
                j1=2*j
                jj(j1-1)=1
                jj(j1)=1
                call dapek(angno,jj,angnoe(j))
                jj(j1-1)=0
                jj(j1)=0
              enddo
              if(ndimf.eq.3) then
                bbcu(ibb,ii) = two * ((emitx*angnoe(1)+emity*angnoe(2))+emitz*angnoe(3))
              else
                bbcu(ibb,ii)=two*(emitx*angnoe(1)+emity*angnoe(2))
              endif
            enddo
          endif
          if(lhc.eq.1) then
            dummy=bbcu(ibb,1)
            bbcu(ibb,1)=bbcu(ibb,2)
            bbcu(ibb,2)=dummy
            dummy=bbcu(ibb,4)
            bbcu(ibb,4)=bbcu(ibb,9)
            bbcu(ibb,9)=dummy
            dummy=bbcu(ibb,5)
            bbcu(ibb,5)=bbcu(ibb,7)
            bbcu(ibb,7)=dummy
            dummy=bbcu(ibb,6)
            bbcu(ibb,6)=bbcu(ibb,10)
            bbcu(ibb,10)=dummy
          endif
          if(lhc.eq.2) then
            bbcu(ibb,1)=bbbx(ix)
            bbcu(ibb,2)=bbby(ix)
            bbcu(ibb,3)=bbbs(ix)
          endif
        else if (beam_expflag .eq. 1) then !New style input
          if(parbe(ix,2).gt.zero) then
             bbcu(ibb,1)=parbe(ix,7)
             bbcu(ibb,4)=parbe(ix,8)
             bbcu(ibb,6)=parbe(ix,9)
             bbcu(ibb,2)=parbe(ix,10)
             bbcu(ibb,9)=parbe(ix,11)
             bbcu(ibb,10)=parbe(ix,12)
             bbcu(ibb,3)=parbe(ix,13)
             bbcu(ibb,5)=parbe(ix,14)
             bbcu(ibb,7)=parbe(ix,15)
             bbcu(ibb,8)=parbe(ix,16)
          endif
          if(parbe(ix,2).eq.zero) then
             bbcu(ibb,1)=parbe(ix,1)
             bbcu(ibb,2)=parbe(ix,3)
             bbcu(ibb,3)=parbe(ix,13)
          endif
        else
           write(lerr,"(a,i0,a)") "UMLAUDA> ERROR beam_expflag was ",beam_expflag,", expected 0 or 1. This is a BUG!"
           call prror
        end if

        if (.not.beam_expfile_open) then
          call f_requestUnit("beam_expert.txt",expertUnit)
          call f_open(unit=expertUnit,file="beam_expert.txt",formatted=.true.,mode="w",status="replace")
          beam_expfile_open = .true.
          !This line will be a comment if copy-pasted into fort.3
          write(expertUnit,"(a,g13.6,a,g13.6,a,g13.6,a)") "/ ******* USING emitx=",emitx,", emity=",emity,", emitz=",emitz," ******"
        endif

        if(parbe(ix,2).eq.0.0) then !4D
          !Note: One should always use the CRLIBM version when converting,
          ! in order to guarantee the exact same results from the converted input file.
          call chr_fromReal(bbcu(ibb,1),tmpStr(1),19,2,rErr)
          call chr_fromReal(bbcu(ibb,2),tmpStr(2),19,2,rErr)
          call chr_fromReal(parbe(ix,5),tmpStr(3),19,2,rErr)
          call chr_fromReal(parbe(ix,6),tmpStr(4),19,2,rErr)
          call chr_fromReal(ptnfac(ix),tmpStr(5),19,2,rErr)
          call chr_fromReal(bbcu(ibb,3),tmpStr(6),19,2,rErr)
          write(expertUnit,"(a48,1x,a1,6(1x,a25))")  bez(ix),"0",tmpStr(1),tmpStr(2),tmpStr(3),tmpStr(4),tmpStr(5),tmpStr(6)
        else                      ! 6D
          call chr_fromReal(parbe(ix,1),tmpStr(1),19,2,rErr)
          call chr_fromReal(parbe(ix,3),tmpStr(2),19,2,rErr)
          call chr_fromReal(parbe(ix,5),tmpStr(3),19,2,rErr)
          call chr_fromReal(parbe(ix,6),tmpStr(4),19,2,rErr)
          write(expertUnit,"(a48,1x,i4,1x,4(1x,a25))") bez(ix),int(parbe(ix,2)),tmpStr(1),tmpStr(2),tmpStr(3),tmpStr(4)

          call chr_fromReal(bbcu(ibb,1),tmpStr(1),19,2,rErr)
          call chr_fromReal(bbcu(ibb,4),tmpStr(2),19,2,rErr)
          call chr_fromReal(bbcu(ibb,6),tmpStr(3),19,2,rErr)
          call chr_fromReal(bbcu(ibb,2),tmpStr(4),19,2,rErr)
          call chr_fromReal(bbcu(ibb,9),tmpStr(5),19,2,rErr)
          write(expertUnit,"(5(a25,1x))") tmpStr(1),tmpStr(2),tmpStr(3),tmpStr(4),tmpStr(5)

          call chr_fromReal(bbcu(ibb,10),tmpStr(1),19,2,rErr)
          call chr_fromReal(bbcu(ibb,3),tmpStr(2),19,2,rErr)
          call chr_fromReal(bbcu(ibb,5),tmpStr(3),19,2,rErr)
          call chr_fromReal(bbcu(ibb,7),tmpStr(4),19,2,rErr)
          call chr_fromReal(bbcu(ibb,8),tmpStr(5),19,2,rErr)
          call chr_fromReal(ptnfac(ix), tmpStr(6),19,2,rErr)
          write(expertUnit,"(6(a25,1x))") tmpStr(1),tmpStr(2),tmpStr(3),tmpStr(4),tmpStr(5),tmpStr(6)
        endif !END if(parbe(ix,2).eq.0.0)

        if((bbcu(ibb,1).le.pieni).or.(bbcu(ibb,2).le.pieni)) then
          goto 9088
        endif
        if(ibbc.eq.1) then
          sfac1=bbcu(ibb,1)+bbcu(ibb,2)
          sfac2=bbcu(ibb,1)-bbcu(ibb,2)
          sfac2s=one
          if(sfac2.lt.zero) sfac2s=-one                            !hr08
          sfac3=sqrt(sfac2**2+(four*bbcu(ibb,3))*bbcu(ibb,3))          !hr03
          if(sfac3 > sfac1) then
            write(lerr,"(a)") "UMLAUDA> ERROR 6D beam-beam with tilt not possible."
            call prror
          end if
          sfac4=(sfac2s*sfac2)/sfac3                                   !hr03
          sfac5=(((-one*sfac2s)*two)*bbcu(ibb,3))/sfac3                !hr03
          sigman(1,ibb)=sqrt(((sfac1+sfac2*sfac4)+(two*bbcu(ibb,3))*sfac5)*half)    !hr03
          sigman(2,ibb)=sqrt(((sfac1-sfac2*sfac4)-(two*bbcu(ibb,3))*sfac5)*half)    !hr03
          bbcu(ibb,11)=sqrt(half*(one+sfac4))
          bbcu(ibb,12)=(-one*sfac2s)*sqrt(half*(one-sfac4))            !hr03
          if(bbcu(ibb,3).lt.zero) bbcu(ibb,12)=-one*bbcu(ibb,12)       !hr03
        else
          bbcu(ibb,11)=one
          sigman(1,ibb)=sqrt(bbcu(ibb,1))
          sigman(2,ibb)=sqrt(bbcu(ibb,2))
        endif
        if(parbe(ix,2).gt.zero) then !6D -> convert units
          do ii=1,10
            bbcu(ibb,ii)=bbcu(ibb,ii)*c1m6
          enddo
        endif
      endif
      goto 440
    endif
    if(kzz.eq.20.and.iqmodc.eq.4) goto 440
    if(kzz.eq.20.and.parbe(ix,2).eq.zero) then                        !hr12
      if(nbeam.ge.1) then
        if(ilinc.eq.0) then
          clobeam(1,imbb(i))=dare(x(1))
          clobeam(2,imbb(i))=dare(x(2))
          clobeam(4,imbb(i))=dare(y(1))*(one+dare(dpda))
          clobeam(5,imbb(i))=dare(y(2))*(one+dare(dpda))
          if(ndimf.eq.3) then
            clobeam(3,imbb(i))=dare(sigmda)
            clobeam(6,imbb(i))=dare(dpda)
          endif
        endif
        if(sigman(1,imbb(i)).eq.sigman(2,imbb(i))) then
          if(ibeco.eq.1) then
            if(ibbc.eq.0) then
              crk=parbe(ix,5)
              cik=parbe(ix,6)
            else
              crk=parbe(ix,5)*bbcu(imbb(i),11) + parbe(ix,6)*bbcu(imbb(i),12)
              cik=parbe(ix,6)*bbcu(imbb(i),11) - parbe(ix,5)*bbcu(imbb(i),12)
            endif
            rho2b=crk**2+cik**2                                          !hr03
            if(rho2b.gt.pieni) then
              if(abs(sigman(1,imbb(i))).lt.pieni) goto 9088
              tkb=rho2b/((two*sigman(1,imbb(i)))*sigman(1,imbb(i)))        !hr03
              beamoff4=(((crad*ptnfac(ix))*crk)/rho2b)*(one-exp_mb(-one*tkb)) !hr03
              beamoff5=(((crad*ptnfac(ix))*cik)/rho2b)*(one-exp_mb(-one*tkb)) !hr03
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
            rkb=((crad*ptnfac(ix))*pisqrt)/rb                            !hr03
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
            tkb=(crk**2/sigman(1,imbb(i))**2+cik**2/sigman(2,imbb(i))**2)*half  !hr03
            xbb=(sigman(2,imbb(i))/sigman(1,imbb(i)))*xrb                !hr03
            zbb=(sigman(1,imbb(i))/sigman(2,imbb(i)))*zrb                !hr03
            call errf(xbb,zbb,cbxb,cbzb)
            beamoff4=(rkb*(crzb-exp_mb(-one*tkb)*cbzb))*sign(one,crk)
            beamoff5=(rkb*(crxb-exp_mb(-one*tkb)*cbxb))*sign(one,cik)
          endif
          if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
          r2bf=two*(sigman(1,imbb(i))**2-sigman(2,imbb(i))**2) !hr08
          rbf=sqrt(r2bf)
          rkbf=((crad*ptnfac(ix))*pisqrt)/rbf                          !hr03
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
            r2b=two*(sigman(2,imbb(i))**2-sigman(1,imbb(i))**2)   !hr08
            rb=sqrt(r2b)
            rkb=((crad*ptnfac(ix))*pisqrt)/rb                            !hr03
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
            tkb=(crk**2/sigman(1,imbb(i))**2+cik**2/sigman(2,imbb(i))**2)*half  !hr03
            xbb=(sigman(2,imbb(i))/sigman(1,imbb(i)))*xrb                !hr03
            zbb=(sigman(1,imbb(i))/sigman(2,imbb(i)))*zrb                !hr03
            call errf(zbb,xbb,cbzb,cbxb)
            beamoff4=(rkb*(crzb-exp_mb(-one*tkb)*cbzb))*sign(one,crk)
            beamoff5=(rkb*(crxb-exp_mb(-one*tkb)*cbxb))*sign(one,cik)
          endif
          if(abs(sigman(1,imbb(i))).lt.pieni.or.abs(sigman(2,imbb(i))).lt.pieni) goto 9088
          r2bf=two*(sigman(2,imbb(i))**2-sigman(1,imbb(i))**2) !hr08
          rbf=sqrt(r2bf)
          rkbf=((crad*ptnfac(ix))*pisqrt)/rbf                          !hr03
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
        goto 440
      endif
      goto 440
    endif
    if(kzz.eq.20.and.parbe(ix,2).gt.zero) then                        !hr12
      if(ilinc.eq.0)then
        clobeam(1,imbb(i))=dare(x(1))
        clobeam(2,imbb(i))=dare(x(2))
        clobeam(4,imbb(i))=dare(y(1))*(one+dare(dpda))
        clobeam(5,imbb(i))=dare(y(2))*(one+dare(dpda))
        if(ndimf.eq.3)then
          clobeam(3,imbb(i))=dare(sigmda)
          clobeam(6,imbb(i))=dare(dpda)
        endif
      endif
#include "include/beam6dfi.f90"
      goto 440
    endif
    if(kzz.eq.41) then
#include "include/alignf.f90"
#include "include/rfmulti_fox.f90"
      goto 440
    endif
    if(kzz.eq.43) then
      temp_angle = ed(ix)
#include "include/xrot_fox.f90"
      goto 440
    endif
    if(kzz.eq.44) then
      temp_angle = ed(ix)
#include "include/yrot_fox.f90"
      goto 440
    endif
    if(kzz.eq.45) then
      temp_angle = ed(ix)
#include "include/srot_fox.f90"
      goto 440
    endif




    if(kzz.eq.23) then
!FOX  CRABAMP=ED(IX)*QQ0 ;

      crabfreq=ek(ix)*c1e3
      crabpht=crabph(ix)


!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT) ;

!FOX  Y(1)=Y(1) - CRABAMP*C1E3/E0F*
!FOX  SIN((SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT))/(1+DPDA) ;

!FOX  EJ1=EJ1-CRABAMP
!FOX  *CRABFREQ*2D0*PI/(CLIGHT)*X(1)*
!FOX  COS((SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT)) ;

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1=(EJF1-E0F)/E0F*C1E3 ;

!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;

      goto 440
  endif
    if(kzz.eq.-23) then
!FOX  CRABAMP=ED(IX)*QQ0 ;
        crabfreq=ek(ix)*c1e3
        crabpht=crabph(ix)
!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT) ;

!FOX  Y(2)=Y(2) - CRABAMP*C1E3/E0F*
!FOX  SIN(KCRABDA)*MTCDA/(ONE+DPDA) ;

!FOX  EJ1=EJ1-CRABAMP
!FOX  *CRABFREQ*2D0*PI/(CLIGHT)*X(2)*
!FOX  COS(KCRABDA) ;

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
!FOX  CRABAMP2=ED(IX)*QQ0 ;

    crabfreq=ek(ix)*c1e3 !JBG Input in MHz changed to kHz
    crabpht2=crabph2(ix)
!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT2) ;

!FOX  Y(1)=Y(1) + (CRABAMP2*CRKVE)*
!FOX  COS(KCRABDA)*MTCDA/(ONE+DPDA);
!FOX  Y(2)=Y(2) - (CRABAMP2*CIKVE)*
!FOX  COS(KCRABDA)*MTCDA/(ONE+DPDA);
!FOX  EJ1=EJ1 - (0.5D0)*(CRABAMP2)*(CRKVE*CRKVE-
!FOX  CIKVE*CIKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*E0F*C1M3*
!FOX  SIN(KCRABDA) ;



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
            goto 440
        endif
      xs=xsi(i) ! JBG change of variables for misal calculations
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP2=ED(IX)*QQ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht2=crabph2(ix)
!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT2) ;

!FOX  Y(2)=Y(2) + (CRABAMP2*CRKVE)*
!FOX  COS(KCRABDA)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=Y(1) + (CRABAMP2*CIKVE)*
!FOX  COS(KCRABDA)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  EJ1=EJ1 -(CRABAMP2)*(CIKVE*CRKVE)
!FOX  *(((CRABFREQ*2D0)*PI)/CLIGHT)*E0F*C1M3*
!FOX  SIN(KCRABDA) ;

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
            goto 440
        endif
      xs=xsi(i)
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP3=ED(IX)*QQ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht3=crabph3(ix)
!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT3) ;

!FOX  Y(1)=Y(1) + 2*(0.5D0)*CRABAMP3*((CRKVE*CRKVE)-
!FOX  (CIKVE*CIKVE))*C1M3*MTCDA/(ONE+DPDA)*
!FOX  COS(KCRABDA);
!FOX  Y(2)=Y(2) - 2*CRABAMP3*(CRKVE*CIKVE)*C1M3*MTCDA/(ONE+DPDA)*
!FOX  COS(KCRABDA);
!FOX  EJ1=EJ1 - 2*(1/6.)*(CRABAMP3)*(CRKVE*CRKVE*CRKVE-
!FOX  3*CIKVE*CIKVE*CRKVE)*(((CRABFREQ*2D0)*PI)/CLIGHT)*
!FOX  C1M6*E0F*
!FOX  SIN(KCRABDA);

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA*MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
      goto 440
      endif
      if(kzz.eq.-27) then
        ! JBG bypass this element if 4D/5D case
        if(iclo6.eq.0) then
            goto 440
        endif
      xs=xsi(i)
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP3=ED(IX)*QQ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht3=crabph3(ix)
!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT3) ;

!FOX  Y(2)=Y(2) - (CRABAMP3*C1M3*
!FOX  COS(KCRABDA)*(MTCDA/(ONE+DPDA))*
!FOX  ((CIKVE*CIKVE)-(CRKVE*CRKVE))) ;
!FOX  Y(1)=Y(1) + 2D0*CRABAMP3*(CRKVE*CIKVE)*C1M3*(MTCDA/(ONE+DPDA))*
!FOX  COS(KCRABDA);
!FOX  EJ1=EJ1 + (ONE/3.0)*(CRABAMP3)*(CIKVE*CIKVE*CIKVE-
!FOX  3.0*CIKVE*CRKVE*CRKVE)*(((CRABFREQ*2.0)*PI)/CLIGHT)*
!FOX  E0F*C1M6*
!FOX  SIN(KCRABDA);

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
            goto 440
        endif
      xs=xsi(i)
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP4=ED(IX)*QQ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht4=crabph4(ix)
!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT) ;

!FOX  Y(1)=Y(1) + 6*(1D0/6D0)*(CRABAMP4)*
!FOX  (CRKVE*CRKVE*CRKVE-(3*CRKVE*CIKVE*CIKVE))*C1M6*
!FOX  COS(KCRABDA)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=Y(2) - 6*(1D0/6D0)*(CRABAMP4)*
!FOX  (3*CIKVE*CRKVE*CRKVE-CIKVE*CIKVE*CIKVE)*C1M6*
!FOX  COS(KCRABDA)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  EJ1=EJ1-6*(1D0/24D0)*(CRABAMP4)*(CRKVE*CRKVE*CRKVE*CRKVE-
!FOX  6*CRKVE*CRKVE*CIKVE*CIKVE+CIKVE*CIKVE*CIKVE*CIKVE)*
!FOX  C1M9*(((CRABFREQ*2D0)*PI)/CLIGHT)*E0F*
!FOX  SIN(KCRABDA) ;


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
      if(kzz.eq.-28) then
        ! JBG bypass this element if 4D/5D case
        if(iclo6.eq.0) then
            goto 440
        endif
      xs=xsi(i)
      zs=zsi(i)
#include "include/alignf.f90"
!FOX  CRABAMP4=ED(IX)*QQ0 ;
          crabfreq=ek(ix)*c1e3
          crabpht4=crabph4(ix)
!FOX  KCRABDA=(SIGMDA/(CLIGHT*(E0F/E0))
!FOX  *CRABFREQ*2D0*PI + CRABPHT) ;

!FOX  Y(1)=Y(1) - (CRABAMP4)*
!FOX  (CIKVE*CIKVE*CIKVE-(3D0*CIKVE*CRKVE*CRKVE))*C1M6*
!FOX  COS(KCRABDA)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=Y(2) - (CRABAMP4)*
!FOX  (3D0*CRKVE*CIKVE*CIKVE-CRKVE*CRKVE*CRKVE)*C1M6*
!FOX  COS(KCRABDA)
!FOX  *MTCDA/(ONE+DPDA) ;
!FOX  EJ1=EJ1-(CRABAMP4)*(CRKVE*CRKVE*CRKVE*CIKVE-
!FOX  CIKVE*CIKVE*CIKVE*CRKVE)*E0F*
!FOX  C1M9*(((CRABFREQ*2D0)*PI)/CLIGHT)*
!FOX  SIN(KCRABDA) ;

!FOX  EJF0=EJF1 ;
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1 = (EJF1-E0F)/E0F*C1E3 ;
!FOX  RV=EJ1/E0*E0F/EJF1 ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=EJF0/EJF1*Y(1) ;
!FOX  Y(2)=EJF0/EJF1*Y(2) ;
      endif
    if(kzz.eq.22) then ! Phase Trombone
#include "include/trombone_fox.f90"
    end if
    if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 440
    if(kzz.eq.15) goto 440
    ipch=0
    if(iqmodc.eq.1) then
      if(ix.eq.iq(1).or.iratioe(ix).eq.iq(1)) then
        ipch=1
      else if(ix.eq.iq(2).or.iratioe(ix).eq.iq(2)) then
        ipch=2
      endif
    endif
    if(ichromc.eq.1) then
      if(ix.eq.crois(1).or.iratioe(ix).eq.crois(1)) then
        ipch=1
      else if(ix.eq.crois(2).or.iratioe(ix).eq.crois(2)) then
        ipch=2
      endif
    endif
    if(ipch.ne.0) then
!FOX  EKK=(SMIDA(IPCH)*RATIOE(IX)+SMIZF(I))*MTCDA/(ONE+DPDA) ;
    else
!FOX  EKK=SMI(I)*MTCDA/(ONE+DPDA) ;
    endif
    xs=xsi(i)
    zs=zsi(i)
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
      if(abs(r0).le.pieni.or.nmz.eq.0) goto 440
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
      goto 440
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

440  continue
  if(ilinc.eq.1) then
    typ=bez(ix)
#include "include/umlalid.f90"
    if(i.eq.nt) goto 470
  endif
430  continue ! END LOOP OVER SINGLE ELEMENTS IN UMLAUDA

!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
!FOX  DELTAS=SIGMDA/RV ;
  if(icav.eq.0.or.ithick.ne.1) then
    if(icoonly.eq.1) then
      xxtr(1,1) = dare(x(1))
      yytr(1,1) = dare(y(1))
      xxtr(1,2) = dare(x(2))
      yytr(1,2) = dare(y(2))
      sigm(1) = dare(sigmda)
      dps(1) = dare(dpda)
    endif
  else
!FOX  CORRAU1(1)=X(1) ;
!FOX  CORRAU1(2)=YP(1) ;
!FOX  CORRAU1(3)=X(2) ;
!FOX  CORRAU1(4)=YP(2) ;
!FOX  CORRAU1(5)=DELTAS ;
!FOX  CORRAU1(6)=DPDA1 ;
    do 435 kkk=1,6
      dpdav2(kkk)=dare(corrau1(kkk))
!FOX  CORRAU1(KKK)=CORRAU1(KKK)-DPDAV2(KKK) ;
435   continue
    if(ivar.gt.ivar1) then
!FOX  CORRAU1(7)=SMIDA(1) ;
!FOX  CORRAU1(8)=SMIDA(2) ;
      dpdav=dare(smida(1))
!FOX  CORRNEW(7)=SMIDA(1)-DPDAV ;
      dpdav=dare(smida(2))
!FOX  CORRNEW(8)=SMIDA(2)-DPDAV ;
    endif
    call dacct(corrau1,nvar,corrnew,nvar,corrau2,nvar)
    do 436 kkk=1,6
!FOX  CORRAU2(KKK)=CORRAU2(KKK)+DPDAV2(KKK) ;
436   continue
!FOX  CORRAU1(2)=CORRAU2(2)/(ONE+CORRAU2(6)) ;
!FOX  CORRAU1(4)=CORRAU2(4)/(ONE+CORRAU2(6)) ;
!FOX  X(1)=CORRAU2(1) ;
!FOX  YP(1)=CORRAU2(2) ;
!FOX  X(2)=CORRAU2(3) ;
!FOX  YP(2)=CORRAU2(4) ;
!FOX  DELTAS=CORRAU2(5) ;
!FOX  DPDA1=CORRAU2(6) ;
!FOX  DPDA=DPDA1*C1M3 ;
!FOX  MOIDA=MTCDA/(ONE+DPDA) ;
!FOX  Y(1)=YP(1)*MTCDA/(ONE+DPDA) ;
!FOX  Y(2)=YP(2)*MTCDA/(ONE+DPDA) ;
!FOX  SIGMDA=DELTAS*RV;
    if(icoonly.eq.1) then
      xxtr(1,1) = dare(x(1))
      yytr(1,1) = dare(y(1))
      xxtr(1,2) = dare(x(2))
      yytr(1,2) = dare(y(2))
      sigm(1) = dare(sigmda)
      dps(1) = dare(dpda)
    endif
  endif
  call dacop(x(1),damap(1))
  call dacop(x(2),damap(3))
  if(ndimf.eq.3) call dacop(sigmda,damap(5))
  if(icoonly.eq.1.or.iqmodc.eq.3) then
    call dacop(y(1),damap(2))
    call dacop(y(2),damap(4))
    if(ndimf.eq.3) call dacop(dpda,damap(6))
    do i=1,nd2
      jj(i)=1
      do ii=1,nd2
        call dapek(damap(ii),jj,aml6(ii,i))
        if(i.eq.6) aml6(ii,i)=aml6(ii,i)*c1e3
      enddo
      jj(i)=0
    enddo
    do i=1,nd2
      aml6(i,i)=aml6(i,i)-one
    enddo
  endif
  call dacop(yp(1),damap(2))
  call dacop(yp(2),damap(4))
  if(ndimf.eq.3) then
    call dacop(deltas,damap(5))
    call dacop(dpda1,damap(6))
  else
    call dacop(dpda1,damap(5))
  endif
  if(iqmodc.eq.2.or.iqmodc.eq.4.or.ilin.ge.2) then
    rewind 18
!Eric
    rewind 26
    call daprid(damap,1,nvar,18)
  endif
!--now do the output
  if(iqmodc.eq.1) call danot(3)
  if(iqmodc.eq.3) call danot(2)
  if(ichromc.eq.1) call danot(4)
  if(ilinc.eq.1.or.ilinc.eq.2.or.iqmodc.eq.1.or.iqmodc.eq.3.or.     &
  &ichromc.eq.1) then
    call mapnorm(damap,f,aa2,a1,xy,h,nord1)
  endif
  if(iqmodc.eq.1.or.iqmodc.eq.3) then
    call gettura(wxys,rrad)
    wxys(3)=abs(wxys(3))
    write(lout,*) (wxys(i),i=1,ndimf)
    do i=1,nd2
      jj(i)=1
      do ii=1,nd2
        call dapek(aa2(ii),jj,tas)
        if(i.eq.6.and.ii.ne.6) tas=tas*c1e3
        if(ii.eq.6.and.i.ne.6) tas=tas*c1m3
        tasm(ii,i)=tas
      enddo
      jj(i)=0
    enddo
  endif
  if(iqmodc.eq.1) then
    call dhdj(h,df)
    do i=1,ndimf
      call dapek(df(ndimf+i),jj,corr(1,i))
    enddo
    corr(1,3)=abs(corr(1,3))
    jj(nd2+1)=1
    call dapek(df(ndimf+1),jj,coefh1)
    call dapek(df(ndimf+2),jj,coefv1)
    jj(nd2+1)=0
    jj(nd2+2)=1
    call dapek(df(ndimf+1),jj,coefh2)
    call dapek(df(ndimf+2),jj,coefv2)
    jj(nd2+2)=0
    det1=coefh1*coefv2-coefv1*coefh2
    if(abs(det1) <= pieni) then
      write(lerr,"(a)") "UMLAUDA> ERROR Quadrupoles are not suited to adjust the tunes."
      call prror
    end if
    corr(2,1)=coefv2/det1
    corr(2,2)=(-one*coefh2)/det1                                     !hr05
    corr(3,1)=(-one*coefv1)/det1                                     !hr05
    corr(3,2)=coefh1/det1
  endif
  if(ichromc.eq.1) then
    call dhdj(h,df)
    jj(nd2+1)=1
    call dapek(df(ndimf+1),jj,corr(1,1))
    call dapek(df(ndimf+2),jj,corr(1,2))
    jj(nd2+2)=1
    call dapek(df(ndimf+1),jj,coefh1)
    call dapek(df(ndimf+2),jj,coefv1)
    jj(nd2+2)=0
    jj(nd2+3)=1
    call dapek(df(ndimf+1),jj,coefh2)
    call dapek(df(ndimf+2),jj,coefv2)
    jj(nd2+3)=0
    jj(nd2+1)=0
    det1=coefh1*coefv2-coefv1*coefh2
    if(abs(det1) <= pieni) then
      write(lerr,"(a)") "UMLAUDA> ERROR Sextupoles are not suited to adjust the chromaticity."
      call prror
    end if
    corr(2,1)=coefv2/det1
    corr(2,2)=(-one*coefh2)/det1                                     !hr05
    corr(3,1)=(-one*coefv1)/det1                                     !hr05
    corr(3,2)=coefh1/det1
  endif
470  continue

  call dadal(damap,6)
  call dadal(damapi,6)
  call dadal(damap1,6)
  call dadal(angno,1)
  call dadal(f,1)
  call dadal(aa2,6)
  call dadal(aa2r,6)
  call dadal(a1,6)
  call dadal(a1r,6)
  call dadal(xy,6)
  call dadal(h,1)
  call dadal(df,6)
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  call comt_daEnd
  return

9088 continue
  write(lerr,"(a)") "UMLAUDA> ERROR Either normalized emittances or the resulting sigma values equal to zero for beam-beam/"
  call prror
  return
!-----------------------------------------------------------------------
10000 format(/t5 ,'---- ENTRY ',i1,'D LINOPT ----')
10010 format(133('-'))
10020 format('  NR     TYP      L-TOTAL    P     PHI          ',        &
  &'BETA         ALFA         GAMMA        DIS        DISP         ',&
  &'CLO        CLOP'/ 1x,                                            &
  &'                    (M)           (2*PI)        ',               &
  &'(M)          (RAD)         (M)         (M)        (RAD)        ',&
  &'(MM)       (MRAD)')
10030 format('|',i6,'|',a8,'|',f12.5,'|','X','|',f12.7,'|',f13.6,'|',   &
  &f13.7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10040 format('|',6x,'|',8x,'|',12x,'|','Y','|',12x,'|',f13.6,'|', f13.7,&
  &'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10050 format('|',6x,'|',a8,'|',12x,'|','S','|',12x,'|',f13.6,'|', f13.  &
  &7,'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10055 format('|',6x,'|',a8,'|',12x,'|','Y','|',12x,'|',f13.6,'|', f13.  &
  &7,'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10060 format('|',6x,'|',8x,'|',12x,'|',103('-'))
10070 format('|',6x,'|',8x,'|',12x,'|','Y','|',f12.7,'|',f13.6,'|', f13.&
  &7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10080 format('|',6x,'|',8x,'|',12x,'|','X','|',12x,'|',f13.6,'|', f13.7,&
  &'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10090 format('|',6x,'|',8x,'|',12x,'|','S','|',12x,'|',f13.6,'|', f13.7,&
  &'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10100 format('|',6x,'|',8x,'|',12x,'|','S','|',f12.7,'|',ES13.6,'|',    &
  &f13.7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10110 format(/t10,'CO-TRACKING ENDED ABNORMALLY'/t10, 'PARTICLE NO. '   &
  &,i7,' AT ELEMENT ',i4/ t10,'HORIZ:  AMPLITUDE = ',ES23.16,        &
  &'   APERTURE = ',f15.3/ t10,'VERT:   AMPLITUDE = ',ES23.16,       &
  &'   APERTURE = ',f15.3/ t10,'ELEMENT - LIST NUMBER ',i4,          &
  &' TYP NUMBER ',i4,' NAME ',a16/)
10120 format(/t10,'CO-TRACKING ENDED ABNORMALLY'/t10, 'PARTICLE NO. '   &
  &,i7,' AT ELEMENT ',i4/ t10,'HORIZ:  AMPLITUDE = ',ES23.16,        &
  &'   APERTURE = ',f15.3/ t10,'VERT:   AMPLITUDE = ',ES23.16,       &
  &'   APERTURE = ',f15.3/ t10,'ELEMENT - LIST NUMBER ',i4,          &
  &' TYP NUMBER ',i4,' NAME ',a16/)
10130 format('  LINEAR OPTICS CALCULATION WITH PRINTOUT ',              &
  &'AFTER EACH BLOCK'/                                               &
  &'   A T T E N T I O N : BETATRON PHASE CALCULATION MIGHT BE WRONG'&
  &,' BY A MULTIPLE OF 0.5 FOR EACH LARGE BLOCK'/)
end subroutine umlauda

!-----------------------------------------------------------------------
!  CALCULATION OF : MOMENTUM-DEPENDING ELEMENT-MATRICES AND
!                   CHANGE OF PATH LENGTHS FOR EACH PARTICLE.
!      SPECIALLY PREPARED FOR NEW D.A. (SIX-DIMENSIONAL VERSION)
!-----------------------------------------------------------------------
subroutine envada
  ! Replaced computed goto with select case, VKBO 27/11/2017
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_commond2
  use mod_lie_dab, only : idao,rscrri,iscrda
  implicit none
  integer i,ien,ih,ip,l,idaa
  real(kind=fPrec) dare,result
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  save
!FOX  B D ;
#include "include/dainicom.f90"
!FOX  D V DA INT FOKQ NORD NVAR ; D V DA INT WFHI NORD NVAR ;
!FOX  D V DA INT DPD NORD NVAR ; D V DA INT DPSQ NORD NVAR ;
!FOX  D V DA INT FOK NORD NVAR ; D V DA INT RHO NORD NVAR ;
!FOX  D V DA INT FOK1 NORD NVAR ; D V DA INT SM1 NORD NVAR ;
!FOX  D V DA INT SM2 NORD NVAR ; D V DA INT SM3 NORD NVAR ;
!FOX  D V DA INT SM4 NORD NVAR ; D V DA INT SM5 NORD NVAR ;
!FOX  D V DA INT SM6 NORD NVAR ; D V DA INT SM12 NORD NVAR ;
!FOX  D V DA INT SM23 NORD NVAR ; D V DA INT AS3 NORD NVAR ;
!FOX  D V DA INT AS4 NORD NVAR ; D V DA INT AS6 NORD NVAR ;
!FOX  D V DA INT SI NORD NVAR ; D V DA INT CO NORD NVAR ;
!FOX  D V DA INT G NORD NVAR ; D V DA INT GL NORD NVAR ;
!FOX  D V DA INT SIQ NORD NVAR ; D V DA INT RHOC NORD NVAR ;
!FOX  D V DA INT HI NORD NVAR ; D V DA INT FI NORD NVAR ;
!FOX  D V DA INT AEK NORD NVAR ; D V DA INT HI1 NORD NVAR ;
!FOX  D V DA INT HP NORD NVAR ; D V DA INT HM NORD NVAR ;
!FOX  D V DA INT HC NORD NVAR ; D V DA INT HS NORD NVAR ;
!FOX  D V DA INT FOKC NORD NVAR ; D V DA INT WF NORD NVAR ;
!FOX  D V DA INT AFOK NORD NVAR ; D V DA INT RHOI NORD NVAR ;
!FOX  D V DA INT WFA NORD NVAR ; D V RE INT RATIOE NELE ;
!FOX  D V RE INT EL NELE ; D V RE INT EK NELE ; D V RE INT ED NELE ;
!FOX  D V RE INT ONE ; D V RE INT ZERO ; D V RE INT TWO ;
!FOX  D V RE INT HALF ; D V RE INT FOUR ; D V RE INT C1E3 ;
!FOX  D V RE INT C2E3 ; D V RE INT C4E3 ;
!FOX  D V IN INT I ; D V IN INT L ; D V IN INT IH ; D V IN INT NE ;
!FOX  D V IN INT NA ; D V IN INT IP ; D V IN INT IPCH ;
!FOX  D F RE DARE 1 ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
!FOX  DPD=ONE+DPDA ;
!FOX  DPSQ=SQRT(DPD) ;
  do i=1,il
    do ih=1,2
      do ip=1,6
!FOX  ALDA(IH,IP)=ZERO ;
!FOX  ASDA(IH,IP)=ZERO ;
      end do
    end do

    !-----------------------------------------------------------------------
    if(abs(el(i)).le.pieni) goto 190
    select case (kz(i))
    case (0)
      goto 20

    !-----------------------------------------------------------------------
    !  RECTANGULAR MAGNET
    !-----------------------------------------------------------------------
    case (1,4)
      if (kz(i) == 1) then
        ih = 1
      else
        ih = 2
      end if
      ! HORIZONTAL
      if(abs(ed(i)).le.pieni) goto 20
!FOX  FOK=EL(I)*ED(I)/DPSQ ;
!FOX  RHO=ONE/ED(I)*DPSQ ;
!FOX  FOK1=SIN(FOK*HALF)/COS(FOK*HALF)/RHO ;
!FOX  SI=SIN(FOK) ;
!FOX  CO=COS(FOK) ;
!FOX  ALDA(IH,1)=ONE ;
!FOX  ALDA(IH,2)=RHO*SI ;
!FOX  ALDA(IH,3)=ZERO ;
!FOX  ALDA(IH,4)=ONE ;
!FOX  ALDA(IH,5)=-DPDA*RHO*(ONE-CO)/DPSQ*C1E3 ;
!FOX  ALDA(IH,6)=-DPDA*TWO*SIN(FOK*HALF)/COS(FOK*HALF)/DPSQ*C1E3 ;
!FOX  SM1=COS(FOK) ;
!FOX  SM2=SIN(FOK)*RHO ;
!FOX  SM3=-SIN(FOK)/RHO ;
!FOX  SM5=-RHO*DPSQ*(ONE-SM1) ;
!FOX  SM6=-SM2*DPSQ/RHO ;
!FOX  SM12=EL(I)-SM1*SM2 ;
!FOX  SM23=SM2*SM3 ;
!FOX  AS3=-RV*(DPDA*RHO/(TWO*DPSQ)*SM23+SM5) ;
!FOX  AS4=-RV*SM23/C2E3 ;
!FOX  AS6=-RV*(EL(I)+SM1*SM2)/C4E3 ;
!FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12+DPDA*(EL(I)-SM2))
!FOX  +EL(I)*(ONE-RV))*C1E3 ;
!FOX  ASDA(IH,2)=-RV*(DPDA/(TWO*RHO*DPSQ)*SM12+SM6)+FOK1*AS3 ;
!FOX  ASDA(IH,3)=AS3 ;
!FOX  ASDA(IH,4)=AS4+TWO*AS6*FOK1 ;
!FOX  ASDA(IH,5)=-RV*SM12/(C4E3*RHO*RHO)+AS6*FOK1*FOK1+FOK1*AS4  ;
!FOX  ASDA(IH,6)=AS6 ;
      ! VERTIKAL
      ih=ih+1
      if(ih.gt.2) ih=1
!FOX  G=SIN(FOK*HALF)/COS(FOK*HALF)/RHO ;
!FOX  GL=EL(I)*G ;
!FOX  ALDA(IH,1)=ONE-GL ;
!FOX  ALDA(IH,2)=EL(I) ;
!FOX  ALDA(IH,3)=-G*(TWO-GL) ;
!FOX  ALDA(IH,4)=ALDA(IH,1) ;
!FOX  AS6=-RV*ALDA(IH,2)/C2E3 ;
!FOX  ASDA(IH,4)=-TWO*AS6*FOK1 ;
!FOX  ASDA(IH,5)=AS6*FOK1*FOK1 ;
!FOX  ASDA(IH,6)=AS6 ;
      goto 190
    case (2)
      goto 100

    !-----------------------------------------------------------------------
    !  SEKTORMAGNET
    !-----------------------------------------------------------------------
    case (3,5)
      if (kz(i) == 3) then
        ih = 1
      else
        ih = 2
      end if
      goto 70

    !-----------------------------------------------------------------------
    !  COMBINED FUNCTION MAGNET
    !-----------------------------------------------------------------------
    case (6,7)
      if (kz(i) == 6) then
        ih = 0
!FOX  FOKQ=EK(I) ;
      else
        ih=1
!FOX  FOKQ=-EK(I) ;
      end if
      if(abs(ek(i)).le.pieni) then
        ih=2
        goto 70
      end if
      if(abs(ed(i)).le.pieni) goto 100
      if(abs(ek(i)-ed(i)**2).le.pieni) goto 20                         !hr08
!FOX  WF=ED(I)/DPSQ ;
!FOX  FOK=FOKQ/DPD-WF*WF ;
!FOX  AFOK=FOK ;
    if(dare(afok).lt.zero) then
!FOX  AFOK=-AFOK ;
    end if
!FOX  HI=SQRT(AFOK) ;
!FOX  FI=HI*EL(I) ;
    if(dare(fok).gt.zero) then
      ! DEFOCUSSING
      ih=ih+1
!FOX  HP=EXP(FI) ;
!FOX  HM=ONE/HP ;
!FOX  HC=(HP+HM)*HALF ;
!FOX  HS=(HP-HM)*HALF ;
!FOX  ALDA(IH,1)=HC ;
!FOX  ALDA(IH,2)=HS/HI ;
!FOX  ALDA(IH,3)=HS*HI ;
!FOX  ALDA(IH,4)=HC ;
!FOX  WFA=WF/AFOK*(ONE-HC)/DPSQ ;
!FOX  WFHI=WF/HI*HS/DPSQ ;
!FOX  ALDA(IH,5)= WFA*DPDA*C1E3 ;
!FOX  ALDA(IH,6)=-WFHI*DPDA*C1E3 ;
!FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;
!FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;
!FOX  ASDA(IH,1)=(RV*(DPDA*DPDA/(FOUR*DPD)*SM12
!FOX  +DPDA*(EL(I)-ALDA(IH,2)))/AFOK*WF*WF+EL(I)*(ONE-RV))*C1E3 ;
!FOX  ASDA(IH,2)=-RV*(DPDA*WF/(TWO*DPSQ)*SM12-DPD*WFHI) ;
!FOX  ASDA(IH,3)=RV*(DPDA*HALF/AFOK/DPD*ED(I)*SM23-DPD*WFA) ;
!FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;
!FOX  ASDA(IH,5)=+RV*SM12*AFOK/C4E3 ;
!FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        ih=ih+1
        if(ih.gt.2) ih=1
!FOX  AEK=EK(I)/DPD ;
        if(dare(aek).lt.zero) then
!FOX  AEK=-AEK ;
        end if
!FOX  HI=SQRT(AEK) ;
!FOX  FI=HI*EL(I) ;
!FOX  SI=SIN(FI) ;
!FOX  CO=COS(FI) ;
!FOX  ALDA(IH,1)=CO ;
!FOX  ALDA(IH,2)=SI/HI ;
!FOX  ALDA(IH,3)=-SI*HI ;
!FOX  ALDA(IH,4)=CO ;
!FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;
!FOX  ASDA(IH,5)=-RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;
!FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
      else
        ih=ih+1
!FOX  SI=SIN(FI) ;
!FOX  CO=COS(FI) ;
!FOX  WFA=WF/AFOK*(ONE-CO)/DPSQ ;
!FOX  WFHI=WF/HI*SI/DPSQ ;
!FOX  ALDA(IH,1)=CO ;
!FOX  ALDA(IH,2)=SI/HI ;
!FOX  ALDA(IH,3)=-SI*HI ;
!FOX  ALDA(IH,4)=CO ;
!FOX  ALDA(IH,5)=-WFA*DPDA*C1E3 ;
!FOX  ALDA(IH,6)=-WFHI*DPDA*C1E3 ;
!FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;
!FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;
!FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12+DPDA
!FOX  *(EL(I)-ALDA(IH,2)))/AFOK*WF*WF+EL(I)*(ONE-RV))*C1E3 ;
!FOX  ASDA(IH,2)=-RV*(DPDA*WF/(TWO*DPSQ)*SM12-DPD*WFHI) ;
!FOX  ASDA(IH,3)=-RV*(DPDA*HALF/AFOK/DPD*ED(I)*SM23-DPD*WFA) ;
!FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;
!FOX  ASDA(IH,5)=-RV*SM12*AFOK/C4E3 ;
!FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        ih=ih+1
        if(ih.gt.2) ih=1
!FOX  AEK=EK(I)/DPD ;
        if(dare(aek).lt.zero) then
!FOX  AEK=-AEK ;
        end if
!FOX  HI=SQRT(AEK) ;
!FOX  FI=HI*EL(I) ;
!FOX  HP=EXP(FI) ;
!FOX  HM=ONE/HP ;
!FOX  HC=(HP+HM)*HALF ;
!FOX  HS=(HP-HM)*HALF ;
!FOX  ALDA(IH,1)=HC ;
!FOX  ALDA(IH,2)=EL(I) ;
!FOX  ALDA(IH,2)=HS/HI ;
!FOX  ALDA(IH,3)=HS*HI ;
!FOX  ALDA(IH,4)=HC ;
!FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;
!FOX  ASDA(IH,5)=+RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;
!FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
      end if
      goto 190

    !-----------------------------------------------------------------------
    !  EDGE FOCUSSING
    !-----------------------------------------------------------------------
    case (8)
!FOX  RHOI=ED(I)/DPSQ ;
!FOX  FOK=RHOI*SIN(EL(I)*RHOI*HALF)/COS(EL(I)*RHOI*HALF) ;
!FOX  ALDA(1,1)=ONE ;
!FOX  ALDA(1,2)=ZERO ;
!FOX  ALDA(1,3)=FOK ;
!FOX  ALDA(1,4)=ONE ;
!FOX  ALDA(2,1)=ONE ;
!FOX  ALDA(2,2)=ZERO ;
!FOX  ALDA(2,3)=-FOK ;
!FOX  ALDA(2,4)=ONE ;
      goto 190
    end select
    !-----------------------------------------------------------------------
    cycle

!-----------------------------------------------------------------------
!  DRIFTLENGTH
!-----------------------------------------------------------------------
   20   do 30 l=1,2
!FOX  ALDA(L,1)=ONE  ;
!FOX  ALDA(L,2)=EL(I) ;
!FOX  ALDA(L,3)=ZERO ;
!FOX  ALDA(L,4)=ONE ;
!FOX  ASDA(L,6)=-RV*ALDA(L,2)/C2E3 ;
   30   continue
!FOX  ASDA(1,1)=EL(I)*(ONE-RV)*C1E3 ;
        goto 190
!-----------------------------------------------------------------------
!  SEKTORMAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
   70   continue
        if(abs(ed(i)).le.pieni) goto 20
!FOX  FOK=EL(I)*ED(I)/DPSQ ;
!FOX  RHO=(ONE/ED(I))*DPSQ ;
!FOX  SI=SIN(FOK) ;
!FOX  CO=COS(FOK) ;
!FOX  RHOC=RHO*(ONE-CO)/DPSQ ;
!FOX  SIQ=SI/DPSQ ;
!FOX  ALDA(IH,1)=CO ;
!FOX  ALDA(IH,2)=RHO*SI ;
!FOX  ALDA(IH,3)=-SI/RHO ;
!FOX  ALDA(IH,4)=CO ;
!FOX  ALDA(IH,5)=-DPDA*RHOC*C1E3 ;
!FOX  ALDA(IH,6)=-DPDA*SIQ*C1E3 ;
!FOX  SM12=EL(I)-ALDA(IH,1)*ALDA(IH,2) ;
!FOX  SM23=ALDA(IH,2)*ALDA(IH,3) ;
!FOX  ASDA(IH,1)=(-RV*(DPDA*DPDA/(FOUR*DPD)*SM12
!FOX  +DPDA*(EL(I)-ALDA(IH,2)))+EL(I)*(ONE-RV))*C1E3 ;
!FOX  ASDA(IH,2)=-RV*(DPDA/(TWO*RHO*DPSQ)*SM12-DPD*SIQ) ;
!FOX  ASDA(IH,3)=-RV*(DPDA*RHO/(TWO*DPSQ)*SM23-DPD*RHOC) ;
!FOX  ASDA(IH,4)=-RV*SM23/C2E3 ;
!FOX  ASDA(IH,5)=-RV*SM12/(C4E3*RHO*RHO) ;
!FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
!--VERTIKAL
        ih=ih+1
        if(ih.gt.2) ih=1
!FOX  ALDA(IH,1)=ONE ;
!FOX  ALDA(IH,2)=EL(I) ;
!FOX  ALDA(IH,3)=ZERO ;
!FOX  ALDA(IH,4)=ONE ;
!FOX  ASDA(IH,6)=-RV*ALDA(IH,2)/C2E3 ;
        goto 190
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
  100   continue
        if(abs(ek(i)).le.pieni) goto 20
!FOX  FOK=EK(I)*MTCDA/(ONE+DPDA) ;
!FOX  AEK=FOK ;
        if(dare(aek).lt.zero) then
!FOX  AEK=-AEK ;
        endif
        ih=0
!FOX  HI=SQRT(AEK) ;
!FOX  FI=EL(I)*HI ;
        if(ek(i).gt.zero) goto 120
  110   ih=ih+1
!FOX  ALDA(IH,1)=COS(FI) ;
!FOX  HI1=SIN(FI) ;
!FOX  ALDA(IH,2)=HI1/HI ;
!FOX  ALDA(IH,3)=-HI1*HI ;
!FOX  ALDA(IH,4)=ALDA(IH,1) ;
!FOX  ASDA(IH,1)=EL(I)*(ONE-RV)*C1E3 ;
!FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;
!FOX  ASDA(IH,5)=-RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;
!FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        if(ih.eq.2) goto 190
!--DEFOCUSSING
  120   ih=ih+1
!FOX  HP=EXP(FI) ;
!FOX  HM=ONE/HP ;
!FOX  HC=(HP+HM)*HALF ;
!FOX  HS=(HP-HM)*HALF ;
!FOX  ALDA(IH,1)=HC ;
!FOX  ALDA(IH,2)=HS/HI ;
!FOX  ALDA(IH,3)=HS*HI ;
!FOX  ALDA(IH,4)=HC ;
!FOX  ASDA(IH,4)=-RV*ALDA(IH,2)*ALDA(IH,3)/C2E3 ;
!FOX  ASDA(IH,5)=+RV*(EL(I)-ALDA(IH,1)*ALDA(IH,2))*AEK/C4E3 ;
!FOX  ASDA(IH,6)=-RV*(EL(I)+ALDA(IH,1)*ALDA(IH,2))/C4E3 ;
        if(ih.eq.1) goto 110
        goto 190
!-----------------------------------------------------------------------
!   NONLINEAR INSERTION
!-----------------------------------------------------------------------
  190 continue
      do ih=1,2
        do ip=1,6
          do ien=1,nord+1
            if (nvar2.eq.5) then
              call dapri6(alda(ih,ip),result,ien,5)
              ald6(i,ih,ip,ien) = result
              call dapri6(asda(ih,ip),result,ien,5)
              asd6(i,ih,ip,ien) = result
            else if (nvar2.eq.6) then
              call dapri6(alda(ih,ip),result,ien,6)
              ald6(i,ih,ip,ien) = result
              call dapri6(asda(ih,ip),result,ien,6)
              asd6(i,ih,ip,ien) = result
            else if (nvar2.eq.4) then
              call dapri6(alda(ih,ip),result,ien,4)
              ald6(i,ih,ip,ien) = result
              call dapri6(asda(ih,ip),result,ien,4)
              asd6(i,ih,ip,ien) = result
            endif
          end do
        end do
      end do
    end do
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION

  return

end subroutine envada

!-----------------------------------------------------------------------
!  CALCULATION OF : MOMENTUM-DEPENDING ELEMENT-MATRICES AND
!                   CHANGE OF PATH LENGTHS FOR EACH PARTICLE.
!      SPECIALLY PREPARED FOR NEW D.A. (SIX-DIMENSIONAL VERSION)
!-----------------------------------------------------------------------
subroutine envquad(i,ipch)
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common
  use mod_common_da
  use mod_lie_dab, only : idao,rscrri,iscrda
  implicit none
  integer i,ih,ipch,idaa
  real(kind=fPrec) dare
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  save
!-----------------------------------------------------------------------
!FOX  B D ;
#include "include/dainicom.f90"
!FOX  D V DA INT FOKQ NORD NVAR ; D V DA INT WFHI NORD NVAR ;
!FOX  D V DA INT DPD NORD NVAR ; D V DA INT DPSQ NORD NVAR ;
!FOX  D V DA INT FOK NORD NVAR ; D V DA INT RHO NORD NVAR ;
!FOX  D V DA INT FOK1 NORD NVAR ; D V DA INT SM1 NORD NVAR ;
!FOX  D V DA INT SM2 NORD NVAR ; D V DA INT SM3 NORD NVAR ;
!FOX  D V DA INT SM4 NORD NVAR ; D V DA INT SM5 NORD NVAR ;
!FOX  D V DA INT SM6 NORD NVAR ; D V DA INT SM12 NORD NVAR ;
!FOX  D V DA INT SM23 NORD NVAR ; D V DA INT AS3 NORD NVAR ;
!FOX  D V DA INT AS4 NORD NVAR ; D V DA INT AS6 NORD NVAR ;
!FOX  D V DA INT SI NORD NVAR ; D V DA INT CO NORD NVAR ;
!FOX  D V DA INT G NORD NVAR ; D V DA INT GL NORD NVAR ;
!FOX  D V DA INT SIQ NORD NVAR ; D V DA INT RHOC NORD NVAR ;
!FOX  D V DA INT HI NORD NVAR ; D V DA INT FI NORD NVAR ;
!FOX  D V DA INT AEK NORD NVAR ; D V DA INT HI1 NORD NVAR ;
!FOX  D V DA INT HP NORD NVAR ; D V DA INT HM NORD NVAR ;
!FOX  D V DA INT HC NORD NVAR ; D V DA INT HS NORD NVAR ;
!FOX  D V DA INT FOKC NORD NVAR ; D V DA INT WF NORD NVAR ;
!FOX  D V DA INT AFOK NORD NVAR ; D V DA INT RHOI NORD NVAR ;
!FOX  D V DA INT WFA NORD NVAR ; D V RE INT RATIOE NELE ;
!FOX  D V RE INT EL NELE ; D V RE INT EK NELE ; D V RE INT ED NELE ;
!FOX  D V RE INT ONE ; D V RE INT ZERO ; D V RE INT TWO ;
!FOX  D V RE INT HALF ; D V RE INT FOUR ; D V RE INT C1E3 ;
!FOX  D V RE INT C2E3 ; D V RE INT C4E3 ;
!FOX  D V IN INT I ; D V IN INT L ; D V IN INT IH ; D V IN INT NE ;
!FOX  D V IN INT NA ; D V IN INT IP ; D V IN INT IPCH ;
!FOX  D F RE DARE 1 ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
!FOX  DPD=ONE+DPDA ;
!FOX  DPSQ=SQRT(DPD) ;
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
  if(abs(ek(i)).le.pieni) goto 100
!FOX  FOK=(SMIDA(IPCH)*RATIOE(I))*MTCDA/(ONE+DPDA) ;
!FOX  AEK=FOK ;
  if(dare(aek).lt.zero) then
!FOX  AEK=-AEK ;
  endif
  ih=0
!FOX  HI=SQRT(AEK) ;
!FOX  FI=EL(I)*HI ;
  if(ek(i).gt.zero) goto 30
20 ih=ih+1
!FOX  ALDAQ(IH,1)=COS(FI) ;
!FOX  HI1=SIN(FI) ;
!FOX  ALDAQ(IH,2)=HI1/HI ;
!FOX  ALDAQ(IH,3)=-HI1*HI ;
!FOX  ALDAQ(IH,4)=ALDAQ(IH,1) ;
!FOX  ASDAQ(IH,1)=EL(I)*(ONE-RV)*C1E3 ;
!FOX  ASDAQ(IH,4)=-RV*ALDAQ(IH,2)*ALDAQ(IH,3)/C2E3 ;
!FOX  ASDAQ(IH,5)=-RV*(EL(I)-ALDAQ(IH,1)*ALDAQ(IH,2))*AEK/C4E3 ;
!FOX  ASDAQ(IH,6)=-RV*(EL(I)+ALDAQ(IH,1)*ALDAQ(IH,2))/C4E3 ;
  if(ih.eq.2) goto 100
!--DEFOCUSSING
30 ih=ih+1
!FOX  HP=EXP(FI) ;
!FOX  HM=ONE/HP ;
!FOX  HC=(HP+HM)*HALF ;
!FOX  HS=(HP-HM)*HALF ;
!FOX  ALDAQ(IH,1)=HC ;
!FOX  ALDAQ(IH,2)=HS/HI ;
!FOX  ALDAQ(IH,3)=HS*HI ;
!FOX  ALDAQ(IH,4)=HC ;
!FOX  ASDAQ(IH,4)=-RV*ALDAQ(IH,2)*ALDAQ(IH,3)/C2E3 ;
!FOX  ASDAQ(IH,5)=+RV*(EL(I)-ALDAQ(IH,1)*ALDAQ(IH,2))*AEK/C4E3 ;
!FOX  ASDAQ(IH,6)=-RV*(EL(I)+ALDAQ(IH,1)*ALDAQ(IH,2))/C4E3 ;
  if(ih.eq.1) goto 20
100 continue
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  return
end subroutine envquad

!-----------------------------------------------------------------------
!  SYNCHROTRON OSCILLATIONS
!        SPECIALLY PREPARED FOR NEW D.A.
!-----------------------------------------------------------------------
subroutine synoda
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_lie_dab, only : idao,iscrri,rscrri,iscrda
  implicit none
  integer ix,idaa,ikz
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  save
!-----------------------------------------------------------------------
!FOX  B D ;
#include "include/dainicom.f90"
!FOX  D V RE INT E0 ; D V RE EXT E0F ;
!FOX  D V RE INT HSY 3 ; D V RE INT PHAS ;
!FOX  D V RE EXT ED NELE ; D V RE EXT HSYC NELE ;
!FOX  D V RE EXT PHASC NELE ;  D V RE INT NUCMDA ;
!FOX  D V RE INT C1E3 ; D V RE INT ONE ; D V IN INT IKZ ;
!FOX  D V IN EXT NELE ; D V IN INT ITION ; D V IN INT IX ;
!FOX  E D ; D V RE INT NUCM0 ; D V RE INT MTCDA ; D V RE INT QQ0 ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
  ix=ixcav

  if(abs(kz(ix)) == 12) then
    ikz = sign(1,kz(ix))
!FOX  EJ1=EJ1+ED(IX)*QQ0*SIN(HSYC(IX)*SIGMDA/C1E3*
!FOX  IKZ+PHASC(IX)) ;
  else
!FOX  EJ1=EJ1+HSY(1)*QQ0*SIN(HSY(3)*SIGMDA/C1E3*ITION+PHAS) ;
  endif
!FOX  EJF1=SQRT(EJ1*EJ1-NUCMDA*NUCMDA) ;
!FOX  DPDA1=(EJF1-E0F)/E0F*C1E3 ;

  return
end subroutine synoda

!----------------------------------------------------------------------*
! PURPOSE:                                                             *
!   MODIFICATION OF WWERF, DOUBLE PRECISION COMPLEX ERROR FUNCTION,    *
!   WRITTEN AT CERN BY K. KOELBIG.                                     *
!   TAKEN FROM MAD8                                                    *
!   VERSION FOR MAP PRODUCTION USING BERZ'S DA PACKAGE                 *
! INPUT:                                                               *
!   XX, YY    (REAL)    ARGUMENT TO CERF.                              *
! OUTPUT:                                                              *
!   WX, WY    (REAL)    FUNCTION RESULT.                               *
!----------------------------------------------------------------------*
subroutine errff(xx,yy,wx,wy)
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_common_track, only : xxtr,yytr,crois,comt_daStart,comt_daEnd
  use mod_common_da
  use mod_lie_dab, only : idao,iscrri,rscrri,iscrda
  implicit none
  integer n,n1,nc,nuu,nuu1,idaa
  real(kind=fPrec) dare,dum
  save
!-----------------------------------------------------------------------
!FOX  B D ;
!FOX  D V DA EXT XX NORD NVAR ; D V DA EXT YY NORD NVAR ;
!FOX  D V DA EXT WX NORD NVAR ; D V DA EXT WY NORD NVAR ;
!FOX  D V DA INT X NORD NVAR ; D V DA INT Y NORD NVAR ;
!FOX  D V DA INT Q NORD NVAR ; D V DA INT H NORD NVAR ;
!FOX  D V DA INT XH NORD NVAR ; D V DA INT YH NORD NVAR ;
!FOX  D V DA INT RX NORD NVAR 33 ; D V DA INT RY NORD NVAR 33 ;
!FOX  D V DA INT TX NORD NVAR ; D V DA INT TN NORD NVAR ;
!FOX  D V DA INT TY NORD NVAR ;D V DA INT SAUX NORD NVAR ;
!FOX  D V DA INT SX NORD NVAR ; D V DA INT SY NORD NVAR ;
!FOX  D V DA INT XL NORD NVAR ;
!FOX  D V RE INT XLIM ; D V RE INT YLIM ; D V RE INT TWO ;
!FOX  D V RE INT ONE ; D V RE INT ZERO ; D V RE INT HALF ;
!FOX  D V RE INT CC ; D V RE INT DUM ;
!FOX  D V IN INT NC ; D V IN INT N ; D V IN INT N1 ; D V IN INT NUU ;
!FOX  D V IN INT NUU1 ; D V IN INT NCC ;
!FOX  D F RE DARE 1 ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
call comt_daStart
!-----------------------------------------------------------------------
!FOX  X=XX ;
!FOX  Y=YY ;
  if(dare(x).lt.zero) then
    write(lout,"(a)") "ERRFF> Problem in DA complex error function: dare(x) < 0"
!FOX    X=-X ;
  endif
  if(dare(y).lt.zero) then
    write(lout,"(a)") "ERRFF> Problem in DA complex error function: dare(y) < 0"
!FOX    Y=-Y ;
  endif
  if(dare(y).lt.ylim.and.dare(x).lt.xlim) then
!FOX    Q=(ONE-Y/YLIM)*SQRT(ONE-X*X/XLIM/XLIM) ;
!FOX    DUM=3.2D0 ;
!FOX    H=ONE/(DUM*Q) ;
    nc=7+int(23.0_fPrec*dare(q))
!FOX    XL=EXP((1-NC)*LOG(H)) ;
!FOX    XH=Y+HALF/H ;
!FOX    YH=X ;
    nuu=10+int(21.0_fPrec*dare(q))
    nuu1=nuu+1
!FOX    RX(NUU1)=ZERO ;
!FOX    RY(NUU1)=ZERO ;
    do 10 n=nuu,1,-1
      n1=n+1
!FOX      TX=XH+N*RX(N1) ;
!FOX      TY=YH-N*RY(N1) ;
!FOX      TN=TX*TX+TY*TY ;
!FOX      RX(N)=HALF*TX/TN ;
!FOX      RY(N)=HALF*TY/TN ;
10   continue
!FOX    SX=ZERO ;
!FOX    SY=ZERO ;
    do 20 n=nc,1,-1
!FOX      SAUX=SX+XL ;
!FOX      SX=RX(N)*SAUX-RY(N)*SY ;
!FOX      SY=RX(N)*SY+RY(N)*SAUX ;
!FOX      XL=H*XL ;
20   continue
!FOX    WX=CC*SX ;
!FOX    WY=CC*SY ;
  else
!FOX    XH=Y ;
!FOX    YH=X ;
!FOX    RX(1)=ZERO ;
!FOX    RY(1)=ZERO ;
    do 30 n=9,1,-1
!FOX      TX=XH+N*RX(1) ;
!FOX      TY=YH-N*RY(1) ;
!FOX      TN=TX*TX+TY*TY ;
!FOX      RX(1)=HALF*TX/TN ;
!FOX      RY(1)=HALF*TY/TN ;
30   continue
!FOX    WX=CC*RX(1) ;
!FOX    WY=CC*RY(1) ;
  endif
!      if(dare(y).eq.0.) then
!!FOX    WX=EXP(-X*X) ;
!      endif
  if(dare(yy).lt.0.d0) then                                          !hr05
!FOX    WX=TWO*EXP(Y*Y-X*X)*COS(TWO*X*Y)-WX ;
!FOX    WY=-TWO*EXP(Y*Y-X*X)*SIN(TWO*X*Y)-WY ;
    if(dare(xx).gt.0.d0) then                                        !hr05
!FOX      WY=-WY ;
    endif
  else
    if(dare(xx).lt.0.d0) then                                        !hr05
!FOX      WY=-WY ;
    endif
  endif
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  call comt_daEnd
  return
end subroutine errff

!-----------------------------------------------------------------------
! WIRE DIFFERENTIAL ALGEBRA
! MODEL OF STRAIGHT CURRENT WIRE
!
!     The model provides a transfer map of a straight current wire.
!     Description:
!     1. Infinitly thin wire with arbitrary orientation.
!     2. Thin element in SixTrack (L)=0
!     3. Parameters:
!     dx, dy: horizontal and vertical distances between wire midpoint
!     and closed orbit [mm]
!     (parameters are given by dx and dy in WIRE block)
!     tx, ty: tilt of the wire w.r.t the closed orbit in the
!     horizontal and vertical planes (in degrees)
!     (parameters are given by tiltx and tilty in WIRE block)
!     L - physical length of the wire element [m]
!     cur - current of the wire [Amperes]
!     embl - embedding drift (integrated length or integration interval) [m]
!     4. The transport map is given for canonical variables (x,px...)
!
! The MAP is constructed out of the following steps:
!     1. Declaration of shifted canonical variables:
!          rx = x+dx; ry = y+dy  in the same way as for the BEAM-BEAM element
!     2. Symplectic Rotation by the tilt angles tx, ty (in 4D space: px, rx, py, ry)
!     3. Wire kick for a longitudinally aligned wire (= kick for tx=ty=0)
!     4. Symplectic Rotation back by the tilt angles -ty, -yx (in 4D space: ...taking only PX, PY)
!--------------------------------------------------------------
!     Normalization factor (in SI) NNORM = (mu0*I*e)/(4*Pi*P0)
!     e -> 1; m0/4Pi -> 1.0e-7; N -> 1.0e-7*I
subroutine wireda(ix,i)

  use floatPrecision
  use mathlib_bouncer
  use physical_constants
  use numerical_constants
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_common_track, only : xxtr,yytr,crois,comt_daStart,comt_daEnd
  use mod_common_da
  use wire
  use mod_lie_dab, only : idao,rscrri,iscrda
  implicit none
  integer ix,idaa,i
  real(kind=fPrec) NNORM_, XCLO, YCLO
  real(kind=fPrec) l,cur,dx,dy,tx,ty,embl,chi
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  save
!-----------------------------------------------------------------------
!FOX  B D ;
#include "include/dainicom.f90"
!FOX  D V DA INT XI NORD NVAR ; D V DA INT YI NORD NVAR ;
!FOX  D V DA INT DXI NORD NVAR ; D V DA INT DYI NORD NVAR ;
!FOX  D V RE INT EMBL ; D V RE INT TX ; D V RE INT TY ;
!FOX  D V RE INT DX ; D V RE INT DY ;
!FOX  D V RE INT XCLO ; D V RE INT YCLO ;
!FOX  D V RE INT CHI ;
!FOX  D V RE INT CUR ;
!FOX  D V RE INT L ; D V RE INT ONE ; D V RE INT TWO ;
!FOX  D V RE INT C1M7 ;
!FOX  D V RE INT C1E3 ; D V RE INT C1M3 ;
!FOX  D V DA INT RTWO_ NORD NVAR ;
!FOX  D V RE INT NNORM_ ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
call comt_daStart
!-----------------------------------------------------------------------
!-- WIRE
!     Normalization factor (in SI) NNORM = (mu0*I*e)/(4*Pi*P0)
!     e -> 1; m0/4Pi -> 1.0e-7; N -> 1.0e-7*I
!     magnetic rigidity

  CHI = (sqrt(e0**2-nucm0**2)*c1e6)/clight
  TX = wire_tiltx(ix) !tilt x [degrees]
  TY = wire_tilty(ix) !tilt y [degrees]
  TX = TX*(pi/c180e0) ![rad]
  TY = TY*(pi/c180e0) ![rad]
  DX = wire_dispx(ix) !displacement x [mm]
  DY = wire_dispy(ix) !displacement y [mm]
  EMBL = wire_lint(ix) !integrated length [m]
  L = wire_lphys(ix) !physical length [m]
  CUR = wire_current(ix)
  XCLO = wire_clo(1,wire_num(i))
  YCLO = wire_clo(2,wire_num(i))
  NNORM_=c1m7/chi

  if (abs(wire_flagco(ix)).ne.1) then
    write(lout,                                                     &
  &fmt='((A,A,/),(A,I0,A,A,/),(A,I0,A,I0,/))')                       &
  &'ERROR: in wirekick -  wire_flagco defined in WIRE block must ',  &
  &'be either 1 or -1!','bez(',ix,') = ',bez(ix),                    &
  &'wire_flagco(',ix,') = ',wire_flagco(ix)
    call prror
  endif

!FOX  YY(1)=YY(1)*C1M3;
!FOX  YY(2)=YY(2)*C1M3;

!!FOX  DXI=DX*C1M3;
!!FOX  DYI=DY*C1M3;
  IF (wire_flagco(ix).eq.1) THEN
!FOX  DXI=(DX+XCLO)*C1M3;
!FOX  DYI=(DY+YCLO)*C1M3;
  ELSE IF (wire_flagco(ix).eq.-1) THEN
!FOX  DXI=DX*C1M3;
!FOX  DYI=DY*C1M3;
  END IF



!-----------------------------------------------------------------------
! X' -> PX'; Y' -> PY
!FOX  YY(1)=YY(1)*(ONE+DPDA)/MTCDA ;
!FOX  YY(2)=YY(2)*(ONE+DPDA)/MTCDA ;

! 1 SHIFT - see the part of the code were wireda is called ....

  IF (wire_flagco(ix).eq.1) THEN
!FOX  XI=(XX(1)+DX)*C1M3;
!FOX  YI=(XX(2)+DY)*C1M3;
  ELSE IF (wire_flagco(ix).eq.-1) THEN
!FOX  XI=(XX(1)+(DX-XCLO))*C1M3;
!FOX  YI=(XX(2)+(DY-YCLO))*C1M3;
  END IF

! ibeco = 0
  if(ibeco.eq.0) then
! 2 symplectic rotation of coordinate system (tx, ty)
!FOX  YI=YI-((XI*SIN(TX))*YY(2))/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(2)*YY(2))/COS(ATAN(YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;
!FOX  XI=XI*(COS(TX)-SIN(TX)*TAN(ATAN(YY(1)/SQRT((ONE+DPDA)*
!FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX)) ;
!FOX  YY(1)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(2)*YY(2))*SIN(ATAN(YY(1)/
!FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;

!FOX  XI=XI-((YI*SIN(TY))*YY(1))/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1))/COS(ATAN(YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;
!FOX  YI=YI*(COS(TY)-SIN(TY)*TAN(ATAN(YY(2)/SQRT((ONE+DPDA)*
!FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY)) ;
!FOX  YY(2)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1))*SIN(ATAN(YY(2)/
!FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;

! 3 apply wire kick
!FOX  RTWO_=XI*XI+YI*YI;
!FOX  YY(1)=YY(1)-(((CUR*NNORM_)*XI)
!FOX  *(SQRT((EMBL+L)*(EMBL+L)+TWO*TWO*RTWO_)
!FOX  -SQRT((EMBL-L)*(EMBL-L)+TWO*TWO*RTWO_)) )/RTWO_;
!FOX  YY(2)=YY(2)-(((CUR*NNORM_)*YI)
!FOX  *(SQRT((EMBL+L)*(EMBL+L)+TWO*TWO*RTWO_)
!FOX  -SQRT((EMBL-L)*(EMBL-L)+TWO*TWO*RTWO_)) )/RTWO_;

! ibeco =1
  elseif(ibeco.eq.1) then

!FOX  DYI=DYI-((DXI*SIN(TX))*YY(2))/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(2)*YY(2))/COS(ATAN(YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;
!FOX  DXI=DXI*(COS(TX)-SIN(TX)*TAN(ATAN(YY(1)/SQRT((ONE+DPDA)*
!FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX)) ;
!FOX  YI=YI-XI*SIN(TX)*YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(2)*YY(2))/COS(ATAN(YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;
!FOX  XI=XI*(COS(TX)-SIN(TX)*TAN(ATAN(YY(1)/SQRT((ONE+DPDA)*
!FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX)) ;
!FOX  YY(1)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(2)*YY(2))*SIN(ATAN(YY(1)/
!FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TX) ;

!FOX  DXI=DXI-((DYI*SIN(TY))*YY(1))/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1))/COS(ATAN(YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;
!FOX  DYI=DYI*(COS(TY)-SIN(TY)*TAN(ATAN(YY(2)/SQRT((ONE+DPDA)*
!FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY)) ;
!FOX  XI=XI-YI*SIN(TY)*YY(1)/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1))/COS(ATAN(YY(2)/SQRT((ONE+DPDA)*(ONE+DPDA)-
!FOX  YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;
!FOX  YI=YI*(COS(TY)-SIN(TY)*TAN(ATAN(YY(2)/SQRT((ONE+DPDA)*
!FOX  (ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY)) ;
!FOX  YY(2)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1))*SIN(ATAN(YY(2)/
!FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))-TY) ;

! 3 apply wire kick
!FOX  RTWO_=XI*XI+YI*YI;
!FOX  YY(1)=YY(1)-(((CUR*NNORM_)*XI)
!FOX  *(SQRT((EMBL+L)*(EMBL+L)+TWO*TWO*RTWO_)
!FOX  -SQRT((EMBL-L)*(EMBL-L)+TWO*TWO*RTWO_)) )/RTWO_;
!FOX  YY(2)=YY(2)-(((CUR*NNORM_)*YI)
!FOX  *(SQRT((EMBL+L)*(EMBL+L)+TWO*TWO*RTWO_)
!FOX  -SQRT((EMBL-L)*(EMBL-L)+TWO*TWO*RTWO_)) )/RTWO_;
! subtract closed orbit kick
! wire kick is negative px -> px - wirekick - (-closed orbit kick)
!FOX  RTWO_=DXI*DXI+DYI*DYI;
!FOX  YY(1)=YY(1)+(((CUR*NNORM_)*DXI)
!FOX  *(SQRT((EMBL+L)*(EMBL+L)+TWO*TWO*RTWO_)
!FOX  -SQRT((EMBL-L)*(EMBL-L)+TWO*TWO*RTWO_)) )/RTWO_;
!FOX  YY(2)=YY(2)+(((CUR*NNORM_)*DYI)
!FOX  *(SQRT((EMBL+L)*(EMBL+L)+TWO*TWO*RTWO_)
!FOX  -SQRT((EMBL-L)*(EMBL-L)+TWO*TWO*RTWO_)) )/RTWO_;
  endif

! 4 symplectic backward rotation of coordinate system (-ty, -tx)
!FOX  YY(2)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1))*SIN(ATAN(YY(2)/
!FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TY) ;
!FOX  YY(1)=SQRT((ONE+DPDA)*(ONE+DPDA)-YY(2)*YY(2))*SIN(ATAN(YY(1)/
!FOX  SQRT((ONE+DPDA)*(ONE+DPDA)-YY(1)*YY(1)-YY(2)*YY(2)))+TX) ;

! PX -> X'; PY -> Y'
!FOX  YY(1)=YY(1)*MTCDA/(ONE+DPDA) ;
!FOX  YY(2)=YY(2)*MTCDA/(ONE+DPDA) ;

!FOX  YY(1)=YY(1)*C1E3;
!FOX  YY(2)=YY(2)*C1E3;

! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  call comt_daEnd
end subroutine wireda

!-----------------------------------------------------------------------
!  CALCULATION OF THE SIX-DIMENSIONAL CLOSED ORBIT
!-----------------------------------------------------------------------
subroutine clorda(nn,idummy,am)
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use matrix_inv
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_common_track
  use mod_common_da
  implicit none
  integer i,i4,icheck,ii,j,j4,k,l,ll,nd2,nn
  real(kind=fPrec) am,cloc,cor,coro,dc,dd,dlo,xx
  integer idummy(nn)
  character(len=6) chp(3),chd(3)
  dimension xx(6),dlo(6),cloc(6),dd(6),dc(6),am(nn,nn)
  integer nerror
  save
!-----------------------------------------------------------------------
  nd2=2*ndimf
  write(lout,10010) nd2
  do l=1,nd2
    xx(l)=zero
    cloc(l)=zero
    dd(l)=zero
    dc(l)=zero
    dlo(l)=zero
    do i=1,nd2
      am(l,i)=zero
    enddo
  enddo
  chp(1)=' CLOX '
  chp(2)=' CLOY '
  chp(3)=' CLOS '
  chd(1)='  D-X '
  chd(2)='  D-Y '
  chd(3)='  D-S '
  cor=zero
  coro=1e38_fPrec
  ii=0
  if(ndimf.eq.3) then
    do l=1,2
      ll=2*l
      cloc(ll-1)=clo6(l)
      cloc(ll)=clop6(l)
    enddo
    cloc(5)=clo6(3)
    cloc(6)=clop6(3)
    if(abs(dppoff).gt.pieni) cloc(6)=dppoff
  else
    do l=1,ndimf
      ll=2*l
      cloc(ll-1)=clo(l)
      cloc(ll)=clop(l)
    enddo
    do l=ndimf+1,3
      ll=2*l
      cloc(ll-1)=zero
      cloc(ll)=zero
    enddo
    cloc(6)=dps(1)
  endif
  do 80 ii=1,itco
    do l=1,2
      ll=2*l
      x(1,l)=cloc(ll-1)
      y(1,l)=cloc(ll)
    enddo
    sigm(1)=cloc(5)
    dps(1)=cloc(6)
    call umlauda
    do i4=1,nd2
      do j4=1,nd2
        am(i4,j4)=aml6(i4,j4)
      enddo
    enddo
    call dinv(nd2,am,nd2,idummy,nerror)
    if(nerror.ne.0) write(lout,*) ' ATTENTION, MATRIX SINGULAR '
    if(ndimf.eq.3) then
      do l=1,2
        ll=2*l
        xx(ll-1)=x(1,l)
        xx(ll)=y(1,l)
      enddo
      xx(5)=sigm(1)
      xx(6)=dps(1)
    else
      do l=1,ndimf
        ll=2*l
        xx(ll-1)=x(1,l)
        xx(ll)=y(1,l)
      enddo
      do l=ndimf+1,3
        ll=2*l
        xx(ll-1)=zero
        xx(ll)=zero
      enddo
    endif
    do l=1,nd2
      dd(l)=cloc(l)-xx(l)
      dc(l)=abs(dd(l))
      if(l.eq.5) dc(5)=dc(5)*c1m2
    enddo
    icheck=0
    do l=1,ndimf
      ll=2*l
      if(dc(ll-1).gt.dma) icheck=1
      if(dc(ll).gt.dmap) icheck=1
    enddo
    if(icheck.eq.0) goto 90
    do k=1,nd2
      dlo(k)=zero
      do j=1,nd2
        dlo(k)=am(k,j)*dd(j)+dlo(k)
      enddo
      if(abs(dppoff).gt.pieni) dlo(6)=zero
    enddo
    write(lout,10020)
    cor=zero
    do l=1,ndimf
      ll=2*l
      write(lout,10060) chp(l),cloc(ll-1),cloc(ll)
      cor=cor+dc(ll-1)**2                                            !hr06
    enddo
    cor=sqrt(cor)
    if(ii.eq.1.or.cor.lt.coro) then
      coro=cor
      do l=1,nd2
        cloc(l)=cloc(l)+dlo(l)
      enddo
      if(ii.ne.itco) then
        write(lout,10030)
        do l=1,ndimf
          ll=2*l
          write(lout,10060) chp(l),cloc(ll-1),cloc(ll)
        enddo
        write(lout,10080) ii,cor
      endif
    else
      write(lout,10040) nd2,ii
      goto 91
    endif
80   continue
  write(lout,10000) itco
  ii=itco
90   continue
  if(ii.ne.itco) then
    do k=1,nd2
      dlo(k)=zero
      do j=1,nd2
        dlo(k)=am(k,j)*dd(j)+dlo(k)
      enddo
      if(abs(dppoff).gt.pieni) dlo(6)=zero
    enddo
    write(lout,10020)
    cor=zero
    do l=1,ndimf
      ll=2*l
      write(lout,10060) chp(l),cloc(ll-1),cloc(ll)
      cor=cor+dc(ll-1)**2                                            !hr06
    enddo
    cor=sqrt(cor)
    if(cor.lt.coro) then
      coro=cor
      do l=1,nd2
        cloc(l)=cloc(l)+dlo(l)
      enddo
      write(lout,10030)
      do l=1,ndimf
        ll=2*l
        write(lout,10060) chp(l),cloc(ll-1),cloc(ll)
      enddo
      write(lout,10080) ii,cor
    else
      write(lout,10040) nd2,ii
      goto 91
    endif
    do l=1,2
      ll=2*l
      x(1,l)=cloc(ll-1)
      y(1,l)=cloc(ll)
    enddo
    sigm(1)=cloc(5)
    dps(1)=cloc(6)
    call umlauda
    do i4=1,nd2
      do j4=1,nd2
        am(i4,j4)=aml6(i4,j4)
      enddo
    enddo
    call dinv(nd2,am,nd2,idummy,nerror)
    if(nerror.ne.0) write(lout,*) ' ATTENTION, MATRIX SINGULAR '
    if(ndimf.eq.3) then
      do l=1,2
        ll=2*l
        xx(ll-1)=x(1,l)
        xx(ll)=y(1,l)
      enddo
      xx(5)=sigm(1)
      xx(6)=dps(1)
    else
      do l=1,ndimf
        ll=2*l
        xx(ll-1)=x(1,l)
        xx(ll)=y(1,l)
      enddo
      do l=ndimf+1,3
        ll=2*l
        xx(ll-1)=zero
        xx(ll)=zero
      enddo
    endif
    do l=1,nd2
      dc(l)=abs(cloc(l)-xx(l))
      if(l.eq.5) dc(5)=dc(5)*c1m2
    enddo
  endif
  write(lout,10050) nd2,ii
  cor=zero
  do l=1,ndimf
    ll=2*l
    write(lout,10070) chp(l),cloc(ll-1),cloc(ll),chd(l),dc(ll-1),dc(ll)
    cor=cor+dc(ll-1)**2                                              !hr06
  enddo
  cor=sqrt(cor)
  write(lout,10080) ii,cor
91   continue
  if(ndimf.eq.3) then
    do l=1,2
      ll=2*l
      clo6(l)=cloc(ll-1)
      clop6(l)=cloc(ll)
    enddo
    clo6(3)=cloc(5)
    clop6(3)=cloc(6)
  else
    do l=1,ndimf
      ll=2*l
      clo(l)=cloc(ll-1)
      clop(l)=cloc(ll)
    enddo
  endif
!-----------------------------------------------------------------------
  return
10000 format(t10,'DA CLOSED ORBIT CALCULATION'/ t10,                    &
  &'MAXIMUM NUMBER OF ITERATIONS ACHIEVED--->',2x,i4/ t10,           &
  &'PROCEDURE MAY NOT HAVE CONVERGED')
10010 format(/131('-')/t10,'ENTERING ',i1,                              &
  &'-D DA CLOSED ORBIT CALCULATION'/)
10020 format(5x,'---- closed orbit before correction----')
10030 format(5x,'---- after DA correction----')
10040 format(/5x,'NO IMPROVEMENT OF ',i1,'-D DA CLOSED ORBIT ',         &
  &'CALCULATION IN ITERATION: ',i4/)
10050 format(t5,'SUCCESSFULL END OF ',i1,'-D DA CLOSED ORBIT ',         &
  &'CALCULATION IN ITERATION: ',i4/)
10060 format(5x,a6,1p,2(1x,g16.9))
10070 format(5x,a6,1p,2(1x,g16.9)/5x,a6,1p,2(1x,g16.9))
10080 format(5x,' ITERAT.=',i3,' ACCURACY=',d13.6/)
end subroutine clorda

!-----------------------------------------------------------------------*
!  FMA                                                                  *
!  M.Fitterer & R. De Maria & K.Sjobak, BE-ABP/HSS                      *
!  last modified: 04-01-2016                                            *
!  purpose: invert the matrix of eigenvecors tas                        *
!           (code copied from postpr only that ta is here fma_tas)      *
!           x(normalized)=fma_tas^-1 x=fma_tas_inv x                    *
!           note: inversion method copied from subroutine postpr        *
!-----------------------------------------------------------------------*
subroutine invert_tas(fma_tas_inv,fma_tas)
  use floatPrecision
  use numerical_constants
  use matrix_inv
  use mod_common_track
  use crcoall
  implicit none

  integer :: i,j            !iterators
  real(kind=fPrec), dimension(6,6), intent(inout) :: fma_tas !tas = normalisation matrix
  real(kind=fPrec), dimension(6,6), intent(out) :: fma_tas_inv !inverse of tas
  integer ierro                   !error messages
!     dummy variables
  real(kind=fPrec), dimension(6,6) :: tdummy !dummy variable for transposing the matrix
  integer, dimension(6) :: idummy !for matrix inversion
!     units: [mm,mrad,mm,mrad,mm,1]
!     invert matrix
!     - set values close to 1 equal to 1
  do i=1,6
      do j=1,6
        fma_tas_inv(i,j)=fma_tas(j,i)
      enddo
  enddo

  if(abs(fma_tas_inv(1,1)).le.pieni.and.abs(fma_tas_inv(2,2)).le.pieni) then
    fma_tas_inv(1,1)=one
    fma_tas_inv(2,2)=one
  endif
  if(abs(fma_tas_inv(3,3)).le.pieni.and.abs(fma_tas_inv(4,4)).le.pieni) then
    fma_tas_inv(3,3)=one
    fma_tas_inv(4,4)=one
  endif
  if(abs(fma_tas_inv(5,5)).le.pieni.and.abs(fma_tas_inv(6,6)).le.pieni) then
    fma_tas_inv(5,5)=one
    fma_tas_inv(6,6)=one
  endif

!     - invert: dinv returns the transposed matrix
  call dinv(6,fma_tas_inv,6,idummy,ierro)
  if (ierro.ne.0) then
      write(lout,*) "Error in INVERT_TAS - Matrix inversion failed!"
      write(lout,*) "Subroutine DINV returned ierro=",ierro
      call prror
  endif

!     - transpose fma_tas_inv
  tdummy=fma_tas_inv
  do i=1,6
    do j=1,6
      fma_tas_inv(i,j)=tdummy(j,i)
    enddo
  enddo
end subroutine invert_tas
