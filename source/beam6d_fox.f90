! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
! ================================================================================================ !
subroutine beaminf(track,param,sigzs,bcu,ibb,ne,ibbc)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use mod_common_da
  use crcoall
  use parpro
  use parbeam, only : beam_expflag,beam_expfile_open
  use mod_lie_dab, only : idao
  implicit none
  integer ibb,ibbc,ne,nsli,idaa
  real(kind=fPrec) alpha,bcu,calpha,cphi,f,param,phi,salpha,sigzs,sphi,star,tphi,phi2,cphi2,sphi2,tphi2
  dimension param(nele,18),bcu(nbb,12),star(3,mbea)
  save
!-----------------------------------------------------------------------
!FOX  B D ;
!FOX  D V DA EXT TRACK NORD NVAR 6 ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
  if (beam_expflag .eq. 0) then
      phi=param(ne,1)
      nsli=param(ne,2)
      alpha=param(ne,3)
      f=param(ne,4)/real(nsli,fPrec)
      phi2=param(ne,18)
  else if(beam_expflag .eq. 1) then
      phi=param(ne,1)
      nsli=param(ne,2)
      alpha=param(ne,3)
      f=param(ne,4)/real(nsli,fPrec)
      !sepax=param(ne,5)     !Not actually used anywhere?
      !sepay=param(ne,6)     !Not actually used anywhere?
      phi2=phi               !Note - phi2 is not a free parameter anymore
  else
      write(lerr,"(a,i0,a)") "ERROR beaminf: beam_expflag was ",beam_expflag," expected 0 or 1. This is a BUG!"
      call prror
  endif
  sphi=sin_mb(phi)
  sphi2=sin_mb(phi2)
  cphi=cos_mb(phi)
  cphi2=cos_mb(phi2)
  tphi=tan_mb(phi)
  tphi2=tan_mb(phi2)
  salpha=sin_mb(alpha)
  calpha=cos_mb(alpha)
!     define slices
  call stsld(star,cphi2,sphi2,sigzs,nsli,calpha,salpha)
  call boostf(sphi,cphi,tphi,salpha,calpha,track)
  call sbcf(star,cphi,cphi2,nsli,f,ibb,bcu,track,ibbc)
  call boostif(sphi,cphi,tphi,salpha,calpha,track)
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  return
end subroutine beaminf

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
! BOOSTF Boost Operation *******************************************
!    P,Q,E are all normalized by P0
! ================================================================================================ !
subroutine boostf(sphi,cphi,tphi,salpha,calpha,track)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common_da
  use mod_lie_dab, only : idao,rscrri,iscrda
  implicit none
  integer idaa
  real(kind=fPrec) calpha,cphi,salpha,sphi,tphi,cphi2,sphi2,tphi2
  save
!-----------------------------------------------------------------------
!FOX  B D ;
!FOX  D V DA EXT TRACK NORD NVAR 6 ; D V DA INT A NORD NVAR ;
!FOX  D V DA INT H NORD NVAR ; D V DA INT SQR1A NORD NVAR ;
!FOX  D V DA INT A1 NORD NVAR ; D V DA INT HD1 NORD NVAR ;
!FOX  D V DA INT H1X NORD NVAR ; D V DA INT H1Y NORD NVAR ;
!FOX  D V DA INT H1Z NORD NVAR ; D V DA INT X1 NORD NVAR ;
!FOX  D V DA INT Y1 NORD NVAR ;
!FOX  D V RE EXT SPHI ; D V RE EXT CPHI ; D V RE EXT TPHI ;
!FOX  D V RE EXT SPHI2 ; D V RE EXT CPHI2 ; D V RE EXT TPHI2 ;
!FOX  D V RE EXT SALPHA ; D V RE EXT CALPHA ;
!FOX  D V RE INT ONE ; D V RE INT C1E3 ;
!FOX  D V DA INT DET NORD NVAR ; D V DA INT H1 NORD NVAR ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
!FOX    H=TRACK(6)+ONE-SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-
!FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)) ;
!FOX    TRACK(6)=TRACK(6)-CALPHA*TPHI*TRACK(2)
!FOX              -TRACK(4)*SALPHA*TPHI+H*TPHI*TPHI ;
!FOX    TRACK(2)=(TRACK(2)-TPHI*H*CALPHA)/CPHI ;
!FOX    TRACK(4)=(TRACK(4)-TPHI*H*SALPHA)/CPHI ;
!FOX    HD1=SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-TRACK(2)*TRACK(2)-
!FOX    TRACK(4)*TRACK(4)) ;
!FOX    H1X=TRACK(2)/HD1 ;
!FOX    H1Y=TRACK(4)/HD1 ;
!FOX    H1Z=ONE-(ONE+TRACK(6))/HD1 ;
!FOX    X1=CALPHA*TPHI*TRACK(5)+(ONE+CALPHA*SPHI*H1X)*TRACK(1)
!FOX       +TRACK(3)*SALPHA*SPHI*H1X ;
!FOX    Y1=SALPHA*TPHI*TRACK(5)+(ONE+SALPHA*SPHI*H1Y)*TRACK(3)
!FOX       +TRACK(1)*CALPHA*SPHI*H1Y ;
!FOX    TRACK(5)=TRACK(5)/CPHI+H1Z*(SPHI*CALPHA*TRACK(1)
!FOX       +SPHI*SALPHA*TRACK(3)) ;
!FOX    TRACK(1)=X1 ;
!FOX    TRACK(3)=Y1 ;
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  return
end subroutine boostf

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
!-----------------------------------------------------------------------
!**SBCF ***Synchro-Beam for headon collision*********************
! ================================================================================================ !
subroutine sbcf(star,cphi,cphi2,nsli,f,ibb,bcu,track,ibbc)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common_da
  use mod_lie_dab, only : idao,rscrri,iscrda
  implicit none
  integer ibb,ibbc,ibbc1,jsli,nsli,idaa
  real(kind=fPrec) bcu,cphi,cphi2,dare,f,sfac,star
  dimension star(3,mbea),bcu(nbb,12)
  save
!-----------------------------------------------------------------------
!FOX  B D ;
!FOX  D V DA EXT TRACK NORD NVAR 6 ;
!FOX  D V DA INT S NORD NVAR ; D V DA INT SP NORD NVAR ;
!FOX  D V DA INT DD NORD NVAR ;
!FOX  D V DA INT DUM NORD NVAR 13 ;
!FOX  D V DA INT SX NORD NVAR ; D V DA INT SY NORD NVAR ;
!FOX  D V DA INT SEPX NORD NVAR ; D V DA INT SEPY NORD NVAR ;
!FOX  D V DA INT SEPX0 NORD NVAR ; D V DA INT SEPY0 NORD NVAR ;
!FOX  D V DA INT COSTH NORD NVAR ; D V DA INT SINTH NORD NVAR ;
!FOX  D V DA INT COSTHP NORD NVAR ; D V DA INT SINTHP NORD NVAR ;
!FOX  D V DA INT BBFX NORD NVAR ; D V DA INT BBFY NORD NVAR ;
!FOX  D V DA INT BBF0 NORD NVAR ; D V DA INT BBGX NORD NVAR ;
!FOX  D V DA INT BBGY NORD NVAR ;
!FOX  D V RE EXT STAR 3 MBEA ; D V RE EXT F ; D V RE EXT BCU NBB 12 ;
!FOX  D V IN EXT IBB ;
!FOX  D V RE INT HALF ; D V RE INT TWO ; D V RE INT FOUR ;
!FOX  D V RE INT ZERO ; D V RE INT ONE ; D V RE INT C1E3 ;
!FOX  D V RE INT SFAC ; D V RE INT CPHI ; D V RE INT CPHI2 ;
!FOX  D V IN INT JSLI ;
!FOX  D F RE DARE 1 ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
  do 2000 jsli=1,nsli
!FOX    S=(TRACK(5)-STAR(3,JSLI))*HALF ;
!FOX    SP=S/CPHI2 ;
!FOX    DUM(1)=BCU(IBB,1)+TWO*BCU(IBB,4)*SP+BCU(IBB,6)*SP*SP ;
!FOX    DUM(2)=BCU(IBB,2)+TWO*BCU(IBB,9)*SP+BCU(IBB,10)*SP*SP ;
!FOX    DUM(3)=BCU(IBB,3)+(BCU(IBB,5)+BCU(IBB,7))*SP+
!FOX    BCU(IBB,8)*SP*SP ;
!FOX    DUM(4)=DUM(1)-DUM(2) ;
!FOX    DUM(5)=DUM(4)*DUM(4)+FOUR*DUM(3)*DUM(3) ;
  if(ibbc.eq.1.and.(abs(dare(dum(4))).gt.pieni.and.               &
  &abs(dare(dum(5))).gt.pieni)) then
      ibbc1=1
!FOX    DUM(5)=SQRT(DUM(5)) ;
  else
      ibbc1=0
  endif
!FOX    SEPX0=TRACK(1)+TRACK(2)*S-STAR(1,JSLI) ;
!FOX    SEPY0=TRACK(3)+TRACK(4)*S-STAR(2,JSLI) ;
  if(ibbc1.eq.1) then
      sfac=one
      if(dare(dum(4)).lt.zero) sfac=(-one*one)
!FOX    DUM(6)=SFAC*DUM(4)/DUM(5) ;
!FOX    DUM(7)=DUM(1)+DUM(2) ;
!FOX    COSTH=HALF*(ONE+DUM(6)) ;
      if(abs(dare(costh)).gt.pieni) then
!FOX    COSTH=SQRT(COSTH) ;
      else
!FOX    COSTH=ZERO ;
      endif
!FOX    SINTH=HALF*(ONE-DUM(6)) ;
      if(abs(dare(sinth)).gt.pieni) then
!FOX    SINTH=SFAC*SQRT(SINTH) ;
      else
!FOX    SINTH=ZERO ;
      endif
      if(dare(dum(3)).lt.zero) then
!FOX    SINTH=-SINTH ;
      endif
!FOX    SY=SFAC*DUM(5) ;
!FOX    SX=(DUM(7)+SY)*HALF ;
!FOX    SY=(DUM(7)-SY)*HALF ;
!FOX    SEPX=SEPX0*COSTH+SEPY0*SINTH ;
!FOX    SEPY=-SEPX0*SINTH+SEPY0*COSTH ;
  else
!FOX    SX=DUM(1) ;
!FOX    SY=DUM(2) ;
!FOX    SEPX=SEPX0 ;
!FOX    SEPY=SEPY0 ;
  endif
  if(dare(sx).gt.dare(sy)) then
      call bbff(sepx,sepy,sx,sy,bbfx,bbfy,bbgx,bbgy)
  else
      call bbff(sepy,sepx,sy,sx,bbfy,bbfx,bbgy,bbgx)
  endif
!FOX    BBFX=F*BBFX ;
!FOX    BBFY=F*BBFY ;
!FOX    BBGX=F*BBGX ;
!FOX    BBGY=F*BBGY ;
  if(ibbc1.eq.1) then
!FOX    DUM(8)=TWO*(BCU(IBB,4)-BCU(IBB,9)+(BCU(IBB,6)-BCU(IBB,10))*SP) ;
!FOX    DUM(9)=BCU(IBB,5)+BCU(IBB,7)+TWO*BCU(IBB,8)*SP ;
!FOX    DUM(10)=(DUM(4)*DUM(8)+FOUR*DUM(3)*DUM(9))/
!FOX    DUM(5)/DUM(5)/DUM(5) ;
!FOX    DUM(11)=SFAC*(DUM(8)/DUM(5)-DUM(4)*DUM(10)) ;
!FOX    DUM(12)=BCU(IBB,4)+BCU(IBB,9)+(BCU(IBB,6)+BCU(IBB,10))*SP ;
!FOX    DUM(13)=SFAC*(DUM(4)*DUM(8)*HALF+TWO*DUM(3)*DUM(9))/DUM(5) ;
      if(abs(dare(costh)).gt.pieni) then
!FOX    COSTHP=DUM(11)/FOUR/COSTH ;
      else
!FOX    COSTHP=ZERO ;
      endif
      if(abs(dare(sinth)).gt.pieni) then
!FOX    SINTHP=-DUM(11)/FOUR/SINTH ;
      else
!FOX    SINTHP=ZERO ;
      endif
!FOX    TRACK(6)=TRACK(6)-
!FOX    (BBFX*(COSTHP*SEPX0+SINTHP*SEPY0)+
!FOX    BBFY*(-SINTHP*SEPX0+COSTHP*SEPY0)+
!FOX    BBGX*(DUM(12)+DUM(13))+BBGY*(DUM(12)-DUM(13)))/
!FOX    CPHI*HALF ;
!FOX    BBF0=BBFX ;
!FOX    BBFX=BBF0*COSTH-BBFY*SINTH ;
!FOX    BBFY=BBF0*SINTH+BBFY*COSTH ;
  else
!FOX    TRACK(6)=TRACK(6)-
!FOX    (BBGX*(BCU(IBB,4)+BCU(IBB,6)*SP)+
!FOX    BBGY*(BCU(IBB,9)+BCU(IBB,10)*SP))/CPHI ;
  endif
!FOX    TRACK(6)=TRACK(6)-(BBFX*(TRACK(2)-BBFX*HALF)+
!FOX    BBFY*(TRACK(4)-BBFY*HALF))*HALF ;
!FOX    TRACK(1)=TRACK(1)+S*BBFX ;
!FOX    TRACK(2)=TRACK(2)-BBFX ;
!FOX    TRACK(3)=TRACK(3)+S*BBFY ;
!FOX    TRACK(4)=TRACK(4)-BBFY ;
2000 continue
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  return
end subroutine sbcf

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
! BOOSTIF **************inverse boost ****************
! ================================================================================================ !
subroutine boostif(sphi,cphi,tphi,salpha,calpha,track)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common_da
  use mod_lie_dab, only : idao,rscrri,iscrda
  implicit none
  integer idaa
  real(kind=fPrec) calpha,cphi,salpha,sphi,tphi
  save
!-----------------------------------------------------------------------
!FOX  B D ;
!FOX  D V DA EXT TRACK NORD NVAR 6 ; D V DA INT A1 NORD NVAR ;
!FOX  D V DA INT H1 NORD NVAR ; D V DA INT SQR1A NORD NVAR ;
!FOX  D V DA INT H1D NORD NVAR ; D V DA INT H1X NORD NVAR ;
!FOX  D V DA INT H1Y NORD NVAR ; D V DA INT H1Z NORD NVAR ;
!FOX  D V DA INT DET NORD NVAR ; D V DA INT X1 NORD NVAR ;
!FOX  D V DA INT Y1 NORD NVAR ; D V DA INT Z1 NORD NVAR ;
!FOX  D V RE INT ONE ; D V RE INT TWO ;
!FOX  D V RE EXT SPHI ; D V RE EXT CPHI ; D V RE EXT TPHI ;
!FOX  D V RE EXT SPHI2 ; D V RE EXT CPHI2 ; D V RE EXT TPHI2 ;
!FOX  D V RE EXT SALPHA ; D V RE EXT CALPHA ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
!FOX    H1D=SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-
!FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)) ;
!FOX    H1X=TRACK(2)/H1D ;
!FOX    H1Y=TRACK(4)/H1D ;
!FOX    H1Z=ONE-(ONE+TRACK(6))/H1D ;
!FOX    H1=(TRACK(6)+ONE-SQRT((ONE+TRACK(6))*(ONE+TRACK(6))-
!FOX    TRACK(2)*TRACK(2)-TRACK(4)*TRACK(4)))*CPHI*CPHI ;
!FOX    DET=ONE/CPHI+TPHI*(H1X*CALPHA+H1Y*SALPHA-H1Z*SPHI) ;
!FOX    X1= TRACK(1)*(ONE/CPHI+SALPHA*(H1Y-H1Z*SALPHA*SPHI)*TPHI)
!FOX       +TRACK(3)*SALPHA*TPHI*(-H1X+H1Z*CALPHA*SPHI)
!FOX       -TRACK(5)*(CALPHA+H1Y*CALPHA*SALPHA*SPHI
!FOX       -H1X*SALPHA*SALPHA*SPHI)*TPHI ;
!FOX    Y1= TRACK(1)*CALPHA*TPHI*(-H1Y+H1Z*SALPHA*SPHI)
!FOX       +TRACK(3)*(ONE/CPHI+CALPHA*(H1X-H1Z*CALPHA*SPHI)*TPHI)
!FOX       -TRACK(5)*(SALPHA-H1Y*CALPHA*CALPHA*SPHI
!FOX       +H1X*CALPHA*SALPHA*SPHI)*TPHI ;
!FOX    Z1=-TRACK(1)*H1Z*CALPHA*SPHI
!FOX       -TRACK(3)*H1Z*SALPHA*SPHI
!FOX       +TRACK(5)*(ONE+H1X*CALPHA*SPHI+H1Y*SALPHA*SPHI) ;
!FOX    TRACK(1)=X1/DET ;
!FOX    TRACK(3)=Y1/DET ;
!FOX    TRACK(5)=Z1/DET ;
!FOX    TRACK(6)=TRACK(6)+CALPHA*SPHI*TRACK(2)
!FOX            +SALPHA*SPHI*TRACK(4) ;
!FOX    TRACK(2)=(TRACK(2)*CPHI+CALPHA*TPHI*H1) ;
!FOX    TRACK(4)=(TRACK(4)*CPHI+SALPHA*TPHI*H1) ;
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  return
end subroutine boostif

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   DA VERSION for SIXTRACK courtesy Peter Leunissen
!   January 1999
!
!**********************************************************************
! BBFF gives transverse (f_x and f_y) and longitudinal(g_x and g_y)
! beam-beam kicks except for the kinematical term (nr_e/\gamma)
! SIGXX is \Sigma
!**********************************************************************
! ================================================================================================ !
subroutine bbff(sepx,sepy,sigxx,sigyy,bbfx,bbfy,bbgx,bbgy)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common_da
  use mod_ranecu
  use mod_lie_dab, only : idao,rscrri,iscrda

  implicit none
  integer idaa
  real(kind=fPrec) dare,hundred,sqrpi2
  parameter(sqrpi2 = 3.544907701811032_fPrec,hundred = c1e2)
  save
!-----------------------------------------------------------------------
!FOX  B D ;
!FOX  D V DA EXT SEPX NORD NVAR ; D V DA EXT SEPY NORD NVAR ;
!FOX  D V DA EXT SIGXX NORD NVAR ; D V DA EXT SIGYY NORD NVAR ;
!FOX  D V DA EXT BBFX NORD NVAR ; D V DA EXT BBFY NORD NVAR ;
!FOX  D V DA EXT BBGX NORD NVAR ; D V DA EXT BBGY NORD NVAR ;
!FOX  D V DA INT SIGXY NORD NVAR ; D V DA INT EXPFAC NORD NVAR ;
!FOX  D V DA INT X NORD NVAR ; D V DA INT FAC NORD NVAR ;
!FOX  D V DA INT FAC2 NORD NVAR ; D V DA INT CONST NORD NVAR ;
!FOX  D V DA INT ARG1X NORD NVAR ; D V DA INT ARG1Y NORD NVAR ;
!FOX  D V DA INT ARG2X NORD NVAR ; D V DA INT ARG2Y NORD NVAR ;
!FOX  D V DA INT WX1 NORD NVAR ; D V DA INT WY1 NORD NVAR ;
!FOX  D V DA INT WX2 NORD NVAR ; D V DA INT WY2 NORD NVAR ;
!FOX  D V DA INT COMFAC NORD NVAR ; D V DA INT COMFAC2 NORD NVAR ;
!FOX  D V RE INT ZERO ; D V RE INT HALF ; D V RE INT ONE ;
!FOX  D V RE INT TWO ; D V RE INT FOUR ; D V RE INT HUNDRED ;
!FOX  D V RE INT SQRPI2 ;
!FOX  D F RE DARE 1 ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------
  if(dare(sigxx).eq.dare(sigyy)) then
!FOX    X=SEPX*SEPX+SEPY*SEPY ;
  if(abs(dare(sigxx)+dare(sigyy)).gt.pieni) then
!FOX      CONST=X/(SIGXX+SIGYY) ;
  else
!FOX      CONST=ZERO ;
  endif
!FOX    EXPFAC=EXP(-CONST) ;
  if(abs(dare(x)).gt.pieni) then
!FOX      BBFX=TWO*SEPX*(ONE-EXPFAC)/X ;
!FOX      BBFY=TWO*SEPY*(ONE-EXPFAC)/X ;
!FOX      COMFAC=-SEPX*BBFX+SEPY*BBFY ;
      if(dare(sigxx).lt.zero) then
!FOX        SIGXX=-SIGXX ;
      endif
      if(dare(sigyy).lt.zero) then
!FOX        SIGYY=-SIGYY ;
      endif
!FOX      COMFAC2=(SIGXX+SIGYY)*(SIGXX+SIGYY) ;
!FOX      BBGX=(COMFAC+FOUR*SEPX*SEPX*CONST/X*EXPFAC)/(TWO*X) ;
!FOX      BBGY=(-COMFAC+FOUR*SEPY*SEPY*CONST/X*EXPFAC)/(TWO*X) ;
  else
!FOX      BBFX=ZERO ;
!FOX      BBFY=ZERO ;
!FOX      BBGX=ZERO ;
!FOX      BBGY=ZERO ;
  endif
  else
!FOX    X=SEPX*SEPX/SIGXX+SEPY*SEPY/SIGYY ;
!FOX    FAC2=TWO*(SIGXX-SIGYY) ;
  if(dare(sigxx).lt.dare(sigyy)) then
!FOX      FAC2=TWO*(SIGYY-SIGXX) ;
  endif
!FOX    FAC=SQRT(FAC2) ;
!FOX    CONST=SQRPI2/FAC ;
!FOX    SIGXY=SQRT(SIGXX/SIGYY) ;
!FOX    ARG1X=(SEPX/FAC) ;
  if(dare(sepx).lt.zero) then
!FOX      ARG1X=-(SEPX/FAC) ;
  endif
!FOX    ARG1Y=(SEPY/FAC) ;
  if(dare(sepy).lt.zero) then
!FOX      ARG1Y=-(SEPY/FAC) ;
  endif
  call errff(arg1x,arg1y,wy1,wx1)
  if(dare(x).lt.hundred) then
!FOX      EXPFAC=EXP(-X*HALF) ;
!FOX      ARG2X=ARG1X/SIGXY ;
!FOX      ARG2Y=ARG1Y*SIGXY ;
      call errff(arg2x,arg2y,wy2,wx2)
!FOX      BBFX=CONST*(WX1-EXPFAC*WX2) ;
!FOX      BBFY=CONST*(WY1-EXPFAC*WY2) ;
      if(dare(sepx).lt.zero) then
!FOX        BBFX=-BBFX ;
      endif
      if(dare(sepy).lt.zero) then
!FOX        BBFY=-BBFY ;
      endif
!FOX      COMFAC=SEPX*BBFX+SEPY*BBFY ;
!FOX      BBGX=-(COMFAC+TWO*(EXPFAC/SIGXY-ONE))/FAC2 ;
!FOX      BBGY= (COMFAC+TWO*(EXPFAC*SIGXY-ONE))/FAC2 ;
  else
!FOX      BBFX=CONST*WX1 ;
!FOX      BBFY=CONST*WY1 ;
      if(dare(sepx).lt.zero) then
!FOX        BBFX=-BBFX ;
      endif
      if(dare(sepy).lt.zero) then
!FOX        BBFY=-BBFY ;
      endif
!FOX      COMFAC=SEPX*BBFX+SEPY*BBFY ;
!FOX      BBGX=-(COMFAC-TWO)/FAC2 ;
!FOX      BBGY= -BBGX ;
  endif
  endif
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION
  return
end subroutine bbff
