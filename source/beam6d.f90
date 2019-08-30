! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
! ================================================================================================ !
subroutine beamint(np,track,param,sigzs,bcu,ibb,ne,ibtyp,ibbc,mtc)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  use crcoall
  use parpro
  use parbeam, only : beam_expflag,beam_expfile_open
  implicit none

  integer ibb,ibbc,ibtyp,ne,np,nsli
  real(kind=fPrec) alpha,calpha,cphi,f,phi,salpha,sigzs,sphi,tphi,phi2,cphi2,sphi2,tphi2

  real(kind=fPrec) :: track(6,npart) !(6,npart)
  real(kind=fPrec) :: mtc(npart)
  real(kind=fPrec) :: param(nele,18) !(nele,18)
  real(kind=fPrec) :: bcu(nbb,12) !(nbb,12)
  real(kind=fPrec) :: star(3,mbea) !(3,mbea)

  save
!-----------------------------------------------------------------------
  if(beam_expflag .eq. 0) then
    phi=param(ne,1)
    nsli=param(ne,2)
    alpha=param(ne,3)
    f=param(ne,4)/real(nsli,fPrec)
    phi2=param(ne,18)
  else if(beam_expflag .eq. 1) then
    alpha=param(ne,3)
    phi=param(ne,1)
    nsli=param(ne,2)
    !sepax=param(ne,4)     !Not actually used anywhere?
    !sepay=param(ne,5)     !Not actually used anywhere?
    f=param(ne,4)/real(nsli,fPrec)
    phi2=phi               !Note - phi2 is not a free parameter anymore
  else
    write(lerr,"(a,i0,a)") "ERROR beamint: beam_expflag was ",beam_expflag," expected 0 or 1. This is a BUG!"
    call prror
  end if

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
  call boost(np,sphi,cphi,tphi,salpha,calpha,track)
  call sbc(np,star,cphi,cphi2,nsli,f,ibtyp,ibb,bcu,track,ibbc,mtc)
  call boosti(np,sphi,cphi,tphi,salpha,calpha,track)
  return
end subroutine beamint

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
! BOOST Boost Operation ********************************************
!    P,Q,E are all normalized by P0
! ================================================================================================ !
subroutine boost(np,sphi,cphi,tphi,salpha,calpha,track)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  implicit none

  integer i,np
  real(kind=fPrec) calpha,cphi,h,h1x,h1y,h1z,hd1,salpha,sphi,tphi,x1,y1

  real(kind=fPrec) :: track(6,npart) !(6,npart)

  save
!-----------------------------------------------------------------------
  do i=1,np
    h=(track(6,i)+one)-sqrt(((one+track(6,i))**2-track(2,i)**2)-track(4,i)**2)

    track(6,i)=((track(6,i)-(calpha*tphi)*track(2,i))-(track(4,i)*salpha)*tphi)+h*tphi**2
    track(2,i)=(track(2,i)-(tphi*h)*calpha)/cphi
    track(4,i)=(track(4,i)-(tphi*h)*salpha)/cphi

    hd1=sqrt(((one+track(6,i))**2-track(2,i)**2)-track(4,i)**2)

    h1x=track(2,i)/hd1
    h1y=track(4,i)/hd1
    h1z=one-(one+track(6,i))/hd1

    x1=((calpha*tphi)*track(5,i)+(one+(calpha*sphi)*h1x)*track(1,i))+((track(3,i)*salpha)*sphi)*h1x
    y1=((salpha*tphi)*track(5,i)+(one+(salpha*sphi)*h1y)*track(3,i))+((track(1,i)*calpha)*sphi)*h1y

    track(5,i)=track(5,i)/cphi+h1z*((sphi*calpha)*track(1,i)+(sphi*salpha)*track(3,i))
    track(1,i)=x1
    track(3,i)=y1
  end do

  return
end subroutine boost

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
!**SBC ***Synchro-Beam for headon collision**********************
!  call BBF  (table) disabled
!****************************************************************
! ================================================================================================ !
subroutine sbc(np,star,cphi,cphi2,nsli,f,ibtyp,ibb,bcu,track,ibbc,mtc)

  use parpro
  use mod_utils
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer i,ibb,ibbc,ibbc1,ibtyp,jsli,np,nsli
  real(kind=fPrec) bbf0,bbfx,bbfy,bbgx,bbgy,costh,costhp,cphi,f,s,sepx,sepx0,sepy,sepy0,sfac,sinth,sinthp,sp,sx,sy,cphi2

  real(kind=fPrec) :: track(6,npart)
  real(kind=fPrec) :: bcu(nbb,12)
  real(kind=fPrec) :: star(3,mbea)
  real(kind=fPrec) :: mtc(npart)
  real(kind=fPrec) :: dum(13)

  save
!-----------------------------------------------------------------------

  do jsli=1,nsli
    do i=1,np
      s=(track(5,i)-star(3,jsli))*half
      sp=s/cphi2
      dum(1)=(bcu(ibb,1)+(two*bcu(ibb,4))*sp)+bcu(ibb,6)*sp**2
      dum(2)=(bcu(ibb,2)+(two*bcu(ibb,9))*sp)+bcu(ibb,10)*sp**2
      dum(3)=(bcu(ibb,3)+(bcu(ibb,5)+bcu(ibb,7))*sp)+bcu(ibb,8)*sp**2
      dum(4)=dum(1)-dum(2)
      dum(5)=dum(4)**2+four*dum(3)**2

      if(ibbc.eq.1.and.(abs(dum(4)).gt.pieni.and.abs(dum(5)).gt.pieni)) then
        ibbc1=1
        dum(5)=sqrt(dum(5))
      else
        ibbc1=0
      end if

  !JBG New set of canonical set of variables at the Col point (CP)
      sepx0=(track(1,i)+track(2,i)*s)-star(1,jsli)
      sepy0=(track(3,i)+track(4,i)*s)-star(2,jsli)
      if(ibbc1.eq.1) then
        sfac=one

        if(dum(4).lt.zero) then
          sfac=-one*one
        end if

        dum(6)=(sfac*dum(4))/dum(5)
        dum(7)=dum(1)+dum(2)
        costh=half*(one+dum(6))
        if(abs(costh).gt.pieni) then
          costh=sqrt(costh)
        else
          costh=zero
        end if

        sinth=half*(one-dum(6))

        if(abs(sinth).gt.pieni) then
          sinth=(sfac)*sqrt(sinth)
        else
          sinth=zero
        end if

        if(dum(3).lt.zero) then
          sinth=-one*sinth
        end if

        sy=sfac*dum(5)
        sx=(dum(7)+sy)*half
        sy=(dum(7)-sy)*half
        sepx=sepx0*costh+sepy0*sinth
        sepy=sepy0*costh-sepx0*sinth
      else
        sx=dum(1)
        sy=dum(2)
        sepx=sepx0
        sepy=sepy0
      end if

      if(sx.gt.sy) then
        call bbf(sepx,sepy,sx,sy,bbfx,bbfy,bbgx,bbgy,ibtyp)
      else
        call bbf(sepy,sepx,sy,sx,bbfy,bbfx,bbgy,bbgx,ibtyp)
      end if

      bbfx=f*bbfx*mtc(i)
      bbfy=f*bbfy*mtc(i)
      bbgx=f*bbgx*mtc(i)
      bbgy=f*bbgy*mtc(i)

      if(ibbc1.eq.1) then
        dum(8)=two*((bcu(ibb,4)-bcu(ibb,9))+(bcu(ibb,6)-bcu(ibb,10))*sp)
        dum(9)=(bcu(ibb,5)+bcu(ibb,7))+(two*bcu(ibb,8))*sp
        dum(10)=(((dum(4)*dum(8)+(four*dum(3))*dum(9))/dum(5))/dum(5))/dum(5)
        dum(11)=sfac*(dum(8)/dum(5)-dum(4)*dum(10))
        dum(12)=(bcu(ibb,4)+bcu(ibb,9))+(bcu(ibb,6)+bcu(ibb,10))*sp
        dum(13)=(sfac*((dum(4)*dum(8))*half+(two*dum(3))*dum(9)))/dum(5)

        if(abs(costh).gt.pieni) then
          costhp=(dum(11)/four)/costh
        else
          costhp=zero
        end if

        if(abs(sinth).gt.pieni) then
          sinthp=((-one*dum(11))/four)/sinth
        else
          sinthp=zero
        end if

        track(6,i)=track(6,i)-((((bbfx*(costhp*sepx0+sinthp*sepy0)+bbfy*(costhp*sepy0-sinthp*sepx0))&
 &                 +bbgx*(dum(12)+dum(13)))+bbgy*(dum(12)-dum(13)))/cphi)*half
        bbf0=bbfx
        bbfx=bbf0*costh-bbfy*sinth
        bbfy=bbf0*sinth+bbfy*costh
      else
        track(6,i)=track(6,i)-(bbgx*(bcu(ibb,4)+bcu(ibb,6)*sp)+bbgy*(bcu(ibb,9)+bcu(ibb,10)*sp))/cphi
      end if

      track(6,i)=track(6,i)-(bbfx*(track(2,i)-bbfx*half)+bbfy*(track(4,i)-bbfy*half))*half
      track(1,i)=track(1,i)+s*bbfx
      track(2,i)=track(2,i)-bbfx
      track(3,i)=track(3,i)+s*bbfy
      track(4,i)=track(4,i)-bbfy
    end do
  end do

  return

end subroutine sbc

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
! BOOSTI **************inverse boost *****************
! ================================================================================================ !
subroutine boosti(np,sphi,cphi,tphi,salpha,calpha,track)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  implicit none

  integer i,np
  real(kind=fPrec) calpha,cphi,det,h1,h1d,h1x,h1y,h1z,salpha,sphi,tphi,x1,y1,z1

  real(kind=fPrec) :: track(6,npart) !(6,npart)

  save
!-----------------------------------------------------------------------
  do i=1,np
    h1d=sqrt(((one+track(6,i))**2-track(2,i)**2)-track(4,i)**2)
    h1x=track(2,i)/h1d
    h1y=track(4,i)/h1d
    h1z=one-(one+track(6,i))/h1d
    h1=((track(6,i)+one)-sqrt(((one+track(6,i))**2-track(2,i)**2)-track(4,i)**2))*cphi**2

    det=one/cphi+tphi*((h1x*calpha+h1y*salpha)-h1z*sphi)

    x1=(track(1,i)*(one/cphi+(salpha*(h1y-(h1z*salpha)*sphi))*tphi)       &
      +((track(3,i)*salpha)*tphi)*((h1z*calpha)*sphi-h1x))                &
      -(track(5,i)*((calpha+((h1y*calpha)*salpha)*sphi)                   &
      -(h1x*salpha**2)*sphi))*tphi

    y1=(((track(1,i)*calpha)*tphi)*((h1z*salpha)*sphi-h1y)                &
      +track(3,i)*(one/cphi+(calpha*(h1x-(h1z*calpha)*sphi))*tphi))       &
      -(track(5,i)*(salpha-(h1y*calpha**2)*sphi                           &
      +((h1x*calpha)*salpha)*sphi))*tphi

    z1=(track(5,i)*((one+(h1x*calpha)*sphi)+(h1y*salpha)*sphi)            &
      -((track(1,i)*h1z)*calpha)*sphi)-((track(3,i)*h1z)*salpha)*sphi

    track(1,i)=x1/det
    track(3,i)=y1/det
    track(5,i)=z1/det
    track(6,i)=(track(6,i)+(calpha*sphi)*track(2,i))+(salpha*sphi)*track(4,i)
    track(2,i)=(track(2,i)*cphi+(calpha*tphi)*h1)
    track(4,i)=(track(4,i)*cphi+(salpha*tphi)*h1)
  end do

  return

end subroutine boosti

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
! ================================================================================================ !
!**BBF   without using table ******************************************
! gives transverse (f_x and f_y) and longitudinal(g_x and g_y)
! beam-beam kicks except for the kinematical term (nr_e/\gamma)
! SIGXX is \Sigma
! ================================================================================================ !
subroutine bbf(sepx,sepy,sigxx,sigyy,bbfx,bbfy,bbgx,bbgy,ibtyp)

  use parpro
  use mod_utils
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  implicit none

  integer ibtyp
  real(kind=fPrec) arg1x,arg1y,arg2x,arg2y,bbfx,bbfy,bbgx,bbgy,comfac,comfac2,const,expfac,fac,fac2,&
    sepx,sepy,sigxx,sigxy,sigyy,sqrpi2,wx1,wx2,wy1,wy2,x,xxyy

  data sqrpi2/3.544907701811032_fPrec/

  save
!-----------------------------------------------------------------------
  if(sigxx.eq.sigyy) then
    x=sepx**2+sepy**2
    xxyy=sigxx+sigyy
    const=zero

    if(abs(xxyy).gt.pieni) then
      const=x/xxyy
    end if

    expfac=exp_mb(-one*const)
    bbfx=zero
    bbfy=zero
    bbgx=zero
    bbgy=zero

    if(abs(x).gt.pieni) then
      bbfx=((two*sepx)*(one-expfac))/x
      bbfy=((two*sepy)*(one-expfac))/x
      comfac=sepy*bbfy-sepx*bbfx
      comfac2=(abs(sigxx)+abs(sigyy))**2
      bbgx=(comfac+(((four*sepx**2)*const)/x)*expfac)/(two*x)
      bbgy=((((four*sepy**2)*const)/x)*expfac-comfac)/(two*x)
    end if
  else
    x=sepx**2/sigxx+sepy**2/sigyy
    fac2=two*abs(sigxx-sigyy)
    fac=sqrt(fac2)
    const=sqrpi2/fac
    sigxy=sqrt(sigxx/sigyy)
    arg1x=abs(sepx/fac)
    arg1y=abs(sepy/fac)

    if(ibtyp.eq.0) call errf(arg1x,arg1y,wy1,wx1)

    if(ibtyp.eq.1) call wzsub(arg1x,arg1y,wy1,wx1)

    if(x.lt.c1e2) then
      expfac=exp_mb(-half*x)
      arg2x=arg1x/sigxy
      arg2y=arg1y*sigxy

      if(ibtyp.eq.0) call errf(arg2x,arg2y,wy2,wx2)

      if(ibtyp.eq.1) call wzsub(arg2x,arg2y,wy2,wx2)

      bbfx=const*(wx1-expfac*wx2)
      bbfy=const*(wy1-expfac*wy2)

      if(sepx.lt.0) bbfx=-one*bbfx

      if(sepy.lt.0) bbfy=-one*bbfy

      comfac=sepx*bbfx+sepy*bbfy
      bbgx=(-one*(comfac+two*(expfac/sigxy -one)))/fac2
      bbgy= (comfac+two*(expfac*sigxy -one))/fac2
    else
      bbfx=const*wx1
      bbfy=const*wy1

      if(sepx.lt.0) bbfx=-one*bbfx

      if(sepy.lt.0) bbfy=-one*bbfy

      comfac=sepx*bbfx+sepy*bbfy
      bbgx=(-one*(comfac-two))/fac2
      bbgy= -one*bbgx
    end if

  end if

  return

end subroutine bbf

! ================================================================================================ !
!
!   Hirata's 6d beam-beam from BBC
!   SIXTRACK version courtesy Peter Leunissen
!   January 1999
!
!*******STSLD*********************************************************
!   makes longitudinal position of the strong slice for all slices
!*********************************************************************
! ================================================================================================ !
subroutine stsld(star,cphi2,sphi2,sigzs,nsli,calpha,salpha)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  implicit none
  integer i,nsli

  real(kind=fPrec) bord,bord1,border,calpha,cphi,cphi2,gauinv,salpha,sigz,sigzs,sphi,sphi2,yy

  real(kind=fPrec) :: star(3,mbea) !(3,mbea)

!-----------------------------------------------------------------------
  data border /eight/
  save
!-----------------------------------------------------------------------
  sigz=sigzs/cphi2
! DEFINE `STARRED' COORDINATES
!  BORD is longitudinal border star(3,mbea) is the barycenter of region
!  divided two borders.
  bord=+border

  do i=nsli,1,-1
    yy=(one/real(nsli,fPrec))*real(i-1,fPrec)

    if(i.ne.1) bord1=gauinv(yy)

    if(i.eq.1) bord1=-one*border

    star(3,i)=(((exp_mb((-one*bord**2)*half)-exp_mb((-one*bord1**2)*half))/sqrt(two*pi))*real(nsli,fPrec))*sigz
    bord=bord1
    !JBG When doing slicing phi=0 for crab crossing
    ! star(1,i)=0.
    ! star(2,i)=0.
    !JBG When doing slicing phi2 different tiltings of the strong beam
    star(1,i)=(star(3,i)*sphi2)*calpha
    star(2,i)=(star(3,i)*sphi2)*salpha
    !star(1,i)=(star(3,i)*sphi)*calpha
    !star(2,i)=(star(3,i)*sphi)*salpha
  end do

  return

end subroutine stsld

! ================================================================================================ !
!  INVERSE OF (INTEGRATED) NORMAL DISTRIBUTION FUNCTION
!              1         X= Y
!     P(Y)=-----------* INTEGRAL EXP(-X**2/2) DX
!          SQRT(2*PI)    X= -INF
!     IF P(Y)=P0, THEN GAUINV(P0)=Y.
!        0 < P0 < 1 ,   -INF < Y < +INF
!  IF THIS ROUTINE IS USED TO CONVERT UNIFORM RANDOM NUMBERS TO
!  GAUSSIAN, MAXIMUM RELATIVE ERROR IN THE DISTRIBUTION FUNCTION
!  DP/DX=EXP(-X**2/2)/SQRT(2*PI) IS LESS THAN 0.640E-3 EVERYWHERE
!  IN THE RANGE  2**(-31) < P0 < 1-2**31.  (MINIMAX APPROXIMATION)
! ================================================================================================ !
real(kind=fPrec) function gauinv(p0)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use crcoall
  implicit none
  real(kind=fPrec) a0,a1,a2,a3,b0,b1,b2,b3,b4,c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,e0,e1,e2,e3,e4,f0,f1,f2,&
    p,p0,p1,p2,pp1,q,qq2,qq3,qq4,qq5,t
!-----------------------------------------------------------------------
  data pp1/0.334624883253_fPrec/, qq2/0.090230446775_fPrec/,            &
       qq3/0.049905685242_fPrec/, qq4/0.027852994157_fPrec/,            &
       qq5/0.015645650215_fPrec/
  data a3/ 4.5585614e+01_fPrec/, a2/ 2.1635544_fPrec/,                  &
       a1/ 2.7724523_fPrec/,     a0/ 2.5050240_fPrec/,                  &
       b4/ 4.0314354e+02_fPrec/, b3/-2.7713713e+02_fPrec/,              &
       b2/ 7.9731883e+01_fPrec/,                                        &
       b1/-1.4946512e+01_fPrec/, b0/ 2.2157257_fPrec/,                  &
       c4/ 4.1394487e+03_fPrec/, c3/-1.5585873e+03_fPrec/,              &
       c2/ 2.4648581e+02_fPrec/,                                        &
       c1/-2.4719139e+01_fPrec/, c0/ 2.4335936_fPrec/,                  &
       d4/ 4.0895693e+04_fPrec/, d3/-8.5400893e+03_fPrec/,              &
       d2/ 7.4942805e+02_fPrec/,                                        &
       d1/-4.1028898e+01_fPrec/, d0/ 2.6346872_fPrec/,                  &
       e4/ 3.9399134e+05_fPrec/, e3/-4.6004775e+04_fPrec/,              &
       e2/ 2.2566998e+03_fPrec/,                                        &
       e1/-6.8317697e+01_fPrec/, e0/ 2.8224654_fPrec/
  data f0/-8.1807613e-02_fPrec/, f1/-2.8358733_fPrec/,                  &
       f2/ 1.4902469_fPrec/
  save
!-----------------------------------------------------------------------
  p=p0-half
  p1=abs(p)
  gauinv=zero ! -Wmaybe-uninitialized
  if(p1.ge.pp1) goto 120
  p2=p**2
  gauinv=(((a3*p2+a2)*p2+a1)*p2+a0)*p
  return

120  q=half-p1
  if(q.le.qq2) goto 140
  gauinv=(((b4*q+b3)*q+b2)*q+b1)*q+b0
  goto 200

140  if(q.le.qq3) goto 150
  gauinv=(((c4*q+c3)*q+c2)*q+c1)*q+c0
  goto 200

150  if(q.le.qq4) goto 160
  gauinv=(((d4*q+d3)*q+d2)*q+d1)*q+d0
  goto 200

160  if(q.le.qq5) goto 170
  gauinv=(((e4*q+e3)*q+e2)*q+e1)*q+e0
  goto 200

170  if(q.le.zero) goto 900
  t=sqrt(-two*log_mb(q))
  gauinv=(t+f0)+f1/(f2+t)

200  if(p.lt.zero) gauinv=-one*gauinv
  return

900  write(lout,910) p0
910  format(' (FUNC.GAUINV) INVALID INPUT ARGUMENT ',1pd20.13)
  call prror
end function gauinv
