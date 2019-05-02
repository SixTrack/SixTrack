module postprocessing

  use floatPrecision
  use parpro

  implicit none

  real(kind=fPrec), private, save :: phase(3,npos+1)
  real(kind=fPrec), private, save :: dani(ninv+1)

contains

subroutine postpr(arg1,arg2)
!-----------------------------------------------------------------------
!  POST PROCESSING
!
!  NFILE   :  FILE UNIT (non-STF) -- always fixed to 90 for STF version.
!  POSI    :  PARTICLE NUMBER
!             (the first particle in pair if ntwin=2, i.e. it is a  pair).
!  NNUML   :  ??
!-----------------------------------------------------------------------
      use mathlib_bouncer
      use numerical_constants
      use matrix_inv
      use crcoall
      use parpro
      use string_tools
      use mod_version
      use mod_time
      use mod_units
      use mod_common_main, only : nnumxv
      use mod_common, only : dpscor,sigcor,icode,idam,its6d,dphix,dphiz,qx0,qz0,&
        dres,dfft,cma1,cma2,nstart,nstop,iskip,iconv,imad,ipos,iav,iwg,ivox,    &
        ivoz,ires,ifh,toptit,kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi,   &
        chromc,tlim,trtime,fort10,fort110,unit10,unit110
#ifdef ROOT
      use root_output
#endif
#ifdef CR
      use checkpoint_restart
#endif
      implicit none

      integer,           intent(in) :: arg1
      integer, optional, intent(in) :: arg2

      integer i,i1,i11,i2,i3,ia,ia0,iaa,iab,iap6,iapx,iapz,ich,idnt,    &
     &ierro,idummy,if1,if2,ife,ife2,ifipa,ifp,ii,ilapa,ilyap,im1,im1s,  &
     &invx,invz,iq,iskc,itopa,iturn,ivo6,iwar6,iwarx,iwarz,j,jm1,jm1s,  &
     &jq,k,k1,nerror,nfft,nfile,nivh,nlost,ntwin,nuex,nuez,nuix,nuiz,   &
     &numl
#ifdef STF
      integer posi,posi1, ia_stf,ifipa_stf,ilapa_stf
#endif
      real tim1,tim2,fxs,fzs
      real(kind=fPrec) const,dle,slope,tle,varlea,wgh
      real(kind=fPrec) alf0,alf04,alf0s2,alf0s3,alf0x2,alf0x3,alf0z2,   &
     &alf0z3,ampx0,ampz0,angi,angii,angiii,ared,ares,armin,armin0,b,b0, &
     &bet0,bet04,bet0s2,bet0s3,bet0x2,bet0x3,bet0z2,bet0z3,biav,bold,c, &
     &c0,c1,c6,clo,cloau,clop,cx,cz,d,d0,d1,dared,dares,di0,di0au,      &
     &di11,dife,dip0,dizu0,dle1,dle1c,dmmac,dnms,dnumlr,dp1,dph6,dphx,  &
     &dphz,dpx,dpxp,dpz,dpzp,dummy,e,e0,e1,emag,emat,emax,emaz,emi,emig,&
     &emii,emiii,emit,emix,emiz,emt,emta,emts,emx,emx0,emxa,emxs,emz,   &
     &emz0,emza,emzs,evt,evt1,evtm,evtma,evtmi,evx,evx1,evx2,evxm,evxma,&
     &evxmi,evz,evz1,evz2,evzm,evzma,evzmi,f,f0,f1,ffx,ffz,finv,g,g0,g1,&
     &gam0s1,gam0s2,gam0s3,gam0x1,gam0x2,gam0x3,gam0z1,gam0z2,gam0z3,h, &
     &h0,h1,p,p1,pcha,pieni2,pinx,pinz,pixr,pizr,pmax,pmin,prec,        &
     &qs0,qwc,ratemx,ratemz,rbeta,s6,sdp6,sdpx,sdpz,sevt,sevx,sevz,     &
     &slopem,sumda,sx,sz,t,ta,ta16,ta26,ta36,ta46,ta56,ta61,ta62,ta63,  &
     &ta64,ta65,tasum,tidnt,tle1,tlo,tph6,tphx,tphz,tpi,txyz,txyz2,x,   &
     &xing,xinv,xp,xp0,xxaux,xxmax,xxmin,xxi,xxr,xyzv,xyzv2,zing,zinv,  &
     &zp,zp0,zzaux,zzmax,zzmin,zzi,zzr

      !The fort.90 file is always with real64, so we need some temps to read it
      ! For the header:
      real(kind=real64) qwc_tmp(3), clo_tmp(3), clop_tmp(3)
      real(kind=real64) di0_tmp(2), dip0_tmp(2)
      real(kind=real64) ta_tmp(6,6)
      real(kind=real64) dmmac_tmp,dnms_tmp,dizu0_tmp,dnumlr_tmp,sigcor_tmp,dpscor_tmp
      real(kind=real64) dummy64

      !For the actual tracking data
      real(kind=real64) b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp
      real(kind=real64) c1_tmp,d1_tmp,e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp

      character(len=80) title(20),chxtit(20),chytit(20)
      character(len=8) cdate,ctime,progrm ! Note: Keep in sync with maincr
      character(len=80) sixtit,commen     ! Note: Keep in sync with mod_common
                                          ! DANGER: If the len changes, CRCHECK will break.

      character(len=11) hvs
      character(len=8192) ch
      character(len=25) ch1
      integer errno,l1,l2,nnuml
      logical rErr
      dimension tle(nlya),dle(nlya)
      dimension wgh(nlya),biav(nlya),slope(nlya),varlea(nlya)
      dimension xinv(ninv),invx(ninv),zinv(ninv),invz(ninv)
      dimension xxr(npos),xxi(npos),zzr(npos),zzi(npos),fxs(npos),fzs(npos)
      dimension bet0(3),alf0(3),t(6,6)
      dimension bet04(2),alf04(2)
      dimension pmin(30),pmax(30)
      dimension idummy(6)
      dimension sumda(60)
      dimension x(2,6),cloau(6),di0au(4)
      dimension qwc(3),clo(3),clop(3),di0(2),dip0(2)
      dimension ta(6,6),txyz(6),txyz2(6),xyzv(6),xyzv2(6),rbeta(6)
      integer itot,ttot
      save

      if(present(arg2)) then
        nnuml = arg2
      else
        nnuml = 0
      end if
#ifdef STF
      posi = arg1
#else
      nfile = arg1
#endif

!----------------------------------------------------------------------
!--TIME START
      pieni2=c1m8
      tlim=c1e7
      call time_timerStart
      tim1=zero
      call time_timerCheck(tim1)

      do i=1,npos
        do j=1,3
          phase(j,i)=zero
        end do
      end do

      do i=1,2
        bet04(i)=zero
        alf04(i)=zero
        di0(i)=zero
        dip0(i)=zero
        di0au(i)=zero
        di0au(i+2)=zero
      end do

      do i=1,3
        bet0(i)=zero
        alf0(i)=zero
        qwc(i)=zero
        clo(i)=zero
        clop(i)=zero
      end do

      do i=1,ninv
        invx(i)=0
        invz(i)=0
        xinv(i)=zero
        zinv(i)=zero
        dani(i)=zero
      end do

      dani(ninv+1)=zero

      do i=1,npos
        xxr(i)=zero
        xxi(i)=zero
        zzr(i)=zero
        zzi(i)=zero
        fxs(i)=0.0
        fzs(i)=0.0
      end do

      do i=1,6
        txyz(i)=zero
        txyz2(i)=zero
        xyzv(i)=zero
        xyzv2(i)=zero
        rbeta(i)=zero
        cloau(i)=zero
        x(1,i)=zero
        x(2,i)=zero
      end do

      do i=1,6
        do j=1,6
          t(i,j)=zero
          ta(i,j)=zero
        end do
      end do

      do i=1,30
        pmax(i)=zero
        pmin(i)=zero
      end do

      do i=1,20
        title(i)=' '
        chxtit(i)=' '
        chytit(i)=' '
      end do

      do i=1,nlya
        tle(i)=zero
        dle(i)=zero
        slope(i)=zero
        varlea(i)=zero
        wgh(i)=zero
        biav(i)=zero
      end do

      do i=1,60
        sumda(i)=zero
      end do

      itot=0
      ttot=0

      do i=1,8
        if (version(i:i).ne.' ') then
          if (version(i:i).ne.'.') then
            itot=itot*10+ichar(version(i:i))-ichar('0')
          else
            ttot=ttot*10**2+itot
            itot=0
          endif
        endif
      enddo

      ttot=ttot*10**2+itot
      sumda(52)=real(ttot,fPrec)                                               !hr06
! and put CPU time for Massimo
! so even if we go to 550 we now get the stats
      sumda(60)=real(trtime,fPrec)
      b0=zero
      nlost=0
      ntwin=1
      nfft=1

      do 120 j=1,npos
        if(nfft.gt.npos/2) goto 130
        nfft=nfft*2
  120 continue
  130 continue
#ifdef STF
      nfile=90
#endif
!----------------------------------------------------------------------
!--READING HEADER
!----------------------------------------------------------------------
      rewind nfile
      ia=0
#ifndef STF
      read(nfile,end=510,iostat=ierro) &
     &     sixtit,commen,cdate,ctime,progrm, &
     &     ifipa,ilapa,itopa,icode,numl, &
     &     qwc_tmp(1),qwc_tmp(2),qwc_tmp(3), &
     &     clo_tmp(1),clop_tmp(1),clo_tmp(2),clop_tmp(2), &
     &     clo_tmp(3),clop_tmp(3), &
     &     di0_tmp(1),dip0_tmp(1),di0_tmp(2),dip0_tmp(2), &
     &     dummy64,dummy64, &
     &     ta_tmp(1,1),ta_tmp(1,2),ta_tmp(1,3), &
     &     ta_tmp(1,4),ta_tmp(1,5),ta_tmp(1,6), &
     &     ta_tmp(2,1),ta_tmp(2,2),ta_tmp(2,3), &
     &     ta_tmp(2,4),ta_tmp(2,5),ta_tmp(2,6), &
     &     ta_tmp(3,1),ta_tmp(3,2),ta_tmp(3,3), &
     &     ta_tmp(3,4),ta_tmp(3,5),ta_tmp(3,6), &
     &     ta_tmp(4,1),ta_tmp(4,2),ta_tmp(4,3), &
     &     ta_tmp(4,4),ta_tmp(4,5),ta_tmp(4,6), &
     &     ta_tmp(5,1),ta_tmp(5,2),ta_tmp(5,3), &
     &     ta_tmp(5,4),ta_tmp(5,5),ta_tmp(5,6), &
     &     ta_tmp(6,1),ta_tmp(6,2),ta_tmp(6,3), &
     &     ta_tmp(6,4),ta_tmp(6,5),ta_tmp(6,6), &
     &     dmmac_tmp,dnms_tmp,dizu0_tmp,dnumlr_tmp, &
     &     sigcor_tmp,dpscor_tmp

      if(ierro.gt.0) then
         write(lout,10320) nfile
#ifdef CR
         goto 551
#else
         goto 550
#endif
      endif

      !Convert it to the current working precission
      do i=1,3
         qwc(i)  = real(qwc_tmp (i), fPrec)
         clo(i)  = real(clo_tmp (i), fPrec)
         clop(i) = real(clop_tmp(i), fPrec)
      end do

      do i=1,2
         di0(i)  = real(di0_tmp (i), fPrec)
         dip0(i) = real(dip0_tmp(i), fPrec)
      enddo

      do i=1,6
         do j=1,6
            ta(j,i) = real(ta_tmp(j,i), fPrec)
         end do
      end do

      dmmac  = real(dmmac_tmp,  fPrec)
      dnms   = real(dnms_tmp,   fPrec)
      dizu0  = real(dizu0_tmp,  fPrec)
      dnumlr = real(dnumlr_tmp, fPrec)
      sigcor = real(sigcor_tmp, fPrec)
      dpscor = real(dpscor_tmp, fPrec)

#ifdef CR
      sumda(1)=nnuml
#else
      sumda(1)=numl
#endif
      idam=1
      if(icode.eq.1.or.icode.eq.2.or.icode.eq.4) idam=1
      if(icode.eq.3.or.icode.eq.5.or.icode.eq.6) idam=2
      if(icode.eq.7) idam=3
      if(ilapa.ne.ifipa) ntwin=2
      if(imad.eq.1.and.progrm.eq.'MAD') then
        imad=0
        rewind nfile
        call join

        read(nfile,end=520,iostat=ierro) &
     &     sixtit,commen,cdate,ctime,progrm, &
     &     ifipa,ilapa,itopa,icode,numl, &
     &     qwc_tmp(1),qwc_tmp(2),qwc_tmp(3), &
     &     clo_tmp(1),clop_tmp(1),clo_tmp(2),clop_tmp(2), &
     &     clo_tmp(3),clop_tmp(3), &
     &     di0_tmp(1),dip0_tmp(1),di0_tmp(2),dip0_tmp(2), &
     &     dummy64,dummy64, &
     &     ta_tmp(1,1),ta_tmp(1,2),ta_tmp(1,3), &
     &     ta_tmp(1,4),ta_tmp(1,5),ta_tmp(1,6), &
     &     ta_tmp(2,1),ta_tmp(2,2),ta_tmp(2,3), &
     &     ta_tmp(2,4),ta_tmp(2,5),ta_tmp(2,6), &
     &     ta_tmp(3,1),ta_tmp(3,2),ta_tmp(3,3), &
     &     ta_tmp(3,4),ta_tmp(3,5),ta_tmp(3,6), &
     &     ta_tmp(4,1),ta_tmp(4,2),ta_tmp(4,3), &
     &     ta_tmp(4,4),ta_tmp(4,5),ta_tmp(4,6), &
     &     ta_tmp(5,1),ta_tmp(5,2),ta_tmp(5,3), &
     &     ta_tmp(5,4),ta_tmp(5,5),ta_tmp(5,6), &
     &     ta_tmp(6,1),ta_tmp(6,2),ta_tmp(6,3), &
     &     ta_tmp(6,4),ta_tmp(6,5),ta_tmp(6,6), &
     &     dmmac_tmp,dnms_tmp,dizu0_tmp,dnumlr_tmp, &
     &     sigcor_tmp,dpscor_tmp

        if(ierro.gt.0) then
           write(lout,10320) nfile
#ifdef CR
           goto 551
#else
           goto 550
#endif
        endif

        !Convert it to the current working precission
        do i=1,3
           qwc(i)  = real(qwc_tmp (i), fPrec)
           clo(i)  = real(clo_tmp (i), fPrec)
           clop(i) = real(clop_tmp(i), fPrec)
        end do

        do i=1,2
           di0(i)  = real(di0_tmp (i), fPrec)
           dip0(i) = real(dip0_tmp(i), fPrec)
        enddo

        do i=1,6
           do j=1,6
              ta(j,i) = real(ta_tmp(j,i), fPrec)
           end do
        end do

        dmmac  = real(dmmac_tmp,  fPrec)
        dnms   = real(dnms_tmp,   fPrec)
        dizu0  = real(dizu0_tmp,  fPrec)
        dnumlr = real(dnumlr_tmp, fPrec)
        sigcor = real(sigcor_tmp, fPrec)
        dpscor = real(dpscor_tmp, fPrec)

        !MadX convention
        ta(1,6)=ta(1,6)*c1e3
        ta(2,6)=ta(2,6)*c1e3
        ta(3,6)=ta(3,6)*c1e3
        ta(4,6)=ta(4,6)*c1e3
        ta(5,6)=ta(5,6)*c1e3
        ta(6,1)=ta(6,1)*c1m3
        ta(6,2)=ta(6,2)*c1m3
        ta(6,3)=ta(6,3)*c1m3
        ta(6,4)=ta(6,4)*c1m3
        ta(6,5)=ta(6,5)*c1m3

      endif
#else
      !Read header lines until a match is found
 555  continue
      read(nfile,end=510,iostat=ierro) &
     &     sixtit,commen,cdate,ctime,progrm, &
     &     ifipa,ilapa,itopa,icode,numl, &
     &     qwc_tmp(1),qwc_tmp(2),qwc_tmp(3), &
     &     clo_tmp(1),clop_tmp(1),clo_tmp(2),clop_tmp(2), &
     &     clo_tmp(3),clop_tmp(3), &
     &     di0_tmp(1),dip0_tmp(1),di0_tmp(2),dip0_tmp(2), &
     &     dummy64,dummy64, &
     &     ta_tmp(1,1),ta_tmp(1,2),ta_tmp(1,3), &
     &     ta_tmp(1,4),ta_tmp(1,5),ta_tmp(1,6), &
     &     ta_tmp(2,1),ta_tmp(2,2),ta_tmp(2,3), &
     &     ta_tmp(2,4),ta_tmp(2,5),ta_tmp(2,6), &
     &     ta_tmp(3,1),ta_tmp(3,2),ta_tmp(3,3), &
     &     ta_tmp(3,4),ta_tmp(3,5),ta_tmp(3,6), &
     &     ta_tmp(4,1),ta_tmp(4,2),ta_tmp(4,3), &
     &     ta_tmp(4,4),ta_tmp(4,5),ta_tmp(4,6), &
     &     ta_tmp(5,1),ta_tmp(5,2),ta_tmp(5,3), &
     &     ta_tmp(5,4),ta_tmp(5,5),ta_tmp(5,6), &
     &     ta_tmp(6,1),ta_tmp(6,2),ta_tmp(6,3), &
     &     ta_tmp(6,4),ta_tmp(6,5),ta_tmp(6,6), &
     &     dmmac_tmp,dnms_tmp,dizu0_tmp,dnumlr_tmp, &
     &     sigcor_tmp,dpscor_tmp

      if(ifipa.ne.posi) then    !IFIPA=first particle, POSI=requested particle
        goto 555                !Get the next header...
      endif
      ! TODO: Protect against no valid headers found,
      ! i.e. posi > itopa.
      if(ierro.gt.0) then
         write(lout,10320) nfile
#ifdef CR
         goto 551
#else
         goto 550
#endif
      endif

      !Convert it to the current working precission
      do i=1,3
         qwc(i)  = real(qwc_tmp (i), fPrec)
         clo(i)  = real(clo_tmp (i), fPrec)
         clop(i) = real(clop_tmp(i), fPrec)
      end do

      do i=1,2
         di0(i)  = real(di0_tmp (i), fPrec)
         dip0(i) = real(dip0_tmp(i), fPrec)
      end do

      do i=1,6
         do j=1,6
            ta(j,i) = real(ta_tmp(j,i), fPrec)
         end do
      end do

      dmmac  = real(dmmac_tmp,  fPrec)
      dnms   = real(dnms_tmp,   fPrec)
      dizu0  = real(dizu0_tmp,  fPrec)
      dnumlr = real(dnumlr_tmp, fPrec)
      sigcor = real(sigcor_tmp, fPrec)
      dpscor = real(dpscor_tmp, fPrec)

#ifdef CR
      sumda(1)=nnuml
#else
      sumda(1)=numl
#endif
      idam=1
      if(icode.eq.1.or.icode.eq.2.or.icode.eq.4) idam=1
      if(icode.eq.3.or.icode.eq.5.or.icode.eq.6) idam=2
      if(icode.eq.7) idam=3
      if(ilapa.ne.ifipa) then !Is first particle != Last particle?
        ntwin=2               !(ntwin=1 is the default in postpr)
!--   binrecs is indexed as 1,2,3,... (=i.e.(91-nfile) in the non-STF version,
!--   while posi values are called as 1,3,5, so using posi1 for crbinrecs index later
      endif
      posi1 = (posi+1)/2 !For both ntwin=1 and 2
#endif

#ifndef STF
!--PREVENT FAULTY POST-PROCESSING
      read(nfile,end=530,iostat=ierro) iaa

      if(ierro.gt.0) then
        write(lout,10320) nfile
        goto 550
      endif

      read(nfile,end=535,iostat=ierro) iab

      if(ierro.gt.0) then
        write(lout,10320) nfile
        goto 550
      endif
#else
!--PREVENT FAULTY POST-PROCESSING
      !--bypass headers
      rewind nfile

      do i=1,itopa,2
         read(nfile)
      enddo

      !--read first track data for particle at posi
      do !The loop is safe, it will anyway end if EOF is reached.
         read(nfile,end=530,iostat=ierro) iaa, j
         if (j.eq.posi) exit
      enddo

      if(ierro.gt.0) then
         write(lout,10320) nfile
         goto 550
      endif

      !--bypass records until the 2nd turn of same particle is reached
      do
         read(nfile,end=535,iostat=ierro) iab, j
         if (j.eq.posi) exit
      enddo

      if(ierro.gt.0) then
        write(lout,10320) nfile
        goto 550
      endif
#endif

      if((((numl+1)/iskip)/(iab-iaa))/iav.gt.nlya) nstop=iav*nlya !hr06

      rewind nfile

!-- Bypassing header to read tracking data later
#ifndef STF
      read(nfile)
#else
      do i=1,itopa,2
         read(nfile) !One header per particle pair.
      enddo
#endif
      sumda(5)=ta(1,1)**2+ta(1,2)**2                                     !hr06
      sumda(6)=ta(3,3)**2+ta(3,4)**2                                     !hr06

      if(iconv.eq.1) then
        cma1=one
        cma2=one
        clo(1)=zero
        clo(2)=zero
        clo(3)=zero
        clop(1)=zero
        clop(2)=zero
        clop(3)=zero
        di0(1)=zero
        di0(2)=zero
        dip0(1)=zero
        dip0(2)=zero

        do i=1,6
          do j=1,6
            if(i.ne.j) then
              ta(i,j)=zero
            else
              ta(i,j)=one
            end if
          end do
       end do
      endif

      cloau(1)= clo(1)
      cloau(2)=clop(1)
      cloau(3)= clo(2)
      cloau(4)=clop(2)
      cloau(5)= clo(3)
      cloau(6)=clop(3)
      di0au(1)= di0(1)
      di0au(2)=dip0(1)
      di0au(3)= di0(2)
      di0au(4)=dip0(2)
      sigcor=cma2
      dpscor=cma1
      sumda(53)= clo(1)
      sumda(54)=clop(1)
      sumda(55)= clo(2)
      sumda(56)=clop(2)
      sumda(57)= clo(3)
      sumda(58)=clop(3)
      if(ifipa.eq.ilapa.and.ndafi.gt.itopa) ndafi=itopa
      if(ilapa-ifipa.eq.1.and.ndafi.gt.itopa/2) ndafi=itopa/2

!-----------------------------------------------------------------------
!  OPTICAL PARAMETERS AT THE STARTING POINT
!-----------------------------------------------------------------------
      ta16=ta(1,6)*c1m3
      ta26=ta(2,6)*c1m3
      ta36=ta(3,6)*c1m3
      ta46=ta(4,6)*c1m3
      ta56=ta(5,6)*c1m3
      ta61=ta(6,1)*c1e3
      ta62=ta(6,2)*c1e3
      ta63=ta(6,3)*c1e3
      ta64=ta(6,4)*c1e3
      ta65=ta(6,5)*c1e3
      bet0(1)=ta(1,1)**2+ta(1,2)**2                                      !hr06
      bet0x2 =ta(1,3)**2+ta(1,4)**2                                      !hr06
      bet0x3 =ta(1,5)**2+ta16**2                                         !hr06
      gam0x1 =ta(2,1)**2+ta(2,2)**2                                      !hr06
      gam0x2 =ta(2,3)**2+ta(2,4)**2                                      !hr06
      gam0x3 =ta(2,5)**2+ta26**2                                         !hr06
      alf0(1)=-one*(ta(1,1)*ta(2,1)+ta(1,2)*ta(2,2))                     !hr06
      alf0x2 =-one*(ta(1,3)*ta(2,3)+ta(1,4)*ta(2,4))                     !hr06
      alf0x3 =-one*(ta(1,5)*ta(2,5)+ta16*ta26)                           !hr06
      bet0(2)=ta(3,3)**2+ta(3,4)**2                                      !hr06
      bet0z2 =ta(3,1)**2+ta(3,2)**2                                      !hr06
      bet0z3 =ta(3,5)**2+ta36**2                                         !hr06
      gam0z1 =ta(4,3)**2+ta(4,4)**2                                      !hr06
      gam0z2 =ta(4,1)**2+ta(4,2)**2                                      !hr06
      gam0z3 =ta(4,5)**2+ta46**2                                         !hr06
      alf0(2)=-one*(ta(3,3)*ta(4,3)+ta(3,4)*ta(4,4))                     !hr06
      alf0z2 =-one*(ta(3,1)*ta(4,1)+ta(3,2)*ta(4,2))                     !hr06
      alf0z3 =-one*(ta(3,5)*ta(4,5)+ta36*ta46)                           !hr06
      bet0(3)=ta(5,5)**2+ta56**2                                         !hr06
      bet0s2 =ta(5,1)**2+ta(5,2)**2                                      !hr06
      bet0s3 =ta(5,3)**2+ta(5,4)**2                                      !hr06
      gam0s1 =ta65**2+ta(6,6)**2                                         !hr06
      gam0s2 =ta61**2+ta62**2                                            !hr06
      gam0s3 =ta63**2+ta64**2                                            !hr06
      alf0(3)=-one*(ta(5,5)*ta65+ta56*ta(6,6))                           !hr06
      alf0s2 =-one*(ta(5,1)*ta61+ta(5,2)*ta62)                           !hr06
      alf0s3 =-one*(ta(5,3)*ta63+ta(5,4)*ta64)                           !hr06
      bet04(1)=bet0(1)
      bet04(2)=bet0(2)
      alf04(1)=alf0(1)
      alf04(2)=alf0(2)

      if(bet0(1).le.pieni.or.bet0(2).le.pieni) then
        write(lout,*) 'WARNING: BETA VALUES ARE ZERO'
        bet0(1)=zero
        bet0(2)=zero
      endif

      do i=1,3
        ii=2*i
        rbeta(ii-1)=sqrt(bet0(i))
        rbeta(ii)=rbeta(ii-1)
        if(abs(rbeta(ii-1)).lt.pieni) rbeta(ii-1)=one
        if(abs(rbeta(ii)).lt.pieni) rbeta(ii)=one
      end do

!----------------------------------------------------------------------
!--SETTING UP OF THE PARAMETERS
!----------------------------------------------------------------------
!--HPLOT TITLES
      if(icode.eq.1) hvs(1:11)='Hor        '
      if(icode.eq.2) hvs(1:11)='    Ver    '
      if(icode.eq.3) hvs(1:11)='Hor Ver    '
      if(icode.eq.4) hvs(1:11)='        Syn'
      if(icode.eq.5) hvs(1:11)='Hor     Syn'
      if(icode.eq.6) hvs(1:11)='    Ver Syn'
      if(icode.eq.7) hvs(1:11)='Hor Ver Syn'
      toptit(2)(1:13)='Particle no. '
      write(toptit(2)(14:16),'(I3)') ifipa !WARNING: Does not work for > 999 particles
      toptit(2)(17:30)=', Phase Space '
      write(toptit(2)(31:41),'(A11)') hvs
      toptit(2)(42:50)=', Dp/p = '
      toptit(2)(62:80)=' '
      toptit(3)(1:5)='Ax = '
      toptit(3)(16:22)=', Ay = '
      toptit(3)(33:80)=' '
      toptit(4)(1:5)='Qx = '
      write(toptit(4)(6:15),10010) qwc(1)
      toptit(4)(16:22)=', Qy = '
      write(toptit(4)(23:32),10010) qwc(2)
      toptit(4)(33:39)=', Qs = '
      write(toptit(4)(39:48),10010) qwc(3)
      toptit(4)(49:80)=' '
      title(1)='Normalized Distance of Phase Space D as a Function of Turn Number N'
      title(2)='Normalized Horizontal Phase Space Projection'
      title(3)='Normalized Vertical Phase Space Projection'
      title(4)='Physical Phase Space Projection'
      title(5)='Synchrotron Phase Space Projection'
      title(6)='Synchrotron versus Horizontal Phase Space Projection'
      title(7)='Synchrotron versus Vertical Phase Space Projection'
      title(8)='Energy E as a Function of Turn Number N'
      title(9)='Stroboscoped Normalized Horizontal Phase Space Projection'
      title(10)='Stroboscoped Normalized Vertical Phase Space Projection'
      title(11)='FFT Analysis of the X Coordinate'
      title(12)='FFT Analysis of the Y Coordinate'
      chxtit(1)='N'
      chytit(1)='D (PI rad)'
      chxtit(2)='X (mm)'
      chytit(2)='Normalized Slope of X (mm)'
      chxtit(3)='Y (mm)'
      chytit(3)='Normalized Slope of Y (MM)'
      chxtit(4)='X (mm)'
      chytit(4)='Y (mm)'
      chxtit(5)='Sigma (mm)'
      chytit(5)='Relative Momentum Deviation'
      chxtit(6)='X (mm)'
      chytit(6)='Relative Momentum Deviation'
      chxtit(7)='Y (mm)'
      chytit(7)='Relative Momentum Deviation'
      chxtit(8)='N'
      chytit(8)='E (MeV)'
      chxtit(9)='X (mm)'
      chytit(9)='Normalized Slope of X (mm)'
      chxtit(10)='Y (mm)'
      chytit(10)='Normalized Slope of Y (mm)'
      chxtit(11)='Horizontal Tune Qx'
      chytit(11)='Normalized FFT Signal'
      chxtit(12)='Vertical Tune Qy'
      chytit(12)='Normalized FFT Signal'
      if(idis.ne.0.or.icow.ne.0.or.istw.ne.0.or.iffw.ne.0) then
        call hplsiz(15.,15.,' ')
        call hplset('VSIZ',.24)
        call hplset('ASIZ',.19)
        call hplset('XLAB',1.5)
        call hplset('YLAB',1.0)
        call hplset('GSIZ',.19)
      endif
      if(iav.lt.1) iav=1
      if(nprint.eq.1) then
        write(lout,10040) sixtit,commen
        write(lout,10050) progrm,ifipa,itopa,hvs,numl,bet0(1),bet0x2,bet0x3,bet0(2),bet0z2,bet0z3,bet0(3),bet0s2,bet0s3,alf0(1), &
     & alf0x2,alf0x3
        write(lout,10060) alf0(2),alf0z2,alf0z3,alf0(3),alf0s2,alf0s3,gam0x1,gam0x2,gam0x3,gam0z1,gam0z2,gam0z3,gam0s1,gam0s2,   &
     & gam0s3,clo(1),clo(2),clo(3),clop(1),clop(2),clop(3),di0(1),di0(2),dip0(1),dip0(2),qwc(1),qwc(2),qwc(3)
        write(lout,10070) iav,nstart,nstop,dphix,dphiz,iwg, qx0,qz0
        write(lout,10080) ivox,ivoz,ires,dres,ifh,dfft
        write(lout,10090) idis,icow,istw,iffw
        write(lout,10100) iskip,iconv,imad,cma1,cma2,nprint,ndafi
      endif ! END if(nprint.eq.1)

!--INITIALISATION
      tpi=eight*atan_mb(one)                                               !hr06
      prec=c1m1
      i1=0
      i11=1
      tlo=zero
      i2=0
      ifp=0
      iwarx=0
      iwarz=0
      iwar6=0
      iapx=0
      iapz=0
      iap6=0
      ivo6=1
      qs0=zero
      armin0=c1e9                                                         !hr06
      armin=armin0
      nivh=ninv/2
      finv=tpi/real(ninv,fPrec)                                                !hr06
      dani(1)=zero
      dani(ninv+1)=tpi

      do i=1,ninv-1
        dani(i+1)=real(i,fPrec)*finv                                             !hr06
      end do

      dle1=zero
      bold=zero
      dle1c=zero
      const=zero
      dphx=zero
      dphz=zero
      dph6=zero
      tphx=zero
      tphz=zero
      tph6=zero
      sdpx=zero
      sdpz=zero
      sdp6=zero
      evx=zero
      evz=zero
      evx2=zero
      evz2=zero
      evt=zero
      sevx=zero
      sevz=zero
      sevt=zero
      emax=zero
      emix=zero
      emaz=zero
      emiz=zero
      emag=zero
      emig=zero
      emxa=zero
      emza=zero
      emta=zero
      emxs=zero
      emzs=zero
      emts=zero
      nuex=0
      nuez=0
      nuix=0
      nuiz=0
      xing=zero
      zing=zero
      pinx=zero
      pinz=zero
      pixr=zero
      pizr=zero
!--INVERTING THE MATRIX OF THE GENERATING VECTORS
!     ta = matrix of eigenvectors already normalized, rotated and ordered, units mm,mrad,mm,mrad,mm,1
!     t  = inverse(ta), units mm,mrad,mm,mrad,mm,1
!
!     This is similar but not exactly the same as the subroutine invert_tas;
!     the "tasum" is missing from there
      do i=1,6
        do j=1,6
          t(i,j)=ta(j,i)
        end do
      end do

      if(abs(t(1,1)).le.pieni.and.abs(t(2,2)).le.pieni) then
        t(1,1)=one
        t(2,2)=one
      endif

      if(abs(t(3,3)).le.pieni.and.abs(t(4,4)).le.pieni) then
        t(3,3)=one
        t(4,4)=one
      endif

      if(abs(t(5,5)).le.pieni.and.abs(t(6,6)).le.pieni) then
        t(5,5)=one
        t(6,6)=one
      endif

      tasum=zero
      its6d=0

      do i=1,6
        tasum=(tasum+abs(t(i,5)))+abs(t(i,6))                            !hr06
      end do

      do i=1,4
        tasum=(tasum+abs(t(5,i)))+abs(t(6,i))                            !hr06
      end do

      tasum=tasum-two

      if(abs(tasum).ge.pieni) its6d=1

      call dinv(6,t,6,idummy,nerror)

      if(nerror.eq.-1) then  !TODO: Using the file number makes no sense in STF case (seen in multiple places)
        write(lout,10290) nfile
        goto 550
      endif

!----------------------------------------------------------------------
!--FIND MINIMUM VALUE OF THE DISTANCE IN PHASESPACE
!----------------------------------------------------------------------
  190 ifipa=0
#ifndef STF
      if(ntwin.eq.1) then
         read(nfile,end=200,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp

         b=real(b_tmp,fPrec)
         c=real(c_tmp,fPrec)
         d=real(d_tmp,fPrec)
         e=real(e_tmp,fPrec)
         f=real(f_tmp,fPrec)
         g=real(g_tmp,fPrec)
         h=real(h_tmp,fPrec)
         p=real(p_tmp,fPrec)

      elseif(ntwin.eq.2) then
         read(nfile,end=200,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp,ilapa,b_tmp,c1_tmp,d1_tmp,    &
     &e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp

         b=real(b_tmp,fPrec)
         c=real(c_tmp,fPrec)
         d=real(d_tmp,fPrec)
         e=real(e_tmp,fPrec)
         f=real(f_tmp,fPrec)
         g=real(g_tmp,fPrec)
         h=real(h_tmp,fPrec)
         p=real(p_tmp,fPrec)

         c1=real(c1_tmp,fPrec)
         d1=real(d1_tmp,fPrec)
         e1=real(e1_tmp,fPrec)
         f1=real(f1_tmp,fPrec)
         g1=real(g1_tmp,fPrec)
         h1=real(h1_tmp,fPrec)
         p1=real(p1_tmp,fPrec)

      endif
#else
!STF case: read tracking data until one reaches right particle.
      if(ntwin.eq.1) then
         read(nfile,end=200,iostat=ierro) ia_stf,ifipa_stf,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp

      elseif(ntwin.eq.2) then
         read(nfile,end=200,iostat=ierro) ia_stf,ifipa_stf,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp,    &
     &        ilapa_stf,b_tmp,c1_tmp,d1_tmp,e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp

      endif

      if(ifipa_stf.ne.posi) then
        goto 190
      endif

      ! Found the right particle; load data in memory (otherwise it's corrupted when EOF is reached)
      ia=ia_stf
      ifipa=ifipa_stf

      b=real(b_tmp,fPrec)
      c=real(c_tmp,fPrec)
      d=real(d_tmp,fPrec)
      e=real(e_tmp,fPrec)
      f=real(f_tmp,fPrec)
      g=real(g_tmp,fPrec)
      h=real(h_tmp,fPrec)
      p=real(p_tmp,fPrec)

      if(ntwin.eq.2) then
         ilapa=ilapa_stf

         c1=real(c1_tmp,fPrec)
         d1=real(d1_tmp,fPrec)
         e1=real(e1_tmp,fPrec)
         f1=real(f1_tmp,fPrec)
         g1=real(g1_tmp,fPrec)
         h1=real(h1_tmp,fPrec)
         p1=real(p1_tmp,fPrec)
      endif
#endif
      if(ierro.gt.0) then
        write(lout,10320) nfile
        goto 550
      endif
      if(ifipa.lt.1) goto 190
      if((ia-nstart).lt.0) goto 190

      if(progrm.eq.'MAD') then
#ifndef STF
        c=c*c1e3
        d=d*c1e3
        e=e*c1e3
        f=f*c1e3
        h=h*c1e3
        p=p*c1e3
        if(ntwin.eq.2) then
          c1=c1*c1e3
          d1=d1*c1e3
          e1=e1*c1e3
          f1=f1*c1e3
          h1=h1*c1e3
          p1=p1*c1e3
        endif
#else
        write(lout,*) "ERROR in postpr: program=MAD not valid for STF."
        call prror(-1)
#endif
      endif ! END if(program.eq.'MAD')

      if(ntwin.eq.2) then
        x(1,1)=c
        x(1,2)=d
        x(1,3)=e
        x(1,4)=f
        x(1,5)=g
        x(1,6)=h
        x(2,1)=c1
        x(2,2)=d1
        x(2,3)=e1
        x(2,4)=f1
        x(2,5)=g1
        x(2,6)=h1
        call distance(x,cloau,di0au,t,b)
      endif
      if(nstop.gt.nstart.and.(ia-nstop).gt.0) goto 200                   !hr06
      if(b.lt.b0.or.abs(b0).le.pieni) b0=b
      goto 190
  200 if(ia.le.0) goto 530

      rewind nfile

!----------------------------------------------------------------------
!--GET FIRST DATA POINT AS A REFERENCE
!----------------------------------------------------------------------
! Skip the header(s)
#ifndef STF
      read(nfile,iostat=ierro)
#else
      do i=1,itopa,2
         read(nfile,iostat=ierro)
      enddo
#endif
      if(ierro.gt.0) then
        write(lout,10320) nfile
#ifdef CR
        goto 551
#else
        goto 550
#endif
      endif

#ifdef CR
!--   Initiate count of binary records
#ifndef STF
      crbinrecs(91-nfile)=1
#else
      crbinrecs(posi1)=1
#endif
#endif

 210  ifipa=0
#ifndef STF
      if(ntwin.eq.1) then
         read(nfile,end=530,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp
         b=real(b_tmp,fPrec)
         c=real(c_tmp,fPrec)
         d=real(d_tmp,fPrec)
         e=real(e_tmp,fPrec)
         f=real(f_tmp,fPrec)
         g=real(g_tmp,fPrec)
         h=real(h_tmp,fPrec)
         p=real(p_tmp,fPrec)

      elseif(ntwin.eq.2) then
         read(nfile,end=530,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp,ilapa,b_tmp,c1_tmp,d1_tmp,    &
     &e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp

         b=real(b_tmp,fPrec)
         c=real(c_tmp,fPrec)
         d=real(d_tmp,fPrec)
         e=real(e_tmp,fPrec)
         f=real(f_tmp,fPrec)
         g=real(g_tmp,fPrec)
         h=real(h_tmp,fPrec)
         p=real(p_tmp,fPrec)

         c1=real(c1_tmp,fPrec)
         d1=real(d1_tmp,fPrec)
         e1=real(e1_tmp,fPrec)
         f1=real(f1_tmp,fPrec)
         g1=real(g1_tmp,fPrec)
         h1=real(h1_tmp,fPrec)
         p1=real(p1_tmp,fPrec)
      endif
#else
!     STF case: read tracking data until one reaches right particle.
      if(ntwin.eq.1) then
         read(nfile,end=530,iostat=ierro) ia_stf,ifipa_stf,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp

      elseif(ntwin.eq.2) then
         read(nfile,end=530,iostat=ierro) ia_stf,ifipa_stf,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp, &
     &        ilapa_stf,b_tmp,c1_tmp,d1_tmp,e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp
      endif

      if(ifipa_stf.ne.posi) then
         goto 210
      endif

      ! Found the right particle; load data in memory (otherwise it's corrupted when EOF is reached)
      ia=ia_stf
      ifipa=ifipa_stf

      b=real(b_tmp,fPrec)
      c=real(c_tmp,fPrec)
      d=real(d_tmp,fPrec)
      e=real(e_tmp,fPrec)
      f=real(f_tmp,fPrec)
      g=real(g_tmp,fPrec)
      h=real(h_tmp,fPrec)
      p=real(p_tmp,fPrec)

      if(ntwin.eq.2) then
         ilapa=ilapa_stf

         c1=real(c1_tmp,fPrec)
         d1=real(d1_tmp,fPrec)
         e1=real(e1_tmp,fPrec)
         f1=real(f1_tmp,fPrec)
         g1=real(g1_tmp,fPrec)
         h1=real(h1_tmp,fPrec)
         p1=real(p1_tmp,fPrec)
      endif
#endif
      if(ierro.gt.0) then
        write(lout,10320) nfile
        goto 550
      endif
#ifdef CR
!     Count one more binary record
#ifndef STF
      crbinrecs(91-nfile)=crbinrecs(91-nfile)+1
#else
      crbinrecs(posi1)=crbinrecs(posi1)+1
#endif
#endif
      if(ifipa.lt.1) goto 210
      if((ia-nstart).lt.0) goto 210
      if(progrm.eq.'MAD') then
#ifndef STF
        c=c*c1e3
        d=d*c1e3
        e=e*c1e3
        f=f*c1e3
        g=g*c1e3
        p=p*c1e3
        if(ntwin.eq.2) then
          c1=c1*c1e3
          d1=d1*c1e3
          e1=e1*c1e3
          f1=f1*c1e3
          g1=g1*c1e3
          p1=p1*c1e3
        endif
#else
        write(lout,*) "ERROR in postpr: program=MAD not valid for STF."
        call prror(-1)
#endif
      endif

      dp1=h
      write(toptit(2)(51:61),10000) dp1-clop(3)
      if(nprint.eq.1.and.ia.eq.0) then
        write(lout,*) 'INITIAL COORDINATES'
        write(lout,*) '       X = ',c
        write(lout,*) '      XP = ',d
        write(lout,*) '       Z = ',e
        write(lout,*) '      ZP = ',f
        write(lout,*) '   SIGMA = ',g
        write(lout,*) '    DP/P = ',h
        write(lout,*) '  ENERGY = ',p
      endif

      if(nstop.gt.nstart.and.(ia-nstop).gt.0) goto 540
      ia=ia-nstart
!--LYAPUNOV
      if(ntwin.eq.2) then
!     first particle
        x(1,1)=c
        x(1,2)=d
        x(1,3)=e
        x(1,4)=f
        x(1,5)=g
        x(1,6)=h
!     twin particle
        x(2,1)=c1
        x(2,2)=d1
        x(2,3)=e1
        x(2,4)=f1
        x(2,5)=g1
        x(2,6)=h1
        call distance(x,cloau,di0au,t,b)
      endif
!--KEEP THE FIRST TURN NUMBER : IA0
      ia0=ia
      xxr(1)=c
      xxi(1)=zero
      zzr(1)=e
      zzi(1)=zero
      c=c-clo(1)
      d=d-clop(1)
      e=e-clo(2)
      f=f-clop(2)
      g=g-clo(3)
      h=h-clop(3)
      c1=c1-clo(1)
      d1=d1-clop(1)
      e1=e1-clo(2)
      f1=f1-clop(2)
      g1=g1-clo(3)
      h1=h1-clop(3)
      if(icode.ge.4) then
        c=c-di0(1)*h
        d=d-dip0(1)*h
        e=e-di0(2)*h
        f=f-dip0(2)*h
        c1=c1-di0(1)*h
        d1=d1-dip0(1)*h
        e1=e1-di0(2)*h
        f1=f1-dip0(2)*h
      endif
!     calculation first particle
!--EMITTANCES
      xp0=bet0(1)*d+alf0(1)*c
      zp0=bet0(2)*f+alf0(2)*e
      emx=(c**2+xp0**2)/bet0(1)                                          !hr06
      emz=(e**2+zp0**2)/bet0(2)                                          !hr06
      if(icode.ge.4.and.its6d.ne.0) then
        c=c+di0(1)*h
        d=d+dip0(1)*h
        e=e+di0(2)*h
        f=f+dip0(2)*h
        c1=c1+di0(1)*h
        d1=d1+dip0(1)*h
        e1=e1+di0(2)*h
        f1=f1+dip0(2)*h
      endif
      emt=emx+emz
      emax=emx
      emix=emx
      emxa=emx
      emaz=emz
      emiz=emz
      emza=emz
      emat=emt
      emit=emt
      emta=emt
      emx0=emx
      emz0=emz

!--COURANT SYNDER
      xyzv(1)=c
      xyzv(2)=d
      xyzv(3)=e
      xyzv(4)=f
      xyzv(5)=g
      xyzv(6)=h

!--CONVERT TO CANONICAL VARIABLES
      if(its6d.eq.1) then
        xyzv(2)=xyzv(2)*((one+xyzv(6))+clop(3))                          !hr06
        xyzv(4)=xyzv(4)*((one+xyzv(6))+clop(3))                          !hr06
      endif

      do iq=1,6
        txyz(iq)=zero
        do jq=1,6
          txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
        end do
      end do

!--INITIAL COORDINATES
      if(nprint.eq.1.and.ia.eq.0) then
        write(lout,*) 'DISTANCE = ',b
      endif

!--EMITTANCES WITH LINEAR COUPLING
      evx=txyz(1)**2+txyz(2)**2                                          !hr06
      evz=txyz(3)**2+txyz(4)**2                                          !hr06
!     calculation second particle
      xyzv2(1)=c1
      xyzv2(2)=d1
      xyzv2(3)=e1
      xyzv2(4)=f1
      xyzv2(5)=g1
      xyzv2(6)=h1

!--CONVERT TO CANONICAL VARIABLES
      if(its6d.eq.1) then
        xyzv2(2)=xyzv2(2)*((one+xyzv2(6))+clop(3))                       !hr06
        xyzv2(4)=xyzv2(4)*((one+xyzv2(6))+clop(3))                       !hr06
      endif

      do iq=1,6
        txyz2(iq)=zero
        do jq=1,6
          txyz2(iq)=txyz2(iq)+t(jq,iq)*xyzv2(jq)
        end do
      end do

      evx2=txyz2(1)**2+txyz2(2)**2                                       !hr06
      evz2=txyz2(3)**2+txyz2(4)**2                                       !hr06
      write(toptit(3)(6:15),10010) sqrt(evx*bet0(1))+sqrt(evz*bet0x2)
      write(toptit(3)(23:32),10010) sqrt(evz*bet0(2))+sqrt(evx*bet0z2)

      if(its6d.eq.1) then
        emiii=txyz(5)**2*cma2**2+txyz(6)**2*cma1**2
      else
        emiii=zero
      endif

!--COURANT SYNDER CONT.
      do iq=1,6
        txyz(iq)=txyz(iq)*rbeta(iq)
      end do

      c0=txyz(1)
      d0=txyz(2)
      e0=txyz(3)
      f0=txyz(4)
      g0=txyz(5)*cma2
      h0=txyz(6)*cma1

!--MIN MAX VALUES
      pmin(2)=b
      pmax(2)=b
      pmin(3)=c0
      pmax(3)=c0
      pmin(4)=d0
      pmax(4)=d0
      pmin(5)=e0
      pmax(5)=e0
      pmin(6)=f0
      pmax(6)=f0
      pmin(9)=g0
      pmax(9)=g0
      pmin(10)=h0
      pmax(10)=h0
      pmin(16)=p
      pmax(16)=p

!--EMITTANCES WITH LINEAR COUPLING CONT.
      emi=evx
      emii=evz
      angi=zero
      angii=zero
      angiii=zero
      if(abs(txyz(1)).gt.pieni.or.abs(txyz(2)).gt.pieni)                &
     &angi=atan2_mb(txyz(2),txyz(1))
      if(abs(txyz(3)).gt.pieni.or.abs(txyz(4)).gt.pieni)                &
     &angii=atan2_mb(txyz(4),txyz(3))
      if(abs(txyz(5)).gt.pieni.or.abs(txyz(6)).gt.pieni)                &
     &angiii=atan2_mb(txyz(6)*cma1,txyz(5)*cma2)
      evt=evx+evz
      evxma=evx
      evzma=evz
      evtma=evt
      evxmi=evx
      evzmi=evz
      evtmi=evt
!--COORDINATE-ANGLE CONVERSION
      call caconv(dpx,d0,c0)
      call caconv(dpz,f0,e0)
      dpxp=tpi+dpx
      dpzp=tpi+dpz
!--INVARIANTS
      call cinvar(dpx,dphix,dpz,dpzp,nuex,emz,zinv,invz)
      call cinvar(dpz,dphiz,dpx,dpxp,nuez,emx,xinv,invx)
!----------------------------------------------------------------------
!--GET DATA POINTS
!----------------------------------------------------------------------
      iskc=0
  240 ifipa=0
#ifndef STF
      if(ntwin.eq.1) then
         read(nfile,end=270,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp

         b=real(b_tmp,fPrec)
         c=real(c_tmp,fPrec)
         d=real(d_tmp,fPrec)
         e=real(e_tmp,fPrec)
         f=real(f_tmp,fPrec)
         g=real(g_tmp,fPrec)
         h=real(h_tmp,fPrec)
         p=real(p_tmp,fPrec)

      elseif(ntwin.eq.2) then
         read(nfile,end=270,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp,ilapa,b_tmp,c1_tmp,d1_tmp,    &
     &e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp

         b=real(b_tmp,fPrec)
         c=real(c_tmp,fPrec)
         d=real(d_tmp,fPrec)
         e=real(e_tmp,fPrec)
         f=real(f_tmp,fPrec)
         g=real(g_tmp,fPrec)
         h=real(h_tmp,fPrec)
         p=real(p_tmp,fPrec)

         c1=real(c1_tmp,fPrec)
         d1=real(d1_tmp,fPrec)
         e1=real(e1_tmp,fPrec)
         f1=real(f1_tmp,fPrec)
         g1=real(g1_tmp,fPrec)
         h1=real(h1_tmp,fPrec)
         p1=real(p1_tmp,fPrec)
      endif
#else
!     STF case: read tracking data until one reaches right particle.
      if(ntwin.eq.1) then
         read(nfile,end=270,iostat=ierro) ia_stf,ifipa_stf,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp

      elseif(ntwin.eq.2) then
         read(nfile,end=270,iostat=ierro) ia_stf,ifipa_stf, &
     &        b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp,ilapa_stf,b_tmp,c1_tmp,d1_tmp,e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp
      endif

      if(ifipa_stf.ne.posi) then
         goto 240
      endif

      ! Found the right particle; load data in memory (otherwise it's corrupted when EOF is reached)
      ia=ia_stf
      ifipa=ifipa_stf

      b=real(b_tmp,fPrec)
      c=real(c_tmp,fPrec)
      d=real(d_tmp,fPrec)
      e=real(e_tmp,fPrec)
      f=real(f_tmp,fPrec)
      g=real(g_tmp,fPrec)
      h=real(h_tmp,fPrec)
      p=real(p_tmp,fPrec)

      if(ntwin.eq.2) then
         ilapa=ilapa_stf

         c1=real(c1_tmp,fPrec)
         d1=real(d1_tmp,fPrec)
         e1=real(e1_tmp,fPrec)
         f1=real(f1_tmp,fPrec)
         g1=real(g1_tmp,fPrec)
         h1=real(h1_tmp,fPrec)
         p1=real(p1_tmp,fPrec)
      endif
#endif
      if(ierro.gt.0) then
        write(lout,10320) nfile
        goto 550
      endif
#ifdef CR
#ifndef STF
!--Increment crbinrecs by 1
      crbinrecs(91-nfile)=crbinrecs(91-nfile)+1
#else
      crbinrecs(posi1)=crbinrecs(posi1)+1
#endif
#endif
      if(ifipa.lt.1) goto 240
      if(progrm.eq.'MAD') then
#ifndef STF
        c=c*c1e3
        d=d*c1e3
        e=e*c1e3
        f=f*c1e3
        g=g*c1e3
        p=p*c1e3
        if(ntwin.eq.2) then
          c1=c1*c1e3
          d1=d1*c1e3
          e1=e1*c1e3
          f1=f1*c1e3
          g1=g1*c1e3
          p1=p1*c1e3
        endif
#else
        write(lout,*) "ERROR in postpr: program=MAD not valid for STF."
        call prror(-1)
#endif
      endif
!--LYAPUNOV
      if(ntwin.eq.2) then
        x(1,1)=c
        x(1,2)=d
        x(1,3)=e
        x(1,4)=f
        x(1,5)=g
        x(1,6)=h
        x(2,1)=c1
        x(2,2)=d1
        x(2,3)=e1
        x(2,4)=f1
        x(2,5)=g1
        x(2,6)=h1
        call distance(x,cloau,di0au,t,b)
      endif
      iskc=iskc+1
      if(mod(iskc,iskip).ne.0) goto 240
      if(nstop.gt.nstart.and.(ia-nstop).gt.0) goto 270
      i1=i1+1
      i11=i1+1
      if(i2.ge.nlya.and.i11.gt.nfft.and.iapx.gt.npos.and.iapz.gt.npos)  &
     &goto 270
      if(i11.le.nfft) then
        xxr(i11)=c
        xxi(i11)=zero
      endif
      if(i11.le.nfft) then
        zzr(i11)=e
        zzi(i11)=zero
      endif
      c=c-clo(1)
      d=d-clop(1)
      e=e-clo(2)
      f=f-clop(2)
      g=g-clo(3)
      h=h-clop(3)
      if(icode.ge.4) then
        c=c-di0(1)*h
        d=d-dip0(1)*h
        e=e-di0(2)*h
        f=f-dip0(2)*h
      endif
!--EMITTANCES
      xp=bet0(1)*d+alf0(1)*c
      zp=bet0(2)*f+alf0(2)*e
      emx=(c**2+xp**2)/bet0(1)                                           !hr06
      emz=(e**2+zp**2)/bet0(2)                                           !hr06
      if(icode.ge.4.and.its6d.ne.0) then
        c=c+di0(1)*h
        d=d+dip0(1)*h
        e=e+di0(2)*h
        f=f+dip0(2)*h
      endif
      emt=emx+emz
      emxa=emxa+emx
      emza=emza+emz
      emta=emta+emt
      emax=max(emx,emax)
      emix=min(emx,emix)
      emaz=max(emz,emaz)
      emiz=min(emz,emiz)
      emat=max(emt,emat)
      emit=min(emt,emit)

!--COURANT SYNDER
      xyzv(1)=c
      xyzv(2)=d
      xyzv(3)=e
      xyzv(4)=f
      xyzv(5)=g
      xyzv(6)=h

!--CONVERT TO CANONICAL VARIABLES
      if(its6d.eq.1) then
        xyzv(2)=xyzv(2)*((one+xyzv(6))+clop(3))                          !hr06
        xyzv(4)=xyzv(4)*((one+xyzv(6))+clop(3))                          !hr06
      endif

      do iq=1,6
        txyz(iq)=zero
        do jq=1,6
          txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
        end do
      end do

!--EMITTANCES WITH LINEAR COUPLING
      evx1=txyz(1)**2+txyz(2)**2                                         !hr06
      evz1=txyz(3)**2+txyz(4)**2                                         !hr06

!--COURANT SYNDER CONT.
      do 260 iq=1,6
        txyz(iq)=txyz(iq)*rbeta(iq)
  260 continue
      c=txyz(1)
      d=txyz(2)
      e=txyz(3)
      f=txyz(4)
      g=txyz(5)*cma2
      h=txyz(6)*cma1

!--MIN MAX VALUES
      pmin(2)= min(pmin(2) ,b)
      pmax(2)= max(pmax(2) ,b)
      pmin(3)= min(pmin(3) ,c)
      pmax(3)= max(pmax(3) ,c)
      pmin(4)= min(pmin(4) ,d)
      pmax(4)= max(pmax(4) ,d)
      pmin(5)= min(pmin(5) ,e)
      pmax(5)= max(pmax(5) ,e)
      pmin(6)= min(pmin(6) ,f)
      pmax(6)= max(pmax(6) ,f)
      pmin(9)= min(pmin(9) ,g)
      pmax(9)= max(pmax(9) ,g)
      pmin(10)=min(pmin(10),h)
      pmax(10)=max(pmax(10),h)
      pmin(16)=min(pmin(16),p)
      pmax(16)=max(pmax(16),p)

!--ADDING (LOG OF) THE DISTANCES OF PHASE SPACE
      ia=ia-nstart

!--GET DIFFERENCE IN THE NUMBER OF TURNS PER DATA ENTRY : IDNT
      if(i1.eq.1) idnt=ia-ia0
      bold=bold+b
      b=b-b0
      dle1c=zero
      if(b.gt.zero) dle1c=log_mb(b)
      if(b.lt.zero) dle1c=-one*log_mb(-one*b)                            !hr06
      dle1=dle1+dle1c

!--EMITTANCES WITH LINEAR COUPLING CONT.
      evt1=evx1+evz1
      evxma=max(evx1,evxma)
      evzma=max(evz1,evzma)
      evtma=max(evt1,evtma)
      evxmi=min(evx1,evxmi)
      evzmi=min(evz1,evzmi)
      evtmi=min(evt1,evtmi)
      evx=evx+evx1
      evz=evz+evz1
      evt=evt+evt1

!--ADDING OF THE PHASE ADVANCES
      sx=c*d0-c0*d
      cx=c0*c+d*d0
      if(iapx.le.npos)                                                  &
     &call cphase(1,dphx,sx,cx,qx0,ivox,iwarx,iapx)
      sz=e*f0-e0*f
      cz=e0*e+f*f0
      if(iapz.le.npos)                                                  &
     &call cphase(2,dphz,sz,cz,qz0,ivoz,iwarz,iapz)
      s6=g*h0-g0*h
      c6=h0*h+g*g0
      if(iap6.le.npos)                                                  &
     &call cphase(3,dph6,s6,c6,qs0,ivo6,iwar6,iap6)

!--AVERAGING AFTER IAV TURNS
      if(mod(i1,iav).eq.0) then
        if(i2.ge.nlya) goto 240
        i2=i2+1
        dle(i2)=dle1/real(iav,fPrec)                                     !hr06
        if(ia.gt.0) then
          tle1=log_mb(real(ia,fPrec))                                          !hr06
          if(i2.gt.1) then
            biav(i2-1)=bold/real(iav,fPrec)                                    !hr06
            if(i2.eq.2) biav(1)=biav(1)*half
            bold=zero
            tle(i2)=(tle1+tlo)*half
            if(abs(tle1-tlo).gt.pieni) then
              wgh(i2)=one/(tle1-tlo)
            else
              write(lout,10310) nfile
              wgh(i2)=zero
            endif
          else
            tle(i2)=tle1*half
            wgh(i2)=one/(tle1)
          endif
        else
          tle(i2)=zero
          wgh(i2)=zero
        endif
        tlo=tle1
        dle1=zero
      endif
!--COORDINATE-ANGLE CONVERSION
      call caconv(dpx,d,c)
      call caconv(dpz,f,e)
      dpxp=tpi+dpx
      dpzp=tpi+dpz
!--INVARIANTS
      call cinvar(dpx,dphix,dpz,dpzp,nuex,emz,zinv,invz)
      call cinvar(dpz,dphiz,dpx,dpxp,nuez,emx,xinv,invx)
!--RESET OF COORDINATES
      c0=c
      d0=d
      e0=e
      f0=f
      g0=g
      h0=h
      goto 240
  270 if(i2.lt.1) i2=1

#ifdef CR
!--Now check that we have correct number of binrecs
!--We can do this only if we know binrecs (NOT post-processing only)
      if (binrec.ne.0) then
#ifndef STF
        if (binrecs(91-nfile).ne.crbinrecs(91-nfile)) then
          write(lout,"(a)") "SIXTRACR> ERROR POSTPR Wrong number of binary records"
          write(lout,"(a,i0,a,3(1x,i0))") "SIXTRACR> Unit ",nfile,", binrec/binrecs/crbinrecs ",&
            binrec,binrecs(91-nfile),crbinrecs(91-nfile)
          write(crlog,"(a)") "SIXTRACR> ERROR POSTPR Wrong number of binary records"
          write(crlog,"(a,i0,a,3(1x,i0))") "SIXTRACR> Unit ",nfile,", binrec/binrecs/crbinrecs ",&
            binrec,binrecs(91-nfile),crbinrecs(91-nfile)
#else
        if (binrecs(posi1).ne.crbinrecs(posi1)) then
          write(lout,"(a)") "SIXTRACR> ERROR POSTPR Wrong number of binary records"
          write(lout,"(a,i0,a,3(1x,i0))") "SIXTRACR> Particle ",posi1,", binrec/binrecs/crbinrecs ",&
            binrec,binrecs(posi1),crbinrecs(posi1)
          write(crlog,"(a)") "SIXTRACR> ERROR POSTPR Wrong number of binary records"
          write(crlog,"(a,i0,a,3(1x,i0))") "SIXTRACR> Particle ",posi1,", binrec/binrecs/crbinrecs ",&
            binrec,binrecs(posi1),crbinrecs(posi1)
#endif
          flush(crlog)
          goto 551
        endif
      endif
#endif

!----------------------------------------------------------------------
!--ANALYSING DATA
!----------------------------------------------------------------------
!--FIT OF DISTANCE IN PHASESPACE + MEAN PHASEADVANCE
      do 280 i=2,i2
        if(iwg.eq.1) call lfitwd(tle,dle,wgh,i,1,slope(i-1),const,varlea(i-1))
        if(iwg.eq.0) call lfitd(tle,dle,i,1,slope(i-1),const,varlea(i-1))
  280 continue
      if(iapx.eq.0) then
        write(lout,*) 'WARNING: IAPX IS ZERO'
        iapx=1
      endif
      if(iapz.eq.0) then
        write(lout,*) 'WARNING: IAPZ IS ZERO'
        iapz=1
      endif
      tphx=dphx/real(iapx,fPrec)                                         !hr06
      tphz=dphz/real(iapz,fPrec)                                         !hr06
      if(iap6.gt.0) tph6=dph6/real(iap6,fPrec)                           !hr06

!--STANDARD DEVIATION OF PHASEADVANCES
      do i=1,iapx
        sdpx=sdpx+(phase(1,i)-tphx)**2                                   !hr06
      end do

      do i=1,iapz
        sdpz=sdpz+(phase(2,i)-tphz)**2                                   !hr06
      end do

      do i=1,iap6
        sdp6=sdp6+(phase(3,i)-tph6)**2                                   !hr06
      end do

      sdpx=sqrt(sdpx)/real(iapx,fPrec)                                   !hr06
      sdpz=sqrt(sdpz)/real(iapz,fPrec)                                   !hr06
      if(iap6.gt.0) sdp6=sqrt(sdp6)/real(iap6,fPrec)                     !hr06

!--AVERAGED EMITTANCES
      di11=i11
      if(i11.eq.0) then
        write(lout,*) '** ERROR ** - I11 IS ZERO'
        goto 550
      endif
      emxa=emxa/di11
      emza=emza/di11
      emta=emta/di11
      evxm=evx/di11
      evzm=evz/di11
      evtm=evt/di11

!--SMEAR CALCULATION AND 4D-SMEAR
      rewind nfile
      !Skip headers
#ifndef STF
      read(nfile,iostat=ierro)
#else
      do i=1,itopa,2
         read(nfile,iostat=ierro)
      enddo
#endif
      if(ierro.gt.0) then
        write(lout,10320) nfile
        goto 550
      endif
      iskc=-1
      do 340 i=1,i11*iskip+nstart
        ifipa=0
        ! Read 1st particle only
#ifndef STF
 315    read(nfile,end=350,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp
        b=real(b_tmp,fPrec)
        c=real(c_tmp,fPrec)
        d=real(d_tmp,fPrec)
        e=real(e_tmp,fPrec)
        f=real(f_tmp,fPrec)
        g=real(g_tmp,fPrec)
        h=real(h_tmp,fPrec)
        p=real(p_tmp,fPrec)
#else
!     STF case: read tracking data until one reaches right particle.
 315    read(nfile,end=350,iostat=ierro) ia_stf,ifipa_stf,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp

        if(ifipa_stf.ne.posi) then
          goto 315
        endif
        ! Found the right particle; load data in memory (otherwise it's corrupted when EOF is reached)
        ia=ia_stf
        ifipa=ifipa_stf

        b=real(b_tmp,fPrec)
        c=real(c_tmp,fPrec)
        d=real(d_tmp,fPrec)
        e=real(e_tmp,fPrec)
        f=real(f_tmp,fPrec)
        g=real(g_tmp,fPrec)
        h=real(h_tmp,fPrec)
        p=real(p_tmp,fPrec)
#endif
        if(ierro.gt.0) then
           write(lout,10320) nfile
           goto 550
        endif

        if(ifipa.lt.1) goto 340
        if(progrm.eq.'MAD') then
#ifndef STF
           c=c*c1e3
           d=d*c1e3
           e=e*c1e3
           f=f*c1e3
           g=g*c1e3
           p=p*c1e3
#else
           write(lout,*) "ERROR in postpr: program=MAD not valid for STF."
           call prror(-1)
#endif
        endif !END if(program.eq.'MAD')

        iskc=iskc+1
        if(mod(iskc,iskip).ne.0) goto 340
        if((ia-nstart).lt.0) goto 340
        c=c-clo(1)
        d=d-clop(1)
        e=e-clo(2)
        f=f-clop(2)
        g=g-clo(3)
        h=h-clop(3)

        if(icode.ge.4) then
          c=c-di0(1)*h
          d=d-dip0(1)*h
          e=e-di0(2)*h
          f=f-dip0(2)*h
        endif

!--MEAN EMITTANCES
        xp=bet0(1)*d+alf0(1)*c
        zp=bet0(2)*f+alf0(2)*e
        emx=(c**2+xp**2)/bet0(1)                                         !hr06
        emz=(e**2+zp**2)/bet0(2)                                         !hr06

        if(icode.ge.4.and.its6d.ne.0) then
          c=c+di0(1)*h
          d=d+dip0(1)*h
          e=e+di0(2)*h
          f=f+dip0(2)*h
        endif

        emt=emx+emz
        emxs=emxs+(emx-emxa)**2                                          !hr06
        emzs=emzs+(emz-emza)**2                                          !hr06
        emts=emts+(emt-emta)**2                                          !hr06

!--COURANT SYNDER
        xyzv(1)=c
        xyzv(2)=d
        xyzv(3)=e
        xyzv(4)=f
        xyzv(5)=g
        xyzv(6)=h

!--CONVERT TO CANONICAL VARIABLES
        if(its6d.eq.1) then
          xyzv(2)=xyzv(2)*((one+xyzv(6))+clop(3))                        !hr06
          xyzv(4)=xyzv(4)*((one+xyzv(6))+clop(3))                        !hr06
        endif

! normalisation with t-matrix = inverse matrix of eigenvectors
        do iq=1,6
          txyz(iq)=zero
          do jq=1,6
            txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
          end do
        end do

!--MEAN EMITTANCES WITH LINEAR COUPLING
        evx=txyz(1)**2+txyz(2)**2                                        !hr06
        evz=txyz(3)**2+txyz(4)**2                                        !hr06

!--COURANT SYNDER CONT.
        do iq=1,6
          txyz(iq)=txyz(iq)*rbeta(iq)
        end do

        c=txyz(1)
        d=txyz(2)
        e=txyz(3)
        f=txyz(4)
        g=txyz(5)
        h=txyz(6)

!--MEAN EMITTANCES WITH LINEAR COUPLING CONT.
        evt=evx+evz
        sevx=sevx+(evx-evxm)**2                                          !hr06
        sevz=sevz+(evz-evzm)**2                                          !hr06
        sevt=sevt+(evt-evtm)**2                                          !hr06
  340 continue
  350 continue

!--SMEAR IN %
      call sinpro(emxa,di11,emxs,emax,emix)
      call sinpro(emza,di11,emzs,emaz,emiz)
      call sinpro(emta,di11,emts,emat,emit)
      call sinpro(evxm,di11,sevx,evxma,evxmi)
      call sinpro(evzm,di11,sevz,evzma,evzmi)
      call sinpro(evtm,di11,sevt,evtma,evtmi)

!----------------------------------------------------------------------
!--PRINTING
!----------------------------------------------------------------------
      if(nstop.lt.ia.and.(ia.lt.numl.or.ia.lt.nint(dnumlr))) nlost=1
      if(nnumxv(ifipa).eq.0.and.nnumxv(ilapa).eq.0) then
        sumda(22)=real(ia,fPrec)                                         !hr06
        sumda(23)=real(ia,fPrec)                                         !hr06
      else
        sumda(22)=real(nnumxv(ifipa),fPrec)                              !hr06
        sumda(23)=real(nnumxv(ilapa),fPrec)                              !hr06
      endif
! #ifdef SIXDA
!       sumda(22)=real(ia,fPrec)                                           !hr06
!       sumda(23)=real(ia,fPrec)                                           !hr06
! #endif
#ifdef CR
! TRY a FIX for nnuml
! should be redumdant now
!     if (nnuml.ne.numl) then
!       if (nint(sumda(22)).eq.numl) sumda(22)=dble(nnuml)
!       if (nint(sumda(23)).eq.numl) sumda(23)=dble(nnuml)
!     endif
#endif
      sumda(2)=real(nlost,fPrec)
      sumda(9)=dp1-clop(3)

!--GET DIFFERENCE IN THE NUMBER OF TURNS PER DATA ENTRY : TIDNT
!--NOW CONSIDERING ONLY TURNS LARGER THAN NSTART
      tidnt=real(((ia-nstart)+idnt)/i11,fPrec)                           !hr06
      if(i2.ge.2) then
        if(nprint.eq.1) write(lout,10110)
        ilyap=0
        slopem=zero

        do i=1,i2-1
          iturn=nint(real((i+1)*iav,fPrec)*tidnt)                        !hr06
          if(nprint.eq.1) write(lout,10120) iturn,biav(i),slope(i),varlea(i)
          if(biav(i).gt.c1m1) ilyap=1
          slopem=max(slopem,slope(i))
        end do

        if(nprint.eq.1) write(lout,10130)
        sumda(10)=biav(i2-1)

        if(ilyap.eq.0) then
         sumda(11)=slope(i2-1)                                           !hr06
        else
         sumda(11)=slopem
        endif
      endif

!--CALCULATION OF AVERAGED PHASEADVANCES
      tph6=abs(tph6)
      if(nprint.eq.1) then
        write(lout,10140)tphx,sdpx,tphz,sdpz,tph6,sdp6,qwc(1),tphx-qwc(1),qwc(2),tphz-qwc(2),qwc(3),tph6-qwc(3),dres,ires
      end if

      sumda(3)=qwc(1)
      sumda(4)=qwc(2)

      if(abs(tphx).gt.pieni) then
        sumda(12)=tphx-qwc(1)
      else
        sumda(12)=zero
      endif

      sumda(13)=sdpx

      if(abs(tphz).gt.pieni) then
        sumda(14)=tphz-qwc(2)
      else
        sumda(14)=zero
      endif

      sumda(15)=sdpz
      sumda(25)=tph6

!--DISTANCE OF Q-VALUES (AVERAGED PHASEADVANCE) TO RESONANCES
      do i=1,21
        do j=1,21
          im1=i-1
          jm1=j-1
          if(im1.eq.0.and.jm1.eq.0) goto 370
          if(im1+jm1.gt.ires) goto 370
          ares=real(im1,fPrec)*tphx+real(jm1,fPrec)*tphz                 !hr06
          dares=anint(ares)
          ares=ares-dares
          if(abs(ares).lt.armin) then
            armin=abs(ares)
            im1s=im1
            jm1s=jm1
          endif
          ared=real(im1,fPrec)*tphx-real(jm1,fPrec)*tphz                 !hr06
          dared=anint(ared)
          ared=ared-dared
          if(abs(ared).lt.armin) then
            armin=abs(ared)
            im1s=im1
            jm1s=-jm1
          endif
          if(abs(ares).lt.dres.and.nprint.eq.1) write(lout,10170) im1,jm1,dares,ares
          if(abs(ared).lt.dres.and.jm1.ne.0.and.im1.ne.0.and.nprint.eq.1) write(lout,10170) im1,-jm1,dared,ared
  370 continue
        end do
      end do


      if(armin.lt.armin0) then
        sumda(16)=real(im1s,fPrec)                                       !hr06
        sumda(17)=real(jm1s,fPrec)                                       !hr06
        sumda(18)=sumda(16)+abs(sumda(17))
      endif
      if(iwarx.eq.1.and.nprint.eq.1) write(lout,10150)
      if(iwarz.eq.1.and.nprint.eq.1) write(lout,10160)

!--Q-VALUES BY AN FFT-ROUTINE
  380 ifp=ifp+1
      ife=2**ifp

      if(ife.le.i11.and.ife.le.nfft) then
        goto 380
      else
        ifp=ifp-1
        ife=ife/2
      endif

      if(ife.eq.0) then
        write(lout,*) '** ERROR ** - IFE IS ZERO'
        goto 550
      endif

      dife=ife

      if(ifp.gt.1) then
        if(nprint.eq.1) write(lout,10180) ife,dfft*100
        call fft(xxr,xxi,ifp,ife)
        call fft(zzr,zzi,ifp,ife)
        xxmax=zero
        zzmax=zero
        xxmin=one
        zzmin=one

        if(ifh.eq.0) then
          if1=1
          if2=ife
          ife2=ife
          pmin(21)=qx0
          pmax(21)=qx0+one
          pmin(23)=qz0
          pmax(23)=qz0+one
        else if(ifh.eq.1) then
          if1=1
          if2=ife/2
          ife2=ife/2
          pmin(21)=qx0
          pmax(21)=qx0+half
          pmin(23)=qz0
          pmax(23)=qz0+half
        else
          if1=ife/2+1
          if2=ife
          ife2=ife/2
          pmin(21)=qx0+half
          pmax(21)=qx0+one
          pmin(23)=qz0+half
          pmax(23)=qz0+one
        endif

        do 390 i=if1,if2
          xxmax=max(xxmax,sqrt(xxr(i)**2+xxi(i)**2))
          zzmax=max(zzmax,sqrt(zzr(i)**2+zzi(i)**2))
          xxmin=min(xxmin,sqrt(xxr(i)**2+xxi(i)**2))
          zzmin=min(zzmin,sqrt(zzr(i)**2+zzi(i)**2))
  390   continue

        if(abs(xxmax).gt.pieni) xxmin=xxmin/xxmax
        if(abs(zzmax).gt.pieni) zzmin=zzmin/zzmax

        if(xxmax.le.pieni) then
          write(lout,*) 'WARNING: XXMAX IS SET TO : ',pieni
          xxmax=pieni
        endif

        if(zzmax.le.pieni) then
          write(lout,*) 'WARNING: ZZMAX IS SET TO : ',pieni
          zzmax=pieni
        endif

        do 400 i=if1,if2
          xxaux=sqrt(xxr(i)**2+xxi(i)**2)
          zzaux=sqrt(zzr(i)**2+zzi(i)**2)
          if(abs(xxaux-xxmax).le.pieni) ffx=(real(i-1,fPrec)/dife)+qx0         !hr06
          if(abs(zzaux-zzmax).le.pieni) ffz=(real(i-1,fPrec)/dife)+qz0         !hr06
          xxaux=xxaux/xxmax
          zzaux=zzaux/zzmax
          if(xxaux.gt.dfft.and.nprint.eq.1) write(lout,10190) real(i-1,fPrec)/dife+qx0,xxaux*c1e2  !hr06
      if(zzaux.gt.dfft.and.nprint.eq.1) write(lout,10200) real(i-1,fPrec)/dife+qz0,zzaux*c1e2      !hr06
  400   continue

        if(nprint.eq.1) write(lout,10210) ffx,ffz,qwc(1),ffx-qwc(1),qwc(2),ffz-qwc(2),dres,ires

!--DISTANCE OF Q-VALUES (FFT) TO RESONANCES
        do i=1,21
          do j=1,21
            im1=i-1
            jm1=j-1
            if(im1.eq.0.and.jm1.eq.0) goto 410
            if(im1+jm1.gt.ires) goto 410
            ares=real(im1,fPrec)*ffx+real(jm1,fPrec)*ffz                             !hr06
            dares=anint(ares)
            ares=ares-dares
            ared=real(im1,fPrec)*ffx-real(jm1,fPrec)*ffz                             !hr06
            dared=anint(ared)
            ared=ared-dared
            if(abs(ares).lt.dres.and.nprint.eq.1) write(lout,10170) im1,jm1,dares,ares
            if(abs(ared).lt.dres.and.jm1.ne.0.and.im1.ne.0.and.nprint.eq.1) write(lout,10170) im1,-jm1,dared,ared
  410   continue
          end do
        end do
      endif

!--PRINT 4-D INVARIANTS WITH LINEAR COUPLING
      if(nprint.eq.1) write(lout,10270) emi,emii,emiii,angi,angii,angiii,evxm,sevx,evxma,evxmi,evzm,sevz,evzma,evzmi,evtm,sevt,  &
     & evtma,evtmi

!--PRINT EMITTANCES AND SMEAR
      ampx0=sqrt(bet0(1)*emx0)
      ampz0=sqrt(bet0(2)*emz0)
      if(nprint.eq.1) write(lout,10220) emx0,ampx0,emz0,ampz0,emxa,emxs,emax,emix,emza,emzs,emaz,emiz,emta,emts,emat,emit
      sumda(46)=emi
      sumda(47)=emii
      sumda(48)=bet0x2
      sumda(49)=bet0z2
      sumda(7)=sqrt(bet0(1)*emi)+sqrt(bet0x2*emii)
      sumda(8)=sqrt(bet0(2)*emii)+sqrt(bet0z2*emi)
      sumda(26)=sqrt(bet0(1)*evx2)+sqrt(bet0x2*evz2)
      sumda(27)=sqrt(bet0(2)*evz2)+sqrt(bet0z2*evx2)
      sumda(19)=sevx
      sumda(20)=sevz
      sumda(21)=sevt
!     sumda(59)=dmmac
!     sumda(60)=dnms
      sumda(59)=dnms
! This place 60 now used for CPU time seconds
! But it is set earlier in case particles are lost very early
      sumda(24)=dizu0
      emax=(emax/c1e2)*emxa+emxa                                        !hr06
      emix=(emix/c1e2)*emxa+emxa                                        !hr06
      emaz=(emaz/c1e2)*emza+emza                                        !hr06
      emiz=(emiz/c1e2)*emza+emza                                        !hr06
      sumda(28)=sqrt(bet0(1)*abs(emix))
      sumda(29)=sqrt(bet0(1)*emxa)
      sumda(30)=sqrt(bet0(1)*emax)
      sumda(31)=sqrt(bet0(2)*abs(emiz))
      sumda(32)=sqrt(bet0(2)*emza)
      sumda(33)=sqrt(bet0(2)*emaz)
      evxma=(evxma/c1e2)*evxm+evxm                                      !hr06
      evxmi=(evxmi/c1e2)*evxm+evxm                                      !hr06
      evzma=(evzma/c1e2)*evzm+evzm                                      !hr06
      evzmi=(evzmi/c1e2)*evzm+evzm                                      !hr06
      sumda(34)=sqrt(bet0(1)*abs(evxmi))
      sumda(35)=sqrt(bet0(1)*evxm)
      sumda(36)=sqrt(bet0(1)*evxma)
      sumda(37)=sqrt(bet0(2)*abs(evzmi))
      sumda(38)=sqrt(bet0(2)*evzm)
      sumda(39)=sqrt(bet0(2)*evzma)
      evtma=(evtma/c1e2)*evtm+evtm                                      !hr06
      evtmi=(evtmi/c1e2)*evtm+evtm                                      !hr06

      if(abs(evxm+evzm).gt.pieni) then
        ratemx=evxm/(evxm+evzm)
        ratemz=evzm/(evxm+evzm)
      else
        ratemx=zero
        ratemz=zero
      endif

      sumda(40)=sqrt((bet0(1)*abs(evtmi))*ratemx)                        !hr06
      sumda(41)=sqrt((bet0(1)*evtm)*ratemx)                              !hr06
      sumda(42)=sqrt((bet0(1)*evtma)*ratemx)                             !hr06
      sumda(43)=sqrt((bet0(2)*abs(evtmi))*ratemz)                        !hr06
      sumda(44)=sqrt((bet0(2)*evtm)*ratemz)                              !hr06
      sumda(45)=sqrt((bet0(2)*evtma)*ratemz)                             !hr06

!--PUT IN THE CHROMATICITY
      sumda(50)=chromc(1)*c1e3
      sumda(51)=chromc(2)*c1e3

!--WRITE DATA FOR THE SUMMARY OF THE POSTPROCESSING ON FILE # 10
! We should really write fort.10 in BINARY!
      write(unit110,iostat=ierro) (sumda(i),i=1,60)
      if(ierro /= 0) then
        write(lerr,"(a,i0)") "POSTPR> ERROR Problems writing to "//trim(fort110)//". iostat = ",ierro
      end if
#ifndef CRLIBM
      write(ch,*,iostat=ierro) (sumda(i),i=1,60)
      write(unit10,"(a)",iostat=ierro) trim(ch)
#else
      l1=1
      do i=1,60
        call chr_fromReal(sumda(i),ch1,19,2,rErr)
        ch(l1:l1+25) = " "//ch1(1:25)
        l1=l1+26
      end do
      write(unit10,"(a)",iostat=ierro) ch(1:l1-1)
#endif
      if(ierro /= 0) then
        write(lerr,"(a,i0)") "POSTPR> ERROR Problems writing to "//trim(fort10)//". iostat = ",ierro
      end if
!--CALCULATION THE INVARIANCES OF THE 4D TRANSVERSAL MOTION
      do 420 i=1,ninv
        if(invx(i).gt.0) then
          nuix=nuix+1
          xing=xing+xinv(i)/real(invx(i),fPrec)                          !hr06
        endif

        if(invz(i).gt.0) then
          nuiz=nuiz+1
          zing=zing+zinv(i)/real(invz(i),fPrec)                          !hr06
        endif
  420 continue

      pinx=real(nuix,fPrec)                                              !hr06
      pinz=real(nuiz,fPrec)                                              !hr06

      if(nuix.ne.0) then
        pixr=real(nuez,fPrec)/real(nuix,fPrec)
        xing=xing/real(nuix,fPrec)                                       !hr06
      endif
      if(nuiz.ne.0) then
        pizr=real(nuex,fPrec)/real(nuiz,fPrec)
        zing=zing/real(nuiz,fPrec)                                       !hr06
      endif
      pinx=(pinx/real(ninv,fPrec))*c1e2                                  !hr06
      pinz=(pinz/real(ninv,fPrec))*c1e2                                  !hr06
      if(nprint.eq.1) write(lout,10230)
      if(nuez.lt.ninv.and.nprint.eq.1) write(lout,10240) nuez,ninv
      if(nuex.lt.ninv.and.nprint.eq.1) write(lout,10250) nuex,ninv
      if(nprint.eq.1) write(lout,10260) nuez,nuix,nuex,nuiz, ninv,pinx,pixr,pinz,pizr,xing,zing

!----------------------------------------------------------------------
!--PLOTTING
!----------------------------------------------------------------------
      pmin(1)=zero
      pmax(1)=real(ia,fPrec)                                             !hr06
      pmin(7)=pmin(3)
      pmax(7)=pmax(3)
      pmin(8)=pmin(5)
      pmax(8)=pmax(5)
      pmin(11)=pmin(3)
      pmax(11)=pmax(3)
      pmin(12)=pmin(10)
      pmax(12)=pmax(10)
      pmin(13)=pmin(5)
      pmax(13)=pmax(5)
      pmin(14)=pmin(10)
      pmax(14)=pmax(10)
      pmax(15)=real(ia,fPrec)                                            !hr06
      pmin(17)=pmin(3)
      pmax(17)=pmax(3)
      pmin(18)=pmin(4)
      pmax(18)=pmax(4)
      pmin(19)=pmin(5)
      pmax(19)=pmax(5)
      pmin(20)=pmin(6)
      pmax(20)=pmax(6)
      pmin(22)=zero
      pmin(24)=zero
      pmax(22)=one
      pmax(24)=one

      do 500 i=1,12
        i2=2*i
        i1=i2-1
        if(pmin(i1).gt.pmax(i1)) pmin(i1)=pmax(i1)
        if(pmin(i2).gt.pmax(i2)) pmin(i2)=pmax(i2)
        if((abs(pmin(i1)-pmax(i1)).le.pieni2) .or.(abs(pmin(i2)-pmax(i2)).le.pieni2)) then
          goto 500
        endif

        do 430 i3=i1,i2
          pcha=(pmax(i3)-pmin(i3))*prec
          pmin(i3)=pmin(i3)-pcha
          pmax(i3)=pmax(i3)+pcha
  430   continue

        if(iffw.eq.2) then
          pmin(22)=xxmin/(one+abs(prec))
          pmin(24)=zzmin/(one+abs(prec))
        endif

        if((i.eq.1.and.idis.eq.1).or. (i.gt.1.and.i.le.8.and.icow.eq.1) &
     &.or. ((i.eq.9.or.i.eq.10).and.istw.eq.1).and. (pmin(i1).ne.pmax   &
     &(i1).and.pmin(i2).ne.pmax(i2))) then

!--HBOOK FRAME
!using http://cds.cern.ch/record/118642
! defines "general title"
          call htitle(title(i))
! hbook 2 (2d hist) id, title, nx, xmin, xmax, ny, ymin, ymax, data_min - "bit precision" ?
          call hbook2(i,' ',2,real(pmin(i1)),real(pmax(i1)), 2,real(pmin(i2)),real(pmax(i2)),0.)
          call hplot(i,' ',' ',0)
          call hplax(chxtit(i),chytit(i))
          call hplsof(4.,14.75,toptit(1),.15,0.,99.,-1)
          call hplsof(4.,14.50,toptit(2),.15,0.,99.,-1)
          call hplsof(4.,14.25,toptit(3),.15,0.,99.,-1)
          call hplsof(4.,14.00,toptit(4),.15,0.,99.,-1)
          call iselnt(10)

          rewind nfile
          !Skip headers
#ifdef STF
          do j=1,itopa,2
#endif
          read(nfile,iostat=ierro)
#ifdef STF
          enddo
#endif
          if(ierro.gt.0) then
            write(lout,10320) nfile
#ifdef CR
            goto 551
#else
            goto 550
#endif
          endif
          iskc=-1
          do 460 j=1,i11*iskip+nstart
 435        ifipa=0
#ifndef STF
            if(ntwin.eq.1) then
               read(nfile,end=470,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp

               b=real(b_tmp,fPrec)
               c=real(c_tmp,fPrec)
               d=real(d_tmp,fPrec)
               e=real(e_tmp,fPrec)
               f=real(f_tmp,fPrec)
               g=real(g_tmp,fPrec)
               h=real(h_tmp,fPrec)
               p=real(p_tmp,fPrec)

            elseif(ntwin.eq.2) then
               read(nfile,end=470,iostat=ierro) ia,ifipa,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp, &
     & ilapa,b_tmp,c1_tmp,d1_tmp,e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp

               b=real(b_tmp,fPrec)
               c=real(c_tmp,fPrec)
               d=real(d_tmp,fPrec)
               e=real(e_tmp,fPrec)
               f=real(f_tmp,fPrec)
               g=real(g_tmp,fPrec)
               h=real(h_tmp,fPrec)
               p=real(p_tmp,fPrec)

               c1=real(c1_tmp,fPrec)
               d1=real(d1_tmp,fPrec)
               e1=real(e1_tmp,fPrec)
               f1=real(f1_tmp,fPrec)
               g1=real(g1_tmp,fPrec)
               h1=real(h1_tmp,fPrec)
               p1=real(p1_tmp,fPrec)
            end if
#else
!     STF case: read tracking data until one reaches right particle.
            if(ntwin.eq.1) then
               read(nfile,end=470,iostat=ierro) ia_stf,ifipa_stf,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp
            elseif(ntwin.eq.2) then
               read(nfile,end=470,iostat=ierro) ia_stf,ifipa_stf,b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp, &
     &              ilapa_stf,b_tmp,c1_tmp,d1_tmp,e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp
            end if

            if(ifipa_stf.ne.posi) then
              goto 435
            end if

            ! Found the right particle; load data in memory (otherwise it's corrupted when EOF is reached)
            ia=ia_stf
            ifipa=ifipa_stf

            b=real(b_tmp,fPrec)
            c=real(c_tmp,fPrec)
            d=real(d_tmp,fPrec)
            e=real(e_tmp,fPrec)
            f=real(f_tmp,fPrec)
            g=real(g_tmp,fPrec)
            h=real(h_tmp,fPrec)
            p=real(p_tmp,fPrec)

            if(ntwin.eq.2) then
               ilapa=ilapa_stf

               c1=real(c1_tmp,fPrec)
               d1=real(d1_tmp,fPrec)
               e1=real(e1_tmp,fPrec)
               f1=real(f1_tmp,fPrec)
               g1=real(g1_tmp,fPrec)
               h1=real(h1_tmp,fPrec)
               p1=real(p1_tmp,fPrec)
            end if
#endif
            if(ierro.gt.0) then
              write(lout,10320) nfile
              goto 550
            end if

            if(ifipa.lt.1) goto 460
            iskc=iskc+1
            if(mod(iskc,iskip).ne.0) goto 460
            if((ia-nstart).lt.0) goto 460
            if(progrm.eq.'MAD') then !NON-STF only
#ifndef STF
              c=c*c1e3
              d=d*c1e3
              e=e*c1e3
              f=f*c1e3
              g=g*c1e3
              p=p*c1e3

              if(ntwin.eq.2) then
                c1=c1*c1e3
                d1=d1*c1e3
                e1=e1*c1e3
                f1=f1*c1e3
                g1=g1*c1e3
                p1=p1*c1e3
              end if
#else
        write(lout,*) "ERROR in postpr: program=MAD not valid for STF."
        call prror(-1)
#endif

            end if
!--LYAPUNOV
            if(ntwin.eq.2) then
              x(1,1)=c
              x(1,2)=d
              x(1,3)=e
              x(1,4)=f
              x(1,5)=g
              x(1,6)=h
              x(2,1)=c1
              x(2,2)=d1
              x(2,3)=e1
              x(2,4)=f1
              x(2,5)=g1
              x(2,6)=h1
              call distance(x,cloau,di0au,t,b)
            endif

            if(icode.ge.4.and.its6d.eq.0) then
              c=c-di0(1)*h
              d=d-dip0(1)*h
              e=e-di0(2)*h
              f=f-dip0(2)*h
            endif

            c=c-clo(1)
            d=d-clop(1)
            e=e-clo(2)
            f=f-clop(2)
            g=g-clo(3)
            h=h-clop(3)
            xyzv(1)=c
            xyzv(2)=d
            xyzv(3)=e
            xyzv(4)=f
            xyzv(5)=g
            xyzv(6)=h
!--CONVERT TO CANONICAL VARIABLES
            if(its6d.eq.1) then
              xyzv(2)=xyzv(2)*((one+xyzv(6))+clop(3))                    !hr06
              xyzv(4)=xyzv(4)*((one+xyzv(6))+clop(3))                    !hr06
            endif

            do iq=1,6
              txyz(iq)=zero
              do jq=1,6
                txyz(iq)=txyz(iq)+t(jq,iq)*xyzv(jq)
              end do
            end do

            do iq=1,6
              txyz(iq)=txyz(iq)*rbeta(iq)
            end do

            c=txyz(1)
            d=txyz(2)
            e=txyz(3)
            f=txyz(4)
            g=txyz(5)*cma2
            h=txyz(6)*cma1
            if(idis.eq.1.and.i.eq.1) call ipm(1,real(ia),real(b))

            if(icow.eq.1.and.i.le.8) then
              if(i.eq.2) call ipm(1,real(c),real(d))
              if(i.eq.3) call ipm(1,real(e),real(f))
              if(i.eq.4) call ipm(1,real(c),real(e))
              if(i.eq.5) call ipm(1,real(g),real(h))
              if(i.eq.6) call ipm(1,real(c),real(h))
              if(i.eq.7) call ipm(1,real(e),real(h))
              if(i.eq.8) call ipm(1,real(ia),real(p))
            endif

            if(istw.eq.1.and.(i.eq.9.or.i.eq.10)) then
              if(i.eq.9) then
                call caconv(dpz,f,e)
                if(abs(dpz).lt.dphiz) call ipm(1,real(c),real(d))
              endif
              if(i.eq.10) then
                call caconv(dpx,d,c)
                if(abs(dpx).lt.dphix) call ipm(1,real(e),real(f))
              endif
            endif
  460     continue
  470     continue
        else if((iffw.eq.1.or.iffw.eq.2).and.(i.eq.11.or.i.eq.12)) then

!--HBOOK FRAME
          call htitle(title(i))
          call hbook2(i,' ',2,real(pmin(i1)),real(pmax(i1)), 2,real(pmin(i2)),real(pmax(i2)),0.)
          if(iffw.eq.2) call hplopt('LOGY',1)
          call hplot(i,' ',' ',0)
          call hplax(chxtit(i),chytit(i))
          call hplsof(4.,14.75,toptit(1),.15,0.,99.,-1)
          call hplsof(4.,14.50,toptit(2),.15,0.,99.,-1)
          call hplsof(4.,14.25,toptit(3),.15,0.,99.,-1)
          call hplsof(4.,14.00,toptit(4),.15,0.,99.,-1)
          call iselnt(10)
          if(i.eq.11) then
            do 480 k=if1,if2
              k1=k-if1+1
              xxaux=sqrt(xxr(k)**2+xxi(k)**2)
              xxaux=xxaux/xxmax
              fxs(k1)=real(dble(k-1)/dife+qx0)                           !hr06
              if(iffw.eq.2) then
                if(abs(xxaux).lt.pieni) then
                  write(lout,*) '* * * ERROR * * *'
                  write(lout,*) 'Apparently horizontal FFT data are corrupted'
                  xxaux=one
                endif
                fzs(k1)=real(log10_mb(xxaux)) !! REAL ???? CHECKME !!!
              else
                fzs(k1)=real(xxaux)
              endif
              if(nprint.eq.1) then
                write(14,10030,iostat=ierro) fxs(k1),fzs(k1)
                if(ierro.ne.0) then
                  write(lout,*)
                  write(lout,*) '*** ERROR ***,PROBLEMS WRITING TO FILE # : ',14
                  write(lout,*) 'ERROR CODE : ',ierro
                  write(lout,*)
                endif
              endif
  480       continue
            call ipl(ife2,fxs,fzs)
          else if(i.eq.12) then
            do 490 k=if1,if2
              k1=k-if1+1
              zzaux=sqrt(zzr(k)**2+zzi(k)**2)
              zzaux=zzaux/zzmax
              fxs(k1)=real(real(k-1,fPrec)/dife+qz0)                           !hr06

              if(iffw.eq.2) then
                if(abs(zzaux).lt.pieni) then
                  write(lout,*) '* * * ERROR * * *'
                  write(lout,*) 'Apparently vertical FFT data are corrupted'
                  zzaux=one
                endif
                fzs(k1)=real(log10_mb(zzaux))
              else
                fzs(k1)=real(zzaux)
              endif

              if(nprint.eq.1) then
                write(15,10030,iostat=ierro) fxs(k1),fzs(k1)
                if(ierro.ne.0) then
                  write(lout,*)
                  write(lout,*) '*** ERROR ***,PROBLEMS WRITING TO FILE # : ',14
                  write(lout,*) 'ERROR CODE : ',ierro
                  write(lout,*)
                endif
              endif
  490       continue
            call ipl(ife2,fxs,fzs)
          endif
        endif
        if(iffw.eq.2) call hplopt('LINY',1)
  500 continue
      if(idis.ne.0.or.icow.ne.0.or.istw.ne.0.or.iffw.ne.0) call hdelet(0)
      goto 560
  510 continue
      write(lout,10300) nfile,'HEADER CORRUPTED'
      goto 550
  520 continue
      write(lout,10300) nfile,'HEADER OF MADFILE CORRUPTED'
      goto 550
  530 continue
      write(lout,10300) nfile,'NO DATA'
      goto 550
  535 continue
      write(lout,10300) nfile,'ONLY START VALUES'
      goto 550
  540 continue
      write(lout,10300) nfile,'WRONG RANGE OF DATA FOR PROCESSING'
      goto 550
#ifdef CR
  551 write(crlog,"(a)") "SIXTRACR> ERROR POSTPR"
      flush(crlog)
! Now we let abend handle the fort.10......
! It will write 0d0 plus CPU time and turn number
! But we empty it as before (if we crash in abend???)
      rewind(unit10)
      endfile(unit10,iostat=ierro)
      call f_close(unit10)
      write(lerr,"(a)") "SIXTRACR> ERROR POSTPR"
      call prror
#endif

 550  continue
!--WRITE DATA FOR THE SUMMARY OF THE POSTPROCESSING ON FILE # 10
!-- Will almost all be zeros but we now have napxto and ttime
! We should really write fort.10 in BINARY!
      write(unit110,iostat=ierro) (sumda(i),i=1,60)
      if(ierro /= 0) then
        write(lerr,"(a,i0)") "POSTPR> ERROR Problems writing to file 110. Error code ",ierro
      endif
#ifndef CRLIBM
      write(ch,*,iostat=ierro) (sumda(i),i=1,60)
      write(unit10,"(a)",iostat=ierro) trim(ch)
#else
      l1=1
      do i=1,60
        call chr_fromReal(sumda(i),ch1,19,2,rErr)
        ch(l1:l1+25) = " "//ch1(1:25)
        l1=l1+26
      enddo
      write(unit10,"(a)",iostat=ierro) ch(1:l1-1)
#endif
      if(ierro /= 0) then
        write(lerr,"(a,i0)") "POSTPR> ERROR Problems writing to "//trim(fort10)//". iostat = ",ierro
      end if
!--REWIND USED FILES
  560 rewind nfile
      if(nprint == 1) then
        rewind 14
        rewind 15
      end if
!--TIME COUNT
      tim2=0.
      call time_timerCheck(tim2)
      if(nprint.eq.1) write(lout,10280) tim2-tim1
!----------------------------------------------------------------------
      return

10000 format(d11.4)
10010 format(f10.6)
10020 format(a80)
10030 format(2f10.6)
10040 format( //131('-')//10x,'OOOOOOOOOOOOOOOOOOOOOO' /10x,            &
     &'OO                  OO' /10x,'OO  POSTPROCESSING  OO' /10x,      &
     &'OO                  OO' /10x,'OOOOOOOOOOOOOOOOOOOOOO'// /10x,    &
     &'TITLE AND COMMENT :'//a80//a80// )
10050 format(10x,'THE FOLLOWING PARAMETERS ARE USED:'//                 &
     &10x,'PROGRAM NAME',t102,a8/                                       &
     &10x,'PARTICLE NUMBER',t102,i7/                                    &
     &10x,'TOTAL NUMBER OF PARTICLES',t102,i7/                          &
     &10x,'PHASE SPACE',t102,a11/                                       &
     &10x,'MAXIMUM NUMBER OF TURNS',t102,i8/                            &
     &10x,'HORIZONTAL BETA',t102,f16.10/                                &
     &10x,'HORIZONTAL BETA-II',t102,f16.10/                             &
     &10x,'HORIZONTAL BETA-III',t102,f16.10/                            &
     &10x,'VERTICAL BETA',t102,f16.10/                                  &
     &10x,'VERTICAL BETA-II',t102,f16.10/                               &
     &10x,'VERTICAL BETA-III',t102,f16.10/                              &
     &10x,'LONGITUDINAL BETA',t98,f20.10/                               &
     &10x,'LONGITUDINAL BETA-II',t102,f16.10/                           &
     &10x,'LONGITUDINAL BETA-III',t102,f16.10/                          &
     &10x,'HORIZONTAL ALFA',t102,f16.10/                                &
     &10x,'HORIZONTAL ALFA-II',t102,f16.10/                             &
     &10x,'HORIZONTAL ALFA-III',t102,f16.10)
10060 format(                                                           &
     &10x,'VERTICAL ALFA',t102,f16.10/                                  &
     &10x,'VERTICAL ALFA-II',t102,f16.10/                               &
     &10x,'VERTICAL ALFA-III',t102,f16.10/                              &
     &10x,'LONGITUDINAL ALFA',t102,f16.10/                              &
     &10x,'LONGITUDINAL ALFA-II',t102,f16.10/                           &
     &10x,'LONGITUDINAL ALFA-III',t102,f16.10/                          &
     &10x,'HORIZONTAL GAMMA',t102,f16.10/                               &
     &10x,'HORIZONTAL GAMMA-II',t102,f16.10/                            &
     &10x,'HORIZONTAL GAMMA-III',t102,f16.10/                           &
     &10x,'VERTICAL GAMMA',t102,f16.10/                                 &
     &10x,'VERTICAL GAMMA-II',t102,f16.10/                              &
     &10x,'VERTICAL GAMMA-III',t102,f16.10/                             &
     &10x,'LONGITUDINAL GAMMA',t102,f16.10/                             &
     &10x,'LONGITUDINAL GAMMA-II',t102,f16.10/                          &
     &10x,'LONGITUDINAL GAMMA-III',t102,f16.10/                         &
     &10x,'HORIZONTAL CLOSED ORBIT',t105,ES17.10/                       &
     &10x,'VERTICAL CLOSED ORBIT',t105,ES17.10/                         &
     &10x,'LONGITUDINAL CLOSED ORBIT',t105,ES17.10/                     &
     &10x,'SLOPE OF HORIZONTAL CLOSED ORBIT',t105,ES17.10/              &
     &10x,'SLOPE OF VERTICAL CLOSED ORBIT',t105,ES17.10/                &
     &10x,'SLOPE OF LONGITUDINAL CLOSED ORBIT',t105,ES17.10/            &
     &10x,'HORIZONTAL DISPERSION',t105,ES17.10/                         &
     &10x,'VERTICAL DISPERSION',t105,ES17.10/                           &
     &10x,'SLOPE OF HORIZONTAL DISPERSION',t105,ES17.10/                &
     &10x,'SLOPE OF VERTICAL DISPERSION',t105,ES17.10/                  &
     &10x,'LINEAR HORIZONTAL TUNE',t105,ES17.10/                        &
     &10x,'LINEAR VERTICAL TUNE',t105,ES17.10/                          &
     &10x,'LINEAR LONGITUDINAL TUNE',t105,ES17.10)
10070 format( 10x,'DATA IS AVERAGED IN SAMPLES OF IAV TURNS',t96,       &
     &'IAV =    ',i7 /10x,'START TURN NUMBER FOR THE ANALYSIS ',t93,    &
     &'NSTART =  ',i9 /10x,'THE ANALYSIS STOPS AFTER TURN NUMBER ',t94, &
     &'NSTOP =  ',i9, /10x,                                             &
     &'HORIZONTAL ANGLE-INTERVAL FOR STROBOSCOPING THE VERTICAL',       &
     &' PHASESPACE PROJECTION',t94,'DPHIX = ',d17.10 /10x,              &
     &'VERTICAL ANGLE-INTERVAL FOR STROBOSCOPING THE HORIZONTAL',       &
     &' PHASESPACE PROJECTION',t94,'DPHIY = ',d17.10 /10x,              &
     &'SWITCH FOR THE WEIGHTING OF THE LINEAR FIT FOR THE ' ,           &
     &'DISTANCE IN PHASESPACE ',t96,'IWG = ',i4 /10x,                   &
     &'INTEGER PART FOR THE HORIZONTAL TUNE ',t96,'QX0 = ' ,f16.10 /10x,&
     &'INTEGER PART FOR THE VERTICAL TUNE ',t96,'QY0 = ' ,f16.10 )
10080 format( 10x,'SWITCH FOR THE QX-VALUE CLOSE TO AN HALF-INTEGER' ,  &
     &' ( INT => 1 ; HALF-INT => 0 )',t95,'IVOX = ',i4 /10x,            &
     &'SWITCH FOR THE QY-VALUE CLOSE TO AN HALF-INTEGER' ,              &
     &' ( INT => 1 ; HALF-INT => 0 )',t95,'IVOY = ',i4 /10x,            &
     &'Q-VALUES ARE CHECKED FOR RESONANCES UP TO ORDER',t95, 'IRES = ', &
     &i4 /10x,'A RESONANCE IS CONSIDERED TO BE STRONG WHEN THE Q-VALUES'&
     &, ' ARE CLOSER TO IT THAN',t95,'DRES = ',d17.10 /10x,             &
     &'SWITCH FOR FFT-RANGE ( IFH=0 => 0-1 ; IFH=1 => 0-.5 ' ,          &
     &'; IFH=2 => .5-1 )',t96,'IFH = ',i4 /10x,                         &
     &'Q-PEAKS OF THE FFT ARE CONSIDERED IF THEY ARE LARGER THAN' ,t95, &
     &'DFFT = ',d17.10 )
10090 format( 10x,                                                      &
     &'SWITCH FOR PRINTING THE DISTANCE IN PHASE SPACE' ,t95,'IDIS = ', &
     &i4 /10x,'SWITCH FOR PRINTING THE COORDINATES',t95,'ICOW = ',i4 /10&
     &x,'SWITCH FOR PRINTING THE STROBOSCOPED PHASESPACES',t95,         &
     &'ISTW = ',i4 /10x,'SWITCH FOR PRINTING THE FFT-SIGNALS',t95,      &
     &'IFFW = ',i4 ,i4 )
10100 format( 10x,'EVERY ISKIP VALUE IS USED FOR THE ANALYSIS' ,t94,    &
     &'ISKIP = ',i4 /10x,                                               &
     &'SWITCH OF COURANT SYNDER TRANSFORMATION (ICONV = 1 => OFF)' ,t93,&
     &' ICONV = ',i4 /10x,                                              &
     &'SWITCH FOR READING MAD DATA ( IMAD = 1 => MAD-DATA ' ,           &
     &'WITH LYAPUNOV ANALYSIS )',t95,'IMAD = ',i4 /10x,                 &
     &'SCALING OF MOMENTUM', ' WITH LYAPUNOV ANALYSIS',t95,'CMA1 = ',f16&
     &.10 /10x,'SCALING OF PATH-LENGTH', ' WITH LYAPUNOV ANALYSIS',t95, &
     &'CMA2 = ',f16.10 /10x,                                            &
     &'SWITCH FOR PRINTING OF THE POSTPROCESSING OUTPUT' ,              &
     &' NPRINT = ( 0 => OFF ; 1 => ON) ',t93,'NPRINT = ',i4 /10x,       &
     &'NUMBER OF BINARY FILES TO BE PROCESSED', ' ( 90 - [90-NDAFI+1] )'&
     &,t94,'NDAFI = ',i4 //)
10110 format(/10x,'ANALYSING THE INCREASE OF THE DISTANCE IN PHASE-' ,  &
     &'SPACE'/10x,53('-')/ //12x,'TURNS',10x,'DISTANCE',20x,            &
     &'SLOPE          RESIDUAL' /10x,70('-'))
10120 format(10x,i7,6x,d17.10,2(2x,f18.11))
10130 format(10x,70('-')//)
10140 format(//10x,'AVERAGED PHASE-ADVANCE' /10x,22('-')/ /10x,         &
     &'X-PHASE :  ',f14.10,'   +/_ ',f14.10 /10x,'Y-PHASE :  ',f14.10,  &
     &'   +/_ ',f14.10/ /10x,'S-PHASE :  ',f14.10,'   +/_ ',f14.10/ /10 &
     &x,'START-QX : ',f14.10,'   CHANGE IN X : ',d17.10 /10x,           &
     &'START-QY : ',f14.10,'   CHANGE IN Y : ',d17.10 /10x,'START-QS : '&
     &,f14.10,'   CHANGE IN S : ',d17.10// /10x,                        &
     &'THE AVERAGED PHASE-ADVANCES ARE CLOSER THEN ',d11.4,' TO ' ,     &
     &'THE FOLLOWING RESONANCES UP TO ',i3,' ORDER'/10x,98('-')/ /10x,  &
     &'NX * QX   +   NY * QY   -      P      =      DELTA'/10x, 52('-'))
10150 format(/10x,'WARNING ! X-PHASE MIGHT NOT BE PRECISE'/)
10160 format(/10x,'WARNING ! Y-PHASE MIGHT NOT BE PRECISE'//)
10170 format(12x,i2,11x,i3,7x,f8.1,9x,d11.4)
10180 format(//10x,'Q-VALUES FROM FFT-ROUTINE' /10x,25('-')/ /10x,      &
     &'THE ANALYSIS WAS DONE WITH ',i7,' ENTRIES.'/ /10x,               &
     &'THE FOLLOWING Q-PEAKS ARE LARGER THEN ',f8.4,' PERCENT.'/ /10x,  &
     &'PLANE          Q-VALUE            SIZE [%]'/10x,43('-'))
10190 format(12x,'X',7x,f14.10,5x,f14.10)
10200 format(12x,'Y',7x,f14.10,5x,f14.10)
10210 format(//10x,'MAXIMUM PEAK'/ /10x,'HORIZONTAL PLANE :  ',f14.10   &
     &/10x,'VERTICAL PLANE   :  ',f14.10/ /10x,'START-QX : ',f14.10,    &
     &'   CHANGE IN X : ',d17.10 /10x,'START-QY : ',f14.10,             &
     &'   CHANGE IN Y : ',d17.10// /10x,                                &
     &'THE MAXIMUM Q-PEAKS ARE CLOSER THEN ',d11.4,' TO ' ,             &
     &'THE FOLLOWING RESONANCES UP TO ',i3,' ORDER'/10x,96('-')/ /10x,  &
     &'NX * QX   +   NY * QY   -      P      =      DELTA'/10x, 52('-'))
10220 format(////10x,'CALCULATION OF THE AVERAGED EMITTANCES' /10x,38(  &
     &'-')// 24x,'START-EMITTANCE           START-AMPLITUDE'// 10x,     &
     &'HORIZONTAL   ',f16.10,9x,f16.10/ 10x,'VERTICAL     ',f16.10,9x,  &
     &f16.10// 14x,'PLANE',10x,'EMITTANCE',16x,'SMEAR',12x,'MAXIMUM',11 &
     &x, 'MINIMUM'/28x,'[PI*MM*MRAD]',15x,'[%]',15x,'[%]',15x, '[%]'/10 &
     &x,86('-')/ 10x,'HORIZONTAL',6x,f16.10,3(6x,f12.6)/ 10x,'VERTICAL',&
     &8x,f16.10,3(6x,f12.6)/ 10x,'SUM',13x,f16.10,3(6x,f12.6)/10x,86('-'&
     &)//)
10230 format(//10x,'INVARIANTS OF THE 4-DIMENSIONAL PHASE-SPACE' /10x,43&
     &('-')//)
10240 format(/10x,'WARNING ! CALCULATION OF THE HORIZONTAL INVARIANT' , &
     &' MIGHT NOT BE PRECISE'/10x,'ONLY ',i5,' ENTRIES COMPARED TO ' ,  &
     &i5,' ANGLE-INTERVALS !'/10x,'INCREASE THE VERTICAL ANGLE-INT' ,   &
     &'ERVAL <DPHIY>'/)
10250 format(/10x,'WARNING ! CALCULATION OF THE VERTICAL INVARIANT' ,   &
     &' MIGHT NOT BE PRECISE'/10x,'ONLY ',i5,' ENTRIES COMPARED TO ' ,  &
     &i5,' ANGLE-INTERVALS !'/10x,'INCREASE THE HORIZONTAL ANGLE-INT' , &
     &'ERVAL <DPHIX>'/)
10260 format(/10x,'THERE ARE ',i5,' ENTRIES FOR THE CALCULATION OF' ,   &
     &' THE HORIZONTAL INVARIANT GROUPED IN ',i5,' ANGLE-INTERVALS' /10 &
     &x,'--------- ',i5,' ENTRIES ----------------------' ,             &
     &'---- VERTICAL   -------------------- ',i5,' ANGLE-INTERVALS' //10&
     &x,'IF THE MOTION IS CLOSE TO FIXPOINTS THE NUMBER OF THOSE' ,     &
     &' ANGLE-INTERVALS WILL BE ONLY A SMALL FRACTION'/10x,             &
     &'OF THE TOTAL NUMBER OF ',i5,' INTERVALS.'/ /25x,                 &
     &'PERCENTAGE OF OCCUPIED INTERVALS     NUMBER OF ENTRIES ' ,       &
     &'PER OCCUPIED INTERVAL' /10x,'HORIZONTAL',16x,f10.6,30x,f12.6 /10 &
     &x,'VERTICAL  ',16x,f10.6,30x,f12.6/ //10x,                        &
     &'THE CALCULATED INVARIANTS ARE IN UNITS OF [PI*MM*MRAD]'/ /10x,   &
     &'HORIZONTAL',10x,f16.10 /10x,'VERTICAL  ',10x,f16.10//)
10270 format(////10x,'LINEARLY DECOUPLED INVARIANTS' /10x,35('-')/ 10x, &
     &'INITIAL EMITTANCE MODE I   :',f16.10/ 10x,                       &
     &'INITIAL EMITTANCE MODE II  :',f16.10/ 10x,                       &
     &'INITIAL EMITTANCE MODE III :',f16.10/ 10x,                       &
     &'INITIAL ANGLE     MODE I   :',f16.10/ 10x,                       &
     &'INITIAL ANGLE     MODE II  :',f16.10/ 10x,                       &
     &'INITIAL ANGLE     MODE III :',f16.10/ /10x,35('-')// 14x,'PLANE',&
     &10x,'EMITTANCE',14x,'4D-SMEAR',11x,'MAXIMUM',11x, 'MINIMUM'/28x,  &
     &'[PI*MM*MRAD]',15x,'[%]',15x,'[%]',15x, '[%]'/10x,86('-')/ 10x,   &
     &'HORIZONTAL',6x,f16.10,3(6x,f12.6)/ 10x,'VERTICAL',8x,f16.10,3    &
     &(6x,f12.6)/ 10x,'SUM',13x,f16.10,3(6x,f12.6)/10x,86('-')//)
10280 format(/10x,'Postprocessing took ',f12.3,' second(s)',            &
     &' of Execution Time'//131('-')//)
10290 format(//10x,'** ERROR ** ----- TRANSFORMATION MATRIX SINGULAR ' ,&
     &'(FILE : ',i2,') -----'//)
10300 format(//10x,'** ERROR ** ----- FILE :',i2,' WITH TRACKING ' ,    &
     &'DATA EMPTY OR CORRUPTED-----'/10x,'PROBLEM : ',a80//)
10310 format(//10x,'** ERROR ** ----- WEIGHTING OF DISTANCE IN PHASE' , &
     &' SPACE (FILE : ',i2,') NOT POSSIBLE-----'//)
10320 format(//10x,'** ERROR ** ----- INPUT DATA CORRUPTED' ,' (FILE : '&
     &,i2,') -----'//)

end subroutine postpr

subroutine fft(ar,ai,m,n)
!---------------------------------------------------------------------
!      A(N) IS A COMPLEX ARRAY. THE INPUT IS A(N)=(X(N),0.0), WHERE
!      X(N) IS THE SEQUENCE TO FOURIER ANALYSE.
!      N=2**M.
!      THIS ROUTINE ONLY WORKS FOR N EQUAL TO A POWER OF TWO.
!      AFTER RETURN A(N)=(..,..) CONTAINS THE COEFFICIENTS OF THE FFT.
!      TO COMPUTE POWER SPECTRA DO   ...=ABS(A(N))
!
!      WRITTEN BY : RUI DILAO
!
!---------------------------------------------------------------------
      use mathlib_bouncer
      use numerical_constants
      implicit none
      integer i,ip,j,k,l,le,le1,m,n,nm1,nv2
      real(kind=fPrec) ar,ai,tr,ti,ui,ur,uur,wr,wi
      dimension ar(n),ai(n)
      save
!-----------------------------------------------------------------------
      n=2**m
      nv2=n/2
      nm1=n-1
      j=1

      do i=1,nm1
        if(i.gt.j)goto 10
        tr=ar(j)
        ti=ai(j)
        ar(j)=ar(i)
        ai(j)=ai(i)
        ar(i)=tr
        ai(i)=ti
   10   k=nv2
   20   if(k.ge.j)goto 30
        j=j-k
        k=k/2
        goto 20
        j=j+k
   30   continue
      end do

      do l=1,m
        le=2**l
        le1=le/2
        ur=one
        ui=zero
        wr=cos_mb(pi/real(le1,fPrec))                                          !hr06
        wi=-one*sin_mb(pi/real(le1,fPrec))                                     !hr06
        do j=1,le1
          do i=j,n,le
            ip=i+le1
            tr=ar(ip)*ur-ai(ip)*ui
            ti=ar(ip)*ui+ai(ip)*ur
            ar(ip)=ar(i)-tr
            ai(ip)=ai(i)-ti
            ar(i)=ar(i)+tr
            ai(i)=ai(i)+ti
          end do

          uur=ur*wr-ui*wi
          ui=ur*wi+ui*wr
          ur=uur
        end do
      end do

      do i=1,n
        ar(i)=(ar(i)/real(n,fPrec))*2                                          !hr06
        ai(i)=(ai(i)/real(n,fPrec))*2                                          !hr06
      end do

      return
end subroutine fft

subroutine caconv(a,b,c)
      use mathlib_bouncer
      use numerical_constants
      implicit none
      real(kind=fPrec) a,b,c
      save
!---------------------------------------------------------------------
      if(abs(b).gt.pieni.or.abs(c).gt.pieni) then
        a=atan2_mb(b,c)
      else
        a=zero
      endif
      return
end subroutine caconv

subroutine cphase(k,a,b,c,d,i,j,ie)
      use mathlib_bouncer
      use numerical_constants
      use parpro
      implicit none
      integer i,ie,j,k
      real(kind=fPrec) a,b,c,d,f,tpi
      save
!---------------------------------------------------------------------
      tpi=eight*atan_mb(one)                                               !hr06
      if(abs(b).gt.pieni.or.abs(c).gt.pieni) then
        f=atan2_mb(b,c)
        ie=ie+1
        phase(k,ie)=f/tpi+d
        if(i.ne.1.and.-f.gt.pieni) phase(k,ie)=phase(k,ie)+one
        a=a+phase(k,ie)
      else
        j=1
      endif
      return
end subroutine cphase

subroutine cinvar(a,b,c,d,j,e,xinv,invx)
      use mathlib_bouncer
      use numerical_constants
      use parpro
      implicit none
      integer i,invx,j
      real(kind=fPrec) a,b,c,d,e,xinv
      dimension xinv(ninv),invx(ninv)
      save
!---------------------------------------------------------------------
      if(abs(a).le.b) then
        do 10 i=1,ninv
          if((c.ge.zero.and.c.ge.dani(i).and.c.lt.dani(i+1)).or. (c.lt.zero.and.d.ge.dani(i).and.d.lt.dani(i+1))) then
            j=j+1
            if(abs(xinv(i)).le.pieni) then
              xinv(i)=e
              invx(i)=1
            else
              xinv(i)=xinv(i)+e
              invx(i)=invx(i)+1
            endif
          endif
   10   continue
      endif
      return
end subroutine cinvar

subroutine sinpro(a,b,c,d,e)
      use mathlib_bouncer
      use numerical_constants
      implicit none
      real(kind=fPrec) a,b,c,d,e
      save
!---------------------------------------------------------------------
      if(abs(a).gt.pieni) then
        if(c.gt.pieni.and.b.gt.pieni) then
          c=(sqrt(c/b)/a)*c1e2                                          !hr06
        else
          c=zero
        endif
        d=((d-a)/a)*c1e2                                                !hr06
        e=((e-a)/a)*c1e2                                                !hr06
      else
        c=zero
        d=zero
        e=zero
      endif
      return
end subroutine sinpro


subroutine join
      use mathlib_bouncer
      use numerical_constants
      use crcoall
      use parpro
      use mod_common, only : dpscor,sigcor,icode,idam,its6d
      implicit none

      integer i,j,ia,idummy,ierro,ifipa,ihalf,ilapa,ipa,ipa1,itopa,numl
      real(kind=fPrec) alf0,bet0,clo,clop,di0,dip0,dps,dummy,e0,qwc,sigm,ta,x,y
      character(len=8) cdate,ctime,progrm ! Note: Keep in sync with maincr
      character(len=80) sixtit,commen     ! Note: Keep in sync with mod_common
                                          ! DANGER: If the len changes, CRCHECK will break.
      dimension bet0(2),alf0(2),ta(6,6)
      dimension qwc(3),clo(3),clop(3)
      dimension x(mpa,2),y(mpa,2),sigm(mpa),dps(mpa)
      dimension di0(2),dip0(2)

      !The fort.90 file is always with real64, so we need some temps to read it
      ! For the header:
      real(kind=real64) qwc_tmp(3), clo_tmp(3), clop_tmp(3)
      real(kind=real64) di0_tmp(2), dip0_tmp(2)
      real(kind=real64) ta_tmp(6,6)
      real(kind=real64) dummy64, dam_tmp
      real(kind=real64) sigcor64, dpscor64, zero64

      !For the actual tracking data
      real(kind=real64) x_tmp(2,2),y_tmp(2,2), sigm_tmp(2),dps_tmp(2),e0_tmp

      save
!-----------------------------------------------------------------------

      sigcor64 = real(one, real64)
      dpscor64 = real(one, real64)
      zero64   = real(zero,real64)

      read(90,iostat=ierro) &
           sixtit,commen,cdate,ctime,progrm, &
           ifipa,ilapa,itopa,icode,numl, &
           qwc_tmp(1),qwc_tmp(2),qwc_tmp(3), &
           clo_tmp(1),clop_tmp(1),clo_tmp(2),clop_tmp(2), &
           clo_tmp(3),clop_tmp(3), &
           di0_tmp(1),dip0_tmp(1),di0_tmp(2),dip0_tmp(2), &
           dummy64,dummy64, &
           ta_tmp(1,1),ta_tmp(1,2),ta_tmp(1,3), &
           ta_tmp(1,4),ta_tmp(1,5),ta_tmp(1,6), &
           ta_tmp(2,1),ta_tmp(2,2),ta_tmp(2,3), &
           ta_tmp(2,4),ta_tmp(2,5),ta_tmp(2,6), &
           ta_tmp(3,1),ta_tmp(3,2),ta_tmp(3,3), &
           ta_tmp(3,4),ta_tmp(3,5),ta_tmp(3,6), &
           ta_tmp(4,1),ta_tmp(4,2),ta_tmp(4,3), &
           ta_tmp(4,4),ta_tmp(4,5),ta_tmp(4,6), &
           ta_tmp(5,1),ta_tmp(5,2),ta_tmp(5,3), &
           ta_tmp(5,4),ta_tmp(5,5),ta_tmp(5,6), &
           ta_tmp(6,1),ta_tmp(6,2),ta_tmp(6,3), &
           ta_tmp(6,4),ta_tmp(6,5),ta_tmp(6,6)

      if (IS_IOSTAT_END(ierro)) then
         write(lout,10000) 90
10000    format(//10x,'** ERROR IN JOIN** ----- INPUT DATA EMPTY (FILE : ',i2,') -----'//)
         return
      endif

      !Convert it to the current working precission
      do i=1,3
         qwc(i)  = real(qwc_tmp (i), fPrec)
         clo(i)  = real(clo_tmp (i), fPrec)
         clop(i) = real(clop_tmp(i), fPrec)
      end do

      do i=1,2
         di0(i)  = real(di0_tmp (i), fPrec)
         dip0(i) = real(dip0_tmp(i), fPrec)
      enddo

      do i=1,6
         do j=1,6
            ta(j,i) = real(ta_tmp(j,i), fPrec)
         end do
      end do

      if(ierro.gt.0) then
         write(lout,10010) 90,ierro
         return
      endif
!-----------------------------------------------------------------------
!  OPTICAL PARAMETERS AT THE STARTING POINT (BECAUSE WHY NOT?!?)
!-----------------------------------------------------------------------
      bet0(1)=ta(1,1)*ta(1,1)+ta(1,2)*ta(1,2)
      alf0(1)=-(ta(1,1)*ta(2,1)+ta(1,2)*ta(2,2))
      bet0(2)=ta(3,3)*ta(3,3)+ta(3,4)*ta(3,4)
      alf0(2)=-(ta(3,3)*ta(4,3)+ta(3,4)*ta(4,4))

      rewind 90

      ihalf=itopa/2

      if(icode.eq.1.or.icode.eq.2.or.icode.eq.4) idam=1
      if(icode.eq.3.or.icode.eq.5.or.icode.eq.6) idam=2
      if(icode.eq.7) idam=3

      !Read pairs of binary files; 90 and 90-ihalf, 89 and 89-ihalf, etc.
      do i=1,ihalf

         !Skip the headers
         read(91-i,iostat=ierro)
         if (IS_IOSTAT_END(ierro)) then
            exit
         endif
         if(ierro.gt.0) then
            write(lout,10010) 91-i,ierro
            exit
         endif

         read(91-i-ihalf,iostat=ierro)
         if (IS_IOSTAT_END(ierro)) then
            exit
         endif
         if(ierro.gt.0) then
            write(lout,10010) 91-i-ihalf,ierro
            exit
         endif

         !Loop to read the actual tracking data and write it to fort.90
         do

            !First file
            read(91-i,iostat=ierro) &
                 ia,ipa,dummy64, &
                 x_tmp(1,1),y_tmp(1,1),x_tmp(1,2),y_tmp(1,2), &
                 sigm_tmp(1),dps_tmp(1),e0_tmp
            if (IS_IOSTAT_END(ierro)) then
               exit
            endif
            !Convert it to the current working precission
            x(1,1)     = real(x_tmp(1,1), fPrec)
            y(1,1)     = real(y_tmp(1,1), fPrec)
            x(1,2)     = real(x_tmp(1,2), fPrec)
            y(1,2)     = real(y_tmp(1,2), fPrec)
            sigm(1)    = real(sigm_tmp(1),fPrec)
            dps(1)     = real(dps_tmp(1), fPrec)
            e0         = real(e0,         fPrec)

            if(ierro.gt.0) then
               write(lout,10010) 91-i,ierro
               exit
            endif

            x(1,1)=x(1,1)*c1e3
            y(1,1)=y(1,1)*c1e3
            x(1,2)=x(1,2)*c1e3
            y(1,2)=y(1,2)*c1e3
            sigm(1)=sigm(1)*c1e3
            e0=e0*c1e3

            !Second file
            read(91-i-ihalf,iostat=ierro) &
                 idummy,idummy,dummy64, &
                 x_tmp(2,1),y_tmp(2,1),x_tmp(2,2),y_tmp(2,2), &
                 sigm_tmp(2),dps_tmp(2)
            if (IS_IOSTAT_END(ierro)) then
               exit
            endif
            !Convert it to the current working precission
            x(2,1)     = real(x_tmp(2,1), fPrec)
            y(2,1)     = real(y_tmp(2,1), fPrec)
            x(2,2)     = real(x_tmp(2,2), fPrec)
            y(2,2)     = real(y_tmp(2,2), fPrec)
            sigm(2)    = real(sigm_tmp(2),fPrec)
            dps(2)     = real(dps_tmp(2), fPrec)
            e0         = real(e0,         fPrec)

            if(ierro.gt.0) then
               write(lout,10010) 91-i-ihalf,ierro
               exit
            endif
            x(2,1)=x(2,1)*c1e3
            y(2,1)=y(2,1)*c1e3
            x(2,2)=x(2,2)*c1e3
            y(2,2)=y(2,2)*c1e3
            sigm(2)=sigm(2)*c1e3

            !OK, we have read a pair, now write it back to fort.90
            ! (possible bug: Should this have been back to fort.(91-i)?)

            !Convert it to real64 before writing
            x_tmp(1,1) = real(x(1,1), real64)
            y_tmp(1,1) = real(y(1,1), real64)
            x_tmp(1,2) = real(x(1,2), real64)
            y_tmp(1,2) = real(x(1,2), real64)

            sigm_tmp(1) = real(sigm(1), real64)
            dps_tmp (1) = real(dps (1), real64)
            e0_tmp      = real(e0,      real64)

            x_tmp(2,1) = real(x(2,1), real64)
            y_tmp(2,1) = real(y(2,1), real64)
            x_tmp(2,2) = real(x(2,2), real64)
            y_tmp(2,2) = real(x(2,2), real64)

            sigm_tmp(2) = real(sigm(2), real64)
            dps_tmp (2) = real(dps (2), real64)

            write(90,iostat=ierro) ia,ipa,dummy64, &
                 x_tmp(1,1),y_tmp(1,1),x_tmp(1,2),y_tmp(1,2), &
                 sigm_tmp(1),dps_tmp(1),e0_tmp, &
                 ipa+1,dummy64,x_tmp(2,1),y_tmp(2,1),x_tmp(2,2),y_tmp(2,2), &
                 sigm_tmp(2),dps_tmp(2),e0_tmp

            if(ierro.ne.0) then
               write(lout,10010) 90,ierro
               exit
            endif

            !Loop back and read the next line...
         end do
         !Ok, we're done reading the pair of files.
         rewind 91-i
         rewind 91-i-ihalf

         !Erase the bottommost file
         write(91-i-ihalf,iostat=ierro)
         if(ierro.ne.0) then
            write(lout,10010) 91-i-ihalf,ierro
         endif
         rewind 90

         !Write "MADTOSIX" into the 'progrm' field of the header of the topmost file
         read(91-i,iostat=ierro) &
           sixtit,commen,cdate,ctime,progrm, &
           ifipa,ilapa,itopa,icode,numl, &
           qwc_tmp(1),qwc_tmp(2),qwc_tmp(3), &
           clo_tmp(1),clop_tmp(1),clo_tmp(2),clop_tmp(2), &
           clo_tmp(3),clop_tmp(3), &
           di0_tmp(1),dip0_tmp(1),di0_tmp(2),dip0_tmp(2), &
           dummy64,dummy64, &
           ta_tmp(1,1),ta_tmp(1,2),ta_tmp(1,3), &
           ta_tmp(1,4),ta_tmp(1,5),ta_tmp(1,6), &
           ta_tmp(2,1),ta_tmp(2,2),ta_tmp(2,3), &
           ta_tmp(2,4),ta_tmp(2,5),ta_tmp(2,6), &
           ta_tmp(3,1),ta_tmp(3,2),ta_tmp(3,3), &
           ta_tmp(3,4),ta_tmp(3,5),ta_tmp(3,6), &
           ta_tmp(4,1),ta_tmp(4,2),ta_tmp(4,3), &
           ta_tmp(4,4),ta_tmp(4,5),ta_tmp(4,6), &
           ta_tmp(5,1),ta_tmp(5,2),ta_tmp(5,3), &
           ta_tmp(5,4),ta_tmp(5,5),ta_tmp(5,6), &
           ta_tmp(6,1),ta_tmp(6,2),ta_tmp(6,3), &
           ta_tmp(6,4),ta_tmp(6,5),ta_tmp(6,6)

         if(ierro.gt.0) then
            write(lout,10010) 91-i,ierro
            goto 40
         endif

         rewind 91-i
         progrm='MADTOSIX'

         write(91-i,iostat=ierro) &
           sixtit,commen,cdate,ctime,progrm, &
           ifipa,ilapa,itopa,icode,numl, &
           qwc_tmp(1),qwc_tmp(2),qwc_tmp(3), &
           clo_tmp(1),clop_tmp(1),clo_tmp(2),clop_tmp(2), &
           clo_tmp(3),clop_tmp(3), &
           di0_tmp(1),dip0_tmp(1),di0_tmp(2),dip0_tmp(2), &
           dummy64,dummy64, &
           ta_tmp(1,1),ta_tmp(1,2),ta_tmp(1,3), &
           ta_tmp(1,4),ta_tmp(1,5),ta_tmp(1,6), &
           ta_tmp(2,1),ta_tmp(2,2),ta_tmp(2,3), &
           ta_tmp(2,4),ta_tmp(2,5),ta_tmp(2,6), &
           ta_tmp(3,1),ta_tmp(3,2),ta_tmp(3,3), &
           ta_tmp(3,4),ta_tmp(3,5),ta_tmp(3,6), &
           ta_tmp(4,1),ta_tmp(4,2),ta_tmp(4,3), &
           ta_tmp(4,4),ta_tmp(4,5),ta_tmp(4,6), &
           ta_tmp(5,1),ta_tmp(5,2),ta_tmp(5,3), &
           ta_tmp(5,4),ta_tmp(5,5),ta_tmp(5,6), &
           ta_tmp(6,1),ta_tmp(6,2),ta_tmp(6,3), &
           ta_tmp(6,4),ta_tmp(6,5),ta_tmp(6,6), &
           zero64,zero64,zero64,zero64, &
           sigcor64,dpscor64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64, &
           zero64,zero64,zero64,zero64

         if(ierro.ne.0) then
            write(lout,10010) 91-i,ierro
            goto 40
         endif

         !copy the particle records from fort.90 to fort.91-i
         do
            read(90,iostat=ierro) &
                 ia,ipa,dam_tmp, &
                 x_tmp(1,1),y_tmp(1,1),x_tmp(1,2),y_tmp(1,2), &
                 sigm_tmp(1),dps_tmp(1),e0_tmp, &
                 ipa1,dam_tmp,x_tmp(2,1),y_tmp(2,1),x_tmp(2,2),y_tmp(2,2),    &
                 &sigm_tmp(2),dps_tmp(2),e0_tmp
            if (IS_IOSTAT_END(ierro)) then
               exit
            endif
            if(ierro.gt.0) then
               write(lout,10010) 90,ierro
               exit
            endif

            write(91-i,iostat=ierro) &
                 ia,ipa,dam_tmp, &
                 x_tmp(1,1),y_tmp(1,1),x_tmp(1,2),y_tmp(1,2), &
                 sigm_tmp(1),dps_tmp(1),e0_tmp, &
                 ipa1,dam_tmp,x_tmp(2,1),y_tmp(2,1),x_tmp(2,2),y_tmp(2,2),    &
                 &sigm_tmp(2),dps_tmp(2),e0_tmp
            if(ierro.ne.0) then
               write(lout,10010) 91-i,ierro
               exit
            endif
         end do

         !Skip the copying of particle records
 40      continue
         rewind 90
         rewind 91-i

      end do ! END of loop over file pairs

      return

10010 format(//10x,'** ERROR IN JOIN** ----- PROBLEMS WITH DATA FILE : ',i2,' ----- ERROR CODE : ',i10//)
end subroutine join


      subroutine writebin_header(ia_p1,ia_p2,fileunit_in, ierro_wbh,    &
     &     cdate,ctime,progrm)
!-------------------------------------------------------------------------
!     Subroutine for writing the header of the binary files (fort.90 etc.)
!     Always converting to real64 before writing to disk
!-----------------------------------------------------------------------
!     K. SJOBAK, October 2017
!-----------------------------------------------------------------------
      use numerical_constants
      use, intrinsic :: iso_fortran_env, only : real64
      use parpro
      use mod_common
      use mod_common_main
      implicit none

      integer, intent(in) :: ia_p1, ia_p2, fileunit_in
      integer, intent(inout) :: ierro_wbh

      character(len=8) cdate,ctime,progrm !Note: Keep in sync with maincr
                                          !DANGER: If the len changes, CRCHECK will break.

      integer i,j

      real(kind=real64) qwcs_tmp(3), clo6v_tmp(3), clop6v_tmp(3)
      real(kind=real64) di0xs_tmp, dip0xs_tmp, di0zs_tmp,dip0zs_tmp
      real(kind=real64) tas_tmp(6,6)
      real(kind=real64) mmac_tmp,nms_tmp,izu0_tmp,numlr_tmp,            &
     &     sigcor_tmp,dpscor_tmp

      real(kind=real64) zero64,one64
      parameter(zero64 = 0.0_real64)
      parameter(one64  = 1.0_real64)

      !Convert from whatever precission is used internally to real64,
      ! which is what should go in the output file
      do i=1,3
         qwcs_tmp  (i) = real(qwcs  (i), real64)
         clo6v_tmp (i) = real(clo6v (i), real64)
         clop6v_tmp(i) = real(clop6v(i), real64)
      enddo

      di0xs_tmp  = real(di0xs, real64)
      dip0xs_tmp = real(dip0xs, real64)
      di0zs_tmp  = real(di0zs, real64)
      dip0zs_tmp = real(dip0zs, real64)

      do i=1,6
         do j=1,6
            tas_tmp(j,i) = real(tas(j,i), real64)
         enddo
      enddo

      mmac_tmp   = 1.0_real64
      nms_tmp    = real(nms(ia_p1), real64)
      izu0_tmp   = real(izu0,       real64)
      numlr_tmp  = real(numlr,      real64)
      sigcor_tmp = real(sigcor,     real64)
      dpscor_tmp = real(dpscor,     real64)

      ! DANGER: IF THE LENGTH OR NUMBER OF THESE FIELDS CHANGE,
      ! CRCHECK WON'T WORK. SEE HOW VARIABLES HBUFF/TBUFF ARE USED.
      ! WE ALSO ASSUME THAT THE INTEGERS ARE ALWAYS 32BIT...

      write(fileunit_in,iostat=ierro_wbh)                               &
     &     sixtit,commen,cdate,ctime,progrm,                            &
     &     ia_p1,ia_p2, napx, icode,numl,                               &
     &     qwcs_tmp(1),qwcs_tmp(2),qwcs_tmp(3),                         &
     &     clo6v_tmp(1),clop6v_tmp(1),clo6v_tmp(2),clop6v_tmp(2),       &
     &     clo6v_tmp(3),clop6v_tmp(3),                                  &
     &     di0xs_tmp,dip0xs_tmp,di0zs_tmp,dip0zs_tmp,                   &
     &     zero64,one64,                                                &
     &     tas_tmp(1,1),tas_tmp(1,2),tas_tmp(1,3),                      &
     &     tas_tmp(1,4),tas_tmp(1,5),tas_tmp(1,6),                      &
     &     tas_tmp(2,1),tas_tmp(2,2),tas_tmp(2,3),                      &
     &     tas_tmp(2,4),tas_tmp(2,5),tas_tmp(2,6),                      &
     &     tas_tmp(3,1),tas_tmp(3,2),tas_tmp(3,3),                      &
     &     tas_tmp(3,4),tas_tmp(3,5),tas_tmp(3,6),                      &
     &     tas_tmp(4,1),tas_tmp(4,2),tas_tmp(4,3),                      &
     &     tas_tmp(4,4),tas_tmp(4,5),tas_tmp(4,6),                      &
     &     tas_tmp(5,1),tas_tmp(5,2),tas_tmp(5,3),                      &
     &     tas_tmp(5,4),tas_tmp(5,5),tas_tmp(5,6),                      &
     &     tas_tmp(6,1),tas_tmp(6,2),tas_tmp(6,3),                      &
     &     tas_tmp(6,4),tas_tmp(6,5),tas_tmp(6,6),                      &
     &     mmac_tmp,nms_tmp,izu0_tmp,numlr_tmp,                         &
     &     sigcor_tmp,dpscor_tmp,                                       &
     &     zero64,zero64,zero64,zero64,zero64,zero64,zero64,zero64,     &
     &     zero64,zero64,zero64,zero64,zero64,zero64,zero64,zero64,     &
     &     zero64,zero64,zero64,zero64,zero64,zero64,zero64,zero64,     &
     &     zero64,zero64,zero64,zero64,zero64,zero64,zero64,zero64,     &
     &     zero64,zero64,zero64,zero64,zero64,zero64,zero64,zero64,     &
     &     zero64,zero64,zero64,zero64

      end subroutine writebin_header

subroutine writebin(nthinerr)
!-------------------------------------------------------------------------
!     Subroutine for writing the the binary files (fort.90 etc.)
!     Always converting to real64 before writing to disk
!-------------------------------------------------------------------------
!     F. SCHMIDT, 3 February 1999
!     K. SJOBAK,    October  2017
!     V.K. Berglyd Olsen, April 2019
!-------------------------------------------------------------------------
      use numerical_constants
      use, intrinsic :: iso_fortran_env, only : real64
      use crcoall
      use parpro
      use mod_common
      use mod_common_main
      use mod_commons
      use mod_common_track
      use mod_common_da
      use mod_settings
#ifdef CR
      use checkpoint_restart
#endif
      implicit none

      integer ia,ia2,ie,nthinerr
#ifdef BOINC
      integer timech
#endif

      real(kind=real64) dam_tmp, xv_tmp(2,2),yv_tmp(2,2),sigmv_tmp(2),dpsv_tmp(2),e0_tmp

      save
!-----------------------------------------------------------------------
#ifdef CR
      if(cr_restart) then
        write(crlog,"(a)") "WRITEBIN> Restarting, so not writing records"
        flush(crlog)
        return
      end if
#endif
         do ia=1,napx-1
!GRD
!     PSTOP=true -> particle lost,
!     partID(ia)=particle ID that is not changing
!     (the particle arrays are compressed to remove lost particles)
!     In the case of no lost particles, all partID(i)=i for 1..npart
            if(.not.pstop(partID(ia)).and..not.pstop(partID(ia)+1).and. &
     &           (mod(partID(ia),2).ne.0)) then !Skip odd particle IDs

               ia2=(partID(ia)+1)/2 !File ID for non-STF & binrecs
               ie=ia+1              !ia = Particle ID 1, ie = Particle ID 2

               if(ntwin.ne.2) then !Write particle partID(ia) only
                  dam_tmp      = real(dam(ia),   real64)
                  xv_tmp(1,1)  = real(xv1(ia),  real64)
                  yv_tmp(1,1)  = real(yv1(ia),  real64)
                  xv_tmp(2,1)  = real(xv2(ia),  real64)
                  yv_tmp(2,1)  = real(yv2(ia),  real64)
                  sigmv_tmp(1) = real(sigmv(ia), real64)
                  dpsv_tmp(1)  = real(dpsv(ia),  real64)
                  e0_tmp       = real(e0,        real64)

                  ! DANGER: IF THE LENGTH OR NUMBER OF THESE FIELDS CHANGE,
                  ! CRCHECK WON'T WORK. SEE HOW VARIABLES HBUFF/TBUFF ARE USED.
                  ! WE ALSO ASSUME THAT THE INTEGERS ARE ALWAYS 32BIT...

#ifndef STF
                  write(91-ia2,iostat=ierro)                            &
#else
                  write(90,iostat=ierro)                                &
#endif
     &               numx,partID(ia),dam_tmp,                           &
     &               xv_tmp(1,1),yv_tmp(1,1),                           &
     &               xv_tmp(2,1),yv_tmp(2,1),                           &
     &               sigmv_tmp(1),dpsv_tmp(1),e0_tmp
#ifndef STF
                  flush(91-ia2)
#else
                  flush(90)
#endif
#ifdef CR
                  binrecs(ia2)=binrecs(ia2)+1
#endif

               else !Write both particles partID(ia) and partID(ia)+1
                    ! Note that dam(ia) (distance in angular phase space)
                    ! is written twice.
                  dam_tmp      = real(dam(ia),   real64)

                  xv_tmp(1,1)  = real(xv1(ia),  real64)
                  yv_tmp(1,1)  = real(yv1(ia),  real64)
                  xv_tmp(2,1)  = real(xv2(ia),  real64)
                  yv_tmp(2,1)  = real(yv2(ia),  real64)
                  sigmv_tmp(1) = real(sigmv(ia), real64)
                  dpsv_tmp(1)  = real(dpsv(ia),  real64)

                  xv_tmp(1,2)  = real(xv1(ie),  real64)
                  yv_tmp(1,2)  = real(yv1(ie),  real64)
                  xv_tmp(2,2)  = real(xv2(ie),  real64)
                  yv_tmp(2,2)  = real(yv2(ie),  real64)
                  sigmv_tmp(2) = real(sigmv(ie), real64)
                  dpsv_tmp(2)  = real(dpsv(ie),  real64)

                  e0_tmp       = real(e0,        real64)

                  ! DANGER: IF THE LENGTH OR NUMBER OF THESE FIELDS CHANGE,
                  ! CRCHECK WON'T WORK. SEE HOW VARIABLES HBUFF/TBUFF ARE USED.
                  ! WE ALSO ASSUME THAT THE INTEGERS ARE ALWAYS 32BIT...
#ifndef STF
                  write(91-ia2,iostat=ierro)                            &
#else
                  write(90,iostat=ierro)                                &
#endif
     &               numx,partID(ia),dam_tmp,                           &
     &               xv_tmp(1,1),yv_tmp(1,1),                           &
     &               xv_tmp(2,1),yv_tmp(2,1),                           &
     &               sigmv_tmp(1),dpsv_tmp(1),e0_tmp,                   &
     &               partID(ia)+1,dam_tmp,                              &
     &               xv_tmp(1,2),yv_tmp(1,2),                           &
     &               xv_tmp(2,2),yv_tmp(2,2),                           &
     &               sigmv_tmp(2),dpsv_tmp(2),e0_tmp
#ifndef STF
                  flush(91-ia2)
#else
                  flush(90)
#endif

#ifdef CR
                  binrecs(ia2)=binrecs(ia2)+1
#endif
               endif
               if(ierro.ne.0) then
                  write(lout,"(2(a,i0))") "WRITEBIN> ERROR Problem writing to file #: ",91-ia2,", error code: ",ierro
#ifdef CR
                  flush(lout)
#else
                  flush(12)
#endif
                  nthinerr=3000
                  return
               endif
            endif
         end do !END "do 10 ia=1,napx-1"
#ifdef CR
      if (lhc.ne.9) then
         binrec=binrec+1
      endif
#endif

end subroutine writebin

! These routines have been moved from sixtrack.s90 because they are only used here. VKBO 2018-05-25
subroutine lfitwd(x,y,w,l,key,a,b,e)
!-----------------------------------------------------------------------
!
!     TO PERFORM A WEIGHTED STRAIGHT LINE FIT
!
!     FOR FORMULAE USED SEE MENZEL, FORMULAS OF PHYSICS P.116
!
!     FIT IS OF Y=AX+B , WITH S**2 ESTIMATOR E. WEIGHTS ARE IN W.
!     IF KEY=0, POINTS WITH Y=0 ARE IGNORED
!     L IS NO. OF POINTS
!
!-----------------------------------------------------------------------
!Eric made DOUBLE PRECISION
  implicit none
  integer icnt,j,key,l
  real(kind=fPrec) a,b,e,x,y,w
  real(kind=fPrec) w2,w2x,w2x2,w2xy,w2y,w2y2,ww,wwf,wwfi
  dimension x(l),y(l),w(l)                                           !hr07
  save
!-----------------------------------------------------------------------
!
!     CALCULATE SUMS
!
!-----------------------------------------------------------------------
  if(l.le.1) goto 1
  w2=0.
  w2x=0.
  w2y=0.
  w2xy=0.
  w2x2=0.
  w2y2=0.
  icnt=0
  do 2 j=1,l
  if(y(j).eq.0..and.key.eq.0) goto 2
  ww=w(j)**2                                                         !hr07
  w2=ww+w2
  wwf=ww*x(j)
  w2x=wwf+w2x
  w2x2=wwf*x(j)+w2x2
  w2xy=wwf*y(j)+w2xy
  wwfi=ww*y(j)
  w2y=wwfi+w2y
  w2y2=wwfi*y(j)+w2y2
  icnt=icnt+1
2 continue
!
!     FIT PARAMETERS
  a=(w2xy-(w2x*w2y)/w2)/(w2x2-w2x**2/w2)
  b=(w2y-a*w2x)/w2
  if(icnt.le.2) goto 3
!Eric
  e=((w2y2-w2y**2/w2)-(w2xy-(w2x*w2y)/w2)**2/(w2x2-w2x**2/w2))/     &!hr07
  &real(icnt-2,fPrec)
  goto 4
!
!     ISUFFICIENT POINTS
1 a=0.
  b=0.
3 e=0.
4 return
end subroutine lfitwd

subroutine lfitd(x,y,l,key,a,b,e)
!-----------------------------------------------------------------------
!
!     TO FIT A STRAIGHT LINE    Y=A*X+B    TO L POINTS WITH ERROR E
!     SEE MENZEL , FORMULAS OF PHYSICS P.116
!     POINTS WITH Y=0 ARE IGNOERD IF KEY=0
!     L IS NO. OF POINTS
!
!-----------------------------------------------------------------------
!Eric made DOUBLE PRECISION
  implicit none
  integer j,key,l
  real(kind=fPrec) a,b,count,e,scartx,scarty
  real(kind=fPrec) sumx,sumxx,sumxy,sumy,sumyy,x,xmed,y,ymed
  dimension x(l),y(l)                                                !hr07
  save
!-----------------------------------------------------------------------
!
!     CALCULATE SUMS
!
!-----------------------------------------------------------------------
  if(l-2.lt.0) goto 25
  if(l-2.ge.0) goto 1
1 count=0.0
  sumx=0.0
  sumy=0.0
  sumxy=0.0
  sumxx=0.0
  sumyy=0.0
  do 10 j=1,l
  if(y(j).eq.0..and.key.eq.0) goto 10
  sumx=sumx+x(j)
  sumy=sumy+y(j)
  count=count+1.0
10 continue
  if(count.le.1.) goto 25
  ymed=sumy/count
  xmed=sumx/count
  do 20 j=1,l
  if(y(j).eq.0..and.key.eq.0) goto 20
  scartx=x(j)-xmed
  scarty=y(j)-ymed
  sumxy=sumxy+scartx   *scarty
  sumxx=sumxx+scartx   *scartx
  sumyy=sumyy+scarty   *scarty
20 continue
!
!     FIT PARAMETERS
  if(sumxx.eq.0.) goto 25
  a=sumxy/sumxx
  b=ymed-a*xmed
  if(count.lt.3.) goto 101
  e=(sumyy-sumxy*a          )/(count-2.0)
  goto 100
!
!     ISUFFICIENT POINTS
25 a=0.0
  b=0.0
101 e=0.0
100 return
end subroutine lfitd

end module postprocessing

subroutine hbook2(i1,c1,i2,r1,r2,i3,r3,r4,r5)
  implicit none
  integer i1,i2,i3
  real r1,r2,r3,r4,r5
  character c1
  save
  return
end subroutine hbook2

subroutine hdelet(i1)
  implicit none
  integer i1
  save
  return
end subroutine hdelet

subroutine hlimit(i1)
  implicit none
  integer i1
  save
  return
end subroutine hlimit

subroutine hplax(c1,c2)
  implicit none
  character c1,c2
  save
  return
end subroutine hplax

subroutine hplcap(i1)
  implicit none
  integer i1
  save
  return
end subroutine hplcap

subroutine hplend()
  implicit none
  save
  return
end subroutine hplend

subroutine hplint(i1)
  implicit none
  integer i1
  save
  return
end subroutine hplint

subroutine hplopt(c1,i1)
  implicit none
  integer i1
  character c1
  save
  return
end subroutine hplopt

subroutine hplot(i1,c1,c2,i2)
  implicit none
  integer i1,i2
  character c1,c2
  save
  return
end subroutine hplot

subroutine hplset(c1,r1)
  implicit none
  real r1
  character c1
  save
  return
end subroutine hplset

subroutine hplsiz(r1,r2,c1)
  implicit none
  real r1,r2
  character c1
  save
  return
end subroutine hplsiz

subroutine hplsof(r1,r2,c1,r3,r4,r5,i1)
  implicit none
  integer i1
  real r1,r2,r3,r4,r5
  character c1
  save
  return
end subroutine hplsof

subroutine htitle(c1)
  implicit none
  character c1
  save
  return
end subroutine htitle

subroutine ipl(i1,r1,r2)
  implicit none
  integer i1
  real r1(*),r2(*)
  save
  return
end subroutine ipl

subroutine ipm(i1,r1,r2)
  implicit none
  integer i1
  real r1,r2
  save
  return
end subroutine ipm

subroutine iselnt(i1)
  implicit none
  integer i1
  save
  return
end subroutine iselnt

subroutine igmeta(i1,i2)
  implicit none
  integer i1,i2
  save
  return
end subroutine igmeta
