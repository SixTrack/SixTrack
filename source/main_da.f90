! ================================================================================================ !
!
!  SIXTRACK
! =========
!  SIXDIMENSIONAL PARTICLE-TRACKING
!
!  DIFFERENTIAL ALGEBRA INCLUDED
!  ONE TURN MAP
!  NO POSTPROCESSING FORSEEN
!
!  DEVELOPPED FROM <RACETRACK> A. WRULICH (DESY 84-026)
! ================================================================================================ !
!  USED DISKS:
!
!  GEOMETRY AND STRENGTH OF THE ACCELERATOR : UNIT  2
!  TRACKING PARAMETER                       : UNIT  3
!  NORMAL PRINTOUT                          : UNIT  6
!  TRACKING DATA                            : UNIT  8
! ================================================================================================ !
program mainda

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use, intrinsic :: iso_fortran_env, only : output_unit
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_units
  use mod_meta
  use mod_time
  use mod_alloc,    only : alloc_init
  use mod_fluc,     only : fluc_randomReport, fluc_errAlign, fluc_errZFZ
  use read_write,   only : readFort33
  use mod_geometry, only : geom_reshuffleLattice
  use mod_version

  implicit none

  integer i,iation,itiono,idate,im,itime,ix,izu,j,k,kpz,kzz,l,ll,ncorruo,ndim,nlino,nlinoo,nmz
  real(kind=fPrec) alf0s1,alf0s2,alf0s3,alf0x2,alf0x3,alf0z2,alf0z3,amp00,bet0s1,bet0s2,bet0s3,     &
    bet0x2,bet0x3,bet0z2,bet0z3,clo0,clop0,dp0,dp10,e0f,eps,epsa,gam0s1,gam0s2,gam0s3,gam0x1,gam0x2,&
    gam0x3,gam0z1,gam0z2,gam0z3,phag,qw,qwc,r0,r0a,rv,&
    tas,tas16,tas26,tas36,tas46,tas56,tas61,tas62,tas63,tas64,tas65

  character(len=8)  cdate,ctime ! Note: Keep in sync with maincr. If the len changes, CRCHECK will break.
  dimension qw(2),qwc(3),clo0(2),clop0(2)
  dimension eps(2),epsa(2)
  dimension tas(6,6)

  ! New Variables
  character(len=:), allocatable :: featList
  character(len=23) timeStamp
  character(len=8)  tsDate
  character(len=10) tsTime
  logical fErr

  ! Features
featList = ""
#ifdef TILT
  featList = featList//" TILT"
#endif
#ifdef FAST
  featList = featList//" FAST"
#endif
#ifdef STF
  featList = featList//" STF"
#endif
#ifdef CRLIBM
  featList = featList//" CRLIBM"
  call disable_xp()
#endif
#ifdef FIO
  featList = featList//" FIO"
#endif

  ! Set to nonzero before calling abend in case of error.
  ! If prror is called, it will be set internally.
  errout = 0

#ifndef CR
  lout=output_unit
#endif
  call f_initUnits
  call meta_initialise ! The meta data file.
  call time_initialise ! The time data file. Need to be as early as possible as it sets cpu time 0.
  call alloc_init      ! Initialise mod_alloc
  call allocate_arrays ! Initial allocation of memory

  ! Open files
  fErr = .false.
  call f_open(unit=12,file="fort.12",formatted=.true.,mode="w", err=fErr)
  call f_open(unit=17,file="fort.17",formatted=.true.,mode="rw",err=fErr) ! DA Files
  call f_open(unit=18,file="fort.18",formatted=.true.,mode="rw",err=fErr) ! DA Files
  call f_open(unit=19,file="fort.19",formatted=.true.,mode="rw",err=fErr) ! DA Files
  call f_open(unit=20,file="fort.20",formatted=.true.,mode="rw",err=fErr) ! DA Files
  call f_open(unit=21,file="fort.21",formatted=.true.,mode="rw",err=fErr) ! DA Files
  call f_open(unit=22,file="fort.22",formatted=.true.,mode="rw",err=fErr) ! DA Files
  call f_open(unit=23,file="fort.23",formatted=.true.,mode="rw",err=fErr) ! DA Files
  call f_open(unit=24,file="fort.24",formatted=.true.,mode="rw",err=fErr) ! DA Files
  call f_open(unit=25,file="fort.25",formatted=.true.,mode="rw",err=fErr) ! DA Files

  call f_open(unit=110,file="fort.110",formatted=.false.,mode="w", err=fErr)
  call f_open(unit=111,file="fort.111",formatted=.false.,mode="rw",err=fErr)

  call time_timeStamp(time_afterFileUnits)

  ! Print Header Info
  call time_timerStart
  call time_timerCheck(time0)

  ! TimeStamp
  call date_and_time(tsDate,tsTime)
  timeStamp = tsDate(1:4)//"-"//tsDate(5:6)//"-"//tsDate(7:8)//" "//&
              tsTime(1:2)//":"//tsTime(3:4)//":"//tsTime(5:10)

  write(lout,"(a)") ""
  write(lout,"(a)") "    SixTrack DA :: Version "//trim(version)//" :: Released "//trim(moddate)
  write(lout,"(a)") "  "//repeat("=",128)
  write(lout,"(a)") "    Git SHA Hash: "//trim(git_revision)
  write(lout,"(a)") "    Built With:   "//trim(adjustl(featList))
  write(lout,"(a)") "    Start Time:   "//timeStamp
  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine

  call meta_write("SixTrackDAVersion", trim(version))
  call meta_write("ReleaseDate",       trim(moddate))
  call meta_write("GitHash",           trim(git_revision))
  call meta_write("StartTime",         timeStamp)

  meta_nPartTurn = 2

  ! Init stuff
  do i=1,2
    eps(i)=zero
    epsa(i)=zero
    qw(i)=zero
    qwc(i)=zero
  end do
  qwc(3)=zero

  call daten
  call time_timeStamp(time_afterDaten)
  if (ithick.eq.1) call allocate_thickarrays
  if(nord <= 0 .or. nvar <= 0) then
    write(lerr,"(a)") "MAINDA> ERROR Order and number of variables have to be larger than 0 to calculate a differential algebra map"
    call prror(-1)
  end if
  if(ithick.eq.1) write(lout,10020)
  if(ithick.eq.0) write(lout,10030)
  call geom_reshuffleLattice
  call ord
  if(allocated(zfz)) call fluc_randomReport
  call clorb(ded)

  do l=1,2
    clo0(l)=clo(l)
    clop0(l)=clop(l)
  end do

  call clorb(zero)

  do l=1,2
    ll=2*l
    di0(l)=(clo0(l)-clo(l))/ded
    dip0(l)=(clop0(l)-clop(l))/ded
  end do

  amp00=amp(1)
  iation=abs(ition)
  call corrorb
  if(irmod2.eq.1) call rmod(dp1)
  if(iqmod.ne.0) call qmod0
  if(ichrom.eq.1.or.ichrom.eq.3) call chroma
  if(iskew.ne.0) call decoup
  dp0=dp1

  !--FOR THE MOMENTUM-SCAN THE MOMENTUM IS THE SAME FOR BOTH PARTICLES
  exz(1,6)=dp1
  exz(2,6)=dp1
  if(ilin.eq.1.or.ilin.eq.3) then
    call linopt(dp1)
  endif
  if(isub.eq.1) call subre(dp1)
  if(ise.eq.1) call search(dp1)
  if(napx.eq.0) goto 160

  ! beam-beam element
  nlino=nlin
  nlin=0
  if(nbeam.ge.1) then
    do i=1,nele
      if((kz(i).eq.20).or.(kz(i).eq.15)) then
        nlin=nlin+1
        if(nlin.gt.nele) then
          write(lerr,"(a)") "MAINDA> ERROR Too many elements for linear optics write-out"
          call prror(-1)
        end if
        bezl(nlin)=bez(i)
      end if
    end do
  end if

  ! MULTIPOLE WITH THEIR RANDOM VALUES ADDED
  izu=0
  do i=1,iu
    ix=ic(i)
    if(ix.le.nblo) cycle
    ix=ix-nblo
    kpz=kp(ix)
    kzz=kz(ix)
    if(kpz.eq.6.or.kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) cycle
    if(kzz.eq.15) cycle
    if(iorg.lt.0) mzu(i)=izu
    izu=mzu(i)+1
    smizf(i)=zfz(izu)*ek(ix)
    smi(i)=sm(ix)+smizf(i)
    izu=izu+1
    xsi(i)=xpl(ix)+zfz(izu)*xrms(ix)
    izu=izu+1
    zsi(i)=zpl(ix)+zfz(izu)*zrms(ix)
    if(mout2.eq.1) then
      if(kzz.eq.11) zfz(izu-2)=zero
      if(abs(ek(ix)).le.pieni) zfz(izu-2)=zero
      if(abs(xrms(ix)).le.pieni) zfz(izu-1)=zero
      if(abs(zrms(ix)).le.pieni) zfz(izu)=zero
      if(icextal(i) > 0) then
        write(31,"(a48,1p,d19.11,2d14.6,d17.9)") bez(ix),zfz(izu-2),zfz(izu-1),zfz(izu),fluc_errAlign(3,icextal(i))
      else if(icextal(i) < 0) then
        write(31,"(a48,1p,d19.11,2d14.6,d17.9)") bez(ix),zfz(izu-2),zfz(izu-1),zfz(izu),fluc_errZFZ(4,-icextal(i))
      else
        write(31,"(a48,1p,d19.11,2d14.6,d17.9)") bez(ix),zfz(izu-2),zfz(izu-1),zfz(izu),zero
      end if
    end if
    if(kzz.eq.11) then
      !Very similar to block "multini"
      r0=ek(ix)
      if(abs(r0).le.pieni) cycle
      nmz=nmu(ix)
      if(nmz.eq.0) then
        izu=izu+2*mmul
        cycle
      end if
      im=irm(ix)
      r0a=one

      do k=1,nmz
        izu=izu+1
        aaiv(k,i)=(ed(ix)*(ak0(im,k)+zfz(izu)*aka(im,k)))/r0a         !hr08
        izu=izu+1
        bbiv(k,i)=(ed(ix)*(bk0(im,k)+zfz(izu)*bka(im,k)))/r0a         !hr08
        r0a=r0a*r0
      end do

      izu=izu+2*mmul-2*nmz
    end if
  end do
  dp10=dp1
  dp1=zero
  if(ichrom > 1) then
    itiono=ition
    ition=0
    call chromda
    ition=itiono
  else
    itiono = 0 ! -Wmaybe-uninitialized
  end if
  dp1=dp10
  if(idp /= 1 .or. iation /= 1) iclo6 = 0
  if(iclo6 == 1 .or. iclo6 == 2) then
    if(iclo6r == 0) then
      clo6(1)  = clo(1)
      clop6(1) = clop(1)
      clo6(2)  = clo(2)
      clop6(2) = clop(2)
      clo6(3)  = zero
      clop6(3) = zero
    else
      write(lout,"(a)") "MAINDA> Reading closed orbit guess from fort.33"
      call readFort33
    end if
    call clorb(zero)
    call betalf(zero,qw)
    call phasad(zero,qwc)
    sigm(1) = clo6(3)
    dps(1)  = clop6(3)
    call qmodda(3,qwc)
    if(ilin >= 2) then
      nlinoo = nlin
      nlin   = nlino
      ilinc  = 1
      call mydaini(2,2,6,3,6,1)
      nlin   = nlinoo
    end if
    dp1 = dp10+clop6(3)
#ifdef DEBUG
!     call dumpbin('ecdclor6',1,3)
!     call abend('end cd clor6                                      ')
#endif
    do i=1,6
      do j=1,6
        tas(i,j)=tasm(i,j)
      end do
    end do
  else
    ncorruo=ncorru
    ncorru=1
    call clorb(zero)
    call betalf(zero,qw)
    call phasad(zero,qwc)
    call clorb(dp1)
    call betalf(dp1,qw)
    call phasad(dp1,qwc)
    ncorru=ncorruo
    dps(1)=dp1
    if(nvar2.le.5) then
      itiono=ition
      ition=0
    end if
#ifdef DEBUG
!       write(*,*) '3rd call qmodda multipole???'
#endif
    call qmodda(2,qwc)
#ifdef DEBUG
!     call dumpbin('aqmodda',2,3)
!     call abend('after  qmodda 2 3                                 ')
#endif
    if(nvar2.le.5) ition=itiono
    if(nvar2.le.4.and.ithick.eq.1) call envar(dp1)
    if(ilin.ge.2) then
      nlinoo=nlin
      nlin=nlino
      iqmodc=2
      call mydaini(1,2,5,2,5,1)
      ilinc=1
      call mydaini(2,2,5,2,5,1)
      nlin=nlinoo
    end if
    do i=1,4
      do j=1,4
        tas(i,j)=tasm(i,j)
      end do
    end do
  end if

  tas16=tas(1,6)*c1m3
  tas26=tas(2,6)*c1m3
  tas36=tas(3,6)*c1m3
  tas46=tas(4,6)*c1m3
  tas56=tas(5,6)*c1m3
  tas61=tas(6,1)*c1e3
  tas62=tas(6,2)*c1e3
  tas63=tas(6,3)*c1e3
  tas64=tas(6,4)*c1e3
  tas65=tas(6,5)*c1e3
  bet0(1)=tas(1,1)**2+tas(1,2)**2                                    !hr08
  bet0x2 =tas(1,3)**2+tas(1,4)**2                                    !hr08
  bet0x3 =tas(1,5)**2+tas16**2                                       !hr08
  gam0x1 =tas(2,1)**2+tas(2,2)**2                                    !hr08
  gam0x2 =tas(2,3)**2+tas(2,4)**2                                    !hr08
  gam0x3 =tas(2,5)**2+tas26**2                                       !hr08
  alf0(1)=-one*(tas(1,1)*tas(2,1)+tas(1,2)*tas(2,2))                 !hr08
  alf0x2 =-one*(tas(1,3)*tas(2,3)+tas(1,4)*tas(2,4))                 !hr08
  alf0x3 =-one*(tas(1,5)*tas(2,5)+tas16*tas26)                       !hr08
  bet0(2)=tas(3,3)**2+tas(3,4)**2                                    !hr08
  bet0z2 =tas(3,1)**2+tas(3,2)**2                                    !hr08
  bet0z3 =tas(3,5)**2+tas36**2                                       !hr08
  gam0z1 =tas(4,3)**2+tas(4,4)**2                                    !hr08
  gam0z2 =tas(4,1)**2+tas(4,2)**2                                    !hr08
  gam0z3 =tas(4,5)**2+tas46**2                                       !hr08
  alf0(2)=-one*(tas(3,3)*tas(4,3)+tas(3,4)*tas(4,4))                 !hr08
  alf0z2 =-one*(tas(3,1)*tas(4,1)+tas(3,2)*tas(4,2))                 !hr08
  alf0z3 =-one*(tas(3,5)*tas(4,5)+tas36*tas46)                       !hr08
  bet0s1 =tas(5,5)**2+tas56**2                                       !hr08
  bet0s2 =tas(5,1)**2+tas(5,2)**2                                    !hr08
  bet0s3 =tas(5,3)**2+tas(5,4)**2                                    !hr08
  gam0s1 =tas65**2+tas(6,6)**2                                       !hr08
  gam0s2 =tas61**2+tas62**2                                          !hr08
  gam0s3 =tas63**2+tas64**2                                          !hr08
  alf0s1 =-one*(tas(5,5)*tas65+tas56*tas(6,6))                       !hr08
  alf0s2 =-one*(tas(5,1)*tas61+tas(5,2)*tas62)                       !hr08
  alf0s3 =-one*(tas(5,3)*tas63+tas(5,4)*tas64)                       !hr08
  if(ierro.eq.0) goto 90
  write(lout,10200) dp1
  goto 160

90 continue
  write(lout,10040)
  phag=(phas*c180e0)/pi                                               !hr08
  if((idp.eq.0).or.(abs(phas).le.pieni.and.ition.eq.0)) then
    write(lout,10140) qwc(1),clo(1),clop(1),bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,qwc(2),&
      clo(2),clop(2),bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2
  end if
  if(idp.eq.1.and.iation.eq.1.and.abs(phas).gt.pieni) then
    if(iclo6.eq.0) then
      write(lout,10120) phag,qwc(1),clo(1),clop(1),bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,qwc(2),&
        clo(2),clop(2),bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2
    else
      write(lout,10130) phag,qwc(1),clo6(1),clop6(1),bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,&
        bet0x3,alf0x3,gam0x3,qwc(2),clo6(2),clop6(2),bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2,&
        bet0z3,alf0z3,gam0z3,qwc(3),clo6(3),clop6(3),bet0s1,alf0s1,gam0s1,bet0s2,alf0s2,gam0s2,bet0s3,alf0s3,gam0s3
    end if
  end if
  if(idp.eq.1.and.ition.eq.0.and.abs(phas).gt.pieni) then
    write(lout,10160) phag,qwc(1),clo(1),clop(1),bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,qwc(2),&
      clo(2),clop(2),bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2
  end if
  if(idp.eq.1.and.abs(phas).le.pieni.and.iation.eq.1) then
    if(iclo6.eq.0) then
      write(lout,10180) qwc(1),clo(1),clop(1),bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,qwc(2),&
        clo(2),clop(2),bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2
    else
      write(lout,10190) qwc(1),clo6(1),clop6(1),bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,&
        bet0x3,alf0x3,gam0x3,qwc(2),clo6(2),clop6(2),bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2,&
        bet0z3,alf0z3,gam0z3,qwc(3),clo6(3),clop6(3),bet0s1,alf0s1,gam0s1,bet0s2,alf0s2,gam0s2,bet0s3,alf0s3,gam0s3
    end if
  end if

  write(lout,10050) dp1
  call anfb(tas)
  if(iclo6.eq.2) then
    x(1,1) = x(1,1) + clo6(1)
    x(1,2) = x(1,2) + clo6(2)
    y(1,1) = y(1,1) + clop6(1)
    y(1,2) = y(1,2) + clop6(2)
    sigm(1) = sigm(1) + clo6(3)
    dps(1) = dps(1) + clop6(3)
    x(2,1) = x(2,1) + clo6(1)
    x(2,2) = x(2,2) + clo6(2)
    y(2,1) = y(2,1) + clop6(1)
    y(2,2) = y(2,2) + clop6(2)
    sigm(2) = sigm(2) + clo6(3)
    dps(2) = dps(2) + clop6(3)
  end if
  do l=1,2
    epsa(l)=amp(l)**2/bet0(l)                                        !hr08
    eps(l)=epsa(l)*c1e6
    x(1,l)=x(1,l)+(clo(l)*real(idz(l),fPrec))*real(1-idfor,fPrec)                !hr08
    y(1,l)=y(1,l)+(clop(l)*real(idz(l),fPrec))*real(1-idfor,fPrec)               !hr08
  end do
  e0f=sqrt(e0**2-pma**2)                                             !hr08
  if(iclo6.eq.0) then
    write(lout,10080) clo(1),clop(1),clo(2),clop(2),idz(1),idz(2),iver, idfor,iclo6,ition
  else
    write(lout,10090) clo6(1),clop6(1),clo6(2),clop6(2),clo6(3),clop6(3), idz(1),idz(2),iver,idfor,iclo6,ition
  endif
  if(idfor.eq.1.and.iclo6.ne.2) goto 110
  ejf(1)=e0f*(one+dps(1))
  ejf(2)=e0f*(one+dps(2))
  ej(1)=sqrt(ejf(1)**2+pma**2)                                       !hr08
  ej(2)=sqrt(ejf(2)**2+pma**2)                                       !hr08
  goto 120
110 continue
  ejf(1)=sqrt(ej(1)**2-pma**2)                                       !hr08
  ejf(2)=sqrt(ej(2)**2-pma**2)                                       !hr08
120 continue
  write(lout,10060) x(1,1),y(1,1),x(1,2),y(1,2),sigm(1),dps(1), x(2,1),y(2,1),x(2,2),y(2,2),sigm(2),dps(2),e0,ej(1),ej(2)
  write(lout,10010) amp,epsa
  write(lout,10170)
  if(e0.gt.pieni) then
    rv=(ej(1)*e0f)/(e0*ejf(1))
    if(ithick == 1) call envars(1,dps(1),rv)
  else
    write(lerr,"(a)") "MAINDA> ERROR Zero or negative energy does not make much sense."
    call prror(-1)
  end if
  if(numl.eq.0.or.numlr.ne.0) then
    write(lout,10070)
    goto 160
  end if
  if(nsix.eq.1.and.nvar2.eq.6) then
    nsix=2
    nvar=nvar-1
    nvar2=5
  end if
  ndim=nvar2/2
  call mydaini(3,nord,nvar,ndim,nvar2,nord1)
  if(inorm.eq.1) call daliesix
! if(icorr.eq.1) then
!   if(nctype.eq.0) call coruord
!   if(nctype.eq.1) call coruglo
! end if
  if(nsix.eq.2) then
    call umschr(19,18)
    nvar2=6
    nvar=nvar+1
    call mydaini(3,nord,nvar,ndim,nvar2,nord1)
  end if
160 continue

!-----------------------------------------------------------------------
! We're done in mainda, no error :)
!-----------------------------------------------------------------------
  call time_timeStamp(time_beforeExit)
  call time_finalise
  call meta_finalise
  call closeUnits
#ifdef CR
  call abend('                                                  ')
#else
  stop
#endif

10010 format(/t10,'UNCOUPLED AMPLITUDES AND EMITTANCES:',&
             /t10,'AMPLITUDE-X = ',f15.3,10x,'AMPLITUDE-Y = ',f15.3, '  MM',&
             /t10,'EMITTANCE-X = ',f15.3,10x,'EMITTANCE-Y =  ',f15.3, '  PI*MRAD*MM')
10020 format(/t10,'STRUCTURE INPUT FILE HAS -THICK- LINEAR ELEMENTS'//)
10030 format(/t10,'STRUCTURE INPUT FILE HAS ONLY -THIN- LINEAR ELEMENTS'//)
10040 format(/131('-'))
10050 format(/t10,'REL. MOMENTUM DEVIATION=',f19.16 &
             /t10,'================================')
10060 format(/5x,'---- INITIAL COORD. OF TWIN-TRAJECTORIES'/ 15(10x,f47.33/))
10070 format(/5x,'NON SENSICAL INPUT: NUML = 0 OR NUMLR NOT 0')
10080 format(/5x,'---- CLOSED ORBIT AND DECOUPLING (1=COU,0=DECOU)'&
             /5x,'/CLX  /',f47.33&
             /5x,'/CLXP /',f47.33&
             /5x,'/CLY  /',f47.33&
             /5x,'/CLYP /',f47.33&
             /5x,'/DCX  /',i13&
             /5x,'/DCY  /',i13&
             /5x,'/IVER /',i13&
             /5x,'/IDFOR/',i13&
             /5x,'/ICLO6/',i13&
             /5x,'/ITION/',i13&
             /5x/)
10090 format(/5x,'---- CLOSED ORBIT AND DECOUPLING (1=COU,0=DECOU)'&
             /5x,'/CLX  /',f47.33&
             /5x,'/CLXP /',f47.33&
             /5x,'/CLY  /',f47.33&
             /5x,'/CLYP /',f47.33&
             /5x,'/CLS  /',f47.33&
             /5x,'/CLSP /',f47.33&
             /5x,'/DCX  /',i13&
             /5x,'/DCY  /',i13&
             /5x,'/IVER /',i13&
             /5x,'/IDFOR/',i13&
             /5x,'/ICLO6/',i13&
             /5x,'/ITION/',i13&
             /5x/)
10120 format(/t10,'TRACKING WITH SYNCHROTRON OSCILLATIONS'// 15x,       &
     &'ACCELERATION WITH PHASE = ',f8.4/ t15,                           &
     &'       TUNE             CLO                CLOP           ',     &
     &'     BET0             ALF0           GAMMA      '//              &
     &t10,'  X  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  Y  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10))
10025 format(/t10,'Run started from binary dump file # 32')
10130 format(/t10,'TRACKING WITH SYNCHROTRON OSCILLATIONS'// 15x,       &
     &'ACCELERATION WITH PHASE = ',f8.4/ t15,                           &
     &'       TUNE             CLO                CLOP           ',     &
     &'     BET0             ALF0           GAMMA      '//              &
     &t10,'  X  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  Y  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  S  ',3(1x,ES17.10),3(1x,ES17.10)/                          &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10))
10140 format(/t10,'TRACKING FOR CONSTANT MOMENTUM DEVIATION'// 15x,     &
     &'------ NO ACCELERATION ------'// t15,                            &
     &'       TUNE         CLO            CLOP           ',             &
     &'   BET0           ALF0           GAMMA      '//                  &
     &t10,'  X  ',f14.10,2(1x,g15.8),1x,f15.9,1x,f15.10,f15.9/          &
     &t60,f15.9,1x,f15.10,f15.9/                                        &
     &t10,'  Y  ',f14.10,2(1x,g15.8),1x,f15.9,1x,f15.10,f15.9/          &
     &t60,f15.9,1x,f15.10,f15.9/)
10150 format(t5//t5,'BACK-TRACKING'/ t5, '============='//)
10160 format(t10,'TRACKING FOR CONSTANT MOMENTUM DEVIATION'// 15x,      &
     &'ACCELERATION WITH PHASE = ',f8.4/ t15,                           &
     &'       TUNE         CLO            CLOP           ',             &
     &'   BET0           ALF0           GAMMA      '//                  &
     &t10,'  X  ',f14.10,2(1x,g15.8),1x,f15.9,1x,f15.10,f15.9/          &
     &t60,f15.9,1x,f15.10,f15.9/                                        &
     &t10,'  Y  ',f14.10,2(1x,g15.8),1x,f15.9,1x,f15.10,f15.9/          &
     &t60,f15.9,1x,f15.10,f15.9/)
10170 format(//131('-')//t10,16('O')/t10,2('O'),12x,2('O')/t10,         &
     &'OO  TRACKING  OO', /t10,2('O'),12x,2('O')/t10,16('O')//131('-')//&
     &)
10180 format(/t10,'TRACKING WITH SYNCHROTRON OSCILLATIONS'// 15x,       &
     &'------ NO ACCELERATION ------'// t15,                            &
     &'       TUNE             CLO                CLOP           ',     &
     &'     BET0             ALF0           GAMMA      '//              &
     &t10,'  X  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  Y  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10))
10190 format(/t10,'TRACKING WITH SYNCHROTRON OSCILLATIONS'// 15x,       &
     &'------ NO ACCELERATION ------'// t15,                            &
     &'       TUNE             CLO                CLOP           ',     &
     &'     BET0             ALF0           GAMMA      '//              &
     &t10,'  X  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  Y  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  S  ',3(1x,ES17.10),3(1x,ES17.10)/                          &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10))
10200 format(t10,'NO OPTICAL SOLUTION FOR',2x,f19.16,2x,'RELATIVE MOMENTUM DEVIATION')

end program mainda
