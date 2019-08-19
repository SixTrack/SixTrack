! ================================================================================================ !
!  Linear Optics Calculations Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Updated: 2019-08-08
! ================================================================================================ !
module mod_linopt

  implicit none

  ! Linear optics files, formerly fort.34 and fort.11
  character(len=15), private, parameter :: linopt_dumpFile   = "linopt_dump.dat"
  character(len=18), private, parameter :: linopt_coupleFile = "linopt_coupled.dat"
  integer,           private, save      :: linopt_dumpUnit   = -1
  integer,           private, save      :: linopt_coupleUnit = -1

contains

! ================================================================================================ !
!  Parse Linear Optics Calculation Line
!  Rewritten from code from DATEN by VKBO
!  Updated: 2019-08-08
! ================================================================================================ !
subroutine linopt_parseInputLine(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common
  use sixtrack_input

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) mode
  integer nSplit,i
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "LINE> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 1) then

    nlin = 0
    ilin = 1

    if(nSplit > 0) mode = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),nt,  iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),ilin,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),ntco,iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),eui, iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),euii,iErr)

    select case(mode)
    case("ELEMENT")
      iprint = 0
    case("BLOCK")
      iprint = 1
    case default
      write(lerr,"(a)") "LINE> ERROR Valid modes are 'BLOCK' or 'ELEMENT'"
      iErr = .true.
    end select

    if(ilin /= 1 .and. ilin /= 2) then
      write(lerr,"(a)") "LINE> ERROR Only 1 (4D) and 2 (6D) are valid options for ilin."
      iErr = .true.
    end if

    if(st_debug) then
      call sixin_echoVal("mode",mode,"LINE",iLine)
      call sixin_echoVal("nt",  nt,  "LINE",iLine)
      call sixin_echoVal("ilin",ilin,"LINE",iLine)
      call sixin_echoVal("ntco",ntco,"LINE",iLine)
      call sixin_echoVal("eui", eui, "LINE",iLine)
      call sixin_echoVal("euii",euii,"LINE",iLine)
    end if
    if(iErr) return

  else

    do i=1,nSplit
      nlin = nlin + 1
      if(nlin > nele) then
        write(lerr,"(2(a,i0))") "LINE> ERROR Too many elements for linear optics write out. Max is ",nele," got ",nlin
        iErr = .true.
        return
      end if
      bezl(nlin) = trim(lnSplit(i))
    end do

  end if

end subroutine linopt_parseInputLine

!-----------------------------------------------------------------------
!  LINEAR PARAMETERS AT THE POSITION OF EVERY ELEMENT OR BLOCK
!-----------------------------------------------------------------------
subroutine linopt(dpp)

  use parpro
  use crcoall
  use mod_units
  use mod_common
  use mod_commons
  use mod_common_track
  use mod_settings
  use string_tools
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use collimation

#ifdef ROOT
  use root_output
#endif

#ifdef HDF5
  use hdf5_output
  use hdf5_linopt
#endif

  integer i,iiii,im,ium,ix,izu,j,jj,jk,jm,k,kpz,kzz,l,l1,ll,nmz,nr,dj
  real(kind=fPrec) aa,aeg,alfa,bb,benkr,beta,bexi,bezii,bl1eg,bl2eg,ci,cikve,clo0,clop0,cr,crkve, &
    crkveuk,di00,dip00,dphi,dpp,dpp1,dppi,dpr,dyy1,dyy2,ekk,etl,phi,phibf,puf,qu,qv,qw,qwc,r0,&
    r0a,t,xl,xs,zl,zs,quz,qvz
#ifdef TILT
  real(kind=fPrec) dyy11,qu1,tiltck,tiltsk
#endif
  character(len=mNameLen) idum

  dimension t(6,4)
  dimension beta(2),alfa(2),phibf(2),phi(2)
  dimension clo0(2),clop0(2),di00(2),dip00(2),qw(2),qwc(3)
  dimension aa(mmul),bb(mmul),dpr(6)
  dimension cr(mmul),ci(mmul)
  dimension aeg(nele,2,6),bl1eg(nblo,2,6),bl2eg(nblo,2,6)
  data dpr/6*zero/

  nhmoni = 0
  nvmoni = 0
  nhcorr = 0
  nvcorr = 0
  ium    = 6

  if(ncorru == 0) then
    write(lout,"(a)") str_divLine
    write(lout,"(a)") ""
    write(lout,"(a)") "    OOOOOOOOOOOOOOOOOOOOO"
    write(lout,"(a)") "    OO                 OO"
    write(lout,"(a)") "    OO  Linear Optics  OO"
    write(lout,"(a)") "    OO                 OO"
    write(lout,"(a)") "    OOOOOOOOOOOOOOOOOOOOO"
    write(lout,"(a)") ""
  end if

  dpr(:)   = zero
  t(:,:)   = zero

  beta(:)  = zero
  alfa(:)  = zero
  phibf(:) = zero
  phi(:)   = zero
  clo0(:)  = zero
  clop0(:) = zero
  di00(:)  = zero
  dip00(:) = zero
  qw(:)    = zero
  qwc(:)   = zero

  aa(:)    = zero
  bb(:)    = zero
  cr(:)    = zero
  ci(:)    = zero

  etl      = zero
  dpr(1)   = dpp*c1e3
  dpr(6)   = one
  dpp1     = dpp+ded

  call clorb(dpp1)
  clo0(1:2)  = clo(1:2)
  clop0(1:2) = clop(1:2)
  call clorb(dpp)
  write(lout,"(a)") ""

  do l=1,2
    ll        = 2*l
    di0(l)    = (clo0(l)-clo(l))/ded
    dip0(l)   = (clop0(l)-clop(l))/ded
    t(6,ll-1) = di0(l)
    t(6,ll)   = dip0(l)
  end do

  if(ncorru == 0) then
    call f_requestUnit(linopt_dumpFile, linopt_dumpUnit)
    call f_open(unit=linopt_dumpUnit,file=linopt_dumpFile,formatted=.true.,mode="w",status="replace")
    write(linopt_dumpUnit,"(a1,1x,a15,1x,a,1x,a4,5(1x,a16))") "#","len_tot",chr_rPad("element",mNameLen),&
      "type","strength","beta(1)","beta(2)","phi(1)","phi(2)"
    write(lout,"(a)") repeat("-",132)
    write(lout,"(a)")                  ""
    write(lout,"(a)")                  "     PLANE |             DISP(MM) |           DISP(MRAD)"
    write(lout,"(a)")                  "    -----------------------------------------------------"
    write(lout,"(a,f20.12,a3,f20.12)") "       X   | ",di0(1)," | ",dip0(1)
    write(lout,"(a,f20.12,a3,f20.12)") "       Y   | ",di0(2)," | ",dip0(2)
    write(lout,"(a)")                  ""
  end if

  call betalf(dpp,qw)
  call phasad(dpp,qwc)

  if(ierro /= 0) then
    write(lerr,"(a)") "LINOPT> ERROR No optical solution"
    call prror
  end if
  if(ncorru == 0) then
    write(lout,"(a)")        ""
    write(lout,"(a,f23.16)") "    Relative energy deviation : ",dpp
    write(lout,"(a,f23.16)") "    Horizontal tune           : ",qwc(1)
    write(lout,"(a,f23.16)") "    Vertical tune             : ",qwc(2)
    write(lout,"(a)")        ""
  end if

  call envar(dpp)

  if(ithick == 1) then
    call envardis(dpp1,aeg,bl1eg,bl2eg)
  end if

!--STARTVALUES OF THE TRAJECTORIES
  do l=1,2
    ll = 2*l
    t(1,ll-1) = clo(l)
    t(1,ll)   = clop(l)
  end do

  do i=1,4
    do j=1,4
      t(i+1,j) = ta(j,i)
      t(i+1,j) = ta(j,i)
    end do
  end do

  if(ncorru == 0 .and. st_quiet == 0) then
    write(lout,"(a)") repeat("-",132)
    write(lout,"(a)") ""
    write(lout,"(a)") "  LINEAR OPTICS CALCULATION WITH PRINTOUT AFTER EACH BLOCK"
    write(lout,"(a)") "  Note: Betatron phase calculation might be wrong by a multiple of 0.5 for each large block"
    write(lout,"(a)") ""
    write(lout,"(a)") repeat("-",132)
    write(lout,"(a)") "|  NR  | ELEM.  | L-TOTAL(M) |P|  PHI(2*PI) |   BETA(M)  |  ALFA(RAD)  |"//&
      "  GAMMA(M) |   DIS(M)  | DISP(RAD) |  CLO(MM)  | CLOP(MRAD)|"
    write(lout,"(a)") repeat("-",132)
  end if

!--START OF THE MACHINE
  idum = "START"
  nr   = 0
  call writelin(nr,idum,etl,phi,t,1,.false.,0)
  if(ntco /= 0) then
    if(mod(nr,ntco) == 0) then
      call cpltwis(idum,t,etl,phi)
    end if
  end if

#ifdef ROOT
  if(root_flag .and. root_Optics == 1) then
    call OpticsRootWrite()
  end if
#endif

!--STRUCTURE ELEMENT LOOP
  if(nt <= 0 .or. nt > iu) then
    nt = iu
  end if
  izu=0

#ifdef HDF5
  if(h5_writeOptics) then
    call h5lin_init
  end if
#endif

  STRUCTLOOP: do k=1,nt
    ix = ic(k)
    if(ix > nblo) goto 220 ! Not a BLOCK
    if(ithick == 1 .and. iprint == 1) goto 160

    jj=0 !initial idx
    dj=1 !step

    if(ix <= 0) then
      ix = -1*ix
      jj = mel(ix)+1 ! initial idx
      dj = -1        ! step
    end if
    jm = mel(ix)

!-- Loop over elements inside the block
    do j=1,jm
      jj = jj+dj       ! Subelement index of current sub=element
      jk = mtyp(ix,jj) ! Single-element index of the current sub-element
      if(ithick == 1 .and. kz(jk) /= 0) goto 120
      if(ithick == 0 .and. kz(jk) /= 0) then
        etl=etl+el(jk)
        write(lerr,"(a)") "LINOPT> ERROR In block '"//trim(bezb(ix))//"': found a thick non-drift element '"//&
          trim(bez(jk))//"' while ithick=1. This should not be possible!"
        call prror
        cycle STRUCTLOOP
      end if

!--IN BLOCK: PURE DRIFTLENGTH (above: If ITHICK=1 and kz!=0, goto 120->MAGNETELEMENT)
      etl = etl+el(jk)

      do l=1,2
        ll = 2*l
        if(abs(t(ll,ll-1)) > pieni) then
          phibf(l) = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
        else
          phibf(l) = pi2
        end if
        do i=1,ium
          t(i,ll-1) = t(i,ll-1)+t(i,ll)*(el(jk))
        end do
      end do

      do l=1,2
        ll=2*l
        if(abs(t(ll,ll-1)) > pieni) then
          dphi = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
        else
          dphi = pi2-phibf(l)
        end if
        if((-one*dphi) > pieni) then
          dphi = dphi+pi
        end if
        phi(l) = phi(l)+dphi/twopi
      end do

      nr = nr+1
      call writelin(nr,bez(jk),etl,phi,t,ix,.true.,k)
      if(ntco /= 0) then
        if(mod(nr,ntco) == 0) then
          call cpltwis(bez(jk),t,etl, phi)
        end if
      end if

#ifdef ROOT
      if(root_flag .and. root_Optics == 1) then
        call OpticsRootWrite()
      end if
#endif

      cycle

!--IN BLOCK: MAGNETELEMENT
120   continue
      if(kz(jk) /= 8) then
        etl = etl+el(jk)
      end if
      do l=1,2
        ll = 2*l

        if(abs(t(ll,ll-1)) > pieni) then
          phibf(l) = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
        else
          phibf(l) = zero
        end if

        puf = t(6,ll-1)
        t(6,ll-1) = (((((aeg(jk,l,1)*(t(1,ll-1)+puf*ded) + aeg(jk,l,2)*(t(1,ll) + t(6,ll)*ded)) &
          + aeg(jk,l,5)*dpp1*c1e3)- a(jk,l,1)*t(1,ll-1))-a(jk,l,2)*t(1,ll))- a(jk,l,5)*dpr(1))/ded
        t(6,ll)   = (((((aeg(jk,l,3)*(t(1,ll-1)+puf*ded) + aeg(jk,l,4)*(t(1,ll) + t(6,ll)*ded)) &
          + aeg(jk,l,6)*dpp1*c1e3)- a(jk,l,3)*t(1,ll-1))-a(jk,l,4)*t(1,ll))- a(jk,l,6)*dpr(1))/ded

        do i=1,ium-1
          puf=t(i,ll-1)
          t(i,ll-1) = (puf*a(jk,l,1)+t(i,ll)*a(jk,l,2))+dpr(i)*a(jk,l,5)
          t(i,ll)   = (puf*a(jk,l,3)+t(i,ll)*a(jk,l,4))+dpr(i)*a(jk,l,6)
        end do
      end do

      do l=1,2
        ll=2*l

        if(abs(t(ll,ll-1)) > pieni) then
          dphi = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
        else
          dphi = -one*phibf(l)
        end if

        if(kz(jk) /= 8 .and. -one*dphi > pieni) then
          dphi=dphi+pi
        end if
        phi(l) = phi(l)+dphi/twopi
      end do

      nr = nr+1
      call writelin(nr,bez(jk),etl,phi,t,ix,.true.,k)
      if(ntco /= 0) then
        if(mod(nr,ntco) == 0) then
          call cpltwis(bez(jk),t,etl, phi)
        end if
      end if

#ifdef ROOT
      if(root_flag .and. root_Optics == 1) then
        call OpticsRootWrite()
      end if
#endif

    end do ! End of loop over elements inside block

    nr = nr+1
    call writelin(nr,bezb(ix),etl,phi,t,ix,.true.,k)
    if(ntco /= 0) then
      if(mod(nr,ntco) == 0) then
        call cpltwis(bezb(ix),t,etl,phi)
      end if
    end if
#ifdef ROOT
    if(root_flag .and. root_Optics == 1) then
      call OpticsRootWrite()
    end if
#endif

    cycle STRUCTLOOP

!--BETACALCULATION FOR SERIES OF BLOCKS (ix >= nblo.and.ithick == 1.and.iprint == 1)
160 continue !if ithick=1 and iprint=1:
    if(ix <= 0) goto 190
!--REGULAR RUN THROUGH BLOCKS
    etl = etl+elbe(ix)

    do l=1,2
      ll=2*l

      if(abs(t(ll,ll-1)) > pieni) then
        phibf(l) = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
      else
        phibf(l) = zero
      end if

      puf = t(6,ll-1)
      t(6,ll-1) = (((((bl1eg(ix,l,1)*(t(1,ll-1)+puf*ded) + bl1eg(ix,l,2)*(t(1,ll)+t(6,ll)*ded)) &
        + bl1eg(ix,l,5)*dpp1*c1e3)- bl1(ix,l,1)*t(1,ll-1))-bl1(ix,l,2)*t(1,ll))- bl1(ix,l,5)*dpr(1))/ded
      t(6,ll)   = (((((bl1eg(ix,l,3)*(t(1,ll-1)+puf*ded) + bl1eg(ix,l,4)*(t(1,ll)+t(6,ll)*ded)) &
        + bl1eg(ix,l,6)*dpp1*c1e3)- bl1(ix,l,3)*t(1,ll-1))-bl1(ix,l,4)*t(1,ll))- bl1(ix,l,6)*dpr(1))/ded

      do i=1,ium-1
        puf = t(i,ll-1)
        t(i,ll-1) = (bl1(ix,l,1)*puf+bl1(ix,l,2)*t(i,ll))+dpr(i)*bl1(ix,l,5)
        t(i,ll)   = (bl1(ix,l,3)*puf+bl1(ix,l,4)*t(i,ll))+dpr(i)*bl1(ix,l,6)
      end do
    end do

    do l=1,2
      ll=2*l
      if(abs(t(ll,ll-1)) > pieni) then
        dphi = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
      else
        dphi = -one*phibf(l)
      end if
      if(-one*dphi > pieni) then
        dphi = dphi+pi
      end if
      phi(l) = phi(l)+dphi/twopi
    end do

    nr = nr+1
    call writelin(nr,bezb(ix),etl,phi,t,ix,.true.,k)
    if(ntco /= 0) then
      if(mod(nr,ntco) == 0) then
        call cpltwis(bezb(ix),t,etl,phi)
      end if
    end if
#ifdef ROOT
    if(root_flag .and. root_Optics == 1) then
      call OpticsRootWrite()
    end if
#endif

    cycle STRUCTLOOP

!--REVERSE RUN THROUGH BLOCKS (ix <= 0)
190 ix  = -ix
    etl = etl+elbe(ix)
    do l=1,2
      ll=2*l

      if(abs(t(ll,ll-1)) > pieni) then
        phibf(l) = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
      else
        phibf(l) = zero
      end if

      puf = t(6,ll-1)
      t(6,ll-1) = (((((bl2eg(ix,l,1)*(t(1,ll-1)+puf*ded) + bl2eg(ix,l,2)*(t(1,ll)+t(6,ll)*ded)) &
        + bl2eg(ix,l,5)*dpp1*c1e3)- bl2(ix,l,1)*t(1,ll-1))-bl2(ix,l,2)*t(1,ll))- bl2(ix,l,5)*dpr(1))/ded
      t(6,ll)   = (((((bl2eg(ix,l,3)*(t(1,ll-1)+puf*ded) + bl2eg(ix,l,4)*(t(1,ll)+t(6,ll)*ded)) &
        + bl2eg(ix,l,6)*dpp1*c1e3)- bl2(ix,l,3)*t(1,ll-1))-bl2(ix,l,4)*t(1,ll))- bl2(ix,l,6)*dpr(1))/ded

      do i=1,ium-1
        puf = t(i,ll-1)
        t(i,ll-1) = (bl2(ix,l,1)*puf+bl2(ix,l,2)*t(i,ll))+dpr(i)*bl2(ix,l,5)
        t(i,ll)   = (bl2(ix,l,3)*puf+bl2(ix,l,4)*t(i,ll))+dpr(i)*bl2(ix,l,6)
      end do
    end do

    do l=1,2
      ll = 2*l

      if(abs(t(ll,ll-1)) > pieni) then
        dphi = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
      else
        dphi = -phibf(l)
      end if

      if(-one*dphi > pieni) then
        dphi = dphi+pi
      end if
      phi(l) = phi(l)+dphi/twopi
    end do

    nr = nr+1
    call writelin(nr,bezb(ix),etl,phi,t,ix,.true.,k)
    if(ntco /= 0) then
      if(mod(nr,ntco) == 0) then
        call cpltwis(bezb(ix),t,etl,phi)
      end if
    end if
#ifdef ROOT
    if(root_flag .and. root_Optics == 1) then
      call OpticsRootWrite()
    end if
#endif

    cycle STRUCTLOOP

    ! NOT A BLOCK / Nonlinear insertion
220 ix   = ix-nblo
    qu   = zero
    qv   = zero
    dyy1 = zero
    dyy2 = zero
    kpz  = kp(ix)
    kzz  = kz(ix)

    ! Cavity
    if(kpz == 6 .or. abs(kzz) == 12) then
      nr = nr+1
      call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
      if(ntco /= 0) then
        if(mod(nr,ntco) == 0) then
          call cpltwis(bez(ix),t,etl,phi)
        end if
      end if
#ifdef ROOT
      if(root_flag .and. root_Optics == 1) then
        call OpticsRootWrite()
      end if
#endif

      cycle STRUCTLOOP
    end if

    ! Beam Beam element .and. fort.3 has BB block
    if(kzz == 20 .and. nbeam >= 1) then
      nbeam = k
      nr    = nr+1
      call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
      if(ntco /= 0) then
        if(mod(nr,ntco) == 0) call cpltwis(bez(ix),t,etl,phi)
      end if
#ifdef ROOT
      if(root_flag .and. root_Optics == 1) then
        call OpticsRootWrite()
      end if
#endif
      cycle STRUCTLOOP
    end if

    ! if kzz==22, starts a do over l; Update t matrix
    if(kzz == 22) then
      do l=1,2
        ll = 2*l
        if(abs(t(ll,ll-1)) > pieni) then
          phibf(l) = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
        else
          phibf(l) = zero
        end if
        do i=1,ium
          puf = t(i,ll-1)
          t(i,ll-1) = (puf*rrtr(imtr(ix),ll-1,ll-1)+t(i,ll)*rrtr(imtr(ix),ll-1,ll))+dpr(i)*rrtr(imtr(ix),ll-1,6)
          t(i,ll)   = (puf*rrtr(imtr(ix),ll,ll-1)+t(i,ll)*rrtr(imtr(ix),ll,ll))+dpr(i)*rrtr(imtr(ix),ll,6)
        end do
        t(1,ll-1) = t(1,ll-1)+cotr(imtr(ix),ll-1)
        t(1,ll)   = t(1,ll)+cotr(imtr(ix),ll)
        if(abs(t(ll,ll-1)) > pieni) then
          dphi = atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
        else
          dphi = -one*phibf(l)
        end if
        if(-one*dphi > pieni) then
          dphi = dphi+pi
        end if
        phi(l) = phi(l)+dphi/twopi
      end do
    end if

    ! Marker, beam-beam, phase-trombone, crab cavity (incl. multipole), or wire
    if(kzz == 0 .or. kzz == 20 .or. kzz == 22 .or. abs(kzz) == 23 .or. abs(kzz) == 26 .or. &
      abs(kzz) == 27 .or. abs(kzz) == 28 .or. abs(kzz) == 15) then

      nr = nr+1
      call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
      if(ntco /= 0) then
        if(mod(nr,ntco) == 0) then
          call cpltwis(bez(ix),t,etl,phi)
        end if
      end if
#ifdef ROOT
      if(root_flag .and. root_Optics == 1) then
        call OpticsRootWrite()
      end if
#endif
      cycle STRUCTLOOP
    end if

    ! Update the matrix etc. for supported blocks
    dyy1 = zero
    dyy2 = zero
    if(iorg < 0) then
      mzu(k) = izu
    end if
    izu = mzu(k)+1
    ekk = (sm(ix)+zfz(izu)*ek(ix))/(one+dpp)
    izu = izu+1
    xs  = xpl(ix)+zfz(izu)*xrms(ix)
    izu = izu+1
    zs  = zpl(ix)+zfz(izu)*zrms(ix)
#include "include/alignl.f90"

    if(kzz >= 0) then
      select case(kzz)

      case (1)
!--HORIZONTAL DIPOLE
        ekk = ekk*c1e3
#include "include/kickl01h.f90"
#include "include/kickq01h.f90"
!--NORMAL QUADRUPOLE
      case(2)
#include "include/kicklxxh.f90"
#include "include/kickq02h.f90"
!--   NORMAL SEXTUPOLE
      case(3)
        ekk = ekk*c1m3
#include "include/kickq03h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL OCTUPOLE
      case(4)
        ekk = ekk*c1m6
#include "include/kicksho.f90"
#include "include/kickq04h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL DECAPOLE
      case(5)
        ekk = ekk*c1m9
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq05h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL DODECAPOLE
      case(6)
        ekk = ekk*c1m12
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq06h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL 14-POLE
      case(7)
        ekk = ekk*c1m15
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq07h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL 16-POLE
      case(8)
        ekk = ekk*c1m18
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq08h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL 18-POLE
      case(9)
        ekk = ekk*c1m21
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq09h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL 20-POLE
      case(10)
        ekk = ekk*c1m24
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq10h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--Multipole block
      case(11)
        r0 = ek(ix)
        if(abs(dki(ix,1)) > pieni) then
          if(abs(dki(ix,3)) > pieni) then
#include "include/multl01.f90"
#include "include/multl08.f90"
            do i=2,ium
#include "include/multl02.f90"
            end do
          else
#include "include/multl03.f90"
#include "include/multl09.f90"
          end if
        end if
        if(abs(dki(ix,2)) > pieni) then
          if(abs(dki(ix,3)) > pieni) then
#include "include/multl04.f90"
#include "include/multl10.f90"
            do i=2,ium
#include "include/multl05.f90"
            end do
          else
#include "include/multl06.f90"
#include "include/multl11.f90"
          end if
        end if
        if(abs(r0) <= pieni) then
          cycle STRUCTLOOP
        end if
        nmz = nmu(ix)
        if(nmz == 0) then
          izu=izu+2*mmul
          nr=nr+1
          call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
          if(ntco /= 0) then
            if(mod(nr,ntco) == 0) then
              call cpltwis(bez(ix),t,etl,phi)
            end if
          end if
#ifdef ROOT
          if(root_flag .and. root_Optics == 1) then
            call OpticsRootWrite()
          end if
#endif

          cycle STRUCTLOOP
        end if
        im    = irm(ix)
        r0a   = one
        benkr = ed(ix)/(one+dpp)
        do l=1,nmz
#include "include/multl07a.f90"
        end do
        if(nmz >= 2) then
#include "include/multl07b.f90"
          do l=3,nmz
#include "include/multl07c.f90"
          end do
        else
#include "include/multl07d.f90"
        end if
#ifdef TILT
#include "include/multl07e.f90"
#endif
        izu = izu+2*mmul-2*nmz

!--Skipped elements
      case(12,13,14,15,16,17,18,19,20,21,22,23)
        cycle STRUCTLOOP

!--DIPEDGE ELEMENT
      case(24)
#include "include/kickldpe.f90"
#include "include/kickqdpe.f90"
!--solenoid
      case(25)
#include "include/kicklso1.f90"
#include "include/kickqso1.f90"

!--Skipped elements
      case(26,27,28)
        cycle STRUCTLOOP

!--Unrecognized element (incl. cav with kp /= 6 for non-collimat/bnlelens)
      case default
        nr = nr+1
        call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
        if(ntco /= 0) then
          if(mod(nr,ntco) == 0) then
            call cpltwis(bez(ix),t,etl,phi)
          end if
        end if
#ifdef ROOT
        if(root_flag .and. root_Optics == 1) then
          call OpticsRootWrite()
        end if
#endif
        cycle STRUCTLOOP
      end select

!--SKEW ELEMENTS
    else if(kzz < 0) then
      kzz = -kzz ! Make it positive
      select case(kzz)
      case(1)
!--VERTICAL DIPOLE
        ekk = ekk*c1e3
#include "include/kickl01v.f90"
#include "include/kickq01v.f90"
!--SKEW QUADRUPOLE
      case(2)
#include "include/kicklxxv.f90"
#include "include/kickq02v.f90"
!--SKEW SEXTUPOLE
      case(3)
        ekk = ekk*c1m3
#include "include/kickq03v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW OCTUPOLE
      case(4)
        ekk = ekk*c1m6
#include "include/kicksho.f90"
#include "include/kickq04v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW DECAPOLE
      case(5)
        ekk = ekk*c1m9
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq05v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW DODECAPOLE
      case(6)
        ekk = ekk*c1m12
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq06v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW 14-POLE
      case(7)
        ekk = ekk*c1m15
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq07v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW 16-POLE
      case(8)
        ekk = ekk*c1m18
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq08v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW 18-POLE
      case(9)
        ekk=ekk*c1m21
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq09v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW 20-POLE
      case(10)
        ekk = ekk*c1m24
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq10v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"

      case default
        ! Unrecognized skew element (including kzz=-12,kp /= 6 for non-collimat/bnlelens)
        nr = nr+1
        call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
        if(ntco /= 0) then
          if(mod(nr,ntco) == 0) then
            call cpltwis(bez(ix),t,etl,phi)
          end if
        end if
#ifdef ROOT
        if(root_flag .and. root_Optics == 1) then
          call OpticsRootWrite()
        end if
#endif
        cycle STRUCTLOOP
      end select
    end if

    ! Done processing an element: go here!
    t(6,2) = t(6,2)-dyy1/(one+dpp)
    t(6,4) = t(6,4)-dyy2/(one+dpp)
    t(1,2) = t(1,2)+dyy1
    t(1,4) = t(1,4)+dyy2
    do i=2,ium
      if(kzz == 24) then
        t(i,2) = (t(i,2)+t(i,1)*qu)-qv*t(i,3)
        t(i,4) = (t(i,4)-t(i,3)*quz)-qvz*t(i,1)
      elseif(kzz == 25) then !--solenoid
#include "include/phassolenoid.f90"
      else
        t(i,4) = (t(i,4)-t(i,3)*qu)-qv*t(i,1)
        t(i,2) = (t(i,2)+t(i,1)*qu)-qv*t(i,3)
      end if
    end do
    bexi  = t(2,1)**2+t(3,1)**2
    bezii = t(4,3)**2+t(5,3)**2
    if(ncorru == 0) then
      if(kz(ix) == 11) then
        if(abs(aa(2)) > pieni .and. nmz > 1) then
          write(linopt_dumpUnit,"(f17.9,1x,a,1x,i4,5(1x,1pe16.9))") etl,bez(ix),-2,aa(2),bexi,bezii,phi
        end if
        do iiii=3,nmz
          if(abs(bb(iiii)) > pieni) then
            write(linopt_dumpUnit,"(f17.9,1x,a,1x,i4,5(1x,1pe16.9))") etl,bez(ix),iiii,bb(iiii),bexi,bezii,phi
          end if
          if(abs(aa(iiii)) > pieni) then
            write(linopt_dumpUnit,"(f17.9,1x,a,1x,i4,5(1x,1pe16.9))") etl,bez(ix),-iiii,aa(iiii),bexi,bezii,phi
          end if
        end do
      else if(abs(ekk) > pieni .and. abs(kz(ix)) >= 3) then
        write(linopt_dumpUnit,"(f17.9,1x,a,1x,i4,5(1x,1pe16.9))") etl,bez(ix),kz(ix),ekk,bexi,bezii,phi
      else if(abs(ekk) > pieni .and. kz(ix) == -2) then
        write(linopt_dumpUnit,"(f17.9,1x,a,1x,i4,5(1x,1pe16.9))") etl,bez(ix),kz(ix),ekk,bexi,bezii,phi
      end if
    end if

    nr = nr+1
    call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
    if(ntco /= 0) then
      if(mod(nr,ntco) == 0) then
        call cpltwis(bez(ix),t,etl,phi)
      end if
    end if
#ifdef ROOT
    if(root_flag .and. root_Optics == 1) then
      call OpticsRootWrite()
    end if
#endif

  end do STRUCTLOOP ! END LOOP OVER ELEMENTS

#ifdef HDF5
  if(h5_writeOptics) call h5lin_saveData
#endif

  call clorb(ded)
  clo0(1:2)  = clo(1:2)
  clop0(1:2) = clop(1:2)
  call clorb(zero)
  do l=1,2
    ll = 2*l
    di0(l)  = (clo0(l)-clo(l))/ded
    dip0(l) = (clop0(l)-clop(l))/ded
  end do
  iiii  = 100
  idum  = "END"
  bexi  = t(2,1)**2+t(3,1)**2
  bezii = t(4,3)**2+t(5,3)**2
  if(ncorru == 0) then
    write(linopt_dumpUnit,"(f17.9,1x,a,1x,i4,5(1x,1pe16.9))") etl,idum,iiii,zero,bexi,bezii,phi
  end if
  if(ncorru == 0) then
    write(lout,"(a)") repeat("-",131)
  end if

  ! Close files
  if(linopt_coupleUnit /= -1) then
    call f_close(linopt_coupleUnit)
  end if
  if(ncorru == 0) then
    call f_close(linopt_dumpUnit)
  end if

end subroutine linopt

! ============================================================================ !
!  Write out linear optics parameters and send to modules that needs it
!  Updated: 2019-07-22
! ============================================================================ !
subroutine writelin(nr,typ,tl,p1,t,ixwl,isBLOC,ielem)

  use parpro
  use crcoall
  use scatter
  use mod_settings
  use mod_common
  use mod_commons
  use mod_common_track
  use collimation
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

#ifdef ROOT
  use iso_c_binding, only: C_NULL_CHAR
  use root_output
#endif

#ifdef HDF5
  use hdf5_output
  use hdf5_linopt
#endif

  integer i,iwrite,ixwl,l,ll,nr
  real(kind=fPrec) al1(2),al2(2),b1(2),b2(2),c(2),cp(2),d(2),dp(2),g1(2),g2(2),p1(2),t(6,4),tl
  character(len=mNameLen) typ
  ! isBLOC == TRUE if ixwl currently refers to a BLOC index, FALSE if it is a SINGLE ELEMENT index
  logical isBLOC
  integer ielem

#ifdef HDF5
  real(kind=fPrec) hdf5Data(17)
#endif

  iwrite = 0
  if(nlin == 0) then
    iwrite = 1
  else
    do i=1,nlin
      if(typ == bezl(i)) iwrite = 1
    end do
  end if
  if(iwrite == 1) then
    do l=1,2
      ll     = 2*l
      b1(l)  = t(ll,ll-1)**2+t(ll+1,ll-1)**2
      b2(l)  = t(6-ll,ll-1)**2+t(7-ll,ll-1)**2
      al1(l) = -one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))
      al2(l) = -one*(t(6-ll,ll-1)*t(6-ll,ll)+t(7-ll,ll-1)*t(7-ll,ll))
      g1(l)  = t(ll,ll)**2+t(ll+1,ll)**2
      g2(l)  = t(6-ll,ll)**2+t(7-ll,ll)**2
      d(l)   = t(6,ll-1)*c1m3
      dp(l)  = t(6,ll)*c1m3
      c(l)   = t(1,ll-1)
      cp(l)  = t(1,ll)
    end do

#ifdef ROOT
    if(root_flag .and. root_Optics == 1) then
      call OpticsRootWriteLin(nr, typ // C_NULL_CHAR,len(typ),tl,c(1),cp(1),c(2),cp(2),&
        b1(1),b1(2),al1(1),al1(2),d(1),d(2),dp(1),dp(2))
    end if
#endif
#ifdef HDF5
    if(h5_writeOptics) then
      hdf5Data(:) = (/tl,&
        p1(1),b1(1),al1(1),g1(1),d(1),dp(1),c(1),cp(1),&
        p1(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),cp(2)/)
      call h5lin_writeLine(nr, typ, hdf5Data)
    end if
#endif

    if(do_coll) then
      tbetax(max(ielem,1))  = b1(1)
      tbetay(max(ielem,1))  = b1(2)
      talphax(max(ielem,1)) = al1(1)
      talphay(max(ielem,1)) = al1(2)
      torbx(max(ielem,1))   = c(1)
      torbxp(max(ielem,1))  = cp(1)
      torby(max(ielem,1))   = c(2)
      torbyp(max(ielem,1))  = cp(2)
      tdispx(max(ielem,1))  = d(1)
      tdispy(max(ielem,1))  = d(2)
    end if

    if(scatter_active .and. .not. isBLOC) then
      call scatter_setLinOpt(ixwl, al1, b1, c, cp, d, dp)
    end if

    if(ncorru == 0) then
      if(st_quiet == 0) then
        write(lout,10000) nr,typ(:8),tl,p1(1),b1(1),al1(1),g1(1),d(1),dp(1),c(1),cp(1)
        write(lout,10010) b2(1),al2(1),g2(1)
        write(lout,10030) typ(9:16)
        write(lout,10020) p1(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),cp(2)
        write(lout,10010) b2(2),al2(2),g2(2)
        write(lout,"(a)") repeat("-",132)
      end if
    else
      if(.not.isBLOC) then
        if(kp(ixwl) == 3) then
          nhmoni = nhmoni + 1
          betam(nhmoni,1)  = b1(1)
          pam(nhmoni,1)    = p1(1)*twopi
          bclorb(nhmoni,1) = c(1)
        else if(kp(ixwl) == 4) then
          nhcorr = nhcorr + 1
          betac(nhcorr,1) = b1(1)
          pac(nhcorr,1)   = p1(1)*twopi
        else if(kp(ixwl) == -3) then
          nvmoni = nvmoni + 1
          betam(nvmoni,2)  = b1(2)
          pam(nvmoni,2)    = p1(2)*twopi
          bclorb(nvmoni,2) = c(2)
        else if(kp(ixwl) == -4) then
          nvcorr = nvcorr + 1
          betac(nvcorr,2) = b1(2)
          pac(nvcorr,2)   = p1(2)*twopi
        end if
      end if
    end if
  end if

  return
10000 format('|',i6,'|',a8,'|',f12.5,'|','X','|',f12.7,'|',f12.6,'|',f13.7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10010 format('|',6x,'|',8x,'|',12x,'|',1x,'|',12x,'|',f12.6,'|', f13.7,'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10020 format('|',6x,'|',8x,'|',12x,'|','Y','|',f12.7,'|',f12.6,'|', f13.7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10030 format('|',6x,'|',a8,'|',12x,'|',102('-'))
end subroutine writelin

!-----------------------------------------------------------------------
!  CALCULATES COUPLED TWISS PARAMETERS AROUND THE RING AND ALSO THE
!  ANGLE OF THE MAJOR AXIS OF A ELLIPSE IN THE X-Y PROJECTION WITH
!  THE X-AXIS. THE 4-D ELLIPSOID IS GIVEN BY THE BOUNDARY OF A
!  DISTRIBUTION OF PARTICLES WITH MAXIMUM EMITANCE OF MODE I AND II,
!  EUI AND EUII RESPECTIVELY.
!  BINARY PRINT ON FILE 11 OF 22 VALUES :
!  POSITION [M],
!  BET(1-4), ALF(1-4), GAM(1-4), COOR-PHI(1-4), COOR-PRIME-PHI(1-4),
!  COUUANGL
!-----------------------------------------------------------------------
subroutine cpltwis(typ,t,etl,phi)

  use parpro
  use mod_units
  use mod_common
  use mod_commons
  use mod_common_track
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
#ifdef ROOT
  use root_output
#endif

  integer i,iwrite
  real(kind=fPrec) alxi,alxii,alzi,alzii,bexi,bexii,bezi,bezii,couuang,etl,gaxi,gaxii,gazi,gazii,   &
    phi(2),phxi,phxii,phxpi,phxpii,phzi,phzii,phzpi,phzpii,t(6,4)
  character(len=mNameLen) typ

  iwrite = 0
  if(nlin == 0) then
    iwrite = 1
  else
    do i=1,nlin
      if(typ == bezl(i)) iwrite = 1
    end do
  end if
  if(iwrite == 1) then
    bexi  = t(2,1)**2+t(3,1)**2
    bexii = t(4,1)**2+t(5,1)**2
    bezi  = t(2,3)**2+t(3,3)**2
    bezii = t(4,3)**2+t(5,3)**2
    alxi  = -one*(t(2,1)*t(2,2)+t(3,1)*t(3,2))
    alxii = -one*(t(4,1)*t(4,2)+t(5,1)*t(5,2))
    alzi  = -one*(t(2,3)*t(2,4)+t(3,3)*t(3,4))
    alzii = -one*(t(4,3)*t(4,4)+t(5,3)*t(5,4))
    gaxi  = t(2,2)**2+t(3,2)**2
    gaxii = t(4,2)**2+t(5,2)**2
    gazi  = t(2,4)**2+t(3,4)**2
    gazii = t(4,4)**2+t(5,4)**2
    if(abs(t(2,1)) > pieni)  phxi   = atan2_mb(t(3,1),t(2,1))
    if(abs(t(4,1)) > pieni)  phxii  = atan2_mb(t(5,1),t(4,1))
    if(abs(t(4,1)) > pieni)  phxii  = atan2_mb(t(5,1),t(4,1))
    if(abs(t(2,3)) > pieni)  phzi   = atan2_mb(t(3,3),t(2,3))
    if(abs(t(4,3)) > pieni)  phzii  = atan2_mb(t(5,3),t(4,3))
    if(abs(t(2,2)) > pieni)  phxpi  = atan2_mb(t(3,2),t(2,2))
    if(abs(t(4,2)) > pieni)  phxpii = atan2_mb(t(5,2),t(4,2))
    if(abs(t(2,4)) > pieni)  phzpi  = atan2_mb(t(3,4),t(2,4))
    if(abs(t(4,4)) > pieni)  phzpii = atan2_mb(t(5,4),t(4,4))
    if(abs(t(2,1)) <= pieni) phxi   = pi2
    if(abs(t(4,1)) <= pieni) then
      if(bexii > pieni)  phxii = pi2
      if(bexii <= pieni) phxii = zero
    end if
    if(abs(t(2,3)) <= pieni) then
      if(bezi > pieni)  phzi = pi2
      if(bezi <= pieni) phzi = zero
    end if
    if(abs(t(4,3)) <= pieni) phzii = pi2
    if(abs(t(2,2)) <= pieni) phxpi = pi2
    if(abs(t(4,2)) <= pieni) then
      if(gaxii > pieni)  phxpii = pi2
      if(gaxii <= pieni) phxpii = zero
    end if
    if(abs(t(2,4)) <= pieni) then
      if(gazi > pieni)  phzpi = pi2
      if(gazi <= pieni) phzpi = zero
    end if
    if(abs(t(4,4)) <= pieni) phzpii = pi2
    if(abs(eui*(bexi-bezi)+euii*(bexii-bezii)) > pieni) then
      couuang = half*atan_mb((two*((eui*sqrt(bexi*bezi))*cos_mb(phxi-phzi)+(euii*sqrt(bexii*bezii)) &
        *cos_mb(phxii-phzii)))/ (eui*(bexi-bezi)+euii*(bexii-bezii)))
    else
      couuang = zero
    end if
    if(linopt_coupleUnit == -1) then
      call f_requestUnit(linopt_coupleFile, linopt_coupleUnit)
      call f_open(unit=linopt_coupleUnit,file=linopt_coupleFile,formatted=.true.,mode="w")
    end if
    write(linopt_coupleUnit,*) typ,etl,phi,bexi,bexii,bezi,bezii,       &
      alxi,alxii,alzi,alzii,gaxi,gaxii,gazi,gazii,                      &
      phxi,phxii,phzi,phzii,phxpi,phxpii,phzpi,phzpii,couuang,          &
      t(6,1),t(6,2),t(6,3),t(6,4),t(1,1),t(1,2),t(1,3),t(1,4)

#ifdef ROOT
    if(root_flag .and. root_Optics == 1) then
      call OpticsRootWriteCpl(phi(1),phi(2),bexi,bexii,bezi,bezii,      &
        alxi,alxii,alzi,alzii,gaxi,gaxii,gazi,gazii,                    &
        phxi,phxii,phzi,phzii,phxpi,phxpii,phzpi,phzpii,couuang,        &
        t(6,1),t(6,2),t(6,3),t(6,4),t(1,1),t(1,2),t(1,3),t(1,4))
    end if
#endif

  end if

end subroutine cpltwis

end module mod_linopt
