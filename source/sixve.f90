! This file contains the routines that are only used by the main version of SixTrack

!-----------------------------------------------------------------------
!  SUBROUTINE TO SUMMARIZE THE RESULTS OF THE POSTPROCESSING
!-----------------------------------------------------------------------
subroutine sumpos

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use crcoall
  use parpro
  use string_tools

  implicit none

  character(len=4) ch
  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn)       :: inLine
  real(kind=fPrec) d(60), dlost
  integer nSplit, ioStat, lineNo, i, j
  logical spErr, fErr
  
  save

  rewind 10
  lineNo = 0
  do i=1,1000
    read(10,"(a)",end=20,iostat=ioStat) inLine
    if(ioStat /= 0) then
      write(lout,"(a,i0)") "SUMPOS> ERROR Failed to read line from 'fort.10'. iostat = ",ioStat
      call prror(-1)
    end if
  
    call chr_split(inLine, lnSplit, nSplit, spErr)
    if(spErr) then
      write(lout,"(a)") "SUMPOS> ERROR Failed to parse line from 'fort.10'"
      call prror(-1)
    end if
    if(nSplit > 60) then
      write(lout,"(a,i0)") "SUMPOS> ERROR Too many elements on a single line of 'fort.10'. Max is 60, got ",nSplit
      call prror(-1)
    end if
    lineNo = lineNo+1

    do j=1,nSplit
      call chr_cast(lnSplit(j),d(j),spErr)
    end do
  
    if(i == 1) write(lout,10000)
    if(abs(d(2)) >= pieni) then
      ch = "LOST"
    else
      ch = "    "
    end if
    if(d(22) >= d(23)) then
      dlost = d(23)
    else
      dlost = d(22)
    end if
    write(lout,10010) nint(dlost),d(3),d(5),d(7),d(9),d(10),d(11),d(12),nint(d(16)),nint(d(18)),    &
      d(19),d(21),ch,d(4),d(6),d(8),d(13),nint(d(17)),d(20),d(25),d(14),d(15)
  end do
  20 continue
  rewind 10

  lineNo = 0
  write(lout,10020)
  do i=1,1000
    read(10,"(a)",end=40,iostat=ioStat) inLine
    if(ioStat /= 0) then
      write(lout,"(a,i0)") "SUMPOS> ERROR Failed to read line from 'fort.10'. iostat = ",ioStat
      call prror(-1)
    end if
  
    call chr_split(inLine, lnSplit, nSplit, spErr)
    if(spErr) then
      write(lout,"(a)") "SUMPOS> ERROR Failed to parse line from 'fort.10'"
      call prror(-1)
    end if
    if(nSplit > 60) then
      write(lout,"(a,i0)") "SUMPOS> ERROR Too many elements on a single line of 'fort.10'. Max is 60, got ",nSplit
      call prror(-1)
    end if
    lineNo = lineNo+1

    do j=1,nSplit
      call chr_cast(lnSplit(j),d(j),spErr)
    end do

    ! Now we are using 60 for CPU in seconds
    ! But note that dnms is now found in word 59.
    ! and we always print the maximum DMMAC as NMAC
    ! or zero which should really be OK I think.
    ! N.B. If particle is lost nms is 0, so we set mmac to zero too
    d(60) = one ! was real(nmac)
    if(nint(d(59)) == 0) d(60) = zero
    write(lout,10030) i,nint(d(59)),nint(d(60)),nint(d(59))*nint(d(24))
  end do

  40 continue
  write(lout,10040)
  return

10000 format(/131('-')/t10,'SUMMARY OF THE POSTPROCESSING' //t1,128(    &
  &'-'), /t1,'|',8x,'|',11x,'|',11x,'|',12x,'|',11x,                 &
  &'|NORMALIZED | SLOPE  |',14x,'|',10x,'|',21x,'|', /t1,            &
  &'|  TURN  |   LINEAR  |   BETA-   | AMPLITUDES | MOMENTUM  |',    &
  &'PHASESPACE | OF THE |  NONLINEAR   |  NEAREST |',7x,'SMEAR OF',6x&
  &,'|', /t1,                                                        &
  &'| NUMBER |   TUNES   | FUNCTIONS |            | DEVIATION |',    &
  &' DISTANCE  |DISTANCE|  DETUNING    | RESONANCE|   THE EMITTANCES'&
  &  ,4x,'|',/t1,128('-'), /t1,                                      &
  &'|        |           |     [M]   |     [MM]   |           |',    &
  &'           |        |              |     |ORD.|',                &
  &'    [%]  |      [%]  |'/t1,128('-'))
10010 format(t1,'|',i8,'|X ',f9.5,'|X ',f9.4,'|X ',f10.6,'|',d11.4, '|',&
  &d11.4,'|',f8.4,'|X ',d12.5,'|X ',i3,'| ',i2,' |X ', f7.3,'|X+Y ', &
  &f7.3,'|' /t1,'|  ',a4,'  |Y ',f9.5,'|Y ',f9.4,'|Y ',f10.6,'|',11x,&
  &'|',11x,'|',8x,'|+/- ',d10.3,'|Y ',i3,'|    |Y ', f7.3,'|    ',7x,&
  &'|' /t1,'|',8x,'|QS ',f8.6,'|  ',9x,'|  ',10x,'|',11x, '|',11x,'|'&
  &,8x,'|Y ',d12.5,'|  ',3x,'|    |  ', 7x,'|    ',7x,'|' /t1,'|',8x,&
  &'|  ',9x,'|  ',9x,'|  ',10x,'|',11x, '|',11x,'|',8x,'|+/- ',d10.3,&
  &'|  ',3x,'|    |  ', 7x,'|    ',7x,'|'/t1,128('-'))
10020 format(/131('-')/t10,'RANDOM SETS USED' //                        &
  &'  CASE  |  # OF RANDOM SET  |  MAX. POSSIBLE SETS   |    ',      &
  &' SEED'/65('-'))
10030 format(3x,i2,13x,i2,19x,i2,13x,i8)
10040 format(65('-')//131('-'))

end subroutine sumpos

!-----------------------------------------------------------------------
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 13-06-2014
!     calculate dcum, as done in linopt and when parsing BLOCs (daten):
!         lengths of thick lens elements are taken on the curvilinear
!         reference system; thus, no difference between the length
!         of SBENDs and the one of RBENDs, as they are both the ARC one;
!     for future needs:
!                ds=two/ed(ix)*asin(el(ix)*ed(ix)/two)
!     always in main code
!-----------------------------------------------------------------------
subroutine cadcum

  use floatPrecision
  use numerical_constants
  use crcoall
  use parpro
  use mod_common
  implicit none

  save

  ! temporary variables
  real(kind=fPrec) tmpdcum
  integer ientry, jentry, kentry, ix

  write(lout,"(a)") "CADCUM> Calculating machine length"

  ! initialise cumulative length
  tmpdcum=zero

  ! loop all over the entries in the accelerator structure
  do ientry=1,iu
    ix=ic(ientry)
    if(ix.gt.nblo) then
      ! SINGLE ELEMENT
      ix=ix-nblo
      if ( el(ix).gt.zero ) tmpdcum=tmpdcum+el(ix)
    else
      ! BLOC: iterate over elements
      do jentry=1,mel(ix)
        kentry=mtyp(ix,jentry)
        if( el(kentry).gt.zero ) tmpdcum=tmpdcum+el(kentry)
      enddo
    endif
!       assign value of dcum
    dcum(ientry)=tmpdcum
!     go to next entry in the acclerator structure
  enddo

!     assign the last value to the closing MARKER:
  dcum(iu+1)=tmpdcum

  if(print_dcum) then
    ! A useful printout. Enabled by the PRINT_DCUM flag in the SETTINGS block
    write(lout,10030) "CADCUM> ","ientry","ix","name"//repeat(" ",44),"    dcum [m]"
    write(lout,10020) "CADCUM> ",0,-1,"START"//repeat(" ",43),dcum(0)
    do ientry=1,iu
      ix=ic(ientry)
      if(ix.gt.nblo) then
        ! SINGLE ELEMENT
        ix=ix-nblo
        write(lout,10020) "CADCUM> ",ientry,ix,bez(ix),dcum(ientry)
      else
        ! BLOC
        write(lout,10020) "CADCUM> ",ientry,ix,bezb(ix),dcum(ientry)
      end if
    end do
    write(lout,10020) "CADCUM> ",iu+1,-1,"END"//repeat(" ",45),dcum(iu+1)
  end if
  write(lout,"(a,f17.10,a)") "CADCUM> Machine length was ",dcum(iu+1)," [m]"

  return

10020 format(a,2(1x,i6),1x,a48,1x,f12.5)
10030 format(a,2(1x,a6),1x,a48,1x,a12)

end subroutine cadcum

subroutine blocksv
  use floatPrecision
  use numerical_constants

  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  implicit none

  integer ia, ikk, j, jm, k, lkk, mkk
  real(kind=fPrec) dpoff
!     integer i,itiono,i1,i2,i3,ia,ia2,iar,iation,ib,ib0,ib1,ib2,ib3,id,&
!    &idate,ie,ig,ii,ikk,im,imonth,iposc,irecuin,itime,ix,izu,j,j2,jj,  &
!    &jm,k,kpz,kzz,l,lkk,ll,m,mkk,napxto,ncorruo,ncrr,nd,nd2,ndafi2,    &
!    &nerror,nlino,nlinoo,nmz,nthinerr
!     double precision alf0s1,alf0s2,alf0s3,alf0x2,alf0x3,alf0z2,alf0z3,&
!    &amp00,bet0s1,bet0s2,bet0s3,bet0x2,bet0x3,bet0z2,bet0z3,chi,coc,   &
!    &dam1,dchi,ddp1,dp0,dp00,dp10,dpoff,dpsic,dps0,dsign,gam0s1,gam0s2,&
!    &gam0s3,gam0x1,gam0x2,gam0x3,gam0z1,gam0z2,gam0z3,phag,r0,r0a,rat0,&
!    &rdev,rmean,rsqsum,rsum,sic,tasia56,tasiar16,tasiar26,tasiar36,    &
!    &tasiar46,tasiar56,tasiar61,tasiar62,tasiar63,tasiar64,tasiar65,   &
!    &taus,x11,x13

save

#ifdef FLUKA
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 11-06-2014
!     entirely re-initialise to 0.0 hv(...) and bl1v(...) arrays
!     inserted in main code by the 'fluka' compilation flag
  do ia=1,npart
    do k=1,nblo
      do lkk=1,2
        do mkk=1,6
          hv(mkk,lkk,ia,k)=zero
          bl1v(mkk,lkk,ia,k)=zero
        end do
      end do
    end do
  end do
#endif
  do 440 k=1,mblo
    jm=mel(k)
    ikk=mtyp(k,1)
    do lkk=1,2
      do mkk=1,6
        do ia=1,napx
          dpoff=dpsv(ia)*c1e3
          if(abs(dpoff).le.pieni) dpoff=one
          hv(mkk,lkk,ia,1)=al(mkk,lkk,ia,ikk)
          if(mkk.eq.5.or.mkk.eq.6) then
            hv(mkk,lkk,ia,1)=hv(mkk,lkk,ia,1)/dpoff
          end if
        end do
      end do
    end do
    if(jm.eq.1) goto 410
    do j=2,jm
      ikk=mtyp(k,j)
      do lkk=1,2
        do ia=1,napx
          dpoff=dpsv(ia)*c1e3
          if(abs(dpoff).le.pieni) dpoff=one
          hv(1,lkk,ia,j)=hv(1,lkk,ia,j-1)*al(1,lkk,ia,ikk)+ hv(3,lkk,ia,j-1)*al(2,lkk,ia,ikk)
          hv(2,lkk,ia,j)=hv(2,lkk,ia,j-1)*al(1,lkk,ia,ikk)+ hv(4,lkk,ia,j-1)*al(2,lkk,ia,ikk)
          hv(3,lkk,ia,j)=hv(1,lkk,ia,j-1)*al(3,lkk,ia,ikk)+ hv(3,lkk,ia,j-1)*al(4,lkk,ia,ikk)
          hv(4,lkk,ia,j)=hv(2,lkk,ia,j-1)*al(3,lkk,ia,ikk)+ hv(4,lkk,ia,j-1)*al(4,lkk,ia,ikk)
!hr05         hv(5,lkk,ia,j)=hv(5,lkk,ia,j-1)*al(1,lkk,ia,ikk)+ hv(6,   &
!hr05&lkk,ia,j-1)*al(2,lkk,ia,ikk)+al(5,lkk,ia,ikk)/dpoff
          hv(5,lkk,ia,j)=(hv(5,lkk,ia,j-1)*al(1,lkk,ia,ikk)+ hv(6,lkk,ia,j-1)*al(2,lkk,ia,ikk))+al(5,lkk,ia,ikk)/dpoff
!hr05         hv(6,lkk,ia,j)=hv(5,lkk,ia,j-1)*al(3,lkk,ia,ikk)+ hv(6,   &
!hr05&lkk,ia,j-1)*al(4,lkk,ia,ikk)+al(6,lkk,ia,ikk)/dpoff
          hv(6,lkk,ia,j)=(hv(5,lkk,ia,j-1)*al(3,lkk,ia,ikk)+ hv(6,lkk,ia,j-1)*al(4,lkk,ia,ikk))+al(6,lkk,ia,ikk)/dpoff
        end do
      end do
    end do
410   do lkk=1,2
      do mkk=1,6
        do ia=1,napx
          bl1v(mkk,lkk,ia,k)=hv(mkk,lkk,ia,jm)
        end do
      end do
    end do
440 continue
end subroutine blocksv

!-----------------------------------------------------------------------
!  CALCULATION OF : MOMENTUM-DEPENDING ELEMENT-MATRICES AND
!                   CHANGE OF PATH LENGTHS FOR EACH PARTICLE.
!  CAUTION:
!          A SPECIAL VERSION FOR VECTORIZATION - AUGUST   1994
!-----------------------------------------------------------------------
subroutine envarsv(dpsv,oidpsv,rvv,ekv)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer

  use parpro
  use mod_common
  use mod_commons
  use mod_common_track
  use mod_common_da

  use mod_alloc

  implicit none
  integer ih1,ih2,j,kz1,l,l1,l2

  !Local version of variables normally found in mod_common_main
  real(kind=fPrec) aek,afok,as3,as4,as6,co,dpd,dpsq,dpsv,fi,    &
        fok,fok1,fokqv,g,gl,hc,hi,hi1,hm,hp,hs,oidpsv,rho,rhoc,rhoi, &
        rvv,si,siq,sm1,sm12,sm2,sm23,sm3,wf,wfa,wfhi

  real(kind=fPrec), allocatable, intent(inout) :: ekv(:,:) !(npart,nele)

  dimension fokqv(npart),dpsv(npart)
  dimension rvv(npart),oidpsv(npart)
  dimension dpd(npart),dpsq(npart),fok(npart),rho(npart)
  dimension fok1(npart),si(npart),co(npart),g(npart),gl(npart)
  dimension sm1(npart),sm2(npart),sm3(npart),sm12(npart)
  dimension as3(npart),as4(npart),as6(npart),sm23(npart)
  dimension rhoc(npart),siq(npart),aek(npart),afok(npart)
  dimension hp(npart),hm(npart),hc(npart),hs(npart),wf(npart)
  dimension wfa(npart),wfhi(npart),rhoi(npart)
  dimension hi(npart),fi(npart),hi1(npart)

  real(kind=fPrec) fokm

!-----------------------------------------------------------------------
  save
!-----------------------------------------------------------------------

  do 10 j=1,napx
    dpd(j)=one+dpsv(j)
    dpsq(j)=sqrt(dpd(j))
10 continue
  do 160 l=1,il
    do l1=1,6
      do j=1,napx
        do l2=1,2
          al(l1,l2,j,l)=zero
          as(l1,l2,j,l)=zero
        enddo
      enddo
    enddo
    if(abs(el(l)).le.pieni) goto 160
    kz1=kz(l)+1
!       goto(20,40,80,60,40,60,100,100,140),kz1
    if (kz1.eq.1) goto 20
    if (kz1.eq.2) goto 40
    if (kz1.eq.3) goto 80
    if (kz1.eq.4) goto 60
    if (kz1.eq.5) goto 40
    if (kz1.eq.6) goto 60
    if (kz1.eq.7) goto 100
    if (kz1.eq.8) goto 100
    if (kz1.eq.9) goto 140
    goto 160
!-----------------------------------------------------------------------
!  DRIFTLENGTH
!-----------------------------------------------------------------------
20   do 30 j=1,napx
      al(1,1,j,l)=one
      al(1,2,j,l)=one
      al(2,1,j,l)=el(l)
      al(2,2,j,l)=el(l)
      al(3,1,j,l)=zero
      al(3,2,j,l)=zero
      al(4,1,j,l)=one
      al(4,2,j,l)=one
      as(6,1,j,l)=((-one*rvv(j))*el(l))/c2e3                         !hr06
      as(6,2,j,l)=as(6,1,j,l)
      as(1,1,j,l)=(el(l)*(one-rvv(j)))*c1e3                          !hr06
30   continue
    goto 160
!-----------------------------------------------------------------------
!  RECTANGULAR MAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
40   fokm=el(l)*ed(l)
    if(abs(fokm).le.pieni) goto 20
    if(kz1.eq.2) then
      ih1=1
      ih2=2
    else
!  RECTANGULAR MAGNET VERTICAL
      ih1=2
      ih2=1
    endif
    do 50 j=1,napx
      fok(j)=fokm/dpsq(j)
      rho(j)=(one/ed(l))*dpsq(j)
      fok1(j)=(tan_mb(fok(j)*half))/rho(j)
      si(j)=sin_mb(fok(j))
      co(j)=cos_mb(fok(j))
      al(1,ih1,j,l)=one
      al(2,ih1,j,l)=rho(j)*si(j)
      al(3,ih1,j,l)=zero
      al(4,ih1,j,l)=one
  al(5,ih1,j,l)=((-one*dpsv(j))*((rho(j)*(one-co(j)))/dpsq(j)))*c1e3 !hr06
  al(6,ih1,j,l)=((-one*dpsv(j))*((two*tan_mb(fok(j)*half))/dpsq(j)))&!hr06
  &*c1e3                                                              !hr06
      sm1(j)=cos_mb(fok(j))
      sm2(j)=sin_mb(fok(j))*rho(j)
      sm3(j)=(-one*sin_mb(fok(j)))/rho(j)                            !hr06
      sm12(j)=el(l)-sm1(j)*sm2(j)
      sm23(j)=sm2(j)*sm3(j)
      as3(j)=(-one*rvv(j))*(((dpsv(j)*rho(j))/(two*dpsq(j)))*sm23(j)-&!hr06
  &(rho(j)*dpsq(j))*(one-sm1(j)))                                     !hr06
      as4(j)=((-one*rvv(j))*sm23(j))/c2e3                            !hr06
      as6(j)=((-one*rvv(j))*(el(l)+sm1(j)*sm2(j)))/c4e3              !hr06
  as(1,ih1,j,l)=(el(l)*(one-rvv(j))-rvv(j)*((dpsv(j)**2/            &!hr06
  &(four*dpd(j)))*sm12(j)+dpsv(j)*(el(l)-sm2(j))))*c1e3               !hr06
  as(2,ih1,j,l)=(-one*rvv(j))*((dpsv(j)/((two*rho(j))*dpsq(j)))*    &!hr06
  &sm12(j)-(sm2(j)*dpsq(j))/rho(j))+fok1(j)*as3(j)                    !hr06
      as(3,ih1,j,l)=as3(j)
      as(4,ih1,j,l)=as4(j)+(two*as6(j))*fok1(j)                      !hr06
      as(5,ih1,j,l)=(as6(j)*fok1(j)**2                              &!hr06
  &-(rvv(j)*sm12(j))/(c4e3*rho(j)**2))+fok1(j)*as4(j)                 !hr06
      as(6,ih1,j,l)=as6(j)
!--VERTIKAL
      g(j)=tan_mb(fok(j)*half)/rho(j)
      gl(j)=el(l)*g(j)
      al(1,ih2,j,l)=one-gl(j)
      al(2,ih2,j,l)=el(l)
      al(3,ih2,j,l)=(-one*g(j))*(two-gl(j))                          !hr06
      al(4,ih2,j,l)=al(1,ih2,j,l)
      as6(j)=((-one*rvv(j))*al(2,ih2,j,l))/c2e3                      !hr06
      as(4,ih2,j,l)=((-one*two)*as6(j))*fok1(j)                      !hr06
      as(5,ih2,j,l)=as6(j)*fok1(j)**2                                !hr06
      as(6,ih2,j,l)=as6(j)
50   continue
    goto 160
!-----------------------------------------------------------------------
!  SEKTORMAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
60   fokm=el(l)*ed(l)
    if(abs(fokm).le.pieni) goto 20
    if(kz1.eq.4) then
      ih1=1
      ih2=2
    else
!  SECTOR MAGNET VERTICAL
      ih1=2
      ih2=1
    endif
    do 70 j=1,napx
      fok(j)=fokm/dpsq(j)
      rho(j)=(one/ed(l))*dpsq(j)
      si(j)=sin_mb(fok(j))
      co(j)=cos_mb(fok(j))
      rhoc(j)=(rho(j)*(one-co(j)))/dpsq(j)                           !hr06
      siq(j)=si(j)/dpsq(j)
      al(1,ih1,j,l)=co(j)
      al(2,ih1,j,l)=rho(j)*si(j)
      al(3,ih1,j,l)=(-one*si(j))/rho(j)                              !hr06
      al(4,ih1,j,l)=co(j)
      al(5,ih1,j,l)=((-one*dpsv(j))*rhoc(j))*c1e3                    !hr06
      al(6,ih1,j,l)=((-one*dpsv(j))*siq(j))*c1e3                     !hr06
      sm12(j)=el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
      sm23(j)=al(2,ih1,j,l)*al(3,ih1,j,l)
    as(1,ih1,j,l)=(el(l)*(one-rvv(j))-rvv(j)*((dpsv(j)**2/           &!hr06
  &(four*dpd(j)))*sm12(j)+dpsv(j)*(el(l)-al(2,ih1,j,l))))*c1e3        !hr06
    as(2,ih1,j,l)=(-one*rvv(j))*((dpsv(j)/((two*rho(j))*dpsq(j)))*   &!hr06
  &sm12(j)-dpd(j)*siq(j))                                             !hr06
      as(3,ih1,j,l)=(-one*rvv(j))*(((dpsv(j)*rho(j))/(two*dpsq(j)))* &!hr06
  &sm23(j)-dpd(j)*rhoc(j))                                            !hr06
      as(4,ih1,j,l)=((-one*rvv(j))*sm23(j))/c2e3                     !hr06
      as(5,ih1,j,l)=((-one*rvv(j))*sm12(j))/(c4e3*rho(j)**2)         !hr06
  as(6,ih1,j,l)=((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l)))/&!hr06
  &c4e3                                                               !hr06
!--VERTIKAL
      al(1,ih2,j,l)=one
      al(2,ih2,j,l)=el(l)
      al(3,ih2,j,l)=zero
      al(4,ih2,j,l)=one
      as(6,ih2,j,l)=((-one*rvv(j))*al(2,ih2,j,l))/c2e3               !hr06
70   continue
    goto 160
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
80   do 90 j=1,napx
      fok(j)=ekv(j,l)*oidpsv(j)
      aek(j)=abs(fok(j))
      hi(j)=sqrt(aek(j))
      fi(j)=el(l)*hi(j)
      if(fok(j).le.zero) then
        al(1,1,j,l)=cos_mb(fi(j))
        hi1(j)=sin_mb(fi(j))
        if(abs(hi(j)).le.pieni) then
          al(2,1,j,l)=el(l)
        else
          al(2,1,j,l)=hi1(j)/hi(j)
        endif
        al(3,1,j,l)=(-one*hi1(j))*hi(j)                              !hr06
        al(4,1,j,l)=al(1,1,j,l)
        as(1,1,j,l)=(el(l)*(one-rvv(j)))*c1e3                        !hr06
        as(4,1,j,l)=(((-one*rvv(j))*al(2,1,j,l))*al(3,1,j,l))/c2e3   !hr06
      as(5,1,j,l)=(((-one*rvv(j))*(el(l)-al(1,1,j,l)*al(2,1,j,l)))* &!hr06
  &aek(j))/c4e3                                                       !hr06
  as(6,1,j,l)=((-one*rvv(j))*(el(l)+al(1,1,j,l)*al(2,1,j,l)))/c4e3   !hr06
!--DEFOCUSSING
        hp(j)=exp_mb(fi(j))
        hm(j)=one/hp(j)
        hc(j)=(hp(j)+hm(j))*half
        hs(j)=(hp(j)-hm(j))*half
        al(1,2,j,l)=hc(j)
        if(abs(hi(j)).le.pieni) then
          al(2,2,j,l)=el(l)
        else
          al(2,2,j,l)=hs(j)/hi(j)
        endif
        al(3,2,j,l)=hs(j)*hi(j)
        al(4,2,j,l)=hc(j)
        as(4,2,j,l)=((-one*rvv(j))*al(2,2,j,l)*al(3,2,j,l))/c2e3     !hr06
      as(5,2,j,l)=((rvv(j)*(el(l)-al(1,2,j,l)*al(2,2,j,l)))*aek(j)) &!hr06
  &/c4e3                                                              !hr06
  as(6,2,j,l)=((-one*rvv(j))*(el(l)+al(1,2,j,l)*al(2,2,j,l)))/c4e3   !hr06
      else
        al(1,2,j,l)=cos_mb(fi(j))
        hi1(j)=sin_mb(fi(j))
        if(abs(hi(j)).le.pieni) then
          al(2,2,j,l)=el(l)
        else
          al(2,2,j,l)=hi1(j)/hi(j)
        endif
        al(3,2,j,l)=(-one*hi1(j))*hi(j)                              !hr06
        al(4,2,j,l)=al(1,2,j,l)
        as(1,2,j,l)=(el(l)*(one-rvv(j)))*c1e3                        !hr06
        as(4,2,j,l)=(((-one*rvv(j))*al(2,2,j,l))*al(3,2,j,l))/c2e3   !hr06
      as(5,2,j,l)=(((-one*rvv(j))*(el(l)-al(1,2,j,l)*al(2,2,j,l)))* &!hr06
  &aek(j))/c4e3                                                       !hr06
  as(6,2,j,l)=((-one*rvv(j))*(el(l)+al(1,2,j,l)*al(2,2,j,l)))/c4e3   !hr06
!--DEFOCUSSING
        hp(j)=exp_mb(fi(j))
        hm(j)=one/hp(j)
        hc(j)=(hp(j)+hm(j))*half
        hs(j)=(hp(j)-hm(j))*half
        al(1,1,j,l)=hc(j)
        if(abs(hi(j)).le.pieni) then
          al(2,1,j,l)=el(l)
        else
          al(2,1,j,l)=hs(j)/hi(j)
        endif
        al(3,1,j,l)=hs(j)*hi(j)
        al(4,1,j,l)=hc(j)
        as(4,1,j,l)=(((-one*rvv(j))*al(2,1,j,l))*al(3,1,j,l))/c2e3   !hr06
      as(5,1,j,l)=((rvv(j)*(el(l)-al(1,1,j,l)*al(2,1,j,l)))*aek(j)) &!hr06
  &/c4e3                                                              !hr06
  as(6,1,j,l)=((-one*rvv(j))*(el(l)+al(1,1,j,l)*al(2,1,j,l)))/c4e3   !hr06
      endif
90   continue
    goto 160
!-----------------------------------------------------------------------
!  COMBINED FUNCTION MAGNET HORIZONTAL
!  FOCUSSING
!-----------------------------------------------------------------------
100   if(kz1.eq.7) then
      do 110 j=1,napx
        fokqv(j)=ekv(j,l)
110     continue
      ih1=1
      ih2=2
    else
!  COMBINED FUNCTION MAGNET VERTICAL
      do 120 j=1,napx
        fokqv(j)=-one*ekv(j,l)                                       !hr06
120     continue
      ih1=2
      ih2=1
    endif
    do 130 j=1,napx
      wf(j)=ed(l)/dpsq(j)
      fok(j)=fokqv(j)/dpd(j)-wf(j)**2                                !hr06
      afok(j)=abs(fok(j))
      hi(j)=sqrt(afok(j))
      fi(j)=hi(j)*el(l)
      if(afok(j).le.pieni) then
        al(1,1,j,l)=one
        al(1,2,j,l)=one
        al(2,1,j,l)=el(l)
        al(2,2,j,l)=el(l)
        al(3,1,j,l)=zero
        al(3,2,j,l)=zero
        al(4,1,j,l)=one
        al(4,2,j,l)=one
        as(6,1,j,l)=((-one*rvv(j))*el(l))/c2e3                       !hr06
        as(6,2,j,l)=as(6,1,j,l)
        as(1,1,j,l)=(el(l)*(one-rvv(j)))*c1e3                        !hr06
      endif
      if(fok(j).lt.(-one*pieni)) then                                !hr06
        si(j)=sin_mb(fi(j))
        co(j)=cos_mb(fi(j))
        wfa(j)=((wf(j)/afok(j))*(one-co(j)))/dpsq(j)                 !hr06
        wfhi(j)=((wf(j)/hi(j))*si(j))/dpsq(j)                        !hr06
        al(1,ih1,j,l)=co(j)
        al(2,ih1,j,l)=si(j)/hi(j)
        al(3,ih1,j,l)=(-one*si(j))*hi(j)                             !hr06
        al(4,ih1,j,l)=co(j)
        al(5,ih1,j,l)=((-one*wfa(j))*dpsv(j))*c1e3                   !hr06
        al(6,ih1,j,l)=((-one*wfhi(j))*dpsv(j))*c1e3                  !hr06
        sm12(j)=el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
        sm23(j)=al(2,ih1,j,l)*al(3,ih1,j,l)
        as(1,ih1,j,l)=(el(l)*(one-rvv(j))-                          &!hr06
  &((rvv(j)*((dpsv(j)**2/(four*dpd(j)))*                             &!hr06
  &sm12(j)+ dpsv(j)*(el(l)-al(2,ih1,j,l))))/afok(j))*wf(j)**2)*c1e3   !hr06
  as(2,ih1,j,l)=(-one*rvv(j))*(((dpsv(j)*wf(j))/(two*dpsq(j)))*     &!hr06
  &sm12(j)-dpd(j)*wfhi(j))                                            !hr06
  as(3,ih1,j,l)=(-one*rvv(j))*(((((dpsv(j)*half)/afok(j))/dpd(j))*  &!hr06
  &ed(l))*sm23(j)-dpd(j)*wfa(j))                                      !hr06
        as(4,ih1,j,l)=((-one*rvv(j))*sm23(j))/c2e3                   !hr06
        as(5,ih1,j,l)=(((-one*rvv(j))*sm12(j))*afok(j))/c4e3         !hr06
  as(6,ih1,j,l)=((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l))) &!hr06
  &/c4e3                                                              !hr06
        aek(j)=abs(ekv(j,l)/dpd(j))
        hi(j)=sqrt(aek(j))
        fi(j)=hi(j)*el(l)
        hp(j)=exp_mb(fi(j))
        hm(j)=one/hp(j)
        hc(j)=(hp(j)+hm(j))*half
        hs(j)=(hp(j)-hm(j))*half
        al(1,ih2,j,l)=hc(j)
        al(2,ih2,j,l)=el(l)
        if(abs(hi(j)).gt.pieni) al(2,ih2,j,l)=hs(j)/hi(j)
        al(3,ih2,j,l)=hs(j)*hi(j)
        al(4,ih2,j,l)=hc(j)
  as(4,ih2,j,l)=(((-one*rvv(j))*al(2,ih2,j,l))*al(3,ih2,j,l))/c2e3   !hr06
      as(5,ih2,j,l)=((rvv(j)*(el(l)-al(1,ih2,j,l)*al(2,ih2,j,l)))*  &!hr06
  &aek(j))/c4e3                                                       !hr06
  as(6,ih2,j,l)=((-one*rvv(j))*(el(l)+al(1,ih2,j,l)*al(2,ih2,j,l))) &!hr06
  &/c4e3                                                              !hr06
      endif
!--DEFOCUSSING
      if(fok(j).gt.pieni) then
        hp(j)=exp_mb(fi(j))
        hm(j)=one/hp(j)
        hc(j)=(hp(j)+hm(j))*half
        hs(j)=(hp(j)-hm(j))*half
        al(1,ih1,j,l)=hc(j)
        al(2,ih1,j,l)=hs(j)/hi(j)
        al(3,ih1,j,l)=hs(j)*hi(j)
        al(4,ih1,j,l)=hc(j)
        wfa(j)=((wf(j)/afok(j))*(one-hc(j)))/dpsq(j)                 !hr06
        wfhi(j)=((wf(j)/hi(j))*hs(j))/dpsq(j)                        !hr06
        al(5,ih1,j,l)= (wfa(j)*dpsv(j))*c1e3                         !hr06
        al(6,ih1,j,l)=((-one*wfhi(j))*dpsv(j))*c1e3                  !hr06
        sm12(j)=el(l)-al(1,ih1,j,l)*al(2,ih1,j,l)
        sm23(j)=al(2,ih1,j,l)*al(3,ih1,j,l)
        as(1,ih1,j,l)=(((rvv(j)*((dpsv(j)**2/(four*dpd(j)))*sm12(j) &
  &+dpsv(j)*(el(l)-al(2,ih1,j,l))))/afok(j))*wf(j)**2+el(l)*         &
  &(one-rvv(j)))*c1e3
      as(2,ih1,j,l)=(-one*rvv(j))*(((dpsv(j)*wf(j))/(two*dpsq(j)))* &!hr06
  &sm12(j)-dpd(j)*wfhi(j))                                            !hr06
    as(3,ih1,j,l)=rvv(j)*(((((dpsv(j)*half)/afok(j))/dpd(j))*ed(l)) &!hr06
  &*sm23(j)-dpd(j)*wfa(j))                                            !hr06
        as(4,ih1,j,l)=((-one*rvv(j))*sm23(j))/c2e3                   !hr06
        as(5,ih1,j,l)=((rvv(j)*sm12(j))*afok(j))/c4e3                !hr06
  as(6,ih1,j,l)=((-one*rvv(j))*(el(l)+al(1,ih1,j,l)*al(2,ih1,j,l))) &!hr06
  &/c4e3                                                              !hr06
        aek(j)=abs(ekv(j,l)/dpd(j))
        hi(j)=sqrt(aek(j))
        fi(j)=hi(j)*el(l)
        si(j)=sin_mb(fi(j))
        co(j)=cos_mb(fi(j))
        al(1,ih2,j,l)=co(j)
        al(2,ih2,j,l)=si(j)/hi(j)
        al(3,ih2,j,l)=(-one*si(j))*hi(j)                             !hr06
        al(4,ih2,j,l)=co(j)
  as(4,ih2,j,l)=(((-one*rvv(j))*al(2,ih2,j,l))*al(3,ih2,j,l))/c2e3   !hr06
  as(5,ih2,j,l)=(((-one*rvv(j))*(el(l)-al(1,ih2,j,l)*al(2,ih2,j,l)))&!hr06
  &*aek(j))/c4e3                                                      !hr06
  as(6,ih2,j,l)=((-one*rvv(j))*(el(l)+al(1,ih2,j,l)*al(2,ih2,j,l))) &!hr06
  &/c4e3                                                              !hr06
      endif
130   continue
    goto 160
!-----------------------------------------------------------------------
!  EDGE FOCUSSING
!-----------------------------------------------------------------------
140   do 150 j=1,napx
      rhoi(j)=ed(l)/dpsq(j)
      fok(j)=rhoi(j)*tan_mb((el(l)*rhoi(j))*half)                    !hr06
      al(1,1,j,l)=one
      al(2,1,j,l)=zero
      al(3,1,j,l)=fok(j)
      al(4,1,j,l)=one
      al(1,2,j,l)=one
      al(2,2,j,l)=zero
      al(3,2,j,l)=-fok(j)
      al(4,2,j,l)=one
150   continue
160 continue

  return
end subroutine envarsv

!-----------------------------------------------------------------------
!  CALCULATION OF THE 4-DIMENSIONAL CLOSED ORBIT INCLUDING DELTA
!-----------------------------------------------------------------------
subroutine mydaini(ncase,nnord,nnvar,nndim,nnvar2,nnord1)

  use floatPrecision
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common, only : ichromc,ilinc,iqmodc
  use mod_common_da
  use mod_lie_dab, only : iscrda,mld_allocArrays
  implicit none
  integer idummy,ncase,ndimfo,ndpt,nis,nndim,nnord,nnord1,nnvar,nnvar2,nord1o,nordo,nvar2o,nvaro
  real(kind=fPrec) am
  dimension am(6,6),idummy(6)
  save
!-----------------------------------------------------------------------
  if(nndim < 2 .or. nndim > 3) then
    write(lout,"(a)") "DAINI> ERROR DA corrections implemented for 4D and 6D only."
    call prror(-1)
  end if
!--------------------
  nordo=nord
  nvaro=nvar
  ndimfo=ndimf
  nvar2o=nvar2
  nord1o=nord1
!--------------------
  nord=nnord
  nvar=nnvar
  ndimf=nndim
  nvar2=nnvar2
  nord1=nnord1
!--------------------
  ndpt=0
  nis=0
!--------------------
  call daeps(preda)
  call idprset(-102)
  call mld_allocArrays(.false.)
  call lieinit(nord,nvar,ndimf,ndpt,0,nis)
#ifdef DEBUG
!     call dumpbin('alieinit',1,11)
!     call abend('alieinit in mydaini                               ')
#endif
  write(lout,10000) nord,nvar,ndimf
  call daall(iscrda,100,'$$IS      ',nord,nvar)
!--closed orbit
#ifdef DEBUG
!     write(*,*) 'ncase=',ncase,' if 1 call clorda'
#endif
  if(ncase.eq.1) call clorda(2*ndimf,idummy,am)
#ifdef DEBUG
!     call dumpbin('aclorda',1,11)
!     call abend('aclorda                                           ')
#endif
!--tune variation
#ifdef DEBUG
!     write(*,*) 'ncase=',ncase,' if 2 call umlauda'
#endif
  if(ncase.eq.2) call umlauda
#ifdef DEBUG
!     if(ncase.eq.2) then
!     call dumpbin('aumlauda',7,77)
!     call abend('aumlauda                                          ')
!     endif
#endif
  iqmodc=0
  ichromc=0
  ilinc=0
  call dadal(iscrda,100)
!--------------------
  nord=nordo
  nvar=nvaro
  nvar2=nvar2o
  ndimf=ndimfo
  nord1=nord1o
!-----------------------------------------------------------------------
10000 format(/131('-')/10x,'DA INITIALIZATION: ORDER = ',i2,            &
  &', # of VARIABLES = ',i2,', DIMENSION = ',i2/)
  return
end subroutine mydaini
