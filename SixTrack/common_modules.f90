!     This file contains modules with variables that are used in all
!     parts of the code, such as lout and crlibm functions, accelerator structure,
!     particle arrays etc.

module parpro
  ! Shared module defining the sizes of key arrays
  ! This is (will be) the "access point" for resizing
  ! the size of the accelerator, the number of particles, etc.
  
  implicit none
  
  integer mbea,mcor,mcop,mmul,mpa,mran,nbb,nblz,ncom,ncor1, &
       nelb,nema,ninv,nlya,nmac,nmon1,npart,nper,nplo,npos,nran, &
       nrco,ntr,nzfz

  !Maximum length of element names
  integer, parameter :: max_name_len = 48
  
  !Max number of particles
#if !defined(BIGNPART) && !defined(HUGENPART)
  parameter(npart = 64,nmac = 1)
#endif
#if defined(BIGNPART) && !defined(HUGENPART)
  !See also:
  ! - subroutine wzsubv
  parameter(npart = 2048,nmac = 1)
#endif
#if !defined(BIGNPART) && defined(HUGENPART)
  !See also:
  ! - subroutine wzsubv
  parameter(npart = 65536,nmac = 1)
#endif
  
  integer, parameter :: nele_initial = 50   ! Must be at least 1
  integer, parameter :: nblo_initial = 2500 ! Setting to 2500 (old value) for now as it blows up when scaling.
  integer :: nele = -1
  integer :: nblo = -1

!Note: nzfz should be = 3*nblz+2*mmul*#MULTIPOLES,
! where #MULTIPOLES are the max number of multipoles in the lattice (up to nblz)
! For now, scale the number of multipoles (from nzfz) as is done in the "no-flag" version:
! 6000/20000 -> 30% multipoles
#ifndef COLLIMAT
#ifdef BIGNBLZ
  parameter(nper=16,nelb=140,nblz=200000,nzfz = 3000000,mmul = 20) !up to 60'000 multipoles
#endif
#ifdef HUGENBLZ
  parameter(nper=16,nelb=280,nblz=400000,nzfz = 6000000,mmul = 20) !up to 120'000 multipoles -> 48MB/nzfz-array (20%)
#endif
#if !defined(BIGNBLZ) && !defined(HUGENBLZ)
  parameter(nper=16,nelb=140,nblz=20000,nzfz = 300000,mmul = 20) !up to 6'000 multipoles
#endif
#else
#ifdef BEAMGAS
  parameter(nper=16,nelb=140,nblz=200000,nzfz = 1920000,mmul = 11) !up to 60'000 multipoles
#else
#ifdef BIGNBLZ
  parameter(nper=16,nelb=140,nblz=200000,nzfz = 1920000,mmul = 11) !up to 60'000 multipoles
#endif
#ifdef HUGENBLZ
  parameter(nper=16,nelb=140,nblz=400000,nzfz = 3840000,mmul = 11) !up to 120'000 multipoles (20%)
#endif
#if !defined(BIGNBLZ) && !defined(HUGENBLZ)
  parameter(nper=16,nelb=140,nblz=15000,nzfz = 144000,mmul = 11) !up to 4500 multipoles
#endif
#endif
#endif
  
  parameter(nran = 2000000,ncom = 100, mran = 500, mpa = 6, nrco = 5, nema = 15)
  parameter(mcor = 10,mcop = mcor+6, mbea = 99)
  parameter(npos = 20000,nlya = 10000,ninv = 1000,nplo = 20000)
  parameter(nmon1 = 600,ncor1 = 600)
  parameter(ntr = 20)
  
  ! Beam-beam lenses
#if !defined(BIGNBLZ) && !defined(HUGENBLZ)
  parameter(nbb = 350)
#endif
#if defined(BIGNBLZ) || defined(HUGENBLZ)
  parameter(nbb = 500)
#endif
  
end module parpro

module mod_common
  use parpro
  use floatPrecision
  use mod_alloc
  use numerical_constants
  
  implicit none

  ! common /erro/
  integer, save :: ierro
  character(len=max_name_len), save :: erbez
  
  ! common /kons/
  real(kind=fPrec), save :: pi,pi2,pisqrt,rad
  
  ! common /str/
  integer, save :: il,mper,mblo,mbloz,msym(nper),kanf,iu,ic(nblz)
  
  ! common /ell/
  real(kind=fPrec), allocatable, save :: ed(:),el(:),ek(:),sm(:) ! (nele)
  integer, allocatable, save :: kz(:),kp(:)  ! (nele)
  
  ! common /bbb/
  real(kind=fPrec), allocatable, save :: bbbx(:), bbby(:), bbbs(:) ! (nele)
  
  ! common /pla/
  real(kind=fPrec), allocatable, save :: xpl(:),xrms(:),zpl(:),zrms(:) ! (nele)
  
  !common /str2/
  integer, allocatable, save :: mel(:)    ! (nblo)
  integer, allocatable, save :: mtyp(:,:) ! (nblo,nelb)
  integer, allocatable, save :: mstr(:)   ! (nblo)

  ! common /mat/
  real(kind=fPrec), allocatable, save :: a(:,:,:)   ! (nele,2,6)
  real(kind=fPrec), allocatable, save :: bl1(:,:,:) ! (nblo,2,6)
  real(kind=fPrec), allocatable, save :: bl2(:,:,:) ! (nblo,2,6)
  
  ! common /syos2/
  real(kind=fPrec), save :: rvf(mpa)
  
  ! common /tra1/
  real(kind=fPrec), save :: rat
  integer, save :: idfor, napx, napxo, numl, niu(2), numlr, nde(2), nwr(4), &
       ird, imc, irew, ntwin, iclo6, iclo6r, iver, ibidu, numlcp, numlmax, nnuml
  
  ! common /syn/
  real(kind=fPrec), save :: qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3), &
       crad,dppoff,sigmoff(nblz),tlen,mtcda
  real(kind=fPrec), allocatable, save ::  hsyc(:), phasc(:) !(nele)
  integer, save :: iicav,ition,idp,ncy,ixcav
  integer, allocatable, save :: itionc(:) !(nele)
  
  ! common /corcom/
  real(kind=fPrec), save :: dpscor,sigcor
  integer, save :: icode,idam,its6d

  ! common /multi/
  real(kind=fPrec), allocatable, save :: bk0(:,:),ak0(:,:), bka(:,:),aka(:,:) !(nele,mmul)

  ! common /mult1/
  real(kind=fPrec), save :: benki
  real(kind=fPrec), allocatable, save :: benkc(:),r00(:) !(nele)
  integer, allocatable, save :: irm(:),nmu(:) !(nele)
  
  ! common /rand0/
  real(kind=fPrec), save :: zfz(nzfz)
  integer, save :: iorg,mzu(nblz),izu0,mmac,mcut
  character(len=:), allocatable, save :: bezr(:,:) !(max_name_len)(3,nele)
  
  ! common /rand1/
  real(kind=fPrec), save :: exterr(nblz,40),extalign(nblz,3),tiltc(nblz),tilts(nblz)
  integer, save :: mout2,icext(nblz),icextal(nblz)
  
  ! common /beo/
  real(kind=fPrec), save :: aper(2),di0(2),dip0(2),ta(6,6)
  
  ! common /clo/
  real(kind=fPrec), save :: dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0
  integer, save :: itco,itcro,itqv,iout,iquiet
  
  ! common /qmodi/
  real(kind=fPrec), save :: qw0(3),amp0
  integer, save :: iq(3),iqmod,iqmod6
  integer, allocatable, save :: kpa(:) !(nele)
  
  ! common /linop/
  character(len=:), allocatable, save :: bez(:), bezl(:) ! (nele)
  character(len=:), allocatable, save :: bezb(:)         ! (nblo)
  real(kind=fPrec), allocatable, save :: elbe(:)         ! (nblo)
  real(kind=fPrec),              save :: eui,euii
  integer,                       save :: ilin,nt,iprint,ntco,nlin
  
  ! common /cororb/
  real(kind=fPrec), save :: betam(nmon1,2),pam(nmon1,2),betac(ncor1,2), &
       pac(ncor1,2),bclorb(nmon1,2)
  integer, save :: nhmoni,nhcorr,nvmoni,nvcorr
  integer, allocatable, save :: ncororb(:) !nele
  
  ! common /clos/
  real(kind=fPrec), save :: sigma0(2)
  integer, save :: iclo, ncorru,ncorrep
  
  ! common /combin/
  real(kind=fPrec), save :: ratio(ncom,20)
  real(kind=fPrec), allocatable, save :: ratioe(:) !(nele)
  integer, save :: icomb0(20),icomb(ncom,20),icoe
  integer, allocatable, save :: iratioe(:) !(nele)
  
  ! common/seacom/m21,m22,m23
  integer, save :: ise,ise1,ise2,ise3,mesa,mp,m21,m22,m23
  integer, allocatable, save :: isea(:) !(nele)
  
  ! common /subres/
  real(kind=fPrec), save :: qxt,qzt,tam1,tam2,totl
  integer, save :: isub,nta,nte,ipt
  
  ! common /secom/
  real(kind=fPrec), save :: rtc(9,18,10,5),rts(9,18,10,5)
  integer, save :: ire(12),ipr(5),irmod2
  
  ! common /secom1/
  real(kind=fPrec), save :: dtr(10)
  integer, save :: nre,nur,nch,nqc,npp,nrr(5),nu(5)
  
  ! common /postr/
  real(kind=fPrec), save :: dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2
  integer, save :: nstart,nstop,iskip,iconv,imad
  
  ! common /posti1/
  integer, save :: ipos,iav,iwg,ivox,ivoz,ires,ifh
  character(len=80), save :: toptit(5) !DANGER: If the len changes, CRCHECK will break.
  
  ! common /posti2/
  integer, save :: kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi
  
  ! common /skew/
  real(kind=fPrec), save :: qwsk(2),betx(2),betz(2),alfx(2),alfz(2)
  integer, save :: iskew,nskew(6)
  
  ! common /pawc/
  real, save :: hmal(nplo)
  
  ! common /tit/
  character (len=80), save :: sixtit,commen !DANGER: If the len changes, CRCHECK will break.
  integer, save :: ithick
  
  ! common/co6d/
  real(kind=fPrec), save :: clo6(3),clop6(3)
  
  ! common /dkic/
  real(kind=fPrec), allocatable, save :: dki(:,:) !(nele,3)

  ! common /beam/
  real(kind=fPrec), save :: sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb), &
       clobeam(6,nbb),beamoff(6,nbb),track6d(6,npart),    &
       sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar
  real(kind=fPrec), allocatable, save :: parbe(:,:) !(nele,18)
  real(kind=fPrec), allocatable, save :: ptnfac(:) !(nele)
  integer, save :: nbeam,ibbc,ibeco,ibtyp,lhc
  
  ! common/trom/
  real(kind=fPrec), save :: cotr(ntr,6),rrtr(ntr,6,6)

  integer, allocatable, save ::  imtr(:) !(nele),
  
  ! common /bb6d/
  real(kind=fPrec), save :: bbcu(nbb,12)
  integer, save :: ibb6d,imbb(nblz)
  
  ! common /acdipco/
  real(kind=fPrec), allocatable, save :: acdipph(:) !(nele)
  integer, allocatable, save :: nturn1(:), nturn2(:), nturn3(:), nturn4(:) !(nele)
  
  ! common /crabco/
  real(kind=fPrec), allocatable, save :: crabph(:),crabph2(:),crabph3(:),crabph4(:) !(nele)
  
  ! common /exact/
  integer, save :: iexact
  
  ! common /sixdim/
  real(kind=fPrec), save :: aml6(6,6),edcor(2)
  
  ! common /postr2/
  integer,          save :: nnumxv(npart)
  
  ! common /correct/
  integer,          save :: ichromc,ilinc,iqmodc
  real(kind=fPrec), save :: corr(3,3),chromc(2),wxys(3),clon(6)
  
  ! common /damp/
  real(kind=fPrec), save :: damp,ampt
  
  ! common /ttime/
  integer,          save :: napxto
  real,             save :: tlim,time0,time1,time2,time3,trtime,pretime,posttime,tottime
  
  ! common /xz/
  real(kind=fPrec), save :: xsi(nblz),zsi(nblz),smi(nblz),smizf(nblz),aai(nblz,mmul),bbi(nblz,mmul)
  
#ifdef FLUKA
  ! This was not in a common block, just defined in +cd commonxz as integer fluka_con
  integer,          save :: fluka_con
#endif
  
contains

subroutine mod_common_allocate_arrays
  
  use crcoall
  use string_tools
  
  implicit none
  
  call alloc(ed,nele,zero,'ed')
  call alloc(el,nele,zero,'el')
  call alloc(ek,nele,zero,'ek')
  call alloc(sm,nele,zero,'sm')
  call alloc(kz,nele,0,'kz')
  call alloc(kp,nele,0,'kp')
  call alloc(bbbx,nele,zero,'bbbx')
  call alloc(bbby,nele,zero,'bbby')
  call alloc(bbbs,nele,zero,'bbbs')
  call alloc(xpl,nele,zero,'xpl')
  call alloc(zpl,nele,zero,'zpl')
  call alloc(xrms,nele,zero,'xrms')
  call alloc(zrms,nele,zero,'zrms')
  call alloc(a,nele,2,6,zero,'a')
  call alloc(hsyc,nele,zero,'hsyc')
  call alloc(phasc,nele,zero,'phasc')
  call alloc(itionc,nele,0,'itionc')
  call alloc(bk0,nele,mmul,zero,'bk0')
  call alloc(ak0,nele,mmul,zero,'ak0')
  call alloc(bka,nele,mmul,zero,'bka')
  call alloc(aka,nele,mmul,zero,'aka')
  call alloc(benkc,nele,zero,'benkc')
  call alloc(r00,nele,zero,'r00')
  call alloc(irm,nele,0,'irm')
  call alloc(nmu,nele,0,'nmu')
  call alloc(bezr, max_name_len, 3, nele, str_nmZeros, 'bezr')
  call alloc(kpa,nele,0,'kpa')
  call alloc(bez, max_name_len, nele, str_nmZeros, 'bez')
! call alloc(bezb, max_name_len, nele, ' ', 'bezb')
  call alloc(bezl, max_name_len, nele, str_nmZeros, 'bezl')
  call alloc(ncororb,nele,0,'ncororb')
  call alloc(ratioe,nele,zero,'ratioe')
  call alloc(iratioe,nele,0,'iratioe')
  call alloc(isea,nele,0,'isea')
  call alloc(dki,nele,3,zero,'dki')
  call alloc(parbe,nele,18,zero,'parbe')
  call alloc(ptnfac,nele,zero,'ptnfac')
  call alloc(imtr,nele,0,'imtr')
  call alloc(acdipph,nele,zero,'acdipph')
  call alloc(nturn1,nele,0,'nturn1')
  call alloc(nturn2,nele,0,'nturn2')
  call alloc(nturn3,nele,0,'nturn3')
  call alloc(nturn4,nele,0,'nturn4')
  call alloc(crabph,nele,zero,'crabph')
  call alloc(crabph2,nele,zero,'crabph2')
  call alloc(crabph3,nele,zero,'crabph3')
  call alloc(crabph4,nele,zero,'crabph4')
  
  ! linopt
  call alloc(bezb,max_name_len,nblo,str_nmZeros,"bezb")
  call alloc(elbe,nblo,zero,"elbe")
  
  ! str2
  call alloc(mel,nblo,0,"mel")
  call alloc(mtyp,nblo,nelb,0,"mtyp")
  call alloc(mstr,nblo,0,"mstr")
  
  ! mat
  call alloc(bl1,nblo,2,6,zero,"bl1")
  call alloc(bl2,nblo,2,6,zero,"bl2")
  
  ! TODO: Certain "COMNUL" functionality should probably be done here as well.
  return !Finished successfully
  
end subroutine mod_common_allocate_arrays

subroutine mod_common_expand_arrays(nele_new, nblo_new)
  
  use string_tools
  
  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: nblo_new
  
  !TODO: Do the thing. This will be long and fugly
  call resize(ed,nele_new,zero,'ed')
  call resize(el,nele_new,zero,'el')
  call resize(ek,nele_new,zero,'ek')
  call resize(sm,nele_new,zero,'sm')
  call resize(kz,nele_new,0,'kz')
  call resize(kp,nele_new,0,'kp')
  call resize(bbbx,nele_new,zero,'bbbx')
  call resize(bbby,nele_new,zero,'bbby')
  call resize(bbbs,nele_new,zero,'bbbs')
  call resize(xpl,nele_new,zero,'xpl')
  call resize(zpl,nele_new,zero,'zpl')
  call resize(xrms,nele_new,zero,'xrms')
  call resize(zrms,nele_new,zero,'zrms')
  call resize(a,nele_new,2,6,zero,'a')
  call resize(hsyc,nele_new,zero,'hsyc')
  call resize(phasc,nele_new,zero,'phasc')
  call resize(itionc,nele_new,0,'itionc')
  call resize(bk0,nele_new,mmul,zero,'bk0')
  call resize(ak0,nele_new,mmul,zero,'ak0')
  call resize(bka,nele_new,mmul,zero,'bka')
  call resize(aka,nele_new,mmul,zero,'aka')
  call resize(benkc,nele_new,zero,'benkc')
  call resize(r00,nele_new,zero,'r00')
  call resize(irm,nele_new,0,'irm')
  call resize(nmu,nele_new,0,'nmu')
  
  call resize(bezr, max_name_len, 3, nele_new, str_nmZeros, 'bezr')
  call resize(kpa,nele_new,0,'kpa')
  
  call resize(bez,  max_name_len, nele_new, str_nmZeros, 'bez')
! call resize(bezb, max_name_len, nele_new, ' ', 'bezb')
  call resize(bezl, max_name_len, nele_new, str_nmZeros, 'bezl')
  
  call resize(ncororb,nele_new,0,'ncororb')
  call resize(ratioe,nele_new,zero,'ratioe')
  call resize(iratioe,nele_new,0,'iratioe')
  call resize(isea,nele_new,0,'isea')
  call resize(dki,nele_new,3,zero,'dki')
  call resize(parbe,nele_new,18,zero,'parbe')
  call resize(ptnfac,nele_new,zero,'ptnfac')
  call resize(imtr,nele_new,0,'imtr')
  call resize(acdipph,nele_new,zero,'acdipph')
  call resize(nturn1,nele_new,0,'nturn1')
  call resize(nturn2,nele_new,0,'nturn2')
  call resize(nturn3,nele_new,0,'nturn3')
  call resize(nturn4,nele_new,0,'nturn4')
  call resize(crabph,nele_new,zero,'crabph')
  call resize(crabph2,nele_new,zero,'crabph2')
  call resize(crabph3,nele_new,zero,'crabph3')
  call resize(crabph4,nele_new,zero,'crabph4')
  
  ! linopt
  call alloc(bezb,max_name_len,nblo_new,str_nmZeros,"bezb")
  call alloc(elbe,nblo_new,zero,"elbe")
  
  ! str2
  call alloc(mel,nblo_new,0,"mel")
  call alloc(mtyp,nblo_new,nelb,0,"mtyp")
  call alloc(mstr,nblo_new,0,"mstr")
  
  ! mat
  call alloc(bl1,nblo_new,2,6,zero,"bl1")
  call alloc(bl2,nblo_new,2,6,zero,"bl2")
  
  ! Note: One should probably here check the allocation status of the arrays, i.e. that they all actually are allocated...
    
end subroutine mod_common_expand_arrays

end module mod_common

module mod_commond
  
  use parpro
  use floatPrecision
  
  implicit none
  
  ! common /dial/
  real(kind=fPrec),            save :: preda
  integer,                     save :: idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)
  
  ! common /norf/
  integer,                     save :: nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2
  
  ! common /tcorr/
  real(kind=fPrec),            save :: weig1,weig2,dpmax
  integer,                     save :: icorr,nctype,namp,nmom,nmom1,nmom2
  character(len=max_name_len), save :: coel(10)
  
end module mod_commond
  
module mod_commont
  
  use parpro
  use floatPrecision
  
  implicit none
  
  ! common /tra/
  real(kind=fPrec), save :: x(mpa,2),y(mpa,2),amp(2),bet0(2),alf0(2),clo(2),clop(2)
  ! common /chrom/
  real(kind=fPrec), save :: cro(2)
  integer,          save :: is(2),ichrom
  
  ! common /tasm/
  real(kind=fPrec), save :: tasm(6,6)
  
  ! common /track/
  integer,          save :: ktrack(nblz),nwri
  real(kind=fPrec), save :: strack(nblz),strackc(nblz),stracks(nblz),strackx(nblz),strackz(nblz),dpsv1(npart)
  
  ! Substitute variables for x,y and is for DA version
  real(kind=fPrec), save :: xxtr(mpa,2),yytr(mpa,2)
  integer,          save :: issss(2)
  
contains

! Copy from actual variables to temp DA variables
subroutine comt_daStart
  xxtr(1:mpa,1:2) = x(1:mpa,1:2)
  yytr(1:mpa,1:2) = y(1:mpa,1:2)
  issss(1:2)      = is(1:2)
end subroutine comt_daStart

! Copy from temp DA variables to actual variables
subroutine comt_daEnd
  x(1:mpa,1:2) = xxtr(1:mpa,1:2)
  y(1:mpa,1:2) = yytr(1:mpa,1:2)
  is(1:2)      = issss(1:2)
end subroutine comt_daEnd

end module mod_commont

module mod_commonmn
  use parpro
  use floatPrecision
  implicit none
  
  ! common /main1/
  real(kind=fPrec), allocatable, save :: ekv(:,:) !(npart,nele)
  real(kind=fPrec), save :: fokqv(npart),aaiv(mmul,nmac,nblz),           &
       bbiv(mmul,nmac,nblz),smiv(nmac,nblz),zsiv(nmac,nblz),             &
       xsiv(nmac,nblz),xsv(npart),zsv(npart),qw(2),qwc(3),clo0(2),       &
       clop0(2),eps(2),epsa(2),ekk(2),cr(mmul),ci(mmul),xv(2,npart),     &
       yv(2,npart),dam(npart),ekkv(npart),sigmv(npart),dpsv(npart),      &
       dp0v(npart),sigmv6(npart),dpsv6(npart),ejv(npart),ejfv(npart),    &
       xlv(npart),zlv(npart),rvv(npart),pttemp,temptr(6),                &
! +if rvet
!        rvet(npart),                                                      &
! +ei
       ejf0v(npart)
  integer, save :: numxv(npart),nms(npart),nlostp(npart),kxxa
  logical, save :: pstop(npart)
  
  ! common /main2/
  real(kind=fPrec), save :: dpd(npart),dpsq(npart),fok(npart),rho(npart),&
       fok1(npart),si(npart),co(npart),g(npart),gl(npart),sm1(npart),    &
       sm2(npart),sm3(npart),sm12(npart),as3(npart),as4(npart),          &
       as6(npart),sm23(npart),rhoc(npart),siq(npart),aek(npart),         &
       afok(npart),hp(npart),hm(npart),hc(npart),hs(npart),wf(npart),    &
       wfa(npart),wfhi(npart),rhoi(npart),hi(npart),fi(npart),hi1(npart),&
       xvl(2,npart),yvl(2,npart),ejvl(npart),dpsvl(npart),oidpsv(npart), &
       sigmvl(npart),aperv(npart,2),clov(2,npart),                       &
       clopv(2,npart),alf0v(npart,2),bet0v(npart,2),ampv(npart)
  integer, save :: iv(npart),ixv(npart)
  
  !common /main3/
  real(kind=fPrec), save :: clo6v(3,npart),clop6v(3,npart),tas(npart,6,6),qwcs(npart,3),di0xs(npart)
  real(kind=fPrec), save :: di0zs(npart),dip0xs(npart),dip0zs(npart),xau(2,6),cloau(6)
  real(kind=fPrec), save :: di0au(4),tau(6,6),tasau(npart,6,6),wx(3),x1(6),x2(6),fake(2,20)
  real(kind=fPrec), allocatable, save :: hv(:,:,:,:),bl1v(:,:,:,:) !(6,2,npart,nblo)
  
  ! common /main4/
  integer,          save :: numx
  real(kind=fPrec), save :: e0f
  logical,          save :: sythckcr ! Only used for CR
  
contains
  subroutine mod_commonmn_allocate_arrays
  use mod_alloc
    implicit none
  end subroutine mod_commonmn_allocate_arrays

  subroutine mod_commonmn_expand_arrays(nele_new)
  use mod_alloc
    implicit none
    integer, intent(in) :: nele_new
    integer stat

  end subroutine mod_commonmn_expand_arrays

  subroutine mod_commonmn_allocate_thickarrays
    use numerical_constants, only : zero
    use crcoall
    use mod_alloc
    implicit none
    
!    integer i1,i2,i3,i
!    integer stat
    
    write(lout,*) "MOD_COMMONMN_ALLOCATE_THICKARRAYS: npart/nele/nblo=", npart, nele, nblo
    
!    allocate(ekv(npart,nele), hv(6,2,npart,nblo), bl1v(6,2,npart,nblo), STAT = stat)
!    if (stat.ne.0) then
!       write(lout,*) "ERROR in mod_commonmn_allocate_thickarrays(); stat=",stat
!       call prror(-1)
!    endif
   
    call alloc(ekv,npart,nele,zero,'ekv')
    call alloc(hv,6,2,npart,nblo,zero,'hv')
    call alloc(bl1v,6,2,npart,nblo,zero,'bl1v')

   !ZERO the newly allocated arrays
!    do i=1,npart
!       do i1=1,nblo
!          do i2=1,2
!             do i3=1,6
!                hv(i3,i2,i,i1)=zero
!                bl1v(i3,i2,i,i1)=zero
!             end do
!          end do
!       end do
!    end do


  end subroutine mod_commonmn_allocate_thickarrays

  subroutine mod_commonmn_expand_thickarrays(nele_new, npart_new, nblo_new)
    use crcoall
    use numerical_constants, only : zero
    use mod_alloc
    implicit none
 
    integer,intent(in) :: nele_new, npart_new, nblo_new

    call resize(ekv,npart_new,nele_new,zero,'ekv')
    call resize(hv,6,2,npart_new,nblo_new,zero,'hv')
    call resize(bl1v,6,2,npart_new,nblo_new,zero,'bl1v')

  end subroutine mod_commonmn_expand_thickarrays

end module mod_commonmn

module mod_commons
  use floatPrecision
  use parpro
  use mod_alloc

  implicit none
  
  ! common /syos/
  real(kind=fPrec), allocatable, save :: as(:,:,:,:),al(:,:,:,:) !(6,2,npart,nele)
  real(kind=fPrec), allocatable, save :: at(:,:,:,:),a2(:,:,:,:) !(6,2,npart,nele)
  real(kind=fPrec), save :: sigm(mpa),dps(mpa)
  integer, save :: idz(2)
  
  ! common /anf/
  real(kind=fPrec), save :: chi0,chid,exz(2,6),dp1
  integer, save :: itra

contains

  subroutine mod_commons_allocate_arrays
    implicit none
  end subroutine mod_commons_allocate_arrays

  subroutine mod_commons_expand_arrays(nele_new)
    implicit none
    integer, intent(in) :: nele_new
    integer stat

  end subroutine mod_commons_expand_arrays
  
  subroutine mod_commons_allocate_thickarrays
    !Subroutine for initial allocation of the big thick tracking arrays 'al' and 'as'
    use numerical_constants, only : zero
    use crcoall
    implicit none
    
!    integer stat
!    integer i1,i3,i4,i
    
    write(lout,*) "MOD_COMMONS_ALLOCATE_THICKARRAYS: npart/nele=", npart, nele
    
!    allocate(al(6,2,npart,nele), as(6,2,npart,nele), STAT = stat)
!    if (stat.ne.0) then
!       write(lout,*) "ERROR in mod_commons_allocate_thickarrays(); stat=",stat
!       call prror(-1)
!    endif
    
!    !ZERO the newly allocated arrays
!    do i=1,nele
!       do i3=1,2
!          do i4=1,6
!             do i1=1,npart
!                al(i4,i3,i1,i)=zero
!                as(i4,i3,i1,i)=zero
!             end do
!          end do
!       end do
!    end do

     call alloc(al,6,2,npart,nele,zero,'al')
     call alloc(as,6,2,npart,nele,zero,'as')
     call alloc(at,6,2,2,nele,zero,'at')
     call alloc(a2,6,2,2,nele,zero,'a2')

  end subroutine mod_commons_allocate_thickarrays

  subroutine mod_commons_expand_thickarrays(nele_new, npart_new)
    use numerical_constants, only : zero
    use crcoall
    implicit none
 
    integer,intent(in) :: nele_new, npart_new

     call resize(al,6,2,npart_new,nele_new,zero,'al')
     call resize(as,6,2,npart_new,nele_new,zero,'as')

  end subroutine mod_commons_expand_thickarrays
end module mod_commons

module mod_commond2
  use floatPrecision
  use parpro, only : nele, nema
  use crcoall
  use numerical_constants
  use mod_alloc

  implicit none

  !Note: These only to be neccessary for thick4d?
  real(kind=fPrec), allocatable, save :: ald6(:,:,:,:), asd6(:,:,:,:) !(nele,2,6,nema)

contains
  subroutine mod_commond2_allocate_arrays
    implicit none
    integer stat

!    allocate(ald6(nele,2,6,nema), asd6(nele,2,6,nema), STAT=stat)
!    if (stat.ne.0) then
!       write(lout,'(A,I8)') "ERROR in SUBROUTINE mod_commond2_allocate_arrays; stat=", stat
!       call prror(-1)
!    endif
    
    call alloc(ald6,nele,2,6,nema,zero,'ald6')
    call alloc(asd6,nele,2,6,nema,zero,'asd6')

  end subroutine mod_commond2_allocate_arrays

  subroutine mod_commond2_expand_arrays(nele_new)
    implicit none
    integer, intent(in) :: nele_new
    integer stat

    call resize(ald6,nele_new,2,6,nema,zero,'ald6')
    call resize(asd6,nele_new,2,6,nema,zero,'asd6')
    
  end subroutine mod_commond2_expand_arrays
  
end module mod_commond2
