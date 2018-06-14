!
!  This file contains all the main modules holding the SixTrack shared variables
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ================================================================================================ !
!  MAIN ARRAYS SIZE COMMON VARIABLES
!  Last modified: 2018-06-13
! ================================================================================================ !
module parpro

  implicit none

! All of the following need a description
  integer :: mbea
  integer :: mcor
  integer :: mcop
  integer :: mmul
  integer :: mpa
  integer :: mran
  integer :: nbb   ! Maximum number of beam-beam elements
  integer :: ncom
  integer :: ncor1
  integer :: nelb
  integer :: nema
  integer :: ninv
  integer :: nlya
  integer :: nmac
  integer :: nmon1
  integer :: npart ! Maximum number of particles
  integer :: nper
  integer :: nplo
  integer :: npos
  integer :: nran
  integer :: nrco
  integer :: ntr
  integer :: nzfz

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
  integer, parameter :: nblo_initial = 50   ! Must be at least 1
  integer, parameter :: nblz_initial = 15000 ! Must be at least 1
  integer :: nele = -1
  integer :: nblo = -1
  integer :: nblz = -1

!Note: nzfz should be = 3*nblz+2*mmul*#MULTIPOLES,
! where #MULTIPOLES are the max number of multipoles in the lattice (up to nblz)
! For now, scale the number of multipoles (from nzfz) as is done in the "no-flag" version:
! 6000/20000 -> 30% multipoles
#ifndef COLLIMAT
#ifdef BIGNBLZ
  parameter(nper=16,nelb=140,nzfz = 3000000,mmul = 20) !up to 60'000 multipoles
#endif
#ifdef HUGENBLZ
  parameter(nper=16,nelb=280,nzfz = 6000000,mmul = 20) !up to 120'000 multipoles -> 48MB/nzfz-array (20%)
#endif
#if !defined(BIGNBLZ) && !defined(HUGENBLZ)
  parameter(nper=16,nelb=140,nzfz = 300000,mmul = 20) !up to 6'000 multipoles
#endif
#else
#ifdef BEAMGAS
  parameter(nper=16,nelb=140,nzfz = 1920000,mmul = 11) !up to 60'000 multipoles
#else
#ifdef BIGNBLZ
  parameter(nper=16,nelb=140,nzfz = 1920000,mmul = 11) !up to 60'000 multipoles
#endif
#ifdef HUGENBLZ
  parameter(nper=16,nelb=140,nzfz = 3840000,mmul = 11) !up to 120'000 multipoles (20%)
#endif
#if !defined(BIGNBLZ) && !defined(HUGENBLZ)
  parameter(nper=16,nelb=140,nzfz = 144000,mmul = 11) !up to 4500 multipoles
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

! ================================================================================================ !
!  PARBEAM COMMON VARIABLES
!  Last modified: 2018-06-12
! ================================================================================================ !
module parbeam

  use floatPrecision

  implicit none

  real(kind=fPrec), parameter :: xcut = 7.77_fPrec
  real(kind=fPrec), parameter :: ycut = 7.46_fPrec
  real(kind=fPrec), parameter :: h    =  1.0_fPrec/63.0_fPrec
  integer,          parameter :: nx   = 490
  integer,          parameter :: ny   = 470
  integer,          parameter :: idim = (nx+2)*(ny+2)

  ! common /wzcom1/
  real(kind=fPrec), save :: hrecip
  integer,          save :: kstep

  ! common /wzcom2/
  real(kind=fPrec), save :: wtreal(idim),wtimag(idim)

  ! common /beam_exp/
  integer,          save :: beam_expflag      ! 0: Old BEAM block, 1: New BEAM::EXPERT
  logical,          save :: beam_expfile_open ! Have we opened the file 'beam_expert.txt'?

end module parbeam

! ================================================================================================ !
!  Global Settings Module
!  Last modiffied: 2018-06-10
!  Holds global settings values and parameters not directly related to the physics
! ================================================================================================ !
module mod_settings

  implicit none

  ! PRINT Flag (fort.3)
  logical, save :: st_print

  ! SETTINGS Block (fort.3)
  integer, save :: st_quiet ! QUIET Level
  logical, save :: st_debug ! Global DEBUG flag

  ! String Stuff
  integer,            parameter :: st_divLen  = 132
  character(len=132), parameter :: st_divLine = repeat("-",132)

end module mod_settings

! ================================================================================================ !
!  THAT BIG COMMON VARIABLES MODULE
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_common

  use parpro
  use floatPrecision
  use numerical_constants, only : c1m6

  implicit none

  ! common /erro/
  integer, save :: ierro
  integer, save :: errout_status
  character(len=max_name_len), save :: erbez

  ! common /kons/
  real(kind=fPrec), save :: pi,pi2,pisqrt,rad

  ! common /str/
  integer,              save :: il,mper,mblo,mbloz,msym(nper),kanf,iu
  integer, allocatable, save :: ic(:) !(nblz)

  ! common /ell/
  real(kind=fPrec), allocatable, save :: ed(:),el(:),ek(:),sm(:)        ! (nele)
  integer,          allocatable, save :: kz(:),kp(:)                    ! (nele)

  ! common /bbb/
  real(kind=fPrec), allocatable, save :: bbbx(:),bbby(:),bbbs(:)        ! (nele)

  ! common /pla/
  real(kind=fPrec), allocatable, save :: xpl(:),xrms(:),zpl(:),zrms(:)  ! (nele)

  ! common /str2/
  integer,          allocatable, save :: mel(:)    ! (nblo)
  integer,          allocatable, save :: mtyp(:,:) ! (nblo,nelb)
  integer,          allocatable, save :: mstr(:)   ! (nblo)

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
  real(kind=fPrec), save :: qs,e0,pma,ej(mpa),ejf(mpa),phas0,phas,hsy(3),crad,dppoff,tlen,mtcda
  integer,          save :: iicav,ition,idp,ncy,ixcav

  real(kind=fPrec), allocatable, save :: sigmoff(:)       ! (nblz)
  real(kind=fPrec), allocatable, save :: hsyc(:),phasc(:) ! (nele)
  integer,          allocatable, save :: itionc(:)        ! (nele)

  ! common /corcom/
  real(kind=fPrec), save :: dpscor,sigcor
  integer,          save :: icode,idam,its6d

  ! Multipole Coefficients
  real(kind=fPrec),              save :: benki
  real(kind=fPrec), allocatable, save :: benkc(:),r00(:)                     ! (nele)
  real(kind=fPrec), allocatable, save :: bk0(:,:),ak0(:,:),bka(:,:),aka(:,:) ! (nele,mmul)
  integer,          allocatable, save :: irm(:),nmu(:)                       ! (nele)

  ! common /rand0/
  real(kind=fPrec), save :: zfz(nzfz)
  integer,          save :: iorg,izu0,mmac,mcut

  character(len=:), allocatable, save :: bezr(:,:) ! (max_name_len)(3,nele)
  integer,          allocatable, save :: mzu(:)    ! (nblz)

  ! common /rand1/
  integer,                       save :: mout2
  integer,          allocatable, save :: icext(:),icextal(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: exterr(:,:)         ! (nblz,40)
  real(kind=fPrec), allocatable, save :: extalign(:,:)       ! (nblz,3)
  real(kind=fPrec), allocatable, save :: tiltc(:),tilts(:)   ! (nblz)

  ! common /beo/
  real(kind=fPrec), save :: aper(2),di0(2),dip0(2),ta(6,6)

  ! common /clo/
  real(kind=fPrec), save :: dma,dmap,dkq,dqq,de0,ded,dsi,dech,dsm0
  integer, save :: itco,itcro,itqv

  ! common /qmodi/
  real(kind=fPrec),              save :: qw0(3),amp0
  integer,                       save :: iq(3),iqmod,iqmod6
  integer,          allocatable, save :: kpa(:) !(nele)

  ! common /linop/
  character(len=:), allocatable, save :: bez(:), bezl(:) ! (nele)
  character(len=:), allocatable, save :: bezb(:)         ! (nblo)
  real(kind=fPrec), allocatable, save :: elbe(:)         ! (nblo)
  real(kind=fPrec),              save :: eui,euii
  integer,                       save :: ilin,nt,iprint,ntco,nlin

  ! common /cororb/
  real(kind=fPrec),              save :: betam(nmon1,2),pam(nmon1,2),betac(ncor1,2),pac(ncor1,2),bclorb(nmon1,2)
  integer,                       save :: nhmoni,nhcorr,nvmoni,nvcorr
  integer,          allocatable, save :: ncororb(:) ! nele

  ! common /clos/
  real(kind=fPrec), save :: sigma0(2)
  integer,          save :: iclo, ncorru,ncorrep

  ! common /combin/
  real(kind=fPrec),              save :: ratio(ncom,20)
  real(kind=fPrec), allocatable, save :: ratioe(:) ! (nele)
  integer,                       save :: icomb0(20),icomb(ncom,20),icoe
  integer,          allocatable, save :: iratioe(:) ! (nele)

  ! common/seacom/m21,m22,m23
  integer,              save :: ise,ise1,ise2,ise3,mesa,mp,m21,m22,m23
  integer, allocatable, save :: isea(:) ! (nele)

  ! common /subres/
  real(kind=fPrec), save :: qxt,qzt,tam1,tam2,totl
  integer,          save :: isub,nta,nte,ipt

  ! common /secom/
  real(kind=fPrec), save :: rtc(9,18,10,5),rts(9,18,10,5)
  integer,          save :: ire(12),ipr(5),irmod2

  ! common /secom1/
  real(kind=fPrec), save :: dtr(10)
  integer,          save :: nre,nur,nch,nqc,npp,nrr(5),nu(5)

  ! common /postr/
  real(kind=fPrec), save :: dphix,dphiz,qx0,qz0,dres,dfft,cma1,cma2
  integer,          save :: nstart,nstop,iskip,iconv,imad

  ! common /posti1/
  integer,           save :: ipos,iav,iwg,ivox,ivoz,ires,ifh
  character(len=80), save :: toptit(5) !DANGER: If the len changes, CRCHECK will break.

  ! common /posti2/
  integer, save :: kwtype,itf,icr,idis,icow,istw,iffw,nprint,ndafi

  ! common /skew/
  real(kind=fPrec), save :: qwsk(2),betx(2),betz(2),alfx(2),alfz(2)
  integer,          save :: iskew,nskew(6)

  ! common /pawc/
  real, save :: hmal(nplo)

  ! common /tit/
  character (len=80), save :: sixtit,commen !DANGER: If the len changes, CRCHECK will break.
  integer,            save :: ithick

  ! common/co6d/
  real(kind=fPrec), save :: clo6(3),clop6(3)

  ! common /dkic/
  real(kind=fPrec), allocatable, save :: dki(:,:) !(nele,3)

  ! common /beam/
  real(kind=fPrec),              save :: sigman(2,nbb),sigman2(2,nbb),sigmanq(2,nbb)
  real(kind=fPrec),              save :: clobeam(6,nbb),beamoff(6,nbb),track6d(6,npart)
  real(kind=fPrec),              save :: sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar
  integer,                       save :: nbeam,ibbc,ibeco,ibtyp,lhc
  real(kind=fPrec), allocatable, save :: parbe(:,:) ! (nele,18)
  real(kind=fPrec), allocatable, save :: ptnfac(:)  ! (nele)

  ! common/trom/
  real(kind=fPrec),              save :: cotr(ntr,6),rrtr(ntr,6,6)
  integer,          allocatable, save :: imtr(:) ! (nele),

  ! common /bb6d/
  real(kind=fPrec), save :: bbcu(nbb,12)
  integer,          save :: ibb6d
  integer,          allocatable, save :: imbb(:) ! (nblz)

  ! common /acdipco/
  real(kind=fPrec), allocatable, save :: acdipph(:) ! (nele)
  integer,          allocatable, save :: nturn1(:),nturn2(:),nturn3(:),nturn4(:) ! (nele)

  ! common /crabco/
  real(kind=fPrec), allocatable, save :: crabph(:),crabph2(:),crabph3(:),crabph4(:) ! (nele)

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

  ! common /mytimes/
  real,             save :: timestart

  ! common /xz/
  real(kind=fPrec), allocatable, save :: xsi(:),zsi(:)     ! (nblz)
  real(kind=fPrec), allocatable, save :: smi(:),smizf(:)   ! (nblz)
  real(kind=fPrec), allocatable, save :: aai(:,:),bbi(:,:) ! (nblz,mmul)

  ! common /dcumdb/
  real(kind=fPrec), allocatable, save :: dcum(:)              ! (0:nblz+1) Machine length in m
  real(kind=fPrec), parameter         :: eps_dcum   = c1m6    ! Tolerance for machine length mismatch [m]
  logical,          parameter         :: print_dcum = .false.

  ! beamdim
  real(kind=fPrec), parameter         :: cc   = 1.12837916709551_fPrec
  real(kind=fPrec), parameter         :: xlim = 5.33_fPrec
  real(kind=fPrec), parameter         :: ylim = 4.29_fPrec

contains

subroutine mod_common_expand_arrays(nele_new, nblo_new, nblz_new)

  use mod_alloc
  use string_tools
  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: nblo_new
  integer, intent(in) :: nblz_new

  call alloc(ed,                       nele_new,       zero,        "ed")
  call alloc(el,                       nele_new,       zero,        "el")
  call alloc(ek,                       nele_new,       zero,        "ek")
  call alloc(sm,                       nele_new,       zero,        "sm")
  call alloc(kz,                       nele_new,       0,           "kz")
  call alloc(kp,                       nele_new,       0,           "kp")
  call alloc(bbbx,                     nele_new,       zero,        "bbbx")
  call alloc(bbby,                     nele_new,       zero,        "bbby")
  call alloc(bbbs,                     nele_new,       zero,        "bbbs")
  call alloc(xpl,                      nele_new,       zero,        "xpl")
  call alloc(zpl,                      nele_new,       zero,        "zpl")
  call alloc(xrms,                     nele_new,       zero,        "xrms")
  call alloc(zrms,                     nele_new,       zero,        "zrms")
  call alloc(a,                        nele_new,2,6,   zero,        "a")
  call alloc(hsyc,                     nele_new,       zero,        "hsyc")
  call alloc(phasc,                    nele_new,       zero,        "phasc")
  call alloc(itionc,                   nele_new,       0,           "itionc")
  call alloc(bk0,                      nele_new, mmul, zero,        "bk0")
  call alloc(ak0,                      nele_new, mmul, zero,        "ak0")
  call alloc(bka,                      nele_new, mmul, zero,        "bka")
  call alloc(aka,                      nele_new, mmul, zero,        "aka")
  call alloc(benkc,                    nele_new,       zero,        "benkc")
  call alloc(r00,                      nele_new,       zero,        "r00")
  call alloc(irm,                      nele_new,       0,           "irm")
  call alloc(nmu,                      nele_new,       0,           "nmu")
  call alloc(bezr,    max_name_len, 3, nele_new,       str_nmZeros, "bezr")
  call alloc(kpa,                      nele_new,       0,           "kpa")
  call alloc(bez,     max_name_len,    nele_new,       str_nmZeros, "bez")
! call alloc(bezb,    max_name_len,    nele_new,       str_nmZeros, "bezb")
  call alloc(bezl,    max_name_len,    nele_new,       str_nmZeros, "bezl")
  call alloc(ncororb,                  nele_new,       0,           "ncororb")
  call alloc(ratioe,                   nele_new,       zero,        "ratioe")
  call alloc(iratioe,                  nele_new,       0,           "iratioe")
  call alloc(isea,                     nele_new,       0,           "isea")
  call alloc(dki,                      nele_new, 3,    zero,        "dki")
  call alloc(parbe,                    nele_new, 18,   zero,        "parbe")
  call alloc(ptnfac,                   nele_new,       zero,        "ptnfac")
  call alloc(imtr,                     nele_new,       0,           "imtr")
  call alloc(acdipph,                  nele_new,       zero,        "acdipph")
  call alloc(nturn1,                   nele_new,       0,           "nturn1")
  call alloc(nturn2,                   nele_new,       0,           "nturn2")
  call alloc(nturn3,                   nele_new,       0,           "nturn3")
  call alloc(nturn4,                   nele_new,       0,           "nturn4")
  call alloc(crabph,                   nele_new,       zero,        "crabph")
  call alloc(crabph2,                  nele_new,       zero,        "crabph2")
  call alloc(crabph3,                  nele_new,       zero,        "crabph3")
  call alloc(crabph4,                  nele_new,       zero,        "crabph4")

  call alloc(bezb,    max_name_len,    nblo_new,       str_nmZeros, "bezb")
  call alloc(elbe,                     nblo_new,       zero,        "elbe")
  call alloc(mel,                      nblo_new,       0,           "mel")
  call alloc(mtyp,                     nblo_new, nelb, 0,           "mtyp")
  call alloc(mstr,                     nblo_new,       0,           "mstr")
  call alloc(bl1,                      nblo_new, 2, 6, zero,        "bl1")
  call alloc(bl2,                      nblo_new, 2, 6, zero,        "bl2")

  call alloc(ic,                       nblz_new,       0,           "ic")
  call alloc(mzu,                      nblz_new,       0,           "mzu")
  call alloc(imbb,                     nblz_new,       0,           "imbb")
  call alloc(icext,                    nblz_new,       0,           "icext")
  call alloc(icextal,                  nblz_new,       0,           "icextal")
  call alloc(exterr,                   nblz_new, 40,   zero,        "exterr")
  call alloc(extalign,                 nblz_new, 3,    zero,        "extalign")
  call alloc(tiltc,                    nblz_new,       zero,        "tiltc")
  call alloc(tilts,                    nblz_new,       zero,        "tilts")
  call alloc(xsi,                      nblz_new,       zero,        "xsi")
  call alloc(zsi,                      nblz_new,       zero,        "zsi")
  call alloc(smi,                      nblz_new,       zero,        "smi")
  call alloc(smizf,                    nblz_new,       zero,        "smizf")
  call alloc(aai,                      nblz_new, mmul, zero,        "aai")
  call alloc(bbi,                      nblz_new, mmul, zero,        "bbi")
  call alloc(dcum,                     nblz_new+1,     zero,        "dcum", 0)
  call alloc(sigmoff,                  nblz_new,       zero,        "sigmoff")

end subroutine mod_common_expand_arrays

end module mod_common

! ================================================================================================ !
!  DA COMMON VARIABLES 1
!  Last modified: 2018-06-12
! ================================================================================================ !
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

! ================================================================================================ !
!  DA COMMON VARIABLES 2
!  Last modified: 2018-06-12
!  These variables are used by FOX and are also declared internally in the FOX code in dainicom.f90
!  The variables are defined here, but there is no direct way to avoid using common/daele/. FixMe
! ================================================================================================ !
module mod_common_da2

  use parpro, only : mcor

  implicit none

  ! common /daele/
  integer, public :: alda(2,6)
  integer, public :: asda(2,6)
  integer, public :: aldaq(2,6)
  integer, public :: asdaq(2,6)
  integer, public :: smida(mcor)
  integer, public :: xx(2)
  integer, public :: yy(2)
  integer, public :: dpda
  integer, public :: dpda1
  integer, public :: sigmda
  integer, public :: ej1
  integer, public :: ejf1
  integer, public :: rv

  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv

end module mod_common_da2

! ================================================================================================ !
!  TRACKING COMMON VARIABLES
!  Last modified: 2018-06-13
! ================================================================================================ !
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
  integer,          save :: nwri
  real(kind=fPrec), save :: dpsv1(npart)

  integer,          allocatable, save :: ktrack(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: strack(:),strackc(:),stracks(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: strackx(:),strackz(:) ! (nblz)

  ! Substitute variables for x,y and is for DA version
  real(kind=fPrec), save :: xxtr(mpa,2),yytr(mpa,2)
  integer,          save :: issss(2)

contains

subroutine mod_commont_expand_arrays(nblz_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: nblz_new

  call alloc(ktrack,  nblz_new, 0,    "ktrack")
  call alloc(strack,  nblz_new, zero, "strack")
  call alloc(strackc, nblz_new, zero, "strackc")
  call alloc(stracks, nblz_new, zero, "stracks")
  call alloc(strackx, nblz_new, zero, "strackx")
  call alloc(strackz, nblz_new, zero, "strackz")

end subroutine mod_commont_expand_arrays

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

! ================================================================================================ !
!  MAIN COMMON VARIABLES
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_commonmn

  use parpro
  use floatPrecision

  implicit none

  ! common /main1/
  real(kind=fPrec), allocatable, save :: ekv(:,:)    ! (npart,nele)
  real(kind=fPrec), allocatable, save :: aaiv(:,:,:) ! (mmul,nmac,nblz)
  real(kind=fPrec), allocatable, save :: bbiv(:,:,:) ! (mmul,nmac,nblz)
  real(kind=fPrec), allocatable, save :: smiv(:,:)   ! (nmac,nblz)
  real(kind=fPrec), allocatable, save :: zsiv(:,:)   ! (nmac,nblz)
  real(kind=fPrec), allocatable, save :: xsiv(:,:)   ! (nmac,nblz)

  real(kind=fPrec), save :: fokqv(npart),xsv(npart),zsv(npart),qw(2),   &
    qwc(3),clo0(2),clop0(2),eps(2),epsa(2),ekk(2),cr(mmul),ci(mmul),    &
    xv(2,npart),yv(2,npart),dam(npart),ekkv(npart),sigmv(npart),        &
    dpsv(npart),dp0v(npart),sigmv6(npart),dpsv6(npart),ejv(npart),      &
    ejfv(npart),xlv(npart),zlv(npart),rvv(npart),pttemp,temptr(6),      &
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

subroutine mod_commonmn_expand_arrays(nblz_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: nblz_new

  call alloc(aaiv, mmul, nmac, nblz_new, zero, "aaiv")
  call alloc(bbiv, mmul, nmac, nblz_new, zero, "bbiv")
  call alloc(smiv,       nmac, nblz_new, zero, "smiv")
  call alloc(zsiv,       nmac, nblz_new, zero, "zsiv")
  call alloc(xsiv,       nmac, nblz_new, zero, "xsiv")

end subroutine mod_commonmn_expand_arrays

subroutine mod_commonmn_allocate_thickarrays

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  call alloc(ekv,npart,nele,zero,'ekv')
  call alloc(hv,6,2,npart,nblo,zero,'hv')
  call alloc(bl1v,6,2,npart,nblo,zero,'bl1v')

end subroutine mod_commonmn_allocate_thickarrays

subroutine mod_commonmn_expand_thickarrays(nele_new, npart_new, nblo_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer,intent(in) :: nele_new, npart_new, nblo_new

  call alloc(ekv,npart_new,nele_new,zero,'ekv')
  call alloc(hv,6,2,npart_new,nblo_new,zero,'hv')
  call alloc(bl1v,6,2,npart_new,nblo_new,zero,'bl1v')

end subroutine mod_commonmn_expand_thickarrays

end module mod_commonmn

! ================================================================================================ !
!  SOMETHING-SOMETHING COMMON VARIABLES
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_commons

  use floatPrecision
  use parpro

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

subroutine mod_commons_allocate_thickarrays

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  call alloc(al,6,2,npart,nele,zero,'al')
  call alloc(as,6,2,npart,nele,zero,'as')
  call alloc(at,6,2,2,nele,zero,'at')
  call alloc(a2,6,2,2,nele,zero,'a2')

end subroutine mod_commons_allocate_thickarrays

subroutine mod_commons_expand_thickarrays(nele_new, npart_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer,intent(in) :: nele_new, npart_new

    call alloc(al,6,2,npart_new,nele_new,zero,'al')
    call alloc(as,6,2,npart_new,nele_new,zero,'as')

end subroutine mod_commons_expand_thickarrays

end module mod_commons

! ================================================================================================ !
!  MAIN COMMON VARIABLES
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_commond2

  use floatPrecision
  use parpro, only : nele, nema
  use crcoall
  use numerical_constants

  implicit none

  !Note: These only to be neccessary for thick4d?
  real(kind=fPrec), allocatable, save :: ald6(:,:,:,:), asd6(:,:,:,:) !(nele,2,6,nema)

contains

subroutine mod_commond2_expand_arrays(nele_new)

  use mod_alloc

  implicit none

  integer, intent(in) :: nele_new

  call alloc(ald6,nele_new,2,6,nema,zero,'ald6')
  call alloc(asd6,nele_new,2,6,nema,zero,'asd6')

end subroutine mod_commond2_expand_arrays

end module mod_commond2

! ================================================================================================ !
!  LIELIB AND DABNEW COMMON VARIABLES
!  Last modified: 2018-06-11
! ================================================================================================ !
module mod_lie_dab

  use floatPrecision

  implicit none

  ! From lielib
  integer, parameter :: ndim  = 3
  integer, parameter :: ndim2 = 6
  integer, parameter :: ntt   = 40
  integer, parameter :: nreso = 20

  integer,          save :: nd,nd2,no,nv
  integer,          save :: ifilt,idpr,iref,itu
  integer,          save :: lienot,iflow,jtune
  integer,          save :: ndc,ndc2,ndpt,ndt
  integer,          save :: nplane(ndim),ista(ndim),idsta(ndim)
  integer,          save :: mx(ndim,nreso),nres
  real(kind=fPrec), save :: epsplane,xplane(ndim)
  real(kind=fPrec), save :: sta(ndim),dsta(ndim),angle(ndim),rad(ndim)
  real(kind=fPrec), save :: ps(ndim),rads(ndim)
  real(kind=fPrec), save :: xintex(0:20)

  ! From dabnew
#ifdef SIXDA
  integer, parameter :: lda = 10000
  integer, parameter :: lst = 20050000
  integer, parameter :: lea = 100000
  integer, parameter :: lia = 5000000
#else
  integer, parameter :: lda = 10000
  integer, parameter :: lst = 200000
  integer, parameter :: lea = 500
  integer, parameter :: lia = 10000
#endif
  integer, parameter :: lno = 120
  integer, parameter :: lnv = 40

  integer,           save :: ndat,nda,ndamaxi,lfi
  integer,           save :: nst,nomax,nvmax,nmmax,nocut
  integer,           save :: idano(lda),idanv(lda),idapo(lda)
  integer,           save :: idalm(lda),idall(lda)
  integer,           save :: i1(lst),i2(lst)
  integer,           save :: ie1(lea),ie2(lea),ieo(lea)
  integer,           save :: ia1(0:lia),ia2(0:lia),ifi(lea)
  logical,           save :: allvec(lda)
  character(len=10), save :: daname(lda)
  real(kind=fPrec),  save :: cc(lst),eps,epsmac,facint(0:lno)

  ! dascr variables
  integer,           save :: idao
  integer,           save :: iscrda(100)
  integer,           save :: iscrri(100)
  real(kind=fPrec),  save :: rscrri(100)

end module mod_lie_dab
