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

  integer, parameter :: mbea  = 99        ! Maximum number of beam-beam slices
  integer, parameter :: mcor  = 10        ! Maximum number of extra parameters (DIFF block)
  integer, parameter :: mcop  = mcor + 6  ! FOX/DA variable
  integer, parameter :: mpa   = 6         ! Maximum number of trajectories
  integer, parameter :: mran  = 500       ! Maximum number of inserted elements (subroutine ord)
  integer, parameter :: ncom  = 100       ! Maximum number of combinations of elements (COMB block)
  integer, parameter :: ncor1 = 600       ! Maximum number of corrections of closed orbit
  integer, parameter :: nema  = 15        ! Maximum order of the one turn map (DIFF block)
  integer, parameter :: ninv  = 1000      ! Number of invariances (postprocessing)
  integer, parameter :: nlya  = 10000     ! Something something postprocessing
  integer, parameter :: nmon1 = 600       ! Maximum number of monitors (closed orbit)
  integer, parameter :: nper  = 16        ! Maximum number of super periods (BLOC list, line 1)
  integer, parameter :: nplo  = 20000     ! Plotting
  integer, parameter :: npos  = 20000     ! Something something postprocessing
  integer, parameter :: nran  = 2000000   ! Maximum size for scaling nzfz
  integer, parameter :: nrco  = 5         ! Maximum order of compensation (RESO block)
  integer, parameter :: mmul  = 20        ! Maximum order of multipoles
  integer, parameter :: nbb   = 500       ! Beam-beam lenses
  integer, parameter :: nelb  = 280       ! Maximum elements per BLOC

  ! Maximum length of element names
  integer, parameter :: mNameLen = 48     ! Maximum length of element names. Keep in sync with MadX
  integer, parameter :: mStrLen  = 161    ! Standard string length
  integer, parameter :: mDivLen  = 132    ! Length of lout output lines
  integer, parameter :: mInputLn = 1600   ! Buffer size for single lines read from input files

  integer :: ntr   = -1   ! Number of phase trombones

  integer :: nzfz  = -1   ! Number of allocated multipole random numbers
  integer :: nele  = -1   ! Number of allocated SINGle elements
  integer :: nblo  = -1   ! Number of allocated BLOCs
  integer :: nblz  = -1   ! Number of allocated STRUcture elements
  integer :: npart = -1   ! Number of allocated particles

  integer, parameter :: nele_initial  = 500
  integer, parameter :: nblo_initial  = 100
  integer, parameter :: nblz_initial  = 1000
  integer, parameter :: npart_initial = 2

  ! Dummy Strings
  character(len=mDivLen),  parameter :: str_divLine = repeat("-",mDivLen)
  character(len=mStrLen),  parameter :: str_dSpace  = repeat(" ",mStrLen)
  character(len=mNameLen), parameter :: str_nmSpace = repeat(" ",mNameLen)

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

  ! SETTINGS Block (fort.3)
  logical, save :: st_print      = .false. ! PRINT flag (fort.3)
  integer, save :: st_quiet      = 0       ! QUIET Level 0=verbose, 1=minimal, 2=quiet
  logical, save :: st_debug      = .false. ! Global DEBUG flag
  logical, save :: st_partsum    = .false. ! Flag to print final particle summary
  integer, save :: st_finalstate = 0       ! Dump particle final state (mod_particles)

  ! Checpoint/Restart Kills Switch Settings
  logical,              save :: st_killswitch = .false. ! Enables the kill on turn number debug feature
  integer, allocatable, save :: st_killturns(:)         ! List of killswitch turns

end module mod_settings

! ================================================================================================ !
!  THAT BIG COMMON VARIABLES MODULE
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_common

  use parpro
  use floatPrecision
  use numerical_constants, only : c1m6, c180e0, pi, two, half, one, zero

  implicit none

  ! common /erro/
  integer, save :: ierro
  integer, save :: errout_status
  character(len=mNameLen), save :: erbez

  ! common /kons/
  real(kind=fPrec),     save :: pi2    = half*pi
  real(kind=fPrec),     save :: twopi  = two*pi
  real(kind=fPrec),     save :: pisqrt = sqrt(pi)
  real(kind=fPrec),     save :: rad    = pi/c180e0

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
  integer,          save :: idfor
  integer,          save :: napx, napxo
  integer,          save :: numl, niu(2), numlr, nde(2), nwr(4)
  integer,          save :: ird, imc, irew, ntwin
  integer,          save :: iclo6, iclo6r, iver
  integer,          save :: numlcp, numlmax, nnuml

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
  real(kind=fPrec), allocatable, save :: benkc(:),r00(:),scalemu(:)         ! (nele)
  real(kind=fPrec), allocatable, save :: bk0(:,:),ak0(:,:),bka(:,:),aka(:,:) ! (nele,mmul)
  integer,          allocatable, save :: irm(:),nmu(:)                       ! (nele)

  ! RF multipoles
  real(kind=fPrec), allocatable, save :: norrfamp(:,:),norrfph(:,:),skrfamp(:,:),skrfph(:,:) ! (nele,mmul)
  integer,          allocatable, save :: nmu_rf(:), irm_rf(:)
  real(kind=fPrec), allocatable, save :: freq_rfm(:)

  ! common /rand0/
  real(kind=fPrec), allocatable, save :: zfz(:) ! (nzfz)
  integer,          save :: iorg,izu0,mcut

  character(len=:), allocatable, save :: bezr(:,:) ! (mNameLen)(3,nele)
  integer,          allocatable, save :: mzu(:)    ! (nblz)

  ! common /rand1/
  integer,                       save :: mout2
  integer,          allocatable, save :: icext(:),icextal(:) ! (nblz)
! real(kind=fPrec), allocatable, save :: exterr(:,:)         ! (nblz,40)
! real(kind=fPrec), allocatable, save :: extalign(:,:)       ! (nblz,3)
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
  real(kind=fPrec),              save :: clobeam(6,nbb),beamoff(6,nbb)
  real(kind=fPrec), allocatable, save :: track6d(:,:) ! (6,npart)
  real(kind=fPrec),              save :: sigz,sige,partnum,parbe14,emitx,emity,emitz
  real(kind=fPrec),              save :: gammar = one
  real(kind=fPrec),              save :: betrel = zero
  integer,                       save :: nbeam,ibbc,ibeco,ibtyp,lhc
  real(kind=fPrec), allocatable, save :: parbe(:,:) ! (nele,18)
  real(kind=fPrec), allocatable, save :: ptnfac(:)  ! (nele)

  ! common/trom/
  real(kind=fPrec), allocatable, save :: cotr(:,:)   ! (ntr,6)
  real(kind=fPrec), allocatable, save :: rrtr(:,:,:) ! (ntr,6,6)
  integer,          allocatable, save :: imtr(:)     ! (nele)

  ! common /bb6d/
  real(kind=fPrec), save :: bbcu(nbb,12)
  integer,          save :: ibb6d
  integer,          allocatable, save :: imbb(:) ! (nblz)

  ! common /acdipco/
  real(kind=fPrec), allocatable, save :: acdipph(:) ! (nele)
  integer,          allocatable, save :: nturn1(:),nturn2(:),nturn3(:),nturn4(:) ! (nele)

  ! common /crabco/
  real(kind=fPrec), allocatable, save :: crabph(:),crabph2(:),crabph3(:),crabph4(:) ! (nele)

  ! common /general-rf multi/
  integer, save :: iord, nordm
  real(kind=fPrec), save :: field_cos(2,mmul), fsddida(2,mmul)
  real(kind=fPrec), save :: field_sin(2,mmul), fcodda(2,mmul)
  ! common /exact/
  integer, save :: iexact
  integer, save :: curveff

  ! common /sixdim/
  real(kind=fPrec), save :: aml6(6,6),edcor(2)

  ! common /postr2/
  integer,          allocatable, save :: nnumxv(:) ! (npart)

  ! common /correct/
  integer,          save :: ichromc,ilinc,iqmodc
  real(kind=fPrec), save :: corr(3,3),chromc(2),wxys(3),clon(6)

  ! common /damp/
  real(kind=fPrec), save :: damp,ampt

  ! common /ttime/
  integer,          save :: napxto
  real,             save :: tlim,time0,time1,time2,time3,trtime,pretime,posttime,tottime

  ! common /xz/
  real(kind=fPrec), allocatable, save :: xsi(:),zsi(:)     ! (nblz)
  real(kind=fPrec), allocatable, save :: smi(:),smizf(:)   ! (nblz)
  real(kind=fPrec), allocatable, save :: aaiv(:,:),bbiv(:,:) ! (nblz,mmul)
  real(kind=fPrec), allocatable, save :: amultip(:,:), bmultip(:,:) ! (nblz,mmul)

  ! common /dcumdb/
  real(kind=fPrec), allocatable, save :: dcum(:)              ! (0:nblz+1) Machine length in m
  real(kind=fPrec), parameter         :: eps_dcum   = c1m6    ! Tolerance for machine length mismatch [m]
  logical,                       save :: print_dcum = .false. ! Set in the SETTINGS block

  ! beamdim
  real(kind=fPrec), parameter         :: cc   = 1.12837916709551_fPrec
  real(kind=fPrec), parameter         :: xlim = 5.33_fPrec
  real(kind=fPrec), parameter         :: ylim = 4.29_fPrec

contains

subroutine mod_common_expand_arrays(nele_new, nblo_new, nblz_new, npart_new)

  use mod_alloc
  use mod_settings
  use numerical_constants, only : zero,one

  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: nblo_new
  integer, intent(in) :: nblz_new
  integer, intent(in) :: npart_new

  call alloc(ed,                   nele_new,       zero,        "ed")
  call alloc(el,                   nele_new,       zero,        "el")
  call alloc(ek,                   nele_new,       zero,        "ek")
  call alloc(sm,                   nele_new,       zero,        "sm")
  call alloc(kz,                   nele_new,       0,           "kz")
  call alloc(kp,                   nele_new,       0,           "kp")
  call alloc(bbbx,                 nele_new,       zero,        "bbbx")
  call alloc(bbby,                 nele_new,       zero,        "bbby")
  call alloc(bbbs,                 nele_new,       zero,        "bbbs")
  call alloc(xpl,                  nele_new,       zero,        "xpl")
  call alloc(zpl,                  nele_new,       zero,        "zpl")
  call alloc(xrms,                 nele_new,       zero,        "xrms")
  call alloc(zrms,                 nele_new,       zero,        "zrms")
  call alloc(a,                    nele_new,2,6,   zero,        "a")
  call alloc(hsyc,                 nele_new,       zero,        "hsyc")
  call alloc(phasc,                nele_new,       zero,        "phasc")
  call alloc(itionc,               nele_new,       0,           "itionc")
  call alloc(bk0,                  nele_new, mmul, zero,        "bk0")
  call alloc(ak0,                  nele_new, mmul, zero,        "ak0")
  call alloc(bka,                  nele_new, mmul, zero,        "bka")
  call alloc(aka,                  nele_new, mmul, zero,        "aka")
  call alloc(benkc,                nele_new,       zero,        "benkc")
  call alloc(norrfamp,             nele_new, mmul, zero,        "norrfamp")
  call alloc(norrfph,              nele_new, mmul, zero,        "norrfph")
  call alloc(skrfamp,              nele_new, mmul, zero,        "skrfamp")
  call alloc(skrfph,               nele_new, mmul, zero,        "skrfph")
  call alloc(freq_rfm,             nele_new,       zero,        "freq_rfm")
  call alloc(r00,                  nele_new,       zero,        "r00")
  call alloc(scalemu,              nele_new,        one,        "scalemu")
  call alloc(irm,                  nele_new,       0,           "irm")
  call alloc(irm_rf,               nele_new,       0,           "irm_rf")
  call alloc(nmu,                  nele_new,       0,           "nmu")
  call alloc(nmu_rf,               nele_new,       0,           "nmu_rf")
  call alloc(bezr,    mNameLen, 3, nele_new,       str_nmSpace, "bezr")
  call alloc(kpa,                  nele_new,       0,           "kpa")
  call alloc(bez,     mNameLen,    nele_new,       str_nmSpace, "bez")
! call alloc(bezb,    mNameLen,    nele_new,       str_nmSpace, "bezb")
  call alloc(bezl,    mNameLen,    nele_new,       str_nmSpace, "bezl")
  call alloc(ncororb,              nele_new,       0,           "ncororb")
  call alloc(ratioe,               nele_new,       one,         "ratioe")
  call alloc(iratioe,              nele_new,       0,           "iratioe")
  call alloc(isea,                 nele_new,       0,           "isea")
  call alloc(dki,                  nele_new, 3,    zero,        "dki")
  call alloc(parbe,                nele_new, 18,   zero,        "parbe")
  call alloc(ptnfac,               nele_new,       zero,        "ptnfac")
  call alloc(imtr,                 nele_new,       0,           "imtr")
  call alloc(acdipph,              nele_new,       zero,        "acdipph")
  call alloc(nturn1,               nele_new,       0,           "nturn1")
  call alloc(nturn2,               nele_new,       0,           "nturn2")
  call alloc(nturn3,               nele_new,       0,           "nturn3")
  call alloc(nturn4,               nele_new,       0,           "nturn4")
  call alloc(crabph,               nele_new,       zero,        "crabph")
  call alloc(crabph2,              nele_new,       zero,        "crabph2")
  call alloc(crabph3,              nele_new,       zero,        "crabph3")
  call alloc(crabph4,              nele_new,       zero,        "crabph4")

  call alloc(bezb,    mNameLen,    nblo_new,       str_nmSpace, "bezb")
  call alloc(elbe,                 nblo_new,       zero,        "elbe")
  call alloc(mel,                  nblo_new,       0,           "mel")
  call alloc(mtyp,                 nblo_new, nelb, 0,           "mtyp")
  call alloc(mstr,                 nblo_new,       0,           "mstr")
  call alloc(bl1,                  nblo_new, 2, 6, zero,        "bl1")
  call alloc(bl2,                  nblo_new, 2, 6, zero,        "bl2")

  call alloc(ic,                   nblz_new,       0,           "ic")
  call alloc(mzu,                  nblz_new,       0,           "mzu")
  call alloc(imbb,                 nblz_new,       0,           "imbb")
  call alloc(icext,                nblz_new,       0,           "icext")
  call alloc(icextal,              nblz_new,       0,           "icextal")
! call alloc(exterr,               nblz_new, 40,   zero,        "exterr")   ! Replaced by compact array in mod_fluc
! call alloc(extalign,             nblz_new, 3,    zero,        "extalign") ! Replaced by compact array in mod_fluc
  call alloc(tiltc,                nblz_new,       one,         "tiltc")
  call alloc(tilts,                nblz_new,       zero,        "tilts")
  call alloc(xsi,                  nblz_new,       zero,        "xsi")
  call alloc(zsi,                  nblz_new,       zero,        "zsi")
  call alloc(smi,                  nblz_new,       zero,        "smi")
  call alloc(smizf,                nblz_new,       zero,        "smizf")
  call alloc(aaiv,          mmul,  nblz_new,       zero,        "aaiv")
  call alloc(bbiv,          mmul,  nblz_new,       zero,        "bbiv")
  call alloc(amultip,       mmul,  nblz_new,       zero,        "amultip")
  call alloc(bmultip,       mmul,  nblz_new,       zero,        "bmultip")
  call alloc(dcum,                 nblz_new+1,     zero,        "dcum", 0)
  call alloc(sigmoff,              nblz_new,       zero,        "sigmoff")

  call alloc(nnumxv,               npart_new,      0,           "nnumxv")
  call alloc(track6d, 6,           npart_new,      zero,        "track6d")

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
  real(kind=fPrec),        save :: preda
  integer,                 save :: idial,nord,nvar,nvar2,nsix,ncor,ipar(mcor)

  ! common /norf/
  integer,                 save :: nordf,nvarf,nord1,ndimf,idptr,inorm,imod1,imod2

end module mod_commond

! ================================================================================================ !
!  TRACKING COMMON VARIABLES
!  Last modified: 2018-06-21
! ================================================================================================ !
module mod_commont

  use parpro
  use floatPrecision

  implicit none

  ! common /tra/
  real(kind=fPrec), save :: x(mpa,2)
  real(kind=fPrec), save :: y(mpa,2)
  real(kind=fPrec), save :: amp(2)
  real(kind=fPrec), save :: bet0(2)
  real(kind=fPrec), save :: alf0(2)
  real(kind=fPrec), save :: clo(2)
  real(kind=fPrec), save :: clop(2)

  ! common /chrom/
  real(kind=fPrec), save :: cro(2)
  integer,          save :: is(2)
  integer,          save :: ichrom

  ! common /tasm/
  real(kind=fPrec), save :: tasm(6,6)

  ! common /track/
  integer,                       save :: nwri
  integer,          allocatable, save :: ktrack(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: strack(:)  ! (nblz)
  real(kind=fPrec), allocatable, save :: strackc(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: stracks(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: strackx(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: strackz(:) ! (nblz)
  real(kind=fPrec), allocatable, save :: dpsv1(:)   ! (npart)

  ! Substitute variables for x,y and is for DA version
  real(kind=fPrec), save :: xxtr(mpa,2)
  real(kind=fPrec), save :: yytr(mpa,2)
  integer,          save :: issss(2)

contains

subroutine mod_commont_expand_arrays(nblz_new,npart_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: nblz_new
  integer, intent(in) :: npart_new

  call alloc(ktrack,  nblz_new,  0,    "ktrack")
  call alloc(strack,  nblz_new,  zero, "strack")
  call alloc(strackc, nblz_new,  zero, "strackc")
  call alloc(stracks, nblz_new,  zero, "stracks")
  call alloc(strackx, nblz_new,  zero, "strackx")
  call alloc(strackz, nblz_new,  zero, "strackz")

  call alloc(dpsv1,   npart_new, zero, "dpsv1")

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
  real(kind=fPrec), allocatable, save :: ekv(:,:)     ! (npart,nele)
  real(kind=fPrec), allocatable, save :: smiv(:)      ! (nblz)
  real(kind=fPrec), allocatable, save :: zsiv(:)      ! (nblz)
  real(kind=fPrec), allocatable, save :: xsiv(:)      ! (nblz)

  real(kind=fPrec), allocatable, save :: fokqv(:)     ! (npart)
  real(kind=fPrec), allocatable, save :: xsv(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: zsv(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: xv1(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: yv1(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: xv2(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: yv2(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: dam(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: ekkv(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: sigmv(:)     ! (npart)
  real(kind=fPrec), allocatable, save :: dpsv(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: dp0v(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: sigmv6(:)    ! (npart)
  real(kind=fPrec), allocatable, save :: dpsv6(:)     ! (npart)
  real(kind=fPrec), allocatable, save :: ejv(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: ejfv(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: xlv(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: zlv(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: rvv(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: ejf0v(:)     ! (npart)

  integer,          allocatable, save :: numxv(:)     ! (npart)
  integer,          allocatable, save :: nms(:)       ! (npart)
  integer,          allocatable, save :: partID(:)    ! (npart)
  integer,          allocatable, save :: parentID(:)  ! (npart)

  logical,          allocatable, save :: pstop(:)     ! (npart)
  logical,          allocatable, save :: llostp(:)    ! (npart)

  real(kind=fPrec),              save :: qw(2)
  real(kind=fPrec),              save :: qwc(3)
  real(kind=fPrec),              save :: clo0(2)
  real(kind=fPrec),              save :: clop0(2)
  real(kind=fPrec),              save :: eps(2)
  real(kind=fPrec),              save :: epsa(2)
  real(kind=fPrec),              save :: ekk(2)
  real(kind=fPrec),              save :: cr(mmul)
  real(kind=fPrec),              save :: ci(mmul)
  real(kind=fPrec),              save :: pttemp
  real(kind=fPrec),              save :: temptr(6)
  integer,                       save :: kxxa

  ! common /main2/
  real(kind=fPrec), allocatable, save :: dpd(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: dpsq(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: fok(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: rho(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: fok1(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: si(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: co(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: g(:)         ! (npart)
  real(kind=fPrec), allocatable, save :: gl(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: sm1(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: sm2(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: sm3(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: sm12(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: as3(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: as4(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: as6(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: sm23(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: rhoc(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: siq(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: aek(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: afok(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: hp(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: hm(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: hc(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: hs(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: wf(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: wfa(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: wfhi(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: rhoi(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: hi(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: fi(:)        ! (npart)
  real(kind=fPrec), allocatable, save :: hi1(:)       ! (npart)
  real(kind=fPrec), allocatable, save :: dpsvl(:)     ! (npart)
  real(kind=fPrec), allocatable, save :: oidpsv(:)    ! (npart)
  real(kind=fPrec), allocatable, save :: sigmvl(:)    ! (npart)
  real(kind=fPrec), allocatable, save :: ejvl(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: ampv(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: xvl(:,:)     ! (2,npart)
  real(kind=fPrec), allocatable, save :: yvl(:,:)     ! (2,npart)
  real(kind=fPrec), allocatable, save :: aperv(:,:)   ! (npart,2)

  integer,          allocatable, save :: iv(:)        ! (npart)
  integer,          allocatable, save :: ixv(:)       ! (npart)

  ! common /main3/
  real(kind=fPrec), allocatable, save :: hv(:,:,:,:)   ! (6,2,npart,nblo)
  real(kind=fPrec), allocatable, save :: bl1v(:,:,:,:) ! (6,2,npart,nblo)
  real(kind=fPrec), allocatable, save :: tasau(:,:,:)  ! (npart,6,6)
  real(kind=fPrec), allocatable, save :: tas(:,:,:)    ! (npart,6,6)
  real(kind=fPrec), allocatable, save :: clo6v(:,:)    ! (3,npart)
  real(kind=fPrec), allocatable, save :: clop6v(:,:)   ! (3,npart)
  real(kind=fPrec), allocatable, save :: qwcs(:,:)     ! (npart,3)
  real(kind=fPrec), allocatable, save :: di0xs(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: di0zs(:)      ! (npart)
  real(kind=fPrec), allocatable, save :: dip0xs(:)     ! (npart)
  real(kind=fPrec), allocatable, save :: dip0zs(:)     ! (npart)
  real(kind=fPrec),              save :: di0au(4)
  real(kind=fPrec),              save :: tau(6,6)
  real(kind=fPrec),              save :: wx(3)
  real(kind=fPrec),              save :: x1(6)
  real(kind=fPrec),              save :: x2(6)
  real(kind=fPrec),              save :: fake(2,20)
  real(kind=fPrec),              save :: xau(2,6)
  real(kind=fPrec),              save :: cloau(6)

  ! common /main4/
  integer,          save :: numx
  real(kind=fPrec), save :: e0f
  logical,          save :: sythckcr ! Only used for CR

contains

subroutine mod_commonmn_expand_arrays(nblz_new,npart_new)

  use mod_alloc
  use numerical_constants, only : zero, one

  implicit none

  integer, intent(in) :: nblz_new
  integer, intent(in) :: npart_new


  call alloc(smiv,             nblz_new,       zero,    "smiv")
  call alloc(zsiv,             nblz_new,       zero,    "zsiv")
  call alloc(xsiv,             nblz_new,       zero,    "xsiv")

  call alloc(fokqv,            npart_new,      zero,    "fokqv")
  call alloc(xsv,              npart_new,      zero,    "xsv")
  call alloc(zsv,              npart_new,      zero,    "zsv")
  call alloc(xv1,              npart_new,      zero,    "xv1")
  call alloc(yv1,              npart_new,      zero,    "yv1")
  call alloc(xv2,              npart_new,      zero,    "xv2")
  call alloc(yv2,              npart_new,      zero,    "yv2")
  call alloc(dam,              npart_new,      zero,    "dam")
  call alloc(ekkv,             npart_new,      zero,    "ekkv")
  call alloc(sigmv,            npart_new,      zero,    "sigmv")
  call alloc(dpsv,             npart_new,      zero,    "dpsv")
  call alloc(dp0v,             npart_new,      zero,    "dp0v")
  call alloc(sigmv6,           npart_new,      zero,    "sigmv6")
  call alloc(dpsv6,            npart_new,      zero,    "dpsv6")
  call alloc(ejv,              npart_new,      zero,    "ejv")
  call alloc(ejfv,             npart_new,      zero,    "ejfv")
  call alloc(xlv,              npart_new,      zero,    "xlv")
  call alloc(zlv,              npart_new,      zero,    "zlv")
  call alloc(rvv,              npart_new,      one,     "rvv")
  call alloc(ejf0v,            npart_new,      zero,    "ejf0v")
  call alloc(numxv,            npart_new,      0,       "numxv")
  call alloc(nms,              npart_new,      0,       "nms")
  call alloc(partID,           npart_new,      0,       "partID")
  call alloc(parentID,         npart_new,      0,       "parentID")
  call alloc(pstop,            npart_new,      .false., "pstop")
  call alloc(llostp,           npart_new,      .false., "llostp")

  call alloc(dpd,              npart_new,      zero,    "dpd")
  call alloc(dpsq,             npart_new,      zero,    "dpsq")
  call alloc(fok,              npart_new,      zero,    "fok")
  call alloc(rho,              npart_new,      zero,    "rho")
  call alloc(fok1,             npart_new,      zero,    "fok1")
  call alloc(si,               npart_new,      zero,    "si")
  call alloc(co,               npart_new,      zero,    "co")
  call alloc(g,                npart_new,      zero,    "g")
  call alloc(gl,               npart_new,      zero,    "gl")
  call alloc(sm1,              npart_new,      zero,    "sm1")
  call alloc(sm2,              npart_new,      zero,    "sm2")
  call alloc(sm3,              npart_new,      zero,    "sm3")
  call alloc(sm12,             npart_new,      zero,    "sm12")
  call alloc(as3,              npart_new,      zero,    "as3")
  call alloc(as4,              npart_new,      zero,    "as4")
  call alloc(as6,              npart_new,      zero,    "as6")
  call alloc(sm23,             npart_new,      zero,    "sm23")
  call alloc(rhoc,             npart_new,      zero,    "rhoc")
  call alloc(siq,              npart_new,      zero,    "siq")
  call alloc(aek,              npart_new,      zero,    "aek")
  call alloc(afok,             npart_new,      zero,    "afok")
  call alloc(hp,               npart_new,      zero,    "hp")
  call alloc(hm,               npart_new,      zero,    "hm")
  call alloc(hc,               npart_new,      zero,    "hc")
  call alloc(hs,               npart_new,      zero,    "hs")
  call alloc(wf,               npart_new,      zero,    "wf")
  call alloc(wfa,              npart_new,      zero,    "wfa")
  call alloc(wfhi,             npart_new,      zero,    "wfhi")
  call alloc(rhoi,             npart_new,      zero,    "rhoi")
  call alloc(hi,               npart_new,      zero,    "hi")
  call alloc(fi,               npart_new,      zero,    "fi")
  call alloc(hi1,              npart_new,      zero,    "hi1")
  call alloc(dpsvl,            npart_new,      zero,    "dpsvl")
  call alloc(oidpsv,           npart_new,      one,     "oidpsv")
  call alloc(sigmvl,           npart_new,      zero,    "sigmvl")
  call alloc(ejvl,             npart_new,      zero,    "ejvl")
  call alloc(ampv,             npart_new,      zero,    "ampv")
  call alloc(xvl,       2,     npart_new,      zero,    "xvl")
  call alloc(yvl,       2,     npart_new,      zero,    "yvl")
  call alloc(aperv,            npart_new, 2,   zero,    "aperv")
  call alloc(iv,               npart_new,      0,       "iv")
  call alloc(ixv,              npart_new,      0,       "ixv")

  call alloc(tasau,            npart_new, 6,6, zero,    "tasau")
  call alloc(tas,              npart_new, 6,6, zero,    "tas")
  call alloc(clo6v,     3,     npart_new,      zero,    "clo6v")
  call alloc(clop6v,    3,     npart_new,      zero,    "clop6v")
  call alloc(qwcs,             npart_new, 3,   zero,    "qwcs")
  call alloc(di0xs,            npart_new,      zero,    "di0xs")
  call alloc(di0zs,            npart_new,      zero,    "di0zs")
  call alloc(dip0xs,           npart_new,      zero,    "dip0xs")
  call alloc(dip0zs,           npart_new,      zero,    "dip0zs")

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
  integer, save      :: lda = -1
  integer, save      :: lst = -1
  integer, save      :: lea = -1
  integer, save      :: lia = -1
  integer, save      :: lno = -1
  integer, parameter :: lnv = 40

  integer,                       save :: ndat,nda,ndamaxi,lfi
  integer,                       save :: nst,nomax,nvmax,nmmax,nocut
  integer,          allocatable, save :: idano(:),idanv(:),idapo(:)
  integer,          allocatable, save :: idalm(:),idall(:)
  integer,          allocatable, save :: i1(:),i2(:)
  integer,          allocatable, save :: ie1(:),ie2(:),ieo(:),ifi(:)
  integer,          allocatable, save :: ia1(:),ia2(:)
  logical,          allocatable, save :: allvec(:)
  character(len=:), allocatable, save :: daname(:)
  real(kind=fPrec), allocatable, save :: cc(:),facint(:)
  real(kind=fPrec),              save :: eps,epsmac

  ! dascr variables
  integer,           save :: idao
  integer,           save :: iscrda(100)
  integer,           save :: iscrri(100)
  real(kind=fPrec),  save :: rscrri(100)

contains

subroutine mld_allocArrays(da_version)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  logical, intent(in) :: da_version

  if(da_version) then
    lda = 10000
    lst = 20050000
    lea = 100000
    lia = 5000000
    lno = 120
  else
    lda = 10000
    lst = 200000
    lea = 500
    lia = 10000
    lno = 120
  end if

  call alloc(idano,lda,0,"idano")
  call alloc(idanv,lda,0,"idanv")
  call alloc(idapo,lda,0,"idapo")
  call alloc(idalm,lda,0,"idalm")
  call alloc(idall,lda,0,"idall")

  call alloc(i1,lst,0,"i1")
  call alloc(i2,lst,0,"i2")

  call alloc(ie1,lea,0,"ie1")
  call alloc(ie2,lea,0,"ie2")
  call alloc(ieo,lea,0,"ieo")
  call alloc(ifi,lea,0,"ifi")

  call alloc(ia1,lia,0,"ia1",0)
  call alloc(ia2,lia,0,"ia2",0)

  call alloc(allvec,   lda,.false.,     "allvec")
  call alloc(daname,10,lda,"          ","daname")

  call alloc(cc,    lst,zero,"cc")
  call alloc(facint,lno,zero,"facint",0)

end subroutine mld_allocArrays

end module mod_lie_dab
