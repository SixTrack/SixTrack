!-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----
!-----                                                                   -----
!-----    NEW BLOCKS PROVIDED FOR THE COLLIMATION STUDIES VIA SIXTRACK   -----
!-----                                                                   -----
!-----        G. ROBERT-DEMOLAIZE, October 27th, 2004                    -----
!-----                                                                   -----
!-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----

! Useful names for future reference in the code:
!GRD Guillaume Robert-Demolaize
!RA  Ralph Assmann
!TW  Thomas Weiler
!CB  Chiara Bracco
!SR  Stefano Redaelli
!RB  Roderik Bruce
!DM  Daniele Mirarchi
!YIL Yngve Inntjore Levinsen
!CT  Claudia Tambasco

module collimation

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use mod_hions
  use mod_alloc
  use mod_units
!  use mod_ranecu
  use mod_ranlux
  use collimation_db

#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif

  implicit none

!  private

!+cd collpara
  integer, parameter :: max_ncoll  = 100
  !integer, parameter :: maxn       = 20000
  integer, parameter :: numeff     = 32
  integer, parameter :: numeffdpop = 29
  integer, parameter :: nc         = 32

!+cd collMatNum
! EQ 2016 added variables for collimator material numbers
  integer, parameter :: nmat  = 14
  integer, parameter :: nrmat = 12

!+cd database
!GRD THIS BLOC IS COMMON TO MAINCR, DATEN, TRAUTHIN AND THIN6D
  logical, save :: do_coll = .false.
  logical, save :: do_select
  logical, save :: do_nominal
  logical, save :: dowrite_dist
  logical, save :: do_oneside
  logical, save :: dowrite_impact
  logical, save :: dowrite_secondary
  logical, save :: dowrite_amplitude
  logical, save :: radial
  logical, save :: systilt_antisymm
  logical, save :: dowritetracks
  logical, save :: cern
  logical, save :: do_mingap

!SEPT2005 for slicing process
  integer, save :: nloop
  integer, save :: rnd_seed
  integer, save :: c_offsettilt_seed
  integer, save :: ibeam
  integer, save :: jobnumber
  integer, save :: do_thisdis
  integer, save :: n_slices
  integer, save :: pencil_distr

  real(kind=fPrec), save :: myenom,mynex,mdex,myney,mdey,                    &
  &nsig_tcp3,nsig_tcsg3,nsig_tcsm3,nsig_tcla3,                       &
  &nsig_tcp7,nsig_tcsg7,nsig_tcsm7,nsig_tcla7,nsig_tclp,nsig_tcli,   &
  &nsig_tcth1,nsig_tcth2,nsig_tcth5,nsig_tcth8,                      &
  &nsig_tctv1,nsig_tctv2,nsig_tctv5,nsig_tctv8,                      &
  &nsig_tcdq,nsig_tcstcdq,nsig_tdi,nsig_tcxrp,nsig_tcryo,            &
!SEPT2005 add these lines for the slicing procedure
  &smin_slices,smax_slices,recenter1,recenter2,                      &
  &fit1_1,fit1_2,fit1_3,fit1_4,fit1_5,fit1_6,ssf1,                   &
  &fit2_1,fit2_2,fit2_3,fit2_4,fit2_5,fit2_6,ssf2,                   &
!SEPT2005,OCT2006 added offset
  &xbeat,xbeatphase,ybeat,ybeatphase,                                &
  &c_rmstilt_prim,c_rmstilt_sec,c_systilt_prim,c_systilt_sec,        &
  &c_rmsoffset_prim,c_rmsoffset_sec,c_sysoffset_prim,                &
  &c_sysoffset_sec,c_rmserror_gap,ndr,                            &
  &driftsx,driftsy,pencil_offset,pencil_rmsx,pencil_rmsy,            &
  &sigsecut3,sigsecut2,enerror,bunchlength

  ! From collimation_comnul
  real(kind=fPrec), public, save :: emitnx0_dist = zero
  real(kind=fPrec), public, save :: emitny0_dist = zero
  real(kind=fPrec), public, save :: emitnx0_collgap = zero
  real(kind=fPrec), public, save :: emitny0_collgap = zero

  real(kind=fPrec), private, save :: nr

  character(len=mNameLen), save :: name_sel
  character(len=80), save :: coll_db
  character(len=16), save :: castordir
  character(len=80), save :: filename_dis

!  common /grd/ myenom,mynex,mdex,myney,mdey,                        &
!  &nsig_tcp3,nsig_tcsg3,nsig_tcsm3,nsig_tcla3,                       &
!  &nsig_tcp7,nsig_tcsg7,nsig_tcsm7,nsig_tcla7,nsig_tclp,nsig_tcli,   &
!  &nsig_tcth1,nsig_tcth2,nsig_tcth5,nsig_tcth8,                      &
!  &nsig_tctv1,nsig_tctv2,nsig_tctv5,nsig_tctv8,                      &
!  &nsig_tcdq,nsig_tcstcdq,nsig_tdi,nsig_tcxrp,nsig_tcryo,            &
!  &smin_slices,smax_slices,recenter1,recenter2,                      &
!  &fit1_1,fit1_2,fit1_3,fit1_4,fit1_5,fit1_6,ssf1,                   &
!  &fit2_1,fit2_2,fit2_3,fit2_4,fit2_5,fit2_6,ssf2,                   &
!  &emitnx0_dist,emitny0_dist,emitnx0_collgap,emitny0_collgap,        &
!  &xbeat,xbeatphase,ybeat,ybeatphase,                                &
!  &c_rmstilt_prim,c_rmstilt_sec,c_systilt_prim,c_systilt_sec,        &
!  &c_rmsoffset_prim,c_rmsoffset_sec,c_sysoffset_prim,                &
!  &c_sysoffset_sec,c_rmserror_gap,nr,                                &
!  &ndr,driftsx,driftsy,pencil_offset,pencil_rmsx,pencil_rmsy,        &
!  &sigsecut3,sigsecut2,enerror,                                      &
!  &bunchlength,coll_db,name_sel,                                     &
!  &castordir,filename_dis,nloop,rnd_seed,c_offsettilt_seed,          &
!  &ibeam,jobnumber,do_thisdis,n_slices,pencil_distr,                 &
!  &do_coll,                                                          &
!  &do_select,do_nominal,dowrite_dist,do_oneside,dowrite_impact,      &
!  &dowrite_secondary,dowrite_amplitude,radial,systilt_antisymm,      &
!  &dowritetracks,cern,cdb_doNSig,do_mingap
!+cd info
  integer, save :: ie, iturn, nabs_total
!  common  /info/ ie,iturn,nabs_total

 !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!+cd dbcommon
! THIS BLOCK IS COMMON TO BOTH THIN6D, BEAMGAS, AND TRAUTHIN SUBROUTINES

  integer ieff,ieffdpop

  ! From collimation_comnul
  real(kind=fPrec), save :: myemitx0_dist = zero
  real(kind=fPrec), save :: myemity0_dist = zero
  real(kind=fPrec), save :: myemitx0_collgap = zero
  real(kind=fPrec), save :: myemity0_collgap = zero

  real(kind=fPrec), save :: myalphay, mybetay, myalphax, mybetax, rselect, myemitx
! myemitx was not saved?
!  common /ralph/ myemitx0_dist,myemity0_dist,myemitx0_collgap,myemity0_collgap,myalphax,myalphay,mybetax,mybetay,rselect

! M. Fiascaris for the collimation team
! variables for global inefficiencies studies
! of normalized and off-momentum halo
! Last modified: July 2016

  real(kind=fPrec), allocatable, save :: neff(:) !(numeff)
  real(kind=fPrec), allocatable, save :: rsig(:) !(numeff)
!  common  /eff/ neff,rsig

  integer, allocatable, save :: counteddpop(:,:) !(npart,numeffdpop)
  integer, allocatable, save :: npartdpop(:) !(numeffdpop)
  integer, allocatable, save :: counted2d(:,:,:) !(npart,numeff,numeffdpop)
  real(kind=fPrec), allocatable, save :: neffdpop(:) !(numeffdpop)
  real(kind=fPrec), allocatable, save :: dpopbins(:) !(numeffdpop)
!  common  /effdpop/ neffdpop,dpopbins,npartdpop,counteddpop

  real(kind=fPrec) dpopmin,dpopmax,mydpop
  real(kind=fPrec), allocatable, save :: neff2d(:,:) !(numeff,numeffdpop)
!  common /eff2d/ neff2d

  integer, allocatable, save :: nimpact(:) !(50)
  real(kind=fPrec), allocatable, save :: sumimpact(:) !(50)
  real(kind=fPrec), allocatable, save :: sqsumimpact(:) !(50)
!  common  /rimpact/ sumimpact,sqsumimpact,nimpact

  character(len=:), allocatable, save :: ename(:) !(mNameLen)(nblz)
  integer, allocatable, save :: nampl(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sum_ax(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sqsum_ax(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sum_ay(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sqsum_ay(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sampl(:) !(nblz)
!  common  /ampl_rev/ sum_ax,sqsum_ax,sum_ay,sqsum_ay,sampl,ename,nampl

  real(kind=fPrec), allocatable, save :: neffx(:) !(numeff)
  real(kind=fPrec), allocatable, save :: neffy(:) !(numeff)
!  common /efficiency/ neffx,neffy

  integer, allocatable, save :: secondary(:) !(npart)
  integer, allocatable, save :: tertiary(:) !(npart)
  integer, allocatable, save :: other(:) !(npart)
  integer, allocatable, save :: scatterhit(:) !(npart)
  integer, allocatable, save :: part_hit_before_pos(:) !(npart)
  integer, allocatable, save :: part_hit_before_turn(:) !(npart)

  real(kind=fPrec), allocatable, save :: part_indiv(:) !(npart)
  real(kind=fPrec), allocatable, save :: part_linteract(:) !(npart)

  integer, allocatable, save :: part_hit_pos(:) !(npart)
  integer, allocatable, save :: part_hit_turn(:) !(npart)
  integer, allocatable, save :: part_abs_pos(:) !(npart)
  integer, allocatable, save :: part_abs_turn(:) !(npart)
  integer, allocatable, save :: part_select(:) !(npart)
  integer, allocatable, save :: nabs_type(:) !(npart)
  integer, save :: n_tot_absorbed
  integer, save :: n_absorbed

  real(kind=fPrec), allocatable, save :: part_impact(:) !(npart)
!  common /stats/ part_impact, part_hit_pos,part_hit_turn, part_hit_before_pos, part_hit_before_turn, &
!  & part_abs_pos,part_abs_turn, nabs_type,part_indiv, part_linteract,secondary,tertiary,other,scatterhit

!  common /n_tot_absorbed/ n_tot_absorbed,n_absorbed
!  common /part_select/ part_select

!  logical firstrun
!  common /firstrun/ firstrun

  integer, save :: nsurvive, nsurvive_end, num_selhit, n_impact
!  common /outcoll/ nsurvive,num_selhit,n_impact,nsurvive_end

  integer, save :: napx00

  real(kind=fPrec), allocatable, save :: db_tilt(:,:) !(max_ncoll,2)

  integer, allocatable, save :: cn_impact(:)  !(max_ncoll)
  integer, allocatable, save :: cn_absorbed(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: caverage(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: csigma(:) !(max_ncoll)
!  common /collsummary/ caverage,csigma,cn_impact,cn_absorbed

! Change the following block to npart
! This is the array that the generated distribution is placed into
  real(kind=fPrec), allocatable, save :: myx(:) !(npart)
  real(kind=fPrec), allocatable, save :: myxp(:) !(npart)
  real(kind=fPrec), allocatable, save :: myy(:) !(npart)
  real(kind=fPrec), allocatable, save :: myyp(:) !(npart)
  real(kind=fPrec), allocatable, save :: myp(:) !(npart)
  real(kind=fPrec), allocatable, save :: mys(:) !(npart)
!  common /coord/ myx,myxp,myy,myyp,myp,mys

  integer, allocatable, save :: counted_r(:,:) !(npart,numeff)
  integer, allocatable, save :: counted_x(:,:) !(npart,numeff)
  integer, allocatable, save :: counted_y(:,:) !(npart,numeff)
!  common /counting/ counted_r,counted_x,counted_y

  character(len=4), save :: smpl
  character(len=80), save :: pfile
!  common /samplenumber/ pfile,smpl,samplenumber
!
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!+cd dblinopt
!
! THIS BLOCK IS COMMON TO WRITELIN,LINOPT,TRAUTHIN,THIN6D AND MAINCR
!
  real(kind=fPrec), allocatable, save :: tbetax(:)  !(nblz)
  real(kind=fPrec), allocatable, save :: tbetay(:)  !(nblz)
  real(kind=fPrec), allocatable, save :: talphax(:) !(nblz)
  real(kind=fPrec), allocatable, save :: talphay(:) !(nblz)
  real(kind=fPrec), allocatable, save :: torbx(:)   !(nblz)
  real(kind=fPrec), allocatable, save :: torbxp(:)  !(nblz)
  real(kind=fPrec), allocatable, save :: torby(:)   !(nblz)
  real(kind=fPrec), allocatable, save :: torbyp(:)  !(nblz)
  real(kind=fPrec), allocatable, save :: tdispx(:)  !(nblz)
  real(kind=fPrec), allocatable, save :: tdispy(:)  !(nblz)

!  common /rtwiss/ tbetax,tbetay,talphax,talphay,torbx,torbxp,torby,torbyp,tdispx,tdispy
!
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
! Variables for finding the collimator with the smallest gap
! and defining, stroring the gap rms error
!
  character(len=mNameLen) :: coll_mingap1, coll_mingap2
  real(kind=fPrec), allocatable, save :: gap_rms_error(:) !(max_ncoll)
  real(kind=fPrec) :: nsig_err, sig_offset
  real(kind=fPrec) :: mingap, gap_h1, gap_h2, gap_h3, gap_h4
  integer :: coll_mingap_id

! common /gap_err/ gap_rms_error
!
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!+cd dbdaten

! IN "+CD DBTRTHIN", "+CD DBDATEN" and "+CD DBTHIN6D"
! logical cut_input
! common /cut/ cut_input

! IN "+CD DBTRTHIN" and "+CD DBDATEN"
  real(kind=fPrec), save :: remitx_dist,remity_dist,remitx_collgap,remity_collgap
! common  /remit/ remitx_dist, remity_dist,remitx_collgap,remity_collgap
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
  logical, save :: coll_found(100)


!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!+cd dbcolcom
  logical, save :: firstcoll,found,onesided
  integer rnd_lux,rnd_k1,rnd_k2

  integer, save :: myix,myktrack

  real(kind=fPrec) nspx,nspy,mux0,muy0
  real(kind=fPrec) ax0,ay0,bx0,by0
  real(kind=fPrec), save :: totals

  ! IN "+CD DBTRTHIN", "+CD DBDATEN" and "+CD DBTHIN6D"
!  logical cut_input
!  common /cut/ cut_input

  real(kind=fPrec), allocatable, save :: xbob(:) !(nblz)
  real(kind=fPrec), allocatable, save :: ybob(:) !(nblz)
  real(kind=fPrec), allocatable, save :: xpbob(:) !(nblz)
  real(kind=fPrec), allocatable, save :: ypbob(:) !(nblz)

  real(kind=fPrec), allocatable, save :: xineff(:) !(npart)
  real(kind=fPrec), allocatable, save :: yineff(:) !(npart)
  real(kind=fPrec), allocatable, save :: xpineff(:) !(npart)
  real(kind=fPrec), allocatable, save :: ypineff(:) !(npart)

!  common /xcheck/ xbob,ybob,xpbob,ypbob,xineff,yineff,xpineff,ypineff

  real(kind=fPrec), allocatable, save :: mux(:) !(nblz)
  real(kind=fPrec), allocatable, save :: muy(:) !(nblz)
!  common /mu/ mux,muy

!  common /collocal/ myix,myktrack,totals,firstcoll,found,onesided

! common /icoll/  icoll
!
!
!  common /materia/mat
!  common /phase/x,xp,z,zp,dpop
!  common /nommom/p0
!  common /cjaw1/zlm
! END BLOCK DBCOLLIM


!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!+cd dbpencil
!
! THIS BLOCK IS COMMON TO THIN6D, TRAUTHIN, COLLIMATE32 AND MAINCR
!
  integer, save :: ipencil
  real(kind=fPrec), save :: xp_pencil0,yp_pencil0
  real(kind=fPrec), allocatable, save :: x_pencil(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: y_pencil(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: pencil_dx(:) !(max_ncoll)
!  common  /pencil/  xp_pencil0,yp_pencil0,pencil_dx,ipencil
!  common  /pencil2/ x_pencil, y_pencil
!
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!+cd dbmkdist
!++ Vectors of coordinates

!  integer :: i,j,mynp,nloop
!  real(kind=fPrec) :: myx(maxn) !(maxn)
!  real(kind=fPrec) :: myxp(maxn) !(maxn)
!  real(kind=fPrec) :: myy(maxn) !(maxn)
!  real(kind=fPrec) :: myyp(maxn) !(maxn)
!  real(kind=fPrec) :: myp(maxn) !(maxn)
!  real(kind=fPrec) :: mys(maxn) !(maxn)
!
!  real(kind=fPrec) myalphax,mybetax,myemitx0,myemitx,mynex,mdex, &
!  &mygammax,myalphay,mybetay,myemity0,myemity,myney,mdey,mygammay,   &
!  &xsigmax,ysigmay,myenom,nr,ndr
!
!

! IN "+CD DBTRTHIN", "+CD DBDATEN", "+CD DBTHIN6D", and "+CD DBMKDIST"
! USED IN MULTIPLE COMMON BLOCKS
  logical, save :: cut_input
!  common /cut/ cut_input

!from +cd interac
!October 2013
!Mean excitation energy (GeV) values added by Claudia for Bethe-Bloch implementation:
!  data (exenergy(i),i=1,5)/ 63.7e-9,166e-9, 322e-9, 727e-9, 823e-9 /
!  data (exenergy(i),i=6,7)/ 78e-9, 78.0e-9 /
!  data (exenergy(i),i=8,nrmat)/ 87.1e-9, 152.9e-9, 424e-9, 320.8e-9, 682.2e-9/
  real(kind=fPrec), parameter :: exenergy(nmat) = &
 & [ 63.7e-9_fPrec, 166e-9_fPrec, 322e-9_fPrec, 727e-9_fPrec, 823e-9_fPrec, 78e-9_fPrec, 78.0e-9_fPrec, 87.1e-9_fPrec, &
 & 152.9e-9_fPrec, 424e-9_fPrec, 320.8e-9_fPrec, 682.2e-9_fPrec, zero, c1e10 ]
! common/meanexen/exenergy(nmat)

!+cd dbtrthin

! Note: no saves needed

!++ Vectors of coordinates

  real(kind=fPrec), private :: mygammax,mygammay


  ! IN "+CD DBTRTHIN" and "+CD DBDATEN"
!  real(kind=fPrec) remitx_dist,remity_dist,
! &     remitx_collgap,remity_collgap
!  common  /remit/ remitx_dist, remity_dist,
! &     remitx_collgap,remity_collgap

  integer, private :: k
!
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!

  logical, public, save :: firstrun
  integer, save :: icoll

!+cd flukavars
! RB DM 2014 added variables for FLUKA output
  real(kind=fPrec), private, save :: xInt,xpInt,yInt,ypInt,sInt
! common/flukaVars/xInt,xpInt,yInt,ypInt,sInt

!+cd funint
  real(kind=fPrec), private, save :: tftot
!  common/funint/tftot

!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!+cd dbthin6d

  integer, save :: num_surhit
  integer, save :: numbin
  integer, save :: ibin
  integer, save :: num_selabs
  integer, save :: iturn_last_hit
  integer, save :: iturn_absorbed
  integer, save :: iturn_survive
  integer, save :: imov
  integer, save :: totalelem
  integer, save :: selelem
  integer, save :: unitnumber
  integer, save :: distnumber
  integer, save :: turnnumber
  integer, private, save :: jb

! SR, 29-08-2005: add the required variable for slicing collimators
  integer, private, save :: jjj

  real(kind=fPrec), save :: zbv

  real(kind=fPrec), save :: c_length    !length in m
  real(kind=fPrec), save :: c_rotation  !rotation angle vs vertical in radian
  real(kind=fPrec), save :: c_aperture  !aperture in m
  real(kind=fPrec), save :: c_offset    !offset in m
  real(kind=fPrec), save :: c_tilt(2)   !tilt in radian
  character(len=4), save :: c_material  !material

  integer, allocatable, save :: ipart(:) !(npart)
  integer, allocatable, save :: flukaname(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: cx(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: cxp(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: cy(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: cyp(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: cp(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: cs(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcx(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcxp(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcy(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcyp(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcp(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcs(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcx0(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcxp0(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcy0(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcyp0(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcp0(:) !(npart)

  real(kind=fPrec), allocatable, save :: xgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: xpgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: ygrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: ypgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: pgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: ejfvgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: sigmvgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: rvvgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: dpsvgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: oidpsvgrd(:) !(npart)
  real(kind=fPrec), allocatable, save :: dpsv1grd(:) !(npart)

  real(kind=fPrec), save :: enom_gev,betax,betay,xmax,ymax
  real(kind=fPrec), save :: nsig,calc_aperture,gammax,gammay,gammax0,gammay0,gammax1,gammay1
  real(kind=fPrec), save :: xj,xpj,yj,ypj,pj
  real(kind=fPrec), save :: arcdx,arcbetax,xdisp,rxjco,ryjco
  real(kind=fPrec), save :: rxpjco,rypjco,c_rmstilt,c_systilt
  real(kind=fPrec), save :: scale_bx, scale_by, scale_bx0, scale_by0, xkick, ykick, bx_dist, by_dist
  real(kind=fPrec), save :: xmax_pencil, ymax_pencil, xmax_nom, ymax_nom, nom_aperture, pencil_aperture

  real(kind=fPrec), allocatable, save :: xp_pencil(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: yp_pencil(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: csum(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: csqsum(:) !(max_ncoll)

  real(kind=fPrec), save :: x_pencil0, y_pencil0, sum, sqsum
  real(kind=fPrec), save :: average, sigma, sigsecut, nspxd, xndisp, zpj

  real(kind=fPrec), save :: dnormx,dnormy,driftx,drifty,xnorm,xpnorm,xangle,ynorm,ypnorm,yangle,grdpiover2,grdpiover4,grd3piover4

!SEPT2005-SR, 29-08-2005 --- add parameter for the array length ---- TW
  real(kind=fPrec), allocatable, save :: x_sl(:) !(100)
  real(kind=fPrec), allocatable, save :: x1_sl(:) !(100)
  real(kind=fPrec), allocatable, save :: x2_sl(:) !(100)
  real(kind=fPrec), allocatable, save :: y1_sl(:) !(100)
  real(kind=fPrec), allocatable, save :: y2_sl(:) !(100)
  real(kind=fPrec), allocatable, save :: angle1(:) !(100)
  real(kind=fPrec), allocatable, save :: angle2(:) !(100)

  real(kind=fPrec), save :: max_tmp, a_tmp1, a_tmp2, ldrift, mynex2, myney2, Nap1pos,Nap2pos,Nap1neg,Nap2neg
  real(kind=fPrec), save :: tiltOffsPos1,tiltOffsPos2,tiltOffsNeg1,tiltOffsNeg2
  real(kind=fPrec), save :: beamsize1, beamsize2,betax1,betax2,betay1,betay2, alphax1, alphax2,alphay1,alphay2,minAmpl
!SEPT2005


!  common /dbthinc/ cx,cxp,cy,cyp,                                   &
!  &cp,cs,rcx,rcxp,rcy,rcyp,                                          &
!  &rcp,rcs,rcx0,rcxp0,rcy0,                                          &
!  &rcyp0,rcp0,enom_gev,betax,betay,xmax,ymax,                        &
!  &nsig,calc_aperture,gammax,gammay,gammax0,gammay0,gammax1,gammay1, &
!  &xj,xpj,yj,ypj,pj,arcdx,arcbetax,xdisp,rxjco,ryjco,                &
!  &rxpjco,rypjco,c_rmstilt,                                          &
!  &c_systilt,scale_bx,scale_by,scale_bx0,scale_by0,xkick,            &
!  &ykick,bx_dist,by_dist,xmax_pencil,ymax_pencil,xmax_nom,ymax_nom,  &
!  &nom_aperture,pencil_aperture,xp_pencil,                           &
!  &yp_pencil,x_pencil0,y_pencil0,sum,sqsum,                          &
!  &csum,csqsum,average,sigma,sigsecut,nspxd,                         &
!  &xndisp,xgrd,xpgrd,ygrd,ypgrd,zpj,                                 &
!  &pgrd,ejfvgrd,sigmvgrd,rvvgrd,                                     &
!  &dpsvgrd,oidpsvgrd,dpsv1grd,                                       &
!  &dnormx,dnormy,driftx,drifty,                                      &
!  &xnorm,xpnorm,xangle,ynorm,ypnorm,yangle,                          &
!  &grdpiover2,grdpiover4,grd3piover4,                                &
!  &x_sl,x1_sl,x2_sl,                                                 &
!  &y1_sl, y2_sl,                                                &
!  &angle1, angle2,                                              &
!  &max_tmp,                                                     &
!  &a_tmp1, a_tmp2, ldrift, mynex2, myney2,                      &
!  &Nap1pos,Nap2pos,Nap1neg,Nap2neg,                             &
!  &tiltOffsPos1,tiltOffsPos2,tiltOffsNeg1,tiltOffsNeg2,         &
!  &beamsize1, beamsize2,betax1,betax2,betay1,betay2,            &
!  &alphax1, alphax2,alphay1,alphay2,minAmpl,                    &
!  &ios,num_surhit,numbin,ibin,                                       &
!  &num_selabs,iturn_last_hit,iturn_absorbed,iturn_survive,imov,      &
!  &ipart,totalelem,selelem,unitnumber,distnumber,turnnumber,         &
!  &jb,flukaname,                                                     &
!  &jjj,ijk,zbv,c_length,c_rotation,                                  &
!  &c_aperture,c_offset,c_tilt,c_material

! myran_gauss,rndm5,

!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!+cd dbcollim
  integer, private, save :: nev
!  real(kind=fPrec), private, save :: c_xmin,c_xmax,c_xpmin,c_xpmax,c_zmin,c_zmax,c_zpmin,c_zpmax,length
!!  common /cmom/xmin,xmax,xpmin,xpmax,zmin,zmax,zpmin,zpmax,length,nev
!
!  real(kind=fPrec), private, save :: c_mybetax,c_mybetaz,mymux,mymuz,atdi
!!  common /other/mybetax,mybetaz,mymux,mymuz,atdi

  real(kind=fPrec), private, save :: length

! Common to interac and dbcollim
  integer, private, save :: mat
  real(kind=fPrec), private, save :: x,xp,z,zp,dpop
  real(kind=fPrec), private, save :: p0
  real(kind=fPrec), private, save :: zlm
!
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!+cd interac
  integer, save :: mcurr
  real(kind=fPrec), save :: xintl(nmat)
  real(kind=fPrec), save :: radl(nmat)
  real(kind=fPrec), save :: zlm1
  real(kind=fPrec), save :: xpsd
  real(kind=fPrec), save :: zpsd
  real(kind=fPrec), save :: psd
  real(kind=fPrec), save :: anuc(nmat)
  real(kind=fPrec), save :: rho(nmat)
  real(kind=fPrec), save :: emr(nmat)
  real(kind=fPrec), parameter :: tlcut = 0.0009982_fPrec
  real(kind=fPrec), save :: hcut(nmat)
  real(kind=fPrec), save :: csect(0:5,nmat)
  real(kind=fPrec), save :: csref(0:5,nmat)
  real(kind=fPrec), save :: bnref(nmat)
  real(kind=fPrec), save :: freep(nmat)
  real(kind=fPrec), save :: cprob(0:5,nmat)
  real(kind=fPrec), save :: bn(nmat)
  real(kind=fPrec), save :: bpp
  real(kind=fPrec), save :: xln15s
  real(kind=fPrec), save :: ecmsq
  real(kind=fPrec), save :: pptot
  real(kind=fPrec), save :: ppel
  real(kind=fPrec), save :: ppsd
  real(kind=fPrec), save :: pptref
  real(kind=fPrec), save :: pperef
  real(kind=fPrec), save :: pref
  real(kind=fPrec), save :: pptco
  real(kind=fPrec), save :: ppeco
  real(kind=fPrec), save :: sdcoe
  real(kind=fPrec), save :: freeco
  real(kind=fPrec), save :: zatom(nmat)
  real(kind=fPrec), save :: dpodx(nmat)

!electron density and plasma energy
  real(kind=fPrec), save :: edens(nmat)
  real(kind=fPrec), save :: pleng(nmat)

! parameter(fnavo=6.02214129e23_fPrec)
  real(kind=fPrec), save :: cgen(200,nmat)
  character(4), save :: mname(nmat)

!  common/mater/anuc(nmat),zatom(nmat),rho(nmat),emr(nmat)
!  common/coul/tlcut,hcut(nmat),cgen(200,nmat),mcurr
!  common/scat/cs(0:5,nmat),csref(0:5,nmat),bnref(nmat),freep(nmat)
!  common/scatu/cprob(0:5,nmat),bn(nmat),bpp,xln15s,ecmsq
!  common/scatu2/xintl(nmat),radl(nmat),mname
!  common/scatpp/pptot,ppel,ppsd
!  common/sppref/pptref,pperef,pref,pptco,ppeco,sdcoe,freeco
!! real(kind=fPrec) exenergy
!! common/meanexen/exenergy(nmat)
!  common/cmcs1/zlm1
!  common/sindif/xpsd,zpsd,psd
!  common/cdpodx/dpodx
!  common/cions/edens(nmat),pleng(nmat)


!  common/materia/mat
!  common/phase/x,xp,z,zp,dpop
!  common/nommom/p0
!  common/cjaw1/zlm

!>
!! block data scdata
!! Cross section inputs and material property database
!! GRD CHANGED ON 2/2003 TO INCLUDE CODE FOR C, C2 from JBJ (rwa)
!<
  integer, private :: i
! Total number of materials are defined in nmat
! Number of real materials are defined in nrmat
! The last materials in nmat are 'vacuum' and 'black',see in sub. SCATIN

! Reference data at pRef=450Gev
  data (mname(i),i=1,nrmat)/ 'Be','Al','Cu','W','Pb','C','C2','MoGR','CuCD', 'Mo', 'Glid', 'Iner'/

  data mname(nmat-1), mname(nmat)/'vacu','blac'/

!GRD IMPLEMENT CHANGES FROM JBJ, 2/2003 RWA
  data (anuc(i),i=1,5)/ 9.01d0,26.98d0,63.55d0,183.85d0,207.19d0/
  data (anuc(i),i=6,7)/12.01d0,12.01d0/
  data (anuc(i),i=8,nrmat)/13.53d0,25.24d0,95.96d0,63.15d0,166.7d0/

  data (zatom(i),i=1,5)/ 4d0, 13d0, 29d0, 74d0, 82d0/
  data (zatom(i),i=6,7)/ 6d0, 6d0/
  data (zatom(i),i=8,nrmat)/ 6.65d0, 11.9d0, 42d0, 28.8d0, 67.7d0/

  data (rho(i),i=1,5)/ 1.848d0, 2.70d0, 8.96d0, 19.3d0, 11.35d0/
  data (rho(i),i=6,7)/ 1.67d0, 4.52d0/
  data (rho(i),i=8,nrmat)/ 2.5d0, 5.4d0, 10.22d0, 8.93d0, 18d0/

  data (radl(i),i=1,5)/ 0.353d0,0.089d0,0.0143d0,0.0035d0,0.0056d0/
  data (radl(i),i=6,7)/ 0.2557d0, 0.094d0/
  data (radl(i),i=8,nrmat)/ 0.1193d0, 0.0316d0, 0.0096d0, 0.0144d0, 0.00385d0/
  data radl(nmat-1),radl(nmat)/ 1.d12, 1.d12 /

!MAY06-GRD value for Tungsten (W) not stated
  data (emr(i),i=1,5)/  0.22d0, 0.302d0, 0.366d0, 0.520d0, 0.542d0/
  data (emr(i),i=6,7)/  0.25d0, 0.25d0/
  data (emr(i),i=8,nrmat)/ 0.25d0, 0.308d0, 0.481d0, 0.418d0, 0.578d0/

  data (hcut(i),i=1,5)/0.02d0, 0.02d0, 3*0.01d0/
  data (hcut(i),i=6,7)/0.02d0, 0.02d0/
  data (hcut(i),i=8,nrmat)/0.02d0, 0.02d0, 0.02d0, 0.02d0, 0.02d0/

  data (dpodx(i),i=1,5)/ .55d0, .81d0, 2.69d0, 5.79d0, 3.4d0 /
  data (dpodx(i),i=6,7)/ .75d0, 1.5d0 /

! All cross-sections are in barns,nuclear values from RPP at 20geV
! Coulomb is integerated above t=tLcut[Gev2] (+-1% out Gauss mcs)

! in Cs and CsRef,1st index: Cross-sections for processes
! 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
! 4:Single Diffractive pp or pn, 5:Coulomb for t above mcs

! Claudia 2013: updated cross section values. Unit: Barn. New 2013:
  data csref(0,1),csref(1,1),csref(5,1)/0.271d0, 0.192d0, 0.0035d-2/
  data csref(0,2),csref(1,2),csref(5,2)/0.643d0, 0.418d0, 0.034d-2/
  data csref(0,3),csref(1,3),csref(5,3)/1.253d0, 0.769d0, 0.153d-2/
  data csref(0,4),csref(1,4),csref(5,4)/2.765d0, 1.591d0, 0.768d-2/
  data csref(0,5),csref(1,5),csref(5,5)/3.016d0, 1.724d0, 0.907d-2/
  data csref(0,6),csref(1,6),csref(5,6)/0.337d0, 0.232d0, 0.0076d-2/
  data csref(0,7),csref(1,7),csref(5,7)/0.337d0, 0.232d0, 0.0076d-2/
  data csref(0,8),csref(1,8),csref(5,8)/0.362d0, 0.247d0, 0.0094d-2/
  data csref(0,9),csref(1,9),csref(5,9)/0.572d0, 0.370d0, 0.0279d-2/
  data csref(0,10),csref(1,10),csref(5,10)/1.713d0,1.023d0,0.265d-2/
  data csref(0,11),csref(1,11),csref(5,11)/1.246d0,0.765d0,0.139d-2/
  data csref(0,12),csref(1,12),csref(5,12)/2.548d0,1.473d0,0.574d-2/

! pp cross-sections and parameters for energy dependence
  data pptref,pperef,sdcoe,pref/0.04d0,0.007d0,0.00068d0,450.0d0/
  data pptco,ppeco,freeco/0.05788d0,0.04792d0,1.618d0/

! Nuclear elastic slope from Schiz et al.,PRD 21(3010)1980
!MAY06-GRD value for Tungsten (W) not stated
  data (bnref(i),i=1,5)/74.7d0,120.3d0,217.8d0,440.3d0,455.3d0/
  data (bnref(i),i=6,7)/70.d0, 70.d0/
  data (bnref(i),i=8,nrmat)/ 76.7d0, 115.0d0, 273.9d0, 208.7d0, 392.1d0/
!GRD LAST 2 ONES INTERPOLATED

! Cprob to choose an interaction in iChoix
  data (cprob(0,i),i=1,nmat)/nmat*zero/
  data (cprob(5,i),i=1,nmat)/nmat*one/
  ! file units
  integer, private, save :: dist0_unit, survival_unit, collgaps_unit, collimator_temp_db_unit
  integer, private, save :: impact_unit, tracks2_unit, pencilbeam_distr_unit, coll_ellipse_unit, all_impacts_unit
  integer, private, save :: FLUKA_impacts_unit, FLUKA_impacts_all_unit, coll_scatter_unit, FirstImpacts_unit, RHIClosses_unit
  integer, private, save :: twisslike_unit, sigmasettings_unit, distsec_unit, efficiency_unit, efficiency_dpop_unit
  integer, private, save :: coll_summary_unit, amplitude_unit, amplitude2_unit, betafunctions_unit, orbitchecking_unit, distn_unit
  integer, private, save :: filename_dis_unit, CollPositions_unit, all_absorptions_unit, efficiency_2d_unit
  integer, private, save :: collsettings_unit, outlun
  ! These are not in use
  !integer, save :: betatron_unit, beta_beat_unit

#ifdef HDF5
  ! Variables to save hdf5 dataset indices
  integer, private, save :: coll_hdf5_survival
  integer, private, save :: coll_hdf5_allImpacts
  integer, private, save :: coll_hdf5_fstImpacts
  integer, private, save :: coll_hdf5_allAbsorb
  integer, private, save :: coll_hdf5_collScatter
#endif

contains

! General routines:
! collimate_init()
! collimate_start_sample()
! collimate_start_turn()
! collimate_start_element()
! collimate_start_collimator()
! collimate_do_collimator()
! collimate_end_collimator()
! collimate_end_element()
! collimate_end_turn()
! collimate_end_sample()
! collimate_exit()
!
! To stop a messy future, each of these should contain calls to
! implementation specific functions: e.g. collimate_init_k2(), etc.
! These should contain the "real" code.

! In addition, these files contain:
! 1: The RNG used in collimation.
! 2: A bunch distribution generator

subroutine collimation_allocate_arrays

  implicit none

  ! Initial allocation handled by expand arrays routine
  call collimation_expand_arrays(npart,nblz,nele)

  ! Fixed allocations follow:
  call alloc(gap_rms_error, max_ncoll, zero, "gap_rms_error") !(max_ncoll)
  call alloc(xp_pencil, max_ncoll, zero, "xp_pencil")
  call alloc(yp_pencil, max_ncoll, zero, "yp_pencil")
  call alloc(csum,      max_ncoll, zero, "csum")
  call alloc(csqsum,    max_ncoll, zero, "csqsum")

  call alloc(x_pencil,  max_ncoll, zero, "x_pencil") !(max_ncoll)
  call alloc(y_pencil,  max_ncoll, zero, "y_pencil") !(max_ncoll)
  call alloc(pencil_dx, max_ncoll, zero, "pencil_dx") !(max_ncoll)

  !SEPT2005-SR, 29-08-2005 --- add parameter for the array length ---- TW
  call alloc(x_sl, 100, zero, "x_sl") !(100)
  call alloc(x1_sl, 100, zero, "x1_sl") !(100)
  call alloc(x2_sl, 100, zero, "x2_sl") !(100)
  call alloc(y1_sl, 100, zero, "y1_sl") !(100)
  call alloc(y2_sl, 100, zero, "y2_sl") !(100)
  call alloc(angle1, 100, zero, "angle1") !(100)
  call alloc(angle2, 100, zero, "angle2") !(100)

  call alloc(npartdpop, numeffdpop, 0, "npartdpop") !(numeffdpop)
  call alloc(neff, numeff, zero, "neff") !(numeff)
  call alloc(rsig, numeff, zero, "rsig") !(numeff)
  call alloc(neffdpop, numeffdpop, zero, "neffdpop") !(numeffdpop)
  call alloc(dpopbins, numeffdpop, zero, "dpopbins") !(numeffdpop)
  call alloc(neff2d, numeff, numeffdpop, zero, "neff2d") !(numeff,numeffdpop)

  call alloc(nimpact, 50, 0, "nimpact") !(50)
  call alloc(sumimpact, 50, zero, "sumimpact") !(50)
  call alloc(sqsumimpact, 50, zero, "sqsumimpact") !(50)

  call alloc(neffx, numeff, zero, "neffx") !(numeff)
  call alloc(neffy, numeff, zero, "neffy") !(numeff)
  call alloc(db_tilt, max_ncoll, 2, zero, "db_tilt") !(max_ncoll,2)

  call alloc(cn_impact, max_ncoll, 0, "cn_impact")  !(max_ncoll)
  call alloc(cn_absorbed, max_ncoll, 0, "cn_absorbed") !(max_ncoll)
  call alloc(caverage, max_ncoll, zero, "caverage") !(max_ncoll)
  call alloc(csigma, max_ncoll, zero, "csigma") !(max_ncoll)

end subroutine collimation_allocate_arrays

subroutine collimation_expand_arrays(npart_new, nblz_new, nele_new)

  implicit none

  integer, intent(in) :: npart_new
  integer, intent(in) :: nblz_new
  integer, intent(in) :: nele_new

  ! Arrays that are always needed
  call alloc(part_abs_turn, npart_new, 0, "part_abs_turn") !(npart_new)

  if(.not. do_coll) return
  ! Arrays that are only needed if Collimation is enabled

  call alloc(tbetax,  nblz_new, zero, 'tbetax')  !(nblz)
  call alloc(tbetay,  nblz_new, zero, 'tbetay')  !(nblz)
  call alloc(talphax, nblz_new, zero, 'talphax') !(nblz)
  call alloc(talphay, nblz_new, zero, 'talphay') !(nblz)
  call alloc(torbx,   nblz_new, zero, 'torbx')   !(nblz)
  call alloc(torbxp,  nblz_new, zero, 'torbxp')  !(nblz)
  call alloc(torby,   nblz_new, zero, 'torby')   !(nblz)
  call alloc(torbyp,  nblz_new, zero, 'torbyp')  !(nblz)
  call alloc(tdispx,  nblz_new, zero, 'tdispx')  !(nblz)
  call alloc(tdispy,  nblz_new, zero, 'tdispy')  !(nblz)

  call alloc(flukaname, npart_new, 0, "flukaname") !(npart)
  call alloc(ipart, npart_new, 0, "ipart") !(npart)
  call alloc(cx,    npart_new, zero, "cx") !(npart)
  call alloc(cxp,   npart_new, zero, "cxp") !(npart)
  call alloc(cy,    npart_new, zero, "cy") !(npart)
  call alloc(cyp,   npart_new, zero, "cyp") !(npart)
  call alloc(cp,    npart_new, zero, "cp") !(npart)
  call alloc(cs,    npart_new, zero, "cs") !(npart)
  call alloc(rcx,   npart_new, zero, "rcx") !(npart)
  call alloc(rcxp,  npart_new, zero, "rcxp") !(npart)
  call alloc(rcy,   npart_new, zero, "rcy") !(npart)
  call alloc(rcyp,  npart_new, zero, "rcyp") !(npart)
  call alloc(rcp,   npart_new, zero, "rcp") !(npart)
  call alloc(rcs,   npart_new, zero, "rcs") !(npart)
  call alloc(rcx0,  npart_new, zero, "rcx0") !(npart)
  call alloc(rcxp0, npart_new, zero, "rcxp0") !(npart)
  call alloc(rcy0,  npart_new, zero, "rcy0") !(npart)
  call alloc(rcyp0, npart_new, zero, "rcyp0") !(npart)
  call alloc(rcp0,  npart_new, zero, "rcp0") !(npart)

  call alloc(xgrd,      npart_new, zero, "xgrd") !(npart)
  call alloc(xpgrd,     npart_new, zero, "xpgrd") !(npart)
  call alloc(ygrd,      npart_new, zero, "ygrd") !(npart)
  call alloc(ypgrd,     npart_new, zero, "ypgrd") !(npart)
  call alloc(pgrd,      npart_new, zero, "pgrd") !(npart)
  call alloc(ejfvgrd,   npart_new, zero, "ejfvgrd") !(npart)
  call alloc(sigmvgrd,  npart_new, zero, "sigmvgrd") !(npart)
  call alloc(rvvgrd,    npart_new, zero, "rvvgrd") !(npart)
  call alloc(dpsvgrd,   npart_new, zero, "dpsvgrd") !(npart)
  call alloc(oidpsvgrd, npart_new, zero, "oidpsvgrd") !(npart)
  call alloc(dpsv1grd,  npart_new, zero, "dpsv1grd") !(npart)

  call alloc(xbob,    nblz_new, zero, "xbob") !(nblz)
  call alloc(ybob,    nblz_new, zero, "ybob") !(nblz)
  call alloc(xpbob,   nblz_new, zero, "xpbob") !(nblz)
  call alloc(ypbob,   nblz_new, zero, "ypbob") !(nblz)

  call alloc(xineff,  npart_new, zero, "xineff") !(npart)
  call alloc(yineff,  npart_new, zero, "yineff") !(npart)
  call alloc(xpineff, npart_new, zero, "xpineff") !(npart)
  call alloc(ypineff, npart_new, zero, "ypineff") !(npart)

  call alloc(mux,     nblz_new, zero, "mux") !(nblz)
  call alloc(muy,     nblz_new, zero, "muy") !(nblz)

  call alloc(counteddpop, npart_new, numeffdpop, 0, "counteddpop") !(npart,numeffdpop)
  call alloc(counted2d,   npart_new, numeff, numeffdpop, 0, "counted2d") !(npart,numeff,numeffdpop)

  call alloc(ename,    mNameLen, nblz_new, ' ', "ename") !(nblz_new)
  call alloc(nampl,    nblz_new, 0, "nampl") !(nblz_new)
  call alloc(sum_ax,   nblz_new, zero, "sum_ax") !(nblz_new)
  call alloc(sqsum_ax, nblz_new, zero, "sqsum_ax") !(nblz_new)
  call alloc(sum_ay,   nblz_new, zero, "sum_ay") !(nblz_new)
  call alloc(sqsum_ay, nblz_new, zero, "sqsum_ay") !(nblz_new)
  call alloc(sampl,    nblz_new, zero, "sampl") !(nblz_new)

  call alloc(secondary,            npart_new, 0, "secondary") !(npart_new)
  call alloc(tertiary,             npart_new, 0, "tertiary") !(npart_new)
  call alloc(other,                npart_new, 0, "other") !(npart_new)
  call alloc(scatterhit,           npart_new, 0, "scatterhit") !(npart_new)
  call alloc(part_hit_before_pos,  npart_new, 0, "part_hit_before_pos") !(npart_new)
  call alloc(part_hit_before_turn, npart_new, 0, "part_hit_before_turn") !(npart_new)
  call alloc(part_hit_pos,         npart_new, 0, "part_hit_pos") !(npart_new)
  call alloc(part_hit_turn,        npart_new, 0, "part_hit_turn") !(npart_new)
  call alloc(part_abs_pos,         npart_new, 0, "part_abs_pos") !(npart_new)
  call alloc(part_select,          npart_new, 0, "part_select") !(npart_new)
  call alloc(nabs_type,            npart_new, 0, "nabs_type") !(npart_new)

  call alloc(part_impact,    npart_new, zero, "part_impact") !(npart_new)
  call alloc(part_indiv,     npart_new, zero, "part_indiv") !(npart_new)
  call alloc(part_linteract, npart_new, zero, "part_linteract") !(npart_new)

  call alloc(counted_r, npart_new, numeff, 0, "counted_r") !(npart_new,numeff)
  call alloc(counted_x, npart_new, numeff, 0, "counted_x") !(npart_new,numeff)
  call alloc(counted_y, npart_new, numeff, 0, "counted_y") !(npart_new,numeff)

  ! Change the following block to npart
  call alloc(myx,  npart_new, zero, "myx")  !(npart_new)
  call alloc(myxp, npart_new, zero, "myxp") !(npart_new)
  call alloc(myy,  npart_new, zero, "myy")  !(npart_new)
  call alloc(myyp, npart_new, zero, "myyp") !(npart_new)
  call alloc(myp,  npart_new, zero, "myp")  !(npart_new)
  call alloc(mys,  npart_new, zero, "mys")  !(npart_new)

end subroutine collimation_expand_arrays

!>
!! collimate_init()
!! This routine is called once at the start of the simulation and
!! can be used to do any initial configuration and/or file loading.
!<
subroutine collimate_init()

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_settings
  use string_tools

  implicit none

#ifdef HDF5
  type(h5_dataField), allocatable :: fldDist0(:)
  integer                         :: fmtDist0, setDist0
#endif
  integer :: i,j

#ifdef HDF5
  if(h5_useForCOLL) call h5_initForCollimation
#endif

#ifdef G4COLLIMAT
! These should be configured in the scatter block when possible/enabled
  real(kind=fPrec) g4_ecut
  integer g4_physics
#endif

  call f_requestUnit('colltrack.out', outlun)
  open(unit=outlun, file='colltrack.out')

  if(st_quiet == 0) then
    write(lout,"(a)") '         -------------------------------'
    write(lout,"(a)")
    write(lout,"(a)") '          Program      C O L L T R A C K '
    write(lout,"(a)")
    write(lout,"(a)") '            R. Assmann           -    AB/ABP'
    write(lout,"(a)") '            C. Bracco            -    AB/ABP'
    write(lout,"(a)") '            V. Previtali         -    AB/ABP'
    write(lout,"(a)") '            S. Redaelli          -    AB/OP'
    write(lout,"(a)") '            G. Robert-Demolaize  -    BNL'
    write(lout,"(a)") '            A. Rossi             -    AB/ABP'
    write(lout,"(a)") '            T. Weiler            -    IEKP'
    write(lout,"(a)") '                 CERN 2001 - 2009'
    write(lout,"(a)")
    write(lout,"(a)") '         -------------------------------'
  else
    write(lout,"(a)") ""
    write(lout,"(a)") " INITIALISING COLLIMATION"
    write(lout,"(a)") ""
  end if

  write(outlun,*)
  write(outlun,*)
  write(outlun,*) '         -------------------------------'
  write(outlun,*)
  write(outlun,*) '         Program      C O L L T R A C K '
  write(outlun,*)
  write(outlun,*) '            R. Assmann       -    AB/ABP'
  write(outlun,*) '             C.Bracco        -    AB/ABP'
  write(outlun,*) '           V. Previtali      -    AB/ABP'
  write(outlun,*) '           S. Redaelli       -    AB/OP'
  write(outlun,*) '      G. Robert-Demolaize    -    BNL'
  write(outlun,*) '             A. Rossi        -    AB/ABP'
  write(outlun,*) '             T. Weiler       -    IEKP'
  write(outlun,*)
  write(outlun,*) '                 CERN 2001 - 2009'
  write(outlun,*)
  write(outlun,*) '         -------------------------------'
  write(outlun,*)
  write(outlun,*)

  if(st_quiet == 0) then
    write(lout,"(a)") '                     R. Assmann, F. Schmidt, CERN'
    write(lout,"(a)") '                           C. Bracco,        CERN'
    write(lout,"(a)") '                           V. Previtali,     CERN'
    write(lout,"(a)") '                           S. Redaelli,      CERN'
    write(lout,"(a)") '                       G. Robert-Demolaize,  BNL'
    write(lout,"(a)") '                           A. Rossi,         CERN'
    write(lout,"(a)") '                           T. Weiler         IEKP'

    write(lout,"(a)")
    write(lout,"(a)") 'Generating particle distribution at FIRST element!'
    write(lout,"(a)") 'Optical functions obtained from Sixtrack internal!'
    write(lout,"(a)") 'Emittance and energy obtained from Sixtrack input!'
    write(lout,"(a)")
    write(lout,"(a)")
  end if

  write(lout,"(a,e15.8)") 'COLL> Info: Betax0 [m]          = ', tbetax(1)
  write(lout,"(a,e15.8)") 'COLL> Info: Betay0 [m]          = ', tbetay(1)
  write(lout,"(a,e15.8)") 'COLL> Info: Alphax0             = ', talphax(1)
  write(lout,"(a,e15.8)") 'COLL> Info: Alphay0             = ', talphay(1)
  write(lout,"(a,e15.8)") 'COLL> Info: Orbitx0 [mm]        = ', torbx(1)
  write(lout,"(a,e15.8)") 'COLL> Info: Orbitxp0 [mrad]     = ', torbxp(1)
  write(lout,"(a,e15.8)") 'COLL> Info: Orbity0 [mm]        = ', torby(1)
  write(lout,"(a,e15.8)") 'COLL> Info: Orbitpy0 [mrad]     = ', torbyp(1)
  write(lout,"(a,e15.8)") 'COLL> Info: Emitx0_dist [um]    = ', remitx_dist
  write(lout,"(a,e15.8)") 'COLL> Info: Emity0_dist [um]    = ', remity_dist
  write(lout,"(a,e15.8)") 'COLL> Info: Emitx0_collgap [um] = ', remitx_collgap
  write(lout,"(a,e15.8)") 'COLL> Info: Emity0_collgap [um] = ', remity_collgap
  write(lout,"(a,e15.8)") 'COLL> Info: E0 [MeV]            = ', e0
  write(lout,"(a)")

  myemitx0_dist = remitx_dist*c1m6
  myemity0_dist = remity_dist*c1m6
  myemitx0_collgap = remitx_collgap*c1m6
  myemity0_collgap = remity_collgap*c1m6

  myalphax = talphax(1)
  myalphay = talphay(1)
  mybetax  = tbetax(1)
  mybetay  = tbetay(1)

!07-2006      myenom   = e0
!      MYENOM   = 1.001*E0
!
  if (myemitx0_dist.le.zero .or. myemity0_dist.le.zero .or. myemitx0_collgap.le.zero .or. myemity0_collgap.le.zero) then
    write(lout,"(a)") "COLL> ERROR Emittances not defined! check collimat block!"
    write(lout,"(a)") "COLL> ERROR Expected format of line 9 in collimat block:"
    write(lout,"(a)") "COLL> ERROR emitnx0_dist  emitny0_dist  emitnx0_collgap  emitny0_collgap"
    write(lout,"(a)") "COLL> ERROR All emittances should be normalized. "//&
      "first put emittance for distribtion generation, then for collimator position etc. units in [mm*mrad]."
    write(lout,"(a)") "COLL> ERROR EXAMPLE: 2.5 2.5 3.5 3.5"
    call prror(-1)
  end if

!++  Calculate the gammas
  mygammax = (one+myalphax**2)/mybetax
  mygammay = (one+myalphay**2)/mybetay

!++  Number of points and generate distribution
!
!GRD SEMI-AUTOMATIC INPUT
!      NLOOP=10
!      MYNEX=6.003
!      MYDEX=0.0015
!      MYNEY=6.003
!      MYDEY=0.0015
!      DO_COLL=1
!      NSIG_PRIM=5.
!      NSIG_SEC=6.
  rselect=64

  write(lout,"(a,i0)")    'COLL> Info: NLOOP               = ', nloop
  write(lout,"(a,i0)")    'COLL> Info: DO_THISDIS          = ', do_thisdis
  write(lout,"(a,e15.8)") 'COLL> Info: MYNEX               = ', mynex
  write(lout,"(a,e15.8)") 'COLL> Info: MYDEX               = ', mdex
  write(lout,"(a,e15.8)") 'COLL> Info: MYNEY               = ', myney
  write(lout,"(a,e15.8)") 'COLL> Info: MYDEY               = ', mdey
  write(lout,"(a,a)")     'COLL> Info: FILENAME_DIS        = ', trim(filename_dis)
  write(lout,"(a,e15.8)") 'COLL> Info: ENERROR             = ', enerror
  write(lout,"(a,e15.8)") 'COLL> Info: BUNCHLENGTH         = ', bunchlength
  write(lout,"(a,i0)")    'COLL> Info: RSELECT             = ', int(rselect)
  write(lout,"(a,l1)")    'COLL> Info: DO_COLL             = ', do_coll
  write(lout,"(a,l1)")    'COLL> Info: DO_NSIG             = ', cdb_doNSig
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCP3           = ', nsig_tcp3
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCSG3          = ', nsig_tcsg3
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCSM3          = ', nsig_tcsm3
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCLA3          = ', nsig_tcla3
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCP7           = ', nsig_tcp7
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCSG7          = ', nsig_tcsg7
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCSM7          = ', nsig_tcsm7
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCLA7          = ', nsig_tcla7
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCLP           = ', nsig_tclp
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCLI           = ', nsig_tcli
! write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTH           = ', nsig_tcth
! write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTV           = ', nsig_tctv
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTH1          = ', nsig_tcth1
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTV1          = ', nsig_tctv1
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTH2          = ', nsig_tcth2
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTV2          = ', nsig_tctv2
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTH5          = ', nsig_tcth5
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTV5          = ', nsig_tctv5
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTH8          = ', nsig_tcth8
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCTV8          = ', nsig_tctv8
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCDQ           = ', nsig_tcdq
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCSTCDQ        = ', nsig_tcstcdq
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TDI            = ', nsig_tdi
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCXRP          = ', nsig_tcxrp
  write(lout,"(a,e15.8)") 'COLL> Info: NSIG_TCRYP          = ', nsig_tcryo

  write(lout,"(a)")
  write(lout,"(a)")       'COLL> INPUT PARAMETERS FOR THE SLICING:'
  write(lout,"(a)")
  write(lout,"(a,i0)")    'COLL> Info: N_SLICES            = ',n_slices
  write(lout,"(a,e15.8)") 'COLL> Info: SMIN_SLICES         = ',smin_slices
  write(lout,"(a,e15.8)") 'COLL> Info: SMAX_SLICES         = ',smax_slices
  write(lout,"(a,e15.8)") 'COLL> Info: RECENTER1           = ',recenter1
  write(lout,"(a,e15.8)") 'COLL> Info: RECENTER2           = ',recenter2
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL> Info: FIT1_1              = ',fit1_1
  write(lout,"(a,e15.8)") 'COLL> Info: FIT1_2              = ',fit1_2
  write(lout,"(a,e15.8)") 'COLL> Info: FIT1_3              = ',fit1_3
  write(lout,"(a,e15.8)") 'COLL> Info: FIT1_4              = ',fit1_4
  write(lout,"(a,e15.8)") 'COLL> Info: FIT1_5              = ',fit1_5
  write(lout,"(a,e15.8)") 'COLL> Info: FIT1_6              = ',fit1_6
  write(lout,"(a,e15.8)") 'COLL> Info: SCALING1            = ',ssf1
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL> Info: FIT2_1              = ',fit2_1
  write(lout,"(a,e15.8)") 'COLL> Info: FIT2_2              = ',fit2_2
  write(lout,"(a,e15.8)") 'COLL> Info: FIT2_3              = ',fit2_3
  write(lout,"(a,e15.8)") 'COLL> Info: FIT2_4              = ',fit2_4
  write(lout,"(a,e15.8)") 'COLL> Info: FIT2_5              = ',fit2_5
  write(lout,"(a,e15.8)") 'COLL> Info: FIT2_6              = ',fit2_6
  write(lout,"(a,e15.8)") 'COLL> Info: SCALING2            = ',ssf2
  write(lout,"(a)")

!SEPT2005
!
! HERE WE CHECK IF THE NEW INPUT IS READ CORRECTLY
!
  write(lout,"(a,e15.8)") 'COLL> Info: EMITXN0_DIST        = ', emitnx0_dist
  write(lout,"(a,e15.8)") 'COLL> Info: EMITYN0_DIST        = ', emitny0_dist
  write(lout,"(a,e15.8)") 'COLL> Info: EMITXN0_COLLGAP     = ', emitnx0_collgap
  write(lout,"(a,e15.8)") 'COLL> Info: EMITYN0_COLLGAP     = ', emitny0_collgap
  write(lout,"(a)")
  write(lout,"(a,l1)")    'COLL> Info: DO_SELECT           = ', do_select
  write(lout,"(a,l1)")    'COLL> Info: DO_NOMINAL          = ', do_nominal
  write(lout,"(a,i0)")    'COLL> Info: RND_SEED            = ', rnd_seed
  write(lout,"(a,l1)")    'COLL> Info: DOWRITE_DIST        = ', dowrite_dist
  write(lout,"(a,a)")     'COLL> Info: NAME_SEL            = ', name_sel
  write(lout,"(a,l1)")    'COLL> Info: DO_ONESIDE          = ', do_oneside
  write(lout,"(a,l1)")    'COLL> Info: DOWRITE_IMPACT      = ', dowrite_impact
  write(lout,"(a,l1)")    'COLL> Info: DOWRITE_SECONDARY   = ', dowrite_secondary
  write(lout,"(a,l1)")    'COLL> Info: DOWRITE_AMPLITUDE   = ', dowrite_amplitude
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL> Info: XBEAT               = ', xbeat
  write(lout,"(a,e15.8)") 'COLL> Info: XBEATPHASE          = ', xbeatphase
  write(lout,"(a,e15.8)") 'COLL> Info: YBEAT               = ', ybeat
  write(lout,"(a,e15.8)") 'COLL> Info: YBEATPHASE          = ', ybeatphase
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL> Info: C_RMSTILT_PRIM      = ', c_rmstilt_prim
  write(lout,"(a,e15.8)") 'COLL> Info: C_RMSTILT_SEC       = ', c_rmstilt_sec
  write(lout,"(a,e15.8)") 'COLL> Info: C_SYSTILT_PRIM      = ', c_systilt_prim
  write(lout,"(a,e15.8)") 'COLL> Info: C_SYSTILT_SEC       = ', c_systilt_sec
  write(lout,"(a,e15.8)") 'COLL> Info: C_RMSOFFSET_PRIM    = ', c_rmsoffset_prim
  write(lout,"(a,e15.8)") 'COLL> Info: C_SYSOFFSET_PRIM    = ', c_sysoffset_prim
  write(lout,"(a,e15.8)") 'COLL> Info: C_RMSOFFSET_SEC     = ', c_rmsoffset_sec
  write(lout,"(a,e15.8)") 'COLL> Info: C_SYSOFFSET_SEC     = ', c_sysoffset_sec
  write(lout,"(a,i0)")    'COLL> Info: C_OFFSETTITLT_SEED  = ', c_offsettilt_seed
  write(lout,"(a,e15.8)") 'COLL> Info: C_RMSERROR_GAP      = ', c_rmserror_gap
  write(lout,"(a,l1)")    'COLL> Info: DO_MINGAP           = ', do_mingap
  write(lout,"(a)")
  write(lout,"(a,l1)")    'COLL> Info: RADIAL              = ', radial
  write(lout,"(a,e15.8)") 'COLL> Info: NR                  = ', nr
  write(lout,"(a,e15.8)") 'COLL> Info: NDR                 = ', ndr
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL: Info: DRIFTSX             = ', driftsx
  write(lout,"(a,e15.8)") 'COLL: Info: DRIFTSY             = ', driftsy
  write(lout,"(a,l1)")    'COLL: Info: CUT_INPUT           = ', cut_input
  write(lout,"(a,l1)")    'COLL: Info: SYSTILT_ANTISYMM    = ', systilt_antisymm
  write(lout,"(a)")
  write(lout,"(a,i0)")    'COLL> Info: IPENCIL             = ', ipencil
  write(lout,"(a,e15.8)") 'COLL> Info: PENCIL_OFFSET       = ', pencil_offset
  write(lout,"(a,e15.8)") 'COLL> Info: PENCIL_RMSX         = ', pencil_rmsx
  write(lout,"(a,e15.8)") 'COLL> Info: PENCIL_RMSY         = ', pencil_rmsy
  write(lout,"(a,i0)")    'COLL> Info: PENCIL_DISTR        = ', pencil_distr
  write(lout,"(a)")
  write(lout,"(a,a)")     'COLL> Info: COLL_DB             = ', coll_db
  write(lout,"(a,i0)")    'COLL> Info: IBEAM               = ', ibeam
  write(lout,"(a)")
  write(lout,"(a,l1)")    'COLL> Info: DOWRITETRACKS       = ', dowritetracks
  write(lout,"(a)")
  write(lout,"(a,l1)")    'COLL> Info: CERN                = ', cern
  write(lout,"(a)")
  write(lout,"(a,a)")     'COLL> Info: CASTORDIR           = ', castordir
  write(lout,"(a)")
  write(lout,"(a,i0)")    'COLL> Info: JOBNUMBER           = ', jobnumber
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL> Info: SIGSECUT2           = ', sigsecut2
  write(lout,"(a,e15.8)") 'COLL> Info: SIGSECUT3           = ', sigsecut3
  write(lout,"(a)")

  napx00 = napx

  write(lout,"(a,i0)")    'COLL> Info: NAPX                = ', napx
  write(lout,"(a,e15.8)") 'COLL> Info: Sigma_x0            = ', sqrt(mybetax*myemitx0_dist)
  write(lout,"(a,e15.8)") 'COLL> Info: Sigma_y0            = ', sqrt(mybetay*myemity0_dist)
  write(lout,"(a)")

! HERE WE SET THE MARKER FOR INITIALIZATION:
  firstrun = .true.

! ...and here is implemented colltrack's beam distribution:

!++  Initialize random number generator
  if (rnd_seed.eq.0) rnd_seed = mclock_liar()
  if (rnd_seed.lt.0) rnd_seed = abs(rnd_seed)
  rnd_lux = 3
  rnd_k1  = 0
  rnd_k2  = 0
  call rluxgo(rnd_lux, rnd_seed, rnd_k1, rnd_k2)
!  call recuin(rnd_seed, 0)
  write(outlun,*) 'INFO>  rnd_seed: ', rnd_seed

!Call distribution routines only if collimation block is in fort.3, otherwise
!the standard sixtrack would be prevented by the 'stop' command
  if(radial) then
    call makedis_radial(myalphax, myalphay, mybetax, mybetay, myemitx0_dist, myemity0_dist, &
                        myenom, nr, ndr, myx, myxp, myy, myyp, myp, mys)
  else
    select case(do_thisdis)
    case(0)
      continue
    case(1)
      call makedis(myalphax, myalphay, mybetax, mybetay, myemitx0_dist, myemity0_dist, &
                   myenom, mynex, mdex, myney, mdey, myx, myxp, myy, myyp, myp, mys)
    case(2)
      call makedis_st(myalphax, myalphay, mybetax, mybetay, myemitx0_dist, myemity0_dist, &
                      myenom, mynex, mdex, myney, mdey, myx, myxp, myy, myyp, myp, mys)
    case(3)
      call makedis_de(myalphax, myalphay, mybetax, mybetay, myemitx0_dist, myemity0_dist, &
                      myenom, mynex, mdex, myney, mdey,myx, myxp, myy, myyp, myp, mys,enerror,bunchlength)
    case(4)
      call readdis(filename_dis, myx, myxp, myy, myyp, myp, mys)
    case(5)
      call makedis_ga(myalphax, myalphay, mybetax, mybetay, myemitx0_dist, myemity0_dist, &
                      myenom, mynex, mdex, myney, mdey, myx, myxp, myy, myyp, myp, mys, enerror, bunchlength )
    case(6)
      call readdis_norm(filename_dis, myalphax, myalphay, mybetax, mybetay, &
                        myemitx0_dist, myemity0_dist, myenom, myx, myxp, myy, myyp, myp, mys, enerror, bunchlength)
    case default
      write(lout,"(a)") "COLL> ERROR Review your distribution parameters!"
      call prror(-1)
    end select
  end if

!++  Reset distribution for pencil beam
!
  if(ipencil.gt.0) then
    write(lout,"(a)") "COLL> WARNING Distributions reset to pencil beam!"
    write(outlun,*) 'WARN>  Distributions reset to pencil beam!'
    do j = 1, napx
      myx(j)  = zero
      myxp(j) = zero
      myy(j)  = zero
      myyp(j) = zero
    end do
  endif

  ! Optionally write the generated particle distribution
  if(dowrite_dist .and. do_thisdis /= 0) then
#ifdef HDF5
    if(h5_useForCOLL) then
      allocate(fldDist0(6))
      fldDist0(1)  = h5_dataField(name="X",  type=h5_typeReal)
      fldDist0(2)  = h5_dataField(name="XP", type=h5_typeReal)
      fldDist0(3)  = h5_dataField(name="Y",  type=h5_typeReal)
      fldDist0(4)  = h5_dataField(name="YP", type=h5_typeReal)
      fldDist0(5)  = h5_dataField(name="S",  type=h5_typeReal)
      fldDist0(6)  = h5_dataField(name="P",  type=h5_typeReal)
      call h5_createFormat("collDist0", fldDist0, fmtDist0)
      call h5_createDataSet("dist0", h5_collID, fmtDist0, setDist0, napx)
      call h5_prepareWrite(setDist0, napx)
      call h5_writeData(setDist0, 1, napx, myx(1:napx))
      call h5_writeData(setDist0, 2, napx, myxp(1:napx))
      call h5_writeData(setDist0, 3, napx, myy(1:napx))
      call h5_writeData(setDist0, 4, napx, myyp(1:napx))
      call h5_writeData(setDist0, 5, napx, mys(1:napx))
      call h5_writeData(setDist0, 6, napx, myp(1:napx))
      call h5_finaliseWrite(setDist0)
      deallocate(fldDist0)
    else
#endif
      call f_requestUnit('dist0.dat', dist0_unit)
      open(unit=dist0_unit,file='dist0.dat') !was 52
      do j = 1, napx
        write(dist0_unit,'(6(1X,E23.15))') myx(j), myxp(j), myy(j), myyp(j), mys(j), myp(j)
      end do
      close(dist0_unit)
#ifdef HDF5
    end if
#endif
  end if

!++  Initialize efficiency array
  do i=1,iu
    sum_ax(i)   = zero
    sqsum_ax(i) = zero
    sum_ay(i)   = zero
    sqsum_ay(i) = zero
    nampl(i)    = zero
    sampl(i)    = zero
  end do

  nspx = zero
  nspy = zero
  ax0  = myalphax
  bx0  = mybetax
  mux0 = mux(1)
  ay0  = myalphay
  by0  = mybetay
  muy0 = muy(1)
  iturn = 1
  ie    = 1
  n_tot_absorbed = 0

  call f_requestUnit('CollPositions.dat', CollPositions_unit)
  open(unit=CollPositions_unit, file='CollPositions.dat')

  ! Read collimator database
  call cdb_readCollDB(coll_db)

!Then do any implementation specific initial loading
#ifdef COLLIMATE_K2
  call collimate_init_k2
#endif

#ifdef MERLINSCATTER
  call collimate_init_merlin
#endif

#ifdef G4COLLIMAT
!! This function lives in the G4Interface.cpp file in the g4collimat folder
!! Accessed by linking libg4collimat.a
!! Set the energy cut at 70% - i.e. 30% energy loss
  g4_ecut = 0.7_fPrec

!! Select the physics engine to use
!! 0 = FTFP_BERT
!! 1 = QGSP_BERT
  g4_physics = 0

  call g4_collimation_init(e0, rnd_seed, g4_ecut, g4_physics)
#endif
  write (lout,"(a)") ""
  write (lout,"(a)") "COLL> Finished collimate initialisation"
  write (lout,"(a)") ""

end subroutine collimate_init

! ================================================================================================ !
!  Parse Input Line
! ================================================================================================ !
subroutine collimate_parseInputLine(inLine, iLine, iErr)

  use string_tools
  use mod_common, only : napx

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  real(kind=fPrec) nSigIn
  integer nSplit, famID
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "COLL> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(iLine)

  case(1)

    if(nSplit /= 1) then
      write(lout,"(a,i0)") "COLL> ERROR Expected 1 value on line 1, got ",nSplit
      iErr = .true.
      return
    end if

    if(nSplit > 0) call chr_cast(lnSPlit(1),do_coll,iErr)

  case(2)

    if(nSplit > 0) call chr_cast(lnSPlit(1),nloop,iErr)
    if(nSplit > 1) call chr_cast(lnSPlit(2),myenom,iErr)

    if(nloop /= 1) then
      write(lout,"(a,i0)") "COLL> ERROR Support for multiple samples is deprecated. nloop must be 1, got ",nloop
      iErr = .true.
      return
    end if

    if(napx*2 > npart) then
      write(lout,"(2(a,i0))") "COLL> ERROR Maximum number of particles is ", npart, ", got ",(napx*2)
      iErr = .true.
      return
   endif

  case(3)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), do_thisdis,  iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), mynex,       iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), mdex,        iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), myney,       iErr)
    if(nSplit > 4)  call chr_cast(lnSPlit(5), mdey,        iErr)
    if(nSplit > 5)  filename_dis = lnSPlit(6)
    if(nSplit > 6)  call chr_cast(lnSPlit(7), enerror,     iErr)
    if(nSplit > 7)  call chr_cast(lnSPlit(8), bunchlength, iErr)

  case(4)
    if(nSplit > 0)  call chr_cast(lnSplit(1), cdb_doNSig,  iErr)
    if(nSplit > 1)  call chr_cast(lnSplit(2), nsig_tcp3,   iErr)
    call cdb_getFamilyID("tcp3",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcp3)
    if(nSplit > 2)  call chr_cast(lnSplit(3), nsig_tcsg3,  iErr)
    call cdb_getFamilyID("tcsg3",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcsg3)
    if(nSplit > 3)  call chr_cast(lnSplit(4), nsig_tcsm3,  iErr)
    call cdb_getFamilyID("tcsm3",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcsm3)
    if(nSplit > 4)  call chr_cast(lnSplit(5), nsig_tcla3,  iErr)
    call cdb_getFamilyID("tcla3",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcla3)
    if(nSplit > 5)  call chr_cast(lnSplit(6), nsig_tcp7,   iErr)
    call cdb_getFamilyID("tcp7",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcp7)
    if(nSplit > 6)  call chr_cast(lnSplit(7), nsig_tcsg7,  iErr)
    call cdb_getFamilyID("tcsg7",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcsg7)
    if(nSplit > 7)  call chr_cast(lnSplit(8), nsig_tcsm7,  iErr)
    call cdb_getFamilyID("tcsm7",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcsm7)
    if(nSplit > 8)  call chr_cast(lnSplit(9), nsig_tcla7,  iErr)
    call cdb_getFamilyID("tcla7",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcla7)
    if(nSplit > 9)  call chr_cast(lnSplit(10),nsig_tclp,   iErr)
    call cdb_getFamilyID("tclp",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tclp)
    if(nSplit > 10) call chr_cast(lnSplit(11),nsig_tcli,   iErr)
    call cdb_getFamilyID("tcli",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcli)
    if(nSplit > 11) call chr_cast(lnSplit(12),nsig_tcdq,   iErr)
    call cdb_getFamilyID("tcdq",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcdq)
    if(nSplit > 12) call chr_cast(lnSplit(13),nsig_tcstcdq,iErr)
    call cdb_getFamilyID("tcstcdq",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcstcdq)
    if(nSplit > 13) call chr_cast(lnSplit(14),nsig_tdi,    iErr)
    call cdb_getFamilyID("tdi",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tdi)

  case(5)
    if(nSplit > 0)  call chr_cast(lnSplit(1), nsig_tcth1,iErr)
    call cdb_getFamilyID("tcth1",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcth1)
    if(nSplit > 1)  call chr_cast(lnSplit(2), nsig_tcth2,iErr)
    call cdb_getFamilyID("tcth2",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcth2)
    if(nSplit > 2)  call chr_cast(lnSplit(3), nsig_tcth5,iErr)
    call cdb_getFamilyID("tcth5",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcth5)
    if(nSplit > 3)  call chr_cast(lnSplit(4), nsig_tcth8,iErr)
    call cdb_getFamilyID("tcth8",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcth8)
    if(nSplit > 4)  call chr_cast(lnSplit(5), nsig_tctv1,iErr)
    call cdb_getFamilyID("tctv1",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tctv1)
    if(nSplit > 5)  call chr_cast(lnSplit(6), nsig_tctv2,iErr)
    call cdb_getFamilyID("tctv2",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tctv2)
    if(nSplit > 6)  call chr_cast(lnSplit(7), nsig_tctv5,iErr)
    call cdb_getFamilyID("tctv5",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tctv5)
    if(nSplit > 7)  call chr_cast(lnSplit(8), nsig_tctv8,iErr)
    call cdb_getFamilyID("tctv8",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tctv8)
    if(nSplit > 8)  call chr_cast(lnSplit(9), nsig_tcxrp,iErr)
    call cdb_getFamilyID("tcxrp",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcxrp)
    if(nSplit > 9)  call chr_cast(lnSplit(10),nsig_tcryo,iErr)
    call cdb_getFamilyID("tcryo",famID,.true.)
    call cdb_setFamilyNSig(famID,nsig_tcryo)

  case(6)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), n_slices,   iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), smin_slices,iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), smax_slices,iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), recenter1,  iErr)
    if(nSplit > 4)  call chr_cast(lnSPlit(5), recenter2,  iErr)

  case(7)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), fit1_1,iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), fit1_2,iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), fit1_3,iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), fit1_4,iErr)
    if(nSplit > 4)  call chr_cast(lnSPlit(5), fit1_5,iErr)
    if(nSplit > 5)  call chr_cast(lnSPlit(6), fit1_6,iErr)
    if(nSplit > 6)  call chr_cast(lnSPlit(7), ssf1,  iErr)

  case(8)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), fit2_1,iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), fit2_2,iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), fit2_3,iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), fit2_4,iErr)
    if(nSplit > 4)  call chr_cast(lnSPlit(5), fit2_5,iErr)
    if(nSplit > 5)  call chr_cast(lnSPlit(6), fit2_6,iErr)
    if(nSplit > 6)  call chr_cast(lnSPlit(7), ssf2,  iErr)

  case(9)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), emitnx0_dist,   iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), emitny0_dist,   iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), emitnx0_collgap,iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), emitny0_collgap,iErr)

  case(10)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), do_select,        iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), do_nominal,       iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), rnd_seed,         iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), dowrite_dist,     iErr)
    if(nSplit > 4)  name_sel = lnSPlit(5)
    if(nSplit > 5)  call chr_cast(lnSPlit(6), do_oneside,       iErr)
    if(nSplit > 6)  call chr_cast(lnSPlit(7), dowrite_impact,   iErr)
    if(nSplit > 7)  call chr_cast(lnSPlit(8), dowrite_secondary,iErr)
    if(nSplit > 8)  call chr_cast(lnSPlit(9), dowrite_amplitude,iErr)

  case(11)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), xbeat,     iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), xbeatphase,iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), ybeat,     iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), ybeatphase,iErr)

  case(12)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), c_rmstilt_prim,   iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), c_rmstilt_sec,    iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), c_systilt_prim,   iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), c_systilt_sec,    iErr)
    if(nSplit > 4)  call chr_cast(lnSPlit(5), c_rmsoffset_prim, iErr)
    if(nSplit > 5)  call chr_cast(lnSPlit(6), c_rmsoffset_sec,  iErr)
    if(nSplit > 6)  call chr_cast(lnSPlit(7), c_sysoffset_prim, iErr)
    if(nSplit > 7)  call chr_cast(lnSPlit(8), c_sysoffset_sec,  iErr)
    if(nSplit > 8)  call chr_cast(lnSPlit(9), c_offsettilt_seed,iErr)
    if(nSplit > 9)  call chr_cast(lnSPlit(10),c_rmserror_gap,   iErr)
    if(nSplit > 10) call chr_cast(lnSPlit(11),do_mingap,        iErr)

  case(13)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), radial,iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), nr,    iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), ndr,   iErr)

  case(14)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), driftsx,         iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), driftsy,         iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), cut_input,       iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), systilt_antisymm,iErr)

  case(15)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), ipencil,      iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), pencil_offset,iErr)
    if(nSplit > 2)  call chr_cast(lnSPlit(3), pencil_rmsx,  iErr)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), pencil_rmsy,  iErr)
    if(nSplit > 4)  call chr_cast(lnSPlit(5), pencil_distr, iErr)
#ifdef G4COLLIMAT
    if(ipencil > 0) then
      write(lout,"(a)") "COLL> ERROR Pencil distribution not supported with geant4"
      iErr = .true.
      return
    endif
#endif

  case(16)
    if(nSplit > 0)  coll_db =  lnSPlit(1)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), ibeam, iErr)

  case(17)
    if(nSplit > 0)  call chr_cast(lnSPlit(1), dowritetracks,iErr)
    if(nSplit > 1)  call chr_cast(lnSPlit(2), cern,         iErr)
    if(nSplit > 2)  castordir = lnSPlit(3)
    if(nSplit > 3)  call chr_cast(lnSPlit(4), jobnumber,    iErr)
    if(nSplit > 4)  call chr_cast(lnSPlit(5), sigsecut2,    iErr)
    if(nSplit > 5)  call chr_cast(lnSPlit(6), sigsecut3,    iErr)

  case default
    write(lout,"(a,i0,a)") "COLL> ERROR Unexpected line ",iLine," encountered."
    iErr = .true.

  end select

end subroutine collimate_parseInputLine

subroutine collimate_postInput(gammar)

  real(kind=fPrec), intent(in) :: gammar

  call collimation_expand_arrays(npart,nblz,nele)

  remitx_dist    = emitnx0_dist*gammar
  remity_dist    = emitny0_dist*gammar
  remitx_collgap = emitnx0_collgap*gammar
  remity_collgap = emitny0_collgap*gammar

end subroutine collimate_postInput

!>
!! collimate_start_sample()
!! This routine is called from trauthin before each sample
!! is injected into thin 6d
!<
subroutine collimate_start_sample(nsample)

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da

  implicit none

  integer j
  integer, intent(in) :: nsample

#ifdef HDF5
  type(h5_dataField), allocatable :: setFields(:)
  integer fmtHdf
#endif

  j = nsample

! HERE WE OPEN ALL THE NEEDED OUTPUT FILES

! Survival Output
#ifdef HDF5
  if(h5_useForCOLL) then
    allocate(setFields(2))
    setFields(1) = h5_dataField(name="TURN",  type=h5_typeInt)
    setFields(2) = h5_dataField(name="NSURV", type=h5_typeInt)
    call h5_createFormat("collSurvival", setFields, fmtHdf)
    call h5_createDataSet("survival", h5_collID, fmtHdf, coll_hdf5_survival, numl)
    deallocate(setFields)
  else
#endif
    call f_requestUnit('survival.dat', survival_unit)
    open(unit=survival_unit, file='survival.dat') ! RB, DM: 2014 bug fix !was 44
    write(survival_unit,*) '# 1=turn 2=n_particle'
#ifdef HDF5
  end if
#endif

  call f_requestUnit("collgaps.dat", collgaps_unit)
  open(unit=collgaps_unit, file="collgaps.dat")
  if(firstrun) write(collgaps_unit,"(a1,1x,a2,1x,a16,4(1x,a19),1x,a4,5(1x,a13),1x,a13)") &
    "#","ID","name            ","angle[rad]","betax[m]","betay[m]","halfgap[m]","mat.",  &
    "length[m]","sigx[m]","sigy[m]","tilt1[rad]","tilt2[rad]","nsig"

  call f_requestUnit('collimator-temp.db', collimator_temp_db_unit)
  open(unit=collimator_temp_db_unit, file='collimator-temp.db') !was 40
!

! TW06/08 added ouputfile for real collimator settings (incluing slicing, ...)
  call f_requestUnit('collsettings.dat', collsettings_unit)
  open(unit=collsettings_unit, file='collsettings.dat') !was 55

  if(firstrun) then
    write(collsettings_unit,*) '# name  slicenumber  halfgap[m]  gap_offset[m] tilt jaw1[rad]  tilt jaw2[rad] length[m] material'
    write(CollPositions_unit,*) '%Ind           Name   Pos[m]'
  end if

  if(dowrite_impact) then
    call f_requestUnit('impact.db', impact_unit)
    open(unit=impact_unit,file='impact.dat') !was 49
    write(impact_unit,*) '# 1=impact 2=divergence'
  endif


  if (dowritetracks) then
!GRD SPECIAL FILE FOR SECONDARY HALO
    if(cern) then
      smpl = '1'

      pfile(1:8) = 'tracks2.'
      pfile(9:9) = smpl
      pfile(10:13) = '.dat'

      call f_requestUnit(pfile(1:13), tracks2_unit)
      open(unit=tracks2_unit,file=pfile(1:13))

    else
      call f_requestUnit('tracks2.dat', tracks2_unit)
      open(unit=tracks2_unit,file='tracks2.dat') !was 38
    end if !end if (cern)

    if(firstrun) write(tracks2_unit,*) '# 1=name 2=turn 3=s 4=x 5=xp 6=y 7=yp 8=DE/E 9=type'

!AUGUST2006:write pencul sheet beam coordiantes to file ---- TW
    call f_requestUnit('pencilbeam_distr.dat', pencilbeam_distr_unit)
    open(unit=pencilbeam_distr_unit, file='pencilbeam_distr.dat') !was 9997
    if(firstrun) write(pencilbeam_distr_unit,*) 'x    xp    y      yp'
#ifdef HDF5
    if(h5_writeTracks2) call h5tr2_init
#endif
  end if !end if (dowritetracks) then

!GRD-SR,09-02-2006 => new series of output controlled by the 'dowrite_impact flag
  if(do_select) then
    call f_requestUnit('coll_ellipse.dat', coll_ellipse_unit)
    open(unit=coll_ellipse_unit, file='coll_ellipse.dat') !was 45
    if(firstrun) then
      write(coll_ellipse_unit,*) '#  1=name 2=x 3=y 4=xp 5=yp 6=E 7=s 8=turn 9=halo 10=nabs_type'
    end if
  end if

  if(dowrite_impact) then
#ifdef HDF5
    if(h5_useForCOLL .and. firstrun) then

      ! All Impacts and All Absorbtions
      allocate(setFields(3))
      setFields(1) = h5_dataField(name="ID",   type=h5_typeInt)
      setFields(2) = h5_dataField(name="TURN", type=h5_typeInt)
      setFields(3) = h5_dataField(name="S",    type=h5_typeReal)
      call h5_createFormat("collAllImpactAbsorb", setFields, fmtHdf)
      call h5_createDataSet("all_impacts",     h5_collID, fmtHdf, coll_hdf5_allImpacts)
      call h5_createDataSet("all_absorptions", h5_collID, fmtHdf, coll_hdf5_allAbsorb)
      deallocate(setFields)

      ! First Impacts
      allocate(setFields(14))
      setFields(1)  = h5_dataField(name="ID",     type=h5_typeInt)
      setFields(2)  = h5_dataField(name="TURN",   type=h5_typeInt)
      setFields(3)  = h5_dataField(name="ICOLL",  type=h5_typeInt)
      setFields(4)  = h5_dataField(name="NABS",   type=h5_typeInt)
      setFields(5)  = h5_dataField(name="S_IMP",  type=h5_typeReal)
      setFields(6)  = h5_dataField(name="S_OUT",  type=h5_typeReal)
      setFields(7)  = h5_dataField(name="X_IN",   type=h5_typeReal)
      setFields(8)  = h5_dataField(name="XP_IN",  type=h5_typeReal)
      setFields(9)  = h5_dataField(name="Y_IN",   type=h5_typeReal)
      setFields(10) = h5_dataField(name="YP_IN",  type=h5_typeReal)
      setFields(11) = h5_dataField(name="X_OUT",  type=h5_typeReal)
      setFields(12) = h5_dataField(name="XP_OUT", type=h5_typeReal)
      setFields(13) = h5_dataField(name="Y_OUT",  type=h5_typeReal)
      setFields(14) = h5_dataField(name="YP_OUT", type=h5_typeReal)
      call h5_createFormat("collFirstImpacts", setFields, fmtHdf)
      call h5_createDataSet("first_impacts", h5_collID, fmtHdf, coll_hdf5_fstImpacts)
      deallocate(setFields)

      ! Coll Scatter
      allocate(setFields(7))
      setFields(1) = h5_dataField(name="ID",    type=h5_typeInt)
      setFields(2) = h5_dataField(name="TURN",  type=h5_typeInt)
      setFields(3) = h5_dataField(name="ICOLL", type=h5_typeInt)
      setFields(4) = h5_dataField(name="NABS",  type=h5_typeInt)
      setFields(5) = h5_dataField(name="DP",    type=h5_typeReal)
      setFields(6) = h5_dataField(name="DX",    type=h5_typeReal)
      setFields(7) = h5_dataField(name="DY",    type=h5_typeReal)
      call h5_createFormat("collScatter", setFields, fmtHdf)
      call h5_createDataSet("coll_scatter", h5_collID, fmtHdf, coll_hdf5_collScatter)
      deallocate(setFields)

    else
#endif
      call f_requestUnit('all_impacts.dat', all_impacts_unit)
      call f_requestUnit('all_absorptions.dat', all_absorptions_unit)
      call f_requestUnit('Coll_Scatter.dat', coll_scatter_unit)
      call f_requestUnit('FirstImpacts.dat', FirstImpacts_unit)

      open(unit=all_impacts_unit, file='all_impacts.dat') !was 46
      open(unit=all_absorptions_unit, file='all_absorptions.dat') !was 47
      open(unit=coll_scatter_unit, file='Coll_Scatter.dat') !was 3998
      open(unit=FirstImpacts_unit, file='FirstImpacts.dat') !was 39

      if (firstrun) then
        write(all_impacts_unit,'(a)') '# 1=name 2=turn 3=s'
        write(all_absorptions_unit,'(a)') '# 1=name 2=turn 3=s'
        write(FirstImpacts_unit,"(a)") "# 1=name, 2=iturn, 3=icoll, 4=nabs, 5=s_imp[m], 6=s_out[m], "//&
          "7=x_in(b!)[m], 8=xp_in, 9=y_in, 10=yp_in, 11=x_out [m], 12=xp_out, 13=y_out, 14=yp_out"
        write(coll_scatter_unit,"(a)") "# 1=icoll, 2=iturn, 3=np, 4=nabs (1:Nuclear-Inelastic,2:Nuclear-Elastic,3:pp-Elastic, "//&
          "4:Single-Diffractive,5:Coulomb), 5=dp, 6=dx', 7=dy'"
      end if ! if (firstrun) then
#ifdef HDF5
    end if
#endif
    call f_requestUnit('FLUKA_impacts.dat', FLUKA_impacts_unit)
    call f_requestUnit('FLUKA_impacts_all.dat', FLUKA_impacts_all_unit)
    open(unit=FLUKA_impacts_unit, file='FLUKA_impacts.dat') !was 48
    open(unit=FLUKA_impacts_all_unit, file='FLUKA_impacts_all.dat') !was 4801
    if (firstrun) then
      write(FLUKA_impacts_unit,'(a)') '# 1=icoll 2=c_rotation 3=s 4=x 5=xp 6=y 7=yp 8=nabs 9=np 10=turn'
      write(FLUKA_impacts_all_unit,'(a)') '# 1=icoll 2=c_rotation 3=s 4=x 5=xp 6=y 7=yp 8=nabs 9=np 10=turn'
    end if ! if (firstrun) then
  end if ! if(dowrite_impact) then

  if(name_sel(1:3).eq.'COL') then
    call f_requestUnit('RHIClosses.dat', RHIClosses_unit)
    open(unit=RHIClosses_unit, file='RHIClosses.dat') !was 555
    if(firstrun) write(RHIClosses_unit,'(a)') '# 1=name 2=turn 3=s 4=x 5=xp 6=y 7=yp 8=dp/p 9=type'
  end if

  ! Copy new particles to tracking arrays. Also add the orbit offset at start of ring!
  if(do_thisdis /= 0) then
    xv1(1:napx00)   = c1e3 *  myx(1:napx00) + torbx(1)
    yv1(1:napx00)   = c1e3 * myxp(1:napx00) + torbxp(1)
    xv2(1:napx00)   = c1e3 *  myy(1:napx00) + torby(1)
    yv2(1:napx00)   = c1e3 * myyp(1:napx00) + torbyp(1)
    sigmv(1:napx00) = mys(1:napx00)
    ejv(1:napx00)   = myp(1:napx00)
  end if

  do i = 1, napx00
    ! FOR NOT FAST TRACKING ONLY
    ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
    rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
    dpsv(j)=(ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
    oidpsv(j)=one/(one+dpsv(j))
    moidpsv(j)=mtc(j)/(one+dpsv(j))
    omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
    dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)

    partID(i)=i
    parentID(i)=i

    do ieff =1, numeff
      counted_r(i,ieff) = 0
      counted_x(i,ieff) = 0
      counted_y(i,ieff) = 0

      do ieffdpop =1, numeffdpop
        counted2d(i,ieff,ieffdpop) = 0
      end do

    end do

    do ieffdpop =1, numeffdpop
      counteddpop(i,ieffdpop) = 0
    end do

  end do

!!!!!!!!!!!!!!!!!!!!!!START THIN6D CUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!++  Some initialization
  do i = 1, numeff
    rsig(i) = (real(i,fPrec)/two - half) + five                           !hr08
  end do

  dpopbins(1) = c1m4

  do i = 2, numeffdpop
    dpopbins(i) = real(i-1,fPrec)*4e-4_fPrec
  end do

  firstcoll = .true.

!GRD HERE WE NEED TO INITIALIZE SOME COLLIMATION PARAMETERS
  napx = napx00

  do j = 1, napx
    part_hit_pos(j)   = 0
    part_hit_turn(j)  = 0
    part_abs_pos(j)   = 0
    part_abs_turn(j)  = 0
    part_select(j)    = 1
    part_indiv(j)     = -c1m6
    part_linteract(j) = zero
    part_impact(j)    = 0
    tertiary(j)       = 0
    secondary(j)      = 0
    other(j)          = 0
    scatterhit(j)     = 0
    nabs_type(j)      = 0
    ipart(j)          = j
    flukaname(j)      = 0
  end do

!++  This we only do once, for the first call to this routine. Numbers
!++  are saved in memory to use exactly the same info for each sample.
!++  COMMON block to decide for first usage and to save coll info.
  if(firstrun) then
  !Reading of collimation database moved to subroutine collimate_init

#ifdef BEAMGAS
!YIL call beam gas initiation routine
  call beamGasInit(myenom)
#endif

  write(lout,"(a)") ""
  write(lout,"(a,i0)") "COLL> Number of collimators: ",cdb_nColl
  do icoll = 1, cdb_nColl
    write(lout,"(a,i5,a)") "COLL> Collimator ",icoll,": "//cdb_cNameUC(icoll)//" "//cdb_cName(icoll)
    coll_found(icoll) = .false.
  end do
  write(lout,"(a)") ""

!******write settings for alignment error in colltrack.out file
  write(outlun,*) ' '
  write(outlun,*) 'Alignment errors settings (tilt, offset,...)'
  write(outlun,*) ' '
  write(outlun,*) 'SETTING> c_rmstilt_prim   : ', c_rmstilt_prim
  write(outlun,*) 'SETTING> c_rmstilt_sec    : ', c_rmstilt_sec
  write(outlun,*) 'SETTING> c_systilt_prim   : ', c_systilt_prim
  write(outlun,*) 'SETTING> c_systilt_sec    : ', c_systilt_sec
  write(outlun,*) 'SETTING> c_rmsoffset_prim : ', c_rmsoffset_prim
  write(outlun,*) 'SETTING> c_rmsoffset_sec  : ', c_rmsoffset_sec
  write(outlun,*) 'SETTING> c_sysoffset_prim : ', c_sysoffset_prim
  write(outlun,*) 'SETTING> c_sysoffset_sec  : ', c_sysoffset_sec
  write(outlun,*) 'SETTING> c_offsettilt seed: ', c_offsettilt_seed
  write(outlun,*) 'SETTING> c_rmserror_gap   : ', c_rmserror_gap
  write(outlun,*) 'SETTING> do_mingap        : ', do_mingap
  write(outlun,*) ' '

!     TW - 01/2007
!     added offset and random_seed for tilt and offset
!*****intialize random generator with offset_seed
  c_offsettilt_seed = abs(c_offsettilt_seed)
  rnd_lux = 3
  rnd_k1  = 0
  rnd_k2  = 0
  call rluxgo(rnd_lux, c_offsettilt_seed, rnd_k1, rnd_k2)
!      write(outlun,*) 'INFO>  c_offsettilt seed: ', c_offsettilt_seed

! reset counter to assure starting at the same position in case of
! using rndm5 somewhere else in the code before
  zbv = rndm5(1)

!++  Generate random tilts (Gaussian distribution plus systematic)
!++  Do this only for the first call of this routine (first sample)
!++  Keep all collimator database info and errors in memeory (COMMON
!++  block) in order to re-use exactly the same information for every
!++  sample.
  if(c_rmstilt_prim.gt.zero .or. c_rmstilt_sec.gt.zero .or. c_systilt_prim.ne.zero .or. c_systilt_sec.ne.zero) then
    do icoll = 1, cdb_nColl
      if(cdb_cNameUC(icoll)(1:3).eq.'TCP') then
        c_rmstilt = c_rmstilt_prim
        c_systilt = c_systilt_prim
      else
        c_rmstilt = c_rmstilt_sec
        c_systilt = c_systilt_sec
      end if

      db_tilt(icoll,1) = c_systilt+c_rmstilt*myran_gauss(three)

      if(systilt_antisymm) then
        db_tilt(icoll,2) = -one*c_systilt+c_rmstilt*myran_gauss(three)
      else
        db_tilt(icoll,2) =      c_systilt+c_rmstilt*myran_gauss(three)
      end if

      write(outlun,*) 'INFO>  Collimator ', cdb_cNameUC(icoll), ' jaw 1 has tilt [rad]: ', db_tilt(icoll,1)
      write(outlun,*) 'INFO>  Collimator ', cdb_cNameUC(icoll), ' jaw 2 has tilt [rad]: ', db_tilt(icoll,2)
    end do
  end if

!++  Generate random offsets (Gaussian distribution plus systematic)
!++  Do this only for the first call of this routine (first sample)
!++  Keep all collimator database info and errors in memeory (COMMON
!++  block) in order to re-use exactly the same information for every
!++  sample and throughout a all run.
 if(c_sysoffset_prim.ne.zero .or. c_sysoffset_sec.ne.zero .or.c_rmsoffset_prim.gt.zero .or.c_rmsoffset_sec.gt.zero) then
   do icoll = 1, cdb_nColl

     if(cdb_cNameUC(icoll)(1:3).eq.'TCP') then
       cdb_cOffset(icoll) = c_sysoffset_prim + c_rmsoffset_prim*myran_gauss(three)
     else
       cdb_cOffset(icoll) = c_sysoffset_sec +  c_rmsoffset_sec*myran_gauss(three)
     end if

     write(outlun,*) 'INFO>  offset: ', cdb_cNameUC(icoll), cdb_cOffset(icoll)
   end do
 endif

!++  Generate random offsets (Gaussian distribution)
!++  Do this only for the first call of this routine (first sample)
!++  Keep all collimator database info and errors in memeory (COMMON
!++  block) in order to re-use exactly the same information for every
!++  sample and throughout a all run.
!         if (c_rmserror_gap.gt.0.) then
!            write(outlun,*) 'INFO> c_rmserror_gap = ',c_rmserror_gap
  do icoll = 1, cdb_nColl
    gap_rms_error(icoll) = c_rmserror_gap * myran_gauss(three)
    write(outlun,*) 'INFO>  gap_rms_error: ', cdb_cNameUC(icoll),gap_rms_error(icoll)
  end do

!---- creating a file with beta-functions at TCP/TCS
  call f_requestUnit('twisslike.out', twisslike_unit)
  open(unit=twisslike_unit, file='twisslike.out') !was 10000
  call f_requestUnit('sigmasettings.out', sigmasettings_unit)
  open(unit=sigmasettings_unit, file='sigmasettings.out') !was 10001
  mingap = 20

  do j=1,iu
! this transformation gives the right marker/name to the corresponding
! beta-dunctions or vice versa ;)
    if(ic(j).le.nblo) then
      do jb=1,mel(ic(j))
        myix=mtyp(ic(j),jb)
      end do
    else
      myix=ic(j)-nblo
    end if

! Using same code-block as below to evalute the collimator opening
! for each collimator, this is needed to get the smallest collimator gap
! in principal only looking for primary and secondary should be enough
! JULY 2008 added changes (V6.503) for names in TCTV -> TCTVA and TCTVB
! both namings before and after V6.503 can be used
    if ( bez(myix)(1:2).eq.'TC'.or. bez(myix)(1:2).eq.'tc'.or. bez(myix)(1:2).eq.'TD'.or. bez(myix)(1:2).eq.'td'&
 &  .or. bez(myix)(1:3).eq.'COL'.or. bez(myix)(1:3).eq.'col') then
      if(bez(myix)(1:3).eq.'TCP' .or. bez(myix)(1:3).eq.'tcp') then
        if(bez(myix)(7:9).eq.'3.B' .or. bez(myix)(7:9).eq.'3.b') then
          nsig = nsig_tcp3
        else
          nsig = nsig_tcp7
        endif
      else if(bez(myix)(1:4).eq.'TCSG' .or. bez(myix)(1:4).eq.'tcsg') then
        if(bez(myix)(8:10).eq.'3.B' .or. bez(myix)(8:10).eq.'3.b' .or. bez(myix)(9:11).eq.'3.B' .or. bez(myix)(9:11).eq.'3.b') then
          nsig = nsig_tcsg3
        else
          nsig = nsig_tcsg7
        endif
        if(bez(myix)(5:6).eq.'.4'.and.bez(myix)(8:9).eq.'6.') then
          nsig = nsig_tcstcdq
        endif
      else if(bez(myix)(1:4).eq.'TCSP' .or. bez(myix)(1:4).eq.'tcsp') then
        if(bez(myix)(9:11).eq.'6.B'.or. bez(myix)(9:11).eq.'6.b') then
          nsig = nsig_tcstcdq
        end if
      else if(bez(myix)(1:4).eq.'TCSM' .or. bez(myix)(1:4).eq.'tcsm') then
        if(bez(myix)(8:10).eq.'3.B' .or. bez(myix)(8:10).eq.'3.b' .or.bez(myix)(9:11).eq.'3.B' .or. bez(myix)(9:11).eq.'3.b') then
          nsig = nsig_tcsm3
        else
          nsig = nsig_tcsm7
        end if
      else if(bez(myix)(1:4).eq.'TCLA' .or. bez(myix)(1:4).eq.'tcla') then
        if(bez(myix)(9:11).eq.'7.B' .or. bez(myix)(9:11).eq.'7.b') then
          nsig = nsig_tcla7
        else
          nsig = nsig_tcla3
        end if
      else if(bez(myix)(1:4).eq.'TCDQ' .or. bez(myix)(1:4).eq.'tcdq') then
         nsig = nsig_tcdq
      ! YIL11: Checking only the IR value for TCT's..
      else if(bez(myix)(1:4).eq.'TCTH'.or.bez(myix)(1:4).eq.'tcth'.or.bez(myix)(1:5).eq.'TCTPH'.or.bez(myix)(1:5).eq.'tctph') then
        if(bez(myix)(8:8).eq.'1' .or. bez(myix)(9:9).eq.'1' ) then
          nsig = nsig_tcth1
        else if(bez(myix)(8:8).eq.'2' .or. bez(myix)(9:9).eq.'2' ) then
          nsig = nsig_tcth2
        else if(bez(myix)(8:8).eq.'5'.or. bez(myix)(9:9).eq.'5' ) then
          nsig = nsig_tcth5
        else if(bez(myix)(8:8).eq.'8' .or.  bez(myix)(9:9).eq.'8' ) then
          nsig = nsig_tcth8
        end if
      else if(bez(myix)(1:4).eq.'TCTV'.or.bez(myix)(1:4).eq.'tctv'.or.bez(myix)(1:5).eq.'TCTPV'.or.bez(myix)(1:5).eq.'tctpv') then
        if(bez(myix)(8:8).eq.'1' .or. bez(myix)(9:9).eq.'1' ) then
           nsig = nsig_tctv1
        else if(bez(myix)(8:8).eq.'2' .or. bez(myix)(9:9).eq.'2' ) then
           nsig = nsig_tctv2
        else if(bez(myix)(8:8).eq.'5' .or. bez(myix)(9:9).eq.'5' ) then
           nsig = nsig_tctv5
        else if(bez(myix)(8:8).eq.'8' .or. bez(myix)(9:9).eq.'8' ) then
           nsig = nsig_tctv8
        end if
      else if(bez(myix)(1:3).eq.'TDI' .or. bez(myix)(1:3).eq.'tdi') then
        nsig = nsig_tdi
      else if(bez(myix)(1:4).eq.'TCLP' .or. bez(myix)(1:4).eq.'tclp' .or.bez(myix)(1:4).eq.'TCL.' .or.bez(myix)(1:4).eq.'tcl.'.or. &
 &            bez(myix)(1:4).eq.'TCLX' .or. bez(myix)(1:4).eq.'tclx') then
        nsig = nsig_tclp
      else if(bez(myix)(1:4).eq.'TCLI' .or. bez(myix)(1:4).eq.'tcli') then
         nsig = nsig_tcli
      else if(bez(myix)(1:4).eq.'TCXR' .or. bez(myix)(1:4).eq.'tcxr') then
        nsig = nsig_tcxrp
      !     TW 04/2008 ---- start adding TCRYO
      else if(bez(myix)(1:5).eq.'TCRYO'.or.bez(myix)(1:5).eq.'tcryo'.or.bez(myix)(1:5).eq.'TCLD.'.or.bez(myix)(1:5).eq.'tcld.') then
        nsig = nsig_tcryo
      !     TW 04/2008 ---- end adding TCRYO
      else if(bez(myix)(1:3).eq.'COL' .or. bez(myix)(1:3).eq.'col') then
        if(bez(myix)(1:4).eq.'COLM'.or.bez(myix)(1:4).eq.'colm'.or.bez(myix)(1:5).eq.'COLH0'.or.bez(myix)(1:5).eq.'colh0') then
          nsig = nsig_tcth1
        else if(bez(myix)(1:5).eq.'COLV0' .or. bez(myix)(1:5).eq.'colv0') then
          nsig = nsig_tcth2
        else if(bez(myix)(1:5).eq.'COLH1' .or. bez(myix)(1:5).eq.'colh1') then
      !     JUNE2005   HERE WE USE NSIG_TCTH2 AS THE OPENING IN THE VERTICAL
      !     JUNE2005   PLANE FOR THE PRIMARY COLLIMATOR OF RHIC; NSIG_TCTH5 STANDS
      !     JUNE2005   FOR THE OPENING OF THE FIRST SECONDARY COLLIMATOR OF RHIC
          nsig = nsig_tcth5
        else if(bez(myix)(1:5).eq.'COLV1' .or. bez(myix)(1:5).eq.'colv1') then
          nsig = nsig_tcth8
        else if(bez(myix)(1:5).eq.'COLH2' .or. bez(myix)(1:5).eq.'colh2') then
          nsig = nsig_tctv1
        end if
!     JUNE2005   END OF DEDICATED TREATMENT OF RHIC OPENINGS
      else
        write(lout,"(a)") "COLL> WARNING Problem detected while writing twisslike.out' and 'sigmasettings.out':"
        write(lout,"(a)") "COLL>         Collimator name '"//trim(adjustl(bez(myix)))//"' was not recognized."//&
          " -> Setting nsig = 1000.0."
        nsig = c1e3
      end if

      do i = 1, cdb_nColl
! start searching minimum gap
        if((cdb_cNameUC(i)(1:mNameLen).eq.bez(myix)(1:mNameLen)).or. &
           (cdb_cName(i)(1:mNameLen).eq.bez(myix)(1:mNameLen))) then
          if( cdb_cLength(i) .gt. zero ) then
            nsig_err = nsig + gap_rms_error(i)

! jaw 1 on positive side x-axis
            gap_h1 = nsig_err - sin_mb(db_tilt(i,1))*cdb_cLength(i)/2
            gap_h2 = nsig_err + sin_mb(db_tilt(i,1))*cdb_cLength(i)/2

! jaw 2 on negative side of x-axis (see change of sign comapred
! to above code lines, alos have a look to setting of tilt angle)
            gap_h3 = nsig_err + sin_mb(db_tilt(i,2))*cdb_cLength(i)/2
            gap_h4 = nsig_err - sin_mb(db_tilt(i,2))*cdb_cLength(i)/2

! find minumum halfgap
! --- searching for smallest halfgap
!! ---scaling for beta beat needed?
!                        if (do_nominal) then
!                           bx_dist = cdb_cBx(icoll) * scale_bx / scale_bx0
!                           by_dist = cdb_cBy(icoll) * scale_by / scale_by0
!                        else
!                           bx_dist = tbetax(j) * scale_bx / scale_bx0
!                           by_dist = tbetay(j) * scale_by / scale_by0
!                        endif
            if (do_nominal) then
              bx_dist = cdb_cBx(icoll)
              by_dist = cdb_cBy(icoll)
            else
              bx_dist = tbetax(j)
              by_dist = tbetay(j)
            end if

            sig_offset = cdb_cOffset(i)/(sqrt(bx_dist**2 * cos_mb(cdb_cRotation(i))**2 + by_dist**2 * sin_mb(cdb_cRotation(i))**2 ))
            write(twisslike_unit,*) bez(myix),tbetax(j),tbetay(j), torbx(j),torby(j), nsig, gap_rms_error(i)
        write(sigmasettings_unit,*) bez(myix), gap_h1, gap_h2, gap_h3, gap_h4, sig_offset, cdb_cOffset(i), nsig, gap_rms_error(i)

            if((gap_h1 + sig_offset) .le. mingap) then
              mingap = gap_h1 + sig_offset
              coll_mingap_id = i
              coll_mingap1 = cdb_cNameUC(i)
              coll_mingap2 = cdb_cName(i)
            else if((gap_h2 + sig_offset) .le. mingap) then
              mingap = gap_h2 + sig_offset
              coll_mingap_id = i
              coll_mingap1 = cdb_cNameUC(i)
              coll_mingap2 = cdb_cName(i)
            else if((gap_h3 - sig_offset) .le. mingap) then
              mingap = gap_h3 - sig_offset
              coll_mingap_id = i
              coll_mingap1 = cdb_cNameUC(i)
              coll_mingap2 = cdb_cName(i)
            else if((gap_h4 - sig_offset) .le. mingap) then
              mingap = gap_h4 - sig_offset
              coll_mingap_id = i
              coll_mingap1 = cdb_cNameUC(i)
              coll_mingap2 = cdb_cName(i)
            end if
          end if
        end if
      end do !do i = 1, cdb_nColl

! could be done more elegant the above code to search the minimum gap
! and should also consider the jaw tilt
    end if
  end do !do j=1,iu

  write(twisslike_unit,*) coll_mingap_id, coll_mingap1, coll_mingap2,  mingap
  write(twisslike_unit,*) 'INFO> IPENCIL initial ', ipencil

! if pencil beam is used and on collimator with smallest gap the
! distribution should be generated, set ipencil to coll_mingap_id
  if (ipencil.gt.0 .and. do_mingap) then
    ipencil = coll_mingap_id
  end if

  write(twisslike_unit,*) 'INFO> IPENCIL new (if do_mingap) ', ipencil
  write(sigmasettings_unit,*) coll_mingap_id, coll_mingap1, coll_mingap2,  mingap

! if pencil beam is used and on collimator with smallest gap the
! distribution should be generated, set ipencil to coll_mingap_id
  write(sigmasettings_unit,*) 'INFO> IPENCIL new (if do_mingap) ',ipencil
  write(sigmasettings_unit,*) 'INFO> rnd_seed is (before reinit)',rnd_seed

  close(twisslike_unit)
  close(sigmasettings_unit)

!****** re-intialize random generator with rnd_seed
!       reinit with initial value used in first call
  rnd_lux = 3
  rnd_k1  = 0
  rnd_k2  = 0
  call rluxgo(rnd_lux, rnd_seed, rnd_k2, rnd_k2)
!  call recuin(rnd_seed, 0)
! TW - 01/2007

!GRD INITIALIZE LOCAL ADDITIVE PARAMETERS, I.E. THE ONE WE DON'T WANT
!GRD TO KEEP OVER EACH LOOP
  do j=1,napx
    tertiary(j)=0
    secondary(j)=0
    other(j)=0
    scatterhit(j)=0
    nabs_type(j) = 0
  end do

  do k = 1, numeff
    neff(k)  = zero
    neffx(k) = zero
    neffy(k) = zero

   do j = 1, numeffdpop
     neff2d(k,j) = zero
   end do
  end do

  do k = 1, numeffdpop
    neffdpop(k)  = zero
    npartdpop(k) = 0
  end do

  do j=1,max_ncoll
    cn_impact(j)   = 0
    cn_absorbed(j) = 0
    csum(j)   = zero
    csqsum(j) = zero
  end do

!++ End of first call stuff (end of first run)
  end if

!GRD NOW WE CAN BEGIN THE LOOPS
end subroutine collimate_start_sample

!>
!! collimate_start_collimator()
!! This routine is called each time we hit a collimator
!<
subroutine collimate_start_collimator(stracki)

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use numerical_constants, only : c5m4

  implicit none

  integer :: j
  real(kind=fPrec), intent(in) :: stracki

  if(bez(myix)(1:3).eq.'TCP' .or. bez(myix)(1:3).eq.'tcp') then
    if(bez(myix)(7:9).eq.'3.B' .or. bez(myix)(7:9).eq.'3.b') then
      nsig = nsig_tcp3
    else
      nsig = nsig_tcp7
    end if

  else if(bez(myix)(1:4).eq.'TCSG' .or.  bez(myix)(1:4).eq.'tcsg') then
    if(bez(myix)(8:10).eq.'3.B'.or.bez(myix)(8:10).eq.'3.b'.or.bez(myix)(9:11).eq.'3.B'.or.bez(myix)(9:11).eq.'3.b') then
      nsig = nsig_tcsg3
    else
      nsig = nsig_tcsg7
    end if
    if((bez(myix)(5:6).eq.'.4'.and.bez(myix)(8:9).eq.'6.')) then
      nsig = nsig_tcstcdq
    end if
  else if(bez(myix)(1:4).eq.'TCSP' .or. bez(myix)(1:4).eq.'tcsp') then
    if(bez(myix)(9:11).eq.'6.B'.or. bez(myix)(9:11).eq.'6.b') then
      nsig = nsig_tcstcdq
    end if
  else if(bez(myix)(1:4).eq.'TCSM' .or. bez(myix)(1:4).eq.'tcsm') then
    if(bez(myix)(8:10).eq.'3.B' .or. bez(myix)(8:10).eq.'3.b' .or. bez(myix)(9:11).eq.'3.B' .or. bez(myix)(9:11).eq.'3.b') then
      nsig = nsig_tcsm3
    else
      nsig = nsig_tcsm7
    end if
  else if(bez(myix)(1:4).eq.'TCLA' .or. bez(myix)(1:4).eq.'tcla') then
    if(bez(myix)(9:11).eq.'7.B' .or. bez(myix)(9:11).eq.'7.b') then
      nsig = nsig_tcla7
    else
      nsig = nsig_tcla3
    endif
  else if(bez(myix)(1:4).eq.'TCDQ' .or. bez(myix)(1:4).eq.'tcdq') then
    nsig = nsig_tcdq
! YIL11: Checking only the IR value for TCT's..
  else if(bez(myix)(1:4).eq.'TCTH' .or. bez(myix)(1:4).eq.'tcth' .or. bez(myix)(1:5).eq.'TCTPH' .or. bez(myix)(1:5).eq.'tctph') then
    if(bez(myix)(8:8).eq.'1' .or. bez(myix)(9:9).eq.'1' ) then
      nsig = nsig_tcth1
    else if(bez(myix)(8:8).eq.'2' .or. bez(myix)(9:9).eq.'2' ) then
      nsig = nsig_tcth2
    else if(bez(myix)(8:8).eq.'5'.or. bez(myix)(9:9).eq.'5' ) then
      nsig = nsig_tcth5
    else if(bez(myix)(8:8).eq.'8' .or. bez(myix)(9:9).eq.'8' ) then
      nsig = nsig_tcth8
    end if
  else if(bez(myix)(1:4).eq.'TCTV' .or.bez(myix)(1:4).eq.'tctv'.or.bez(myix)(1:5).eq.'TCTPV' .or.bez(myix)(1:5).eq.'tctpv' ) then
    if(bez(myix)(8:8).eq.'1' .or. bez(myix)(9:9).eq.'1' ) then
       nsig = nsig_tctv1
    else if(bez(myix)(8:8).eq.'2' .or. bez(myix)(9:9).eq.'2' ) then
       nsig = nsig_tctv2
    else if(bez(myix)(8:8).eq.'5' .or. bez(myix)(9:9).eq.'5' ) then
       nsig = nsig_tctv5
    else if(bez(myix)(8:8).eq.'8' .or. bez(myix)(9:9).eq.'8' ) then
       nsig = nsig_tctv8
    end if
  else if(bez(myix)(1:3).eq.'TDI' .or. bez(myix)(1:3).eq.'tdi') then
    nsig = nsig_tdi
  else if(bez(myix)(1:4).eq.'TCLP' .or.bez(myix)(1:4).eq.'tclp'.or.bez(myix)(1:4).eq.'TCL.'.or.bez(myix)(1:4).eq.'tcl.'.or. &
&         bez(myix)(1:4).eq.'TCLX' .or.bez(myix)(1:4).eq.'tclx') then
    nsig = nsig_tclp
  else if(bez(myix)(1:4).eq.'TCLI' .or. bez(myix)(1:4).eq.'tcli') then
    nsig = nsig_tcli
  else if(bez(myix)(1:4).eq.'TCXR' .or. bez(myix)(1:4).eq.'tcxr') then
    nsig = nsig_tcxrp
  else if(bez(myix)(1:5).eq.'TCRYO'.or.bez(myix)(1:5).eq.'tcryo'.or.bez(myix)(1:5).eq.'TCLD.' .or. bez(myix)(1:5).eq.'tcld.') then
    nsig = nsig_tcryo
  else if(bez(myix)(1:3).eq.'COL' .or. bez(myix)(1:3).eq.'col') then
    if(bez(myix)(1:4).eq.'COLM' .or. bez(myix)(1:4).eq.'colm' .or. bez(myix)(1:5).eq.'COLH0' .or. bez(myix)(1:5).eq.'colh0') then
      nsig = nsig_tcth1
    elseif(bez(myix)(1:5).eq.'COLV0' .or. bez(myix)(1:5).eq.'colv0') then
      nsig = nsig_tcth2
    else if(bez(myix)(1:5).eq.'COLH1' .or. bez(myix)(1:5).eq.'colh1') then
!     JUNE2005   HERE WE USE NSIG_TCTH2 AS THE OPENING IN THE VERTICAL
!     JUNE2005   PLANE FOR THE PRIMARY COLLIMATOR OF RHIC; NSIG_TCTH5 STANDS
!     JUNE2005   FOR THE OPENING OF THE FIRST SECONDARY COLLIMATOR OF RHIC
      nsig = nsig_tcth5
    else if(bez(myix)(1:5).eq.'COLV1' .or. bez(myix)(1:5).eq.'colv1') then
      nsig = nsig_tcth8
    else if(bez(myix)(1:5).eq.'COLH2' .or. bez(myix)(1:5).eq.'colh2') then
      nsig = nsig_tctv1
    end if
  else
    if(firstrun.and.iturn.eq.1) then
      write(lout,"(a)") "COLL> WARNING When setting opening for collimator '"//trim(adjustl(bez(myix)))//&
        "' from fort.3. Name not recognized. Setting nsig = 1000.0"
    end if
  nsig=c1e3
!JUNE2005   END OF DEDICATED TREATMENT OF RHIC OPENINGS
  end if

!++  Write trajectory for any selected particle
  c_length = zero

! SR, 23-11-2005: To avoid binary entries in 'amplitude.dat'
  if( firstrun ) then
    if(rselect.gt.0 .and. rselect.lt.65) then
      do j = 1, napx
        xj  = (xv1(j)-torbx(ie))/c1e3
        xpj = (yv1(j)-torbxp(ie))/c1e3
        yj  = (xv2(j)-torby(ie))/c1e3
        ypj = (yv2(j)-torbyp(ie))/c1e3
        pj  = ejv(j)/c1e3

        if(iturn.eq.1.and.j.eq.1) then
          sum_ax(ie)=zero
          sum_ay(ie)=zero
        end if

!-- DRIFT PART
        if(stracki.eq.0.) then
          if(iexact.eq.0) then
            xj  = xj + half*c_length*xpj
            yj  = yj + half*c_length*ypj
          else
            zpj = sqrt(one-xpj**2-ypj**2)
            xj  = xj + half*c_length*(xpj/zpj)
            yj  = yj + half*c_length*(ypj/zpj)
          end if
        end if

        gammax = (one + talphax(ie)**2)/tbetax(ie)
        gammay = (one + talphay(ie)**2)/tbetay(ie)

        if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
          nspx = sqrt(abs( gammax*(xj)**2 + two*talphax(ie)*xj*xpj +tbetax(ie)*xpj**2 )/myemitx0_collgap)
          nspy = sqrt(abs( gammay*(yj)**2 + two*talphay(ie)*yj*ypj +tbetay(ie)*ypj**2 )/myemity0_collgap)
          sum_ax(ie)   = sum_ax(ie) + nspx
          sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
          sum_ay(ie)   = sum_ay(ie) + nspy
          sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
          nampl(ie)    = nampl(ie) + 1
        else
          nspx = zero
          nspy = zero
        end if

          sampl(ie)    = totals
          ename(ie)    = bez(myix)(1:mNameLen)
      end do !do j = 1, napx
    end if !if(rselect.gt.0 .and. rselect.lt.65) then
  end if !if( firstrun ) then

!GRD HERE WE LOOK FOR ADEQUATE DATABASE INFORMATION
  found = .false.

!     SR, 01-09-2005: to set found = .TRUE., add the condition L>0!!
  do j = 1, cdb_nColl
    if((cdb_cNameUC(j)(1:mNameLen).eq.bez(myix)(1:mNameLen)) .or. &
       (cdb_cName(j)(1:mNameLen).eq.bez(myix)(1:mNameLen))) then
      if( cdb_cLength(j) .gt. zero ) then
        found = .true.
        icoll = j
        if(firstrun) then
          coll_found(j) = .TRUE.
          write(CollPositions_unit,*) j, cdb_cNameUC(j), totals
        end if
      end if
    end if
  end do

  if(.not. found .and. firstrun .and. iturn.eq.1) then
    write(lout,"(a)") "COLL> WARNING Collimator not found in colldb: '"//trim(bez(myix))//"'"
  end if

end subroutine collimate_start_collimator

!>
!! collimate_do_collimator()
!! This routine is calls the actual scattering functions
!<
subroutine collimate_do_collimator(stracki)

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use numerical_constants, only : c5m4

  implicit none

  integer :: j

  real(kind=fPrec), intent(in) :: stracki

#ifdef G4COLLIMAT
  integer g4_lostc
  integer :: part_hit_flag = 0
  integer :: part_abs_flag = 0
  real(kind=fPrec) x_tmp,y_tmp,xp_tmp,yp_tmp
#endif

!-----------------------------------------------------------------------
!GRD NEW COLLIMATION PARAMETERS
!-----------------------------------------------------------------------
!++  Get the aperture from the beta functions and emittance
!++  A simple estimate of beta beating can be included that
!++  has twice the betatron phase advance
  if(.not. cdb_doNSig) nsig = cdb_cNSig(icoll)

  scale_bx = (one + xbeat*sin_mb(four*pi*mux(ie)+xbeatphase) )
  scale_by = (one + ybeat*sin_mb(four*pi*muy(ie)+ybeatphase) )

  if(firstcoll) then
    scale_bx0 = scale_bx
    scale_by0 = scale_by
    firstcoll = .false.
  end if
!-------------------------------------------------------------------
!++  Assign nominal OR design beta functions for later
  if(do_nominal) then
    bx_dist = cdb_cBx(icoll) * scale_bx / scale_bx0
    by_dist = cdb_cBy(icoll) * scale_by / scale_by0
  else
    bx_dist = tbetax(ie) * scale_bx / scale_bx0
    by_dist = tbetay(ie) * scale_by / scale_by0
  end if

!++  Write beam ellipse at selected collimator
  if (((cdb_cNameUC(icoll).eq.name_sel(1:mNameLen)) .or. (cdb_cName(icoll).eq.name_sel(1:mNameLen))) .and. do_select) then
    do j = 1, napx
      write(coll_ellipse_unit,'(1X,I8,6(1X,E15.7),3(1X,I4,1X,I4))') ipart(j),xv1(j), xv2(j), yv1(j), yv2(j), &
     &        ejv(j), mys(j),iturn,secondary(j)+tertiary(j)+other(j)+scatterhit(j),nabs_type(j)
    end do
  end if

!-------------------------------------------------------------------
!++  Output to temporary database and screen
  if(iturn.eq.1.and.firstrun) then
    write(collimator_temp_db_unit,*) '# '
    write(collimator_temp_db_unit,*) cdb_cNameUC(icoll)!(1:11)
    write(collimator_temp_db_unit,*) cdb_cMaterial(icoll)
    write(collimator_temp_db_unit,*) cdb_cLength(icoll)
    write(collimator_temp_db_unit,*) cdb_cRotation(icoll)
    write(collimator_temp_db_unit,*) cdb_cOffset(icoll)
    write(collimator_temp_db_unit,*) tbetax(ie)
    write(collimator_temp_db_unit,*) tbetay(ie)

    write(outlun,*) ' '
    write(outlun,*)   'Collimator information: '
    write(outlun,*) ' '
    write(outlun,*) 'Name:                ', cdb_cNameUC(icoll)!(1:11)
    write(outlun,*) 'Material:            ', cdb_cMaterial(icoll)
    write(outlun,*) 'Length [m]:          ', cdb_cLength(icoll)
    write(outlun,*) 'Rotation [rad]:      ', cdb_cRotation(icoll)
    write(outlun,*) 'Offset [m]:          ', cdb_cOffset(icoll)
    write(outlun,*) 'Design beta x [m]:   ', cdb_cBx(icoll)
    write(outlun,*) 'Design beta y [m]:   ', cdb_cBy(icoll)
    write(outlun,*) 'Optics beta x [m]:   ', tbetax(ie)
    write(outlun,*) 'Optics beta y [m]:   ', tbetay(ie)
  end if

!-------------------------------------------------------------------
!++  Calculate aperture of collimator
!JUNE2005   HERE ONE HAS TO HAVE PARTICULAR TREATMENT OF THE OPENING OF
!JUNE2005   THE PRIMARY COLLIMATOR OF RHIC
  if(cdb_cNameUC(icoll)(1:4).ne.'COLM') then
    nsig = nsig + gap_rms_error(icoll)
    xmax = nsig*sqrt(bx_dist*myemitx0_collgap)
    ymax = nsig*sqrt(by_dist*myemity0_collgap)
    xmax_pencil = (nsig+pencil_offset)*sqrt(bx_dist*myemitx0_collgap)
    ymax_pencil = (nsig+pencil_offset)*sqrt(by_dist*myemity0_collgap)
    xmax_nom   = cdb_cNSig(icoll)*sqrt(cdb_cBx(icoll)*myemitx0_collgap)
    ymax_nom   = cdb_cNSig(icoll)*sqrt(cdb_cBy(icoll)*myemity0_collgap)
    c_rotation = cdb_cRotation(icoll)
    c_length   = cdb_cLength(icoll)
    c_material = cdb_cMaterial(icoll)
    c_offset   = cdb_cOffset(icoll)
    c_tilt(1)  = db_tilt(icoll,1)
    c_tilt(2)  = db_tilt(icoll,2)

    calc_aperture   = sqrt( xmax**2 * cos_mb(c_rotation)**2 + ymax**2 * sin_mb(c_rotation)**2 )
    nom_aperture    = sqrt( xmax_nom**2 * cos_mb(c_rotation)**2 + ymax_nom**2 * sin_mb(c_rotation)**2 )
    pencil_aperture = sqrt( xmax_pencil**2 * cos_mb(c_rotation)**2+ ymax_pencil**2 * sin_mb(c_rotation)**2 )

!++  Get x and y offsets at collimator center point
    x_pencil(icoll) = xmax_pencil * (cos_mb(c_rotation))
    y_pencil(icoll) = ymax_pencil * (sin_mb(c_rotation))

!++  Get corresponding beam angles (uses xp_max)
    xp_pencil(icoll) = -one * sqrt(myemitx0_collgap/tbetax(ie))*talphax(ie)* xmax / sqrt(myemitx0_collgap*tbetax(ie))
    yp_pencil(icoll) = -one * sqrt(myemity0_collgap/tbetay(ie))*talphay(ie)* ymax / sqrt(myemity0_collgap*tbetay(ie))
    xp_pencil0 = xp_pencil(icoll)
    yp_pencil0 = yp_pencil(icoll)

    pencil_dx(icoll) = sqrt(xmax_pencil**2 * cos_mb(c_rotation)**2 + ymax_pencil**2 * sin_mb(c_rotation)**2)-calc_aperture

!++ TW -- tilt for of jaw for pencil beam
!++ as in Ralphs orig routine, but not in collimate subroutine itself
!            nprim = 3
!            if ( (icoll.eq.ipencil) &
!     &           icoll.le.nprim .and. (j.ge.(icoll-1)*nev/nprim)        &
!     &           .and. (j.le.(icoll)*nev/nprim))) then
! this is done for every bunch (64 particle bucket)
! important: Sixtrack calculates in "mm" and collimate2 in "m"
! therefore 1E-3 is used to

! RB: added condition that pencil_distr.ne.3 in order to do the tilt
    if((icoll.eq.ipencil).and.(iturn.eq.1).and. (pencil_distr.ne.3)) then
!!               write(*,*) " ************************************** "
!!               write(*,*) " * INFO> seting tilt for pencil beam  * "
!!               write(*,*) " ************************************** "

!! respects if the tilt symmetric or not, for systilt_antiymm c_tilt is
!! -systilt + rmstilt otherwise +systilt + rmstilt
!!               if (systilt_antisymm) then
!! to align the jaw/pencil to the beam always use the minus regardless which
!! orientation of the jaws was used (symmetric/antisymmetric)
      c_tilt(1) = c_tilt(1) +    (xp_pencil0*cos_mb(c_rotation) + sin_mb(c_rotation)*yp_pencil0)
      c_tilt(2) = c_tilt(2) -one*(xp_pencil0*cos_mb(c_rotation) + sin_mb(c_rotation)*yp_pencil0)
      write(lout,*) "INFO> Changed tilt1  ICOLL  to  ANGLE  ", icoll, c_tilt(1)
      write(lout,*) "INFO> Changed tilt2  ICOLL  to  ANGLE  ", icoll, c_tilt(2)
    end if
!++ TW -- tilt angle changed (added to genetated on if spec. in fort.3)

!JUNE2005   HERE IS THE SPECIAL TREATMENT...
  else if(cdb_cNameUC(icoll)(1:4).eq.'COLM') then
    xmax = nsig_tcth1*sqrt(bx_dist*myemitx0_collgap)
    ymax = nsig_tcth2*sqrt(by_dist*myemity0_collgap)

    c_rotation = cdb_cRotation(icoll)
    c_length   = cdb_cLength(icoll)
    c_material = cdb_cMaterial(icoll)
    c_offset   = cdb_cOffset(icoll)
    c_tilt(1)  = db_tilt(icoll,1)
    c_tilt(2)  = db_tilt(icoll,2)
    calc_aperture = xmax
    nom_aperture = ymax
  end if

!-------------------------------------------------------------------
!++  Further output
  if(firstrun) then
    if(iturn.eq.1) then
      write(outlun,*) xp_pencil(icoll), yp_pencil(icoll), pencil_dx(icoll)
      write(outlun,'(a,i4)') 'Collimator number:   ', icoll
      write(outlun,*) 'Beam size x [m]:     ', sqrt(tbetax(ie)*myemitx0_collgap), "(from collgap emittance)"
      write(outlun,*) 'Beam size y [m]:     ', sqrt(tbetay(ie)*myemity0_collgap), "(from collgap emittance)"
      write(outlun,*) 'Divergence x [urad]:     ', c1e6*xp_pencil(icoll)
      write(outlun,*) 'Divergence y [urad]:     ', c1e6*yp_pencil(icoll)
      write(outlun,*) 'Aperture (nom) [m]:  ', nom_aperture
      write(outlun,*) 'Aperture (cal) [m]:  ', calc_aperture
      write(outlun,*) 'Collimator halfgap [sigma]:  ', nsig
      write(outlun,*) 'RMS error on halfgap [sigma]:  ', gap_rms_error(icoll)
      write(outlun,*) ' '

      write(collgaps_unit,"(i4,1x,a16,4(1x,e19.10),1x,a4,5(1x,e13.5),1x,f13.6)") &
        icoll,cdb_cName(icoll)(1:16),cdb_cRotation(icoll),tbetax(ie),tbetay(ie),calc_aperture, &
        cdb_cMaterial(icoll),cdb_cLength(icoll),sqrt(tbetax(ie)*myemitx0_collgap), &
        sqrt(tbetay(ie)*myemity0_collgap),db_tilt(icoll,1),db_tilt(icoll,2),nsig

! coll settings file
      if(n_slices.le.1) then
        write(collsettings_unit,'(a,1x,i10,5(1x,e13.5),1x,a)')          &
     &cdb_cNameUC(icoll)(1:12),                                            &
     &n_slices,calc_aperture,                                           &
     &cdb_cOffset(icoll),                                                 &
     &db_tilt(icoll,1),                                                 &
     &db_tilt(icoll,2),                                                 &
     &cdb_cLength(icoll),                                                 &
     &cdb_cMaterial(icoll)
      end if !if(n_slices.le.1) then
    end if !if(iturn.eq.1) then
  end if !if(firstrun) then

!++  Assign aperture which we define as the FULL width (factor 2)!!!
!JUNE2005 AGAIN, SOME SPECIFIC STUFF FOR RHIC
  if(cdb_cNameUC(icoll)(1:4).eq.'COLM') then
    c_aperture = two*calc_aperture
    nom_aperture = two*nom_aperture
  else if(cdb_cNameUC(icoll)(1:4).ne.'COLM') then
    c_aperture = two*calc_aperture
  end if

  c_aperture = two*calc_aperture
!          IF(IPENCIL.GT.zero) THEN
!          C_APERTURE = 2.*pencil_aperture

  if(firstrun.and.iturn.eq.1.and.icoll.eq.7) then
    call f_requestUnit('distsec', distsec_unit)
    open(unit=distsec_unit,file='distsec') !was 99
    do j=1,napx
      write(distsec_unit,'(4(1X,E15.7))') xv1(j),yv1(j),xv2(j),yv2(j)
    end do
    close(distsec_unit)
  end if

! RB: addition matched halo sampled directly on the TCP using pencil beam flag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((iturn.eq.1).and.(ipencil.eq.icoll).and.(pencil_distr.eq.3)) then

!     create distribution where the normalized distance between jaw and beam is the smallest
!     - this is where particles will first impact:
!     without imperfections, it is:
!              -- at the face of the collimator for the case of beta'<0 (POSITIVE alpha - beam converging) and
!              -- at the exit of the collimator for the case of beta'>0 (NEGATIVE alpha beam diverging)

!     with imperfections: include errors on gap, tilt and offset. We have to calculate the normalized distance
!     to each corner separately!

!     First: calculate optical parameters at start and end of collimator (half a collimator length upstream and
!     downstream of present s-position)
!     Assuming a purely vertical or horizontal halo - need to add more conditions for other cases!

!     Using standard twiss transfer matrix for a drift : ( new_halo_model_checks.nb )
!     at start of collimator:
    ldrift = -c_length / two !Assign the drift length over which the optics functions are propagated
    betax1 = tbetax(ie) - two*ldrift*talphax(ie) + (ldrift**2 * (one+talphax(ie)**2))/tbetax(ie)
    betay1 = tbetay(ie) - two*ldrift*talphay(ie) + (ldrift**2 * (one+talphay(ie)**2))/tbetay(ie)

    alphax1 = talphax(ie) - (ldrift*(1+talphax(ie)**2))/tbetax(ie)
    alphay1 = talphay(ie) - (ldrift*(1+talphay(ie)**2))/tbetay(ie)

!   at end of collimator:
    ldrift = c_length / two
    betax2 = tbetax(ie) - two*ldrift*talphax(ie) + (ldrift**2 * (one+talphax(ie)**2))/tbetax(ie)
    betay2 = tbetay(ie) - two*ldrift*talphay(ie) + (ldrift**2 * (one+talphay(ie)**2))/tbetay(ie)

    alphax2 = talphax(ie) - (ldrift*(1+talphax(ie)**2))/tbetax(ie)
    alphay2 = talphay(ie) - (ldrift*(1+talphay(ie)**2))/tbetay(ie)

!   calculate beam size at start and end of collimator. account for collimation plane
    if((mynex.gt.0).and.(myney.eq.zero)) then  ! horizontal halo
      beamsize1 = sqrt(betax1 * myemitx0_collgap)
      beamsize2 = sqrt(betax2 * myemitx0_collgap)
    else if((mynex.eq.0).and.(myney.gt.zero)) then   ! vertical halo
      beamsize1 = sqrt(betay1 * myemity0_collgap)
      beamsize2 = sqrt(betay2 * myemity0_collgap)
    else
      write(lout,*) "attempting to use a halo not purely in the horizontal or vertical plane with pencil_dist=3 - abort."
      call prror(-1)
    end if

!   calculate offset from tilt of positive and negative jaws, at start and end
!   remember: tilt angle is defined such that one corner stays at nominal position, the other corner is more open

!   jaw in positive x (or y):
    if(c_tilt(1).ge.0) then
      tiltOffsPos1 = zero
      tiltOffsPos2 = abs(sin_mb(c_tilt(1))) * c_length
    else
      tiltOffsPos1 = abs(sin_mb(c_tilt(1))) * c_length
      tiltOffsPos2 = zero
    end if

!   jaw in negative x (or y):
    if(c_tilt(2).ge.0) then
      tiltOffsNeg1 = abs(sin_mb(c_tilt(2))) * c_length
      tiltOffsNeg2 = zero
    else
      tiltOffsNeg1 = zero
      tiltOffsNeg2 = abs(sin_mb(c_tilt(2))) * c_length
    end if

!   calculate half distance from jaws to beam center (in units of beam sigma) at the beginning of the collimator,
!     positive and neg jaws.
    Nap1pos=((c_aperture/two + c_offset) + tiltOffsPos1)/beamsize1
    Nap2pos=((c_aperture/two + c_offset) + tiltOffsPos2)/beamsize2
    Nap1neg=((c_aperture/two - c_offset) + tiltOffsNeg1)/beamsize1
    Nap2neg=((c_aperture/two - c_offset) + tiltOffsNeg2)/beamsize2

!   Minimum normalized distance from jaw to beam center - this is the n_sigma at which the halo should be generated
    minAmpl = min(Nap1pos,Nap2pos,Nap1neg,Nap2neg)

!   Assign amplitudes in x and y for the halo generation function
    if((mynex.gt.0).and.(myney.eq.zero)) then ! horizontal halo
       mynex2 = minAmpl
    else if((mynex.eq.0).and.(myney.gt.zero)) then ! vertical halo
       myney2 = minAmpl
    end if               ! other cases taken care of above - in these cases, program has already stopped

!   assign optics parameters to use for the generation of the starting halo - at start or end of collimator
    if((minAmpl.eq.Nap1pos).or.(minAmpl.eq.Nap1neg)) then ! min normalized distance occurs at start of collimator
      mybetax=betax1
      mybetay=betay1
      myalphax=alphax1
      myalphay=alphay1
      ldrift = -c_length / two
    else               ! min normalized distance occurs at end of collimator
      mybetax=betax2
      mybetay=betay2
      myalphax=alphax2
      myalphay=alphay2
      ldrift = c_length / two
    end if

!   create new pencil beam distribution with spread at start or end of collimator at the minAmpl
!   note: if imperfections are active, equal amounts of particles are still generated on the two jaws.
!   but it might be then that only one jaw is hit on the first turn, thus only by half of the particles
!   the particle generated on the other side will then hit the same jaw several turns later, possibly smearing the impact parameter
!   This could possibly be improved in the future.
    call makedis_coll(myalphax,myalphay, mybetax, mybetay, myemitx0_collgap, myemity0_collgap, &
 &                    myenom, mynex2, mdex, myney2, mdey, myx, myxp, myy, myyp, myp, mys)

    do j = 1, napx
      xv1(j)  = c1e3*myx(j)  + torbx(ie)
      yv1(j)  = c1e3*myxp(j) + torbxp(ie)
      xv2(j)  = c1e3*myy(j)  + torby(ie)
      yv2(j)  = c1e3*myyp(j) + torbyp(ie)
      sigmv(j) = mys(j)
      ejv(j)   = myp(j)

!      as main routine will track particles back half a collimator length (to start of jaw),
!      track them now forward (if generated at face) or backward (if generated at end)
!      1/2 collimator length to center of collimator (ldrift pos or neg)
       xv1(j)  = xv1(j) - ldrift*yv1(j)
       xv2(j)  = xv2(j) - ldrift*yv2(j)

!      write out distribution - generated either at the BEGINNING or END of the collimator
!       write(4997,'(6(1X,E15.7))') myx(j), myxp(j), myy(j), myyp(j), mys(j), myp(j)
    end do
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end RB addition

!++  Copy particle data to 1-dim array and go back to meters
  do j = 1, napx
    rcx(j)  = (xv1(j)-torbx(ie)) /c1e3
    rcxp(j) = (yv1(j)-torbxp(ie))/c1e3
    rcy(j)  = (xv2(j)-torby(ie)) /c1e3
    rcyp(j) = (yv2(j)-torbyp(ie))/c1e3
    rcp(j)  = ejv(j)/c1e3
    rcs(j)  = zero
    part_hit_before_turn(j) = part_hit_turn(j)
    part_hit_before_pos(j)  = part_hit_pos(j)
    rcx0(j)  = rcx(j)
    rcxp0(j) = rcxp(j)
    rcy0(j)  = rcy(j)
    rcyp0(j) = rcyp(j)
    rcp0(j)  = rcp(j)
    ejf0v(j) = ejfv(j)

!++  For zero length element track back half collimator length
!  DRIFT PART
    if (stracki.eq.0.) then
      if(iexact.eq.0) then
        rcx(j)  = rcx(j) - half*c_length*rcxp(j)
        rcy(j)  = rcy(j) - half*c_length*rcyp(j)
      else
        zpj=sqrt(one-rcxp(j)**2-rcyp(j)**2)
        rcx(j) = rcx(j) - half*c_length*(rcxp(j)/zpj)
        rcy(j) = rcy(j) - half*c_length*(rcyp(j)/zpj)
      end if
    else
      write(lout,"(a,f13.6)") "COLL> ERROR Non-zero length collimator: '"//trim(cdb_cNameUC(icoll))//"' length = ",stracki
      call prror
    end if

    flukaname(j) = ipart(j)
  end do

!++  Do the collimation tracking
  enom_gev = myenom*c1m3

!++  Allow primaries to be one-sided, if requested
  if ((cdb_cNameUC(icoll)(1:3).eq.'TCP' .or. cdb_cNameUC(icoll)(1:3).eq.'COL') .and. do_oneside) then
    onesided = .true.
  else
    onesided = .false.
  end if

!GRD HERE IS THE MAJOR CHANGE TO THE CODE: IN ORDER TO TRACK PROPERLY THE
!GRD SPECIAL RHIC PRIMARY COLLIMATOR, IMPLEMENTATION OF A DEDICATED ROUTINE
  if(found) then
    if(cdb_cNameUC(icoll)(1:4).eq.'COLM') then
      call collimaterhic(c_material,                                    &
     &              c_length, c_rotation,                               &
     &              c_aperture, nom_aperture,                           &
     &              c_offset, c_tilt,                                   &
     &              rcx, rcxp, rcy, rcyp,                               &
     &              rcp, rcs, napx, enom_gev,                           &
     &              part_hit_pos,part_hit_turn,                         &
     &              part_abs_pos,part_abs_turn,                         &
     &              part_impact, part_indiv, part_linteract,            &
     &              onesided,                                           &
!GRD let's also add the FLUKA possibility
     &              flukaname)
    else

!GRD-SR, 09-02-2006
!Force the treatment of the TCDQ equipment as a onsided collimator.
!Both for Beam 1 and Beam 2, the TCDQ is at positive x side.
!              if(cdb_cNameUC(icoll)(1:4).eq.'TCDQ' ) onesided = .true.
! to treat all collimators onesided
! -> only for worst case TCDQ studies
      if(cdb_cNameUC(icoll)(1:4).eq.'TCDQ') onesided = .true.
      if(cdb_cNameUC(icoll)(1:5).eq.'TCXRP') onesided = .true.

!==> SLICE here is possible
!
!     SR, 29-08-2005: Slice the collimator jaws in 'n_slices' pieces
!     using two 4th-order polynomial fits. For each slices, the new
!     gaps and centre are calculates
!     It is assumed that the jaw point closer to the beam defines the
!     nominal aperture.
!
!     SR, 01-09-2005: new official version - input assigned through
!     the 'fort.3' file.
!               if (n_slices.gt.1d0 .and.                                &
!     &              totals.gt.smin_slices .and.                         &
!     &              totals.lt.smax_slices .and.                         &
!     &              cdb_cNameUC(icoll)(1:4).eq.'TCSG' ) then
!                  if (firstrun) then
!                  write(*,*) 'INFOslice - Collimator ',
!     &              cdb_cNameUC(icoll), ' sliced in ',n_slices,
!     &              ' pieces!'
!                  endif
!CB

      if(n_slices.gt.one .and. totals.gt.smin_slices .and. totals.lt.smax_slices .and. &
 &      (cdb_cNameUC(icoll)(1:4).eq.'TCSG' .or. cdb_cNameUC(icoll)(1:3).eq.'TCP' .or. cdb_cNameUC(icoll)(1:4).eq.'TCLA'.or. &
 &       cdb_cNameUC(icoll)(1:3).eq.'TCT' .or. cdb_cNameUC(icoll)(1:4).eq.'TCLI'.or. cdb_cNameUC(icoll)(1:4).eq.'TCL.'.or.  &
!     RB: added slicing of TCRYO as well
 &       cdb_cNameUC(icoll)(1:5).eq.'TCRYO')) then

        if(firstrun) then
          write(lout,*) 'INFO> slice - Collimator ', cdb_cNameUC(icoll), ' sliced in ',n_slices, ' pieces !'
        end if

!!     In this preliminary try, all secondary collimators are sliced.
!!     Slice only collimators with finite length!!
!               if (cdb_cNameUC(icoll)(1:4).eq.'TCSG' .and.
!     &              c_length.gt.0d0 ) then
!!     Slice the primaries, to have more statistics faster!
!!               if (cdb_cNameUC(icoll)(1:3).eq.'TCP' .and.
!!     +              c_length.gt.0d0 ) then
!!
!!
!!     Calculate longitudinal positions of slices and corresponding heights
!!     and angles from the fit parameters.
!!     -> MY NOTATION: y1_sl: jaw at x > 0; y2_sl: jaw at x < 0;
!!     Note: here, take (n_slices+1) points in order to calculate the
!!           tilt angle of the last slice!!

!     CB:10-2007 deformation of the jaws scaled with length
        do jjj=1,n_slices+1
          x_sl(jjj) = (jjj-1) * c_length / real(n_slices,fPrec)

          y1_sl(jjj) = fit1_1 + fit1_2*x_sl(jjj) + fit1_3/c_length*(x_sl(jjj)**2) +           &
 &                           fit1_4*(x_sl(jjj)**3) + fit1_5*(x_sl(jjj)**4) + fit1_6*(x_sl(jjj)**5)

          y2_sl(jjj) = -one * (fit2_1 + fit2_2*x_sl(jjj) + fit2_3/c_length*(x_sl(jjj)**2) +   &
 &                           fit2_4*(x_sl(jjj)**3) + fit2_5*(x_sl(jjj)**4) + fit2_6*(x_sl(jjj)**5))
        end do

!       Apply the slicing scaling factors (ssf's):
!       CB:10-2007 coordinates rotated of the tilt
        do jjj=1,n_slices+1
          y1_sl(jjj) = ssf1 * y1_sl(jjj)
          y2_sl(jjj) = ssf2 * y2_sl(jjj)
! CB code
          x1_sl(jjj) = x_sl(jjj) *cos_mb(db_tilt(icoll,1))-y1_sl(jjj)*sin_mb(db_tilt(icoll,1))
          x2_sl(jjj) = x_sl(jjj) *cos_mb(db_tilt(icoll,2))-y2_sl(jjj)*sin_mb(db_tilt(icoll,2))
          y1_sl(jjj) = y1_sl(jjj)*cos_mb(db_tilt(icoll,1))+x_sl(jjj) *sin_mb(db_tilt(icoll,1))
          y2_sl(jjj) = y2_sl(jjj)*cos_mb(db_tilt(icoll,2))+x_sl(jjj) *sin_mb(db_tilt(icoll,2))
        end do

!       Sign of the angle defined differently for the two jaws!
        do jjj=1,n_slices
          angle1(jjj) = (( y1_sl(jjj+1) - y1_sl(jjj) ) / ( x1_sl(jjj+1)-x1_sl(jjj) ))
          angle2(jjj) = (( y2_sl(jjj+1) - y2_sl(jjj) ) / ( x2_sl(jjj+1)-x2_sl(jjj) ))
        end do

!       Sign of the angle defined differently for the two jaws!
!                    do jjj=1,n_slices
!                       angle1(jjj) = ( y1_sl(jjj+1) - y1_sl(jjj) ) /     &
!       &                    (c_length / dble(n_slices) )
!                       angle2(jjj) = ( y2_sl(jjj+1) - y2_sl(jjj) ) /     &
!       &                    (c_length / dble(n_slices) )
!                    enddo
!       For both jaws, look for the 'deepest' point (closest point to beam)
!       Then, shift the vectors such that this closest point defines
!       the nominal aperture
!       Index here must go up to (n_slices+1) in case the last point is the
!       closest (and also for the later calculation of 'a_tmp1' and 'a_tmp2')

!       SR, 01-09-2005: add the recentring flag, as given in 'fort.3' to
!       choose whether recentre the deepest point or not
        max_tmp = c1e6
        do jjj=1, n_slices+1
          if( y1_sl(jjj).lt.max_tmp ) then
            max_tmp = y1_sl(jjj)
          end if
        end do

        do jjj=1, n_slices+1
          y1_sl(jjj) = y1_sl(jjj) - (max_tmp * recenter1) + (half*c_aperture)
        end do

        max_tmp = -c1e6

        do jjj=1, n_slices+1
          if( y2_sl(jjj).gt.max_tmp ) then
            max_tmp = y2_sl(jjj)
          end if
        end do

        do jjj=1, n_slices+1
          y2_sl(jjj) = y2_sl(jjj) - (max_tmp * recenter2) - (half*c_aperture)
        end do

!!      Check the collimator jaw surfaces (beam frame, before taking into
!!      account the azimuthal angle of the collimator)
        if(firstrun) then
          write(lout,*) 'Slicing collimator ',cdb_cNameUC(icoll)
           do jjj=1,n_slices
             write(lout,*) x_sl(jjj), y1_sl(jjj), y2_sl(jjj), angle1(jjj), angle2(jjj), db_tilt(icoll,1), db_tilt(icoll,2)
           end do
        end if
!
!!     Check the calculation of slice gap and centre
!                  if (firstrun) then
!                     write(*,*) 'Verify centre and gap!'
!                     do jjj=1,n_slices
!                        if ( angle1(jjj).gt.0d0 ) then
!                           a_tmp1 = y1_sl(jjj)
!                        else
!                           a_tmp1 = y1_sl(jjj+1)
!                        endif
!                        if ( angle2(jjj).lt.0d0 ) then
!                           a_tmp2 = y2_sl(jjj)
!                        else
!                           a_tmp2 = y2_sl(jjj+1)
!                        endif
!                        write(*,*) a_tmp1 - a_tmp2,
!     +                       0.5 * ( a_tmp1 + a_tmp2 )
!                     enddo
!                  endif
!
!       Now, loop over the number of slices and call collimate2 each time!
!       For each slice, the corresponding offset and angle are to be used.
        do jjj=1,n_slices

!         First calculate aperture and centre of the slice
!         Note that:
!         (1)due to our notation for the angle sign,
!         the rotation point of the slice (index j or j+1)
!         DEPENDS on the angle value!!
!         (2) New version of 'collimate2' is required: one must pass
!         the slice number in order the calculate correctly the 's'
!         coordinate in the impact files.

!         Here, 'a_tmp1' and 'a_tmp2' are, for each slice, the closest
!         corners to the beam
          if( angle1(jjj).gt.zero ) then
            a_tmp1 = y1_sl(jjj)
          else
            a_tmp1 = y1_sl(jjj+1)
          end if

          if( angle2(jjj).lt.zero ) then
            a_tmp2 = y2_sl(jjj)
          else
            a_tmp2 = y2_sl(jjj+1)
          end if

!!     Write down the information on slice centre and offset
!                     if (firstrun) then
!                        write(*,*) 'Processing slice number ',jjj,
!     &                       ' of ',n_slices,' for the collimator ',
!     &                       cdb_cNameUC(icoll)
!                        write(*,*) 'Aperture [m]= ',
!     &                       a_tmp1 - a_tmp2
!                        write(*,*) 'Offset [m]  = ',
!     &                       0.5 * ( a_tmp1 + a_tmp2 )
!                     endif
!!
!     Be careful! the initial tilt must be added!
!     We leave it like this for the moment (no initial tilt)
!         c_tilt(1) = c_tilt(1) + angle1(jjj)
!         c_tilt(2) = c_tilt(2) + angle2(jjj)
          c_tilt(1) = angle1(jjj)
          c_tilt(2) = angle2(jjj)
!     New version of 'collimate2' is required: one must pass the
!     slice number in order the calculate correctly the 's'
!     coordinate in the impact files.
!     +                    a_tmp1 - a_tmp2,
!     +                    0.5 * ( a_tmp1 + a_tmp2 ),
! -- TW SEP07 added compatility for tilt, gap and ofset errors to slicing
! -- TW gaprms error is already included in the c_aperture used above
! -- TW tilt error is added to y1_sl and y2_sl therfore included in
! -- TW angle1 and angle2 no additinal changes needed
! -- TW offset error directly added to call of collimate2

! --- TW JUNE08
          if (firstrun) then
            write(collsettings_unit,'(a,1x,i10,5(1x,e13.5),1x,a)')      &
     &                       cdb_cNameUC(icoll)(1:12),                     &
     &                       jjj,                                       &
     &                       (a_tmp1 - a_tmp2)/two,                     &
     &                       half * (a_tmp1 + a_tmp2) + c_offset,       &
     &                       c_tilt(1),                                 &
     &                       c_tilt(2),                                 &
     &                       c_length / real(n_slices,fPrec),           &
     &                       cdb_cMaterial(icoll)
          end if
! --- TW JUNE08
                     call collimate2(c_material,                        &
     &                    c_length / real(n_slices,fPrec),              &
     &                    c_rotation,                                   &
     &                    a_tmp1 - a_tmp2,                              &
     &                    half * ( a_tmp1 + a_tmp2 ) + c_offset,        &
     &                    c_tilt,                                       &
     &                    rcx, rcxp, rcy, rcyp,                         &
     &                    rcp, rcs, napx, enom_gev,                     &
     &                    part_hit_pos, part_hit_turn,                  &
     &                    part_abs_pos, part_abs_turn,                  &
     &                    part_impact, part_indiv,                      &
     &                    part_linteract, onesided, flukaname,          &
     &                    secondary,                                    &
     &                    jjj, nabs_type)
        end do !do jjj=1,n_slices
      else !if(n_slices.gt.one .and. totals.gt.smin_slices .and. totals.lt.smax_slices .and.
!     Treatment of non-sliced collimators

#ifdef G4COLLIMAT
!! Add the geant4 geometry
        if(firstrun.and.iturn.eq.1) then
          call g4_add_collimator(cdb_cNameUC(icoll), c_material, c_length, c_aperture, c_rotation, c_offset)
        endif

!! Here we do the real collimation
!! First set the correct collimator
        call g4_set_collimator(cdb_cNameUC(icoll))
        flush(lout)

!! Loop over all our particles
        g4_lostc = 0
        do j = 1, napx
          if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
!! Rotate particles in the frame of the collimator
!! There is more precision if we do it here rather
!! than in the g4 geometry
            x_tmp = rcx(j)
            y_tmp = rcy(j)
            xp_tmp = rcxp(j)
            yp_tmp = rcyp(j)
            rcx(j) = x_tmp*cos_mb(c_rotation) +sin_mb(c_rotation)*y_tmp
            rcy(j) = y_tmp*cos_mb(c_rotation) -sin_mb(c_rotation)*x_tmp
            rcxp(j) = xp_tmp*cos_mb(c_rotation)+sin_mb(c_rotation)*yp_tmp
            rcyp(j) = yp_tmp*cos_mb(c_rotation)-sin_mb(c_rotation)*xp_tmp

!! Call the geant4 collimation function
            call g4_collimate(rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j))

!! Get the particle back + information
            call g4_collimate_return(rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j), part_hit_flag, part_abs_flag, &
 &                                   part_impact(j), part_indiv(j), part_linteract(j))

!! Rotate back into the accelerator frame
            x_tmp   = rcx(j)
            y_tmp   = rcy(j)
            xp_tmp  = rcxp(j)
            yp_tmp  = rcyp(j)
            rcx(j)  = x_tmp *cos_mb(-one*c_rotation) + sin_mb(-one*c_rotation)*y_tmp
            rcy(j)  = y_tmp *cos_mb(-one*c_rotation) - sin_mb(-one*c_rotation)*x_tmp
            rcxp(j) = xp_tmp*cos_mb(-one*c_rotation) + sin_mb(-one*c_rotation)*yp_tmp
            rcyp(j) = yp_tmp*cos_mb(-one*c_rotation) - sin_mb(-one*c_rotation)*xp_tmp

!           If a particle hit
            if(part_hit_flag.ne.0) then
              part_hit_pos(j) = ie
              part_hit_turn(j) = iturn
            end if

!           If a particle died (the checking if it is already dead is at the start of the loop)
!           Geant just has a general inelastic process that single diffraction is part of
!           Therefore we can not know if this interaction was SD or some other inelastic type
            if(part_abs_flag.ne.0) then
              if(dowrite_impact) then
!! FLUKA_impacts.dat
                write(FLUKA_impacts_unit,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))') &
 &                    icoll,c_rotation,zero,zero,zero,zero,zero,part_abs_flag,flukaname(j),iturn
              end if

              part_abs_pos(j)  = ie
              part_abs_turn(j) = iturn
              rcx(j) = 99.99e-3_fPrec
              rcy(j) = 99.99e-3_fPrec
              g4_lostc = g4_lostc + 1
            end if
          flush(lout)
          end if !part_abs_pos(j) .ne. 0 .and. part_abs_turn(j) .ne. 0
        end do   !do j = 1, napx
!      write(lout,*) 'COLLIMATOR LOSSES ', cdb_cNameUC(icoll), g4_lostc
#endif
#ifndef G4COLLIMAT
! This is what is called in a normal collimation run
                  call collimate2(c_material, c_length, c_rotation,     &
     &                 c_aperture, c_offset, c_tilt,                    &
     &                 rcx, rcxp, rcy, rcyp,                            &
     &                 rcp, rcs, napx, enom_gev,                        &
     &                 part_hit_pos,part_hit_turn,                      &
     &                 part_abs_pos, part_abs_turn,                     &
     &                 part_impact, part_indiv, part_linteract,         &
     &                 onesided, flukaname, secondary, 1, nabs_type)
#endif
      end if !if (n_slices.gt.one .and.
    end if !if(cdb_cNameUC(icoll)(1:4).eq.'COLM') then
  end if !if (found) then
end subroutine collimate_do_collimator

!>
!! collimate_end_collimator()
!! This routine is called at the exit of a collimator
!<
subroutine collimate_end_collimator()

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use numerical_constants, only : c5m4

  implicit none

  integer :: j

#ifdef HDF5
  ! For tracks2
  integer hdfturn,hdfpid,hdftyp
  real(kind=fPrec) hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdfs
#endif

  ! real(kind=fPrec) stracki ! stracki makes no sense here

!++  Output information:
!++
!++  PART_HIT_POS (MAX_NPART)  Hit flag for last hit
!++  PART_HIT_TURN(MAX_NPART)  Hit flag for last hit
!++  PART_ABS_POS (MAX_NPART)  Abs flag
!++  PART_ABS_TURN(MAX_NPART)  Abs flag
!++  PART_IMPACT  (MAX_NPART)  Impact parameter (0 for inner face)
!++  PART_INDIV   (MAX_NPART)  Divergence of impacting particles
!------------------------------------------------------------------------------
!++  Calculate average impact parameter and save info for all
!++  collimators. Copy information back and do negative drift.
  n_impact = 0
  n_absorbed = 0
  sum      = zero
  sqsum    = zero

!++  Copy particle data back and do path length stuff; check for absorption
!++  Add orbit offset back.
  do j = 1, napx

!APRIL2005 IN ORDER TO GET RID OF NUMERICAL ERRORS, JUST DO THE TREATMENT FOR
!APRIL2005 IMPACTING PARTICLES...
    if(part_hit_pos(j) .eq.ie .and. part_hit_turn(j).eq.iturn) then
!++  For zero length element track back half collimator length
! DRIFT PART
      ! if (stracki.eq.0.) then ! stracki makes no sense here
        if(iexact.eq.0) then
          rcx(j)  = rcx(j) - half*c_length*rcxp(j)
          rcy(j)  = rcy(j) - half*c_length*rcyp(j)
        else
          zpj=sqrt(one-rcxp(j)**2-rcyp(j)**2)
          rcx(j) = rcx(j) - half*c_length*(rcxp(j)/zpj)
          rcy(j) = rcy(j) - half*c_length*(rcyp(j)/zpj)
        end if
      ! end if ! stracki makes no sense here

!++  Now copy data back to original verctor
      xv1(j) = rcx(j)  * c1e3 + torbx(ie)
      yv1(j) = rcxp(j) * c1e3 + torbxp(ie)
      xv2(j) = rcy(j)  * c1e3 + torby(ie)
      yv2(j) = rcyp(j) * c1e3 + torbyp(ie)
      ejv(j)  = rcp(j)  * c1e3

!++  Energy update, as recommended by Frank
      ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
      rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
      dpsv(j)=(ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
      oidpsv(j)=one/(one+dpsv(j))
      moidpsv(j)=mtc(j)/(one+dpsv(j))
      omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
      dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)
      yv1(j)   = ejf0v(j)/ejfv(j)*yv1(j)
      yv2(j)   = ejf0v(j)/ejfv(j)*yv2(j)

!++   For absorbed particles set all coordinates to zero. Also
!++   include very large offsets, let's say above 100mm or
!++   100mrad.
      if( (part_abs_pos(j).ne.0 .and. part_abs_turn(j).ne.0) .or.&
 &      xv1(j).gt.c1e2 .or. yv1(j).gt.c1e2 .or. xv2(j).gt.c1e2 .or. yv2(j).gt.c1e2) then
        xv1(j) = zero
        yv1(j) = zero
        xv2(j) = zero
        yv2(j) = zero
        ejv(j)  = myenom
        sigmv(j)= zero
        part_abs_pos(j)=ie
        part_abs_turn(j)=iturn
        secondary(j) = 0
        tertiary(j)  = 0
        other(j)     = 0
        scatterhit(j)= 0
        nabs_type(j) = 0
      end if

!APRIL2005 ...OTHERWISE JUST GET BACK FORMER COORDINATES
    else
      xv1(j) = rcx0(j)  * c1e3 + torbx(ie)
      yv1(j) = rcxp0(j) * c1e3 + torbxp(ie)
      xv2(j) = rcy0(j)  * c1e3 + torby(ie)
      yv2(j) = rcyp0(j) * c1e3 + torbyp(ie)
      ejv(j)  = rcp0(j)  * c1e3
    end if

!++  First check for particle interaction at this collimator and this turn
    if(part_hit_pos (j).eq.ie .and. part_hit_turn(j).eq.iturn) then

!++  Fill the change in particle angle into histogram
      if(dowrite_impact) then
#ifdef HDF5
        if(h5_useForCOLL) then
          call h5_prepareWrite(coll_hdf5_allImpacts, 1)
          call h5_writeData(coll_hdf5_allImpacts, 1, 1, ipart(j))
          call h5_writeData(coll_hdf5_allImpacts, 2, 1, iturn)
          call h5_writeData(coll_hdf5_allImpacts, 3, 1, sampl(ie))
          call h5_finaliseWrite(coll_hdf5_allImpacts)
        else
#endif
          write(all_impacts_unit,'(i8,1x,i4,1x,f8.2)') ipart(j),iturn,sampl(ie)
#ifdef HDF5
        end if
#endif
      end if

      ! Particle has impacted
      if(part_abs_pos(j) .ne.0 .and. part_abs_turn(j).ne.0) then
        if(dowrite_impact) then
#ifdef HDF5
          if(h5_useForCOLL) then
            call h5_prepareWrite(coll_hdf5_allAbsorb, 1)
            call h5_writeData(coll_hdf5_allAbsorb, 1, 1, ipart(j))
            call h5_writeData(coll_hdf5_allAbsorb, 2, 1, iturn)
            call h5_writeData(coll_hdf5_allAbsorb, 3, 1, sampl(ie))
            call h5_finaliseWrite(coll_hdf5_allAbsorb)
          else
#endif
            write(all_absorptions_unit,'(i8,1x,i4,1x,f8.2)') ipart(j),iturn,sampl(ie)
#ifdef HDF5
          end if
#endif
        end if

      !Here we've found a newly hit particle
      else if(part_abs_pos (j).eq.0 .and.  part_abs_turn(j).eq.0) then
        xkick = rcxp(j) - rcxp0(j)
        ykick = rcyp(j) - rcyp0(j)

        ! Indicate wether this is a secondary / tertiary / other particle;
        !  note that 'scatterhit' (equals 8 when set) is set in SCATTER.
        if(cdb_cNameUC(icoll)(1:3).eq.'TCP'   .or. &
           cdb_cNameUC(icoll)(1:4).eq.'COLM'  .or. &
           cdb_cNameUC(icoll)(1:5).eq.'COLH0' .or. &
           cdb_cNameUC(icoll)(1:5).eq.'COLV0'       ) then
          secondary(j) = 1
        else if(cdb_cNameUC(icoll)(1:3).eq.'TCS'   .or. &
                cdb_cNameUC(icoll)(1:4).eq.'COLH1' .or. &
                cdb_cNameUC(icoll)(1:4).eq.'COLV1' .or. &
                cdb_cNameUC(icoll)(1:4).eq.'COLH2'       ) then
          tertiary(j)  = 2
       else if((cdb_cNameUC(icoll)(1:3).eq.'TCL') .or. &
               (cdb_cNameUC(icoll)(1:3).eq.'TCT') .or. &
               (cdb_cNameUC(icoll)(1:3).eq.'TCD') .or. &
               (cdb_cNameUC(icoll)(1:3).eq.'TDI')       ) then
          other(j)     = 4
        end if
      else
        write(lout,*) "Error in collimate_end_collimator"
        write(lout,*) "Particle cannot be both absorbed and not absorbed."
        write(lout,*) part_abs_pos (j),  part_abs_turn(j)
        call prror(-1)
      end if

!GRD THIS LOOP MUST NOT BE WRITTEN INTO THE "IF(FIRSTRUN)" LOOP !!!!!
      if(dowritetracks) then
        if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
          if((secondary(j) .eq. 1 .or. &
              tertiary(j)  .eq. 2 .or. &
              other(j)     .eq. 4 .or. &
              scatterhit(j).eq.8         ) .and. &
             (xv1(j).lt.99.0_fPrec .and. xv2(j).lt.99.0_fPrec).and.&
!GRD HERE WE APPLY THE SAME KIND OF CUT THAN THE SIGSECUT PARAMETER
             ((((xv1(j)*c1m3)**2 / (tbetax(ie)*myemitx0_collgap)) .ge. real(sigsecut2,fPrec)) .or. &
             (((xv2(j)*c1m3)**2  / (tbetay(ie)*myemity0_collgap)) .ge. real(sigsecut2,fPrec)) .or. &
             (((xv1(j)*c1m3)**2  / (tbetax(ie)*myemitx0_collgap)) + &
             ((xv2(j)*c1m3)**2   / (tbetay(ie)*myemity0_collgap)) .ge. sigsecut3)) ) &
             then

            xj  = (xv1(j)-torbx(ie))  /c1e3
            xpj = (yv1(j)-torbxp(ie)) /c1e3
            yj  = (xv2(j)-torby(ie))  /c1e3
            ypj = (yv2(j)-torbyp(ie)) /c1e3

#ifdef HDF5
            if(h5_writeTracks2) then
              ! We write trajectories before and after element in this case.
              hdfpid  = ipart(j)
              hdfturn = iturn
              hdfs    = sampl(ie)-half*c_length
              hdfx    = (rcx0(j)*c1e3+torbx(ie)) - half*c_length*(rcxp0(j)*c1e3+torbxp(ie))
              hdfxp   = rcxp0(j)*c1e3+torbxp(ie)
              hdfy    = (rcy0(j)*c1e3+torby(ie)) - half*c_length*(rcyp0(j)*c1e3+torbyp(ie))
              hdfyp   = rcyp0(j)*c1e3+torbyp(ie)
              hdfdee  = (ejv(j)-myenom)/myenom
              hdftyp  = secondary(j)+tertiary(j)+other(j)+scatterhit(j)
              call h5tr2_writeLine(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdftyp)

              hdfs  = sampl(ie)+half*c_length
              hdfx  = xv1(j) + half*c_length*yv1(j)
              hdfxp = yv1(j)
              hdfy  = xv2(j) + half*c_length*yv2(j)
              hdfyp = yv2(j)
              call h5tr2_writeLine(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdftyp)
            else
#endif
              write(tracks2_unit,'(1x,i8,1x,i4,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)') &
                ipart(j),iturn,sampl(ie)-half*c_length,           &
                (rcx0(j)*c1e3+torbx(ie))-half*c_length*(rcxp0(j)*c1e3+torbxp(ie)), &
                rcxp0(j)*c1e3+torbxp(ie),                                          &
                (rcy0(j)*c1e3+torby(ie))-half*c_length*(rcyp0(j)*c1e3+torbyp(ie)), &
                rcyp0(j)*c1e3+torbyp(ie),                                          &
                (ejv(j)-myenom)/myenom,secondary(j)+tertiary(j)+other(j)+scatterhit(j)

              write(tracks2_unit,'(1x,i8,1x,i4,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)') &
                ipart(j),iturn,sampl(ie)+half*c_length,           &
                xv1(j)+half*c_length*yv1(j),yv1(j),                             &
                xv2(j)+half*c_length*yv2(j),yv2(j),(ejv(j)-myenom)/myenom,      &
                secondary(j)+tertiary(j)+other(j)+scatterhit(j)
#ifdef HDF5
            end if
#endif
          end if ! if((secondary(j).eq.1.or.tertiary(j).eq.2.or.other(j).eq.4)
          ! .and.(xv1(j).lt.99.0_fPrec.and.xv2(j).lt.99.0_fPrec) and.
        end if !if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
      end if !if(dowritetracks) then

!++  Calculate impact observables, fill histograms, save collimator info, ...
      n_impact = n_impact + 1
      sum = sum + part_impact(j)
      sqsum = sqsum + part_impact(j)**2
      cn_impact(icoll) = cn_impact(icoll) + 1
      csum(icoll) = csum(icoll) + part_impact(j)
      csqsum(icoll) = csqsum(icoll) + part_impact(j)**2

!++  If the interacting particle was lost, add-up counters for absorption
!++  Note: a particle with x/y >= 99. never hits anything any more in
!++        the logic of this program. Be careful to always fulfill this!
      if(part_abs_pos(j).ne.0 .and. part_abs_turn(j).ne.0) then
        n_absorbed = n_absorbed + 1
        cn_absorbed(icoll) = cn_absorbed(icoll) + 1
        n_tot_absorbed = n_tot_absorbed + 1
        iturn_last_hit = part_hit_before_turn(j)
        iturn_absorbed = part_hit_turn(j)
        if(iturn_last_hit.eq.0) then
          iturn_last_hit = iturn_absorbed
          iturn_survive  = iturn_absorbed - iturn_last_hit
        end if
      end if

!++  End of check for hit this turn and element
    end if
  end do ! end do j = 1, napx

!++  Calculate statistical observables and save into files...
  if (n_impact.gt.0) then
    average = sum/n_impact

    if (sqsum/n_impact.ge.average**2) then
      sigma = sqrt(sqsum/n_impact - average**2)
    else
      sigma = zero
    end if
  else
    average = zero
    sigma   = zero
  end if

  if(cn_impact(icoll).gt.0) then
    caverage(icoll) = csum(icoll)/cn_impact(icoll)

    if((caverage(icoll)**2).gt.(csqsum(icoll)/cn_impact(icoll))) then
      csigma(icoll) = 0
    else
      csigma(icoll) = sqrt(csqsum(icoll)/cn_impact(icoll) - caverage(icoll)**2)
    end if
  end if

!-----------------------------------------------------------------
!++  For a  S E L E C T E D  collimator only consider particles that
!++  were scattered on this selected collimator at the first turn. All
!++  other particles are discarded.
!++  - This is switched on with the DO_SELECT flag in the input file.
!++  - Note that the part_select(j) flag defaults to 1 for all particles.

! should name_sel(1:11) extended to allow longer names as done for
! coll the coll_ellipse.dat file !!!!!!!!
  if(((cdb_cNameUC(icoll).eq.name_sel(1:mNameLen)).or.&
      (cdb_cName(icoll).eq.name_sel(1:mNameLen))) .and. iturn.eq.1  ) then
    num_selhit = 0
    num_surhit = 0
    num_selabs = 0

    do j = 1, napx
      if( part_hit_pos (j).eq.ie .and. part_hit_turn(j).eq.iturn ) then

      num_selhit = num_selhit+1

      if(part_abs_pos(j) .eq.0 .and. part_abs_turn(j).eq.0) then
        num_surhit = num_surhit+1
      else
        num_selabs = num_selabs + 1
      end if

!++  If we want to select only partciles interacting at the specified
!++  collimator then remove all other particles and reset the number
!++  of the absorbed particles to the selected collimator.
      else if(do_select.and.firstrun) then
        part_select(j) = 0
        n_tot_absorbed = num_selabs
      end if
    end do

!++  Calculate average impact parameter and save distribution into file
!++  only for selected collimator
    n_impact = 0
    sum      = zero
    sqsum    = zero

    do j = 1, napx
      if( part_hit_pos (j).eq.ie .and. part_hit_turn(j).eq.iturn ) then
        if(part_impact(j).lt.-half) then
          write(lout,*) 'ERR>  Found invalid impact parameter!', part_impact(j)
          write(outlun,*) 'ERR>  Invalid impact parameter!', part_impact(j)
          call prror(-1)
        end if

        n_impact = n_impact + 1
        sum = sum + part_impact(j)
        sqsum = sqsum + part_impact(j)**2
        if(part_hit_pos (j).ne.0 .and. part_hit_turn(j).ne.0 .and.dowrite_impact ) then
          write(impact_unit,*) part_impact(j), part_indiv(j)
        end if
      end if
    end do

    if(n_impact.gt.0) then
      average = sum/n_impact
      if(sqsum/n_impact.ge.average**2) then
        sigma = sqrt(sqsum/n_impact - average**2)
      else
        sigma = zero
      end if
    end if

!++  Some information
    write(lout,"(a,i8)")    'COLL> Selected collimator had N hits. N: ', num_selhit
    write(lout,"(a,i8)")    'COLL> Number of impacts                : ', n_impact
    write(lout,"(a,i8)")    'COLL> Number of escaped protons        : ', num_surhit
    write(lout,"(a,e15.8)") 'COLL> Average impact parameter [m]     : ', average
    write(lout,"(a,e15.8)") 'COLL> Sigma impact parameter [m]       : ', sigma

    if (dowrite_impact) then
      close(impact_unit)
    end if

!++  End of    S E L E C T E D   collimator
  end if

end subroutine collimate_end_collimator

!>
!! collimate_end_sample()
!! This routine is called from trauthin after each sample
!! has been tracked by thin6d
!<
subroutine collimate_end_sample(j)

  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use crcoall
#ifdef ROOT
  use root_output
#endif

  implicit none

  integer, intent(in) :: j

#ifdef HDF5
  type(h5_dataField), allocatable :: fldHdf(:)
  integer fmtHdf, setHdf
#endif

!++  Save particle offsets to a file
  ! close(beta_beat_unit)
  close(survival_unit)

  if(dowrite_impact) close(impact_unit)

  if(dowritetracks) then
    if(cern) close(tracks2_unit)
#ifdef HDF5
    if(cern .and. h5_writeTracks2) call h5tr2_finalise
#endif
  end if

!------------------------------------------------------------------------
!++  Write the number of absorbed particles
  write(outlun,*) 'INFO>  Number of impacts             : ', n_tot_absorbed+nsurvive_end
  write(outlun,*) 'INFO>  Number of impacts at selected : ', num_selhit
  write(outlun,*) 'INFO>  Number of surviving particles : ', nsurvive_end
  write(outlun,*) 'INFO>  Number of absorbed particles  : ', n_tot_absorbed
  write(outlun,*)

  if(n_tot_absorbed.ne.0) then                                       !hr08
    write(outlun,*) ' INFO>  Eff_r @  8 sigma    [e-4] : ', (neff(5)/real(n_tot_absorbed,fPrec))/c1m4              !hr08
    write(outlun,*) ' INFO>  Eff_r @ 10 sigma    [e-4] : ', (neff(9)/real(n_tot_absorbed,fPrec))/c1m4              !hr08
    write(outlun,*) ' INFO>  Eff_r @ 10-20 sigma [e-4] : ', ((neff(9)-neff(19))/(real(n_tot_absorbed,fPrec)))/c1m4 !hr08
    write(outlun,*)
    write(outlun,*) neff(5)/real(n_tot_absorbed,fPrec), neff(9)/real(n_tot_absorbed,fPrec), &
 & (neff(9)-neff(19))/(real(n_tot_absorbed,fPrec)), ' !eff'
    write(outlun,*)
  else
    write(lout,*) 'NO PARTICLE ABSORBED'
  endif

  write(lout,"(a)")
  write(lout,"(a,i8)") 'COLL> Number of impacts             : ', n_tot_absorbed+nsurvive_end
  write(lout,"(a,i8)") 'COLL> Number of impacts at selected : ', num_selhit
  write(lout,"(a,i8)") 'COLL> Number of surviving particles : ', nsurvive_end
  write(lout,"(a,i8)") 'COLL> Number of absorbed particles  : ', n_tot_absorbed
  write(lout,"(a)")

  if(n_tot_absorbed.ne.0) then
    write(lout,"(a,f20.12)") 'COLL> Eff_r @  8 sigma    [e-4] : ', (neff(5)/real(n_tot_absorbed,fPrec))/c1m4
    write(lout,"(a,f20.12)") 'COLL> Eff_r @ 10 sigma    [e-4] : ', (neff(9)/real(n_tot_absorbed,fPrec))/c1m4
    write(lout,"(a,f20.12)") 'COLL> Eff_r @ 10-20 sigma [e-4] : ', ((neff(9)-neff(19))/real(n_tot_absorbed,fPrec))/c1m4
  else
    write(lout,"(a)") 'COLL> No particle absorbed'
  endif
  write(lout,"(a)")

! Write efficiency file
#ifdef HDF5
  if(h5_useForCOLL .and. n_tot_absorbed /= 0) then
    allocate(fldHdf(8))
    fldHdf(1) = h5_dataField(name="RAD_SIGMA",  type=h5_typeReal)
    fldHdf(2) = h5_dataField(name="NEFFX/NTOT", type=h5_typeReal)
    fldHdf(3) = h5_dataField(name="NEFFY/NTOT", type=h5_typeReal)
    fldHdf(4) = h5_dataField(name="NEFF/NTOT",  type=h5_typeReal)
    fldHdf(5) = h5_dataField(name="NEFFX",      type=h5_typeReal)
    fldHdf(6) = h5_dataField(name="NEFFY",      type=h5_typeReal)
    fldHdf(7) = h5_dataField(name="NEFF",       type=h5_typeReal)
    fldHdf(8) = h5_dataField(name="NTOT",       type=h5_typeInt)
    call h5_createFormat("collEfficiency", fldHdf, fmtHdf)
    call h5_createDataSet("efficiency", h5_collID, fmtHdf, setHdf, numeff)
    call h5_prepareWrite(setHdf, numeff)
    call h5_writeData(setHdf, 1, numeff, rsig(1:numeff))
    call h5_writeData(setHdf, 2, numeff, neffx(1:numeff)/real(n_tot_absorbed,real64))
    call h5_writeData(setHdf, 3, numeff, neffy(1:numeff)/real(n_tot_absorbed,real64))
    call h5_writeData(setHdf, 4, numeff, neff(1:numeff)/real(n_tot_absorbed,real64))
    call h5_writeData(setHdf, 5, numeff, neffx(1:numeff))
    call h5_writeData(setHdf, 6, numeff, neffy(1:numeff))
    call h5_writeData(setHdf, 7, numeff, neff(1:numeff))
    call h5_writeData(setHdf, 8, numeff, n_tot_absorbed)
    call h5_finaliseWrite(setHdf)
    deallocate(fldHdf)
  else
#endif
    call f_requestUnit('efficiency.dat', efficiency_unit)
    open(unit=efficiency_unit, file='efficiency.dat') !was 1991
    if(n_tot_absorbed /= 0) then
      write(efficiency_unit,*) '# 1=rad_sigma 2=frac_x 3=frac_y 4=frac_r' ! This is not correct?
      do k=1,numeff
        write(efficiency_unit,'(7(1x,e15.7),1x,I5)') rsig(k), neffx(k)/real(n_tot_absorbed,fPrec), &
 & neffy(k)/real(n_tot_absorbed,fPrec), neff(k)/real(n_tot_absorbed,fPrec), neffx(k), neffy(k), neff(k), n_tot_absorbed
      end do
    else
      write(lout,*) 'NO PARTICLE ABSORBED'
    end if
    close(efficiency_unit)
#ifdef HDF5
  end if
#endif

! Write efficiency vs dp/p file
#ifdef HDF5
  if(h5_useForCOLL .and. n_tot_absorbed /= 0) then
    allocate(fldHdf(5))
    fldHdf(1) = h5_dataField(name="DP/P",        type=h5_typeReal)
    fldHdf(2) = h5_dataField(name="NDPOP/TNABS", type=h5_typeReal)
    fldHdf(3) = h5_dataField(name="NDPOP",       type=h5_typeReal)
    fldHdf(4) = h5_dataField(name="TNABS",       type=h5_typeInt)
    fldHdf(5) = h5_dataField(name="NPART",       type=h5_typeInt)
    call h5_createFormat("collEfficiencyDPOP", fldHdf, fmtHdf)
    call h5_createDataSet("efficiency_dpop", h5_collID, fmtHdf, setHdf, numeffdpop)
    call h5_prepareWrite(setHdf, numeffdpop)
    call h5_writeData(setHdf, 1, numeffdpop, dpopbins(1:numeffdpop))
    call h5_writeData(setHdf, 2, numeffdpop, neffdpop(1:numeffdpop)/real(n_tot_absorbed,real64))
    call h5_writeData(setHdf, 3, numeffdpop, neffdpop(1:numeffdpop))
    call h5_writeData(setHdf, 4, numeffdpop, n_tot_absorbed)
    call h5_writeData(setHdf, 5, numeffdpop, npartdpop(1:numeffdpop))
    call h5_finaliseWrite(setHdf)
    deallocate(fldHdf)
  else
#endif
    call f_requestUnit('efficiency_dpop.dat', efficiency_dpop_unit)
    open(unit=efficiency_dpop_unit, file='efficiency_dpop.dat') !was 1992
    if(n_tot_absorbed /= 0) then
      write(efficiency_dpop_unit,*) '# 1=dp/p 2=n_dpop/tot_nabs 3=n_dpop 4=tot_nabs 5=npart'
      do k=1,numeffdpop
        write(efficiency_dpop_unit,'(3(1x,e15.7),2(1x,I5))') dpopbins(k), neffdpop(k)/real(n_tot_absorbed,fPrec), neffdpop(k), &
            n_tot_absorbed, npartdpop(k)
      end do
    else
      write(lout,*) 'NO PARTICLE ABSORBED'
    end if
    close(efficiency_dpop_unit)
#ifdef HDF5
  end if
#endif

! Write 2D efficiency file (eff vs. A_r and dp/p)
#ifdef HDF5
  if(h5_useForCOLL .and. n_tot_absorbed /= 0) then
    allocate(fldHdf(5))
    fldHdf(1) = h5_dataField(name="RAD_SIGMA", type=h5_typeReal)
    fldHdf(2) = h5_dataField(name="DP/P",      type=h5_typeReal)
    fldHdf(3) = h5_dataField(name="N/TNABS",   type=h5_typeReal)
    fldHdf(4) = h5_dataField(name="N",         type=h5_typeReal)
    fldHdf(5) = h5_dataField(name="TNABS",     type=h5_typeInt)
    call h5_createFormat("collEfficiency2D", fldHdf, fmtHdf)
    call h5_createDataSet("efficiency_2d", h5_collID, fmtHdf, setHdf, numeffdpop)
    do i=1,numeff
      call h5_prepareWrite(setHdf, numeffdpop)
      call h5_writeData(setHdf, 1, numeffdpop, rsig(i))
      call h5_writeData(setHdf, 2, numeffdpop, dpopbins(1:numeffdpop))
      call h5_writeData(setHdf, 3, numeffdpop, neff2d(i,1:numeffdpop)/real(n_tot_absorbed,fPrec))
      call h5_writeData(setHdf, 4, numeffdpop, neff2d(i,1:numeffdpop))
      call h5_writeData(setHdf, 5, numeffdpop, n_tot_absorbed)
      call h5_finaliseWrite(setHdf)
    end do
    deallocate(fldHdf)
  else
#endif
    call f_requestUnit('efficiency_2d.dat', efficiency_2d_unit)
    open(unit=efficiency_2d_unit, file='efficiency_2d.dat') !was 1993
    if(n_tot_absorbed /= 0) then
      write(efficiency_2d_unit,*) '# 1=rad_sigma 2=dp/p 3=n/tot_nabs 4=n 5=tot_nabs'
      do i=1,numeff
        do k=1,numeffdpop
          write(efficiency_2d_unit,'(4(1x,e15.7),1(1x,I5))') rsig(i), dpopbins(k),neff2d(i,k)/real(n_tot_absorbed,fPrec), &
                neff2d(i,k), n_tot_absorbed
        end do
      end do
    else
      write(lout,*) 'NO PARTICLE ABSORBED'
    end if
    close(efficiency_2d_unit)
#ifdef HDF5
  end if
#endif

! Write collimation summary file
#ifdef HDF5
  if(h5_useForCOLL) then
    allocate(fldHdf(7))
    fldHdf(1) = h5_dataField(name="ICOLL",    type=h5_typeInt)
    fldHdf(2) = h5_dataField(name="COLLNAME", type=h5_typeChar, size=mNameLen)
    fldHdf(3) = h5_dataField(name="NIMP",     type=h5_typeInt)
    fldHdf(4) = h5_dataField(name="NABS",     type=h5_typeInt)
    fldHdf(5) = h5_dataField(name="IMP_AV",   type=h5_typeReal)
    fldHdf(6) = h5_dataField(name="IMP_SIG",  type=h5_typeReal)
    fldHdf(7) = h5_dataField(name="LENGTH",   type=h5_typeReal)
    call h5_createFormat("collSummary", fldHdf, fmtHdf)
    call h5_createDataSet("coll_summary", h5_collID, fmtHdf, setHdf)
    ! There is a lot of overhead in writing line by line, but this is a small log file anyway.
    do i=1, cdb_nColl
      if(cdb_cLength(i) > zero .and. coll_found(i)) then
        call h5_prepareWrite(setHdf, 1)
        call h5_writeData(setHdf, 1, 1, i)
        call h5_writeData(setHdf, 2, 1, cdb_cNameUC(i))
        call h5_writeData(setHdf, 3, 1, cn_impact(i))
        call h5_writeData(setHdf, 4, 1, cn_absorbed(i))
        call h5_writeData(setHdf, 5, 1, caverage(i))
        call h5_writeData(setHdf, 6, 1, csigma(i))
        call h5_writeData(setHdf, 7, 1, cdb_cLength(i))
        call h5_finaliseWrite(setHdf)
      end if
    end do
    deallocate(fldHdf)
  else
#endif
    call f_requestUnit('coll_summary.dat', coll_summary_unit)
    open(unit=coll_summary_unit, file='coll_summary.dat') !was 50
    write(coll_summary_unit,*) '# 1=icoll 2=collname 3=nimp 4=nabs 5=imp_av 6=imp_sig 7=length'
    do icoll = 1, cdb_nColl
      if(cdb_cLength(icoll) > zero .and. coll_found(icoll)) then
        write(coll_summary_unit,'(i4,1x,a,2(1x,i5),2(1x,e15.7),3x,f4.1)') icoll, cdb_cNameUC(icoll), cn_impact(icoll), &
          cn_absorbed(icoll), caverage(icoll), csigma(icoll),cdb_cLength(icoll)
      end if
    end do
    close(coll_summary_unit)
#ifdef HDF5
  end if
#endif

#ifdef ROOT
  if(root_flag .and. root_Collimation.eq.1) then
    do icoll = 1, cdb_nColl
      if(cdb_cLength(icoll).gt.zero) then
        call CollimatorLossRootWrite(icoll, cdb_cNameUC(icoll), len(cdb_cNameUC(icoll)), cn_impact(icoll), cn_absorbed(icoll), &
          caverage(icoll), csigma(icoll), cdb_cLength(icoll))
      end if
    end do
  end if

  ! flush the root file
!  call SixTrackRootWrite()
#endif

end subroutine collimate_end_sample

!>
!! collimate_exit()
!! This routine is called once at the end of the simulation and
!! can be used to do any final postrocessing and/or file saving.
!<
subroutine collimate_exit()

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da

  implicit none

  integer :: i,j

  close(outlun)
  close(collgaps_unit)

  if(dowritetracks) then
    if(.not. cern) close(tracks2_unit)
#ifdef HDF5
    if(.not. cern .and. h5_writeTracks2) call h5tr2_finalise
#endif
    if(name_sel(1:3).eq.'COL') close(RHIClosses_unit)
  endif

  if(do_select) then
    close(coll_ellipse_unit)
  endif

  if(dowrite_impact) then
    close(all_impacts_unit)
    close(all_absorptions_unit)
    close(FLUKA_impacts_unit)
    close(FLUKA_impacts_all_unit)
    close(coll_scatter_unit)
    close(FirstImpacts_unit)
  endif

  call f_requestUnit('amplitude.dat', amplitude_unit)
  call f_requestUnit('amplitude2.dat', amplitude2_unit)
  call f_requestUnit('betafunctions.dat', betafunctions_unit)
  open(unit=amplitude_unit, file='amplitude.dat') !was 56
  open(unit=amplitude2_unit, file='amplitude2.dat') !was 51
  open(unit=betafunctions_unit, file='betafunctions.dat') !was 57

  if(dowrite_amplitude) then
    write(amplitude_unit,"(a)")                                         &
      "# 1=ielem 2=name 3=s 4=AX_AV 5=AX_RMS 6=AY_AV 7=AY_RMS "//       &
      "8=alphax 9=alphay 10=betax 11=betay 12=orbitx "//                &
      "13=orbity 14=tdispx 15=tdispy 16=xbob 17=ybob 18=xpbob 19=ypbob"

    do i=1,iu
       write(amplitude_unit,'(i4, (1x,a16), 17(1x,e20.13))')             &!hr08
      &i, ename(i), sampl(i),                                            &!hr08
      &sum_ax(i)/real(max(nampl(i),1),fPrec),                            &!hr08
      &sqrt(abs((sqsum_ax(i)/real(max(nampl(i),1),fPrec))-               &!hr08
      &(sum_ax(i)/real(max(nampl(i),1),fPrec))**2)),                     &!hr08
      &sum_ay(i)/real(max(nampl(i),1),fPrec),                            &!hr08
      &sqrt(abs((sqsum_ay(i)/real(max(nampl(i),1),fPrec))-               &!hr08
      &(sum_ay(i)/real(max(nampl(i),1),fPrec))**2)),                     &!hr08
      &talphax(i), talphay(i),                                           &!hr08
      &tbetax(i), tbetay(i), torbx(i), torby(i),                         &!hr08
      &tdispx(i), tdispy(i),                                             &!hr08
      &xbob(i),ybob(i),xpbob(i),ypbob(i)                                  !hr08
    end do

    write(amplitude2_unit,"(a)") "# 1=ielem 2=name 3=s 4=ORBITX 5=orbity 6=tdispx 7=tdispy 8=xbob 9=ybob 10=xpbob 11=ypbob"

    do i=1,iu
      write(amplitude2_unit,'(i4, (1x,a16), 9(1x,e15.7))') i, ename(i), sampl(i), torbx(i), torby(i), tdispx(i), tdispy(i), &
            xbob(i), ybob(i), xpbob(i), ypbob(i)
    end do

    write(betafunctions_unit,"(a)") "# 1=ielem 2=name       3=s             4=TBETAX(m)     5=TBETAY(m)     6=TORBX(mm)"// &
      "    7=TORBY(mm) 8=TORBXP(mrad)   9=TORBYP(mrad)  10=TDISPX(m)  11=MUX()    12=MUY()"


    do i=1,iu
!     RB: added printout of closed orbit and angle
      write(betafunctions_unit,'(i5, (1x,a16), 10(1x,e15.7))') i, ename(i), sampl(i), tbetax(i), tbetay(i), torbx(i), torby(i), &
        torbxp(i), torbyp(i), tdispx(i), mux(i), muy(i)
    end do
  endif

  close(amplitude_unit)
  close(amplitude2_unit)
  close(betafunctions_unit)

!GRD
!      DO J=1,iu
!        DO I=1,numl
!        xaveragesumoverturns(j)  = xaverage(j,i)
!     &                             + xaverage(j,MAX((i-1),1))
!        yaveragesumoverturns(j)  = yaverage(j,i)
!     &                             + yaverage(j,MAX((i-1),1))
!        xpaveragesumoverturns(j) = xpaverage(j,i)
!     &                             + xpaverage(j,MAX((i-1),1))
!        ypaveragesumoverturns(j) = ypaverage(j,i)
!     &                             + ypaverage(j,MAX((i-1),1))
!        END DO
!        xclosedorbitcheck(j)=(xaveragesumoverturns(j)
!     &                        +xaverage(j,numl))/(2*numl)
!        yclosedorbitcheck(j)=(yaveragesumoverturns(j)
!     &                        +yaverage(j,numl))/(2*numl)
!        xpclosedorbitcheck(j)=(xpaveragesumoverturns(j)
!     &                        +xpaverage(j,numl))/(2*numl)
!        ypclosedorbitcheck(j)=(ypaveragesumoverturns(j)
!     &                        +ypaverage(j,numl))/(2*numl)
!      END DO
!
!      OPEN(unit=99, file='xchecking.dat')
!      WRITE(99,*) '# 1=s 2=x 3=xp 4=y 5=yp'
!      DO J=1,iu
!      WRITE(99,'(i, 5(1x,e15.7))')
!     &     j, SAMPL(j),
!     &     xclosedorbitcheck(j), xpclosedorbitcheck(j),
!     &     yclosedorbitcheck(j), ypclosedorbitcheck(j)
!      END DO
!      CLOSE(99)
!GRD
!GRD WE CAN ALSO MAKE AN ORBIT CHECKING
!GRD

  call f_requestUnit('orbitchecking.dat', orbitchecking_unit)
  open(unit=orbitchecking_unit, file='orbitchecking.dat') !was 99
  write(orbitchecking_unit,*) '# 1=s 2=torbitx 3=torbity'

  do j=1,iu
    write(orbitchecking_unit,'(i5, 3(1x,e15.7))') j, sampl(j),torbx(j), torby(j)
  end do

  close(orbitchecking_unit)
  close(CollPositions_unit)

#ifdef G4COLLIMAT
  call g4_terminate()
#endif

end subroutine collimate_exit

!>
!! This routine is called at the start of each tracking turn
!<
subroutine collimate_start_turn(n)

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da

  implicit none

  integer, intent(in) :: n

  iturn=n
  totals=zero !This keeps track of the s position of the current element, which is also done by cadcum
end subroutine collimate_start_turn

!>
!! This routine is called at the start of every element
!<
subroutine collimate_start_element(i)

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da

  implicit none

  integer, intent(in) :: i
  integer j

  ie=i
!++  For absorbed particles set all coordinates to zero. Also
!++  include very large offsets, let's say above 100mm or
!++  100mrad.
  do j = 1, napx
    if( (part_abs_pos(j).ne.0 .and. part_abs_turn(j).ne.0) .or.&
 &  xv1(j).gt.c1e2 .or. yv1(j).gt.c1e2 .or. xv2(j).gt.c1e2 .or. yv2(j).gt.c1e2) then
      xv1(j) = zero
      yv1(j) = zero
      xv2(j) = zero
      yv2(j) = zero
      ejv(j)  = myenom
      sigmv(j)= zero
      part_abs_pos(j)=ie
      part_abs_turn(j)=iturn
      secondary(j) = 0
      tertiary(j)  = 0
      other(j)     = 0
      scatterhit(j)= 0
      nabs_type(j) = 0
    end if
  end do

!GRD SAVE COORDINATES OF PARTICLE 1 TO CHECK ORBIT
  if(firstrun) then
    xbob(ie)=xv1(1)
    ybob(ie)=xv2(1)
    xpbob(ie)=yv1(1)
    ypbob(ie)=yv2(1)
  end if

!++  Here comes sixtrack stuff
  if(ic(i).le.nblo) then
    do jb=1,mel(ic(i))
      myix=mtyp(ic(i),jb)
    end do
  else
    myix=ic(i)-nblo
  end if

!++  Make sure we go into collimation routine for any definition
!++  of collimator element, relying on element name instead.
  if (                                                          &
!GRD HERE ARE SOME CHANGES TO MAKE RHIC TRAKING AVAILABLE
!APRIL2005
     &(bez(myix)(1:3).eq.'TCP'.or.bez(myix)(1:3).eq.'tcp') .or.         &
     &(bez(myix)(1:3).eq.'TCS'.or.bez(myix)(1:3).eq.'tcs') .or.         &
!UPGRADE January 2005
     &(bez(myix)(1:3).eq.'TCL'.or.bez(myix)(1:3).eq.'tcl') .or.         &
     &(bez(myix)(1:3).eq.'TCT'.or.bez(myix)(1:3).eq.'tct') .or.         &
     &(bez(myix)(1:3).eq.'TCD'.or.bez(myix)(1:3).eq.'tcd') .or.         &
     &(bez(myix)(1:3).eq.'TDI'.or.bez(myix)(1:3).eq.'tdi') .or.         &
! UPGRADE MAI 2006 -> TOTEM
     &(bez(myix)(1:3).eq.'TCX'.or.bez(myix)(1:3).eq.'tcx') .or.         &
! TW 04/2008 adding TCRYO
     &(bez(myix)(1:3).eq.'TCR'.or.bez(myix)(1:3).eq.'tcr') .or.         &
!RHIC
     &(bez(myix)(1:3).eq.'COL'.or.bez(myix)(1:3).eq.'col') ) then

    myktrack = 1
  else
    myktrack = ktrack(i)
  endif

end subroutine collimate_start_element

!>
!! collimate_end_element()
!! This routine is called at the end of every element
!<
subroutine collimate_end_element

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da

  implicit none

  integer j

#ifdef HDF5
  integer hdfturn,hdfpid,hdftyp
  real(kind=fPrec) hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdfs
#endif

  if(firstrun) then
    if(rselect.gt.0 .and. rselect.lt.65) then
      do j = 1, napx
        xj  = (xv1(j)-torbx(ie)) /c1e3
        xpj = (yv1(j)-torbxp(ie))/c1e3
        yj  = (xv2(j)-torby(ie)) /c1e3
        ypj = (yv2(j)-torbyp(ie))/c1e3
        pj  = ejv(j)/c1e3

        if(iturn.eq.1.and.j.eq.1) then
          sum_ax(ie) = zero
          sum_ay(ie) = zero
        endif

        if(tbetax(ie).gt.zero) then
          gammax = (one + talphax(ie)**2)/tbetax(ie)
          gammay = (one + talphay(ie)**2)/tbetay(ie)
        else
          gammax = (one + talphax(ie-1)**2)/tbetax(ie-1)
          gammay = (one + talphay(ie-1)**2)/tbetay(ie-1)
        endif

        if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
          if(tbetax(ie).gt.0.) then
            nspx = sqrt(abs( gammax*(xj)**2 + two*talphax(ie)*xj*xpj +   tbetax(ie)*xpj**2 )/myemitx0_collgap)
            nspy = sqrt(abs( gammay*(yj)**2 + two*talphay(ie)*yj*ypj +   tbetay(ie)*ypj**2 )/myemity0_collgap)
          else
            nspx = sqrt(abs( gammax*(xj)**2 + two*talphax(ie-1)*xj*xpj + tbetax(ie-1)*xpj**2 )/myemitx0_collgap)
            nspy = sqrt(abs( gammay*(yj)**2 + two*talphay(ie-1)*yj*ypj + tbetay(ie-1)*ypj**2 )/myemity0_collgap)
          end if

          sum_ax(ie)   = sum_ax(ie) + nspx
          sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
          sum_ay(ie)   = sum_ay(ie) + nspy
          sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
          nampl(ie)    = nampl(ie) + 1
        else
          nspx = zero
          nspy = zero
        end if

        sampl(ie) = totals
        ename(ie) = bez(myix)(1:mNameLen)
      end do
    end if
  end if

!GRD THIS LOOP MUST NOT BE WRITTEN INTO THE "IF(FIRSTRUN)" LOOP !!!!
  if (dowritetracks) then
    do j = 1, napx
      xj     = (xv1(j)-torbx(ie)) /c1e3
      xpj    = (yv1(j)-torbxp(ie))/c1e3
      yj     = (xv2(j)-torby(ie)) /c1e3
      ypj    = (yv2(j)-torbyp(ie))/c1e3

      arcdx = 2.5_fPrec
      arcbetax = c180e0

      if (xj.le.zero) then
        xdisp = xj + (pj-myenom)/myenom * arcdx* sqrt(tbetax(ie)/arcbetax)
      else
        xdisp = xj - (pj-myenom)/myenom * arcdx* sqrt(tbetax(ie)/arcbetax)
      end if

      xndisp = xj

      nspxd = sqrt(abs(gammax*xdisp**2 +  two*talphax(ie)*xdisp*xpj +  tbetax(ie)*xpj**2)/myemitx0_collgap)
      nspx  = sqrt(abs(gammax*xndisp**2 + two*talphax(ie)*xndisp*xpj + tbetax(ie)*xpj**2)/myemitx0_collgap)
      nspy  = sqrt(abs(gammay*yj**2 +     two*talphay(ie)*yj*ypj +     tbetay(ie)*ypj**2)/myemity0_collgap)

      if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then

!GRD HERE WE APPLY THE SAME KIND OF CUT THAN THE SIGSECUT PARAMETER
         if((secondary(j) .eq. 1 .or. &
             tertiary(j)  .eq. 2 .or. &
             other(j)     .eq. 4 .or. &
             scatterhit(j).eq. 8       ) .and. &
             (xv1(j).lt.99.0_fPrec .and. xv2(j).lt.99.0_fPrec) .and. &
             ((((xv1(j)*c1m3)**2 / (tbetax(ie)*myemitx0_collgap)) .ge. real(sigsecut2,fPrec)).or. &
             (((xv2(j)*c1m3)**2  / (tbetay(ie)*myemity0_collgap)) .ge. real(sigsecut2,fPrec)).or. &
             (((xv1(j)*c1m3)**2  / (tbetax(ie)*myemitx0_collgap)) + &
             ((xv2(j)*c1m3)**2  /  (tbetay(ie)*myemity0_collgap)) .ge. sigsecut3)) ) &
             then

          xj  = (xv1(j)-torbx(ie)) /c1e3
          xpj = (yv1(j)-torbxp(ie))/c1e3
          yj  = (xv2(j)-torby(ie)) /c1e3
          ypj = (yv2(j)-torbyp(ie))/c1e3
#ifdef HDF5
          if(h5_writeTracks2) then
            hdfpid=ipart(j)
            hdfturn=iturn
            hdfs=sampl(ie)
            hdfx=xv1(j)
            hdfxp=yv1(j)
            hdfy=xv2(j)
            hdfyp=yv2(j)
            hdfdee=(ejv(j)-myenom)/myenom
            hdftyp=secondary(j)+tertiary(j)+other(j)+scatterhit(j)
            call h5tr2_writeLine(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdftyp)
          else
#endif
            write(tracks2_unit,'(1x,i8,1x,i4,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)') ipart(j), iturn, sampl(ie), &
              xv1(j), yv1(j), xv2(j), yv2(j), (ejv(j)-myenom)/myenom, secondary(j)+tertiary(j)+other(j)+scatterhit(j)
#ifdef HDF5
          end if
#endif
        end if
      end if
    end do
  end if !!JUNE2005 here I close the "if(dowritetracks)" outside of the firstrun flag

end subroutine collimate_end_element

!>
!! collimate_end_turn()
!! This routine is called at the end of every turn
!<
subroutine collimate_end_turn

  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use crcoall

#ifdef ROOT
  use root_output
#endif

  implicit none

  integer j
  integer napx_pre

#ifdef HDF5
  ! For tracks2
  integer hdfturn,hdfpid,hdftyp
  real(kind=fPrec) hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdfs
  ! For other output
  type(h5_dataField), allocatable :: fldHdf(:)
  integer fmtHdf, setHdf
#endif

!__________________________________________________________________
!++  Now do analysis at selected elements...

!++  Save twiss functions of present element
  ax0  = talphax(ie)
  bx0  = tbetax(ie)
  mux0 = mux(ie)
  ay0  = talphay(ie)
  by0  = tbetay(ie)
  muy0 = muy(ie)

!GRD GET THE COORDINATES OF THE PARTICLES AT THE IEth ELEMENT:
  do j = 1,napx
    xgrd(j)  = xv1(j)
    xpgrd(j) = yv1(j)
    ygrd(j)  = xv2(j)
    ypgrd(j) = yv2(j)

    xineff(j)  = xv1(j) - torbx (ie)
    xpineff(j) = yv1(j) - torbxp(ie)
    yineff(j)  = xv2(j) - torby (ie)
    ypineff(j) = yv2(j) - torbyp(ie)

    pgrd(j)  = ejv(j)
    ejfvgrd(j) = ejfv(j)
    sigmvgrd(j) = sigmv(j)
    rvvgrd(j) = rvv(j)
    dpsvgrd(j) = dpsv(j)
    oidpsvgrd(j) = oidpsv(j)
    dpsv1grd(j) = dpsv1(j)

!GRD IMPORTANT: ALL PARTICLES ABSORBED ARE CONSIDERED TO BE LOST,
!GRD SO WE GIVE THEM A LARGE OFFSET
    if(part_abs_pos(j).ne.0 .and. part_abs_turn(j).ne.0) then
      xgrd(j) = 99.5_fPrec
      ygrd(j) = 99.5_fPrec
    end if
  end do

!++  For LAST ELEMENT in the ring calculate the number of surviving
!++  particles and save into file versus turn number
  if(ie.eq.iu) then
    nsurvive = 0

    do j = 1, napx
      if (xgrd(j).lt.99.0_fPrec .and. ygrd(j).lt.99.0_fPrec) then
        nsurvive = nsurvive + 1
      end if
    end do

#ifdef HDF5
    if(h5_useForCOLL) then
      call h5_prepareWrite(coll_hdf5_survival, 1)
      call h5_writeData(coll_hdf5_survival, 1, 1, iturn)
      call h5_writeData(coll_hdf5_survival, 2, 1, nsurvive)
      call h5_finaliseWrite(coll_hdf5_survival)
    else
#endif
      write(survival_unit,'(2i7)') iturn, nsurvive
#ifdef HDF5
    end if
#endif

#ifdef ROOT
    if(root_flag .and. root_Collimation.eq.1) then
      call SurvivalRootWrite(iturn, nsurvive)
    end if
#endif

    if (iturn.eq.numl) then
      nsurvive_end = nsurvive_end + nsurvive
    end if
  end if

!=======================================================================
!++  Do collimation analysis at element 20 ("zero" turn) or LAST
!++  ring element.

!++  If selecting, look at number of scattered particles at selected
!++  collimator. For the "zero" turn consider the information at element
!++  20 (before collimation), otherwise take information at last ring
!++  element.
  if (do_coll .and. (  (iturn.eq.1 .and. ie.eq.20) .or. (ie.eq.iu) ) ) then

!++  Calculate gammas
!------------------------------------------------------------------------
    gammax = (1 + talphax(ie)**2)/tbetax(ie)
    gammay = (1 + talphay(ie)**2)/tbetay(ie)

!________________________________________________________________________
!++  Loop over all particles.
    do j = 1, napx
!
!------------------------------------------------------------------------
!++  Save initial distribution of particles that were scattered on
!++  the first turn at the selected primary collimator
!
!            IF (DOWRITE_DIST .AND. DO_SELECT .AND. ITURN.EQ.1 .AND.
!     &          PART_SELECT(j).EQ.1) THEN
!              WRITE(987,'(4(1X,E15.7))') X00(J), XP00(J),
!     &                                        Y00(J), YP00(J)
!            ENDIF
!------------------------------------------------------------------------
!++  Do the binning in amplitude, only considering particles that were
!++  not absorbed before.

      if (xgrd(j).lt.99.0_fPrec .and. ygrd(j) .lt.99.0_fPrec .and. (part_select(j).eq.1 .or. ie.eq.20)) then

!++  Normalized amplitudes are calculated

!++  Allow to apply some dispersive offset. Take arc dispersion (2m) and
!++  normalize with arc beta_x function (180m).
        arcdx    = 2.5_fPrec
        arcbetax = c180e0
        xdisp = abs(xgrd(j)*c1m3) + abs((pgrd(j)-myenom)/myenom)*arcdx * sqrt(tbetax(ie)/arcbetax)
        nspx = sqrt(abs(gammax*xdisp**2 +two*talphax(ie)*xdisp*(xpgrd(j)*c1m3)+tbetax(ie)*(xpgrd(j)*c1m3)**2 )/myemitx0_collgap)
        nspy = sqrt(abs(gammay*(ygrd(j)*c1m3)**2 + two*talphay(ie)*(ygrd(j)*c1m3*ypgrd(j)*c1m3)+ tbetay(ie)*(ypgrd(j)*c1m3)**2 )&
   &           /myemity0_collgap)

!++  Populate the efficiency arrays at the end of each turn...
! Modified by M.Fiascaris, July 2016
        if(ie.eq.iu) then
          do ieff = 1, numeff
            if(counted_r(j,ieff).eq.0 .and. sqrt( &
            &((xineff(j)*c1m3)**2 + (talphax(ie)*xineff(j)*c1m3 + tbetax(ie)*xpineff(j)*c1m3)**2)/(tbetax(ie)*myemitx0_collgap)+&
            &((yineff(j)*c1m3)**2 + (talphay(ie)*yineff(j)*c1m3 + tbetay(ie)*ypineff(j)*c1m3)**2)/(tbetay(ie)*myemity0_collgap))&
            &.ge.rsig(ieff)) then
              neff(ieff) = neff(ieff)+one
              counted_r(j,ieff)=1
            end if

!++ 2D eff
            do ieffdpop = 1, numeffdpop
              if(counted2d(j,ieff,ieffdpop).eq.0 .and.abs((ejv(j)-myenom)/myenom).ge.dpopbins(ieffdpop)) then
                neff2d(ieff,ieffdpop) = neff2d(ieff,ieffdpop)+one
                counted2d(j,ieff,ieffdpop)=1
              end if
            end do

            if(counted_x(j,ieff).eq.0 .and.sqrt(((xineff(j)*c1m3)**2 + &
            &(talphax(ie)*xineff(j)*c1m3 + tbetax(ie)*xpineff(j)*c1m3)**2)/(tbetax(ie)*myemitx0_collgap)).ge.rsig(ieff)) then
              neffx(ieff) = neffx(ieff) + one
              counted_x(j,ieff)=1
            end if

            if(counted_y(j,ieff).eq.0 .and. &
            &sqrt(((yineff(j)*c1m3)**2 + (talphay(ie)*yineff(j)*c1m3 + tbetay(ie)*ypineff(j)*c1m3)**2)/ &
            &tbetay(ie)*myemity0_collgap).ge.rsig(ieff)) then
              neffy(ieff) = neffy(ieff) + one
              counted_y(j,ieff)=1
            end if
          end do !do ieff = 1, numeff

          do ieffdpop = 1, numeffdpop
            if(counteddpop(j,ieffdpop).eq.0) then
              dpopmin = zero
              mydpop = abs((ejv(j)-myenom)/myenom)
              if(ieffdpop.gt.1) dpopmin = dpopbins(ieffdpop-1)

              dpopmax = dpopbins(ieffdpop)
              if(mydpop.ge.dpopmin .and. mydpop.lt.mydpop) then
                npartdpop(ieffdpop)=npartdpop(ieffdpop)+1
              end if
            end if

            if(counteddpop(j,ieffdpop).eq.0 .and.abs((ejv(j)-myenom)/myenom).ge.dpopbins(ieffdpop)) then
              neffdpop(ieffdpop) = neffdpop(ieffdpop)+one
              counteddpop(j,ieffdpop)=1
            end if
          end do !do ieffdpop = 1, numeffdpop
        end if !if(ie.eq.iu) then

!++  Do an emittance drift
        driftx = driftsx*sqrt(tbetax(ie)*myemitx0_collgap)
        drifty = driftsy*sqrt(tbetay(ie)*myemity0_collgap)

        if(ie.eq.iu) then
          dnormx = driftx / sqrt(tbetax(ie)*myemitx0_collgap)
          dnormy = drifty / sqrt(tbetay(ie)*myemity0_collgap)
          xnorm  = (xgrd(j)*c1m3) / sqrt(tbetax(ie)*myemitx0_collgap)
          xpnorm = (talphax(ie)*(xgrd(j)*c1m3)+ tbetax(ie)*(xpgrd(j)*c1m3)) / sqrt(tbetax(ie)*myemitx0_collgap)
          xangle = atan2_mb(xnorm,xpnorm)
          xnorm  = xnorm  + dnormx*sin_mb(xangle)
          xpnorm = xpnorm + dnormx*cos_mb(xangle)
          xgrd(j)  = c1e3 * (xnorm * sqrt(tbetax(ie)*myemitx0_collgap))
          xpgrd(j) = c1e3 * ((xpnorm*sqrt(tbetax(ie)*myemitx0_collgap)-talphax(ie)*xgrd(j)*c1m3)/tbetax(ie))

          ynorm  = (ygrd(j)*c1m3)/ sqrt(tbetay(ie)*myemity0_collgap)
          ypnorm = (talphay(ie)*(ygrd(j)*c1m3)+tbetay(ie)*(ypgrd(j)*c1m3)) / sqrt(tbetay(ie)*myemity0_collgap)
          yangle = atan2_mb(ynorm,ypnorm)
          ynorm  = ynorm  + dnormy*sin_mb(yangle)
          ypnorm = ypnorm + dnormy*cos_mb(yangle)
          ygrd(j)  = c1e3 * (ynorm * sqrt(tbetay(ie)*myemity0_collgap))
          ypgrd(j) = c1e3 * ((ypnorm*sqrt(tbetay(ie)*myemity0_collgap)-talphay(ie)*ygrd(j)*c1m3)/tbetay(ie))
        end if

!------------------------------------------------------------------------
!++  End of check for selection flag and absorption
      end if
!++  End of do loop over particles
    end do
!_________________________________________________________________
!++  End of collimation efficiency analysis for selected particles
  end if

!------------------------------------------------------------------
!++  For LAST ELEMENT in the ring compact the arrays by moving all
!++  lost particles to the end of the array.
  napx_pre = napx
  if(ie.eq.iu) then
    imov = 0
    do j = 1, napx
      if(xgrd(j).lt.99.0_fPrec .and. ygrd(j).lt.99.0_fPrec) then
        llostp(j) = .false.
      else
        llostp(j) = .true.
      end if
    end do

    ! Move the lost particles to the end of the arrays
    call shuffleLostParticles

    write(lout,"(3(a,i0))") "COLL> Compacted the particle distributions: ",napx_pre," --> ",napx,", turn = ",iturn
    flush(lout)

! napx gets updated by shuffleLostParticles
!    napx = imov
  endif

  ! Write final distribution
  if(dowrite_dist .and. ie == iu .and. iturn == numl) then
#ifdef HDF5
    if(h5_useForCOLL) then
      allocate(fldHdf(6))
      fldHdf(1) = h5_dataField(name="X",  type=h5_typeReal)
      fldHdf(2) = h5_dataField(name="XP", type=h5_typeReal)
      fldHdf(3) = h5_dataField(name="Y",  type=h5_typeReal)
      fldHdf(4) = h5_dataField(name="YP", type=h5_typeReal)
      fldHdf(5) = h5_dataField(name="Z",  type=h5_typeReal)
      fldHdf(6) = h5_dataField(name="E",  type=h5_typeReal)
      call h5_createFormat("collDistN", fldHdf, fmtHdf)
      call h5_createDataSet("distn", h5_collID, fmtHdf, setHdf, napx)
      call h5_prepareWrite(setHdf, napx)
      call h5_writeData(setHdf, 1, napx, (xgrd(1:napx) -torbx(1)) /c1e3)
      call h5_writeData(setHdf, 2, napx, (xpgrd(1:napx)-torbxp(1))/c1e3)
      call h5_writeData(setHdf, 3, napx, (ygrd(1:napx) -torby(1)) /c1e3)
      call h5_writeData(setHdf, 4, napx, (ypgrd(1:napx)-torbyp(1))/c1e3)
      call h5_writeData(setHdf, 5, napx, sigmvgrd(1:napx))
      call h5_writeData(setHdf, 6, napx, ejfvgrd(1:napx))
      call h5_finaliseWrite(setHdf)
      deallocate(fldHdf)
    else
#endif
      call f_requestUnit('distn.dat', distn_unit)
      open(unit=distn_unit, file='distn.dat') !was 9998
      write(distn_unit,*) '# 1=x 2=xp 3=y 4=yp 5=z 6 =E'
      do j = 1, napx
        write(distn_unit,'(6(1X,E23.15))') (xgrd(j)-torbx(1))/c1e3, (xpgrd(j)-torbxp(1))/c1e3, (ygrd(j)-torby(1))/c1e3, &
          (ypgrd(j)-torbyp(1))/c1e3, sigmvgrd(j), ejfvgrd(j)
      end do
      close(distn_unit)
#ifdef HDF5
    end if
#endif
  end if

!GRD NOW ONE HAS TO COPY BACK THE NEW DISTRIBUTION TO ITS "ORIGINAL NAME"
!GRD AT THE END OF EACH TURN
  if(ie.eq.iu) then
    do j = 1,napx
      xv1(j) = xgrd(j)
      yv1(j) = xpgrd(j)
      xv2(j) = ygrd(j)
      yv2(j) = ypgrd(j)
      ejv(j)  = pgrd(j)
      ejfv(j)   = ejfvgrd(j)
      sigmv(j)  = sigmvgrd(j)
      rvv(j)    = rvvgrd(j)
      dpsv(j)   = dpsvgrd(j)
      oidpsv(j) = oidpsvgrd(j)
      dpsv1(j)  = dpsv1grd(j)
    end do
  end if

  if(firstrun) then
    if(rselect.gt.0 .and. rselect.lt.65) then
      do j = 1, napx
        xj  = (xv1(j)-torbx(ie)) /c1e3
        xpj = (yv1(j)-torbxp(ie))/c1e3
        yj  = (xv2(j)-torby(ie)) /c1e3
        ypj = (yv2(j)-torbyp(ie))/c1e3
        pj  = ejv(j)/c1e3

        if(iturn.eq.1.and.j.eq.1) then
          sum_ax(ie)=zero
          sum_ay(ie)=zero
        end if

        if(tbetax(ie).gt.0.) then
          gammax = (one + talphax(ie)**2)/tbetax(ie)
          gammay = (one + talphay(ie)**2)/tbetay(ie)
        else
          gammax = (one + talphax(ie-1)**2)/tbetax(ie-1)
          gammay = (one + talphay(ie-1)**2)/tbetay(ie-1)
        end if

        if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
          if(tbetax(ie).gt.0.) then
            nspx = sqrt(abs( gammax*(xj)**2 + two*talphax(ie)*xj*xpj + tbetax(ie)*xpj**2 )/myemitx0_collgap)
            nspy = sqrt(abs( gammay*(yj)**2 + two*talphay(ie)*yj*ypj + tbetay(ie)*ypj**2 )/myemity0_collgap)
          else
            nspx = sqrt(abs( gammax*(xj)**2 + two*talphax(ie-1)*xj*xpj +tbetax(ie-1)*xpj**2 )/myemitx0_collgap)
            nspy = sqrt(abs( gammay*(yj)**2 + two*talphay(ie-1)*yj*ypj +tbetay(ie-1)*ypj**2 )/myemity0_collgap)
          end if

          sum_ax(ie)   = sum_ax(ie) + nspx
          sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
          sum_ay(ie)   = sum_ay(ie) + nspy
          sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
          nampl(ie)    = nampl(ie) + 1
          sampl(ie)    = totals
          ename(ie)    = bez(myix)(1:mNameLen)
        else
          nspx = zero
          nspy = zero
        end if !if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
      end do !do j = 1, napx
    end if !if(rselect.gt.0 .and. rselect.lt.65) then
  end if !if(firstrun) then

!GRD THIS LOOP MUST NOT BE WRITTEN INTO THE "IF(FIRSTRUN)" LOOP !!!!
  if(dowritetracks) then
    do j = 1, napx
      xj    = (xv1(j)-torbx(ie))/c1e3
      xpj   = (yv1(j)-torbxp(ie))/c1e3
      yj    = (xv2(j)-torby(ie))/c1e3
      ypj   = (yv2(j)-torbyp(ie))/c1e3
      arcdx = 2.5_fPrec
      arcbetax = c180e0

      if(xj.le.0.) then
        xdisp = xj + (pj-myenom)/myenom * arcdx * sqrt(tbetax(ie)/arcbetax)
      else
        xdisp = xj - (pj-myenom)/myenom * arcdx * sqrt(tbetax(ie)/arcbetax)
      end if

      xndisp = xj
      nspxd  = sqrt(abs(gammax*xdisp**2  + two*talphax(ie)*xdisp*xpj  + tbetax(ie)*xpj**2)/myemitx0_collgap)
      nspx   = sqrt(abs(gammax*xndisp**2 + two*talphax(ie)*xndisp*xpj + tbetax(ie)*xpj**2)/myemitx0_collgap)
      nspy   = sqrt(abs( gammay*yj**2    + two*talphay(ie)*yj*ypj     + tbetay(ie)*ypj**2)/myemity0_collgap)

      if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
!GRD HERE WE APPLY THE SAME KIND OF CUT THAN THE SIGSECUT PARAMETER
        if((secondary(j) .eq. 1 .or. &
            tertiary(j)  .eq. 2 .or. &
            other(j)     .eq. 4 .or. &
            scatterhit(j).eq. 8        ) .and. &
            (xv1(j).lt.99.0_fPrec .and. xv2(j).lt.99.0_fPrec) .and. &
            ((((xv1(j)*c1m3)**2 / (tbetax(ie)*myemitx0_collgap)) .ge. real(sigsecut2,fPrec)).or. &
            (((xv2(j)*c1m3)**2  / (tbetay(ie)*myemity0_collgap)) .ge. real(sigsecut2,fPrec)).or. &
            (((xv1(j)*c1m3)**2  / (tbetax(ie)*myemitx0_collgap)) + &
            ((xv2(j)*c1m3)**2  / (tbetay(ie)*myemity0_collgap)) .ge. sigsecut3)) ) &
            then

          xj     = (xv1(j)-torbx(ie))/c1e3
          xpj    = (yv1(j)-torbxp(ie))/c1e3
          yj     = (xv2(j)-torby(ie))/c1e3
          ypj    = (yv2(j)-torbyp(ie))/c1e3
#ifdef HDF5
          if(h5_writeTracks2) then
            hdfpid=ipart(j)
            hdfturn=iturn
            hdfs=sampl(ie)
            hdfx=xv1(j)
            hdfxp=yv1(j)
            hdfy=xv2(j)
            hdfyp=yv2(j)
            hdfdee=(ejv(j)-myenom)/myenom
            hdftyp=secondary(j)+tertiary(j)+other(j)+scatterhit(j)
            call h5tr2_writeLine(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdftyp)
          else
#endif
            write(tracks2_unit,'(1x,i8,1x,i4,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)') ipart(j),iturn,sampl(ie), &
              xv1(j),yv1(j),xv2(j),yv2(j),(ejv(j)-myenom)/myenom,secondary(j)+tertiary(j)+other(j)+scatterhit(j)
#ifdef HDF5
          end if
#endif
        end if !if ((secondary(j).eq.1.or.tertiary(j).eq.2.or.other(j).eq.4.or.scatterhit(j).eq.8
      end if !if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
    end do ! do j = 1, napx
  end if !if(dowritetracks) then
!=======================================================================
end subroutine collimate_end_turn

!>
!! "Merlin" scattering collimation configuration
!! This routine pre-calcuates some varibles for
!! the nuclear properties
!<
subroutine collimate_init_merlin()

  implicit none

  integer i

! compute the electron densnity and plasma energy for each material
  do i=1, nmat
    edens(i) = CalcElectronDensity(zatom(i),rho(i),anuc(i))
    pleng(i) = CalcPlasmaEnergy(edens(i))
  end do

end subroutine collimate_init_merlin

!>
!! K2 scattering collimation configuration
!<
subroutine collimate_init_k2()
!nothing currently
end subroutine collimate_init_k2

!>
!! subroutine collimate2(c_material, c_length, c_rotation,           &
!!-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----
!!----                                                                    -----
!!-----  NEW ROUTINES PROVIDED FOR THE COLLIMATION STUDIES VIA SIXTRACK   -----
!!-----                                                                   -----
!!-----          G. ROBERT-DEMOLAIZE, November 1st, 2004                  -----
!!-----                                                                   -----
!!-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----
!!++  Based on routines by JBJ. Changed by RA 2001.
!!GRD
!!GRD MODIFIED VERSION FOR COLLIMATION SYSTEM: G. ROBERT-DEMOLAIZE
!!GRD
!!
!!++  - Deleted all HBOOK stuff.
!!++  - Deleted optics routine and all parser routines.
!!++  - Replaced RANMAR call by RANLUX call
!!++  - Included RANLUX code from CERNLIB into source
!!++  - Changed dimensions from CGen(100,nmat) to CGen(200,nmat)
!!++  - Replaced FUNPRE with FUNLXP
!!++  - Replaced FUNRAN with FUNLUX
!!++  - Included all CERNLIB code into source: RANLUX, FUNLXP, FUNLUX,
!!++                                         FUNPCT, FUNLZ, RADAPT,
!!++                                           RGS56P
!!++    with additional entries:             RLUXIN, RLUXUT, RLUXAT,
!!++                                           RLUXGO
!!++
!!++  - Changed program so that Nev is total number of particles
!!++    (scattered and not-scattered)
!!++  - Added debug comments
!!++  - Put real dp/dx
!<
subroutine collimate2(c_material, c_length, c_rotation,           &
     &c_aperture, c_offset, c_tilt,x_in, xp_in, y_in,yp_in,p_in, s_in,  &
     &     np, enom,                                                    &
     &     lhit_pos, lhit_turn,                                         &
     &     part_abs_pos_local, part_abs_turn_local,                     &
     &     impact, indiv, lint, onesided,name,                          &
     &     flagsec, j_slices, nabs_type)

  use crcoall
  use parpro
  use mod_common, only : iexact
  implicit none

! BLOCK DBCOLLIM
! This block is common to collimaterhic and collimate2
! It is NOT compatible with block DBCOMMON, as some variable names overlap...


  logical onesided,hit
! integer nprim,filel,mat,nev,j,nabs,nhit,np,icoll,nabs_tmp
  integer nprim,j,nabs,nhit,np

  integer, allocatable :: lhit_pos(:) !(npart)
  integer, allocatable :: lhit_turn(:) !(npart)
  integer, allocatable :: part_abs_pos_local(:) !(npart)
  integer, allocatable :: part_abs_turn_local(:) !(npart)
  integer, allocatable :: name(:) !(npart)
  integer, allocatable :: nabs_type(:) !(npart)
!MAY2005

  real(kind=fPrec), allocatable :: x_in(:) !(npart)
  real(kind=fPrec), allocatable :: xp_in(:) !(npart)
  real(kind=fPrec), allocatable :: y_in(:) !(npart)
  real(kind=fPrec), allocatable :: yp_in(:) !(npart)
  real(kind=fPrec), allocatable :: p_in(:) !(npart)
  real(kind=fPrec), allocatable :: s_in(:) !(npart)
  real(kind=fPrec), allocatable :: indiv(:) !(npart)
  real(kind=fPrec), allocatable :: lint(:) !(npart)
  real(kind=fPrec), allocatable :: impact(:) !(npart)
  real(kind=fPrec) keeps,fracab,drift_length,mirror,tiltangle

  real(kind=fPrec) c_length    !length in m
  real(kind=fPrec) c_rotation  !rotation angle vs vertical in radian
  real(kind=fPrec) c_aperture  !aperture in m
  real(kind=fPrec) c_offset    !offset in m
  real(kind=fPrec) c_tilt(2)   !tilt in radian
  character(len=4) c_material  !material

  real(kind=fPrec) x00,z00,p,sp,s,enom

!AUGUST2006 Added ran_gauss for generation of pencil/     ------- TW
!           sheet beam distribution  (smear in x and y)
!

  real(kind=fPrec) x_flk,xp_flk,y_flk,yp_flk,zpj

  real(kind=fPrec) s_impact
  integer flagsec(npart)

!     SR, 18-08-2005: add temporary variable to write in FirstImpacts
!     the initial distribution of the impacting particles in the
!     collimator frame.
  real(kind=fPrec) xinn,xpinn,yinn,ypinn

!     SR, 29-08-2005: add the slice number to calculate the impact
!     location within the collimator.
!     j_slices = 1 for the a non sliced collimator!
  integer j_slices

  save
!      write(lout,*) 'In col2 ', c_material, c_length, c_aperture,       &
!     &c_offset, c_tilt, x_in, xp_in, y_in,p_in, np, enom
!=======================================================================
! Be=1 Al=2 Cu=3 W=4 Pb=5
! LHC uses:    Al, 0.2 m
!              Cu, 1.0 m

  if(c_material.eq.'BE') then
    mat = 1
  else if(c_material.eq.'AL') then
    mat = 2
  else if(c_material.eq.'CU') then
    mat = 3
  else if(c_material.eq.'W') then
    mat = 4
  else if(c_material.eq.'PB') then
    mat = 5
  else if(c_material.eq.'C') then
    mat = 6
  else if(c_material.eq.'C2') then
    mat = 7
  else if(c_material.eq.'MoGR') then
    mat = 8
  else if(c_material.eq.'CuCD') then
    mat = 9
  else if(c_material.eq.'Mo') then
    mat = 10
  else if(c_material.eq.'Glid') then
    mat = 11
  else if(c_material.eq.'Iner') then
    mat = 12
!02/2008 TW added vacuum and black absorber (was missing)
  else if(c_material.eq.'VA') then
    mat = nmat-1
  else if(c_material.eq.'BL') then
    mat = nmat
  else
    write(lout,*)
    write(lout,*) 'ERR>  In subroutine collimate2:'
    write(lout,*) 'ERR>  Material "', c_material, '" not found.'
    write(lout,*) 'ERR>  Check your CollDB! Stopping now.'
    call prror(-1)
  end if

  length  = c_length
  nev = np
  p0  = enom

!++  Initialize scattering processes
  call scatin(p0)

! EVENT LOOP,  initial distribution is here a flat distribution with
! xmin=x-, xmax=x+, etc. from the input file

  nhit    = 0
  fracab  = zero
  mirror  = one

!==> SLICE here

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  do j = 1, nev
! SR-GRD (04-08-2005):
!   Don't do scattering process for particles already absorbed
    if ( part_abs_pos_local(j) .ne. 0 .and. part_abs_turn_local(j) .ne. 0) goto 777

    impact(j) = -one
    lint(j)   = -one
    indiv(j)  = -one

    x   = x_in(j)
    xp  = xp_in(j)
    z   = y_in(j)
    zp  = yp_in(j)
    p   = p_in(j)
    sp   = zero
    dpop = (p - p0)/p0
    x_flk  = zero
    y_flk  = zero
    xp_flk = zero
    yp_flk = zero

!++  Transform particle coordinates to get into collimator coordinate
!++  system
!
!++  First check whether particle was lost before
!        if (x.lt.99d-3 .and. z.lt.99d-3) then
!++  First do rotation into collimator frame
    x  = x_in(j) *cos_mb(c_rotation)+sin_mb(c_rotation)*y_in(j)
    z  = y_in(j) *cos_mb(c_rotation)-sin_mb(c_rotation)*x_in(j)
    xp = xp_in(j)*cos_mb(c_rotation)+sin_mb(c_rotation)*yp_in(j)
    zp = yp_in(j)*cos_mb(c_rotation)-sin_mb(c_rotation)*xp_in(j)
!
!++  For one-sided collimators consider only positive X. For negative
!++  X jump to the next particle


! RB: adding exception from goto if it's
    if ((onesided .and. x.lt.zero).and. ((icoll.ne.ipencil) .or. (iturn.ne.1))) goto 777

!++  Now mirror at the horizontal axis for negative X offset
    if(x.lt.zero) then                                             !hr09
      mirror = -one
      tiltangle = -one*c_tilt(2)
    end if

    if(x.ge.zero) then                                             !hr09
      mirror = one
      tiltangle = c_tilt(1)
    end if

    x  = mirror * x
    xp = mirror * xp
!
!          if (j.eq.1) then
!             write(*,*) 'INFOtilt',
!     &            icoll, j_slices, c_tilt(1), c_tilt(2),
!     &            mirror, tiltangle, c_offset, c_aperture/2
!          endif
!
!++  Shift with opening and offset
!
    x  = (x - c_aperture/two) - mirror*c_offset                    !hr09
!
!++  Include collimator tilt
!
    if(tiltangle.gt.zero) then                                    !hr09
      xp = xp - tiltangle
    end if

    if(tiltangle.lt.zero) then
      x  = x + sin_mb(tiltangle) * c_length
      xp = xp - tiltangle
    end if

!++  For selected collimator, first turn reset particle distribution
!++  to simple pencil beam
!
! -- TW why did I set this to 0, seems to be needed for getting
!       right amplitude => no "tilt" of jaw for the first turn !!!!
!          c_tilt(1) = 0d0
!          c_tilt(2) = 0d0

    nprim = 3

    if(( (icoll.eq.ipencil .and. iturn.eq.1) .or. (iturn.eq.1.and. ipencil.eq.999 .and. icoll.le.nprim .and. &
 &    (j.ge.(icoll-1)*nev/nprim) .and. (j.le.(icoll)*nev/nprim))).and.(pencil_distr.ne.3)) then
! RB addition : don't go in this if-statement if pencil_distr=3. This distribution is generated in main loop instead

! -- TW why did I set this to 0, seems to be needed for getting
!       right amplitude => no "tilt" of jaw for the first turn !!!!
      c_tilt(1) = zero
      c_tilt(2) = zero

!AUGUST2006: Standard pencil beam as implemented by GRD ------- TW
      if(pencil_rmsx.eq.zero .and. pencil_rmsy.eq.zero) then     !hr09
        x  = pencil_dx(icoll)
        xp = zero
        z  = zero
        zp = zero
      end if

!AUGUST2006: Rectangular (pencil-beam) sheet-beam with  ------ TW
!            pencil_offset is the rectangulars center
!            pencil_rmsx defines spread of impact parameter
!            pencil_rmsy defines spread parallel to jaw surface

      if(pencil_distr.eq.0 .and.(pencil_rmsx.ne.0..or.pencil_rmsy.ne.0.)) then
! how to assure that all generated particles are on the jaw ?!
        x  = pencil_dx(icoll)+pencil_rmsx*(real(rndm4(),fPrec)-half)
        xp = zero
        z  = pencil_rmsy*(real(rndm4(),fPrec)-half)
        zp = zero
      end if

!AUGUST2006: Gaussian (pencil-beam) sheet-beam with ------- TW
!            pencil_offset is the mean  gaussian distribution
!            pencil_rmsx defines spread of impact parameter
!            pencil_rmsy defines spread parallel to jaw surface

      if(pencil_distr.eq.1 .and.(pencil_rmsx.ne.zero.or.pencil_rmsy.ne.zero )) then
        x  = pencil_dx(icoll) + pencil_rmsx*ran_gauss(two)
! all generated particles are on the jaw now
        x  = sqrt(x**2)
        xp = zero
        z  = pencil_rmsy*ran_gauss(two)
        zp = zero
      end if

!AUGUST2006: Gaussian (pencil-beam) sheet-beam with ------- TW
!            pencil_offset is the mean  gaussian distribution
!            pencil_rmsx defines spread of impact parameter
!                        here pencil_rmsx is not gaussian!!!
!            pencil_rmsy defines spread parallel to jaw surface

      if(pencil_distr.eq.2 .and.(pencil_rmsx.ne.zero.or.pencil_rmsy.ne.zero )) then   !hr09
        x  = pencil_dx(icoll) + pencil_rmsx*(real(rndm4(),fPrec)-half)                  !hr09
! all generated particles are on the jaw now
        x  = sqrt(x**2)
        xp = zero
        z  = pencil_rmsy*ran_gauss(two)
        zp = zero
      end if

!JULY2007: Selection of pos./neg. jaw  implemented by GRD ---- TW

! ensure that for onesided only particles on pos. jaw are created
      if(onesided) then
        mirror = one
      else
!     if(rndm4().lt.0.5) mirror = -1d0
!     if(rndm4().ge.0.5) mirror = 1d0  => using two different random
        if(rndm4().lt.half) then
          mirror = -one
        else
          mirror = one
        end if
      end if

! -- TW SEP07 if c_tilt is set to zero before entering pencil beam
!             section the assigning of the tilt will result in
!             assigning zeros
      if(mirror.lt.zero) then                                     !hr09
!!     tiltangle = -one*c_tilt(2)
        tiltangle = c_tilt(2)
      else
        tiltangle = c_tilt(1)
      end if
!!!!--- commented this out since particle is tilted after leaving
!!!!--- collimator -> remove this  code fragment in final verion
!!             x  = mirror * x
!!             xp = mirror * xp


!++  Include collimator tilt
! this is propably not correct
!
!             xp =  (xp_pencil0*cos_mb(c_rotation)+                         &
!     &            sin_mb(c_rotation)*yp_pencil0)
!             if (tiltangle.gt.0.) then
!                xp = xp - tiltangle
!!             endif
!!             elseif (tiltangle.lt.0.) then
!             else
!               x  = x + sin_mb(tiltangle) * c_length
!               xp = xp - tiltangle
!             endif
!
      write(pencilbeam_distr_unit,'(f10.8,(2x,f10.8),(2x,f10.8),(2x,f10.8),(2x,f10.8))') x, xp, z, zp, tiltangle

    end if !if(( (icoll.eq.ipencil .and. iturn.eq.1) .or. (itu

!          if(rndm4().lt.0.5) mirror = -abs(mirror)
!          if(rndm4().ge.0.5) mirror = abs(mirror)
!        endif
!
!     SR, 18-08-2005: after finishing the coordinate transformation,
!     or the coordinate manipulations in case of pencil beams,
!     write down the initial coordinates of the impacting particles
    xinn  = x
    xpinn = xp
    yinn  = z
    ypinn = zp
!
!++  Possibility to slice here (RA,SR: 29-08-2005)
!
!++  particle passing above the jaw are discarded => take new event
!++  entering by the face, shorten the length (zlm) and keep track of
!++  entrance longitudinal coordinate (keeps) for histograms
!
!++  The definition is that the collimator jaw is at x>=0.
!
!++  1) Check whether particle hits the collimator
    hit   = .false.
    s     = zero                                                !hr09
    keeps = zero                                                !hr09
    zlm   = -one * length

    if(x.ge.zero) then                                            !hr09

!++  Particle hits collimator and we assume interaction length ZLM equal
!++  to collimator length (what if it would leave collimator after
!++  small length due to angle???)
      zlm = length
      impact(j) = x
      indiv(j) = xp
    else if(xp.le.zero) then                                      !hr09

!++  Particle does not hit collimator. Interaction length ZLM is zero.
      zlm = zero
    else
!
!++  Calculate s-coordinate of interaction point
      s = (-one*x) / xp
      if(s.le.0) then
        write(lout,*) 'S.LE.0 -> This should not happen'
        call prror(-1)
      end if

      if(s .lt. length) then
        zlm = length - s
        impact(j) = zero
        indiv(j) = xp
      else
        zlm = zero
      end if
    end if !if(x.ge.0.d0) then

!++  First do the drift part
! DRIFT PART
    drift_length = length - zlm
    if(drift_length.gt.zero) then                                 !hr09
      if(iexact.eq.0) then
        x  = x + xp* drift_length
        z  = z + zp * drift_length
        sp = sp + drift_length
      else
        zpj = sqrt(one-xp**2-zp**2)
        x = x + drift_length*(xp/zpj)
        z = z + drift_length*(zp/zpj)
        sp = sp + drift_length
      end if
    end if
!
!++  Now do the scattering part
!
    if (zlm.gt.zero) then
!JUNE2005
      s_impact = sp
!JUNE2005
      nhit = nhit + 1
!            WRITE(*,*) J,X,XP,Z,ZP,SP,DPOP
!     RB: add new input arguments to jaw icoll,iturn,ipart for writeout
      call jaw(s, nabs, icoll,iturn,name(j),dowrite_impact)

      nabs_type(j) = nabs
!JUNE2005
!JUNE2005 SR+GRD: CREATE A FILE TO CHECK THE VALUES OF IMPACT PARAMETERS
!JUNE2005
!     SR, 29-08-2005: Add to the longitudinal coordinates the position
!     of the slice beginning

      if(dowrite_impact) then
        if(flagsec(j).eq.0) then
#ifdef HDF5
          if(h5_useForCOLL) then
            call h5_prepareWrite(coll_hdf5_fstImpacts, 1)
            call h5_writeData(coll_hdf5_fstImpacts, 1,  1, name(j))
            call h5_writeData(coll_hdf5_fstImpacts, 2,  1, iturn)
            call h5_writeData(coll_hdf5_fstImpacts, 3,  1, icoll)
            call h5_writeData(coll_hdf5_fstImpacts, 4,  1, nabs)
            call h5_writeData(coll_hdf5_fstImpacts, 5,  1, s_impact + (real(j_slices,fPrec)-one) * c_length)
            call h5_writeData(coll_hdf5_fstImpacts, 6,  1, s+sp + (real(j_slices,fPrec)-one) * c_length)
            call h5_writeData(coll_hdf5_fstImpacts, 7,  1, xinn)
            call h5_writeData(coll_hdf5_fstImpacts, 8,  1, xpinn)
            call h5_writeData(coll_hdf5_fstImpacts, 9,  1, yinn)
            call h5_writeData(coll_hdf5_fstImpacts, 10, 1, ypinn)
            call h5_writeData(coll_hdf5_fstImpacts, 11, 1, x)
            call h5_writeData(coll_hdf5_fstImpacts, 12, 1, xp)
            call h5_writeData(coll_hdf5_fstImpacts, 13, 1, z)
            call h5_writeData(coll_hdf5_fstImpacts, 14, 1, zp)
            call h5_finaliseWrite(coll_hdf5_fstImpacts)
          else
#endif
            write(FirstImpacts_unit,'(i5,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))') &
                name(j),iturn,icoll,nabs,                               &
                s_impact + (real(j_slices,fPrec)-one) * c_length,       &
                s+sp + (real(j_slices,fPrec)-one) * c_length,           &
                xinn,xpinn,yinn,ypinn,                                  &
                x,xp,z,zp
#ifdef HDF5
          end if
#endif
        end if
      end if
!!     SR, 18-08-2005: add also the initial coordinates of the
!!                     impacting particles!
!            if(flagsec(j).eq.0) then
!              write(333,'(i5,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))')&
!     +              name(j),iturn,icoll,nabs,s_impact,s+sp,
!     +              xinn,xpinn,yinn,ypinn,
!     +              x,xp,z,zp
!            endif
!     !Old format...
!            if(flagsec(j).eq.0) then
!              write(333,'(i5,1x,i4,1x,i2,1x,i1,2(1x,f5.3),2(1x,e16.7))')
!     &name(j),iturn,icoll,nabs,s_impact,s+sp,impact(j),x
!            endif
!JUNE2005

      lhit_pos(j)  = ie
      lhit_turn(j) = iturn

!-- September2006  TW added from Ralphs code
!--------------------------------------------------------------
!++ Change tilt for pencil beam impact
!
!            if ( (icoll.eq.ipencil                                      &
!     &           .and. iturn.eq.1)   .or.                               &
!     &           (iturn.eq.1 .and. ipencil.eq.999 .and.                 &
!     &                             icoll.le.nprim .and.                 &
!     &            (j.ge.(icoll-1)*nev/nprim) .and.                      &
!     &            (j.le.(icoll)*nev/nprim)                              &
!     &           )  ) then
!
!               if (.not. changed_tilt1(icoll) .and. mirror.gt.0.) then
! ----- Maybe a warning would be nice that c_tilt is overwritten !!!!!
! changed xp_pencil0(icoll) to xp_pencil0 due to definition mismatch
! this has to be solved if necassary and understood
!                 c_tilt(1) = xp_pencil0(icoll)*cos_mb(c_rotation)+         &
!     &                       sin_mb(c_rotation)*yp_pencil0(icoll)
!                 c_tilt(1) = xp_pencil0*cos_mb(c_rotation)+                &
!     &                       sin_mb(c_rotation)*yp_pencil0
!                 write(*,*) "INFO> Changed tilt1  ICOLL  to  ANGLE  ",  &
!     &                   icoll, c_tilt(1), j
!                 changed_tilt1(icoll) = .true.
!               elseif (.not. changed_tilt2(icoll)                       &
!     &                                   .and. mirror.lt.0.) then
! changed xp_pencil0(icoll) to xp_pencil0 due to definition mismatch
! this has to be solved if necassary and understood
!                 c_tilt(2) = -1.*(xp_pencil0(icoll)*cos_mb(c_rotation)+    &
!     &                       sin_mb(c_rotation)*yp_pencil0(icoll))
!                 c_tilt(2) = -1.*(xp_pencil0*cos_mb(c_rotation)+           &
!     &                       sin_mb(c_rotation)*yp_pencil0)
!                 write(*,*) "INFO> Changed tilt2  ICOLL  to  ANGLE  ",  &
!     &                   icoll, c_tilt(2), j
!                 changed_tilt2(icoll) = .true.
!               endif
!            endif
!
!----------------------------------------------------------------
!-- September 2006
!
!++  If particle is absorbed then set x and y to 99.99 mm
!     SR: before assigning new (x,y) for nabs=1, write the
!     inelastic impact file .

!     RB: writeout should be done for both inelastic and single diffractive. doing all transformations
!       in x_flk and making the set to 99.99 mm conditional for nabs=1
!!! /* start RB fix */

! transform back to lab system for writeout.
! keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

      x_flk = xInt
      xp_flk = xpInt

      if(tiltangle.gt.zero) then
        x_flk  = x_flk  + tiltangle*(sInt+sp)
        xp_flk = xp_flk + tiltangle
      else if(tiltangle.lt.zero) then
        xp_flk = xp_flk + tiltangle
        x_flk  = x_flk - sin_mb(tiltangle) * ( length -(sInt+sp) )
      end if

      x_flk  = (x_flk + c_aperture/two) + mirror*c_offset
      x_flk  = mirror * x_flk
      xp_flk = mirror * xp_flk
      y_flk  = yInt   * cos_mb(-one*c_rotation) - x_flk  * sin_mb(-one*c_rotation)
      yp_flk = ypInt  * cos_mb(-one*c_rotation) - xp_flk * sin_mb(-one*c_rotation)
      x_flk  = x_flk  * cos_mb(-one*c_rotation) + yInt   * sin_mb(-one*c_rotation)
      xp_flk = xp_flk * cos_mb(-one*c_rotation) + ypInt  * sin_mb(-one*c_rotation)

! write out all impacts to all_impacts.dat
      if(dowrite_impact) then
        write(FLUKA_impacts_all_unit,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))') &
     &              icoll,c_rotation,                                   &
     &              sInt + sp + (real(j_slices,fPrec)-one) * c_length,  &
     &              x_flk*c1e3, xp_flk*c1e3, y_flk*c1e3, yp_flk*c1e3,   &
     &              nabs,name(j),iturn
      end if

! standard FLUKA_impacts writeout of inelastic and single diffractive
      if((nabs.eq.1).OR.(nabs.eq.4)) then

!     SR, 29-08-2005: Include the slice numer!
        if(dowrite_impact) then
          write(FLUKA_impacts_unit,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))') &
     &icoll,c_rotation,                                                 &
     &sInt + sp + (real(j_slices,fPrec)-one) * c_length,                &!hr09
     &x_flk*c1e3, xp_flk*c1e3, y_flk*c1e3, yp_flk*c1e3,                 &
     &nabs,name(j),iturn
        end if
!
!     Finally, the actual coordinate change to 99 mm
        if(nabs.eq.1) then
          fracab = fracab + 1
          x = 99.99e-3_fPrec
          z = 99.99e-3_fPrec
          part_abs_pos_local(j) = ie
          part_abs_turn_local(j) = iturn
          lint(j) = zlm
        end if
      end if !if((nabs.eq.1).OR.(nabs.eq.4)) then
    end if !if (zlm.gt.0.) then
!!! /* end RB fix */

!++  Do the rest drift, if particle left collimator early
!  DRIFT PART
    if(nabs.ne.1 .and. zlm.gt.zero) then
      drift_length = (length-(s+sp))
      if(drift_length.gt.c1m15) then
        if(iexact.eq.0) then
          x  = x + xp * drift_length
          z  = z + zp * drift_length
          sp = sp + drift_length
        else
          zpj = sqrt(one-xp**2-zp**2)
          x = x + drift_length*(xp/zpj)
          z = z + drift_length*(zp/zpj)
          sp = sp + drift_length
        end if
      end if
      lint(j) = zlm - drift_length
    end if

!++  Transform back to particle coordinates with opening and offset
    if(x.lt.99.0d-3) then

!++  Include collimator tilt
      if(tiltangle.gt.zero) then                                 !hr09
        x  = x  + tiltangle*c_length
        xp = xp + tiltangle
      else if(tiltangle.lt.zero) then                             !hr09
        x  = x + tiltangle*c_length
        xp = xp + tiltangle
        x  = x - sin_mb(tiltangle) * c_length
      end if

!++  Transform back to particle coordinates with opening and offset
      z00 = z
      x00 = x + mirror*c_offset
      x = (x + c_aperture/two) + mirror*c_offset                   !hr09

!++  Now mirror at the horizontal axis for negative X offset
      x  = mirror * x
      xp = mirror * xp

!++  Last do rotation into collimator frame
      x_in(j)  = x  *cos_mb(-one*c_rotation) + z  *sin_mb(-one*c_rotation)
      y_in(j)  = z  *cos_mb(-one*c_rotation) - x  *sin_mb(-one*c_rotation)
      xp_in(j) = xp *cos_mb(-one*c_rotation) + zp *sin_mb(-one*c_rotation)
      yp_in(j) = zp *cos_mb(-one*c_rotation) - xp *sin_mb(-one*c_rotation)

      if(( (icoll.eq.ipencil.and. iturn.eq.1).or. &
  &        (iturn.eq.1 .and.ipencil.eq.999 .and.icoll.le.nprim .and.(j.ge.(icoll-1)*nev/nprim) .and.(j.le.(icoll)*nev/nprim)))&
  &             .and.(pencil_distr.ne.3)) then    ! RB: adding condition that this shouldn't be done if pencil_distr=3

        x00  = mirror * x00
        x_in(j)  = x00  *cos_mb(-one*c_rotation) + z00  *sin_mb(-one*c_rotation)
        y_in(j)  = z00  *cos_mb(-one*c_rotation) - x00  *sin_mb(-one*c_rotation)

        xp_in(j) = xp_in(j) + mirror*xp_pencil0
        yp_in(j) = yp_in(j) + mirror*yp_pencil0
        x_in(j)  = x_in(j)  + mirror*x_pencil(icoll)
        y_in(j)  = y_in(j)  + mirror*y_pencil(icoll)
      end if

      p_in(j) = (one + dpop) * p0
!     SR, 30-08-2005: add the initial position of the slice
      s_in(j) = sp + (real(j_slices,fPrec)-one) * c_length               !hr09

    else
      x_in(j) = x
      y_in(j) = z
    end if !if(x.lt.99.0d-3) then

! output for comparing the particle in accelerator frame
!
!c$$$          if(dowrite_impact) then
!c$$$             write(9996,'(i5,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))')  &
!c$$$     &            name(j),iturn,icoll,nabs,                             &
!c$$$     &            s_in(j),                                              &
!c$$$     &            s+sp + (dble(j_slices)-1d0) * c_length,               &!hr09
!c$$$     &            x_in(j),xp_in(j),y_in(j),yp_in(j),                    &
!c$$$     &            x,xp,z,zp
!c$$$          endif
!
!++  End of check for particles not being lost before
!
!        endif
!
!        IF (X.GT.99.00) WRITE(*,*) 'After : ', X, X_IN(J)
!
!++  End of loop over all particles
!
 777  continue
  end do
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!      WRITE(*,*) 'Number of particles:            ', Nev
!      WRITE(*,*) 'Number of particle hits:        ', Nhit
!      WRITE(*,*) 'Number of absorped particles:   ', fracab
!      WRITE(*,*) 'Number of escaped particles:    ', Nhit-fracab
!      WRITE(*,*) 'Fraction of absorped particles: ', 100.*fracab/Nhit
!
end subroutine collimate2

!>
!! collimaterhic()
!! Collimation for RHIC
!<
subroutine collimaterhic(c_material, c_length, c_rotation,        &
     &     c_aperture, n_aperture,                                      &
     &     c_offset, c_tilt,                                            &
     &     x_in, xp_in, y_in,                                           &
     &     yp_in, p_in, s_in, np, enom,                                 &
     &     lhit_pos,lhit_turn,                                          &
     &     part_abs_pos_local, part_abs_turn_local,                     &
     &     impact, indiv, lint, onesided,                               &
     &     name)
!
!++  Based on routines by JBJ. Changed by RA 2001.
!
!++  - Deleted all HBOOK stuff.
!++  - Deleted optics routine and all parser routines.
!++  - Replaced RANMAR call by RANLUX call
!++  - Included RANLUX code from CERNLIB into source
!++  - Changed dimensions from CGen(100,nmat) to CGen(200,nmat)
!++  - Replaced FUNPRE with FUNLXP
!++  - Replaced FUNRAN with FUNLUX
!++  - Included all CERNLIB code into source: RANLUX, FUNLXP, FUNLUX,
!++                                         FUNPCT, FUNLZ, RADAPT,
!++                                           RGS56P
!++    with additional entries:             RLUXIN, RLUXUT, RLUXAT,
!++                                           RLUXGO
!++
!++  - Changed program so that Nev is total number of particles
!++    (scattered and not-scattered)
!++  - Added debug comments
!++  - Put real dp/dx
!

  use crcoall
  use parpro
  implicit none

!
! BLOCK DBCOLLIM
! This block is common to collimaterhic and collimate2
! It is NOT compatible with block DBCOMMON, as some variable names overlap...


  logical onesided
! integer nprim,filel,mat,nev,j,nabs,nhit,np,icoll,nabs_tmp
  integer np

  integer, allocatable :: lhit_pos(:) !(npart)
  integer, allocatable :: lhit_turn(:) !(npart)
  integer, allocatable :: part_abs_pos_local(:) !(npart)
  integer, allocatable :: part_abs_turn_local(:) !(npart)
  integer, allocatable :: name(:) !(npart)
!MAY2005

  real(kind=fPrec), allocatable :: x_in(:) !(npart)
  real(kind=fPrec), allocatable :: xp_in(:) !(npart)
  real(kind=fPrec), allocatable :: y_in(:) !(npart)
  real(kind=fPrec), allocatable :: yp_in(:) !(npart)
  real(kind=fPrec), allocatable :: p_in(:) !(npart)
  real(kind=fPrec), allocatable :: s_in(:) !(npart)
  real(kind=fPrec), allocatable :: indiv(:) !(npart)
  real(kind=fPrec), allocatable :: lint(:) !(npart)
  real(kind=fPrec), allocatable :: impact(:) !(npart)

  real(kind=fPrec) c_length    !length in m
  real(kind=fPrec) c_rotation  !rotation angle vs vertical in radian
  real(kind=fPrec) c_aperture  !aperture in m
  real(kind=fPrec) c_offset    !offset in m
  real(kind=fPrec) c_tilt(2)   !tilt in radian
  character(len=4) c_material  !material

  real(kind=fPrec) enom
  real(kind=fPrec) n_aperture  !aperture in m for the vertical plane
  save
!=======================================================================
  write(lout,"(a)") "COLL> ERROR collimateRHIC is no longer supported!"
  call prror(-1)
end subroutine collimaterhic
!
!-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----GRD-----
!! END collimaterhic()

!>
!! ichoix(ma)
!! Select a scattering type (elastic, sd, inelastic, ...)
!<
function ichoix(ma)
  implicit none
  integer ma,i,ichoix
  real(kind=fPrec) aran
  aran=real(rndm4(),fPrec)
  i=1
10  if( aran.gt.cprob(i,ma) ) then
      i=i+1
      goto 10
    end if

    ichoix=i
    return
end function ichoix

!>
!! gettran(inter,xmat,p)
!! This function determines: GETTRAN - rms transverse momentum transfer
!! Note: For single-diffractive scattering the vector p of momentum
!! is modified (energy loss is applied)
!<
real(kind=fPrec) function gettran(inter,xmat,p)

  implicit none

  integer, intent(in) :: inter,xmat
  real(kind=fPrec) :: p

  integer :: length
  real(kind=fPrec) :: t,xm2,bsd
  real(kind=fPrec) :: truth,xran(1)

  ! Neither if-statements below have an else, so defaultingfuction return to zero.
  gettran = zero ! -Wmaybe-uninitialized

! inter=2: Nuclear Elastic, 3: pp Elastic, 4: Single Diffractive, 5:Coulomb
#ifndef MERLINSCATTER
  if( inter.eq.2 ) then
    gettran = (-one*log_mb(real(rndm4(),fPrec)))/bn(xmat)                  !hr09

  else if( inter .eq. 3 ) then
    gettran = (-one*log_mb(real(rndm4(),fPrec)))/bpp                       !hr09

  else if( inter .eq. 4 ) then
    xm2 = exp_mb( real(rndm4(),fPrec) * xln15s )
    p = p  * (one - xm2/ecmsq)
    if( xm2 .lt. two ) then
      bsd = two * bpp
    else if (( xm2 .ge. two ).and. ( xm2 .le. five )) then
      bsd = ((106.0_fPrec-17.0_fPrec*xm2) *  bpp )/ 36.0_fPrec             !hr09
!    else if ( xm2 .gt. five ) then
    else !makes the compiler more happy
      bsd = (seven * bpp) / 12.0_fPrec                                     !hr09
    end if
      gettran = (-one*log_mb(real(rndm4(),fPrec)))/bsd                     !hr09

  else if( inter.eq.5 ) then
    length=1
    call funlux( cgen(1,mat), xran, length)
    truth=xran(1)
    t=real(truth,fPrec)                                                    !hr09
    gettran = t
  end if
#else

  if( inter.eq.2 ) then
    gettran = (-one*log_mb(real(rndm4(),fPrec)))/bn(xmat)                  !hr09

  else if( inter .eq. 3 ) then
    call merlinscatter_get_elastic_t(gettran)

  else if( inter .eq. 4 ) then
    call merlinscatter_get_sd_xi(xm2)
    call merlinscatter_get_sd_t(gettran)
    p = p  * (one - (xm2/ecmsq))

  else if ( inter.eq.5 ) then
    length=1
    call funlux( cgen(1,mat) , xran, length)
    truth=xran(1)
    t=real(truth,fPrec)                                                 !hr09
    gettran = t
  end if

#endif
  return
end function gettran

!>
!! tetat(t,p,tx,tz)
!! ???
!!
!<
subroutine tetat(t,p,tx,tz)

  implicit none

  real(kind=fPrec) t,p,tx,tz,va,vb,va2,vb2,r2,teta
  teta = sqrt(t)/p

! Generate sine and cosine of an angle uniform in [0,2pi](see RPP)
10 va  =(two*real(rndm4(),fPrec))-one                                      !hr09
  vb = real(rndm4(),fPrec)
  va2 = va**2
  vb2 = vb**2
  r2 = va2 + vb2
  if ( r2.gt.one) go to 10
  tx = teta * ((two*va)*vb) / r2                                    !hr09
  tz = teta * (va2 - vb2) / r2
  return
end subroutine tetat

!>
!! ruth(t)
!! Calculate the rutherford scattering cross section
!<
function ruth(t)

  implicit none

  real(kind=fPrec) ruth,t
  real(kind=fPrec) cnorm,cnform
  parameter(cnorm=2.607e-5_fPrec,cnform=0.8561e3_fPrec) ! DM: changed 2.607d-4 to 2.607d-5 to fix Rutherford bug

  ruth=(cnorm*exp_mb(((-one*real(t,fPrec))*cnform)*emr(mcurr)**2))*(zatom(mcurr)/real(t,fPrec))**2
end function ruth


!>
!! scatin(plab)
!! Configure the K2 scattering routine cross sections
!!
!<
subroutine scatin(plab)
  use physical_constants

  implicit none

  integer ma,i
  real(kind=fPrec) plab
  real(kind=fPrec) tlow,thigh

  ecmsq = (two * pmap) * plab                                   !hr09
#ifndef MERLINSCATTER
  xln15s=log_mb(0.15_fPrec*ecmsq)                                           !hr09

!Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
  pptot=0.041084_fPrec-0.0023302_fPrec*log_mb(ecmsq)+0.00031514_fPrec*log_mb(ecmsq)**2

!Claudia used the fit from TOTEM for ppel (in barn)
  ppel=(11.7_fPrec-1.59_fPrec*log_mb(ecmsq)+0.134_fPrec*log_mb(ecmsq)**2)/c1e3

!Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
  ppsd=(4.3_fPrec+0.3_fPrec*log_mb(ecmsq))/c1e3
#endif

#ifdef MERLINSCATTER
!No crlibm...
  call merlinscatter_setup(plab,rnd_seed)
  call merlinscatter_setdata(pptot,ppel,ppsd)
#endif

!Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
  bpp=7.156_fPrec+1.439_fPrec*log_mb(sqrt(ecmsq))

! unmeasured tungsten data,computed with lead data and power laws
  bnref(4) = bnref(5)*(anuc(4) / anuc(5))**(two/three)
  emr(4) = emr(5) * (anuc(4)/anuc(5))**(one/three)

! Compute cross-sections (CS) and probabilities + Interaction length
! Last two material treated below statement number 100

  tlow=tlcut                                                   !hr09
  do ma=1,nrmat
    mcurr=ma
! prepare for Rutherford differential distribution
    thigh=hcut(ma)                                             !hr09
    call funlxp ( ruth , cgen(1,ma) ,tlow, thigh )

! freep: number of nucleons involved in single scattering
    freep(ma) = freeco * anuc(ma)**(one/three)

! compute pp and pn el+single diff contributions to cross-section
! (both added : quasi-elastic or qel later)
    csect(3,ma) = freep(ma) * ppel
    csect(4,ma) = freep(ma) * ppsd

! correct TOT-CSec for energy dependence of qel
! TOT CS is here without a Coulomb contribution
    csect(0,ma) = csref(0,ma) + freep(ma) * (pptot - pptref)
    bn(ma) = (bnref(ma) * csect(0,ma)) / csref(0,ma)                    !hr09
! also correct inel-CS
    csect(1,ma) = (csref(1,ma) * csect(0,ma)) / csref(0,ma)                !hr09
!
! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect(2,ma) = ((csect(0,ma) - csect(1,ma)) - csect(3,ma)) - csect(4,ma)         !hr09
    csect(5,ma) = csref(5,ma)
! Now add Coulomb
    csect(0,ma) = csect(0,ma) + csect(5,ma)
! Interaction length in meter
  xintl(ma) = (c1m2*anuc(ma))/(((fnavo * rho(ma))*csect(0,ma))*1d-24) !hr09

! Filling CProb with cumulated normalised Cross-sections
    do i=1,4
      cprob(i,ma)=cprob(i-1,ma)+csect(i,ma)/csect(0,ma)
    end do
  end do

! Last two materials for 'vaccum' (nmat-1) and 'full black' (nmat)
  cprob(1,nmat-1) = one
  cprob(1,nmat)   = one
  xintl(nmat-1)   = c1e12
  xintl(nmat)     = zero
  return
end subroutine scatin

!>
!! jaw(s,nabs,icoll,iturn,ipart,dowrite_impact)
!! ???
!!     RB: adding as input arguments to jaw variables icoll,iturn,ipart
!!         these are only used for the writeout of particle histories
!!
!!++  Input:   ZLM is interaction length
!!++           MAT is choice of material
!!
!!++  Output:  nabs = 1   Particle is absorped
!!++           nabs = 4   Single-diffractive scattering
!!++           dpop       Adjusted for momentum loss (dE/dx)
!!++           s          Exit longitudinal position
!!
!!++  Physics:  If monte carlo interaction length greater than input
!!++            interaction length, then use input interaction length
!!++            Is that justified???
!!
!!     nabs=1....absorption
!!
!<
subroutine jaw(s,nabs,icoll,iturn,ipart,dowrite_impact)

  implicit none

  integer nabs,inter,iturn,icoll,ipart,nabs_tmp ! RB: added variables icoll,iturn,ipart for writeout
  logical dowrite_impact
  real(kind=fPrec) m_dpodx     !CT, RB, DM
  real(kind=fPrec) p,rlen,s,t,dxp,dzp,p1,zpBef,xpBef,pBef
!...cne=1/(sqrt(b))
!...dpodx=dE/(dx*c)

!++  Note that the input parameter is dpop. Here the momentum p is
!++  constructed out of this input.

  p=p0*(one+dpop)
  nabs=0
  nabs_tmp=nabs

  if(mat.eq.nmat) then
!++  Collimator treated as black absorber
    nabs=1
    nabs_tmp=nabs
    s=zero

    if(dowrite_impact) then
      ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
      if(h5_useForCOLL) then
        call coll_hdf5_writeCollScatter(icoll, iturn, ipart, nabs_tmp, -one, zero, zero)
      else
#endif
      write(coll_scatter_unit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))') icoll, iturn, ipart, nabs_tmp, -one, zero, zero
#ifdef HDF5
      endif
#endif
    end if
    return
  else if(mat.eq.nmat-1) then
!++  Collimator treated as drift
    s=zlm
    x=x+s*xp
    z=z+s*zp

    if(dowrite_impact) then
      ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
      if(h5_useForCOLL) then
        call coll_hdf5_writeCollScatter(icoll, iturn, ipart, nabs_tmp, -one, zero, zero)
      else
#endif
      write(coll_scatter_unit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))') icoll, iturn, ipart, nabs_tmp, -one, zero, zero
#ifdef HDF5
      endif
#endif
    end if

    return
  end if

!++  Initialize the interaction length to input interaction length
  rlen=zlm

!++  Do a step for a point-like interaction. This is a loop with
!++  label 10!!!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!++  Get monte-carlo interaction length.

10  zlm1=(-one*xintl(mat))*log_mb(real(rndm4(),fPrec))                          !hr09
  nabs_tmp=0 !! type of interaction reset before following scattering process
  xpBef=xp ! save angles and momentum before scattering
  zpBef=zp
  pBef=p

!++  If the monte-carlo interaction length is longer than the
!++  remaining collimator length, then put it to the remaining
!++  length, do multiple coulomb scattering and return.
!++  LAST STEP IN ITERATION LOOP
  if(zlm1.gt.rlen) then
    zlm1=rlen
    call mcs(s)
    s=(zlm-rlen)+s                                                    !hr09
#ifndef MERLINSCATTER
    call calc_ion_loss(mat,p,rlen,m_dpodx)  ! DM routine to include tail
    p=p-m_dpodx*s
#endif
#ifdef MERLINSCATTER
!void calc_ion_loss_merlin_(double* p, double* ElectronDensity, double* PlasmaEnergy, double* MeanIonisationEnergy, double* result)
    call merlinscatter_calc_ion_loss(p,edens(mat), pleng(mat),exenergy(mat),s,m_dpodx)
    p=p-m_dpodx
#endif

    dpop=(p-p0)/p0
    if(dowrite_impact) then
      ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
      if(h5_useForCOLL) then
        call coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef)
      else
#endif
      write(coll_scatter_unit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))') icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
#ifdef HDF5
      endif
#endif
    end if
    return
  end if

!++  Otherwise do multi-coulomb scattering.
!++  REGULAR STEP IN ITERATION LOOP
  call mcs(s)

!++  Check if particle is outside of collimator (X.LT.0) after
!++  MCS. If yes, calculate output longitudinal position (s),
!++  reduce momentum (output as dpop) and return.
!++  PARTICLE LEFT COLLIMATOR BEFORE ITS END.
  if(x.le.zero) then
    s=(zlm-rlen)+s                                                    !hr09

#ifndef MERLINSCATTER
    call calc_ion_loss(mat,p,rlen,m_dpodx)
    p=p-m_dpodx*s
#endif
#ifdef MERLINSCATTER
    call merlinscatter_calc_ion_loss(p,edens(mat),pleng(mat),exenergy(mat),s,m_dpodx)
    p=p-m_dpodx
#endif
    dpop=(p-p0)/p0

    if(dowrite_impact) then
      ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
      if(h5_useForCOLL) then
        call coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef)
      else
#endif
      write(coll_scatter_unit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))') icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
#ifdef HDF5
      endif
#endif
    end if

    return
  end if

!++  Check whether particle is absorbed. If yes, calculate output
!++  longitudinal position (s), reduce momentum (output as dpop)
!++  and return.
!++  PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.

  inter=ichoix(mat)
  nabs=inter
  nabs_tmp=nabs

! RB, DM: save coordinates before interaction for writeout to FLUKA_impacts.dat
  xInt=x
  xpInt=xp
  yInt=z
  ypInt=zp
  sInt=(zlm-rlen)+zlm1                                                 !hr09

  if(inter.eq.1) then
    s=(zlm-rlen)+zlm1                                                 !hr09

#ifndef MERLINSCATTER
    call calc_ion_loss(mat,p,rlen,m_dpodx)
    p=p-m_dpodx*s
#endif
#ifdef MERLINSCATTER
    call merlinscatter_calc_ion_loss(p,edens(mat),pleng(mat),exenergy(mat),s,m_dpodx)
    p=p-m_dpodx
#endif

    dpop=(p-p0)/p0

    ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
    if(h5_useForCOLL) then
      call coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs_tmp,-one,zero,zero)
    else
#endif
    write(coll_scatter_unit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))') icoll,iturn,ipart,nabs_tmp,-one,zero,zero
#ifdef HDF5
    endif
#endif
    return
  end if

!++  Now treat the other types of interaction, as determined by ICHOIX:

!++      Nuclear-Elastic:          inter = 2
!++      pp Elastic:               inter = 3
!++      Single-Diffractive:       inter = 4    (changes momentum p)
!++      Coulomb:                  inter = 5

!++  As the single-diffractive interaction changes the momentum, save
!++  input momentum in p1.
  p1 = p

!++  Gettran returns some monte carlo number, that, as I believe, gives
!++  the rms transverse momentum transfer.
  t = gettran(inter,mat,p)

!++  Tetat calculates from the rms transverse momentum transfer in
!++  monte-carlo fashion the angle changes for x and z planes. The
!++  angle change is proportional to SQRT(t) and 1/p, as expected.
  call tetat(t,p,dxp,dzp)

!++  Apply angle changes
  xp=xp+dxp
  zp=zp+dzp

!++  Treat single-diffractive scattering.
  if(inter.eq.4) then

!++ added update for s
    s=(zlm-rlen)+zlm1                                                !hr09
    xpsd=dxp
    zpsd=dzp
    psd=p1
!
!++  Add this code to get the momentum transfer also in the calling
!++  routine...
    dpop=(p-p0)/p0
  end if

  if(dowrite_impact) then
    ! write coll_scatter.dat for complete scattering histories.
    ! Includes changes in angle from both
#ifdef HDF5
    if(h5_useForCOLL) then
      call coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef)
    else
#endif
    write(coll_scatter_unit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))') icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
#ifdef HDF5
    endif
#endif
  end if

!++  Calculate the remaining interaction length and close the iteration
!++  loop.
  rlen=rlen-zlm1
  goto 10

end subroutine jaw

#ifdef HDF5
subroutine coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs,dp,dx,dy)

  integer,          intent(in) :: icoll,iturn,ipart,nabs
  real(kind=fPrec), intent(in) :: dp,dx,dy

  call h5_prepareWrite(coll_hdf5_collScatter, 1)
  call h5_writeData(coll_hdf5_collScatter, 1, 1, ipart)
  call h5_writeData(coll_hdf5_collScatter, 2, 1, iturn)
  call h5_writeData(coll_hdf5_collScatter, 3, 1, icoll)
  call h5_writeData(coll_hdf5_collScatter, 4, 1, nabs)
  call h5_writeData(coll_hdf5_collScatter, 5, 1, dp)
  call h5_writeData(coll_hdf5_collScatter, 6, 1, dx)
  call h5_writeData(coll_hdf5_collScatter, 7, 1, dy)
  call h5_finaliseWrite(coll_hdf5_collScatter)

end subroutine coll_hdf5_writeCollScatter
#endif

!>
!! mcs(s)
!!++  Input:   zlm1   Monte-carlo interaction length
!!
!!++  Output:  s      Longitudinal position
!!++           p0     Reference momentum
!!++           dpop   Relative momentum offset
!!
!!     collimator: x>0 and y<zlm1
!<
subroutine mcs(s)
      implicit none
!      save h,dh,bn
      real(kind=fPrec) h,dh,theta,rlen0,rlen,ae,be,bn0,s
      real(kind=fPrec) radl_mat,rad_len ! Claudia 2013 added variables


!   bn=sqrt(3)/(number of sigmas for s-determination(=4))
      data h/.001d0/dh/.0001d0/bn0/.4330127019d0/

      radl_mat=radl(mat)
      theta=13.6d-3/(p0*(1.d0+dpop))      !Claudia added log part
      rad_len=radl(mat)                    !Claudia

      x=(x/theta)/radl(mat)                                              !hr09
      xp=xp/theta
      z=(z/theta)/radl(mat)                                              !hr09
      zp=zp/theta
      rlen0=zlm1/radl(mat)
      rlen=rlen0
10    ae=bn0*x
      be=bn0*xp
      call soln3(ae,be,dh,rlen,s)
      if(s.lt.h) s=h
      call scamcs(x,xp,s,radl_mat)
      if(x.le.0.d0) then
       s=(rlen0-rlen)+s                                                  !hr09
       goto 20
      end if
      if(s+dh.ge.rlen) then
       s=rlen0
       goto 20
      end if
      rlen=rlen-s
      goto 10
20    call scamcs(z,zp,s,radl_mat)
      s=s*radl(mat)
      x=(x*theta)*radl(mat)                                              !hr09
      xp=xp*theta
      z=(z*theta)*radl(mat)                                              !hr09
      zp=zp*theta
end subroutine mcs

!>
!! scamcs(xx,xxp,s,radl_mat)
!! ???
!<
subroutine scamcs(xx,xxp,s,radl_mat)

  implicit none

  real(kind=fPrec) v1,v2,r2,a,z1,z2,ss,s,xx,xxp,x0,xp0
  real(kind=fPrec) radl_mat

  x0=xx
  xp0=xxp

5 v1=2d0*real(rndm4(),fPrec)-1d0
  v2=2d0*real(rndm4(),fPrec)-1d0
  r2=v1**2+v2**2                                                     !hr09
  if(r2.ge.1.d0) goto 5

  a=sqrt((-2.d0*log_mb(r2))/r2)                                         !hr09
  z1=v1*a
  z2=v2*a
  ss=sqrt(s)
  xx=x0+s*(xp0+(half*ss)*(one+0.038_fPrec*log_mb(s))*(z2+z1*0.577350269_fPrec)) !Claudia: added logarithmic part in mcs formula
  xxp=xp0+ss*z2*(one+0.038_fPrec*log_mb(s))
end subroutine scamcs

!>
!! soln3(a,b,dh,smax,s)
!! ???
!<
subroutine soln3(a,b,dh,smax,s)

  implicit none

  real(kind=fPrec) b,a,s,smax,c,dh
  if(b.eq.zero) then
    s=a**0.6666666666666667_fPrec
!      s=a**(two/three)
    if(s.gt.smax) s=smax
    return
  end if

  if(a.eq.zero) then
    if(b.gt.zero) then
      s=b**2
    else
      s=zero
    end if
    if(s.gt.smax) s=smax
    return
  end if

  if(b.gt.zero) then
    if(smax**3.le.(a+b*smax)**2) then
      s=smax
      return
    else
      s=smax*half
      call iterat(a,b,dh,s)
    end if
  else
    c=(-one*a)/b
    if(smax.lt.c) then
      if(smax**3.le.(a+b*smax)**2) then
        s=smax
        return
      else
        s=smax*half
        call iterat(a,b,dh,s)
      end if
    else
      s=c*half
      call iterat(a,b,dh,s)
    end if
  end if

end subroutine soln3

subroutine iterat(a,b,dh,s)

  implicit none

  real(kind=fPrec) ds,s,a,b,dh

  ds=s
10 ds=ds*half

  if(s**3.lt.(a+b*s)**2) then
    s=s+ds
  else
    s=s-ds
  end if

  if(ds.lt.dh) then
    return
  else
    goto 10
  end if

end subroutine iterat

!>
!! get_dpodx(p,mat_i)
!! calculate mean ionization energy loss according to Bethe-Bloch
!<
function get_dpodx(p, mat_i)          !Claudia
  use physical_constants

  implicit none

  real(kind=fPrec), intent(in) :: p
  integer, intent(in) :: mat_i

  real(kind=fPrec) PE,K,gamma_p
  real(kind=fPrec) beta_p,gamma_s,beta_s,me2,mp2,T,part_1,part_2,I_s,delta
  parameter(K=0.307075)
  real(kind=fPrec) get_dpodx

  mp2       = pmap**2
  me2       = pmae**2
  beta_p    = one
  gamma_p   = p/pmap
  beta_s    = beta_p**2
  gamma_s   = gamma_p**2
  T         = (2*pmae*beta_s*gamma_s)/(1+(2*gamma_p*pmae/pmap)+me2/mp2)
  PE        = sqrt(rho(mat_i)*zatom(mat_i)/anuc(mat_i))*28.816e-9_fPrec
  I_s       = exenergy(mat_i)**2
  part_1    = K*zatom(mat_i)/(anuc(mat_i)*beta_s)
  delta     = log_mb(PE/exenergy(mat_i))+log_mb(beta_p*gamma_p)-half
  part_2    = half*log_mb((two*pmae*beta_s*gamma_s*T)/I_s)
  get_dpodx = part_1*(part_2-beta_s-delta)*rho(mat_i)*c1m1
  return
end function get_dpodx

!>
!! CalcElectronDensity(AtomicNumber, Density, AtomicMass)
!! Function to calculate the electron density in a material
!! Should give the number per cubic meter
!<
function CalcElectronDensity(AtomicNumber, Density, AtomicMass)
  implicit none

  real(kind=fPrec) AtomicNumber, Density, AtomicMass
  real(kind=fPrec) Avogadro
  real(kind=fPrec) CalcElectronDensity
  real(kind=fPrec) PartA, PartB
  parameter (Avogadro = 6.022140857e23_fPrec)
  PartA = AtomicNumber * Avogadro * Density
  !1e-6 factor converts to n/m^-3
  PartB = AtomicMass * c1m6
  CalcElectronDensity = PartA/PartB
  return
end function CalcElectronDensity

!>
!! CalcPlasmaEnergy(ElectronDensity)
!! Function to calculate the plasma energy in a material
!! CalculatePlasmaEnergy = (PlanckConstantBar * sqrt((ElectronDensity *(ElectronCharge**2)) / &
!!& (ElectronMass * FreeSpacePermittivity)))/ElectronCharge*eV;
!<
function CalcPlasmaEnergy(ElectronDensity)

  implicit none

  real(kind=fPrec) ElectronDensity
  real(kind=fPrec) CalcPlasmaEnergy
  real(kind=fPrec) sqrtAB,PartA,PartB,FSPC2

  !Values from the 2016 PDG
  real(kind=fPrec) PlanckConstantBar,ElectronCharge,ElectronMass
  real(kind=fPrec) ElectronCharge2
  real(kind=fPrec) FreeSpacePermittivity,FreeSpacePermeability
  real(kind=fPrec) SpeedOfLight,SpeedOfLight2

  parameter (PlanckConstantBar = 1.054571800e-34_fPrec)
  parameter (ElectronCharge = 1.6021766208e-19_fPrec)
  parameter (ElectronCharge2 = ElectronCharge*ElectronCharge)
  parameter (ElectronMass = 9.10938356e-31_fPrec)
  parameter (SpeedOfLight = 299792458.0_fPrec)
  parameter (SpeedOfLight2 = SpeedOfLight*SpeedOfLight)

  parameter (FreeSpacePermeability = 16.0e-7_fPrec*atan(one)) ! Henry per meter
  parameter (FSPC2 = FreeSpacePermeability*SpeedOfLight2)
  parameter (FreeSpacePermittivity = one/FSPC2)
  parameter (PartB = ElectronMass * FreeSpacePermittivity)

  PartA = ElectronDensity * ElectronCharge2

  sqrtAB = sqrt(PartA/PartB)
  CalcPlasmaEnergy=PlanckConstantBar*sqrtAB/ElectronCharge*c1m9
  return
end function CalcPlasmaEnergy

!>
!! calc_ion_loss(IS,PC,DZ,EnLo)
!! subroutine for the calculazion of the energy loss by ionization
!! Either mean energy loss from Bethe-Bloch, or higher energy loss, according to finite probability from cross section
!! written by DM for crystals, introduced in main code by RB
!<
subroutine calc_ion_loss(IS, PC, DZ, EnLo)

! IS material ID
! PC momentum in GeV
! DZ length traversed in material (meters)
! EnLo energy loss in GeV/meter

  use physical_constants

  implicit none

  integer IS

  real(kind=fPrec) PC,DZ,EnLo,exEn
  real(kind=fPrec) k !Daniele: parameters for dE/dX calculation (const,electron radius,el. mass, prot.mass)
  real(kind=fPrec) enr,mom,betar,gammar,bgr !Daniele: energy,momentum,beta relativistic, gamma relativistic
  real(kind=fPrec) Tmax,plen !Daniele: maximum energy tranfer in single collision, plasma energy (see pdg)
  real(kind=fPrec) thl,Tt,cs_tail,prob_tail
  real(kind=fPrec) ranc

  data k/0.307075_fPrec/      !constant in front bethe-bloch [MeV g^-1 cm^2]
! The following values are now taken from physical_constants
!  data re/2.818d-15/    !electron radius [m]
!  data me/0.510998910/  !electron mass [MeV/c^2]
!  data mp/938.272013/   !proton mass [MeV/c^2]

  mom    = PC*c1e3                    ! [GeV/c] -> [MeV/c]
  enr    = (mom*mom+pmap*pmap)**half  ! [MeV]
  gammar = enr/pmap
  betar  = mom/enr
  bgr    = betar*gammar

! mean excitation energy - convert to MeV
  exEn=exenergy(IS)*c1e3

! Tmax is max energy loss from kinematics
  Tmax=(two*pmae*bgr**2)/(one+two*gammar*pmae/pmap+(pmae/pmap)**2) ![MeV]

! plasma energy - see PDG 2010 table 27.1
  plen = ((rho(IS)*zatom(IS)/anuc(IS))**half)*28.816e-6_fPrec ![MeV]

! calculate threshold energy
! Above this threshold, the cross section for high energy loss is calculated and then
! a random number is generated to determine if tail energy loss should be applied, or only mean from Bethe-Bloch
! below threshold, only the standard bethe-bloch is used (all particles get average energy loss)

! thl is 2* width of landau distribution (as in fig 27.7 in PDG 2010). See Alfredo's presentation for derivation
  thl = four*k*zatom(IS)*DZ*c1e2*rho(IS)/(anuc(IS)*betar**2) ![MeV]
!     write(3456,*) thl     ! should typically be >0.06MeV for approximations to be valid - check!

! Bethe Bloch mean energy loss
  EnLo = ((k*zatom(IS))/(anuc(IS)*betar**2))*(half*log_mb((two*pmae*bgr*bgr*Tmax)/(exEn*exEn))-betar**two-&
& log_mb(plen/exEn)-log_mb(bgr)+half)

  EnLo = EnLo*rho(IS)*c1m1*DZ  ![GeV]

! threshold Tt is bethe bloch + 2*width of Landau distribution
  Tt = EnLo*c1e3+thl      ![MeV]

! cross section - see Alfredo's presentation for derivation
  cs_tail = ((k*zatom(IS))/(anuc(IS)*betar**2))*((half*((one/Tt)-(one/Tmax)))-(log_mb(Tmax/Tt)*(betar**2) &
 &        /(two*Tmax))+((Tmax-Tt)/(four*(gammar**2)*(pmap**2))))

! probability of being in tail: cross section * density * path length
  prob_tail = cs_tail*rho(IS)*DZ*c1e2;

  ranc = real(rndm4(),fPrec)

! determine based on random number if tail energy loss occurs.
  if(ranc.lt.prob_tail) then
    EnLo = ((k*zatom(IS))/(anuc(IS)*betar**2))*(half*log_mb((two*pmae*bgr*bgr*Tmax)/(exEn*exEn))-betar**two- &
 &       log_mb(plen/exEn)-log_mb(bgr)+half+(TMax**2)/(eight*(gammar**2)*(pmap**2)))

    EnLo = EnLo*rho(IS)*c1m1 ![GeV/m]
  else
    ! if tial energy loss does not occur, just use the standard Bethe Bloch
    EnLo = EnLo/DZ  ![GeV/m]
  endif

  RETURN

end subroutine calc_ion_loss

subroutine makedis(myalphax, myalphay, mybetax, mybetay,    &
     &myemitx0, myemity0, myenom, mynex, mdex, myney, mdey,             &
     &myx, myxp, myy, myyp, myp, mys)

!  Generate distribution

  use crcoall
  use mod_common, only : napx
  implicit none

  integer :: j
  real(kind=fPrec), allocatable :: myx(:) !(npart)
  real(kind=fPrec), allocatable :: myxp(:) !(npart)
  real(kind=fPrec), allocatable :: myy(:) !(npart)
  real(kind=fPrec), allocatable :: myyp(:) !(npart)
  real(kind=fPrec), allocatable :: myp(:) !(npart)
  real(kind=fPrec), allocatable :: mys(:) !(npart)

  real(kind=fPrec) myalphax,mybetax,myemitx0,myemitx,mynex,mdex, &
  &mygammax,myalphay,mybetay,myemity0,myemity,myney,mdey,mygammay,   &
  &xsigmax,ysigmay,myenom

  save
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!
!++  Calculate the gammas
!
  mygammax = (one+myalphax**2)/mybetax
  mygammay = (one+myalphay**2)/mybetay
!++TW 11/07 reset j, helps if subroutine is called twice
! was done during try to reset distribution, still needed
! will this subroutine ever called twice?
  j = 0
!
!++  Number of points and generate distribution
  write(lout,"(a)")
  write(lout,"(a)") 'COLL> Generation of particle distribution Version 1:'
  write(lout,"(a)") 'COLL> This routine generates particles in phase space X/XP and Y/YP ellipses, as defined in the input '
  write(lout,"(a)") 'COLL> parameters. Distribution is flat in the band. X and Y are fully uncorrelated.'
  write(lout,"(a)")

  write(outlun,*)
  write(outlun,*) 'Generation of particle distribution Version 1'
  write(outlun,*)
  write(outlun,*) 'This routine generates particles in phase space'
  write(outlun,*) 'X/XP and Y/YP ellipses, as defined in the input'
  write(outlun,*) 'parameters. Distribution is flat in the band.'
  write(outlun,*) 'X and Y are fully uncorrelated.'
  write(outlun,*)
  write(outlun,*) 'INFO>  Number of particles   = ', napx
  write(outlun,*) 'INFO>  Av number of x sigmas = ', mynex
  write(outlun,*) 'INFO>  +- spread in x sigmas = ', mdex
  write(outlun,*) 'INFO>  Av number of y sigmas = ', myney
  write(outlun,*) 'INFO>  +- spread in y sigmas = ', mdey
  write(outlun,*) 'INFO>  Nominal beam energy   = ', myenom
  write(outlun,*) 'INFO>  Sigma_x0 = ', sqrt(mybetax*myemitx0)
  write(outlun,*) 'INFO>  Sigma_y0 = ', sqrt(mybetay*myemity0)
  write(outlun,*) 'INFO>  Beta x   = ', mybetax
  write(outlun,*) 'INFO>  Beta y   = ', mybetay
  write(outlun,*) 'INFO>  Alpha x  = ', myalphax
  write(outlun,*) 'INFO>  Alpha y  = ', myalphay
  write(outlun,*)

  do while (j.lt.napx)
    j = j + 1
    myemitx = myemitx0*(mynex + ((two*real(rndm4()-half,fPrec))*mdex) )**2
    xsigmax = sqrt(mybetax*myemitx)
    myx(j)  = xsigmax * sin_mb((two*pi)*real(rndm4(),fPrec))
    if(rndm4().gt.half) then
      myxp(j) = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax
    else
      myxp(j) = -one*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax
    end if

    myemity = myemity0*(myney + ((two*real(rndm4()-half,fPrec))*mdey) )**2
    ysigmay = sqrt(mybetay*myemity)
    myy(j)  = ysigmay * sin_mb((two*pi)*real(rndm4(),fPrec))
    if(rndm4().gt.half) then
      myyp(j) = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay
    else
      myyp(j) = -one*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay
    end if

    myp(j) = myenom
    mys(j) = zero

!++  Dangerous stuff, just for the moment
    if (cut_input) then
      !0.1d-3 -> c1m4
      if((.not. (myy(j).lt.-0.008e-3_fPrec .and. myyp(j).lt. c1m4 .and.myyp(j).gt.zero) ) .and. &
&        (.not. (myy(j).gt. 0.008e-3_fPrec .and. myyp(j).gt.-c1m4 .and.myyp(j).lt.zero) ) ) then
        j = j - 1
      end if
    end if
  end do
  return
end subroutine makedis

!========================================================================
! SR, 08-05-2005: Add the finite beam size in the othe dimension
subroutine makedis_st(myalphax, myalphay, mybetax, mybetay, &
     &     myemitx0, myemity0, myenom, mynex, mdex, myney, mdey,  &
     &     myx, myxp, myy, myyp, myp, mys)

!     Uses the old routine 'MAKEDIS' for the halo plane and adds the
!     transverse beam size in the other plane (matched distrubutions
!     are generated starting from thetwiss functions).
!     If 'mynex' and 'myney' are BOTH set to zero, nominal bunches
!     centred in the aperture centre are generated. (SR, 08-05-2005)

  use crcoall
  use mod_common, only : napx
  implicit none

  integer :: j
  real(kind=fPrec), allocatable :: myx(:) !(npart)
  real(kind=fPrec), allocatable :: myxp(:) !(npart)
  real(kind=fPrec), allocatable :: myy(:) !(npart)
  real(kind=fPrec), allocatable :: myyp(:) !(npart)
  real(kind=fPrec), allocatable :: myp(:) !(npart)
  real(kind=fPrec), allocatable :: mys(:) !(npart)

  real(kind=fPrec) myalphax,mybetax,myemitx0,myemitx,mynex,mdex, &
  &mygammax,myalphay,mybetay,myemity0,myemity,myney,mdey,mygammay,   &
  &xsigmax,ysigmay,myenom

  real(kind=fPrec) iix, iiy, phix, phiy
  save

!-----------------------------------------------------------------------
!++  Generate particle distribution
!++  Generate random distribution, assuming optical parameters at IP1
!++  Calculate the gammas
  mygammax = (one+myalphax**2)/mybetax
  mygammay = (one+myalphay**2)/mybetay
  do j=1, napx
    if((mynex.gt.zero).and.(myney.eq.zero)) then
      myemitx = myemitx0*(mynex+((two*real(rndm4()-half,fPrec))*mdex))**2
      xsigmax = sqrt(mybetax*myemitx)
      myx(j)  = xsigmax * sin_mb((two*pi)*real(rndm4(),fPrec))

      if (rndm4().gt.half) then
        myxp(j) =     sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax
      else
       myxp(j) = -one*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax
      end if

      phiy = (two*pi)*real(rndm4(),fPrec)
      iiy = (-one*myemity0) * log_mb( real(rndm4(),fPrec) )
      myy(j) = sqrt((two*iiy)*mybetay) * cos_mb(phiy)
      myyp(j) = (-one*sqrt((two*iiy)/mybetay)) * (sin_mb(phiy) + myalphay * cos_mb(phiy))

    else if ( mynex.eq.zero.and.myney.gt.zero ) then
      myemity = myemity0*(myney+((two*real(rndm4()-half,fPrec))*mdey))**2
      ysigmay = sqrt(mybetay*myemity)
      myy(j)  = ysigmay * sin_mb((two*pi)*real(rndm4(),fPrec))

      if (rndm4().gt.half) then
        myyp(j) =      sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay
      else
        myyp(j) = -one*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay
      end if

      phix = (two*pi)*real(rndm4(),fPrec)
      iix = (-one* myemitx0) * log_mb( real(rndm4(),fPrec) )
      myx(j) = sqrt((two*iix)*mybetax) * cos_mb(phix)
      myxp(j) = (-one*sqrt((two*iix)/mybetax)) * (sin_mb(phix) +myalphax * cos_mb(phix))

    else if( mynex.eq.zero.and.myney.eq.zero ) then
      phix = (two*pi)*real(rndm4(),fPrec)
      iix = (-one*myemitx0) * log_mb( real(rndm4(),fPrec) )
      myx(j) = sqrt((two*iix)*mybetax) * cos_mb(phix)
      myxp(j) = (-one*sqrt((two*iix)/mybetax)) * (sin_mb(phix) +myalphax * cos_mb(phix))
      phiy = (two*pi)*real(rndm4(),fPrec)
      iiy = (-one*myemity0) * log_mb( real(rndm4(),fPrec) )
      myy(j) = sqrt((two*iiy)*mybetay) * cos_mb(phiy)
      myyp(j) = (-one*sqrt((two*iiy)/mybetay)) * (sin_mb(phiy) + myalphay * cos_mb(phiy))
    else
      write(lout,"(a)") "COLL> ERROR Bbeam parameters not correctly set!"
    end if
    myp(j) = myenom
    mys(j) = zero
  end do
  return
end subroutine makedis_st

!========================================================================
!
!     RB: new routine to sample part of matched phase ellipse which is outside
!     the cut of the jaws
!     Assuming cut of the jaw at mynex for hor plane.
!     largest amplitude outside of jaw is mynex + mdex.  Analog for vertical plane.

!     same routine as makedis_st, but rejection sampling to get
!     only particles hitting the collimator on the same turn.

!     Treat as a pencil beam in main routine.

subroutine makedis_coll(myalphax, myalphay, mybetax, mybetay,  myemitx0, myemity0, &
 &                        myenom, mynex, mdex, myney, mdey, myx, myxp, myy, myyp, myp, mys)

  use crcoall
  use mod_common, only : napx
  implicit none

  integer :: j
  real(kind=fPrec), allocatable :: myx(:) !(npart)
  real(kind=fPrec), allocatable :: myxp(:) !(npart)
  real(kind=fPrec), allocatable :: myy(:) !(npart)
  real(kind=fPrec), allocatable :: myyp(:) !(npart)
  real(kind=fPrec), allocatable :: myp(:) !(npart)
  real(kind=fPrec), allocatable :: mys(:) !(npart)

  real(kind=fPrec) myalphax,mybetax,myemitx0,myemitx,mynex,mdex, &
  &mygammax,myalphay,mybetay,myemity0,myemity,myney,mdey,mygammay,   &
  &xsigmax,ysigmay,myenom


  real(kind=fPrec) iix, iiy, phix,phiy,cutoff

  save
!
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!++  Calculate the gammas

  mygammax = (one+myalphax**2)/mybetax
  mygammay = (one+myalphay**2)/mybetay

! calculate cutoff in x or y from the collimator jaws.
  if((mynex.gt.zero).and.(myney.eq.zero)) then
    cutoff=mynex*sqrt(mybetax*myemitx0)
  else
    cutoff=myney*sqrt(mybetay*myemity0)
  end if

      do j=1, napx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if((mynex.gt.zero).and.(myney.eq.zero)) then  ! halo in x
 887        continue
            myemitx = myemitx0*(mynex+(real(rndm4(),fPrec)*mdex))**2
            xsigmax = sqrt(mybetax*myemitx)
            myx(j)  = xsigmax * sin_mb((two*pi)*real(rndm4(),fPrec))
            if(abs(myx(j)).lt.cutoff) goto 887

            if(rndm4().gt.half) then
              myxp(j) = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax
            else
              myxp(j) = -one*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax
            end if

            phiy = (two*pi)*real(rndm4(),fPrec)
            iiy = (-one*myemity0) * log_mb( real(rndm4(),fPrec) )
            myy(j) = sqrt((two*iiy)*mybetay) * cos_mb(phiy)
            myyp(j) = (-one*sqrt((two*iiy)/mybetay)) * (sin_mb(phiy) + myalphay * cos_mb(phiy))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         else if( mynex.eq.zero.and.myney.gt.zero ) then  ! halo in y
 886        continue
            myemity = myemity0*(myney+(real(rndm4(),fPrec)*mdey))**2
            ysigmay = sqrt(mybetay*myemity)
            myy(j)   = ysigmay * sin_mb((two*pi)*real(rndm4(),fPrec))
            if(abs(myy(j)).lt.cutoff) goto 886

            if(rndm4().gt.half) then
              myyp(j) = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay
            else
              myyp(j) = -one*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay
            end if

            phix = (two*pi)*real(rndm4(),fPrec)
            iix = (-one* myemitx0) * log_mb( real(rndm4(),fPrec) )
            myx(j) = sqrt((two*iix)*mybetax) * cos_mb(phix)
            myxp(j) = (-one*sqrt((two*iix)/mybetax)) * (sin_mb(phix) + myalphax * cos_mb(phix))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nominal bunches centered in the aperture - can't apply rejection sampling. return with error
         else if( mynex.eq.zero.and.myney.eq.zero ) then
           write(lout,*) "Stop in makedis_coll. attempting to use halo type 3 with Gaussian dist. "
           call prror(-1)
         else
           write(lout,*) "Error - beam parameters not correctly set!"
         end if

         myp(j) = myenom
         mys(j) = zero

      end do

      return
end subroutine makedis_coll

!========================================================================
!
! SR, 09-05-2005: Add the energy spread and the finite bunch length.
!                 Gaussian distributions assumed
subroutine makedis_de( myalphax, myalphay, mybetax, mybetay, &
     &     myemitx0, myemity0, myenom, mynex, mdex, myney, mdey,        &
     &     myx, myxp, myy, myyp, myp, mys,                              &
     &     enerror,bunchlength)

!     Uses the old routine 'MAKEDIS' for the halo plane and adds the
!     transverse beam size in the other plane (matched distrubutions
!     are generated starting from thetwiss functions).
!     If 'mynex' and 'myney' are BOTH set to zero, nominal bunches
!     centred in the aperture centre are generated. (SR, 08-05-2005)

  use crcoall
  use mod_common, only : napx
  implicit none

  integer :: j
  real(kind=fPrec), allocatable :: myx(:) !(npart)
  real(kind=fPrec), allocatable :: myxp(:) !(npart)
  real(kind=fPrec), allocatable :: myy(:) !(npart)
  real(kind=fPrec), allocatable :: myyp(:) !(npart)
  real(kind=fPrec), allocatable :: myp(:) !(npart)
  real(kind=fPrec), allocatable :: mys(:) !(npart)

  real(kind=fPrec) myalphax,mybetax,myemitx0,myemitx,mynex,mdex, &
  &mygammax,myalphay,mybetay,myemity0,myemity,myney,mdey,mygammay,   &
  &xsigmax,ysigmay,myenom

  real(kind=fPrec) iix, iiy, phix, phiy
  real(kind=fPrec) enerror, bunchlength
  real(kind=fPrec) en_error, bunch_length
  real(kind=fPrec) long_cut
  real(kind=fPrec) a_st, b_st
  save
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!++  Calculate the gammas

  mygammax = (one+myalphax**2)/mybetax
  mygammay = (one+myalphay**2)/mybetay

!     Assign bunch length and dp/p depending on the energy
!     Check if the units in metres are correct!
!GRD      if ( myenom.eq.7e6 ) then
!GRD         en_error     = 1.129e-4
!GRD         bunch_length = 7.55e-2
!GRD      elseif ( myenom.eq.4.5e5 ) then
!GRD         en_error     = 3.06e-4
!GRD         bunch_length = 11.24e-2
!GRD      else
  en_error = enerror
  bunch_length = bunchlength

!GRD         write(lout,*)"Warning-Energy different from LHC inj or top!"
!GRD         write(lout,*)"  => 7TeV values of dp/p and bunch length used!"
!GRD      endif

  write(lout,*) "Generation of bunch with dp/p and length:"
  write(lout,*) "  RMS bunch length  = ", bunch_length
  write(lout,*) "  RMS energy spread = ", en_error

  do j=1, napx
    if((mynex.gt.zero).and.(myney.eq.zero)) then
      myemitx = myemitx0*(mynex+((two*real(rndm4()-half,fPrec))*mdex))**2
      xsigmax = sqrt(mybetax*myemitx)
      myx(j)  = xsigmax * sin_mb((two*pi)*real(rndm4(),fPrec))

      if (rndm4().gt.half) then
        myxp(j) = sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax
      else
        myxp(j) = -one*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax
      end if

      phiy = (two*pi)*real(rndm4(),fPrec)
      iiy = (-one*myemity0) * log_mb( real(rndm4(),fPrec) )
      myy(j) = sqrt((two*iiy)*mybetay) * cos_mb(phiy)
      myyp(j) = (-one*sqrt((two*iiy)/mybetay)) * (sin_mb(phiy) + myalphay * cos_mb(phiy))

    else if( mynex.eq.zero.and.myney.gt.zero ) then
      myemity = myemity0*(myney+((two*real(rndm4()-half,fPrec))*mdey))**2
      ysigmay = sqrt(mybetay*myemity)
      myy(j)   = ysigmay * sin_mb((two*pi)*real(rndm4(),fPrec))

      if(rndm4().gt.half) then
        myyp(j) = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay
      else
        myyp(j) = -one*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay
      end if

      phix = (two*pi)*real(rndm4(),fPrec)
      iix = (-one*myemitx0) * log_mb( real(rndm4(),fPrec) )
      myx(j) = sqrt((two*iix)*mybetax) * cos_mb(phix)
      myxp(j) = (-one*sqrt((two*iix)/mybetax)) * (sin_mb(phix) + myalphax * cos_mb(phix))

    else if( mynex.eq.zero.and.myney.eq.zero ) then
      phix = (two*pi)*real(rndm4(),fPrec)
      iix = (-one*myemitx0) * log_mb( real(rndm4(),fPrec) )
      myx(j) = sqrt((two*iix)*mybetax) * cos_mb(phix)
      myxp(j) = (-one*sqrt((two*iix)/mybetax)) * (sin_mb(phix) + myalphax * cos_mb(phix))
      phiy = (two*pi)*real(rndm4(),fPrec)
      iiy = (-one*myemity0) * log_mb( real(rndm4(),fPrec) )
      myy(j) = sqrt((two*iiy)*mybetay) * cos_mb(phiy)
      myyp(j) = (-one*sqrt((two*iiy)/mybetay)) * (sin_mb(phiy) + myalphay * cos_mb(phiy))
    else
      write(lout,*) "Error - beam parameters not correctly set!"
    end if
  end do

! SR, 11-08-2005 For longitudinal phase-space, add a cut at 2 sigma
!++   1st: generate napxnumbers within the chose cut
  long_cut = 2
  j = 1
  do while (j.le.napx)
    a_st = ran_gauss(five)
    b_st = ran_gauss(five)

    do while ((a_st**2+b_st**2).gt.long_cut**2)
      a_st = ran_gauss(five)
      b_st = ran_gauss(five)
    end do

    mys(j) = a_st
    myp(j) = b_st
    j = j + 1
  end do

!++   2nd: give the correct values
  do j=1,napx
    myp(j) = myenom * (one + myp(j) * en_error)
    mys(j) = bunch_length * mys(j)
  end do

  return
end subroutine makedis_de


!========================================================================
subroutine readdis(filename_dis,myx,myxp,myy,myyp,myp,mys)
!
!     SR, 09-08-2005
!     Format for the input file:
!               x, y   -> [ m ]
!               xp, yp -> [ rad ]
!               s      -> [ mm ]
!               DE     -> [ MeV ]

  use crcoall
  use parpro
  use string_tools
  use mod_common, only : napx
  implicit none

  integer :: j
  real(kind=fPrec), allocatable :: myx(:) !(npart)
  real(kind=fPrec), allocatable :: myxp(:) !(npart)
  real(kind=fPrec), allocatable :: myy(:) !(npart)
  real(kind=fPrec), allocatable :: myyp(:) !(npart)
  real(kind=fPrec), allocatable :: myp(:) !(npart)
  real(kind=fPrec), allocatable :: mys(:) !(npart)

  character(len=80)   filename_dis

  integer stat

  character(len=mInputLn) inLine
  character(len=:), allocatable :: lnSplit(:)
  integer nSplit
  logical spErr

  save

  write(lout,"(a)") "COLL> Reading input bunch from file '"//filename_dis//"'"

  call f_requestUnit(filename_dis, filename_dis_unit)
  open(unit=filename_dis_unit, file=filename_dis, iostat=stat,status="OLD",action="read") !was 53
  if(stat.ne.0)then
    write(lout,"(a)")    "COLL> ERROR Subroutine readdis: Could not open the file."
    write(lout,"(a,i0)") "COLL>       Got iostat=",stat
    goto 20
  end if

  do j=1,napx
    read(filename_dis_unit,"(a)",end=10,err=20) inLine
    call chr_split(inLine, lnSplit, nSplit, spErr)
    if(spErr) then
      write(lout,"(a)") "COLL> ERROR Failed to parse input line from particle distribution file."
      goto 20
    end if
    if(nSplit /= 6) then
      write(lout,"(a)") "COLL> ERROR Expected 6 values per line in particle distribution file."
      goto 20
    end if
    call chr_cast(lnSplit(1),myx(j), spErr)
    call chr_cast(lnSplit(2),myxp(j),spErr)
    call chr_cast(lnSplit(3),myy(j), spErr)
    call chr_cast(lnSplit(4),myyp(j),spErr)
    call chr_cast(lnSplit(5),mys(j), spErr)
    call chr_cast(lnSplit(6),myp(j), spErr)
    if(spErr) then
      write(lout,"(a)") "COLL> ERROR Failed to parse value from particle distribution file."
      goto 20
    end if
  end do

  !TODO: Double-check that this is an OK way of dealing with reading less-than-expected particles from the file
 10   napx = j - 1
  write(lout,"(a,i0)") "COLL> Number of particles read from the file = ",napx

  close(filename_dis_unit)

  return

 20   continue

  write(lout,"(a)") "COLL> I/O Error on Unit 53 in subroutine readdis"
  call prror(-1)

end subroutine readdis

!========================================================================
!
subroutine readdis_norm(filename_dis,  myalphax, myalphay, mybetax, mybetay, &
 &           myemitx, myemity, myenom, myx, myxp, myy, myyp, myp, mys, enerror, bunchlength)
!     Format for the input file:
!               x, y   -> [ sigma ]
!               xp, yp -> [ sigma ]
!               s      -> [ sigma ]
!               DE     -> [ sigma ]

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use string_tools
  implicit none

  integer :: j
  real(kind=fPrec), allocatable :: myx(:) !(npart)
  real(kind=fPrec), allocatable :: myxp(:) !(npart)
  real(kind=fPrec), allocatable :: myy(:) !(npart)
  real(kind=fPrec), allocatable :: myyp(:) !(npart)
  real(kind=fPrec), allocatable :: myp(:) !(npart)
  real(kind=fPrec), allocatable :: mys(:) !(npart)

  real(kind=fPrec) myalphax, myalphay
  real(kind=fPrec) mybetax,myemitx
  real(kind=fPrec) mybetay,myemity,myenom

  character(len=80)   filename_dis
  real(kind=fPrec) enerror, bunchlength

  integer stat

  real(kind=fPrec) normx, normy, normxp, normyp, normp, norms
  real(kind=fPrec) myemitz

  character(len=mInputLn) inLine
  character(len=:), allocatable :: lnSplit(:)
  integer nSplit
  logical spErr

  if (iclo6.eq.0) then
    write(lout,"(a)") "COLL> ERROR DETECTED: Incompatible flag           "
    write(lout,"(a)") "COLL> in line 2 of the TRACKING block             "
    write(lout,"(a)") "COLL> of fort.3 for calculating the closed orbit  "
    write(lout,"(a)") "COLL> (iclo6 must not be =0). When using an input "
    write(lout,"(a)") "COLL> distribution in normalized coordinates for  "
    write(lout,"(a)") "COLL> collimation the closed orbit is needed for a"
    write(lout,"(a)") "COLL> correct TAS matrix for coordinate transform."
    call prror(-1)
  endif

  write(lout,"(a)") "COLL> Reading input bunch from file '"//filename_dis//"'"

  call f_requestUnit(filename_dis, filename_dis_unit)
  open(unit=filename_dis_unit, file=filename_dis, iostat=stat, status="OLD",action="read") !was 53
  if(stat.ne.0)then
    write(lout,"(a)")    "COLL> ERROR Subroutine readdis: Could not open the file."
    write(lout,"(a,i0)") "COLL>       Got iostat=",stat
    goto 20
  end if

  do j=1,napx
    read(filename_dis_unit,"(a)",end=10,err=20) inLine
    call chr_split(inLine, lnSplit, nSplit, spErr)
    if(spErr) then
      write(lout,"(a)") "COLL> ERROR Failed to parse input line from particle distribution file."
      goto 20
    end if
    if(nSplit /= 6) then
      write(lout,"(a)") "COLL> ERROR Expected 6 values per line in particle distribution file."
      goto 20
    end if
    call chr_cast(lnSplit(1),normx, spErr)
    call chr_cast(lnSplit(2),normxp,spErr)
    call chr_cast(lnSplit(3),normy, spErr)
    call chr_cast(lnSplit(4),normyp,spErr)
    call chr_cast(lnSplit(5),norms, spErr)
    call chr_cast(lnSplit(6),normp, spErr)
    if(spErr) then
      write(lout,"(a)") "COLL> ERROR Failed to parse value from particle distribution file."
      goto 20
    end if
! A normalized distribution with x,xp,y,yp,z,zp is read and
! transformed with the TAS matrix T , which is the transformation matrix
! from normalized to physical coordinates it is scaled with the geometric
! emittances in diag matrix S. x = T*S*normx
! units of TAS matrix # m,rad,m,rad,m,1
! The collimation coordinates/units are
! x[m], x'[rad], y[m], y'[rad]$, sig[mm], dE [MeV].

!         write(lout,*) " myenom [MeV]= ",myenom
!         write(lout,*) " myemitx [m]= ",myemitx
!         write(lout,*) " myemity [m]= ",myemity
!         write(lout,*) " bunchlength [mm]= ",bunchlength
!         write(lout,*) " enerror = ",enerror

         !convert bunchlength from [mm] to [m]
         ! enerror is the energy spread
    myemitz  = bunchlength * c1m3 * enerror


! scaling the TAS matrix entries of the longitudinal coordinate. tas(ia,j,k)  ia=the particle for which the tas was written

    myx(j)   =                            &
     &     normx  * sqrt(myemitx)*tas(1,1,1) + &
     &     normxp * sqrt(myemitx)*tas(1,1,2) + &
     &     normy  * sqrt(myemity)*tas(1,1,3) + &
     &     normyp * sqrt(myemity)*tas(1,1,4) + &
     &     norms  * sqrt(myemitz)*tas(1,1,5) + &
     &     normp  * sqrt(myemitz)*c1m3*tas(1,1,6)

    myxp(j)  =                            &
     &     normx  * sqrt(myemitx)*tas(1,2,1) + &
     &     normxp * sqrt(myemitx)*tas(1,2,2) + &
     &     normy  * sqrt(myemity)*tas(1,2,3) + &
     &     normyp * sqrt(myemity)*tas(1,2,4) + &
     &     norms  * sqrt(myemitz)*tas(1,2,5) + &
     &     normp  * sqrt(myemitz)*c1m3*tas(1,2,6)

    myy(j)   =                            &
     &     normx  * sqrt(myemitx)*tas(1,3,1) + &
     &     normxp * sqrt(myemitx)*tas(1,3,2) + &
     &     normy  * sqrt(myemity)*tas(1,3,3) + &
     &     normyp * sqrt(myemity)*tas(1,3,4) + &
     &     norms  * sqrt(myemitz)*tas(1,3,5) + &
     &     normp  * sqrt(myemitz)*c1m3*tas(1,3,6)

    myyp(j)  =                            &
     &     normx  * sqrt(myemitx)*tas(1,4,1) + &
     &     normxp * sqrt(myemitx)*tas(1,4,2) + &
     &     normy  * sqrt(myemity)*tas(1,4,3) + &
     &     normyp * sqrt(myemity)*tas(1,4,4) + &
     &     norms  * sqrt(myemitz)*tas(1,4,5) + &
     &     normp  * sqrt(myemitz)*c1m3*tas(1,4,6)

    mys(j)   =                            &
     &     normx  * sqrt(myemitx)*tas(1,5,1) + &
     &     normxp * sqrt(myemitx)*tas(1,5,2) + &
     &     normy  * sqrt(myemity)*tas(1,5,3) + &
     &     normyp * sqrt(myemity)*tas(1,5,4) + &
     &     norms  * sqrt(myemitz)*tas(1,5,5) + &
     &     normp  * sqrt(myemitz)*c1m3*tas(1,5,6)

    myp(j)   =                                    &
     &     normx  * sqrt(myemitx)*c1e3*tas(1,6,1) + &
     &     normxp * sqrt(myemitx)*c1e3*tas(1,6,2) + &
     &     normy  * sqrt(myemity)*c1e3*tas(1,6,3) + &
     &     normyp * sqrt(myemity)*c1e3*tas(1,6,4) + &
     &     norms  * sqrt(myemitz)*c1e3*tas(1,6,5) + &
     &     normp  * sqrt(myemitz)*tas(1,6,6)

! add the momentum
! convert to canonical variables
! dE/E with unit [1] from the closed orbit is added
!For the 4D coordinates the closed orbit
! will be added by SixTrack itself later on.
     myxp(j)  = myxp(j)*(one+myp(j)+clop6v(3,1))
     myyp(j)  = myyp(j)*(one+myp(j)+clop6v(3,1))
! unit conversion for collimation [m] to [mm]
     mys(j)   = mys(j)*c1e3
     myp(j)   = myenom*(one+myp(j))

  end do

  !TODO: Double-check that this is an OK way of dealing with reading less-than-expected particles from the file
10   napx = j - 1
  write(lout,"(a,i0)") "COLL> Number of particles read from the file = ",napx

  close(filename_dis_unit)
  return

20 continue
   write(lout,"(a)") "COLL> I/O Error on Unit 53 in subroutine readdis"
   call prror(-1)

end subroutine readdis_norm


!========================================================================
!
subroutine makedis_radial( myalphax, myalphay, mybetax,      &
     &mybetay, myemitx0, myemity0, myenom, nr, ndr, myx, myxp, myy, myyp, myp, mys)

  use mod_common, only : napx
  use crcoall
  implicit none

  integer :: j
  real(kind=fPrec), allocatable :: myx(:) !(npart)
  real(kind=fPrec), allocatable :: myxp(:) !(npart)
  real(kind=fPrec), allocatable :: myy(:) !(npart)
  real(kind=fPrec), allocatable :: myyp(:) !(npart)
  real(kind=fPrec), allocatable :: myp(:) !(npart)
  real(kind=fPrec), allocatable :: mys(:) !(npart)

  real(kind=fPrec) myalphax,mybetax,myemitx0,myemitx,mynex,mdex, &
  &mygammax,myalphay,mybetay,myemity0,myemity,myney,mdey,mygammay,   &
  &xsigmax,ysigmay,myenom,nr,ndr

  save
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!++  Calculate the gammas

  mygammax = (one+myalphax**2)/mybetax
  mygammay = (one+myalphay**2)/mybetay

!++  Number of points and generate distribution

  mynex = nr/sqrt(two)
  mdex = ndr/sqrt(two)
  myney = nr/sqrt(two)
  mdey = ndr/sqrt(two)

  write(lout,*)
  write(lout,*) 'Generation of particle distribution Version 2'
  write(lout,*)
  write(lout,*) 'This routine generates particles in that are fully'
  write(lout,*) 'correlated between X and Y.'
  write(lout,*)

  write(outlun,*)
  write(outlun,*) 'Generation of particle distribution Version 2'
  write(outlun,*)
  write(outlun,*) 'This routine generates particles in that are fully'
  write(outlun,*) 'correlated between X and Y.'
  write(outlun,*)
  write(outlun,*)
  write(outlun,*) 'INFO>  Number of particles   = ', napx
  write(outlun,*) 'INFO>  Av number of x sigmas = ', mynex
  write(outlun,*) 'INFO>  +- spread in x sigmas = ', mdex
  write(outlun,*) 'INFO>  Av number of y sigmas = ', myney
  write(outlun,*) 'INFO>  +- spread in y sigmas = ', mdey
  write(outlun,*) 'INFO>  Nominal beam energy   = ', myenom
  write(outlun,*) 'INFO>  Sigma_x0 = ', sqrt(mybetax*myemitx0)
  write(outlun,*) 'INFO>  Sigma_y0 = ', sqrt(mybetay*myemity0)
  write(outlun,*)

  do while (j.lt.napx)

    j = j + 1
    myemitx = myemitx0*(mynex + ((two*real(rndm4()-half,fPrec))*mdex) )**2  !hr09
    xsigmax = sqrt(mybetax*myemitx)
    myx(j)  = xsigmax * sin_mb((two*pi)*real(rndm4(),fPrec))              !hr09

    if (rndm4().gt.half) then
      myxp(j) =      sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax !hr09
    else
      myxp(j) = -one*sqrt(myemitx/mybetax-myx(j)**2/mybetax**2)-(myalphax*myx(j))/mybetax !hr09
    endif

    myemity = myemity0*(myney + ((two*real(rndm4()-half,fPrec))*mdey) )**2  !hr09
    ysigmay = sqrt(mybetay*myemity)
    myy(j)  = ysigmay * sin_mb((two*pi)*real(rndm4(),fPrec))          !hr09

    if (rndm4().gt.half) then
      myyp(j)  = sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay      !hr09
    else
      myyp(j)  = -one*sqrt(myemity/mybetay-myy(j)**2/mybetay**2)-(myalphay*myy(j))/mybetay !hr09
    endif

!APRIL2005
    myp(j)   = myenom
!        if(j.eq.1) then
!          myp(j)   = myenom*(1-0.05)
!!       do j=2,mynp
!        else
!          myp(j) = myp(1) + (j-1)*2d0*0.05*myenom/(mynp-1)
!        endif
!APRIL2005
    mys(j)   = zero

!++  Dangerous stuff, just for the moment
!
!        IF ( (.NOT. (Y(j).LT.-.008e-3 .AND. YP(j).LT.0.1e-3 .AND.
!     1               YP(j).GT.0.0) ) .AND.
!     2       (.NOT. (Y(j).GT..008e-3 .AND. YP(j).GT.-0.1e-3 .AND.
!     3               YP(j).LT.0.0) ) ) THEN
!          J = J - 1
!        ENDIF
!
  end do

  return
end subroutine makedis_radial

!>
!! \brief The routine makes an initial Gaussian distribution
!!
!!     Uses the old routine 'MAKEDIS' for the halo plane and adds the\n
!!     transverse beam size in the other plane (matched distrubutions\n
!!     are generated starting from the twiss functions).\n
!!     If 'mynex' and 'myney' are BOTH set to zero, nominal bunches\n
!!     centred in the aperture centre are generated. (SR, 08-05-2005)
!!
!!     YIL EDIT 2010: particle 0 is always on orbit...
!!
!! @author Javier Barranco <jbarranc@cern.ch>
!! @param myalphax
!! @param myalphay
!! @param mybetax
!! @param mybetay
!! @param myemitx0
!! @param myemity0
!! @param myenom
!! @param mynex
!! @param mdex
!! @param myney
!! @param mdey
!! @param myx
!! @param myxp
!! @param myy
!! @param myyp
!! @param myp
!! @param mys
!! @param enerror
!! @param bunchlength
!!
!! @date Last modified: 06. August 2009
!! @see ran_gauss
!!
!<
subroutine makedis_ga( myalphax, myalphay, mybetax, mybetay, myemitx0, myemity0, myenom, mynex, mdex, myney, mdey, &
 &  myx, myxp, myy, myyp, myp, mys, enerror, bunchlength )

  use crcoall
  use parpro
  use mod_common_track
  use mod_common, only : napx
  implicit none

  integer :: j
  real(kind=fPrec), allocatable :: myx(:) !(npart)
  real(kind=fPrec), allocatable :: myxp(:) !(npart)
  real(kind=fPrec), allocatable :: myy(:) !(npart)
  real(kind=fPrec), allocatable :: myyp(:) !(npart)
  real(kind=fPrec), allocatable :: myp(:) !(npart)
  real(kind=fPrec), allocatable :: mys(:) !(npart)

  real(kind=fPrec) myalphax,mybetax,myemitx0,myemitx,mynex,mdex, &
  &mygammax,myalphay,mybetay,myemity0,myemity,myney,mdey,mygammay,   &
  &xsigmax,ysigmay,myenom

  real(kind=fPrec) enerror, bunchlength
  real(kind=fPrec) en_error, bunch_length

  real(kind=fPrec) long_cut
  real(kind=fPrec) a_st, b_st
  integer startpar

  save

!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!++  Generate random distribution, assuming optical parameters at IP1
!
!++  Calculate the gammas

  mygammax = (one+myalphax**2)/mybetax
  mygammay = (one+myalphay**2)/mybetay
  en_error = enerror
  bunch_length = bunchlength

  write (lout,*) "Generation of bunch with dp/p and length:"
  write (lout,*) "  RMS bunch length  = ", bunch_length
  write (lout,*) "  RMS energy spread = ", en_error
! JBG August 2007
  write (lout,*)
  write (lout,*) "   ***STEP 1 for Gaussian Beam***"
  write (lout,*)
  write (lout,*) "   Beam generated with 5 sigma cut"
  write (lout,*)
  write (lout,*) "  Parameters used for Distribution Generation"
  write (lout,*) "  BetaX =", mybetax
  write (lout,*) "  BetaY =", mybetay
  write (lout,*) "  EmittanceX =", myemitx0
  write (lout,*) "  EmittanceY =", myemity0
  write (lout,*)

  startpar=1

#ifdef BEAMGAS
! YIL July 2010 first particle on orbit
!  initial xangle (if any) is not
!  yet applied at this point...
!  so we can set all to 0.
  startpar=2
  myx(1)  = zero     !hr13
  myy(1)  = zero     !hr13
  myxp(1) = zero     !hr13
  myyp(1) = zero     !hr13
  myp(1)  = myenom
  mys(1)  = zero
!YIL end edit July 2010
#endif

  do j=startpar, napx
! JBG July 2007
! Option added for septum studies

    myemitx = myemitx0
    xsigmax = sqrt(mybetax*myemitx)
    myx(j)  = xsigmax * ran_gauss(mynex)
    myxp(j) = ran_gauss(mynex)*sqrt(myemitx/mybetax)-((myalphax*myx(j))/mybetax)  !hr13
    myemity = myemity0
    ysigmay = sqrt(mybetay*myemity)
    myy(j)  = ysigmay * ran_gauss(myney)
    myyp(j) = ran_gauss(myney)*sqrt(myemity/mybetay)-((myalphay*myy(j))/mybetay)  !hr13
  end do

! SR, 11-08-2005 For longitudinal phase-space, add a cut at 2 sigma
!++   1st: generate napxnumbers within the chosen cut

  long_cut = 2
  j = startpar

  do while (j.le.napx)
    a_st = ran_gauss(five)
    b_st = ran_gauss(five)

    do while ((a_st*a_st+b_st*b_st).gt.long_cut*long_cut)
      a_st = ran_gauss(five)
      b_st = ran_gauss(five)
    end do

    mys(j) = a_st
    myp(j) = b_st
    j = j + 1
  end do

!++   2nd: give the correct values
  do j=startpar,napx
    myp(j) = myenom * (one + myp(j) * en_error)
    mys(j) = bunch_length * mys(j)
  end do

  return
end subroutine makedis_ga

function rndm4()

  implicit none

  integer len, in
  real(kind=fPrec) rndm4, a
  save IN,a
  parameter ( len =  30000 )
  dimension a(len)
  data in/1/

  if( in.eq.1 ) then
    call ranlux(a,len)
!    call ranecu(a,len,-1)
    rndm4=a(1)
    in=2
  else
    rndm4=a(in)
    in=in+1
    if(in.eq.len+1)in=1
  endif

  return

end function rndm4


!ccccccccccccccccccccccccccccccccccccccc
!-TW-01/2007
! function rndm5(irnd) , irnd = 1 will reset
! inn counter => enables reproducible set of
! random unmbers
!cccccccccccccccccccccccccccccccccc
!
function rndm5(irnd)

  implicit none

  integer len, inn, irnd
  real(kind=fPrec) rndm5,a
  save

  parameter( len =  30000 )
  dimension a(len)
  data inn/1/
!
! reset inn to 1 enable reproducible random numbers
  if( irnd .eq. 1) inn = 1

  if( inn.eq.1 ) then
    call ranlux(a,len)
!     call ranecu(a,len,-1)
    rndm5=a(1)
    inn=2
  else
    rndm5=a(inn)
    inn=inn+1
    if(inn.eq.len+1)inn=1
  end if

  return
end function rndm5

!ccccccccccccccccccccccccccccccccccccccc
real(kind=fPrec) function myran_gauss(cut)
!*********************************************************************
!
! myran_gauss - will generate a normal distribution from a uniform
!     distribution between [0,1].
!     See "Communications of the ACM", V. 15 (1972), p. 873.
!
!     cut - real(kind=fPrec) - cut for distribution in units of sigma
!     the cut must be greater than 0.5
!
!     changed rndm4 to rndm5(irnd) and defined flag as true
!
!*********************************************************************

  implicit none

  logical flag
  real(kind=fPrec) x, u1, u2, twopi, r,cut
  save

  flag = .true. !Does this initialize only once, or is it executed every pass?
                !See ran_gauss(cut)

  twopi=eight*atan_mb(one)
1 if (flag) then
    r = real(rndm5(0),fPrec)
    r = max(r, half**32)
    r = min(r, one-half**32)
    u1 = sqrt(-two*log_mb( r ))
    u2 = real(rndm5(0),fPrec)
    x = u1 * cos_mb(twopi*u2)
  else
     x = u1 * sin_mb(twopi*u2)
  endif

  flag = .not. flag

!     cut the distribution if cut > 0.5
  if (cut .gt. half .and. abs(x) .gt. cut) goto 1

  myran_gauss = x
  return
end function myran_gauss


!cccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine funlxp (func,xfcum,x2low,x2high)
!         F. JAMES,   Sept, 1994
!
!         Prepares the user function FUNC for FUNLUX
!         Inspired by and mostly copied from FUNPRE and FUNRAN
!         except that
!    1. FUNLUX uses RANLUX underneath,
!    2. FUNLXP expands the first and last bins to cater for
!              functions with long tails on left and/or right,
!    3. FUNLXP calls FUNPCT to do the actual finding of percentiles.
!    4. both FUNLXP and FUNPCT use RADAPT for Gaussian integration.
!
      use crcoall
      implicit none
      external func
      integer ifunc,ierr
      real(kind=fPrec) x2high,x2low,xfcum,rteps,xhigh,xlow,xrange,uncert,x2,tftot1,x3,tftot2,func
      dimension xfcum(200)
      parameter (rteps=0.0002)
      save ifunc
      data ifunc/0/
      ifunc = ifunc + 1
!         FIND RANGE WHERE FUNCTION IS NON-ZERO.
      call funlz(func,x2low,x2high,xlow,xhigh)
      xrange = xhigh-xlow
      if(xrange .le. 0.)  then
        write(lout,'(A,2G15.5)') ' FUNLXP finds function range .LE.0',xlow,xhigh
        go to 900
      endif
      call radapt(func,xlow,xhigh,1,rteps,zero,tftot ,uncert)
!      WRITE(6,1003) IFUNC,XLOW,XHIGH,TFTOT
! 1003 format(' FUNLXP: integral of USER FUNCTION', i3,' from ',e12.5,' to ',e12.5,' is ',e14.6)
!
!      WRITE (6,'(A,A)') ' FUNLXP preparing ',
!     + 'first the whole range, then left tail, then right tail.'
      call funpct(func,ifunc,xlow,xhigh,xfcum,1,99,tftot,ierr)
      if (ierr .gt. 0)  go to 900
      x2 = xfcum(3)
      call radapt(func,xlow,x2,1,rteps,zero,tftot1 ,uncert)
      call funpct(func,ifunc,xlow,x2 ,xfcum,101,49,tftot1,ierr)
      if (ierr .gt. 0)  go to 900
      x3 = xfcum(98)
      call radapt(func,x3,xhigh,1,rteps,zero,tftot2 ,uncert)
      call funpct(func,ifunc,x3,xhigh,xfcum,151,49,tftot2,ierr)
      if (ierr .gt. 0)  go to 900
!      WRITE(6,1001) IFUNC,XLOW,XHIGH
! 1001 format(' FUNLXP has prepared USER FUNCTION', i3, ' between',g12.3,' and',g12.3,' for FUNLUX')

      return
  900 continue
      write(lout,*) ' Fatal error in FUNLXP. FUNLUX will not work.'
end subroutine funlxp

subroutine funpct(func,ifunc,xlow,xhigh,xfcum,nlo,nbins,tftot,ierr)
!        Array XFCUM is filled from NLO to NLO+NBINS, which makes
!        the number of values NBINS+1, or the number of bins NBINS
      use crcoall
      implicit none
      external func
      integer ierr,nbins,nlo,ifunc,nz,ibin,maxz,iz,nitmax,ihome
      real(kind=fPrec) tftot,xhigh,xlow,func,xfcum,rteps,tpctil,tz,tzmax,x,f,tcum,  &
     &x1,f1,dxmax,fmin,fminz,xincr,tincr,xbest,dtbest,tpart,x2,precis,  &
     &refx,uncert,tpart2,dtpar2,dtabs,aberr
      dimension xfcum(*)
      parameter (rteps=0.005, nz=10, maxz=20, nitmax=6,precis=1e-6)
!      DOUBLE PRECISION TPCTIL, TZ, TCUM, XINCR, DTABS,
!     &  TINCR, TZMAX, XBEST, DTBEST, DTPAR2
!
      ierr = 0
      if (tftot .le. 0.) go to 900
      tpctil = tftot/real(nbins)                                         !hr09
      tz = tpctil/real(nz)
      tzmax = tz * 2.
      xfcum(nlo) = xlow
      xfcum(nlo+nbins) = xhigh
      x = xlow
      f = func(x)
      if (f .lt. 0.) go to 900
!         Loop over percentile bins
      do 600 ibin = nlo, nlo+nbins-2
      tcum = 0.
      x1 = x
      f1 = f
      dxmax = (xhigh -x) / nz
      fmin = tz/dxmax
      fminz = fmin
!         Loop over trapezoids within a supposed percentil
      do 500 iz= 1, maxz
      xincr = tz/max(f1,fmin,fminz)
  350 x = x1 + xincr
      f = func(x)
      if (f .lt. 0.) go to 900
      tincr = ((x-x1) * 0.5) * (f+f1)                                    !hr09
      if (tincr .lt. tzmax) go to 370
      xincr = xincr * 0.5
      go to 350
  370 continue
      tcum = tcum + tincr
      if (tcum .ge. tpctil*0.99) go to 520
      fminz = (tz*f)/ (tpctil-tcum)                                      !hr09
      f1 = f
      x1 = x
  500 continue
      write(lout,*) ' FUNLUX:  WARNING. FUNPCT fails trapezoid.'
!         END OF TRAPEZOID LOOP
!         Adjust interval using Gaussian integration with
!             Newton corrections since F is the derivative
  520 continue
      x1 = xfcum(ibin)
      xbest = x
      dtbest = tpctil
      tpart = tpctil
!         Allow for maximum NITMAX more iterations on RADAPT
      do 550 ihome= 1, nitmax
  535 xincr = (tpctil-tpart) / max(f,fmin)
      x = xbest + xincr
      x2 = x
        if (ihome .gt. 1 .and. x2 .eq. xbest) then
        write(lout,'(A,G12.3)') ' FUNLUX: WARNING from FUNPCT: insufficient precision at X=',x
        go to 580
        endif
      refx = abs(x)+precis
      call radapt(func,x1,x2,1,rteps,zero,tpart2,uncert)
      dtpar2 = tpart2-tpctil
      dtabs = abs(dtpar2)
      if(abs(xincr)/refx .lt. precis) goto 545
      if(dtabs .lt. dtbest) goto 545
      xincr = xincr * 0.5
      goto 535
  545 dtbest = dtabs
      xbest = x
      tpart = tpart2
      f = func(x)
      if(f .lt. 0.) goto 900
      if(dtabs .lt. rteps*tpctil) goto 580
  550 continue
      write(lout,'(A,I4)') ' FUNLUX: WARNING from FUNPCT: cannot converge, bin',ibin

  580 continue
      xincr = (tpctil-tpart) / max(f,fmin)
      x = xbest + xincr
      xfcum(ibin+1) = x
      f = func(x)
      if(f .lt. 0.) goto 900
  600 continue
!         END OF LOOP OVER BINS
      x1 = xfcum((nlo+nbins)-1)                                          !hr09
      x2 = xhigh
      call radapt(func,x1,x2,1,rteps,zero,tpart ,uncert)
      aberr = abs(tpart-tpctil)/tftot
!      WRITE(6,1001) IFUNC,XLOW,XHIGH
      if(aberr .gt. rteps)  write(lout,1002) aberr
      return
  900 write(lout,1000) x,f
      ierr = 1
      return
 1000 format(/' FUNLUX fatal error in FUNPCT: function negative:'/      &
&,' at X=',e15.6,', F=',e15.6/)
! 1001 FORMAT(' FUNPCT has prepared USER FUNCTION',I3,
!     + ' between',G12.3,' and',G12.3,' for FUNLUX.')
 1002 format(' WARNING: Relative error in cumulative distribution may be as big as',f10.7)

end subroutine funpct

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine funlux(array,xran,len)
!         Generation of LEN random numbers in any given distribution,
!         by 4-point interpolation in the inverse cumulative distr.
!         which was previously generated by FUNLXP
      implicit none
      integer len,ibuf,j,j1
      real(kind=fPrec) array,xran,gap,gapinv,tleft,bright,gaps,gapins,x,p,a,b
      dimension array(200)
      dimension xran(len)
!  Bin width for main sequence, and its inverse
      parameter (gap= 1./99.,  gapinv=99.)
!  Top of left tail, bottom of right tail (each tail replaces 2 bins)
      parameter (tleft= 2./99.,bright=97./99.)
!  Bin width for minor sequences (tails), and its inverse
      parameter (gaps=tleft/49.,  gapins=1./gaps)
!
!   The array ARRAY is assumed to have the following structure:
!        ARRAY(1-100) contains the 99 bins of the inverse cumulative
!                     distribution of the entire function.
!        ARRAY(101-150) contains the 49-bin blowup of main bins
!                       1 and 2 (left tail of distribution)
!        ARRAY(151-200) contains the 49-bin blowup of main bins
!                       98 and 99 (right tail of distribution)
!
      x = zero ! -Wmaybe-uninitialized
      call ranlux(xran,len)
!      call ranecu(xran,len,-1)

      do 500 ibuf= 1, len
      x = xran(ibuf)
      j = int(  x    *gapinv) + 1
      if (j .lt. 3)  then
         j1 = int( x *gapins)
             j = j1 + 101
             j = max(j,102)
             j = min(j,148)
         p = (   x -gaps*real(j1-1)) * gapins                            !hr09
         a = (p+1.0) * array(j+2) - (p-2.0)*array(j-1)
         b = (p-1.0) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.0))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5  !hr09
      else if (j .gt. 97)  then
         j1 = int((x-bright)*gapins)
             j = j1 + 151
             j = max(j,152)
             j = min(j,198)
         p = ((x -bright) -gaps*(j1-1)) * gapins                         !hr09
         a = (p+1.0) * array(j+2) - (p-2.0)*array(j-1)
         b = (p-1.0) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.0))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5  !hr09
      else
!      J = MAX(J,2)
!      J = MIN(J,98)
         p = (   x -gap*real(j-1)) * gapinv                              !hr09
         a = (p+1.) * array(j+2) - (p-2.)*array(j-1)
         b = (p-1.) * array(j) - p * array(j+1)
      xran(ibuf) = ((a*p)*(p-1.))*0.16666667 + ((b*(p+1.))*(p-2.))*0.5   !hr09
      endif
  500 continue
      tftot = x
      return
end subroutine funlux

subroutine funlz(func,x2low,x2high,xlow,xhigh)
! FIND RANGE WHERE FUNC IS NON-ZERO.
! WRITTEN 1980, F. JAMES
! MODIFIED, NOV. 1985, TO FIX BUG AND GENERALIZE
! TO FIND SIMPLY-CONNECTED NON-ZERO REGION (XLOW,XHIGH)
! ANYWHERE WITHIN THE GIVEN REGION (X2LOW,H2HIGH).
!    WHERE 'ANYWHERE' MEANS EITHER AT THE LOWER OR UPPER
!    EDGE OF THE GIVEN REGION, OR, IF IN THE MIDDLE,
!    COVERING AT LEAST 1% OF THE GIVEN REGION.
! OTHERWISE IT IS NOT GUARANTEED TO FIND THE NON-ZERO REGION.
! IF FUNCTION EVERYWHERE ZERO, FUNLZ SETS XLOW=XHIGH=0.
      use crcoall
      implicit none
      external func
      integer logn,nslice,i,k
      real(kind=fPrec) xhigh,xlow,x2high,x2low,func,xmid,xh,xl,xnew
      xlow = x2low
      xhigh = x2high
!         FIND OUT IF FUNCTION IS ZERO AT ONE END OR BOTH
      xmid = xlow
      if (func(xlow) .gt. 0.) go to 120
      xmid = xhigh
      if (func(xhigh) .gt. 0.)  go to 50
!         FUNCTION IS ZERO AT BOTH ENDS,
!         LOOK FOR PLACE WHERE IT IS NON-ZERO.
      do 30 logn= 1, 7
      nslice = 2**logn
      do 20 i= 1, nslice, 2
      xmid = xlow + (real(i) * (xhigh-xlow)) / real(nslice)              !hr09
      if (func(xmid) .gt. 0.)  go to 50
   20 continue
   30 continue
!         FALLING THROUGH LOOP MEANS CANNOT FIND NON-ZERO VALUE
      write(lout,554)
      write(lout,555) xlow, xhigh
      xlow = 0.
      xhigh = 0.
      go to 220
!
   50 continue
!         DELETE 'LEADING' ZERO RANGE
      xh = xmid
      xl = xlow
      do 70 k= 1, 20
      xnew = 0.5*(xh+xl)
      if (func(xnew) .eq. 0.) go to 68
      xh = xnew
      go to 70
   68 xl = xnew
   70 continue
      xlow = xl
      write(lout,555) x2low,xlow
  120 continue
      if (func(xhigh) .gt. 0.) go to 220
!         DELETE 'TRAILING' RANGE OF ZEROES
      xl = xmid
      xh = xhigh
      do 170 k= 1, 20
      xnew = 0.5*(xh+xl)
      if (func(xnew) .eq. 0.) go to 168
      xl = xnew
      go to 170
  168 xh = xnew
  170 continue
      xhigh = xh
      write(lout,555) xhigh, x2high
!
  220 continue
      return
  554 format('0CANNOT FIND NON-ZERO FUNCTION VALUE')
  555 format(' FUNCTION IS ZERO FROM X=',e12.5,' TO ',e12.5)
end subroutine funlz

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine radapt(f,a,b,nseg,reltol,abstol,res,err)

! RES = Estimated Integral of F from A to B,
! ERR = Estimated absolute error on RES.
! NSEG  specifies how the adaptation is to be done:
!    =0   means use previous binning,
!    =1   means fully automatic, adapt until tolerance attained.
!    =n>1 means first split interval into n equal segments,
!         then adapt as necessary to attain tolerance.
! The specified tolerances are:
!        relative: RELTOL ;  absolute: ABSTOL.
!    It stop s when one OR the other is satisfied, or number of
!    segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
!    can be set to zero, in which case only the other is used.

      implicit none

      external f
      integer nseg,ndim,nter,nsegd,i,iter,ibig
      real(kind=fPrec) err,res,abstol,reltol,b,a,xlo,xhi,tval,ters,te,root,xhib,bin,xlob,bige,hf,xnew,r1,f
      real(kind=fPrec) tvals,terss

      parameter (ndim=100)
      parameter (r1 = 1., hf = r1/2.)

      dimension xlo(ndim),xhi(ndim),tval(ndim),ters(ndim)
      save xlo,xhi,tval,ters,nter
      data nter /0/

      if(nseg .le. 0)  then
       if(nter .eq. 0) then
        nsegd=1
        go to 2
       endif
       tvals=zero
       terss=zero
       do 1 i = 1,nter
       call rgs56p(f,xlo(i),xhi(i),tval(i),te)
       ters(i)=te**2
       tvals=tvals+tval(i)                                         !hr09
       terss=terss+ters(i)
    1  continue
       root= sqrt(two*terss)                                      !hr09
       go to 9
      endif
      nsegd=min(nseg,ndim)
    2 xhib=a
      bin=(b-a)/real(nsegd,fPrec)                                              !hr09
      do 3 i = 1,nsegd
      xlo(i)=xhib
      xlob=xlo(i)
      xhi(i)=xhib+bin
      if(i .eq. nsegd) xhi(i)=b
      xhib=xhi(i)
      call rgs56p(f,xlob,xhib,tval(i),te)
      ters(i)=te**2
    3 continue
      nter=nsegd
      do 4 iter = 1,ndim
      tvals=tval(1)                                                !hr09
      terss=ters(1)                                                !hr09
      do 5 i = 2,nter
      tvals=tvals+tval(i)                                          !hr09
      terss=terss+ters(i)                                          !hr09
    5 continue
      root=sqrt(two*terss)                                       !hr09

      if(root .le. abstol .or. root .le. reltol*abs(tvals)) then
        goto 9
      end if

      if(nter .eq. ndim) go to 9
      bige=ters(1)
      ibig=1
      do 6 i = 2,nter
      if(ters(i) .gt. bige) then
       bige=ters(i)
       ibig=i
      endif
    6 continue
      nter=nter+1
      xhi(nter)=xhi(ibig)
      xnew=hf*(xlo(ibig)+xhi(ibig))
      xhi(ibig)=xnew
      xlo(nter)=xnew
      call rgs56p(f,xlo(ibig),xhi(ibig),tval(ibig),te)
      ters(ibig)=te**2
      call rgs56p(f,xlo(nter),xhi(nter),tval(nter),te)
      ters(nter)=te**2
    4 continue
    9 res=tvals                                                    !hr09
      err=root
      return
end subroutine radapt

!cccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine rgs56p(f,a,b,res,err)

  implicit none

  integer i
  real(kind=fPrec) err,res,b,a,f,w6,x6,w5,x5,rang,r1,hf
  real(kind=fPrec) e5,e6

  parameter (r1 = 1., hf = r1/2.)
  dimension x5(5),w5(5),x6(6),w6(6)

  data (x5(i),w5(i),i=1,5)                                          &
 &/4.6910077030668004e-02, 1.1846344252809454e-01,                  &
 &2.3076534494715846e-01, 2.3931433524968324e-01,                   &
 &5.0000000000000000e-01, 2.8444444444444444e-01,                   &
 &7.6923465505284154e-01, 2.3931433524968324e-01,                   &
 &9.5308992296933200e-01, 1.1846344252809454e-01/

  data (x6(i),w6(i),i=1,6)                                          &
 &/3.3765242898423989e-02, 8.5662246189585178e-02,                  &
 &1.6939530676686775e-01, 1.8038078652406930e-01,                   &
 &3.8069040695840155e-01, 2.3395696728634552e-01,                   &
 &6.1930959304159845e-01, 2.3395696728634552e-01,                   &
 &8.3060469323313225e-01, 1.8038078652406930e-01,                   &
 &9.6623475710157601e-01, 8.5662246189585178e-02/

  rang=b-a
  e5=zero
  e6=zero
  do i = 1,5
    e5=e5+dble(w5(i)*f(a+rang*x5(i)))                                  !hr09
    e6=e6+dble(w6(i)*f(a+rang*x6(i)))                                  !hr09
  end do

  e6=e6+dble(w6(6)*f(a+rang*x6(6)))
  res=real((dble(hf)*(e6+e5))*dble(rang))                            !hr09
  err=real(abs((e6-e5)*dble(rang)))                                  !hr09
  return
end subroutine rgs56p

!*********************************************************************
!
! Define INTEGER function MCLOCK that can differ from system to system
! For re-initializtion of random generator
!
!*********************************************************************
integer function mclock_liar( )
  use crcoall
  implicit none

  save

  integer    mclock
  integer    count_rate, count_max
  logical    clock_ok

!        MCLOCK_LIAR = MCLOCK()

  clock_ok = .true.

  if (clock_ok) then
    call system_clock( mclock, count_rate, count_max )
    if ( count_max .eq. 0 ) then
      clock_ok = .false.
      write(lout,"(a)") 'COLL> System Clock not present or not Responding'
      write(lout,"(a)") 'COLL> R.N.G. Reseed operation disabled.'
    endif
  endif

  mclock_liar = mclock

  return
end function mclock_liar


real(kind=fPrec) function ran_gauss(cut)
!*********************************************************************
!
! RAN_GAUSS - will generate a normal distribution from a uniform
!   distribution between [0,1].
!   See "Communications of the ACM", V. 15 (1972), p. 873.
!
! cut - real(kind=fPrec) - cut for distribution in units of sigma
!                the cut must be greater than 0.5
!
!*********************************************************************

  use crcoall
  use parpro
  implicit none

  logical flag
  DATA flag/.TRUE./
  real(kind=fPrec) x, u1, u2, twopi, r,cut

  save

  twopi=eight*atan_mb(one) !Why not 2*pi, where pi is in block "common"?
1 if (flag) then
    r = real(rndm4(),fPrec)
    r = max(r, half**32)
    r = min(r, one-half**32)
    u1 = sqrt(-two*log_mb( r ))
    u2 = real(rndm4(),fPrec)
    x = u1 * cos_mb(twopi*u2)
  else
    x = u1 * sin_mb(twopi*u2)
  endif

  flag = .not. flag

!  cut the distribution if cut > 0.5
  if (cut .gt. half .and. abs(x) .gt. cut) goto 1

  ran_gauss = x
  return
end function ran_gauss

end module collimation
