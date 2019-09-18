! ================================================================================================ !
!
!  SixTrack Collimation Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  G. Robert-Demolaize, R. Assmann, T. Weiler, C. Bracco, S. Redaelli, R. Bruce, D. Mirarchi,
!  Y.I. Levinsen, C. Tambasco, J. Molson, K.N. Sjobak, V.K. Berglyd Olsen
!
!  Created: 2004-10-27
!  Updated: 2019-09-12
!
! ================================================================================================ !
module collimation

  use parpro
  use floatPrecision
  use numerical_constants
  use coll_common, only : max_ncoll

  implicit none

  integer, parameter :: numeff     = 32
  integer, parameter :: numeffdpop = 29
  integer, parameter :: nc         = 32

  ! Logical Flags
  logical, public,  save :: do_coll          = .false.
  logical, public,  save :: coll_oldBlock    = .false.
  logical, private, save :: do_select        = .false.
  logical, private, save :: do_nominal       = .false.
  logical, private, save :: do_oneside       = .false.
  logical, private, save :: systilt_antisymm = .false.
  logical, private, save :: cern             = .false.
  logical, private, save :: do_mingap        = .false.
  logical, public,  save :: firstrun         = .true.
  logical, private, save :: cut_input        = .false. ! Not in use?

  integer, private, save :: icoll      = 0
  integer, private, save :: nloop      = 1
  integer, private, save :: ibeam      = 1
  integer, private, save :: jobnumber  = 0

  ! Distribution
  integer,          private, save :: do_thisdis   = 0
  real(kind=fPrec), public,  save :: myenom       = zero

  ! Jaw Slicing
  integer,          private, save :: n_slices     = 0
  real(kind=fPrec), private, save :: smin_slices  = zero
  real(kind=fPrec), private, save :: smax_slices  = zero
  real(kind=fPrec), private, save :: recenter1    = zero
  real(kind=fPrec), private, save :: recenter2    = zero
  real(kind=fPrec), private, save :: jaw_fit(2,6) = zero
  real(kind=fPrec), private, save :: jaw_ssf(2)   = zero

  ! Beta-beating
  real(kind=fPrec), private, save :: xbeat        = zero
  real(kind=fPrec), private, save :: xbeatphase   = zero
  real(kind=fPrec), private, save :: ybeat        = zero
  real(kind=fPrec), private, save :: ybeatphase   = zero

  real(kind=fPrec), private, save :: c_rmstilt_prim    = zero
  real(kind=fPrec), private, save :: c_rmstilt_sec     = zero
  real(kind=fPrec), private, save :: c_systilt_prim    = zero
  real(kind=fPrec), private, save :: c_systilt_sec     = zero
  real(kind=fPrec), private, save :: c_rmsoffset_prim  = zero
  real(kind=fPrec), private, save :: c_rmsoffset_sec   = zero
  real(kind=fPrec), private, save :: c_sysoffset_prim  = zero
  real(kind=fPrec), private, save :: c_sysoffset_sec   = zero
  real(kind=fPrec), private, save :: c_rmserror_gap    = zero
  integer,          private, save :: c_offsettilt_seed = 0

  ! Radial Dist
  logical,          private, save :: radial  = .false.

  ! Emittance Drift
  real(kind=fPrec), private, save :: driftsx = zero
  real(kind=fPrec), private, save :: driftsy = zero

  real(kind=fPrec), private, save :: sigsecut3 = one
  real(kind=fPrec), private, save :: sigsecut2 = one

  real(kind=fPrec), private, save :: emitnx0_dist = zero
  real(kind=fPrec), private, save :: emitny0_dist = zero
  real(kind=fPrec), private, save :: emitnx0_collgap = zero
  real(kind=fPrec), private, save :: emitny0_collgap = zero

  character(len=mNameLen),  private, save :: name_sel  = " "
  character(len=16),        private, save :: castordir = " "

  integer, save :: ie, iturn, nabs_total

  integer ieff,ieffdpop

  real(kind=fPrec), private, save :: myemitx0_dist    = zero
  real(kind=fPrec), private, save :: myemity0_dist    = zero
  real(kind=fPrec), public,  save :: myemitx0_collgap = zero
  real(kind=fPrec), public,  save :: myemity0_collgap = zero

  real(kind=fPrec), private, save :: myalphay
  real(kind=fPrec), private, save :: mybetay
  real(kind=fPrec), private, save :: myalphax
  real(kind=fPrec), private, save :: mybetax
  real(kind=fPrec), private, save :: rselect = 64.0 ! Not set anywhere, but used in if statements
  real(kind=fPrec), private, save :: myemitx

  ! M. Fiascaris for the collimation team
  ! variables for global inefficiencies studies
  ! of normalized and off-momentum halo
  ! Last modified: July 2016

  real(kind=fPrec), allocatable, save :: neff(:) !(numeff)
  real(kind=fPrec), allocatable, save :: rsig(:) !(numeff)

  integer, allocatable, save :: counteddpop(:,:) !(npart,numeffdpop)
  integer, allocatable, save :: npartdpop(:) !(numeffdpop)
  integer, allocatable, save :: counted2d(:,:,:) !(npart,numeff,numeffdpop)
  real(kind=fPrec), allocatable, save :: neffdpop(:) !(numeffdpop)
  real(kind=fPrec), allocatable, save :: dpopbins(:) !(numeffdpop)

  real(kind=fPrec) dpopmin,dpopmax,mydpop
  real(kind=fPrec), allocatable, save :: neff2d(:,:) !(numeff,numeffdpop)

  integer, allocatable, save :: nimpact(:) !(50)
  real(kind=fPrec), allocatable, save :: sumimpact(:) !(50)
  real(kind=fPrec), allocatable, save :: sqsumimpact(:) !(50)

  integer, allocatable, save :: nampl(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sum_ax(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sqsum_ax(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sum_ay(:) !(nblz)
  real(kind=fPrec), allocatable, save :: sqsum_ay(:) !(nblz)

  real(kind=fPrec), allocatable, save :: neffx(:) !(numeff)
  real(kind=fPrec), allocatable, save :: neffy(:) !(numeff)

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

  integer, save :: nsurvive, nsurvive_end, num_selhit, n_impact

  integer, allocatable, save :: cn_impact(:)  !(max_ncoll)
  integer, allocatable, save :: cn_absorbed(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: caverage(:) !(max_ncoll)
  real(kind=fPrec), allocatable, save :: csigma(:) !(max_ncoll)

  integer, allocatable, save :: counted_r(:,:) !(npart,numeff)
  integer, allocatable, save :: counted_x(:,:) !(npart,numeff)
  integer, allocatable, save :: counted_y(:,:) !(npart,numeff)

  character(len=4), save :: smpl
  character(len=80), save :: pfile

  ! Variables for finding the collimator with the smallest gap
  ! and defining, stroring the gap rms error

  character(len=mNameLen) :: coll_mingap2
  real(kind=fPrec), allocatable, save :: gap_rms_error(:) !(max_ncoll)
  real(kind=fPrec) :: nsig_err, sig_offset
  real(kind=fPrec) :: mingap, gap_h1, gap_h2, gap_h3, gap_h4
  integer :: coll_mingap_id

  real(kind=fPrec), save :: remitx_dist,remity_dist,remitx_collgap,remity_collgap

  logical, save :: firstcoll
  integer rnd_lux,rnd_k1,rnd_k2

  integer, save :: myix

  real(kind=fPrec), public  :: nspx,nspy,mux0,muy0
  real(kind=fPrec), private :: ax0,ay0,bx0,by0     ! These are set, but never used

  real(kind=fPrec), allocatable, save :: xbob(:) !(nblz)
  real(kind=fPrec), allocatable, save :: ybob(:) !(nblz)
  real(kind=fPrec), allocatable, save :: xpbob(:) !(nblz)
  real(kind=fPrec), allocatable, save :: ypbob(:) !(nblz)

  real(kind=fPrec), allocatable, save :: xineff(:) !(npart)
  real(kind=fPrec), allocatable, save :: yineff(:) !(npart)
  real(kind=fPrec), allocatable, save :: xpineff(:) !(npart)
  real(kind=fPrec), allocatable, save :: ypineff(:) !(npart)

  real(kind=fPrec), allocatable, save :: mux(:) !(nblz)
  real(kind=fPrec), allocatable, save :: muy(:) !(nblz)

  integer, save :: num_surhit
  integer, save :: numbin
  integer, save :: ibin
  integer, save :: num_selabs
  integer, save :: iturn_last_hit
  integer, save :: iturn_absorbed
  integer, save :: iturn_survive
  integer, save :: totalelem
  integer, save :: selelem
  integer, save :: unitnumber
  integer, save :: distnumber
  integer, save :: turnnumber

  real(kind=fPrec), save :: c_length    !length in m
  real(kind=fPrec), save :: c_rotation  !rotation angle vs vertical in radian
  real(kind=fPrec), save :: c_aperture  !aperture in m
  real(kind=fPrec), save :: c_offset    !offset in m
  real(kind=fPrec), save :: c_tilt(2)   !tilt in radian
  character(len=4), save :: c_material  !material

  real(kind=fPrec), allocatable, private, save :: rcx0(:)  !(npart)
  real(kind=fPrec), allocatable, private, save :: rcxp0(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcy0(:)  !(npart)
  real(kind=fPrec), allocatable, private, save :: rcyp0(:) !(npart)
  real(kind=fPrec), allocatable, private, save :: rcp0(:)  !(npart)

  real(kind=fPrec), private, save :: xj, xpj, yj, ypj, pj

  real(kind=fPrec), save :: enom_gev,betax,betay,xmax,ymax
  real(kind=fPrec), save :: nsig,calc_aperture,gammax,gammay,gammax0,gammay0,gammax1,gammay1
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

  real(kind=fPrec), save :: max_tmp, a_tmp1, a_tmp2, ldrift, mynex2, myney2, Nap1pos,Nap2pos,Nap1neg,Nap2neg
  real(kind=fPrec), save :: tiltOffsPos1,tiltOffsPos2,tiltOffsNeg1,tiltOffsNeg2
  real(kind=fPrec), save :: beamsize1, beamsize2,betax1,betax2,betay1,betay2, alphax1, alphax2,alphay1,alphay2,minAmpl

#ifdef HDF5
  ! Variables to save hdf5 dataset indices
  integer, private, save :: coll_hdf5_survival
  integer, private, save :: coll_hdf5_allImpacts
  integer, private, save :: coll_hdf5_fstImpacts
  integer, private, save :: coll_hdf5_allAbsorb
  integer, private, save :: coll_hdf5_collScatter
#endif

contains

subroutine collimation_allocate_arrays

  use mod_alloc

  implicit none

  ! Initial allocation handled by expand arrays routine
  call collimation_expand_arrays(npart,nblz)

  ! Fixed allocations follow:
  call alloc(gap_rms_error, max_ncoll, zero, "gap_rms_error") !(max_ncoll)
  call alloc(xp_pencil, max_ncoll, zero, "xp_pencil")
  call alloc(yp_pencil, max_ncoll, zero, "yp_pencil")
  call alloc(csum,      max_ncoll, zero, "csum")
  call alloc(csqsum,    max_ncoll, zero, "csqsum")

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

  call alloc(cn_impact, max_ncoll, 0, "cn_impact")  !(max_ncoll)
  call alloc(cn_absorbed, max_ncoll, 0, "cn_absorbed") !(max_ncoll)
  call alloc(caverage, max_ncoll, zero, "caverage") !(max_ncoll)
  call alloc(csigma, max_ncoll, zero, "csigma") !(max_ncoll)

end subroutine collimation_allocate_arrays

subroutine collimation_expand_arrays(npart_new, nblz_new)

  use mod_alloc
  use coll_common

  integer, intent(in) :: npart_new
  integer, intent(in) :: nblz_new

  ! Arrays that are always needed
  call alloc(part_abs_turn, npart_new, 0, "part_abs_turn") !(npart_new)

  if(.not. do_coll) return
  ! Arrays that are only needed if Collimation is enabled

  ! Allocate Common Variables
  call coll_expandArrays(npart_new, nblz_new)

  call alloc(rcx0,  npart_new, zero, "rcx0") !(npart)
  call alloc(rcxp0, npart_new, zero, "rcxp0") !(npart)
  call alloc(rcy0,  npart_new, zero, "rcy0") !(npart)
  call alloc(rcyp0, npart_new, zero, "rcyp0") !(npart)
  call alloc(rcp0,  npart_new, zero, "rcp0") !(npart)

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

  call alloc(nampl,    nblz_new, 0, "nampl") !(nblz_new)
  call alloc(sum_ax,   nblz_new, zero, "sum_ax") !(nblz_new)
  call alloc(sqsum_ax, nblz_new, zero, "sqsum_ax") !(nblz_new)
  call alloc(sum_ay,   nblz_new, zero, "sum_ay") !(nblz_new)
  call alloc(sqsum_ay, nblz_new, zero, "sqsum_ay") !(nblz_new)

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

end subroutine collimation_expand_arrays

! ================================================================================================ !
!  Collimation Init
!  This routine is called once at the start of the simulation and can be used to do any initial
!  configuration and/or file loading.
! ================================================================================================ !
subroutine collimate_init

  use crcoall
  use parpro
  use coll_common
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_settings
  use mod_time
  use string_tools
  use coll_k2
  use coll_db
  use coll_dist
  use mod_units
  use mod_ranlux
  use mod_particles
#ifdef HDF5
  use hdf5_output
#endif
#ifdef G4COLLIMATION
  use geant4
#endif

  implicit none

#ifdef HDF5
  type(h5_dataField), allocatable :: fldDist0(:)
  integer                         :: fmtDist0, setDist0
#endif
  integer :: i,j,fUnit

#ifdef HDF5
  if(h5_useForCOLL) call h5_initForCollimation
#endif

  call f_requestUnit("colltrack.out", outlun)
  call f_open(unit=outlun,file="colltrack.out",formatted=.true.,mode="w")

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

  write(lout,"(a)")       ""
  write(lout,"(a)")       str_divLine
  write(lout,"(a)")       " INITIALISING COLLIMATION"
  write(lout,"(a)")       str_divLine
  write(lout,"(a)")       ""
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

  myemitx0_dist    = remitx_dist*c1m6
  myemity0_dist    = remity_dist*c1m6
  myemitx0_collgap = remitx_collgap*c1m6
  myemity0_collgap = remity_collgap*c1m6

  myalphax = talphax(1)
  myalphay = talphay(1)
  mybetax  = tbetax(1)
  mybetay  = tbetay(1)

  if(myemitx0_dist <= zero .or. myemity0_dist <= zero .or. myemitx0_collgap <= zero .or. myemity0_collgap <= zero) then
    write(lerr,"(a)") "COLL> ERROR Emittances not defined! check collimat block!"
    write(lerr,"(a)") "COLL> ERROR Expected format of line 9 in collimat block:"
    write(lerr,"(a)") "COLL> ERROR emitnx0_dist  emitny0_dist  emitnx0_collgap  emitny0_collgap"
    write(lerr,"(a)") "COLL> ERROR All emittances should be normalized. "//&
      "first put emittance for distribtion generation, then for collimator position etc. units in [mm*mrad]."
    write(lerr,"(a)") "COLL> ERROR EXAMPLE: 2.5 2.5 3.5 3.5"
    call prror
  end if

  write(lout,"(a,i0)")    'COLL> Info: NLOOP               = ', nloop
  write(lout,"(a,i0)")    'COLL> Info: DIST_TYPES          = ', do_thisdis
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_NEX            = ', cdist_ampX
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_DEX            = ', cdist_smearX
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_NEY            = ', cdist_ampY
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_DEY            = ', cdist_smearY
  write(lout,"(a,a)")     'COLL> Info: DIST_FILE           = ', trim(cdist_fileName)
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_EN_ERROR       = ', cdist_spreadE
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_BUNCHLENGTH    = ', cdist_bunchLen
  write(lout,"(a,i0)")    'COLL> Info: RSELECT             = ', int(rselect)
  write(lout,"(a,l1)")    'COLL> Info: DO_COLL             = ', do_coll
  write(lout,"(a,l1)")    'COLL> Info: DO_NSIG             = ', cdb_doNSig
  do i=1,cdb_nFam
    write(lout,"(a,a19,a3,f13.6)") "COLL> Info: ",chr_rPad("NSIG_"//trim(chr_toUpper(cdb_famName(i))),19)," = ",cdb_famNSig(i)
  end do
  write(lout,"(a)")
  write(lout,"(a)")       'COLL> INPUT PARAMETERS FOR THE SLICING:'
  write(lout,"(a)")
  write(lout,"(a,i0)")    'COLL> Info: N_SLICES            = ', n_slices
  write(lout,"(a,e15.8)") 'COLL> Info: SMIN_SLICES         = ', smin_slices
  write(lout,"(a,e15.8)") 'COLL> Info: SMAX_SLICES         = ', smax_slices
  write(lout,"(a,e15.8)") 'COLL> Info: RECENTER1           = ', recenter1
  write(lout,"(a,e15.8)") 'COLL> Info: RECENTER2           = ', recenter2
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(1,1)        = ', jaw_fit(1,1)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(1,2)        = ', jaw_fit(1,2)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(1,3)        = ', jaw_fit(1,3)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(1,4)        = ', jaw_fit(1,4)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(1,5)        = ', jaw_fit(1,5)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(1,6)        = ', jaw_fit(1,6)
  write(lout,"(a,e15.8)") 'COLL> Info: SCALING1            = ', jaw_ssf(1)
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(2,1)        = ', jaw_fit(2,1)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(2,2)        = ', jaw_fit(2,2)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(2,3)        = ', jaw_fit(2,3)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(2,4)        = ', jaw_fit(2,4)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(2,5)        = ', jaw_fit(2,5)
  write(lout,"(a,e15.8)") 'COLL> Info: JAW_FIT(2,6)        = ', jaw_fit(2,6)
  write(lout,"(a,e15.8)") 'COLL> Info: SCALING2            = ', jaw_ssf(2)
  write(lout,"(a)")
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
  write(lout,"(a,e15.8)") 'COLL> Info: NR                  = ', cdist_ampR
  write(lout,"(a,e15.8)") 'COLL> Info: NDR                 = ', cdist_smearR
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
  write(lout,"(a,a)")     'COLL> Info: COLL_DB             = ', cdb_fileName
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
  write(lout,"(a,i0)")    'COLL> Info: NAPX                = ', napx
  write(lout,"(a,e15.8)") 'COLL> Info: Sigma_x0            = ', sqrt(mybetax*myemitx0_dist)
  write(lout,"(a,e15.8)") 'COLL> Info: Sigma_y0            = ', sqrt(mybetay*myemity0_dist)
  write(lout,"(a)")

  ! Initialize random number generator
  if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  rnd_lux = 3
  rnd_k1  = 0
  rnd_k2  = 0
  call rluxgo(rnd_lux, rnd_seed, rnd_k1, rnd_k2)
  write(outlun,*) 'INFO>  rnd_seed: ', rnd_seed

  ! Call distribution routines only if collimation block is in fort.3
  cdist_energy    = myenom
  cdist_alphaX    = myalphax
  cdist_alphaY    = myalphay
  cdist_betaX     = mybetax
  cdist_betaY     = mybetay
  cdist_emitX     = myemitx0_dist
  cdist_emitY     = myemity0_dist
  cdist_emitXColl = myemitx0_collgap
  cdist_emitYColl = myemity0_collgap
  if(radial) then
    call cdist_makeRadial
  else
    call cdist_makeDist(do_thisdis)
  end if

  ! Reset distribution for pencil beam
  if(ipencil > 0) then
    write(lout,"(a)") "COLL> WARNING Distributions reset to pencil beam!"
    write(outlun,*)   "WARN> Distributions reset to pencil beam!"
    xv1(1:napx) = zero
    yv1(1:napx) = zero
    xv2(1:napx) = zero
    yv2(1:napx) = zero
  endif

  ! Optionally write the generated particle distribution
  if(dowrite_dist) then
#ifdef HDF5
    if(h5_useForCOLL) then
      allocate(fldDist0(6))
      fldDist0(1) = h5_dataField(name="X",     type=h5_typeReal)
      fldDist0(2) = h5_dataField(name="XP",    type=h5_typeReal)
      fldDist0(3) = h5_dataField(name="Y",     type=h5_typeReal)
      fldDist0(4) = h5_dataField(name="YP",    type=h5_typeReal)
      fldDist0(5) = h5_dataField(name="SIGMA", type=h5_typeReal)
      fldDist0(6) = h5_dataField(name="E",     type=h5_typeReal)
      call h5_createFormat("collDist0", fldDist0, fmtDist0)
      call h5_createDataSet("dist0", h5_collID, fmtDist0, setDist0, napx)
      call h5_prepareWrite(setDist0, napx)
      call h5_writeData(setDist0, 1, napx, xv1(1:napx))
      call h5_writeData(setDist0, 2, napx, yv1(1:napx))
      call h5_writeData(setDist0, 3, napx, xv2(1:napx))
      call h5_writeData(setDist0, 4, napx, yv2(1:napx))
      call h5_writeData(setDist0, 5, napx, sigmv(1:napx))
      call h5_writeData(setDist0, 6, napx, ejv(1:napx))
      call h5_finaliseWrite(setDist0)
      deallocate(fldDist0)
    else
#endif
      call f_requestUnit("dist0.dat", fUnit)
      call f_open(unit=fUnit,file="dist0.dat",formatted=.true.,mode="w")
      do j=1, napx
        write(fUnit,"(6(1x,e23.15))") xv1(j), yv1(j), xv2(j), yv2(j), sigmv(j), ejv(j)
      end do
      call f_close(fUnit)
#ifdef HDF5
    end if
#endif
  end if

  ! Initialise efficiency array
  do i=1,iu
    sum_ax(i)   = zero
    sqsum_ax(i) = zero
    sum_ay(i)   = zero
    sqsum_ay(i) = zero
    nampl(i)    = zero
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

  ! Collimator Database
  call cdb_readCollDB                 ! Read the collimator DB
  call cdb_setLHCOnesided(do_oneside) ! Set LHC onesided collimators
  call cdb_writeDB_newFromOld         ! Write a copy of the db in new format, if provided in old format

  ! Then do any implementation specific initial loading
#ifdef COLLIMATE_K2
  call k2coll_init
#endif
#ifdef MERLINSCATTER
  call k2coll_merlinInit
#endif

#ifdef G4COLLIMATION

  if(n_slices /= 0) then
    write(lerr,"(a)") "COLL> ERROR Cannot use jaw fit in G4COLLIMATION version of SixTrack"
    call prror
  end if

! Open the edep file
  call f_requestUnit(fort208,unit208)
  call f_open(unit=unit208,file=fort208,formatted=.true.,mode="w")

!! This function lives in the G4Interface.cpp file in the g4collimat folder
!! Accessed by linking libg4collimat.a
!! Set the energy cut at 70% - i.e. 30% energy loss
!  g4_ecut = 0.7_fPrec
!  g4_ecut = zero

!! Select the physics engine to use
!! 0 = FTFP_BERT
!! 1 = QGSP_BERT
!  g4_physics = 0

  call g4_collimation_init(e0, rnd_seed, g4_recut, g4_aecut, g4_rcut, g4_rangecut_mm, g4_v0, trim(g4_phys_str), &
    g4_debug, g4_keep_stable, g4_edep)
#endif

  write (lout,"(a)") ""
  write (lout,"(a)") "COLL> Finished collimate initialisation"
  write (lout,"(a)") ""

  ! Adding the orbit offset at start of ring
  if(do_thisdis /= 0 .or. radial) then
    xv1(1:napx) = c1e3 * xv1(1:napx) + torbx(1)
    yv1(1:napx) = c1e3 * yv1(1:napx) + torbxp(1)
    xv2(1:napx) = c1e3 * xv2(1:napx) + torby(1)
    yv2(1:napx) = c1e3 * yv2(1:napx) + torbyp(1)
  end if

  call part_updatePartEnergy(1,.false.)

  call collimate_openFiles
  call collimate_start

end subroutine collimate_init

! ================================================================================================ !
!  Parse Input Line
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-04-16
! ================================================================================================ !
subroutine collimate_parseInputLine(inLine, iLine, iErr)

  use crcoall
  use coll_db
  use coll_dist
  use string_tools
  use coll_common
  use mod_common, only : napx

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  real(kind=fPrec) nSigIn(23), rTmp
  integer nSplit, famID
  logical spErr, fErr

  nSigIn(:) = cdb_defColGap

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "COLL> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  if(nSplit == 1 .and. iLine == 1) then
    ! This is the old block format
    coll_oldBlock = .true.
  end if
  if(coll_oldBlock) goto 10

  !  Parse new style COLL block
  ! ============================
  select case(lnSplit(1))

  case("DO_COLL")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR DO_COLL expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       DO_COLL true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), do_coll, iErr)

  case("ENERGY")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR ENERGY expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       ENERGY energy[MeV]"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), myenom, iErr)

  case("DIST_TYPE")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR DIST_TYPE expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       DIST_TYPE 0-6"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), do_thisdis, iErr)
    if(do_thisdis < 0 .or. do_thisdis > 6) then
      write(lerr,"(a,i0)") "COLL> ERROR DIST_TYPE must be between 0 and 6, got ",do_thisdis
      iErr = .true.
      return
    end if

  case("DIST_PARAM")
    if(nSplit /= 5 .and. nSplit /= 7) then
      write(lerr,"(a,i0)") "COLL> ERROR DIST_PARAM expects 4 or 6 values, got ",nSplit-1
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), cdist_ampX,   iErr)
    call chr_cast(lnSplit(3), cdist_smearX, iErr)
    call chr_cast(lnSplit(4), cdist_ampY,   iErr)
    call chr_cast(lnSplit(5), cdist_smearY, iErr)
    if(nSplit == 7) then
      call chr_cast(lnSplit(6), cdist_spreadE,  iErr)
      call chr_cast(lnSplit(7), cdist_bunchLen, iErr)
    end if

  case("DIST_FILE")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR DIST_FILE expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       DIST_FILE filename"
      iErr = .true.
      return
    end if
    cdist_fileName = trim(lnSplit(2))

  case("NSIG_FAM")
    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "COLL> ERROR NSIG_FAM expects 2 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       NSIG_FAM name nsig"
      iErr = .true.
      return
    end if
    if(len_trim(lnSplit(2)) > cdb_fNameLen) then
      write(lerr,"(2(a,i0))") "COLL> ERROR NSIG_FAM family name can be maximum ",cdb_fNameLen,&
        " characters, got ",len_trim(lnSplit(2))
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(3),rTmp,iErr)
    call cdb_addFamily(lnSplit(2),rTmp,famID,fErr)
    if(fErr) then
      write(lerr,"(a,i0)") "COLL> ERROR NSIG_FAM family '"//trim(lnSplit(2))//"' defined more than once"
      iErr = .true.
      return
    end if

  case("DO_NSIG")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR DO_NSIG expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       DO_NSIG true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), cdb_doNSig, iErr)

  case("JAW_SLICE")
    write(lerr,"(a)") "COLL> ERROR The new COLLIMATION block no longer supports the JAW_SLICE flag"
    write(lerr,"(a)") "COLL>       The feature has been moved the the collimator database"
    iErr = .true.
    return

  case("JAW_FIT1")
    if(nSplit /= 8) then
      write(lerr,"(a,i0)") "COLL> ERROR JAW_FIT1 expects 7 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       JAW_FIT1 fit1.1 fit1.2 fit1.3 fit1.4 fit1.5 fit1.6 scale"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), jaw_fit(1,1),iErr)
    call chr_cast(lnSplit(3), jaw_fit(1,2),iErr)
    call chr_cast(lnSplit(4), jaw_fit(1,3),iErr)
    call chr_cast(lnSplit(5), jaw_fit(1,4),iErr)
    call chr_cast(lnSplit(6), jaw_fit(1,5),iErr)
    call chr_cast(lnSplit(7), jaw_fit(1,6),iErr)
    call chr_cast(lnSplit(8), jaw_ssf(1),  iErr)

  case("JAW_FIT2")
    if(nSplit /= 8) then
      write(lerr,"(a,i0)") "COLL> ERROR JAW_FIT2 expects 7 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       JAW_FIT2 fit2.1 fit2.2 fit2.3 fit2.4 fit2.5 fit2.6 scale"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), jaw_fit(2,1),iErr)
    call chr_cast(lnSplit(3), jaw_fit(2,2),iErr)
    call chr_cast(lnSplit(4), jaw_fit(2,3),iErr)
    call chr_cast(lnSplit(5), jaw_fit(2,4),iErr)
    call chr_cast(lnSplit(6), jaw_fit(2,5),iErr)
    call chr_cast(lnSplit(7), jaw_fit(2,6),iErr)
    call chr_cast(lnSplit(8), jaw_ssf(2),  iErr)

  case("EMIT","EMITTANCE")
    if(nSplit /= 5) then
      write(lerr,"(a,i0)") "COLL> ERROR EMIT expects 4 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       EMIT ex_dist ey_dist ex_colgap ey_colgap"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), emitnx0_dist,   iErr)
    call chr_cast(lnSplit(3), emitny0_dist,   iErr)
    call chr_cast(lnSplit(4), emitnx0_collgap,iErr)
    call chr_cast(lnSplit(5), emitny0_collgap,iErr)

  case("DO_SELECT")
    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "COLL> ERROR DO_SELECT expects 2 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       DO_SELECT true|false name"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), do_select, iErr)
    name_sel = chr_toLower(trim(lnSplit(3)))

  case("DO_NOMINAL")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR DO_NOMINAL expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       DO_NOMINAL true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), do_nominal, iErr)

  case("SEED")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR SEED expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       SEED rnd_seed"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), rnd_seed, iErr)

  case("DO_ONESIDE")
    write(lerr,"(a)") "COLL> ERROR The new COLLIMATION block no longer supports the DO_ONESIDE flag"
    write(lerr,"(a)") "COLL>       The feature has been moved the the collimator database"
    iErr = .true.
    return

  case("WRITE_DIST")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR WRITE_DIST expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       WRITE_DIST true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dowrite_dist, iErr)

  case("WRITE_IMPACT")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR WRITE_IMPACT expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       WRITE_IMPACT true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dowrite_impact, iErr)

  case("WRITE_SECOND","WRITE_SECONDARY")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR WRITE_SECOND expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       WRITE_SECOND true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dowrite_secondary, iErr)

  case("WRITE_AMPL","WRITE_AMPLITUDE")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR WRITE_AMPL expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       WRITE_AMPL true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dowrite_amplitude, iErr)

  case("BETA_BEAT")
    if(nSplit /= 5) then
      write(lerr,"(a,i0)") "COLL> ERROR BETA_BEAT expects 4 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       BETA_BEAT xbeat xbeat_phase ybeat ybeat_phase"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), xbeat,     iErr)
    call chr_cast(lnSplit(3), xbeatphase,iErr)
    call chr_cast(lnSplit(4), ybeat,     iErr)
    call chr_cast(lnSplit(5), ybeatphase,iErr)

  case("ALIGNERR_PRIM")
    if(nSplit /= 5) then
      write(lerr,"(a,i0)") "COLL> ERROR ALIGNERR_PRIM expects 4 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       ALIGNERR_PRIM rms_tilt sys_tilt rms_offset sys_offset"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), c_rmstilt_prim,   iErr)
    call chr_cast(lnSplit(3), c_systilt_prim,   iErr)
    call chr_cast(lnSplit(4), c_rmsoffset_prim, iErr)
    call chr_cast(lnSplit(5), c_sysoffset_prim, iErr)

  case("ALIGNERR_SEC")
    if(nSplit /= 5) then
      write(lerr,"(a,i0)") "COLL> ERROR ALIGNERR_SEC expects 4 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       ALIGNERR_SEC rms_tilt sys_tilt rms_offset sys_offset"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), c_rmstilt_sec,   iErr)
    call chr_cast(lnSplit(3), c_systilt_sec,   iErr)
    call chr_cast(lnSplit(4), c_rmsoffset_sec, iErr)
    call chr_cast(lnSplit(5), c_sysoffset_sec, iErr)

  case("ALIGNERR_GAP")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR ALIGNERR_GAP expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       ALIGNERR_GAP rmserror_gap"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2),c_rmserror_gap, iErr)

  case("ALIGNERR_SEED")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR ALIGNERR_SEED expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       ALIGNERR_SEED seed"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2),c_offsettilt_seed, iErr)

  case("DO_MINGAP")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR DO_MINGAP expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       DO_MINGAP true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), do_mingap, iErr)

  case("DO_RADIAL")
    if(nSplit /= 4) then
      write(lerr,"(a,i0)") "COLL> ERROR DO_RADIAL expects 3 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       DO_RADIAL true|false amp smear"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), radial,       iErr)
    call chr_cast(lnSplit(3), cdist_ampR,   iErr)
    call chr_cast(lnSplit(4), cdist_smearR, iErr)

  case("EMIT_DRIFT")
    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "COLL> ERROR EMIT_DRIFT expects 2 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       EMIT_DRIFT driftx drifty"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), driftsx, iErr)
    call chr_cast(lnSplit(3), driftsy, iErr)

  case("SYSTILT_ANTI")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR SYSTILT_ANTI expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       SYSTILT_ANTI true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), systilt_antisymm, iErr)

  case("PENCIL")
    if(nSplit /= 6) then
      write(lerr,"(a,i0)") "COLL> ERROR PENCIL expects 5 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       PENCIL ipencil offset rmsx rmsy distr"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), ipencil,      iErr)
    call chr_cast(lnSplit(3), pencil_offset,iErr)
    call chr_cast(lnSplit(4), pencil_rmsx,  iErr)
    call chr_cast(lnSplit(5), pencil_rmsy,  iErr)
    call chr_cast(lnSplit(6), pencil_distr, iErr)

  case("COLLDB")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR COLLDB expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       COLLDB filename"
      iErr = .true.
      return
    end if
    cdb_fileName = trim(lnSplit(2))

  case("BEAM_NUM")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR BEAM_NUM expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       BEAM_NUM 1|2"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), ibeam, iErr)
    if(ibeam /= 1 .and. ibeam /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR BEAM_NUM must be 1 or 2, got ",ibeam
      iErr = .true.
      return
    end if

  case("WRITE_TRACKS")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR WRITE_TRACKS expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       WRITE_TRACKS true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dowritetracks, iErr)

  case("CERN")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR CERN expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       CERN true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), cern, iErr)

  case("SIGSECUT")
    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "COLL> ERROR SIGSECUT expects 2 values, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       SIGSECUT sigma_xy sigma_r"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), sigsecut2, iErr)
    call chr_cast(lnSplit(3), sigsecut3, iErr)

  case default
    write(lerr,"(a)") "COLL> ERROR Unknown keyword '"//trim(lnSplit(1))//"'"
    write(lerr,"(a)") "COLL>       The parser is assuming you want keyword/value block parsing as the first line had one value."
    write(lerr,"(a)") "COLL>       The two formats cannot be mixed."
    iErr = .true.
    return

  end select

  return

  !  Parse old style COLL block
  ! ============================
10 continue

  select case(iLine)

  case(1)
    if(nSplit /= 1) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 1 value on line 1, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1),do_coll,iErr)

  case(2)
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 2 values on line 2, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1),nloop,iErr)
    call chr_cast(lnSplit(2),myenom,iErr)

    if(nloop /= 1) then
      write(lerr,"(a,i0)") "COLL> ERROR Multiple samples is no longer supported. nloop must be 1, got ",nloop
      iErr = .true.
      return
    end if

  case(3)
    if(nSplit /= 8) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 8 values on line 3, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), do_thisdis,     iErr)
    call chr_cast(lnSplit(2), cdist_ampX,     iErr)
    call chr_cast(lnSplit(3), cdist_smearX,   iErr)
    call chr_cast(lnSplit(4), cdist_ampY,     iErr)
    call chr_cast(lnSplit(5), cdist_smearY,   iErr)
    cdist_fileName = lnSplit(6)
    call chr_cast(lnSplit(7), cdist_spreadE,  iErr)
    call chr_cast(lnSplit(8), cdist_bunchLen, iErr)

  case(4)
    if(nSplit /= 14) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 14 values on line 4, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), cdb_doNSig,iErr)
    call chr_cast(lnSplit(2), nSigIn(1), iErr)
    call chr_cast(lnSplit(3), nSigIn(2), iErr)
    call chr_cast(lnSplit(4), nSigIn(3), iErr)
    call chr_cast(lnSplit(5), nSigIn(4), iErr)
    call chr_cast(lnSplit(6), nSigIn(5), iErr)
    call chr_cast(lnSplit(7), nSigIn(6), iErr)
    call chr_cast(lnSplit(8), nSigIn(7), iErr)
    call chr_cast(lnSplit(9), nSigIn(8), iErr)
    call chr_cast(lnSplit(10),nSigIn(9), iErr)
    call chr_cast(lnSplit(11),nSigIn(10),iErr)
    call chr_cast(lnSplit(12),nSigIn(11),iErr)
    call chr_cast(lnSplit(13),nSigIn(12),iErr)
    call chr_cast(lnSplit(14),nSigIn(13),iErr)
    call cdb_addFamily("tcp3",   nSigIn(1), famID,fErr)
    call cdb_addFamily("tcsg3",  nSigIn(2), famID,fErr)
    call cdb_addFamily("tcsm3",  nSigIn(3), famID,fErr)
    call cdb_addFamily("tcla3",  nSigIn(4), famID,fErr)
    call cdb_addFamily("tcp7",   nSigIn(5), famID,fErr)
    call cdb_addFamily("tcsg7",  nSigIn(6), famID,fErr)
    call cdb_addFamily("tcsm7",  nSigIn(7), famID,fErr)
    call cdb_addFamily("tcla7",  nSigIn(8), famID,fErr)
    call cdb_addFamily("tclp",   nSigIn(9), famID,fErr)
    call cdb_addFamily("tcli",   nSigIn(10),famID,fErr)
    call cdb_addFamily("tcdq",   nSigIn(11),famID,fErr)
    call cdb_addFamily("tcstcdq",nSigIn(12),famID,fErr)
    call cdb_addFamily("tdi",    nSigIn(13),famID,fErr)

  case(5)
    if(nSplit < 8 .or. nSplit > 10) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 8-10 values on line 5, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), nSigIn(14),iErr)
    call chr_cast(lnSplit(2), nSigIn(15),iErr)
    call chr_cast(lnSplit(3), nSigIn(16),iErr)
    call chr_cast(lnSplit(4), nSigIn(17),iErr)
    call chr_cast(lnSplit(5), nSigIn(18),iErr)
    call chr_cast(lnSplit(6), nSigIn(19),iErr)
    call chr_cast(lnSplit(7), nSigIn(20),iErr)
    call chr_cast(lnSplit(8), nSigIn(21),iErr)
    if(nSplit > 8) then
      call chr_cast(lnSplit(9), nSigIn(22),iErr)
    end if
    if(nSplit > 9) then
      call chr_cast(lnSplit(10),nSigIn(23),iErr)
    end if
    call cdb_addFamily("tcth1",nSigIn(14),famID,fErr)
    call cdb_addFamily("tcth2",nSigIn(15),famID,fErr)
    call cdb_addFamily("tcth5",nSigIn(16),famID,fErr)
    call cdb_addFamily("tcth8",nSigIn(17),famID,fErr)
    call cdb_addFamily("tctv1",nSigIn(18),famID,fErr)
    call cdb_addFamily("tctv2",nSigIn(19),famID,fErr)
    call cdb_addFamily("tctv5",nSigIn(20),famID,fErr)
    call cdb_addFamily("tctv8",nSigIn(21),famID,fErr)
    call cdb_addFamily("tcxrp",nSigIn(22),famID,fErr)
    call cdb_addFamily("tcryo",nSigIn(23),famID,fErr)

  case(6)
    if(nSplit /= 5) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 5 values on line 6, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), n_slices,   iErr)
    call chr_cast(lnSplit(2), smin_slices,iErr)
    call chr_cast(lnSplit(3), smax_slices,iErr)
    call chr_cast(lnSplit(4), recenter1,  iErr)
    call chr_cast(lnSplit(5), recenter2,  iErr)

  case(7)
    if(nSplit /= 7) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 7 values on line 7, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), jaw_fit(1,1),iErr)
    call chr_cast(lnSplit(2), jaw_fit(1,2),iErr)
    call chr_cast(lnSplit(3), jaw_fit(1,3),iErr)
    call chr_cast(lnSplit(4), jaw_fit(1,4),iErr)
    call chr_cast(lnSplit(5), jaw_fit(1,5),iErr)
    call chr_cast(lnSplit(6), jaw_fit(1,6),iErr)
    call chr_cast(lnSplit(7), jaw_ssf(1),  iErr)

  case(8)
    if(nSplit /= 7) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 7 values on line 8, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), jaw_fit(2,1),iErr)
    call chr_cast(lnSplit(2), jaw_fit(2,2),iErr)
    call chr_cast(lnSplit(3), jaw_fit(2,3),iErr)
    call chr_cast(lnSplit(4), jaw_fit(2,4),iErr)
    call chr_cast(lnSplit(5), jaw_fit(2,5),iErr)
    call chr_cast(lnSplit(6), jaw_fit(2,6),iErr)
    call chr_cast(lnSplit(7), jaw_ssf(2),  iErr)

  case(9)
    if(nSplit /= 4) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 4 values on line 9, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), emitnx0_dist,   iErr)
    call chr_cast(lnSplit(2), emitny0_dist,   iErr)
    call chr_cast(lnSplit(3), emitnx0_collgap,iErr)
    call chr_cast(lnSplit(4), emitny0_collgap,iErr)

  case(10)
    if(nSplit /= 9) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 9 values on line 10, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), do_select,        iErr)
    call chr_cast(lnSplit(2), do_nominal,       iErr)
    call chr_cast(lnSplit(3), rnd_seed,         iErr)
    call chr_cast(lnSplit(4), dowrite_dist,     iErr)
    name_sel = chr_toLower(lnSplit(5))
    call chr_cast(lnSplit(6), do_oneside,       iErr)
    call chr_cast(lnSplit(7), dowrite_impact,   iErr)
    call chr_cast(lnSplit(8), dowrite_secondary,iErr)
    call chr_cast(lnSplit(9), dowrite_amplitude,iErr)

  case(11)
    if(nSplit /= 4) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 4 values on line 11, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), xbeat,     iErr)
    call chr_cast(lnSplit(2), xbeatphase,iErr)
    call chr_cast(lnSplit(3), ybeat,     iErr)
    call chr_cast(lnSplit(4), ybeatphase,iErr)

  case(12)
    if(nSplit /= 11) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 11 values on line 12, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), c_rmstilt_prim,   iErr)
    call chr_cast(lnSplit(2), c_rmstilt_sec,    iErr)
    call chr_cast(lnSplit(3), c_systilt_prim,   iErr)
    call chr_cast(lnSplit(4), c_systilt_sec,    iErr)
    call chr_cast(lnSplit(5), c_rmsoffset_prim, iErr)
    call chr_cast(lnSplit(6), c_rmsoffset_sec,  iErr)
    call chr_cast(lnSplit(7), c_sysoffset_prim, iErr)
    call chr_cast(lnSplit(8), c_sysoffset_sec,  iErr)
    call chr_cast(lnSplit(9), c_offsettilt_seed,iErr)
    call chr_cast(lnSplit(10),c_rmserror_gap,   iErr)
    call chr_cast(lnSplit(11),do_mingap,        iErr)

  case(13)
    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 3 values on line 13, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), radial,       iErr)
    call chr_cast(lnSplit(2), cdist_ampR,   iErr)
    call chr_cast(lnSplit(3), cdist_smearR, iErr)

  case(14)
    if(nSplit /= 4) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 4 values on line 14, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), driftsx,         iErr)
    call chr_cast(lnSplit(2), driftsy,         iErr)
    call chr_cast(lnSplit(3), cut_input,       iErr)
    call chr_cast(lnSplit(4), systilt_antisymm,iErr)

  case(15)
    if(nSplit /= 5) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 5 values on line 15, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), ipencil,      iErr)
    call chr_cast(lnSplit(2), pencil_offset,iErr)
    call chr_cast(lnSplit(3), pencil_rmsx,  iErr)
    call chr_cast(lnSplit(4), pencil_rmsy,  iErr)
    call chr_cast(lnSplit(5), pencil_distr, iErr)
#ifdef G4COLLIMATION
    if(ipencil > 0) then
      write(lerr,"(a)") "COLL> ERROR Pencil distribution not supported with geant4"
      iErr = .true.
      return
    end if
#endif

  case(16)
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 2 values on line 16, got ",nSplit
      iErr = .true.
      return
    end if
    cdb_fileName = lnSPlit(1)
    call chr_cast(lnSPlit(2), ibeam, iErr)

  case(17)
    if(nSplit /= 6) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 6 values on line 17, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), dowritetracks,iErr)
    call chr_cast(lnSplit(2), cern,         iErr)
    castordir = lnSplit(3)
    call chr_cast(lnSplit(4), jobnumber,    iErr)
    call chr_cast(lnSplit(5), sigsecut2,    iErr)
    call chr_cast(lnSplit(6), sigsecut3,    iErr)

  case default
    write(lerr,"(a,i0,a)") "COLL> ERROR Unexpected line ",iLine," encountered."
    iErr = .true.

  end select

end subroutine collimate_parseInputLine

subroutine collimate_postInput(gammar)

  use crcoall
  use coll_db

  real(kind=fPrec), intent(in) :: gammar

  call collimation_expand_arrays(npart,nblz)

  remitx_dist    = emitnx0_dist*gammar
  remity_dist    = emitny0_dist*gammar
  remitx_collgap = emitnx0_collgap*gammar
  remity_collgap = emitny0_collgap*gammar

  if(myenom == zero) then
    write(lerr,"(a)") "COLL> ERROR Beam energy cannot be zero"
    call prror
  end if

  if(cdb_fileName == " ") then
    write(lerr,"(a)") "COLL> ERROR No collimator database file specified"
    call prror
  end if

end subroutine collimate_postInput

subroutine collimate_openFiles

  use mod_units
  use string_tools
  use mod_common, only : numl
  use coll_common
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif

#ifdef HDF5
  type(h5_dataField), allocatable :: setFields(:)
  integer fmtHdf
#endif

  ! Survival File
  call f_requestUnit(coll_survivalFile,coll_survivalUnit)
  call f_open(unit=coll_survivalUnit,file=coll_survivalFile,formatted=.true.,mode="w")
  write(coll_survivalUnit,"(a7,1x,a9)") "#  turn","n_part"

  ! Collimator Gaps File
  call f_requestUnit(coll_gapsFile,coll_gapsUnit)
  call f_open(unit=coll_gapsUnit,file=coll_gapsFile,formatted=.true.,mode="w")
  write(coll_gapsUnit,"(a1,1x,a2,1x,a16,4(1x,a19),1x,a4,5(1x,a13),1x,a13)")     &
    "#","ID","name            ","angle[rad]","betax[m]","betay[m]","halfgap[m]",&
    "mat.","length[m]","sigx[m]","sigy[m]","tilt1[rad]","tilt2[rad]","nsig"

  ! Collimator Settings (Jaw Slices)
  call f_requestUnit(coll_settingsFile,coll_settingsUnit)
  call f_open(unit=coll_settingsUnit,file=coll_settingsFile,formatted=.true.,mode="w")
  write(coll_settingsUnit,"(a20,1x,a10,5(1x,a13),1x,a4)") chr_rPad("# name",20),"slice","halfgap[m]","gapoffset[m]",&
    "tiltjaw1[rad]","tiltjaw2[rad]","length[m]","mat."

  ! Positions
  call f_requestUnit(coll_positionsFile,coll_positionsUnit)
  call f_open(unit=coll_positionsUnit,file=coll_positionsFile,formatted=.true.,mode="w")
  write(coll_positionsUnit,"(a)") "# Ind           Name   Pos[m]"

  ! Twiss-Like File
  call f_requestUnit(coll_twissLikeFile,coll_twissLikeUnit)
  call f_open(unit=coll_twissLikeUnit,file=coll_twissLikeFile,formatted=.true.,mode="w")

  ! Sigma Settings File
  call f_requestUnit(coll_sigmaSetFile,coll_sigmaSetUnit)
  call f_open(unit=coll_sigmaSetUnit,file=coll_sigmaSetFile,formatted=.true.,mode="w")

  ! Tracks Files
  if(dowritetracks) then
    call f_requestUnit(coll_tracksFile,coll_tracksUnit)
    call f_open(unit=coll_tracksUnit,file=coll_tracksFile,formatted=.true.,mode="w")
    write(coll_tracksUnit,"(a)") "# name turn s x xp y yp DE/E type"

    call f_requestUnit(coll_pencilFile,coll_pencilUnit)
    call f_open(unit=coll_pencilUnit, file=coll_pencilFile,formatted=.true.,mode="w")
    write(coll_pencilUnit,"(a)") "# x xp y yp"
  end if

  if(do_select) then
    call f_requestUnit(coll_ellipseFile,coll_ellipseUnit)
    call f_open(unit=coll_ellipseUnit,file=coll_ellipseFile,formatted=.true.,mode="w")
    write(coll_ellipseUnit,"(a)") "# name x y xp yp E s turn halo nabs_type"
  end if

  if(dowrite_impact) then
    call f_requestUnit(coll_allImpactFile, coll_allImpactUnit)
    call f_requestUnit(coll_allAbsorbFile, coll_allAbsorbUnit)
    call f_requestUnit(coll_scatterFile,   coll_scatterUnit)
    call f_requestUnit(coll_fstImpactFile, coll_fstImpactUnit)
    call f_requestUnit(coll_impactFile,    coll_impactUnit)
    call f_requestUnit(coll_flukImpFile,   coll_flukImpUnit)
    call f_requestUnit(coll_flukImpAllFile,coll_flukImpAllUnit)
    call f_requestUnit(coll_jawProfileFile,coll_jawProfileUnit)

    call f_open(unit=coll_allImpactUnit, file=coll_allImpactFile, formatted=.true.,mode="w")
    call f_open(unit=coll_allAbsorbUnit, file=coll_allAbsorbFile, formatted=.true.,mode="w")
    call f_open(unit=coll_scatterUnit,   file=coll_scatterFile,   formatted=.true.,mode="w")
    call f_open(unit=coll_fstImpactUnit, file=coll_fstImpactFile, formatted=.true.,mode="w")
    call f_open(unit=coll_impactUnit,    file=coll_impactFile,    formatted=.true.,mode="w")
    call f_open(unit=coll_flukImpUnit,   file=coll_flukImpFile,   formatted=.true.,mode="w")
    call f_open(unit=coll_flukImpAllUnit,file=coll_flukImpAllFile,formatted=.true.,mode="w")
    call f_open(unit=coll_jawProfileUnit,file=coll_jawProfileFile,formatted=.true.,mode="w")

    write(coll_allImpactUnit,"(a)") "# 1=name 2=turn 3=s"
    write(coll_allAbsorbUnit,"(a)") "# 1=name 2=turn 3=s"
    write(coll_fstImpactUnit,"(a)") "# 1=name, 2=iturn, 3=icoll, 4=nabs, 5=s_imp[m], 6=s_out[m], "//&
      "7=x_in(b!)[m], 8=xp_in, 9=y_in, 10=yp_in, 11=x_out [m], 12=xp_out, 13=y_out, 14=yp_out"
    write(coll_scatterUnit,"(a)") "# 1=icoll, 2=iturn, 3=np, 4=nabs (1:Nuclear-Inelastic,2:Nuclear-Elastic,3:pp-Elastic, "//&
      "4:Single-Diffractive,5:Coulomb), 5=dp, 6=dx', 7=dy'"
    write(coll_impactUnit,"(a)") "# impact divergence"
    write(coll_flukImpUnit,"(a)") "# 1=icoll 2=c_rotation 3=s 4=x 5=xp 6=y 7=yp 8=nabs 9=np 10=turn"
    write(coll_flukImpAllUnit,"(a)") "# 1=icoll 2=c_rotation 3=s 4=x 5=xp 6=y 7=yp 8=nabs 9=np 10=turn"
    write(coll_jawProfileUnit,"(a1,1x,a6,1x,2(a7,1x),5(a17,1x),a12)") "#", "icoll", "iturn", "np", "x[m]", "xp[]", "y[m]", "yp[]",&
      "s[m]", "[1:in,2:out]"
  end if

#ifdef HDF5

  !  HDF5 Initialisation for Collimation
  ! =====================================

  if(h5_useForCOLL .eqv. .false.) return

  allocate(setFields(2))
  setFields(1) = h5_dataField(name="TURN",  type=h5_typeInt)
  setFields(2) = h5_dataField(name="NSURV", type=h5_typeInt)
  call h5_createFormat("collSurvival", setFields, fmtHdf)
  call h5_createDataSet("survival", h5_collID, fmtHdf, coll_hdf5_survival, numl)
  deallocate(setFields)

  if(dowritetracks) then
    if(h5_writeTracks2) call h5tr2_init
  end if

  if(dowrite_impact) then

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

  end if

#endif

end subroutine collimate_openFiles

!>
!! collimate_start_sample()
!! This routine is called from trauthin before each sample
!! is injected into thin 6d
!<
subroutine collimate_start

  use parpro
  use crcoall
  use coll_common
  use mod_common
  use mod_common_main
  use mod_common_track
  use coll_db
  use mod_ranlux
  use mathlib_bouncer
  use mod_units

  integer i,j,k
  real(kind=fPrec) dummy

  do i=1,napx
    do ieff=1,numeff
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
    rsig(i) = (real(i,fPrec)/two - half) + five
  end do

  dpopbins(1) = c1m4

  do i = 2, numeffdpop
    dpopbins(i) = real(i-1,fPrec)*4e-4_fPrec
  end do

  firstcoll = .true.

!GRD HERE WE NEED TO INITIALIZE SOME COLLIMATION PARAMETERS

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
  end do

#ifdef BEAMGAS
!YIL call beam gas initiation routine
  call beamGasInit(myenom)
#endif

  write(lout,"(a)") ""
  write(lout,"(a,i0)") "COLL> Number of collimators: ",cdb_nColl
  do icoll = 1, cdb_nColl
    write(lout,"(a,i5,a)") "COLL> Collimator ",icoll,": "//cdb_cName(icoll)
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

!++  Generate random tilts (Gaussian distribution plus systematic)
!++  Do this only for the first call of this routine (first sample)
!++  Keep all collimator database info and errors in memeory (COMMON
!++  block) in order to re-use exactly the same information for every
!++  sample.
  if(c_rmstilt_prim.gt.zero .or. c_rmstilt_sec.gt.zero .or. c_systilt_prim.ne.zero .or. c_systilt_sec.ne.zero) then
    do icoll = 1, cdb_nColl
      if(cdb_cName(icoll)(1:3) == "tcp") then
        c_rmstilt = c_rmstilt_prim
        c_systilt = c_systilt_prim
      else
        c_rmstilt = c_rmstilt_sec
        c_systilt = c_systilt_sec
      end if

      cdb_cTilt(1,icoll) = c_systilt+c_rmstilt*ran_gauss2(three)

      if(systilt_antisymm) then
        cdb_cTilt(2,icoll) = -one*c_systilt+c_rmstilt*ran_gauss2(three)
      else
        cdb_cTilt(2,icoll) =      c_systilt+c_rmstilt*ran_gauss2(three)
      end if

      write(outlun,*) 'INFO>  Collimator ', cdb_cName(icoll), ' jaw 1 has tilt [rad]: ', cdb_cTilt(1,icoll)
      write(outlun,*) 'INFO>  Collimator ', cdb_cName(icoll), ' jaw 2 has tilt [rad]: ', cdb_cTilt(2,icoll)
    end do
  end if

  ! In case we're using old type jaw fit, this is where we generate the parameters for the new method
  ! After this, the number of slices is also stored per collimator, and can be extracted again later
  call cdb_setMasterJawFit(n_slices, smin_slices, smax_slices, recenter1, recenter2, jaw_fit, jaw_ssf)

!++  Generate random offsets (Gaussian distribution plus systematic)
!++  Do this only for the first call of this routine (first sample)
!++  Keep all collimator database info and errors in memeory (COMMON
!++  block) in order to re-use exactly the same information for every
!++  sample and throughout a all run.
 if(c_sysoffset_prim.ne.zero .or. c_sysoffset_sec.ne.zero .or.c_rmsoffset_prim.gt.zero .or.c_rmsoffset_sec.gt.zero) then
   do icoll = 1, cdb_nColl

     if(cdb_cName(icoll)(1:3) == "tcp") then
       cdb_cOffset(icoll) = c_sysoffset_prim + c_rmsoffset_prim*ran_gauss2(three)
     else
       cdb_cOffset(icoll) = c_sysoffset_sec +  c_rmsoffset_sec*ran_gauss2(three)
     end if

     write(outlun,*) 'INFO>  offset: ', cdb_cName(icoll), cdb_cOffset(icoll)
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
    gap_rms_error(icoll) = c_rmserror_gap * ran_gauss2(three)
    write(outlun,*) 'INFO>  gap_rms_error: ', cdb_cName(icoll),gap_rms_error(icoll)
  end do

!---- creating a file with beta-functions at TCP/TCS
  mingap = 20

  do j=1,iu
! this transformation gives the right marker/name to the corresponding
! beta-dunctions or vice versa ;)
    if(ic(j) <= nblo) then
      myix = mtyp(ic(j),mel(ic(j)))
    else
      myix = ic(j)-nblo
    end if

    if(cdb_elemMap(myix) > 0) then
      nsig = cdb_cNSig(cdb_elemMap(myix))
    else
      nsig = cdb_defColGap
    end if

    do i=1,cdb_nColl
! start searching minimum gap
      if(cdb_cName(i) == bez(myix)) then
        if( cdb_cLength(i) > zero ) then
          nsig_err = nsig + gap_rms_error(i)

! jaw 1 on positive side x-axis
          gap_h1 = nsig_err - sin_mb(cdb_cTilt(1,i))*cdb_cLength(i)/2
          gap_h2 = nsig_err + sin_mb(cdb_cTilt(1,i))*cdb_cLength(i)/2

! jaw 2 on negative side of x-axis (see change of sign comapred
! to above code lines, alos have a look to setting of tilt angle)
          gap_h3 = nsig_err + sin_mb(cdb_cTilt(2,i))*cdb_cLength(i)/2
          gap_h4 = nsig_err - sin_mb(cdb_cTilt(2,i))*cdb_cLength(i)/2

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
          write(coll_twissLikeUnit,*) bez(myix),tbetax(j),tbetay(j), torbx(j),torby(j), nsig, gap_rms_error(i)
          write(coll_sigmaSetUnit,*) bez(myix), gap_h1, gap_h2, gap_h3, gap_h4, sig_offset, cdb_cOffset(i), nsig, gap_rms_error(i)

          if((gap_h1 + sig_offset) .le. mingap) then
            mingap = gap_h1 + sig_offset
            coll_mingap_id = i
            coll_mingap2 = cdb_cName(i)
          else if((gap_h2 + sig_offset) .le. mingap) then
            mingap = gap_h2 + sig_offset
            coll_mingap_id = i
            coll_mingap2 = cdb_cName(i)
          else if((gap_h3 - sig_offset) .le. mingap) then
            mingap = gap_h3 - sig_offset
            coll_mingap_id = i
            coll_mingap2 = cdb_cName(i)
          else if((gap_h4 - sig_offset) .le. mingap) then
            mingap = gap_h4 - sig_offset
            coll_mingap_id = i
            coll_mingap2 = cdb_cName(i)
          end if
        end if
      end if
    end do !do i = 1, cdb_nColl

  end do !do j=1,iu

  write(coll_twissLikeUnit,*) coll_mingap_id, coll_mingap2,  mingap
  write(coll_twissLikeUnit,*) 'INFO> IPENCIL initial ', ipencil

! if pencil beam is used and on collimator with smallest gap the
! distribution should be generated, set ipencil to coll_mingap_id
  if (ipencil.gt.0 .and. do_mingap) then
    ipencil = coll_mingap_id
  end if

  write(coll_twissLikeUnit,*) 'INFO> IPENCIL new (if do_mingap) ', ipencil
  write(coll_sigmaSetUnit,*) coll_mingap_id, coll_mingap2,  mingap

! if pencil beam is used and on collimator with smallest gap the
! distribution should be generated, set ipencil to coll_mingap_id
  write(coll_sigmaSetUnit,*) 'INFO> IPENCIL new (if do_mingap) ',ipencil
  write(coll_sigmaSetUnit,*) 'INFO> rnd_seed is (before reinit)',rnd_seed

  call f_close(coll_twissLikeUnit)
  call f_close(coll_sigmaSetUnit)

!****** re-intialize random generator with rnd_seed
!       reinit with initial value used in first call

  ! This sets the random geenrator back to the default seed rather than the one used for coll gaps.
  ! However, this doesn't actually restore the random generator to the state it would have been in without the
  ! coll gaps errors being generated as rndm5() will extract 30000 new random numbers from ranlux and discard
  ! the ones it already holds and would have used.
  ! Alternatively, we can use ranecu instead, which is capable of continuing a chain of random numbers from
  ! a given set of seeds.
  ! It is probably unnecessary to use different random seeds here in the first place.
  rnd_lux = 3
  rnd_k1  = 0
  rnd_k2  = 0
  call rluxgo(rnd_lux, rnd_seed, rnd_k2, rnd_k2)
! call recuin(rnd_seed, 0)
  dummy = rndm5(1) ! Reset rndm5 too

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

end subroutine collimate_start

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-11
!  Updated: 2019-09-11
!  Wrapper routine for main collimation call in tracking.
! ================================================================================================ !
subroutine collimate_trackThin(stracki, isColl)

  use mod_time
  use mod_common
  use mod_common_main
  use mod_common_track

  real(kind=fPrec), intent(in) :: stracki
  logical,          intent(in) :: isColl

  integer j

  if(isColl) then

    call time_startClock(time_clockCOLL)
    call collimate_start_collimator(stracki)
    call collimate_do_collimator(stracki)
    call collimate_end_collimator(stracki)
    call time_stopClock(time_clockCOLL)

  else

    gammax = (one + talphax(ie)**2)/tbetax(ie)
    gammay = (one + talphay(ie)**2)/tbetay(ie)

    do j=1,napx
      xj  = (xv1(j)-torbx(ie))/c1e3
      xpj = (yv1(j)-torbxp(ie))/c1e3
      yj  = (xv2(j)-torby(ie))/c1e3
      ypj = (yv2(j)-torbyp(ie))/c1e3
      pj  = ejv(j)/c1e3

      if(firstrun) then
        if(iturn == 1 .and. j == 1) then
          sum_ax(ie) = zero
          sum_ay(ie) = zero
        end if
      end if

      if(part_abs_pos(j) == 0 .and. part_abs_turn(j) == 0) then
        nspx         = sqrt(abs(gammax*(xj)**2 + two*talphax(ie)*xj*xpj + tbetax(ie)*xpj**2)/myemitx0_collgap)
        nspy         = sqrt(abs(gammay*(yj)**2 + two*talphay(ie)*yj*ypj + tbetay(ie)*ypj**2)/myemity0_collgap)
        sum_ax(ie)   = sum_ax(ie) + nspx
        sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
        sum_ay(ie)   = sum_ay(ie) + nspy
        sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
        nampl(ie)    = nampl(ie) + 1
      else
        nspx = zero
        nspy = zero
      end if
    end do

  end if

end subroutine collimate_trackThin

! ================================================================================================ !
!  Updated: 2019-09-12
!  Called when we enter a collimator.
!  The routine doesn't check that the collimator exists! This should be checked before calling it.
! ================================================================================================ !
subroutine collimate_start_collimator(stracki)

  use coll_db
  use mod_common
  use mod_common_main
  use mod_common_track

  real(kind=fPrec), intent(in) :: stracki

  integer j

  icoll    = cdb_elemMap(myix)
  nsig     = cdb_cNSig(icoll)
  c_length = zero

  ! SR, 23-11-2005: To avoid binary entries in 'amplitude.dat'
  if(firstrun) then

    gammax = (one + talphax(ie)**2)/tbetax(ie)
    gammay = (one + talphay(ie)**2)/tbetay(ie)

    do j=1,napx
      xj  = (xv1(j)-torbx(ie))/c1e3
      xpj = (yv1(j)-torbxp(ie))/c1e3
      yj  = (xv2(j)-torby(ie))/c1e3
      ypj = (yv2(j)-torbyp(ie))/c1e3
      pj  = ejv(j)/c1e3

      if(iturn == 1 .and. j == 1) then
        sum_ax(ie) = zero
        sum_ay(ie) = zero
      end if

      ! DRIFT PART
      if(stracki == zero) then
        if(iexact) then
          zpj = sqrt(one-xpj**2-ypj**2)
          xj  = xj + half*c_length*(xpj/zpj)
          yj  = yj + half*c_length*(ypj/zpj)
        else
          xj  = xj + half*c_length*xpj
          yj  = yj + half*c_length*ypj
        end if
      end if

      if(part_abs_pos(j) == 0 .and. part_abs_turn(j) == 0) then
        nspx = sqrt(abs(gammax*xj**2 + two*talphax(ie)*xj*xpj +tbetax(ie)*xpj**2)/myemitx0_collgap)
        nspy = sqrt(abs(gammay*yj**2 + two*talphay(ie)*yj*ypj +tbetay(ie)*ypj**2)/myemity0_collgap)
        sum_ax(ie)   = sum_ax(ie) + nspx
        sqsum_ax(ie) = sqsum_ax(ie) + nspx**2
        sum_ay(ie)   = sum_ay(ie) + nspy
        sqsum_ay(ie) = sqsum_ay(ie) + nspy**2
        nampl(ie)    = nampl(ie) + 1
      else
        nspx = zero
        nspy = zero
      end if
    end do
  end if

end subroutine collimate_start_collimator

!>
!! collimate_do_collimator()
!! This routine is calls the actual scattering functions
!<
subroutine collimate_do_collimator(stracki)

  use crcoall
  use parpro
  use coll_common
  use mod_common
  use mod_commons
  use mod_common_da
  use mod_common_main
  use mod_common_track
  use numerical_constants, only : c5m4
  use coll_db
  use coll_k2
  use coll_jawfit
  use coll_dist
  use mod_units
  use mathlib_bouncer
  use mod_alloc
  use string_tools
#ifdef ROOT
  use root_output
#endif
#ifdef G4COLLIMATION
  use geant4
#endif

  implicit none

  real(kind=fPrec), intent(in) :: stracki

  integer j, iSlice, nSlices
  logical onesided, linside(napx)
  real(kind=fPrec) jawLength, jawAperture, jawOffset, jawTilt(2)
  real(kind=fPrec) x_Dump,xpDump,y_Dump,ypDump,s_Dump

#ifdef G4COLLIMATION
  integer :: g4_lostc
  integer :: g4_npart
  integer :: part_hit_flag = 0
  integer :: part_abs_flag = 0
  real(kind=fPrec) x_tmp,y_tmp,xp_tmp,yp_tmp

  ! ien0,ien1: ion energy entering/leaving the collimator
  ! energy in MeV
  real(kind=fPrec)    :: ien0, ien1
  integer(kind=int16) :: nnuc0,nnuc1
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
  if(chr_toLower(cdb_cName(icoll)) == name_sel .and. do_select) then
    do j=1,napx
      write(coll_ellipseUnit,"(1x,i8,6(1x,e15.7),3(1x,i4,1x,i4))") partID(j),xv1(j),xv2(j),yv1(j),yv2(j), &
        ejv(j),sigmv(j),iturn,secondary(j)+tertiary(j)+other(j)+scatterhit(j),nabs_type(j)
    end do
  end if

  if(iturn == 1 .and. firstrun) then
    write(outlun,*) ' '
    write(outlun,*)   'Collimator information: '
    write(outlun,*) ' '
    write(outlun,*) 'Name:                ', cdb_cName(icoll)
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
  c_tilt(1)  = cdb_cTilt(1,icoll)
  c_tilt(2)  = cdb_cTilt(2,icoll)

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
! important: Sixtrack calculates in "mm" and k2coll_collimate in "m"
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

  ! Extract number of jaw fit slices for this collimator
  if(cdb_cSliced(icoll) > 0) then ! Collimator is sliced
    nSlices = jaw_getSliceCount(cdb_cSliced(icoll))
    if(nSlices < 0) then
      write(lerr,"(a)")    "COLL> ERROR Invalid entry in jaw fit database for collimator '"//trim(cdb_cName(icoll))//"'"
      write(lerr,"(a,i0)") "COLL>       Value returned for number of slices is ",nSlices
      call prror
    end if
  else
    nSlices = 1
  end if

  ! Further output
  if(firstrun) then
    if(iturn == 1) then
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

      write(coll_gapsUnit,"(i4,1x,a16,4(1x,e19.10),1x,a4,5(1x,e13.5),1x,f13.6)") &
        icoll,cdb_cName(icoll)(1:16),cdb_cRotation(icoll),tbetax(ie),tbetay(ie),calc_aperture, &
        cdb_cMaterial(icoll),cdb_cLength(icoll),sqrt(tbetax(ie)*myemitx0_collgap), &
        sqrt(tbetay(ie)*myemity0_collgap),cdb_cTilt(1,icoll),cdb_cTilt(2,icoll),nsig

      ! Write to coll settings file if we have 0 or 1 slices
      if(nSlices <= 1) then
        write(coll_settingsUnit,"(a20,1x,i10,5(1x,1pe13.6),1x,a)") cdb_cName(icoll)(1:20), nSlices, calc_aperture, &
          cdb_cOffset(icoll), cdb_cTilt(1,icoll), cdb_cTilt(2,icoll), cdb_cLength(icoll), cdb_cMaterial(icoll)
      end if
    end if
  end if

  c_aperture = two*calc_aperture
!          IF(IPENCIL.GT.zero) THEN
!          C_APERTURE = 2.*pencil_aperture

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
    if((cdist_ampX.gt.0).and.(cdist_ampY.eq.zero)) then  ! horizontal halo
      beamsize1 = sqrt(betax1 * myemitx0_collgap)
      beamsize2 = sqrt(betax2 * myemitx0_collgap)
    else if((cdist_ampX.eq.0).and.(cdist_ampY.gt.zero)) then   ! vertical halo
      beamsize1 = sqrt(betay1 * myemity0_collgap)
      beamsize2 = sqrt(betay2 * myemity0_collgap)
    else
      write(lerr,"(a)") "COLL> ERROR Attempting to use a halo not purely in the horizontal "//&
        "or vertical plane with pencil_dist=3 - abort."
      call prror
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
    if((cdist_ampX.gt.0).and.(cdist_ampY.eq.zero)) then ! horizontal halo
      mynex2 = minAmpl
    else if((cdist_ampX.eq.0).and.(cdist_ampY.gt.zero)) then ! vertical halo
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
    call cdist_makeDist_coll(myalphax,myalphay,mybetax,mybetay,mynex2,myney2)

    do j = 1, napx
      xv1(j)  = c1e3*xv1(j) + torbx(ie)
      yv1(j)  = c1e3*yv1(j) + torbxp(ie)
      xv2(j)  = c1e3*xv2(j) + torby(ie)
      yv2(j)  = c1e3*yv2(j) + torbyp(ie)

!      as main routine will track particles back half a collimator length (to start of jaw),
!      track them now forward (if generated at face) or backward (if generated at end)
!      1/2 collimator length to center of collimator (ldrift pos or neg)
       xv1(j)  = xv1(j) - ldrift*yv1(j)
       xv2(j)  = xv2(j) - ldrift*yv2(j)
    end do
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end RB addition

!++  Copy particle data to 1-dim array and go back to meters
  do j=1,napx
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

    ! For zero length element track back half collimator length
    ! DRIFT PART
    if(stracki == 0) then
      if(iexact) then
        zpj=sqrt(one-rcxp(j)**2-rcyp(j)**2)
        rcx(j) = rcx(j) - (half*c_length)*(rcxp(j)/zpj)
        rcy(j) = rcy(j) - (half*c_length)*(rcyp(j)/zpj)
      else
        rcx(j) = rcx(j) - (half*c_length)*rcxp(j)
        rcy(j) = rcy(j) - (half*c_length)*rcyp(j)
      end if
    else
      write(lerr,"(a,f13.6)") "COLL> ERROR Non-zero length collimator '"//trim(cdb_cName(icoll))//"' with length = ",stracki
      call prror
    end if
  end do

!++  Do the collimation tracking
  enom_gev = myenom*c1m3

  ! Allow treatment of collimators as one-sided
  if(cdb_cSides(icoll) == 1) then
    onesided = .true.
  else if(cdb_cSides(icoll) == 2) then
    onesided = .true.
  else
    onesided = .false.
  end if

  linside(:) = .false.

  if(cdb_cSliced(icoll) > 0) then ! Treatment of sliced collimators
    ! Now, loop over the number of slices and call k2coll_collimate each time.
    ! For each slice, the corresponding offset and angle are to be used.
    do iSlice=1,nSlices
      jawAperture = c_aperture
      jawOffset   = c_offset
      jawTilt     = c_tilt
      call jaw_getFitSliceValues(cdb_cSliced(icoll), iSlice, jawLength, jawAperture, jawOffset, jawTilt)
      if(firstrun) then
        write(coll_settingsUnit,"(a20,1x,i10,5(1x,1pe13.6),1x,a)") cdb_cName(icoll)(1:20), iSlice,  &
          jawAperture/two, jawOffset, jawTilt(1), jawTilt(2), jawLength, cdb_cMaterial(icoll)
      end if
      call k2coll_collimate(icoll, iturn, ie, jawLength, c_rotation, jawAperture,       &
        jawOffset, jawTilt, rcx, rcxp, rcy, rcyp, rcp, rcs, enom_gev, part_hit_pos,           &
        part_hit_turn, part_abs_pos, part_abs_turn, part_impact, part_indiv, part_linteract,        &
        onesided, secondary, iSlice, nabs_type, linside)
    end do

  else ! Treatment of non-sliced collimators

#ifndef G4COLLIMATION

    call k2coll_collimate(icoll, iturn, ie, c_length, c_rotation, c_aperture, c_offset, &
      c_tilt, rcx, rcxp, rcy, rcyp, rcp, rcs, enom_gev, part_hit_pos,part_hit_turn,           &
      part_abs_pos, part_abs_turn, part_impact, part_indiv, part_linteract, onesided, secondary, 1, &
      nabs_type, linside)

#else

    !! Add the geant4 geometry
    if(firstrun .and. iturn == 1) then
      call g4_add_collimator(cdb_cName(icoll), c_material, c_length, c_aperture, c_rotation, torbx(ie), torby(ie))
    end if

!! Here we do the real collimation
!! First set the correct collimator
    call g4_set_collimator(cdb_cName(icoll))
    flush(lout)

!! Loop over all our particles
    g4_lostc = 0
    nnuc0 = 0
    ien0  = zero
    nnuc1 = 0
    ien1  = zero

    if(g4_debug .eqv. .true.) then
      write(lout,"(2a)") 'COLLIMATOR:', cdb_cName(icoll)
      write(lout,"(12a)") chr_lpad('id',33), chr_lpad('pdgid',12), chr_lpad('mass',25), chr_lpad('x',25), chr_lpad('y',25), &
                          chr_lpad('xp',25), chr_lpad('yp',25), chr_lpad('p',25), chr_lpad('spin_x',25), chr_lpad('spin_y',25),&
                          chr_lpad('spin_z',25)
      flush(lout)
    end if

    do j = 1, napx
!!!!          if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
!! Rotate particles in the frame of the collimator
!! There is more precision if we do it here rather
!! than in the g4 geometry

      if(g4_debug .eqv. .true.) then
        write(lout,"(a,2(1X,I11),10(1X,E24.16))") 'g4 sending particle: ', j, pdgid(j), nucm(j), rcx(j), rcy(j), rcxp(j), &
          rcyp(j), rcp(j), spin_x(j), spin_y(j), spin_z(j), sigmv(j)
      end if

      x_tmp = rcx(j)
      y_tmp = rcy(j)
      xp_tmp = rcxp(j)
      yp_tmp = rcyp(j)
      rcx(j) =  x_tmp *cos_mb(c_rotation) + sin_mb(c_rotation)*y_tmp
      rcy(j) =  y_tmp *cos_mb(c_rotation) - sin_mb(c_rotation)*x_tmp
      rcxp(j) = xp_tmp*cos_mb(c_rotation) + sin_mb(c_rotation)*yp_tmp
      rcyp(j) = yp_tmp*cos_mb(c_rotation) - sin_mb(c_rotation)*xp_tmp

!! Add all particles
      call g4_add_particle(rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j), pdgid(j), nzz(j), naa(j), nqq(j), nucm(j), &
        sigmv(j), spin_x(j), spin_y(j), spin_z(j))

! Log input energy + nucleons as per the FLUKA coupling
      nnuc0   = nnuc0 + naa(j)
      ien0    = ien0 + rcp(j) * c1e3
    end do

!! Call the geant4 collimation function
!            call g4_collimate(rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j))
    call g4_collimate()

!! Get the particle number back
    call g4_get_particle_count(g4_npart)

!! resize arrays
    call expand_arrays(nele, g4_npart, nblz, nblo, nbb)

!! Reset napx to the correct value
    napx = g4_npart

    if(g4_debug .eqv. .true.) then
      write(lout,"(12a)") chr_lpad('id',33), chr_lpad('pdgid',12), chr_lpad('mass',25), chr_lpad('x',25), chr_lpad('y',25), &
                          chr_lpad('xp',25), chr_lpad('yp',25), chr_lpad('p',25), chr_lpad('spin_x',25), chr_lpad('spin_y',25),&
                          chr_lpad('spin_z',25)
      flush(lout)
    end if

    do j = 1, napx
!! Get the particle back + information
!! Remember C arrays start at 0, fortran at 1 here.
      call g4_collimate_return(j-1, rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j), pdgid(j), nucm(j), nzz(j), naa(j), nqq(j), &
        sigmv(j), part_hit_flag, part_abs_flag, part_impact(j), part_indiv(j), part_linteract(j), spin_x(j), spin_y(j), spin_z(j))

      partID(j) = j
      pstop (j) = .false.

!! Rotate back into the accelerator frame
      x_tmp   = rcx(j)
      y_tmp   = rcy(j)
      xp_tmp  = rcxp(j)
      yp_tmp  = rcyp(j)
      rcx(j)  = x_tmp *cos_mb(-one*c_rotation) + sin_mb(-one*c_rotation)*y_tmp
      rcy(j)  = y_tmp *cos_mb(-one*c_rotation) - sin_mb(-one*c_rotation)*x_tmp
      rcxp(j) = xp_tmp*cos_mb(-one*c_rotation) + sin_mb(-one*c_rotation)*yp_tmp
      rcyp(j) = yp_tmp*cos_mb(-one*c_rotation) - sin_mb(-one*c_rotation)*xp_tmp

! This needs fixing - FIXME
!            sigmv(j) = zero
!            sigmv(j) = s - (g4_v0*g4_time)
      part_impact(j) = 0
      part_indiv(j) = 0
      part_linteract(j) = 0

! Log output energy + nucleons as per the FLUKA coupling
      nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
      ien1        = ien1  + rcp(j) * c1e3                   ! outcoming energy

! Fix hits
! if(part_hit_pos(j) .eq.ie .and. part_hit_turn(j).eq.iturn)
      part_hit_pos(j)  = ie
      part_hit_turn(j) = iturn
      part_abs_pos(j) = 0
      part_abs_turn(j) = 0

!!           If a particle hit
!            if(part_hit_flag.ne.0) then
!              part_hit_pos(j) = ie
!              part_hit_turn(j) = iturn
!            end if
!
!!           If a particle died (the checking if it is already dead is at the start of the loop)
!!           Geant just has a general inelastic process that single diffraction is part of
!!           Therefore we can not know if this interaction was SD or some other inelastic type
!            if(part_abs_flag.ne.0) then
!              if(dowrite_impact) then
!!! FLUKA_impacts.dat
!                write(coll_flukImpUnit,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))') &
! &                    icoll,c_rotation,zero,zero,zero,zero,zero,part_abs_flag,flukaname(j),iturn
!              end if
!
!              part_abs_pos(j)  = ie
!              part_abs_turn(j) = iturn
!              rcx(j) = 99.99e-3_fPrec
!              rcy(j) = 99.99e-3_fPrec
!              g4_lostc = g4_lostc + 1
!            end if

      if(g4_debug .eqv. .true.) then
        write(lout,"(a,2(1X,I11),10(1X,E24.16))") 'g4 return particle:  ', j, pdgid(j), nucm(j), rcx(j), rcy(j), rcxp(j), &
              rcyp(j), rcp(j), spin_x(j), spin_y(j), spin_z(j), sigmv(j)
      end if

      flush(lout)
!!!!          end if !part_abs_pos(j) .ne. 0 .and. part_abs_turn(j) .ne. 0
    end do   !do j = 1, napx

    call g4_collimation_clear()

    if((ien0-ien1) > one) then
#ifdef ROOT
      if(root_flag .and. root_Collimation == 1) then
        call root_EnergyDeposition(icoll, nnuc0-nnuc1,c1m3*(ien0-ien1))
      end if
#endif
      write(unit208,"(2(i5,1x),e24.16)") icoll, (nnuc0-nnuc1), c1m3*(ien0-ien1)
      flush(unit208)
    end if

#endif
  end if

  if(dowrite_impact) then
    ! update writeout of jaw profiles
    do j=1,napx
      if(linside(j) .and. sqrt(rcx(j)**2 + rcy(j)**2) < 99.0e-3_fPrec) then
        x_Dump = rcx (j)*cos_mb(c_rotation)+sin_mb(c_rotation)*rcy (j)
        xpDump = rcxp(j)*cos_mb(c_rotation)+sin_mb(c_rotation)*rcyp(j)
        y_Dump = rcy (j)*cos_mb(c_rotation)-sin_mb(c_rotation)*rcx (j)
        ypDump = rcyp(j)*cos_mb(c_rotation)-sin_mb(c_rotation)*rcxp(j)
        s_Dump = c_length
        write(coll_jawProfileUnit,"(3(1x,i7),5(1x,e17.9),1x,i1)") &
          icoll,iturn,partID(j),x_Dump,xpDump,y_Dump,ypDump,s_Dump,2
      end if
    end do
  end if

end subroutine collimate_do_collimator

!>
!! collimate_end_collimator()
!! This routine is called at the exit of a collimator
!<
subroutine collimate_end_collimator(stracki)

  use crcoall
  use parpro
  use coll_common
  use mod_common
  use mod_commons
  use mod_common_da
  use mod_common_main
  use mod_common_track
  use numerical_constants, only : c5m4
  use coll_db
  use mod_units
  use string_tools
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif
#ifdef G4COLLIMATION
  use geant4
#endif

  implicit none

  integer :: j

#ifdef HDF5
  ! For tracks2
  integer hdfturn,hdfpid,hdftyp
  real(kind=fPrec) hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdfs
#endif

  real(kind=fPrec), intent(in) :: stracki

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

#ifdef G4COLLIMATION
  do j = 1, napx
      if (stracki.eq.0.) then
        if(iexact .eqv. .false.) then
!          write(lout,*) 'iexact = 0', rcxp(j), rcyp(j)
          rcx(j)  = rcx(j) - half*c_length*rcxp(j)
          rcy(j)  = rcy(j) - half*c_length*rcyp(j)
        else
!          write(lout,*) 'iexact = 1', rcxp(j), rcyp(j)
          zpj=sqrt(one-rcxp(j)**2-rcyp(j)**2)
          rcx(j) = rcx(j) - half*c_length*(rcxp(j)/zpj)
          rcy(j) = rcy(j) - half*c_length*(rcyp(j)/zpj)
        end if
      end if

!++  Now copy data back to original verctor
      xv1(j) = rcx(j)  * c1e3 + torbx(ie)
      yv1(j) = rcxp(j) * c1e3 + torbxp(ie)
      xv2(j) = rcy(j)  * c1e3 + torby(ie)
      yv2(j) = rcyp(j) * c1e3 + torbyp(ie)
      ejv(j) = rcp(j)  * c1e3


! Update mtc and other arrays.
!      ejv(j)    = rcp(j)  * c1e3
!!!  write(lout,*) 'ejfv', ejv(j), nucm(j)
!      ejfv  (j) = sqrt((ejv(j)-nucm(j))*(ejv(j)+nucm(j)))   ! hisix: ion mass
      ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
!!!  write(lout,*) 'ejfv ok', ejfv(j)
      rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))                 ! hisix: remains unchanged
      dpsv  (j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f         ! hisix: new delta
      oidpsv(j) = one/(one+dpsv(j))
      dpsv1 (j) = (dpsv(j)*c1e3)*oidpsv(j)
      mtc     (j) = (nqq(j)*nucm0)/(qq0*nucm(j))            ! hisix: mass to charge
      moidpsv (j) = mtc(j)*oidpsv(j)                        ! hisix
      omoidpsv(j) = c1e3*((one-mtc(j))*oidpsv(j))           ! hisix

!++  Energy update, as recommended by Frank
!      ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
!      rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
!      dpsv(j)=(ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
!      oidpsv(j)=one/(one+dpsv(j))
!      moidpsv(j)=mtc(j)/(one+dpsv(j))
!      omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
!      dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)
      yv1(j)   = ejf0v(j)/ejfv(j)*yv1(j)
      yv2(j)   = ejf0v(j)/ejfv(j)*yv2(j)
!!!  write(lout,*) 'Coordinate loop end'
end do
#endif

#ifndef G4COLLIMATION
!++  Copy particle data back and do path length stuff; check for absorption
!++  Add orbit offset back.
  do j = 1, napx

!APRIL2005 IN ORDER TO GET RID OF NUMERICAL ERRORS, JUST DO THE TREATMENT FOR
!APRIL2005 IMPACTING PARTICLES...
    if(part_hit_pos(j) .eq.ie .and. part_hit_turn(j).eq.iturn) then
!++  For zero length element track back half collimator length
! DRIFT PART
      if (stracki.eq.0.) then
        if(iexact) then
          zpj=sqrt(one-rcxp(j)**2-rcyp(j)**2)
          rcx(j) = rcx(j) - half*c_length*(rcxp(j)/zpj)
          rcy(j) = rcy(j) - half*c_length*(rcyp(j)/zpj)
        else
          rcx(j)  = rcx(j) - half*c_length*rcxp(j)
          rcy(j)  = rcy(j) - half*c_length*rcyp(j)
        end if
      end if

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
      mtc(j) = (nqq(j)*nucm0)/(qq0*nucm(j))            ! hisix: mass to charge
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
        ejv(j) = myenom
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
          call h5_writeData(coll_hdf5_allImpacts, 1, 1, partID(j))
          call h5_writeData(coll_hdf5_allImpacts, 2, 1, iturn)
          call h5_writeData(coll_hdf5_allImpacts, 3, 1, dcum(ie))
          call h5_finaliseWrite(coll_hdf5_allImpacts)
        else
#endif
          write(coll_allImpactUnit,"(i8,1x,i4,1x,f8.2)") partID(j),iturn,dcum(ie)
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
            call h5_writeData(coll_hdf5_allAbsorb, 1, 1, partID(j))
            call h5_writeData(coll_hdf5_allAbsorb, 2, 1, iturn)
            call h5_writeData(coll_hdf5_allAbsorb, 3, 1, dcum(ie))
            call h5_finaliseWrite(coll_hdf5_allAbsorb)
          else
#endif
            write(coll_allAbsorbUnit,"(i8,1x,i4,1x,f8.2)") partID(j),iturn,dcum(ie)
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
        if(cdb_cName(icoll)(1:3) == "tcp") then
          secondary(j) = 1
        else if(cdb_cName(icoll)(1:3) == "tcs") then
          tertiary(j)  = 2
        else if((cdb_cName(icoll)(1:3) == "tcl") .or. (cdb_cName(icoll)(1:3) == "tct") .or. &
                (cdb_cName(icoll)(1:3) == "tcd") .or. (cdb_cName(icoll)(1:3) == "tdi")) then
          other(j)     = 4
        end if
      else
        write(lerr,"(a)")          "COLL> ERROR Particle cannot be both absorbed and not absorbed"
        write(lerr,"(a,2(1x,i0))") "COLL>      ",part_abs_pos (j),part_abs_turn(j)
        call prror
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
             ((((xv1(j)*c1m3)**2 / (tbetax(ie)*myemitx0_collgap)) .ge. sigsecut2) .or. &
             (((xv2(j)*c1m3)**2  / (tbetay(ie)*myemity0_collgap)) .ge. sigsecut2) .or. &
             (((xv1(j)*c1m3)**2  / (tbetax(ie)*myemitx0_collgap)) + &
             ((xv2(j)*c1m3)**2   / (tbetay(ie)*myemity0_collgap)) .ge. sigsecut3)) ) then

            xj  = (xv1(j)-torbx(ie))  /c1e3
            xpj = (yv1(j)-torbxp(ie)) /c1e3
            yj  = (xv2(j)-torby(ie))  /c1e3
            ypj = (yv2(j)-torbyp(ie)) /c1e3

#ifdef HDF5
            if(h5_writeTracks2) then
              ! We write trajectories before and after element in this case.
              hdfpid  = partID(j)
              hdfturn = iturn
              hdfs    = dcum(ie)-half*c_length
              hdfx    = (rcx0(j)*c1e3+torbx(ie)) - half*c_length*(rcxp0(j)*c1e3+torbxp(ie))
              hdfxp   = rcxp0(j)*c1e3+torbxp(ie)
              hdfy    = (rcy0(j)*c1e3+torby(ie)) - half*c_length*(rcyp0(j)*c1e3+torbyp(ie))
              hdfyp   = rcyp0(j)*c1e3+torbyp(ie)
              hdfdee  = (ejv(j)-myenom)/myenom
              hdftyp  = secondary(j)+tertiary(j)+other(j)+scatterhit(j)
              call h5tr2_writeLine(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdftyp)

              hdfs  = dcum(ie)+half*c_length
              hdfx  = xv1(j) + half*c_length*yv1(j)
              hdfxp = yv1(j)
              hdfy  = xv2(j) + half*c_length*yv2(j)
              hdfyp = yv2(j)
              call h5tr2_writeLine(hdfpid,hdfturn,hdfs,hdfx,hdfxp,hdfy,hdfyp,hdfdee,hdftyp)
            else
#endif
              write(coll_tracksUnit,"(1x,i8,1x,i4,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)") &
                partID(j),iturn,dcum(ie)-half*c_length,                            &
                (rcx0(j)*c1e3+torbx(ie))-half*c_length*(rcxp0(j)*c1e3+torbxp(ie)), &
                rcxp0(j)*c1e3+torbxp(ie),                                          &
                (rcy0(j)*c1e3+torby(ie))-half*c_length*(rcyp0(j)*c1e3+torbyp(ie)), &
                rcyp0(j)*c1e3+torbyp(ie),                                          &
                (ejv(j)-myenom)/myenom,secondary(j)+tertiary(j)+other(j)+scatterhit(j)

              write(coll_tracksUnit,"(1x,i8,1x,i4,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)") &
                partID(j),iturn,dcum(ie)+half*c_length,                         &
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
  if(chr_toLower(cdb_cName(icoll)) == name_sel .and. iturn == 1) then
    num_selhit = 0
    num_surhit = 0
    num_selabs = 0

    do j=1,napx
      if(part_hit_pos (j) == ie .and. part_hit_turn(j) == iturn) then
        num_selhit = num_selhit+1
        if(part_abs_pos(j)  == 0 .and. part_abs_turn(j) == 0) then
          num_surhit = num_surhit+1
        else
          num_selabs = num_selabs + 1
        end if
!++  If we want to select only partciles interacting at the specified
!++  collimator then remove all other particles and reset the number
!++  of the absorbed particles to the selected collimator.
      else if(do_select .and. firstrun) then
        part_select(j) = 0
        n_tot_absorbed = num_selabs
      end if
    end do

!++  Calculate average impact parameter and save distribution into file
!++  only for selected collimator
    n_impact = 0
    sum      = zero
    sqsum    = zero

    do j=1,napx
      if(part_hit_pos(j) == ie .and. part_hit_turn(j) == iturn) then
        if(part_impact(j) < -half) then
          write(lerr,"(a,i0)") "COLL> ERROR Found invalid impact parameter ", part_impact(j)
          write(outlun,*) 'ERR>  Invalid impact parameter!', part_impact(j)
          call prror
        end if

        n_impact = n_impact + 1
        sum = sum + part_impact(j)
        sqsum = sqsum + part_impact(j)**2
        if(part_hit_pos (j).ne.0 .and. part_hit_turn(j).ne.0 .and.dowrite_impact ) then
          write(coll_impactUnit,*) part_impact(j), part_indiv(j)
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

    if(dowrite_impact) then
      call f_close(coll_impactUnit)
    end if

!++  End of    S E L E C T E D   collimator
  end if
#endif


end subroutine collimate_end_collimator

! ================================================================================================ !
!  Collimate Exit
! ================================================================================================ !
subroutine collimate_exit

  use parpro
  use coll_common
  use mod_common
  use mod_commons
  use mod_common_da
  use mod_common_main
  use mod_common_track
  use crcoall
  use coll_db
  use mod_units
  use string_tools
#ifdef ROOT
  use root_output
#endif
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif
#ifdef G4COLLIMATION
  use geant4
#endif

  implicit none

  ! integer, intent(in) :: j

#ifdef HDF5
  type(h5_dataField), allocatable :: fldHdf(:)
  integer fmtHdf, setHdf
#endif
  integer i,k,ix

!++  Save particle offsets to a file
  call f_close(coll_survivalUnit)

  if(dowrite_impact) then
    call f_close(coll_impactUnit)
  end if

  if(dowritetracks) then
    call f_close(coll_tracksUnit)
#ifdef HDF5
    if(h5_writeTracks2) call h5tr2_finalise
#endif
  end if

!------------------------------------------------------------------------
!++  Write the number of absorbed particles
  write(outlun,*) 'INFO>  Number of impacts             : ', n_tot_absorbed+nsurvive_end
  write(outlun,*) 'INFO>  Number of impacts at selected : ', num_selhit
  write(outlun,*) 'INFO>  Number of surviving particles : ', nsurvive_end
  write(outlun,*) 'INFO>  Number of absorbed particles  : ', n_tot_absorbed
  write(outlun,*)

  if(n_tot_absorbed.ne.0) then
    write(outlun,*) ' INFO>  Eff_r @  8 sigma    [e-4] : ', (neff(5)/real(n_tot_absorbed,fPrec))/c1m4
    write(outlun,*) ' INFO>  Eff_r @ 10 sigma    [e-4] : ', (neff(9)/real(n_tot_absorbed,fPrec))/c1m4
    write(outlun,*) ' INFO>  Eff_r @ 10-20 sigma [e-4] : ', ((neff(9)-neff(19))/(real(n_tot_absorbed,fPrec)))/c1m4
    write(outlun,*)
    write(outlun,*) neff(5)/real(n_tot_absorbed,fPrec), neff(9)/real(n_tot_absorbed,fPrec), &
      (neff(9)-neff(19))/(real(n_tot_absorbed,fPrec)), ' !eff'
    write(outlun,*)
  else
    write(lout,"(a)") "COLL> No particles absorbed"
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
    call f_requestUnit(coll_efficFile,coll_efficUnit)
    call f_open(unit=coll_efficUnit,file=coll_efficFile,formatted=.true.,mode="w")
    if(n_tot_absorbed /= 0) then
      write(coll_efficUnit,"(a1,1x,a13,6(1x,a15),1x,a8)") "#","rad_sigma",&
        "frac_x","frac_y","frac_r","eff_x","eff_y","eff_r","n_abs"
      do k=1,numeff
        write(coll_efficUnit,"(7(1x,e15.7),1x,i8)") rsig(k), neffx(k)/real(n_tot_absorbed,fPrec), &
          neffy(k)/real(n_tot_absorbed,fPrec), neff(k)/real(n_tot_absorbed,fPrec), neffx(k), neffy(k), neff(k), n_tot_absorbed
      end do
    else
      write(coll_efficUnit,"(a)") "No particles absorbed"
    end if
    call f_close(coll_efficUnit)
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
    call f_requestUnit(coll_efficDPFile,coll_efficDPUnit)
    call f_open(unit=coll_efficDPUnit,file=coll_efficDPFile,formatted=.true.,mode="w")
    if(n_tot_absorbed /= 0) then
      write(coll_efficDPUnit,"(a1,1x,a13,2(1x,a15),2(1x,a8))") "#","dp/p","n_dpop/tot_nabs","n_dpop","tot_nabs","npart"
      do k=1,numeffdpop
        write(coll_efficDPUnit,"(e15.7,2(1x,e15.7),2(1x,i8))") dpopbins(k), neffdpop(k)/real(n_tot_absorbed,fPrec), neffdpop(k), &
          n_tot_absorbed, npartdpop(k)
      end do
    else
      write(coll_efficDPUnit,"(a)") "No particles absorbed"
    end if
    call f_close(coll_efficDPUnit)
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
    call f_requestUnit(coll_effic2DFile,coll_effic2DUnit)
    call f_open(unit=coll_effic2DUnit,file=coll_effic2DFile,formatted=.true.,mode="w")
    if(n_tot_absorbed /= 0) then
      write(coll_effic2DUnit,"(a1,1x,a13,3(1x,a15),1x,a8)") "#","rad_sigma","dp/p","n/tot_nabs","n","tot_nabs"
      do i=1,numeff
        do k=1,numeffdpop
          write(coll_effic2DUnit,"(e15.7,3(1x,e15.7),1x,i8)") rsig(i), dpopbins(k), &
            neff2d(i,k)/real(n_tot_absorbed,fPrec), neff2d(i,k), n_tot_absorbed
        end do
      end do
    else
      write(coll_effic2DUnit,"(a)") "No particles absorbed"
    end if
    call f_close(coll_effic2DUnit)
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
      if(cdb_cLength(i) > zero .and. cdb_cFound(i)) then
        call h5_prepareWrite(setHdf, 1)
        call h5_writeData(setHdf, 1, 1, i)
        call h5_writeData(setHdf, 2, 1, cdb_cName(i))
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
    call f_requestUnit(coll_summaryFile,coll_summaryUnit)
    call f_open(unit=coll_summaryUnit,file=coll_summaryFile,formatted=.true.,mode="w")
    write(coll_summaryUnit,"(a1,1x,a5,1x,a20,2(1x,a8),2(1x,a15),1x,a6)") "#","icoll",chr_rPad("collname",20),&
      "nimp","nabs","imp_av","imp_sig","length"
    do icoll = 1, cdb_nColl
      if(cdb_cLength(icoll) > zero .and. cdb_cFound(icoll)) then
        write(coll_summaryUnit,"(i7,1x,a20,2(1x,i8),2(1x,e15.7),1x,f6.2)") icoll, cdb_cName(icoll)(1:20), cn_impact(icoll), &
          cn_absorbed(icoll), caverage(icoll), csigma(icoll), cdb_cLength(icoll)
      end if
    end do
    call f_close(coll_summaryUnit)
#ifdef HDF5
  end if
#endif

#ifdef ROOT
  if(root_flag .and. root_Collimation.eq.1) then
    do icoll = 1, cdb_nColl
      if(cdb_cLength(icoll).gt.zero) then
        call CollimatorLossRootWrite(icoll, cdb_cName(icoll), len(cdb_cName(icoll)), cn_impact(icoll), cn_absorbed(icoll), &
          caverage(icoll), csigma(icoll), cdb_cLength(icoll))
      end if
    end do
  end if
#endif

  call f_close(outlun)
  call f_close(coll_gapsUnit)
  call f_close(coll_positionsUnit)

  if(dowritetracks) then
    call f_close(coll_tracksUnit)
#ifdef HDF5
    if(h5_writeTracks2) call h5tr2_finalise
#endif
  end if

  if(do_select) then
    call f_close(coll_ellipseUnit)
  end if

  if(dowrite_impact) then
    call f_close(coll_allImpactUnit)
    call f_close(coll_allAbsorbUnit)
    call f_close(coll_flukImpUnit)
    call f_close(coll_flukImpAllUnit)
    call f_close(coll_scatterUnit)
    call f_close(coll_fstImpactUnit)
    call f_close(coll_jawProfileUnit)
  end if

  if(dowrite_amplitude) then
    ! Write amplitude.dat
    call f_requestUnit(coll_ampFile,coll_ampUnit)
    call f_open(unit=coll_ampUnit,file=coll_ampFile,formatted=.true.,mode="w")
    write(coll_ampUnit,"(a1,1x,a6,1x,a20,19(1x,a20))") "#","ielem",chr_rPad("name",20),"s","AX_AV","AX_RMS","AY_AV","AY_RMS",&
      "alphax","alphay","betax","betay","orbitx","orbity","dispx","dispy","xbob","ybob","xpbob","ypbob","mux","muy"
    do i=1,iu
      if(ic(i) <= nblo) then
        ix = mtyp(ic(i),mel(ic(i)))
      else
        ix = ic(i)-nblo
      end if
      write(coll_ampUnit,"(i8,1x,a20,19(1x,1pe20.13))") i, bez(ix)(1:20), dcum(i),                       &
        sum_ax(i)/real(max(nampl(i),1),fPrec),                                                           &
        sqrt(abs((sqsum_ax(i)/real(max(nampl(i),1),fPrec))-(sum_ax(i)/real(max(nampl(i),1),fPrec))**2)), &
        sum_ay(i)/real(max(nampl(i),1),fPrec),                                                           &
        sqrt(abs((sqsum_ay(i)/real(max(nampl(i),1),fPrec))-(sum_ay(i)/real(max(nampl(i),1),fPrec))**2)), &
        talphax(i), talphay(i), tbetax(i), tbetay(i), torbx(i), torby(i), tdispx(i), tdispy(i),          &
        xbob(i), ybob(i), xpbob(i), ypbob(i), mux(i), muy(i)
    end do
    call f_close(coll_ampUnit)
  end if

  ! Write orbitchecking.dat
  call f_requestUnit(coll_orbitCheckFile,coll_orbitCheckUnit)
  call f_open(unit=coll_orbitCheckUnit,file=coll_orbitCheckFile,formatted=.true.,mode="w")
  write(coll_orbitCheckUnit,"(a1,1x,a6,3(1x,a15))") "#","s","torbitx","torbity"
  do i=1,iu
    write(coll_orbitCheckUnit,"(i8,3(1x,1pe15.7))") i, dcum(i), torbx(i), torby(i)
  end do
  call f_close(coll_orbitCheckUnit)

#ifdef G4COLLIMATION
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

end subroutine collimate_start_turn

!>
!! This routine is called at the start of every element
!<
subroutine collimate_start_element(i)

  use mod_common
  use mod_common_main

  integer, intent(in) :: i
  integer j

  ie=i
#ifndef G4COLLIMATION
  ! For absorbed particles set all coordinates to zero. Also
  ! include very large offsets, let's say above 100mm or 100mrad.
  do j=1,napx
    if((part_abs_pos(j) /= 0 .and. part_abs_turn(j) /= 0) .or. &
      xv1(j) > c1e2 .or. yv1(j) > c1e2 .or. xv2(j) > c1e2 .or. yv2(j) > c1e2) then
      xv1(j)   = zero
      yv1(j)   = zero
      xv2(j)   = zero
      yv2(j)   = zero
      ejv(j)   = myenom
      sigmv(j) = zero
      secondary(j)  = 0
      tertiary(j)   = 0
      other(j)      = 0
      scatterhit(j) = 0
      nabs_type(j)  = 0
      part_abs_pos(j)  = ie
      part_abs_turn(j) = iturn
    end if
  end do
#endif

!GRD SAVE COORDINATES OF PARTICLE 1 TO CHECK ORBIT
  if(firstrun) then
    xbob(ie)=xv1(1)
    ybob(ie)=xv2(1)
    xpbob(ie)=yv1(1)
    ypbob(ie)=yv2(1)
  end if

!++  Here comes sixtrack stuff
  if(ic(i) <= nblo) then
    myix = mtyp(ic(i),mel(ic(i)))
  else
    myix = ic(i)-nblo
  end if

end subroutine collimate_start_element

!>
!! collimate_end_element()
!! This routine is called at the end of every element
!<
subroutine collimate_end_element

  use crcoall
  use parpro
  use coll_common
  use mod_common
  use mod_commons
  use mod_common_da
  use mod_common_main
  use mod_common_track
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif

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
             ((((xv1(j)*c1m3)**2 / (tbetax(ie)*myemitx0_collgap)) .ge. sigsecut2).or. &
             (((xv2(j)*c1m3)**2  / (tbetay(ie)*myemity0_collgap)) .ge. sigsecut2).or. &
             (((xv1(j)*c1m3)**2  / (tbetax(ie)*myemitx0_collgap)) + &
             ((xv2(j)*c1m3)**2  /  (tbetay(ie)*myemity0_collgap)) .ge. sigsecut3)) ) then

          xj  = (xv1(j)-torbx(ie)) /c1e3
          xpj = (yv1(j)-torbxp(ie))/c1e3
          yj  = (xv2(j)-torby(ie)) /c1e3
          ypj = (yv2(j)-torbyp(ie))/c1e3
#ifdef HDF5
          if(h5_writeTracks2) then
            call h5tr2_writeLine(partID(j),iturn,dcum(ie),xv1(j),yv1(j),xv2(j),yv2(j),&
              (ejv(j)-myenom)/myenom,secondary(j)+tertiary(j)+other(j)+scatterhit(j))
          else
#endif
            write(coll_tracksUnit,"(1x,i8,1x,i4,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)") partID(j), iturn, dcum(ie), &
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
  use coll_common
  use mod_common
  use mod_commons
  use mod_common_da
  use mod_common_main
  use mod_common_track
  use crcoall
  use mod_units
  use mathlib_bouncer

#ifdef ROOT
  use root_output
#endif
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif

  implicit none

  integer j, fUnit

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

  do j=1,napx
    xineff(j)  = xv1(j) - torbx (ie)
    xpineff(j) = yv1(j) - torbxp(ie)
    yineff(j)  = xv2(j) - torby (ie)
    ypineff(j) = yv2(j) - torbyp(ie)
    ! All particles absorbed are considered to be lost, so we give them a large offset
    if(part_abs_pos(j) /= 0 .and. part_abs_turn(j) /= 0) then
      xv1(j) = 99.5_fPrec
      xv2(j) = 99.5_fPrec
    end if
  end do

  ! For LAST ELEMENT in the ring calculate the number of surviving
  ! particles and save into file versus turn number
  if(ie == iu) then
    nsurvive = 0

    do j=1,napx
      if(xv1(j) < 99.0_fPrec .and. xv2(j) < 99.0_fPrec) then
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
      write(coll_survivalUnit,"(i7,1x,i9)") iturn, nsurvive
#ifdef HDF5
    end if
#endif

#ifdef ROOT
    if(root_flag .and. root_Collimation == 1) then
      call SurvivalRootWrite(iturn, nsurvive)
    end if
#endif

    if(iturn == numl) then
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
  if(do_coll .and. ((iturn == 1 .and. ie == 20) .or. (ie == iu))) then
    gammax = (1 + talphax(ie)**2)/tbetax(ie)
    gammay = (1 + talphay(ie)**2)/tbetay(ie)

    do j = 1, napx
      ! Do the binning in amplitude, only considering particles that were not absorbed before.
      if(xv1(j) < 99.0_fPrec .and. xv2(j) < 99.0_fPrec .and. (part_select(j) == 1 .or. ie == 20)) then
        ! Normalized amplitudes are calculated
        ! Allow to apply some dispersive offset. Take arc dispersion (2m) and normalize with arc beta_x function (180m).
        arcdx     = 2.5_fPrec
        arcbetax  = c180e0
        xdisp     = abs(xv1(j)*c1m3) + (abs((ejv(j)-myenom)/myenom)*arcdx) * sqrt(tbetax(ie)/arcbetax)
        nspx      = sqrt(                                                       &
                      abs(gammax*xdisp**2 +                                     &
                        ((two*talphax(ie))*xdisp)*(yv1(j)*c1m3) +               &
                        tbetax(ie)*(yv1(j)*c1m3)**2                             &
                      )/myemitx0_collgap                                        &
                    )
        nspy      = sqrt(                                                       &
                      abs(gammay*(xv2(j)*c1m3)**2 +                             &
                        ((two*talphay(ie))*(xv2(j)*c1m3))*(yv2(j)*c1m3) +       &
                        tbetay(ie)*(yv2(j)*c1m3)**2                             &
                      )/myemity0_collgap                                        &
                    )

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
          xnorm  = (xv1(j)*c1m3) / sqrt(tbetax(ie)*myemitx0_collgap)
          xpnorm = (talphax(ie)*(xv1(j)*c1m3)+ tbetax(ie)*(yv1(j)*c1m3)) / sqrt(tbetax(ie)*myemitx0_collgap)
          xangle = atan2_mb(xnorm,xpnorm)
          xnorm  = xnorm  + dnormx*sin_mb(xangle)
          xpnorm = xpnorm + dnormx*cos_mb(xangle)
          xv1(j) = c1e3 * (xnorm * sqrt(tbetax(ie)*myemitx0_collgap))
          yv1(j) = c1e3 * ((xpnorm*sqrt(tbetax(ie)*myemitx0_collgap)-talphax(ie)*xv1(j)*c1m3)/tbetax(ie))

          ynorm  = (xv2(j)*c1m3)/ sqrt(tbetay(ie)*myemity0_collgap)
          ypnorm = (talphay(ie)*(xv2(j)*c1m3)+tbetay(ie)*(yv2(j)*c1m3)) / sqrt(tbetay(ie)*myemity0_collgap)
          yangle = atan2_mb(ynorm,ypnorm)
          ynorm  = ynorm  + dnormy*sin_mb(yangle)
          ypnorm = ypnorm + dnormy*cos_mb(yangle)
          xv2(j) = c1e3 * (ynorm * sqrt(tbetay(ie)*myemity0_collgap))
          yv2(j) = c1e3 * ((ypnorm*sqrt(tbetay(ie)*myemity0_collgap)-talphay(ie)*xv2(j)*c1m3)/tbetay(ie))
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
  if(ie == iu) then
    do j = 1, napx
      if(xv1(j) < 99.0_fPrec .and. xv2(j) < 99.0_fPrec) then
        llostp(j) = .false.
      else
        llostp(j) = .true.
      end if
    end do

    ! Move the lost particles to the end of the arrays
    call shuffleLostParticles
  end if

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
      call h5_writeData(setHdf, 1, napx, (xv1(1:napx)-torbx(1)) /c1e3)
      call h5_writeData(setHdf, 2, napx, (yv1(1:napx)-torbxp(1))/c1e3)
      call h5_writeData(setHdf, 3, napx, (xv2(1:napx)-torby(1)) /c1e3)
      call h5_writeData(setHdf, 4, napx, (yv2(1:napx)-torbyp(1))/c1e3)
      call h5_writeData(setHdf, 5, napx, sigmv(1:napx))
      call h5_writeData(setHdf, 6, napx, ejfv(1:napx))
      call h5_finaliseWrite(setHdf)
      deallocate(fldHdf)
    else
#endif
      call f_requestUnit("distn.dat",fUnit)
      call f_open(unit=fUnit,file="distn.dat",formatted=.true.,mode="w",status="replace")
      write(fUnit,"(a)") "# 1=x 2=xp 3=y 4=yp 5=z 6=E"
      do j=1,napx
        write(fUnit,"(6(1x,e23.15))") (xv1(j)-torbx(1))/c1e3, (yv1(j)-torbxp(1))/c1e3, (xv2(j)-torby(1))/c1e3, &
          (yv2(j)-torbyp(1))/c1e3, sigmv(j), ejfv(j)
      end do
      call f_close(fUnit)
#ifdef HDF5
    end if
#endif
  end if

  if(firstrun) then
    if(rselect > 0 .and. rselect < 65) then ! The value rselect is fixed to 64, this if is probably redundant.
      do j=1,napx
        xj  = (xv1(j)-torbx(ie)) /c1e3
        xpj = (yv1(j)-torbxp(ie))/c1e3
        yj  = (xv2(j)-torby(ie)) /c1e3
        ypj = (yv2(j)-torbyp(ie))/c1e3
        pj  = ejv(j)/c1e3

        if(iturn == 1 .and. j == 1) then
          sum_ax(ie) = zero
          sum_ay(ie) = zero
        end if

        if(tbetax(ie) > 0.) then
          gammax = (one + talphax(ie)**2)/tbetax(ie)
          gammay = (one + talphay(ie)**2)/tbetay(ie)
        else
          gammax = (one + talphax(ie-1)**2)/tbetax(ie-1)
          gammay = (one + talphay(ie-1)**2)/tbetay(ie-1)
        end if

        if(part_abs_pos(j) == 0 .and. part_abs_turn(j) == 0) then
          if(tbetax(ie) > 0.) then
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
        else
          nspx = zero
          nspy = zero
        end if !if(part_abs_pos(j).eq.0 .and. part_abs_turn(j).eq.0) then
      end do !do j = 1, napx
    end if !if(rselect.gt.0 .and. rselect.lt.65) then
  end if !if(firstrun) then

!GRD THIS LOOP MUST NOT BE WRITTEN INTO THE "IF(FIRSTRUN)" LOOP !!!!
  if(dowritetracks) then
    do j=1, napx
      xj    = (xv1(j)-torbx(ie))/c1e3
      xpj   = (yv1(j)-torbxp(ie))/c1e3
      yj    = (xv2(j)-torby(ie))/c1e3
      ypj   = (yv2(j)-torbyp(ie))/c1e3
      arcdx = 2.5_fPrec
      arcbetax = c180e0

      if(xj <= 0.) then
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
            ((((xv1(j)*c1m3)**2 / (tbetax(ie)*myemitx0_collgap)) .ge. sigsecut2).or. &
            (((xv2(j)*c1m3)**2  / (tbetay(ie)*myemity0_collgap)) .ge. sigsecut2).or. &
            (((xv1(j)*c1m3)**2  / (tbetax(ie)*myemitx0_collgap)) + &
            ((xv2(j)*c1m3)**2  / (tbetay(ie)*myemity0_collgap)) .ge. sigsecut3)) ) then

          xj  = (xv1(j)-torbx(ie))/c1e3
          xpj = (yv1(j)-torbxp(ie))/c1e3
          yj  = (xv2(j)-torby(ie))/c1e3
          ypj = (yv2(j)-torbyp(ie))/c1e3
#ifdef HDF5
          if(h5_writeTracks2) then
            call h5tr2_writeLine(partID(j),iturn,dcum(ie),xv1(j),yv1(j),xv2(j),yv2(j),&
              (ejv(j)-myenom)/myenom,secondary(j)+tertiary(j)+other(j)+scatterhit(j))
          else
#endif
            write(coll_tracksUnit,"(1x,i8,1x,i4,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)") partID(j),iturn,dcum(ie), &
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

#ifdef HDF5
subroutine coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs,dp,dx,dy)

  use hdf5_output

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

end module collimation
