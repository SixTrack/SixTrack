! ================================================================================================ !
!
!  SixTrack Collimation Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  G. Robert-Demolaize, R. Assmann, T. Weiler, C. Bracco, S. Redaelli, R. Bruce, D. Mirarchi,
!  Y.I. Levinsen, C. Tambasco, J. Molson, K.N. Sjobak, V.K. Berglyd Olsen
!
!  Created: 2004-10-27
!  Updated: 2019-10-21
!
! ================================================================================================ !
module collimation

  use parpro
  use floatPrecision
  use numerical_constants

  implicit none

  integer, parameter :: numeff     = 32
  integer, parameter :: numeffdpop = 29

  ! Logical Flags
  logical, public,  save :: do_coll          = .false.
  logical, private, save :: coll_oldBlock    = .false.
  logical, private, save :: do_select        = .false.
  logical, private, save :: do_nominal       = .false.
  logical, private, save :: do_oneside       = .false.
  logical, private, save :: systilt_antisymm = .false.
  logical, private, save :: do_mingap        = .false.
  logical, private, save :: firstcoll        = .true.

  integer, private, save :: icoll            = 0
  integer, private, save :: ie               = 1
  integer, private, save :: iturn            = 1
  integer, private, save :: c_ix             = 0

  ! Distribution
  integer,          private, save :: do_thisdis   = 0
  real(kind=fPrec), public,  save :: c_enom       = zero
  logical,          private, save :: radial       = .false.

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

  ! Emittance Drift
  real(kind=fPrec), private, save :: driftsx = zero
  real(kind=fPrec), private, save :: driftsy = zero

  real(kind=fPrec), private, save :: sigsecut3 = one
  real(kind=fPrec), private, save :: sigsecut2 = one

  ! Normalised emittances from input
  real(kind=fPrec), private, save :: emitnx0_dist    = zero
  real(kind=fPrec), private, save :: emitny0_dist    = zero
  real(kind=fPrec), private, save :: emitnx0_collgap = zero
  real(kind=fPrec), private, save :: emitny0_collgap = zero

  ! Geometric emittances
  real(kind=fPrec), private, save :: c_emitx0_dist    = zero
  real(kind=fPrec), private, save :: c_emity0_dist    = zero
  real(kind=fPrec), private, save :: c_emitx0_collgap = zero
  real(kind=fPrec), private, save :: c_emity0_collgap = zero

  character(len=mNameLen), private, save :: name_sel = " "

  ! M. Fiascaris for the collimation team
  ! Variables for global inefficiencies studies of normalized and off-momentum halo
  integer,          allocatable, private, save :: counteddpop(:,:) ! (npart,numeffdpop)
  integer,          allocatable, private, save :: npartdpop(:)     ! (numeffdpop)
  integer,          allocatable, private, save :: counted2d(:,:,:) ! (npart,numeff,numeffdpop)
  integer,          allocatable, private, save :: counted_r(:,:)   ! (npart,numeff)
  integer,          allocatable, private, save :: counted_x(:,:)   ! (npart,numeff)
  integer,          allocatable, private, save :: counted_y(:,:)   ! (npart,numeff)
  real(kind=fPrec), allocatable, private, save :: neffx(:)         ! (numeff)
  real(kind=fPrec), allocatable, private, save :: neffy(:)         ! (numeff)
  real(kind=fPrec), allocatable, private, save :: neff(:)          ! (numeff)
  real(kind=fPrec), allocatable, private, save :: rsig(:)          ! (numeff)
  real(kind=fPrec), allocatable, private, save :: neffdpop(:)      ! (numeffdpop)
  real(kind=fPrec), allocatable, private, save :: dpopbins(:)      ! (numeffdpop)
  real(kind=fPrec), allocatable, private, save :: neff2d(:,:)      ! (numeff,numeffdpop)

  ! Arrays allocated to nblz
  integer,          allocatable, private, save :: nampl(:)
  real(kind=fPrec), allocatable, private, save :: sum_ax(:)
  real(kind=fPrec), allocatable, private, save :: sqsum_ax(:)
  real(kind=fPrec), allocatable, private, save :: sum_ay(:)
  real(kind=fPrec), allocatable, private, save :: sqsum_ay(:)

  ! Arrays allocated to npart
  integer,          allocatable, private, save :: part_hit_pos(:)   ! Hit flag for last hit
  integer,          allocatable, private, save :: part_hit_turn(:)  ! Hit flag for last hit
  integer,          allocatable, public,  save :: part_abs_pos(:)   ! Absorbed in element
  integer,          allocatable, public,  save :: part_abs_turn(:)  ! Absorbed in turn
  integer,          allocatable, private, save :: part_select(:)
  integer,          allocatable, private, save :: nabs_type(:)
  integer,          allocatable, public,  save :: nhit_stage(:)
  real(kind=fPrec), allocatable, private, save :: part_linteract(:)
  real(kind=fPrec), allocatable, private, save :: part_indiv(:)     ! Divergence of impacting particles
  real(kind=fPrec), allocatable, private, save :: part_impact(:)    ! Impact parameter (0 for inner face)
  real(kind=fPrec), allocatable, private, save :: rcx0(:)
  real(kind=fPrec), allocatable, private, save :: rcxp0(:)
  real(kind=fPrec), allocatable, private, save :: rcy0(:)
  real(kind=fPrec), allocatable, private, save :: rcyp0(:)
  real(kind=fPrec), allocatable, private, save :: rcp0(:)

  integer, private, save :: n_tot_absorbed = 0
  integer, private, save :: n_absorbed     = 0
  integer, private, save :: nabs_total     = 0
  integer, private, save :: nsurvive       = 0
  integer, private, save :: nsurvive_end   = 0
  integer, private, save :: num_selhit     = 0

  ! Variables for finding the collimator with the smallest gap and defining, storing the gap rms error
  real(kind=fPrec), allocatable, private, save :: xbob(:)    ! (nblz)
  real(kind=fPrec), allocatable, private, save :: ybob(:)    ! (nblz)
  real(kind=fPrec), allocatable, private, save :: xpbob(:)   ! (nblz)
  real(kind=fPrec), allocatable, private, save :: ypbob(:)   ! (nblz)
  real(kind=fPrec), allocatable, private, save :: xineff(:)  ! (npart)
  real(kind=fPrec), allocatable, private, save :: yineff(:)  ! (npart)
  real(kind=fPrec), allocatable, private, save :: xpineff(:) ! (npart)
  real(kind=fPrec), allocatable, private, save :: ypineff(:) ! (npart)

  real(kind=fPrec), private, save :: scale_bx0 = zero
  real(kind=fPrec), private, save :: scale_by0 = zero
  real(kind=fPrec), private, save :: bx_dist   = zero
  real(kind=fPrec), private, save :: by_dist   = zero

#ifdef CR
! CR arrays
  integer,          allocatable, private, save :: counteddpop_cr(:,:) ! (npart,numeffdpop)
  integer,          allocatable, private, save :: npartdpop_cr(:)     ! (numeffdpop)
  integer,          allocatable, private, save :: counted2d_cr(:,:,:) ! (npart,numeff,numeffdpop)
  integer,          allocatable, private, save :: counted_r_cr(:,:)   ! (npart,numeff)
  integer,          allocatable, private, save :: counted_x_cr(:,:)   ! (npart,numeff)
  integer,          allocatable, private, save :: counted_y_cr(:,:)   ! (npart,numeff)
  real(kind=fPrec), allocatable, private, save :: neffx_cr(:)         ! (numeff)
  real(kind=fPrec), allocatable, private, save :: neffy_cr(:)         ! (numeff)
  real(kind=fPrec), allocatable, private, save :: neff_cr(:)          ! (numeff)
  real(kind=fPrec), allocatable, private, save :: rsig_cr(:)          ! (numeff)
  real(kind=fPrec), allocatable, private, save :: neffdpop_cr(:)      ! (numeffdpop)
  real(kind=fPrec), allocatable, private, save :: dpopbins_cr(:)      ! (numeffdpop)
  real(kind=fPrec), allocatable, private, save :: neff2d_cr(:,:)      ! (numeff,numeffdpop)

  ! Arrays allocated to nblz
  integer,          allocatable, private, save :: nampl_cr(:)
  real(kind=fPrec), allocatable, private, save :: sum_ax_cr(:)
  real(kind=fPrec), allocatable, private, save :: sqsum_ax_cr(:)
  real(kind=fPrec), allocatable, private, save :: sum_ay_cr(:)
  real(kind=fPrec), allocatable, private, save :: sqsum_ay_cr(:)
  real(kind=fPrec), allocatable, private, save :: xbob_cr(:)
  real(kind=fPrec), allocatable, private, save :: ybob_cr(:)
  real(kind=fPrec), allocatable, private, save :: xpbob_cr(:)
  real(kind=fPrec), allocatable, private, save :: ypbob_cr(:)

  ! Arrays allocated to npart
  integer,          allocatable, private, save :: part_hit_pos_cr(:)   ! Hit flag for last hit
  integer,          allocatable, private, save :: part_hit_turn_cr(:)  ! Hit flag for last hit
  integer,          allocatable, public,  save :: part_abs_pos_cr(:)   ! Absorbed in element
  integer,          allocatable, public,  save :: part_abs_turn_cr(:)  ! Absorbed in turn
  integer,          allocatable, private, save :: part_select_cr(:)
  integer,          allocatable, private, save :: nabs_type_cr(:)
  integer,          allocatable, public,  save :: nhit_stage_cr(:)
  real(kind=fPrec), allocatable, private, save :: part_linteract_cr(:)
  real(kind=fPrec), allocatable, private, save :: part_indiv_cr(:)     ! Divergence of impacting particles
  real(kind=fPrec), allocatable, private, save :: part_impact_cr(:)    ! Impact parameter (0 for inner face)

  real(kind=fPrec), allocatable, private, save :: xineff_cr(:)  ! (npart)
  real(kind=fPrec), allocatable, private, save :: yineff_cr(:)  ! (npart)
  real(kind=fPrec), allocatable, private, save :: xpineff_cr(:) ! (npart)
  real(kind=fPrec), allocatable, private, save :: ypineff_cr(:) ! (npart)

  integer, allocatable, save :: cry_proc_cr(:)
  integer, allocatable, save :: cry_proc_prev_cr(:)
  integer, allocatable, save :: cry_proc_tmp_cr(:)

  integer, private, save :: n_tot_absorbed_cr = 0
  integer, private, save :: n_absorbed_cr     = 0
  integer, private, save :: nabs_total_cr     = 0
  integer, private, save :: nsurvive_cr       = 0
  integer, private, save :: nsurvive_end_cr   = 0
  integer, private, save :: num_selhit_cr     = 0
#endif

contains

subroutine collimation_expand_arrays(npart_new, nblz_new)

  use mod_alloc
  use coll_common

  integer, intent(in) :: npart_new
  integer, intent(in) :: nblz_new

  logical :: efficAlloc = .false.

  ! Arrays that are always needed
  call alloc(part_abs_turn, npart_new, 0, "part_abs_turn")

  if(.not. do_coll) return

  ! Arrays below are only needed if collimation is enabled

  ! Allocate Common Variables
  call coll_expandArrays(npart_new)

  call alloc(rcx0,           npart_new,  zero, "rcx0")
  call alloc(rcxp0,          npart_new,  zero, "rcxp0")
  call alloc(rcy0,           npart_new,  zero, "rcy0")
  call alloc(rcyp0,          npart_new,  zero, "rcyp0")
  call alloc(rcp0,           npart_new,  zero, "rcp0")

  call alloc(part_hit_pos,   npart_new,  0,    "part_hit_pos")
  call alloc(part_hit_turn,  npart_new,  0,    "part_hit_turn")
  call alloc(part_abs_pos,   npart_new,  0,    "part_abs_pos")
  call alloc(part_select,    npart_new,  1,    "part_select")
  call alloc(nabs_type,      npart_new,  0,    "nabs_type")
  call alloc(nhit_stage,     npart_new,  0,    "nhit_stage")

  call alloc(part_impact,    npart_new,  zero, "part_impact")
  call alloc(part_indiv,     npart_new, -c1m6, "part_indiv")
  call alloc(part_linteract, npart_new,  zero, "part_linteract")

  if(dowrite_amplitude) then
    call alloc(nampl,    nblz_new,  0,    "nampl")
    call alloc(sum_ax,   nblz_new,  zero, "sum_ax")
    call alloc(sqsum_ax, nblz_new,  zero, "sqsum_ax")
    call alloc(sum_ay,   nblz_new,  zero, "sum_ay")
    call alloc(sqsum_ay, nblz_new,  zero, "sqsum_ay")
    call alloc(xbob,     nblz_new,  zero, "xbob")
    call alloc(ybob,     nblz_new,  zero, "ybob")
    call alloc(xpbob,    nblz_new,  zero, "xpbob")
    call alloc(ypbob,    nblz_new,  zero, "ypbob")
  end if

  if(dowrite_efficiency) then
    call alloc(xineff,      npart_new,                     zero, "xineff")
    call alloc(yineff,      npart_new,                     zero, "yineff")
    call alloc(xpineff,     npart_new,                     zero, "xpineff")
    call alloc(ypineff,     npart_new,                     zero, "ypineff")
    call alloc(counted_r,   npart_new, numeff,             0,    "counted_r")
    call alloc(counted_x,   npart_new, numeff,             0,    "counted_x")
    call alloc(counted_y,   npart_new, numeff,             0,    "counted_y")
    call alloc(counteddpop, npart_new,         numeffdpop, 0,    "counteddpop")
    call alloc(counted2d,   npart_new, numeff, numeffdpop, 0,    "counted2d")
  end if

  if(dowrite_efficiency .and. .not. efficAlloc) then
    ! These are fixed size, so only need to be allocated once
    call alloc(npartdpop,         numeffdpop, 0,    "npartdpop")
    call alloc(neff,      numeff,             zero, "neff")
    call alloc(rsig,      numeff,             zero, "rsig")
    call alloc(neffdpop,          numeffdpop, zero, "neffdpop")
    call alloc(dpopbins,          numeffdpop, zero, "dpopbins")
    call alloc(neff2d,    numeff, numeffdpop, zero, "neff2d")
    call alloc(neffx,     numeff,             zero, "neffx")
    call alloc(neffy,     numeff,             zero, "neffy")
    efficAlloc = .true.
  end if

end subroutine collimation_expand_arrays

subroutine coll_shuffleLostPart

  use mod_common
  use mod_common_main
  use coll_common

  integer j, tnapx

  tnapx = napx
  do j=napx,1,-1
    if(llostp(j) .eqv. .false.) cycle

    part_hit_pos(j:tnapx)   = cshift(part_hit_pos(j:tnapx),   1)
    part_hit_turn(j:tnapx)  = cshift(part_hit_turn(j:tnapx),  1)
    part_abs_pos(j:tnapx)   = cshift(part_abs_pos(j:tnapx),   1)
    part_abs_turn(j:tnapx)  = cshift(part_abs_turn(j:tnapx),  1)
    part_select(j:tnapx)    = cshift(part_select(j:tnapx),    1)
    part_impact(j:tnapx)    = cshift(part_impact(j:tnapx),    1)
    part_indiv(j:tnapx)     = cshift(part_indiv(j:tnapx),     1)
    part_linteract(j:tnapx) = cshift(part_linteract(j:tnapx), 1)

    nabs_type(j:tnapx)      = cshift(nabs_type(j:tnapx),      1)
    nhit_stage(j:tnapx)     = cshift(nhit_stage(j:tnapx),     1)

    cry_proc_prev(j:tnapx)  = cshift(cry_proc_prev(j:tnapx),  1)
    cry_proc_tmp(j:tnapx)   = cshift(cry_proc_tmp(j:tnapx),   1)

    tnapx = tnapx - 1
  end do

  if(dowrite_efficiency) then
    tnapx = napx
    do j=napx,1,-1
      if(llostp(j) .eqv. .false.) cycle

      counted_r(j:tnapx,:) = cshift(counted_r(j:tnapx,:), 1, 1)
      counted_x(j:tnapx,:) = cshift(counted_x(j:tnapx,:), 1, 1)
      counted_y(j:tnapx,:) = cshift(counted_y(j:tnapx,:), 1, 1)

      tnapx = tnapx - 1
    end do
  end if

end subroutine coll_shuffleLostPart

! ================================================================================================ !
!  Collimation Init
!  This routine is called once at the start of the simulation and can be used to do any initial
!  configuration and/or file loading.
! ================================================================================================ !
subroutine coll_init

  use crcoall
  use mod_common
  use mod_common_main
  use mod_common_track
  use coll_k2
  use coll_db
  use coll_dist
  use coll_common
  use coll_crystal
  use coll_materials
  use mod_time
  use mod_units
  use mod_ranlux
  use mod_particles
#ifdef HDF5
  use hdf5_output
#endif
#ifdef G4COLLIMATION
  use geant4
#endif

#ifdef HDF5
  type(h5_dataField), allocatable :: fldDist0(:)
  integer                         :: fmtDist0, setDist0
#endif
  integer i,j,fUnit,minGapID
  real(kind=fPrec) dummy,c_rmstilt,c_systilt

#ifdef HDF5
  if(h5_useForCOLL) call h5_initForCollimation
#endif

  c_emitx0_dist    = (emitnx0_dist*gammar)*c1m6
  c_emity0_dist    = (emitny0_dist*gammar)*c1m6
  c_emitx0_collgap = (emitnx0_collgap*gammar)*c1m6
  c_emity0_collgap = (emitny0_collgap*gammar)*c1m6

  if(c_emitx0_dist <= zero .or. c_emity0_dist <= zero .or. c_emitx0_collgap <= zero .or. c_emity0_collgap <= zero) then
    write(lerr,"(a)") "COLL> ERROR Emittances not defined! check collimat block!"
    write(lerr,"(a)") "COLL> ERROR Expected format of line 9 in collimat block:"
    write(lerr,"(a)") "COLL> ERROR emitnx0_dist  emitny0_dist  emitnx0_collgap  emitny0_collgap"
    write(lerr,"(a)") "COLL> ERROR All emittances should be normalized. "//&
      "first put emittance for distribtion generation, then for collimator position etc. units in [mm*mrad]."
    write(lerr,"(a)") "COLL> ERROR EXAMPLE: 2.5 2.5 3.5 3.5"
    call prror
  end if

  call coll_echoSettings

  ! Initialize random number generator
  if(rnd_seed == 0) rnd_seed = time_getSysClock()
  if(rnd_seed <  0) rnd_seed = abs(rnd_seed)
  call rluxgo(3, rnd_seed, 0, 0)
  write(outlun,*) 'INFO>  rnd_seed: ', rnd_seed
  flush(outlun)
#ifdef CR
  coll_trackoutPos = coll_trackoutPos + 1
#endif

  ! Call distribution routines
  call cdist_init(c_enom,c_emitx0_dist,c_emity0_dist,c_emitx0_collgap,c_emity0_collgap)
  if(radial) then
    call cdist_makeRadial
  else
    call cdist_makeDist(do_thisdis)
  end if

  ! Reset distribution for pencil beam
  if(ipencil > 0) then
    write(lout,"(a)") "COLL> WARNING Distributions reset to pencil beam!"
    write(outlun,*)   "WARN> Distributions reset to pencil beam!"
    flush(outlun)
#ifdef CR
    coll_trackoutPos = coll_trackoutPos + 1
#endif
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
      call f_open(unit=fUnit,file="dist0.dat",formatted=.true.,mode="w",status="replace")
      do j=1, napx
        write(fUnit,"(6(1x,e23.15))") xv1(j), yv1(j), xv2(j), yv2(j), sigmv(j), ejv(j)
      end do
      call f_close(fUnit)
#ifdef HDF5
    end if
#endif
  end if

  ! Collimator Database and Materials
  call collmat_init                   ! Set default values for collimator materials
  call cdb_readCollDB                 ! Read the collimator DB
  call cdb_setLHCOnesided(do_oneside) ! Set LHC onesided collimators
  call cdb_writeDB_newFromOld         ! Write a copy of the db in new format, if provided in old format

  ! Then do any implementation specific initial loading
  call k2coll_init
  if(coll_hasCrystal) then
    call cry_init
  end if
#ifdef MERLINSCATTER
  call k2coll_merlinInit
#endif

  ! Open the edep file
  if(unit208 == -1) then
    call f_requestUnit(fort208,unit208)
    call f_open(unit=unit208,file=fort208,formatted=.true.,mode="w",status="replace")
#ifdef CR
    fort208Pos = 0
#endif
  end if

#ifdef G4COLLIMATION

  if(n_slices /= 0) then
    write(lerr,"(a)") "COLL> ERROR Cannot use jaw fit in G4COLLIMATION version of SixTrack"
    call prror
  end if

  ! This function lives in the G4Interface.cpp file in the g4collimat folder
  ! Accessed by linking libg4collimat.a
  ! Set the energy cut at 70% - i.e. 30% energy loss
  ! g4_ecut = 0.7_fPrec
  ! g4_ecut = zero

  ! Select the physics engine to use
  ! 0 = FTFP_BERT
  ! 1 = QGSP_BERT
  ! g4_physics = 0

  call g4_collimation_init(e0, rnd_seed, g4_recut, g4_aecut, g4_rcut, g4_rangecut_mm, g4_v0, trim(g4_phys_str), &
    g4_debug, g4_keep_stable, g4_edep, g4_neutral)
#endif

  write (lout,"(a)") ""
  write (lout,"(a)") "COLL> Finished collimate initialisation"
  write (lout,"(a)") ""

  ! Adding the orbit offset at start of ring
  if(do_thisdis /= 0 .or. radial) then
    xv1(1:napx) = c1e3*xv1(1:napx) + torbx(1)
    yv1(1:napx) = c1e3*yv1(1:napx) + torbxp(1)
    xv2(1:napx) = c1e3*xv2(1:napx) + torby(1)
    yv2(1:napx) = c1e3*yv2(1:napx) + torbyp(1)
  end if

  call part_updatePartEnergy(1,.false.)
  call coll_openFiles

  ! Initialisation
  if(dowrite_efficiency) then
    do i=1,numeff
      rsig(i) = (real(i,fPrec)/two - half) + five
    end do
    dpopbins(1) = c1m4
    do i=2,numeffdpop
      dpopbins(i) = real(i-1,fPrec)*4.0e-4_fPrec
    end do
  end if

#ifdef BEAMGAS
  call beamGasInit(c_enom)
#endif

  write(lout,"(a)") ""
  write(lout,"(a,i0)") "COLL> Number of collimators: ",cdb_nColl
  do i=1,cdb_nColl
    if(cdb_cFound(i)) then
      write(lout,"(a,i5,a)") "COLL> Found Collimator   ",i,": "//trim(cdb_cName(i))
    else
      write(lout,"(a,i5,a)") "COLL> Missing Collimator ",i,": "//trim(cdb_cName(i))
    end if
  end do
  write(lout,"(a)") ""

  ! Write settings for alignment error in colltrack.out file
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
  flush(outlun)
#ifdef CR
  coll_trackoutPos = coll_trackoutPos + 15
#endif

  ! Intialise random generator with offset_seed
  c_offsettilt_seed = abs(c_offsettilt_seed)
  call rluxgo(3, c_offsettilt_seed, 0, 0)

  ! Generate random tilts (Gaussian distribution plus systematic)
  if(c_rmstilt_prim > zero .or. c_rmstilt_sec > zero .or. c_systilt_prim /= zero .or. c_systilt_sec /= zero) then
    do i=1,cdb_nColl
      if(cdb_cStage(i) == cdb_stgPrimary) then
        c_rmstilt = c_rmstilt_prim
        c_systilt = c_systilt_prim
      else
        c_rmstilt = c_rmstilt_sec
        c_systilt = c_systilt_sec
      end if
      cdb_cTilt(1,i) = c_systilt + c_rmstilt*ran_gauss2(three)
      if(systilt_antisymm) then
        cdb_cTilt(2,i) = -one*c_systilt + c_rmstilt*ran_gauss2(three)
      else
        cdb_cTilt(2,i) =      c_systilt + c_rmstilt*ran_gauss2(three)
      end if
      write(outlun,*) 'INFO>  Collimator ',trim(cdb_cName(i)),' jaw 1 has tilt [rad]: ',cdb_cTilt(1,i)
      write(outlun,*) 'INFO>  Collimator ',trim(cdb_cName(i)),' jaw 2 has tilt [rad]: ',cdb_cTilt(2,i)
    end do

    do i=1,cdb_nColl
      if(cdb_cStage(i) == cdb_stgPrimary) then
        cdb_cOffset(i) = c_sysoffset_prim + c_rmsoffset_prim*ran_gauss2(three)
      else
        cdb_cOffset(i) = c_sysoffset_sec  + c_rmsoffset_sec*ran_gauss2(three)
      end if
      write(outlun,*) 'INFO>  Offset: ',trim(cdb_cName(i)),cdb_cOffset(i)
    end do
  end if
  if(c_rmserror_gap > zero) then
    do i=1,cdb_nColl
      gap_rms_error(i) = c_rmserror_gap * ran_gauss2(three)
      write(outlun,*) 'INFO>  Gap RMS error: ',trim(cdb_cName(i)),gap_rms_error(i)
    end do
  end if

  ! In case we're using old type jaw fit, this is where we generate the parameters for the new method
  ! After this, the number of slices is also stored per collimator, and can be extracted again later
  call cdb_setMasterJawFit(n_slices, smin_slices, smax_slices, recenter1, recenter2, jaw_fit, jaw_ssf)

  call coll_getMinGapID(minGapID)

  ! if pencil beam is used and on collimator with smallest gap the
  ! distribution should be generated, set ipencil to minGapID
  if(ipencil > 0 .and. do_mingap) then
    ipencil = minGapID
  end if

  ! This sets the random geenrator back to the default seed rather than the one used for coll gaps.
  ! However, this doesn't actually restore the random generator to the state it would have been in without the
  ! coll gaps errors being generated as rndm5() will extract 30000 new random numbers from ranlux and discard
  ! the ones it already holds and would have used.
  ! Alternatively, we can use ranecu instead, which is capable of continuing a chain of random numbers from
  ! a given set of seeds.
  ! It is probably unnecessary to use different random seeds here in the first place.
  call rluxgo(3, rnd_seed, 0, 0)
!  dummy = coll_rand() ! Reset rndm5 too

end subroutine coll_init

! ================================================================================================ !
!  Parse Input Line
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-04-16
! ================================================================================================ !
subroutine coll_parseInputLine(inLine, iLine, iErr)

  use crcoall
  use coll_db
  use coll_dist
  use string_tools
  use coll_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  real(kind=fPrec) nSigIn(23), rTmp
  integer nSplit, famID, iDum
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
    call chr_cast(lnSplit(2), c_enom, iErr)

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

  case("WRITE_EFFIC","WRITE_EFFICIENCY")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR WRITE_EFFIC expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       WRITE_EFFIC true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dowrite_efficiency, iErr)

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
    write(lerr,"(a)") "COLL> ERROR The BEAM_NUM flag has been removed"
    iErr = .true.
    return

  case("WRITE_TRACKS")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR WRITE_TRACKS expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       WRITE_TRACKS true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dowrite_tracks, iErr)

  case("WRITE_CRYCOORDS")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "COLL> ERROR WRITE_CRYCOORDS expects 1 value, got ",nSplit-1
      write(lerr,"(a)")    "COLL>       WRITE_CRYCOORDS true|false"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), dowrite_crycoord, iErr)

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
    call chr_cast(lnSplit(1),iDum,iErr)
    call chr_cast(lnSplit(2),c_enom,iErr)

    if(iDum /= 1) then
      write(lerr,"(a,i0)") "COLL> ERROR Multiple samples is no longer supported. nloop must be 1, got ",iDum
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
    if(nSplit /= 9 .and. nSplit /= 10) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 9 or 10 values on line 10, got ",nSplit
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
    if(nSplit > 9) then
      ! This one is optional as it's been added later
      call chr_cast(lnSplit(10), dowrite_efficiency,iErr)
    end if

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
    ! Value 3 is ignored
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
    cdb_fileName = lnSplit(1)
    ! The second value is ignored

  case(17)
    if(nSplit /= 6) then
      write(lerr,"(a,i0)") "COLL> ERROR Expected 6 values on line 17, got ",nSplit
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(1), dowrite_tracks,iErr)
    ! The second value is ignored
    ! The third value is ignored
    ! The fourth value is ignored
    call chr_cast(lnSplit(5), sigsecut2, iErr)
    call chr_cast(lnSplit(6), sigsecut3, iErr)

  case default
    write(lerr,"(a,i0,a)") "COLL> ERROR Unexpected line ",iLine," encountered."
    iErr = .true.

  end select

end subroutine coll_parseInputLine

! ================================================================================================ !
!  Post-input checks for the sanity of parameters
! ================================================================================================ !
subroutine coll_postInput()

  use crcoall
  use coll_db

  ! Call one extra time as some arrays depend on input values
  call collimation_expand_arrays(npart,nblz)

  if(c_enom <= zero) then
    write(lerr,"(a)") "COLL> ERROR Beam energy must bea larger than zero"
    call prror
  end if

  if(cdb_fileName == " ") then
    write(lerr,"(a)") "COLL> ERROR No collimator database file specified"
    call prror
  end if

end subroutine coll_postInput

! ================================================================================================ !
!  Open the main output files for collimation
! ================================================================================================ !
subroutine coll_openFiles

  use mod_units
  use string_tools
  use mod_common, only : numl
  use mod_settings
  use coll_common
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif

#ifdef HDF5
  type(h5_dataField), allocatable :: setFields(:)
  integer fmtHdf
#endif

!open files, in CR mode perform additional checks

  ! Survival File
#ifdef CR
  if(coll_survivalFilePos == -1) then
#endif
    call f_requestUnit(coll_survivalFile,coll_survivalUnit)
    call f_open(unit=coll_survivalUnit,file=coll_survivalFile,formatted=.true.,mode="w",status="replace")
    write(coll_survivalUnit,"(a7,1x,a9)") "#  turn","n_part"
    flush(coll_survivalUnit)
#ifdef CR
    coll_survivalFilePos = 1
  end if
#endif

  ! Collimator Gaps File
#ifdef CR
  if(coll_gapsFilePos == -1) then
#endif
    call f_requestUnit(coll_gapsFile,coll_gapsUnit)
    call f_open(unit=coll_gapsUnit,file=coll_gapsFile,formatted=.true.,mode="w",status="replace")
    write(coll_gapsUnit,"(a1,1x,a2,1x,a16,4(1x,a19),1x,a4,5(1x,a13),1x,a13)")     &
      "#","ID","name            ","angle[rad]","betax[m]","betay[m]","halfgap[m]",&
      "mat.","length[m]","sigx[m]","sigy[m]","tilt1[rad]","tilt2[rad]","nsig"
    flush(coll_gapsUnit)
#ifdef CR
    coll_gapsFilePos = 1
  end if
#endif

  ! Collimator Settings (Jaw Slices)
#ifdef CR
  if(coll_settingsFilePos == -1) then
#endif
    call f_requestUnit(coll_settingsFile,coll_settingsUnit)
    call f_open(unit=coll_settingsUnit,file=coll_settingsFile,formatted=.true.,mode="w",status="replace")
    write(coll_settingsUnit,"(a20,1x,a10,5(1x,a13),1x,a4)") chr_rPad("# name",20),"slice","halfgap[m]","gapoffset[m]",&
      "tiltjaw1[rad]","tiltjaw2[rad]","length[m]","mat."
    flush(coll_settingsUnit)
#ifdef CR
    coll_settingsFilePos = 1
  end if
#endif

  ! Positions
#ifdef CR
  if(coll_positionsFilePos == -1) then
#endif
    call f_requestUnit(coll_positionsFile,coll_positionsUnit)
    call f_open(unit=coll_positionsUnit,file=coll_positionsFile,formatted=.true.,mode="w",status="replace")
    write(coll_positionsUnit,"(a)") "# Ind           Name   Pos[m]"
    flush(coll_positionsUnit)
#ifdef CR
    coll_positionsFilePos = 1
  end if
#endif

  ! Tracks Files
  if(dowrite_tracks) then
#ifdef CR
    if(coll_tracksFilePos == -1) then
#endif
      call f_requestUnit(coll_tracksFile,coll_tracksUnit)
      call f_open(unit=coll_tracksUnit,file=coll_tracksFile,formatted=.true.,mode="w",status="replace")
      write(coll_tracksUnit,"(a)") "# name turn s x xp y yp DE/E type"
      flush(coll_tracksUnit)
#ifdef CR
      coll_tracksFilePos = 1
    end if
#endif
  end if

#ifdef CR
  if(coll_pencilFilePos == -1) then
#endif
    call f_requestUnit(coll_pencilFile,coll_pencilUnit)
    call f_open(unit=coll_pencilUnit, file=coll_pencilFile,formatted=.true.,mode="w",status="replace")
    write(coll_pencilUnit,"(a)") "# x xp y yp"
    flush(coll_pencilUnit)
#ifdef CR
    coll_pencilFilePos = 1
  end if
#endif

  ! Crystal Files
  if(coll_hasCrystal .and. dowrite_crycoord) then
    if(st_debug) then
#ifdef CR
      if(coll_cryEntFilePos == -1) then
#endif
        call f_requestUnit(coll_cryEntFile,coll_cryEntUnit)
        call f_open(unit=coll_cryEntUnit,file=coll_cryEntFile,formatted=.true.,mode="w",status="replace")
        write(coll_cryEntUnit,"(a1,1x,a6,1x,a8,1x,a20,1x,a4,2(1x,a3),5(1x,a15))") &
          "#","partID","turn",chr_rPad("collimator",20),"mat.","hit","abs","x","xp","y","yp","p"
        flush(coll_cryEntUnit)
#ifdef CR
        coll_cryEntFilePos = 1
      end if
#endif

#ifdef CR
      if(coll_cryExitFilePos == -1) then
#endif
        call f_requestUnit(coll_cryExitFile,coll_cryExitUnit)
        call f_open(unit=coll_cryExitUnit,file=coll_cryExitFile,formatted=.true.,mode="w",status="replace")
        write(coll_cryExitUnit,"(a1,1x,a6,1x,a8,1x,a20,1x,a4,2(1x,a3),5(1x,a15))") &
          "#","partID","turn",chr_rPad("collimator",20),"mat.","hit","abs","x","xp","y","yp","p"
        flush(coll_cryExitUnit)
#ifdef CR
        coll_cryExitFilePos = 1
      end if
#endif
    end if

#ifdef CR
    if(coll_cryInterFilePos == -1) then
#endif
      call f_requestUnit(coll_cryInterFile,coll_cryInterUnit)
      call f_open(unit=coll_cryInterUnit,file=coll_cryInterFile,formatted=.true.,mode="w",status="replace")
      write(coll_cryInterUnit,"(a1,1x,a6,1x,a8,1x,a20,1x,a4,1x,a4,10(1x,a15))") &
          "#","partID","turn",chr_rPad("collimator",20),"prev","proc","kickx","kicky","Ein","Eout", &
          "xpin","ypin","cryangle","xin","yin"
      flush(coll_cryInterUnit)
#ifdef CR
      coll_cryInterFilePos = 1
    end if
#endif
  end if

  if(do_select) then
#ifdef CR
    if(coll_ellipseFilePos == -1) then
#endif
      call f_requestUnit(coll_ellipseFile,coll_ellipseUnit)
      call f_open(unit=coll_ellipseUnit,file=coll_ellipseFile,formatted=.true.,mode="w",status="replace")
      write(coll_ellipseUnit,"(a)") "# name x y xp yp E s turn halo nabs_type"
      flush(coll_ellipseUnit)
#ifdef CR
      coll_ellipseFilePos = 1
    end if
#endif
  end if

  if(dowrite_impact) then

#ifdef CR
    if(coll_allImpactFilePos == -1) then
#endif
      call f_requestUnit(coll_allImpactFile, coll_allImpactUnit)
      call f_open(unit=coll_allImpactUnit, file=coll_allImpactFile, formatted=.true.,mode="w",status="replace")
      write(coll_allImpactUnit,"(a)") "# 1=name 2=turn 3=s"
      flush(coll_allImpactUnit)
#ifdef CR
      coll_allImpactFilePos = 1
    end if
#endif

#ifdef CR
    if(coll_allAbsorbFilePos == -1) then
#endif
      call f_requestUnit(coll_allAbsorbFile, coll_allAbsorbUnit)
      call f_open(unit=coll_allAbsorbUnit, file=coll_allAbsorbFile, formatted=.true.,mode="w",status="replace")
      write(coll_allAbsorbUnit,"(a)") "# 1=name 2=turn 3=s"
      flush(coll_allAbsorbUnit)
#ifdef CR
      coll_allAbsorbFilePos = 1
    end if
#endif

#ifdef CR
    if(coll_fstImpactFilePos == -1) then
#endif
      call f_requestUnit(coll_fstImpactFile, coll_fstImpactUnit)
      call f_open(unit=coll_fstImpactUnit, file=coll_fstImpactFile, formatted=.true.,mode="w",status="replace")
      write(coll_fstImpactUnit,"(a)") "# 1=name, 2=iturn, 3=icoll, 4=nabs, 5=s_imp[m], 6=s_out[m], "//&
        "7=x_in(b!)[m], 8=xp_in, 9=y_in, 10=yp_in, 11=x_out [m], 12=xp_out, 13=y_out, 14=yp_out"
      flush(coll_fstImpactUnit)
#ifdef CR
      coll_fstImpactFilePos = 1
    end if
#endif

#ifdef CR
    if(coll_scatterFilePos == -1) then
#endif
      call f_requestUnit(coll_scatterFile,   coll_scatterUnit)
      call f_open(unit=coll_scatterUnit,   file=coll_scatterFile,   formatted=.true.,mode="w",status="replace")
      write(coll_scatterUnit,"(a)") "# 1=icoll, 2=iturn, 3=np, 4=nabs (1:Nuclear-Inelastic,2:Nuclear-Elastic,3:pp-Elastic, "//&
        "4:Single-Diffractive,5:Coulomb), 5=dp, 6=dx', 7=dy'"
      flush(coll_scatterUnit)
#ifdef CR
      coll_scatterFilePos = 1
    end if
#endif

#ifdef CR
    if(coll_flukImpFilePos == -1) then
#endif
      call f_requestUnit(coll_flukImpFile,   coll_flukImpUnit)
      call f_open(unit=coll_flukImpUnit,   file=coll_flukImpFile,   formatted=.true.,mode="w",status="replace")
      write(coll_flukImpUnit,"(a)") "# 1=icoll 2=c_rotation 3=s 4=x 5=xp 6=y 7=yp 8=nabs 9=np 10=turn"
      flush(coll_flukImpUnit)
#ifdef CR
      coll_flukImpFilePos = 1
    end if
#endif

#ifdef CR
    if(coll_flukImpAllFilePos == -1) then
#endif
      call f_requestUnit(coll_flukImpAllFile,coll_flukImpAllUnit)
      call f_open(unit=coll_flukImpAllUnit,file=coll_flukImpAllFile,formatted=.true.,mode="w",status="replace")
      write(coll_flukImpAllUnit,"(a)") "# 1=icoll 2=c_rotation 3=s 4=x 5=xp 6=y 7=yp 8=nabs 9=np 10=turn"
      flush(coll_flukImpAllUnit)
#ifdef CR
      coll_flukImpAllFilePos = 1
    end if
#endif

#ifdef CR
    if(coll_jawProfileFilePos == -1) then
#endif
      call f_requestUnit(coll_jawProfileFile,coll_jawProfileUnit)
      call f_open(unit=coll_jawProfileUnit,file=coll_jawProfileFile,formatted=.true.,mode="w",status="replace")
      write(coll_jawProfileUnit,"(a1,1x,a6,1x,2(a7,1x),5(a17,1x),a12)") "#", "icoll", "iturn", "np", "x[m]", "xp[]", "y[m]", &
        "yp[]", "s[m]", "[1:in,2:out]"
      flush(coll_jawProfileUnit)
#ifdef CR
      coll_jawProfileFilePos = 1
    end if
#endif
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

  if(dowrite_tracks) then
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

end subroutine coll_openFiles

! ================================================================================================ !
!  Statistics for the amplitude.dat file
! ================================================================================================ !
subroutine coll_computeStats

  use mod_common
  use mod_common_main
  use mod_common_track
  use coll_common

  integer j
  real(kind=fPrec) gammax,gammay,xj,xpj,yj,ypj,pj,nspx,nspy

  if(dowrite_amplitude .eqv. .false.) then
    ! We're not writing the amplitude.dat file anyway, so return here.
    ! It shouldn't be called anyway, but just to make sure we don't segfault if it does.
    return
  end if

  gammax = (one + talphax(ie)**2)/tbetax(ie)
  gammay = (one + talphay(ie)**2)/tbetay(ie)

!  if(firstrun .and. iturn == 1) then
  if(iturn == 1) then
    sum_ax(ie) = zero
    sum_ay(ie) = zero
  end if

  do j=1,napx
    xj  = (xv1(j)-torbx(ie))/c1e3
    xpj = (yv1(j)-torbxp(ie))/c1e3
    yj  = (xv2(j)-torby(ie))/c1e3
    ypj = (yv2(j)-torbyp(ie))/c1e3
    pj  = ejv(j)/c1e3

    if(part_abs_pos(j) == 0 .and. part_abs_turn(j) == 0) then
      nspx         = sqrt(abs((gammax*(xj)**2 + ((two*talphax(ie))*xj)*xpj) + tbetax(ie)*xpj**2)/c_emitx0_collgap)
      nspy         = sqrt(abs((gammay*(yj)**2 + ((two*talphay(ie))*yj)*ypj) + tbetay(ie)*ypj**2)/c_emity0_collgap)
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

end subroutine coll_computeStats

! ================================================================================================ !
!  The routine that calls the actual scattering modules for a collimator
! ================================================================================================ !
subroutine coll_doCollimator(stracki)

  use crcoall
  use mod_time
  use mod_common
  use coll_common
  use mod_settings
  use mod_particles
  use mod_common_main
  use mod_common_track
  use coll_db
  use coll_k2
  use coll_jawfit
  use coll_dist
  use coll_crystal
  use mod_units
  use mathlib_bouncer
  use string_tools
#ifdef G4COLLIMATION
  use geant4
#endif

  real(kind=fPrec), intent(in) :: stracki

  integer j, iSlice, nSlices
  logical onesided, linside(napx), isAbs, isHit
  real(kind=fPrec) nsig,c_length,jawLength,jawAperture,jawOffset,jawTilt(2),x_Dump,xpDump,y_Dump,   &
    ypDump,s_Dump,xmax,ymax,calc_aperture,zpj,xmax_pencil,ymax_pencil,xmax_nom,ymax_nom,            &
    nom_aperture,scale_bx,scale_by,c_tilt(2),c_offset,c_aperture,c_rotation,cry_tilt, &
    cRot,sRot

  call time_startClock(time_clockCOLL)

  ! Get collimator values from DB
  icoll      = cdb_elemMap(c_ix)
  nsig       = cdb_cNSig(icoll)
  c_rotation = cdb_cRotation(icoll)
  c_length   = cdb_cLength(icoll)
  c_offset   = cdb_cOffset(icoll)
  c_tilt(:)  = cdb_cTilt(:,icoll)
  cRot       = cos_mb(c_rotation)
  sRot       = sin_mb(c_rotation)
  if(iturn == 1.and. dowrite_amplitude) then
    call coll_computeStats
  end if

  ! Get the aperture from the beta functions and emittance
  ! A simple estimate of beta beating can be included that has twice the betatron phase advance
  if(.not. cdb_doNSig) then
    nsig = cdb_cNSig(icoll)
  end if

  scale_bx = one + xbeat*sin_mb(xbeatphase)
  scale_by = one + ybeat*sin_mb(ybeatphase)

  if(firstcoll) then
    scale_bx0 = scale_bx
    scale_by0 = scale_by
    firstcoll = .false.
  end if

  ! Assign nominal OR design beta functions for later
  if(do_nominal) then
    bx_dist = (cdb_cBx(icoll)*scale_bx)/scale_bx0
    by_dist = (cdb_cBy(icoll)*scale_by)/scale_by0
  else
    bx_dist = (tbetax(ie)*scale_bx)/scale_bx0
    by_dist = (tbetay(ie)*scale_by)/scale_by0
  end if

  ! Write beam ellipse at selected collimator
  ! Checking lower case of collimator name against name_sel (which is already lower case)
  ! This is to ensure compatibility with old style COLL block which used upper case collimator name
  if(chr_toLower(cdb_cName(icoll)) == name_sel .and. do_select) then
    do j=1,napx
      write(coll_ellipseUnit,"(1x,i8,6(1x,e15.7),3(1x,i4,1x,i4))") partID(j),xv1(j),xv2(j),yv1(j),yv2(j), &
        ejv(j),sigmv(j),iturn,nhit_stage(j),nabs_type(j)
#ifdef CR
      flush(coll_ellipseUnit)
      coll_ellipseFilePos = coll_ellipseFilePos + 1
#endif
    end do
  end if

  if(iturn == 1) then
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
    flush(outlun)
#ifdef CR
    coll_trackoutPos = coll_trackoutPos + 12
#endif
  end if

  ! Calculate aperture of collimator
  nsig        = nsig + gap_rms_error(icoll)
  xmax        = nsig*sqrt(bx_dist*c_emitx0_collgap)
  ymax        = nsig*sqrt(by_dist*c_emity0_collgap)
  xmax_pencil = (nsig + pencil_offset)*sqrt(bx_dist*c_emitx0_collgap)
  ymax_pencil = (nsig + pencil_offset)*sqrt(by_dist*c_emity0_collgap)
  xmax_nom    = cdb_cNSig(icoll)*sqrt(cdb_cBx(icoll)*c_emitx0_collgap)
  ymax_nom    = cdb_cNSig(icoll)*sqrt(cdb_cBy(icoll)*c_emity0_collgap)

  cry_proc(:) = -1

  calc_aperture = sqrt(    xmax**2*cRot**2 +     ymax**2*sRot**2)
  nom_aperture  = sqrt(xmax_nom**2*cRot**2 + ymax_nom**2*sRot**2)

  ! Get x and y offsets at collimator center point
  x_pencil(icoll) = xmax_pencil * cRot
  y_pencil(icoll) = ymax_pencil * sRot

  ! Get corresponding beam angles (uses xp_max)
  xp_pencil(icoll) = (((-one*sqrt(c_emitx0_collgap/tbetax(ie)))*talphax(ie))*xmax)/sqrt(c_emitx0_collgap*tbetax(ie))
  yp_pencil(icoll) = (((-one*sqrt(c_emity0_collgap/tbetay(ie)))*talphay(ie))*ymax)/sqrt(c_emity0_collgap*tbetay(ie))
  xp_pencil0 = xp_pencil(icoll)
  yp_pencil0 = yp_pencil(icoll)

  pencil_dx(icoll) = sqrt(xmax_pencil**2 * cRot**2 + ymax_pencil**2 * sRot**2) - calc_aperture

  ! Added condition that pencil_distr /= 3 in order to do the tilt
  if(icoll == ipencil .and. iturn == 1 .and. pencil_distr /= 3) then
    ! To align the jaw/pencil to the beam always use the minus regardless which
    ! orientation of the jaws was used (symmetric/antisymmetric)
    c_tilt(1) = c_tilt(1) + (xp_pencil0 * cRot + sRot * yp_pencil0)
    c_tilt(2) = c_tilt(2) - (xp_pencil0 * cRot + sRot * yp_pencil0)
    write(lout,*) "INFO> Changed tilt1  ICOLL  to  ANGLE  ", icoll, c_tilt(1)
    write(lout,*) "INFO> Changed tilt2  ICOLL  to  ANGLE  ", icoll, c_tilt(2)
  end if

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
  if(iturn == 1) then
    if(iturn == 1) then
      write(outlun,*) xp_pencil(icoll), yp_pencil(icoll), pencil_dx(icoll)
      write(outlun,'(a,i4)') 'Collimator number:   ', icoll
      write(outlun,*) 'Beam size x [m]:     ', sqrt(tbetax(ie)*c_emitx0_collgap), "(from collgap emittance)"
      write(outlun,*) 'Beam size y [m]:     ', sqrt(tbetay(ie)*c_emity0_collgap), "(from collgap emittance)"
      write(outlun,*) 'Divergence x [urad]:     ', c1e6*xp_pencil(icoll)
      write(outlun,*) 'Divergence y [urad]:     ', c1e6*yp_pencil(icoll)
      write(outlun,*) 'Aperture (nom) [m]:  ', nom_aperture
      write(outlun,*) 'Aperture (cal) [m]:  ', calc_aperture
      write(outlun,*) 'Collimator halfgap [sigma]:  ', nsig
      write(outlun,*) 'RMS error on halfgap [sigma]:  ', gap_rms_error(icoll)
      write(outlun,*) ' '
      flush(outlun)
#ifdef CR
      coll_trackoutPos = coll_trackoutPos + 11
#endif

      write(coll_gapsUnit,"(i4,1x,a16,4(1x,e19.10),1x,a4,5(1x,e13.5),1x,f13.6)") &
        icoll,cdb_cName(icoll)(1:16),cdb_cRotation(icoll),tbetax(ie),tbetay(ie),calc_aperture, &
        cdb_cMaterial(icoll),cdb_cLength(icoll),sqrt(tbetax(ie)*c_emitx0_collgap), &
        sqrt(tbetay(ie)*c_emity0_collgap),cdb_cTilt(1,icoll),cdb_cTilt(2,icoll),nsig
#ifdef CR
      flush(coll_gapsUnit)
      coll_gapsFilePos = coll_gapsFilePos + 1
#endif

      ! Write to coll settings file if we have 0 or 1 slices
      if(nSlices <= 1) then
        write(coll_settingsUnit,"(a20,1x,i10,5(1x,1pe13.6),1x,a)") cdb_cName(icoll)(1:20), nSlices, calc_aperture, &
          cdb_cOffset(icoll), cdb_cTilt(1,icoll), cdb_cTilt(2,icoll), cdb_cLength(icoll), cdb_cMaterial(icoll)
#ifdef CR
        flush(coll_settingsUnit)
        coll_settingsFilePos = coll_settingsFilePos + 1
#endif
      end if
    end if
  end if !if(firstrun)

  c_aperture = two*calc_aperture

  ! Addition matched halo sampled directly on the TCP using pencil beam flag
  if(iturn == 1 .and. ipencil == icoll .and. pencil_distr == 3) then
    call coll_matchedHalo(c_tilt,c_offset,c_aperture,c_length)
    call part_updatePartEnergy(1,.true.)
    if(st_debug) then
      call part_writeState("pencilbeam_distr_type3.dat", .true., .false.)
    end if
  end if

  ! Copy particle data to 1-dim array and go back to meters
  do j=1,napx
    rcx(j)   = (xv1(j)-torbx(ie)) /c1e3
    rcxp(j)  = (yv1(j)-torbxp(ie))/c1e3
    rcy(j)   = (xv2(j)-torby(ie)) /c1e3
    rcyp(j)  = (yv2(j)-torbyp(ie))/c1e3
    rcp(j)   = ejv(j)/c1e3
    rcs(j)   = zero

    rcx0(j)  = rcx(j)
    rcxp0(j) = rcxp(j)
    rcy0(j)  = rcy(j)
    rcyp0(j) = rcyp(j)
    rcp0(j)  = rcp(j)
    ejf0v(j) = ejfv(j)

    ! For zero length element track back half collimator length
    ! DRIFT PART
    if(stracki == zero) then
      if(iexact) then
        zpj    = sqrt(one - rcxp(j)**2 - rcyp(j)**2)
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

  ! Allow treatment of collimators as one-sided
  onesided = cdb_cSides(icoll) == 1 .or. cdb_cSides(icoll) == 2

  linside(:) = .false.

  if(cdb_cSliced(icoll) > 0) then ! Treatment of sliced collimators

    ! Now, loop over the number of slices and call k2coll_collimate each time.
    ! For each slice, the corresponding offset and angle are to be used.
    if(cdb_isCrystal(icoll)) then
      write(lerr,"(a)") "COLL> ERROR A crystal collimator cannot be sliced"
      call prror
    end if
    do iSlice=1,nSlices
      jawAperture = c_aperture
      jawOffset   = c_offset
      jawTilt     = c_tilt
      call jaw_getFitSliceValues(cdb_cSliced(icoll), iSlice, jawLength, jawAperture, jawOffset, jawTilt)
      if(iturn == 1) then
        write(coll_settingsUnit,"(a20,1x,i10,5(1x,1pe13.6),1x,a)") cdb_cName(icoll)(1:20), iSlice, &
          jawAperture/two, jawOffset, jawTilt(1), jawTilt(2), jawLength, cdb_cMaterial(icoll)
#ifdef CR
        flush(coll_settingsUnit)
        coll_settingsFilePos = coll_settingsFilePos + 1
#endif
      end if
      call k2coll_collimate(icoll, iturn, ie, jawLength, c_rotation, jawAperture,            &
        jawOffset, jawTilt, rcx, rcxp, rcy, rcyp, rcp, rcs, c_enom*c1m3, part_hit_pos,       &
        part_hit_turn, part_abs_pos, part_abs_turn, part_impact, part_indiv, part_linteract, &
        onesided, nhit_stage, iSlice, nabs_type, linside)
    end do

  else ! Treatment of non-sliced collimators

#ifndef G4COLLIMATION
    if(cdb_isCrystal(icoll)) then
      call cry_startElement(icoll,ie,c_emitx0_dist,c_emity0_dist,cry_tilt,c_length)
    end if
    call k2coll_collimate(icoll, iturn, ie, c_length, c_rotation, c_aperture, c_offset, &
      c_tilt, rcx, rcxp, rcy, rcyp, rcp, rcs, c_enom*c1m3, part_hit_pos, part_hit_turn, &
      part_abs_pos, part_abs_turn, part_impact, part_indiv, part_linteract,             &
      onesided, nhit_stage, 1, nabs_type, linside)

    if(cdb_isCrystal(icoll) .and. dowrite_crycoord) then
      if(st_debug) then
        do j=1,napx
          isHit = part_hit_pos(j) == ie .and. part_hit_turn(j) == iturn
          isAbs = part_abs_pos(j) == ie .and. part_abs_turn(j) == iturn
          write(coll_cryEntUnit,"(i8,1x,i8,1x,a20,1x,a4,2(3x,l1),5(1x,1pe15.8))")    &
            partID(j),iturn,cdb_cName(icoll)(1:20),cdb_cMaterial(icoll),isHit,isAbs, &
            rcx0(j),rcxp0(j),rcy0(j),rcyp0(j),rcp0(j)
#ifdef CR
          flush(coll_cryEntUnit)
          coll_cryEntFilePos = coll_cryEntFilePos + 1
#endif
          write(coll_cryExitUnit,"(i8,1x,i8,1x,a20,1x,a4,2(3x,l1),5(1x,1pe15.8))")   &
            partID(j),iturn,cdb_cName(icoll)(1:20),cdb_cMaterial(icoll),isHit,isAbs, &
            rcx(j),rcxp(j),rcy(j),rcyp(j),rcp(j)
#ifdef CR
          flush(coll_cryExitUnit)
          coll_cryExitFilePos = coll_cryExitFilePos + 1
#endif
        end do
      end if
      do j=1,napx
        if(cry_proc(j) > 0) then
          write(coll_cryInterUnit,"(i8,1x,i8,1x,a20,1x,i4,1x,i4,10(1x,1pe15.8))")        &
            partID(j), iturn, cdb_cName(icoll)(1:20),cry_proc_prev(j),cry_proc(j),rcxp(j)-rcxp0(j), &
            rcyp(j)-rcyp0(j),rcp0(j),rcp(j),rcxp0(j),rcyp0(j),cry_tilt,rcx0(j),rcy0(j)
#ifdef CR
          flush(coll_cryInterUnit)
          coll_cryInterFilePos = coll_cryInterFilePos + 1
#endif
        end if
      end do
    end if
#else
    call coll_doCollimator_Geant4(c_aperture,c_rotation,c_length,onesided)
#endif

  end if

#ifndef G4COLLIMATION
  ! Calculate average impact parameter and save info for all collimators.
  ! Copy information back and do negative drift.

  if(dowrite_impact) then
    ! Update writeout of jaw profiles
    do j=1,napx
      if(linside(j) .and. sqrt(rcx(j)**2 + rcy(j)**2) < 99.0e-3_fPrec) then
        x_Dump =  rcx(j)*cRot + sRot*rcy(j)
        xpDump = rcxp(j)*cRot + sRot*rcyp(j)
        y_Dump =  rcy(j)*cRot - sRot*rcx(j)
        ypDump = rcyp(j)*cRot - sRot*rcxp(j)
        s_Dump = c_length
        write(coll_jawProfileUnit,"(3(1x,i7),5(1x,e17.9),1x,i1)") &
          icoll,iturn,partID(j),x_Dump,xpDump,y_Dump,ypDump,s_Dump,2
#ifdef CR
        flush(coll_jawProfileUnit)
        coll_jawProfileFilePos = coll_jawProfileFilePos + 1
#endif
      end if
    end do
  end if
#endif

#ifdef G4COLLIMATION
  do j=1,napx
    if(stracki == zero) then
      if(iexact) then
        zpj    = sqrt(one-rcxp(j)**2-rcyp(j)**2)
        rcx(j) = rcx(j) - (half*c_length)*(rcxp(j)/zpj)
        rcy(j) = rcy(j) - (half*c_length)*(rcyp(j)/zpj)
      else
        rcx(j) = rcx(j) - (half*c_length)*rcxp(j)
        rcy(j) = rcy(j) - (half*c_length)*rcyp(j)
      end if
    end if

    ! Copy data back to the original vector
    xv1(j) =  rcx(j)*c1e3 + torbx(ie)
    yv1(j) = rcxp(j)*c1e3 + torbxp(ie)
    xv2(j) =  rcy(j)*c1e3 + torby(ie)
    yv2(j) = rcyp(j)*c1e3 + torbyp(ie)
    ejv(j) =  rcp(j)*c1e3
  end do
  call part_updatePartEnergy(1,.true.)
#else
  ! Copy particle data back and do path length stuff; check for absorption
  ! Add orbit offset back
  do j=1,napx
    if(part_hit_pos(j) == ie .and. part_hit_turn(j) == iturn) then
      ! In order to get rid of numerical errors, just do the treatment for impacting particles
      if(stracki == zero) then
        ! For zero length element track back half collimator length
        if(iexact) then
          zpj    = sqrt(one-rcxp(j)**2-rcyp(j)**2)
          rcx(j) = rcx(j) - (half*c_length)*(rcxp(j)/zpj)
          rcy(j) = rcy(j) - (half*c_length)*(rcyp(j)/zpj)
        else
          rcx(j) = rcx(j) - (half*c_length)*rcxp(j)
          rcy(j) = rcy(j) - (half*c_length)*rcyp(j)
        end if
      end if

      ! Copy data back to the original vector
      xv1(j) =  rcx(j)*c1e3 + torbx(ie)
      yv1(j) = rcxp(j)*c1e3 + torbxp(ie)
      xv2(j) =  rcy(j)*c1e3 + torby(ie)
      yv2(j) = rcyp(j)*c1e3 + torbyp(ie)
      ejv(j) =  rcp(j)*c1e3
    end if
  end do

  call part_updatePartEnergy(1,.true.)

  ! The aperture check in this do loop should be reviewed and possibly removed
  do j=1,napx
    if(part_hit_pos(j) == ie .and. part_hit_turn(j) == iturn) then
      ! For absorbed particles set all coordinates to zero. Also include very
      ! large offsets, let's say above 100mm or 100mrad.
      if((part_abs_pos(j) /= 0 .and. part_abs_turn(j) /= 0) .or. &
        xv1(j) > c1e2 .or. yv1(j) > c1e2 .or. xv2(j) > c1e2 .or. yv2(j) > c1e2) then
        xv1(j)           = zero
        yv1(j)           = zero
        xv2(j)           = zero
        yv2(j)           = zero
        ejv(j)           = c_enom
        sigmv(j)         = zero
        part_abs_pos(j)  = ie
        part_abs_turn(j) = iturn
        pstop(j)         = .true.
        numxv(j)         = numx
        nabs_type(j)     = 0
        nhit_stage(j)    = 0
      end if
    end if
  end do

  if(dowrite_impact) then
    call coll_writeImpactAbsorb
  end if

  n_absorbed = 0
  do j=1,napx
    if(part_hit_pos(j) == ie .and. part_hit_turn(j) == iturn) then

      ! Calculate impact observables, fill histograms, save collimator info
      cn_impact(icoll) = cn_impact(icoll) + 1
      csum(icoll)      = csum(icoll)   + part_impact(j)
      csqsum(icoll)    = csqsum(icoll) + part_impact(j)**2

      if(part_abs_pos(j) /= 0 .and. part_abs_turn(j) /= 0) then
        ! If the interacting particle was lost, add-up counters for absorption
        ! Note: a particle with x/y >= 99. never hits anything any more in the logic of this program.
        ! Be careful to always fulfill this!
        n_absorbed         = n_absorbed + 1
        cn_absorbed(icoll) = cn_absorbed(icoll) + 1
        n_tot_absorbed     = n_tot_absorbed + 1

      elseif(part_abs_pos (j) == 0 .and. part_abs_turn(j) == 0) then
        nhit_stage(j) = ior(nhit_stage(j),cdb_cStage(icoll)) ! Record the hit type
      else
        write(lerr,"(a)")          "COLL> ERROR Particle cannot be both absorbed and not absorbed"
        write(lerr,"(a,2(1x,i0))") "COLL>      ",part_abs_pos (j),part_abs_turn(j)
        call prror
      end if

    end if
  end do

  if(dowrite_tracks) then
    call coll_writeTracks2(1)
  end if

  if(cn_impact(icoll) > 0) then
    caverage(icoll) = csum(icoll)/cn_impact(icoll)

    if(caverage(icoll)**2 > csqsum(icoll)/cn_impact(icoll)) then
      csigma(icoll) = 0
    else
      csigma(icoll) = sqrt(csqsum(icoll)/cn_impact(icoll) - caverage(icoll)**2)
    end if
  end if

  ! Checking lower case of collimator name against name_sel (which is already lower case)
  ! This is to ensure compatibility with old style COLL block which used upper case collimator name
  if(chr_toLower(cdb_cName(icoll)) == name_sel .and. iturn == 1) then
    call coll_writeSelectedCollimator
  end if
#endif

  call time_stopClock(time_clockCOLL)

end subroutine coll_doCollimator

! ================================================================================================ !
!  Collimate Exit
! ================================================================================================ !
subroutine coll_exitCollimation

  use crcoall
  use mod_units
  use mod_common
  use mod_common_track
  use string_tools
  use coll_common
  use coll_db
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif
#ifdef ROOT
  use root_output
#endif
#ifdef G4COLLIMATION
  use geant4
#endif

#ifdef HDF5
  type(h5_dataField), allocatable :: fldHdf(:)
  integer fmtHdf, setHdf
#endif
  integer i,k,ix

  ! Write the number of absorbed particles
  write(outlun,*) "INFO>  Number of impacts             : ", n_tot_absorbed+nsurvive_end
  write(outlun,*) "INFO>  Number of impacts at selected : ", num_selhit
  write(outlun,*) "INFO>  Number of surviving particles : ", nsurvive_end
  write(outlun,*) "INFO>  Number of absorbed particles  : ", n_tot_absorbed
  write(outlun,*)
  flush(outlun)
#ifdef CR
  coll_trackoutPos = coll_trackoutPos + 5
#endif

  if(n_tot_absorbed /= 0) then
    if(dowrite_efficiency) then
      write(outlun,*) " INFO>  Eff_r @  8 sigma    [e-4] : ", (neff(5)/real(n_tot_absorbed,fPrec))/c1m4
      write(outlun,*) " INFO>  Eff_r @ 10 sigma    [e-4] : ", (neff(9)/real(n_tot_absorbed,fPrec))/c1m4
      write(outlun,*) " INFO>  Eff_r @ 10-20 sigma [e-4] : ", ((neff(9)-neff(19))/(real(n_tot_absorbed,fPrec)))/c1m4
      write(outlun,*)
      write(outlun,*) neff(5)/real(n_tot_absorbed,fPrec), neff(9)/real(n_tot_absorbed,fPrec), &
        (neff(9)-neff(19))/(real(n_tot_absorbed,fPrec)), " !eff"
      write(outlun,*)
      flush(outlun)
#ifdef CR
      coll_trackoutPos = coll_trackoutPos + 6
#endif
    else
      write(outlun,*) "INFO> Efficiency calculations not enabled"
      flush(outlun)
#ifdef CR
      coll_trackoutPos = coll_trackoutPos + 1
#endif
    end if
  else
    write(lout,"(a)") "COLL> No particles absorbed"
  end if

  write(lout,"(a)")
  write(lout,"(a,i8)") "COLL> Number of impacts             : ", n_tot_absorbed+nsurvive_end
  write(lout,"(a,i8)") "COLL> Number of impacts at selected : ", num_selhit
  write(lout,"(a,i8)") "COLL> Number of surviving particles : ", nsurvive_end
  write(lout,"(a,i8)") "COLL> Number of absorbed particles  : ", n_tot_absorbed
  write(lout,"(a)")

  if(n_tot_absorbed /= 0) then
    if(dowrite_efficiency) then
      write(lout,"(a,f20.12)") "COLL> Eff_r @  8 sigma    [e-4] : ", (neff(5)/real(n_tot_absorbed,fPrec))/c1m4
      write(lout,"(a,f20.12)") "COLL> Eff_r @ 10 sigma    [e-4] : ", (neff(9)/real(n_tot_absorbed,fPrec))/c1m4
      write(lout,"(a,f20.12)") "COLL> Eff_r @ 10-20 sigma [e-4] : ", ((neff(9)-neff(19))/real(n_tot_absorbed,fPrec))/c1m4
    else
      write(lout,"(a)") "COLL> Efficiency calculations not enabled"
    end if
  else
    write(lout,"(a)") "COLL> No particle absorbed"
  end if
  write(lout,"(a)")

  if(dowrite_efficiency) then
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
      call f_open(unit=coll_efficUnit,file=coll_efficFile,formatted=.true.,mode="w",status="replace")
      if(n_tot_absorbed /= 0) then
        write(coll_efficUnit,"(a1,1x,a13,6(1x,a15),1x,a8)") "#","rad_sigma",&
          "frac_x","frac_y","frac_r","eff_x","eff_y","eff_r","n_abs"
        flush(coll_efficUnit)
#ifdef CR
        coll_efficFilePos = 1
#endif
        do k=1,numeff
          write(coll_efficUnit,"(7(e15.7,1x),i8)") rsig(k), neffx(k)/real(n_tot_absorbed,fPrec), &
            neffy(k)/real(n_tot_absorbed,fPrec), neff(k)/real(n_tot_absorbed,fPrec), neffx(k), neffy(k), neff(k), n_tot_absorbed
          flush(coll_efficUnit)
#ifdef CR
          coll_efficFilePos = coll_efficFilePos + 1
#endif
        end do
      else
        write(coll_efficUnit,"(a)") "No particles absorbed"
        flush(coll_efficUnit)
#ifdef CR
        coll_efficFilePos = coll_efficFilePos + 1
#endif
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
      call f_open(unit=coll_efficDPUnit,file=coll_efficDPFile,formatted=.true.,mode="w",status="replace")
      if(n_tot_absorbed /= 0) then
        write(coll_efficDPUnit,"(a1,1x,a13,2(1x,a15),2(1x,a8))") "#","dp/p","n_dpop/tot_nabs","n_dpop","tot_nabs","npart"
        flush(coll_efficDPUnit)
#ifdef CR
        coll_efficDPFilePos = 1
#endif
        do k=1,numeffdpop
          write(coll_efficDPUnit,"(e15.7,2(1x,e15.7),2(1x,i8))") dpopbins(k), neffdpop(k)/real(n_tot_absorbed,fPrec), neffdpop(k), &
            n_tot_absorbed, npartdpop(k)
          flush(coll_efficDPUnit)
#ifdef CR
          coll_efficDPFilePos = coll_efficDPFilePos + 1
#endif
        end do
      else
        write(coll_efficDPUnit,"(a)") "No particles absorbed"
        flush(coll_efficDPUnit)
#ifdef CR
        coll_efficDPFilePos = coll_efficDPFilePos + 1
#endif
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
      call f_open(unit=coll_effic2DUnit,file=coll_effic2DFile,formatted=.true.,mode="w",status="replace")
      if(n_tot_absorbed /= 0) then
        write(coll_effic2DUnit,"(a1,1x,a13,3(1x,a15),1x,a8)") "#","rad_sigma","dp/p","n/tot_nabs","n","tot_nabs"
        flush(coll_effic2DUnit)
#ifdef CR
        coll_effic2DFilePos = 1
#endif
        do i=1,numeff
          do k=1,numeffdpop
            write(coll_effic2DUnit,"(e15.7,3(1x,e15.7),1x,i8)") rsig(i), dpopbins(k), &
              neff2d(i,k)/real(n_tot_absorbed,fPrec), neff2d(i,k), n_tot_absorbed
            flush(coll_effic2DUnit)
#ifdef CR
            coll_effic2DFilePos = coll_effic2DFilePos + 1
#endif
          end do
        end do
      else
        write(coll_effic2DUnit,"(a)") "No particles absorbed"
        flush(coll_effic2DUnit)
#ifdef CR
        coll_effic2DFilePos = coll_effic2DFilePos + 1
#endif
      end if
      call f_close(coll_effic2DUnit)
#ifdef HDF5
    end if
#endif
  end if

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
    call f_open(unit=coll_summaryUnit,file=coll_summaryFile,formatted=.true.,mode="w",status="replace")
    write(coll_summaryUnit,"(a1,1x,a5,1x,a20,2(1x,a8),2(1x,a15),1x,a6)") "#","icoll",chr_rPad("collname",20),&
      "nimp","nabs","imp_av","imp_sig","length"
      flush(coll_summaryUnit)
#ifdef CR
      coll_summaryFilePos = 1
#endif
    do icoll = 1, cdb_nColl
      if(cdb_cLength(icoll) > zero .and. cdb_cFound(icoll)) then
        write(coll_summaryUnit,"(i7,1x,a20,2(1x,i8),2(1x,e15.7),1x,f6.3)") icoll, cdb_cName(icoll)(1:20), cn_impact(icoll), &
          cn_absorbed(icoll), caverage(icoll), csigma(icoll), cdb_cLength(icoll)
        flush(coll_summaryUnit)
#ifdef CR
        coll_summaryFilePos = coll_summaryFilePos + 1
#endif
      end if
    end do
    call f_close(coll_summaryUnit)
#ifdef HDF5
  end if
#endif

#ifdef ROOT
  if(root_flag .and. root_Collimation == 1) then
    do icoll = 1, cdb_nColl
      if(cdb_cLength(icoll) > zero) then
        call CollimatorLossRootWrite(icoll, cdb_cName(icoll), len(cdb_cName(icoll)), cn_impact(icoll), cn_absorbed(icoll), &
          caverage(icoll), csigma(icoll), cdb_cLength(icoll))
      end if
    end do
  end if
#endif

  call f_close(outlun)
  call f_close(coll_gapsUnit)
  call f_close(coll_positionsUnit)
  call f_close(coll_survivalUnit)

  if(dowrite_tracks) then
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
    call f_open(unit=coll_ampUnit,file=coll_ampFile,formatted=.true.,mode="w",status="replace")
    write(coll_ampUnit,"(a1,1x,a6,1x,a20,17(1x,a20))") "#","ielem",chr_rPad("name",20),"s","AX_AV","AX_RMS","AY_AV","AY_RMS",&
      "alphax","alphay","betax","betay","orbitx","orbity","dispx","dispy","xbob","ybob","xpbob","ypbob"
    flush(coll_ampUnit)
#ifdef CR
    coll_ampFilePos = 1
#endif
    do i=1,iu
      if(ic(i) <= nblo) then
        ix = mtyp(ic(i),mel(ic(i)))
      else
        ix = ic(i)-nblo
      end if
      write(coll_ampUnit,"(i8,1x,a20,17(1x,1pe20.13))") i, bez(ix)(1:20), dcum(i),                       &
        sum_ax(i)/real(max(nampl(i),1),fPrec),                                                           &
        sqrt(abs((sqsum_ax(i)/real(max(nampl(i),1),fPrec))-(sum_ax(i)/real(max(nampl(i),1),fPrec))**2)), &
        sum_ay(i)/real(max(nampl(i),1),fPrec),                                                           &
        sqrt(abs((sqsum_ay(i)/real(max(nampl(i),1),fPrec))-(sum_ay(i)/real(max(nampl(i),1),fPrec))**2)), &
        talphax(i), talphay(i), tbetax(i), tbetay(i), torbx(i), torby(i), tdispx(i), tdispy(i),          &
        xbob(i), ybob(i), xpbob(i), ypbob(i)
      flush(coll_ampUnit)
#ifdef CR
      coll_ampFilePos = coll_ampFilePos + 1
#endif
    end do
    call f_close(coll_ampUnit)
  end if

  ! Write orbitchecking.dat
  call f_requestUnit(coll_orbitCheckFile,coll_orbitCheckUnit)
  call f_open(unit=coll_orbitCheckUnit,file=coll_orbitCheckFile,formatted=.true.,mode="w",status="replace")
  write(coll_orbitCheckUnit,"(a1,1x,a6,3(1x,a15))") "#","s","torbitx","torbity"
  flush(coll_orbitCheckUnit)
#ifdef CR
  coll_orbitCheckFilePos = 1
#endif
  do i=1,iu
    write(coll_orbitCheckUnit,"(i8,3(1x,1pe15.7))") i, dcum(i), torbx(i), torby(i)
    flush(coll_orbitCheckUnit)
#ifdef CR
    coll_orbitCheckFilePos = coll_orbitCheckFilePos + 1
#endif
  end do
  call f_close(coll_orbitCheckUnit)

#ifdef G4COLLIMATION
  call g4_terminate()
#endif

end subroutine coll_exitCollimation

! ================================================================================================ !
!  This routine is called at the start of each tracking turn
! ================================================================================================ !
subroutine coll_startTurn(n)
  integer, intent(in) :: n
  iturn = n
end subroutine coll_startTurn

! ================================================================================================ !
!  This routine is called at the start of every element
! ================================================================================================ !
subroutine coll_startElement(iStru, iSing)

  use mod_common
  use mod_common_main
  use coll_common

  integer, intent(in) :: iStru
  integer, intent(in) :: iSing
  integer j

  ie   = iStru
  c_ix = iSing

#ifndef G4COLLIMATION
  ! For absorbed particles set all coordinates to zero.
  ! Also include very large offsets, let's say above 100mm or 100mrad.
  do j=1,napx
    if((part_abs_pos(j) /= 0 .and. part_abs_turn(j) /= 0) &
      .or. xv1(j) > c1e2 .or. yv1(j) > c1e2 .or. xv2(j) > c1e2 .or. yv2(j) > c1e2) then
      xv1(j)           = zero
      yv1(j)           = zero
      xv2(j)           = zero
      yv2(j)           = zero
      ejv(j)           = c_enom
      sigmv(j)         = zero
      nabs_type(j)     = 0
      nhit_stage(j)    = 0
      part_abs_pos(j)  = ie
      part_abs_turn(j) = iturn
      pstop(j)         = .true.
      numxv(j)         = numx
    end if
  end do
#endif

  ! Save coordinates of particle 1 to check orbit for amplitude file
  if(iturn == 1.and. dowrite_amplitude) then
    xbob(ie)  = xv1(1)
    ybob(ie)  = xv2(1)
    xpbob(ie) = yv1(1)
    ypbob(ie) = yv2(1)
  end if

end subroutine coll_startElement

! ================================================================================================ !
!  This routine is called at the end of every element
! ================================================================================================ !
subroutine coll_endElement

  use coll_common
  use mod_common

  if(iturn == 1.and. dowrite_amplitude) then
    call coll_computeStats
  end if

  ! Note: Not in firstrun
  if(dowrite_tracks) then
    call coll_writeTracks2(2)
  end if

end subroutine coll_endElement

! ================================================================================================ !
!  This routine is called at the end of every turn
! ================================================================================================ !
subroutine coll_endTurn

  use mod_time
  use mod_units
  use mod_common
  use mod_common_main
  use mod_common_track
  use coll_common
#ifdef ROOT
  use root_output
#endif
#ifdef HDF5
  use hdf5_output
#endif

  integer j,fUnit
#ifdef HDF5
  ! For other output
  type(h5_dataField), allocatable :: fldHdf(:)
  integer fmtHdf, setHdf
#endif

  call time_startClock(time_clockCOLL)

  do j=1,napx
    if(dowrite_efficiency) then
      xineff(j)  = xv1(j) - torbx (ie)
      xpineff(j) = yv1(j) - torbxp(ie)
      yineff(j)  = xv2(j) - torby (ie)
      ypineff(j) = yv2(j) - torbyp(ie)
    end if
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
#ifdef CR
      flush(coll_survivalUnit)
      coll_survivalFilePos = coll_survivalFilePos + 1
#endif
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

  if(iturn == 1 .or. ie == iu) then
    call coll_doEfficiency
  end if

  ! For LAST ELEMENT in the ring compact the arrays by moving all
  ! lost particles to the end of the array.
  if(ie == iu) then
    do j=1,napx
      if(xv1(j) < 99.0_fPrec .and. xv2(j) < 99.0_fPrec) then
        llostp(j) = .false.
      else
        llostp(j) = .true.
      end if
    end do
#ifndef G4COLLIMATION
    ! Move the lost particles to the end of the arrays
    call shuffleLostParticles
#endif
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

  ! Call end element one extra time
  call coll_endElement

  call time_stopClock(time_clockCOLL)

end subroutine coll_endTurn

! ================================================================================================ !
!  R. Bruce
!
!  Create distribution where the normalized distance between jaw and beam is the smallest.
!  This is where particles will first impact:
!  Without imperfections, it is:
!   - at the face of the collimator for the case of beta'<0 (POSITIVE alpha - beam converging) and
!   - at the exit of the collimator for the case of beta'>0 (NEGATIVE alpha beam diverging)
!  With imperfections:
!   - include errors on gap, tilt and offset. We have to calculate the normalized distance to each
!     corner separately!
!
! Calculate optical parameters at start and end of collimator (half a collimator length upstream
! and downstream of present s-position). Assuming a purely vertical or horizontal halo - need to
! add more conditions for other cases!
!
! Using standard twiss transfer matrix for a drift at start of collimator
! ================================================================================================ !
subroutine coll_matchedHalo(c_tilt,c_offset,c_aperture,c_length)

  use crcoall
  use coll_dist
  use mod_common
  use mod_common_main
  use mod_common_track
  use mathlib_bouncer

  real(kind=fPrec), intent(in) :: c_tilt(2)
  real(kind=fPrec), intent(in) :: c_offset
  real(kind=fPrec), intent(in) :: c_aperture
  real(kind=fPrec), intent(in) :: c_length

  integer j
  real(kind=fPrec) Nap1pos,Nap2pos,Nap1neg,Nap2neg,tiltOffsPos1,tiltOffsPos2,tiltOffsNeg1,     &
    tiltOffsNeg2,beamsize1,beamsize2,minAmpl,ldrift,c_nex2,c_ney2,betax1,betax2,betay1,betay2, &
    alphax1,alphax2,alphay1,alphay2,c_alphax,c_alphay,c_betax,c_betay

  ! Assign the drift length over which the optics functions are propagated
  ldrift = -c_length/two
  betax1 = tbetax(ie) - two*ldrift*talphax(ie) + (ldrift**2 * (one+talphax(ie)**2))/tbetax(ie)
  betay1 = tbetay(ie) - two*ldrift*talphay(ie) + (ldrift**2 * (one+talphay(ie)**2))/tbetay(ie)

  alphax1 = talphax(ie) - (ldrift*(1+talphax(ie)**2))/tbetax(ie)
  alphay1 = talphay(ie) - (ldrift*(1+talphay(ie)**2))/tbetay(ie)

  ! At end of collimator:
  ldrift = c_length/two
  betax2 = tbetax(ie) - two*ldrift*talphax(ie) + (ldrift**2 * (one+talphax(ie)**2))/tbetax(ie)
  betay2 = tbetay(ie) - two*ldrift*talphay(ie) + (ldrift**2 * (one+talphay(ie)**2))/tbetay(ie)

  alphax2 = talphax(ie) - (ldrift*(1+talphax(ie)**2))/tbetax(ie)
  alphay2 = talphay(ie) - (ldrift*(1+talphay(ie)**2))/tbetay(ie)

  ! Calculate beam size at start and end of collimator. account for collimation plane
  if(cdist_ampX > zero .and. cdist_ampY == zero) then ! Horizontal halo
    beamsize1 = sqrt(betax1 * c_emitx0_collgap)
    beamsize2 = sqrt(betax2 * c_emitx0_collgap)
  else if(cdist_ampX == zero .and. cdist_ampY > zero) then ! Vertical halo
    beamsize1 = sqrt(betay1 * c_emity0_collgap)
    beamsize2 = sqrt(betay2 * c_emity0_collgap)
  else
    write(lerr,"(a)") "COLL> ERROR Attempting to use a halo not purely in the horizontal "//&
      "or vertical plane with pencil_dist=3 - abort."
    call prror
  end if

  ! Calculate offset from tilt of positive and negative jaws, at start and end
  ! Remember: tilt angle is defined such that one corner stays at nominal position, the other corner is more open

  ! Jaw in positive x (or y):
  if(c_tilt(1) >= zero) then
    tiltOffsPos1 = zero
    tiltOffsPos2 = abs(sin_mb(c_tilt(1))) * c_length
  else
    tiltOffsPos1 = abs(sin_mb(c_tilt(1))) * c_length
    tiltOffsPos2 = zero
  end if

  ! Jaw in negative x (or y):
  if(c_tilt(2) >= zero) then
    tiltOffsNeg1 = abs(sin_mb(c_tilt(2))) * c_length
    tiltOffsNeg2 = zero
  else
    tiltOffsNeg1 = zero
    tiltOffsNeg2 = abs(sin_mb(c_tilt(2))) * c_length
  end if

  ! Calculate half distance from jaws to beam center (in units of beam sigma) at the beginning of the collimator,
  ! Positive and neg jaws.
  Nap1pos = ((c_aperture/two + c_offset) + tiltOffsPos1)/beamsize1
  Nap2pos = ((c_aperture/two + c_offset) + tiltOffsPos2)/beamsize2
  Nap1neg = ((c_aperture/two - c_offset) + tiltOffsNeg1)/beamsize1
  Nap2neg = ((c_aperture/two - c_offset) + tiltOffsNeg2)/beamsize2

  ! Minimum normalized distance from jaw to beam center - this is the n_sigma at which the halo should be generated
  minAmpl = min(Nap1pos,Nap2pos,Nap1neg,Nap2neg)

  ! Assign amplitudes in x and y for the halo generation function
  if(cdist_ampX > zero .and. cdist_ampY == zero) then ! Horizontal halo
    c_nex2 = minAmpl
    c_ney2 = zero
  else if(cdist_ampX == zero .and. cdist_ampY > zero) then ! Vertical halo
    c_ney2 = minAmpl
    c_nex2 = zero
  end if ! Other cases taken care of above - in these cases, program has already stopped

  ! Assign optics parameters to use for the generation of the starting halo - at start or end of collimator
  if(minAmpl == Nap1pos .or. minAmpl == Nap1neg) then ! min normalized distance occurs at start of collimator
    ldrift = -c_length/two
    c_alphax = alphax1
    c_alphay = alphay1
    c_betax  = betax1
    c_betay  = betay1
  else ! Min normalized distance occurs at end of collimator
    ldrift = c_length/two
    c_alphax = alphax2
    c_alphay = alphay2
    c_betax  = betax2
    c_betay  = betay2
  end if

  ! create new pencil beam distribution with spread at start or end of collimator at the minAmpl
  ! note: if imperfections are active, equal amounts of particles are still generated on the two jaws.
  ! but it might be then that only one jaw is hit on the first turn, thus only by half of the particles
  ! the particle generated on the other side will then hit the same jaw several turns later, possibly smearing the impact parameter
  ! This could possibly be improved in the future.
  call cdist_makeDist_coll(c_alphax,c_alphay,c_betax,c_betay,c_nex2,c_ney2)

  do j=1,napx
    xv1(j) = c1e3*xv1(j) + torbx(ie)
    yv1(j) = c1e3*yv1(j) + torbxp(ie)
    xv2(j) = c1e3*xv2(j) + torby(ie)
    yv2(j) = c1e3*yv2(j) + torbyp(ie)

    ! as main routine will track particles back half a collimator length (to start of jaw),
    ! track them now forward (if generated at face) or backward (if generated at end)
    ! 1/2 collimator length to center of collimator (ldrift pos or neg)
    xv1(j) = xv1(j) - ldrift*yv1(j)
    xv2(j) = xv2(j) - ldrift*yv2(j)
  end do

end subroutine coll_matchedHalo

! ================================================================================================ !
!  Find the smallest gap, and also write sigmasettings.out
!  Updated: 2019-10-10
! ================================================================================================ !
subroutine coll_getMinGapID(minGapID)

  use coll_db
  use coll_common
  use mod_units
  use string_tools
  use mathlib_bouncer
  use mod_common_track

  integer, intent(out) :: minGapID

  integer i,j
  real(kind=fPrec) gapH1,gapH2,gapH3,gapH4,minGap,nSigErr,sigOffset

#ifdef CR
  if(coll_sigmaSetFilePos == -1) then
#endif
    call f_requestUnit(coll_sigmaSetFile,coll_sigmaSetUnit)
    call f_open(unit=coll_sigmaSetUnit,file=coll_sigmaSetFile,formatted=.true.,mode="w",status="replace")
    write(coll_sigmaSetUnit,"(a1,1x,a18,12(1x,a16))") "#",chr_rPad("collimator",18),&
      "gap_h1","gap_h2","gap_h3","gap_h4","sig_offset","coll_offset","nsig",        &
      "gap_rms_error","beta_x","beta_y","orb_x","orb_y"
    flush(coll_sigmaSetUnit)
#ifdef CR
    coll_sigmaSetFilePos = 1
#endif

  minGap = 20.0_fPrec
  do i=1,cdb_nColl
    ! Start searching minimum gap
    if(cdb_cFound(i) .and. cdb_cLength(i) > zero) then
      nSigErr = cdb_cNSig(i) + gap_rms_error(i)

      ! Jaw 1 on positive side x-axis
      gapH1 = nSigErr - sin_mb(cdb_cTilt(1,i))*cdb_cLength(i)/2
      gapH2 = nSigErr + sin_mb(cdb_cTilt(1,i))*cdb_cLength(i)/2

      ! Jaw 2 on negative side of x-axis (see change of sign compared
      ! to above code lines, alos have a look to setting of tilt angle)
      gapH3 = nSigErr + sin_mb(cdb_cTilt(2,i))*cdb_cLength(i)/2
      gapH4 = nSigErr - sin_mb(cdb_cTilt(2,i))*cdb_cLength(i)/2

      j = cdb_struMap(i) ! The structure index of the collimator
      if(do_nominal) then
        bx_dist = cdb_cBx(i)
        by_dist = cdb_cBy(i)
      else
        bx_dist = tbetax(j)
        by_dist = tbetay(j)
      end if

      sigOffset = cdb_cOffset(i)/(sqrt(bx_dist**2 * cos_mb(cdb_cRotation(i))**2 + by_dist**2 * sin_mb(cdb_cRotation(i))**2))
      write(coll_sigmaSetUnit,"(a20,7(1x,f16.6),1x,1pe16.9,4(1x,f16.6))") cdb_cName(i),&
        gapH1,gapH2,gapH3,gapH4,sigOffset,cdb_cOffset(i),cdb_cNSig(i),gap_rms_error(i),&
        tbetax(j),tbetay(j),torbx(j),torby(j)
      flush(coll_sigmaSetUnit)
#ifdef CR
      coll_sigmaSetFilePos = coll_sigmaSetFilePos + 1
#endif

      if((gapH1 + sigOffset) <= minGap) then
        minGap   = gapH1 + sigOffset
        minGapID = i
      else if((gapH2 + sigOffset) <= minGap) then
        minGap   = gapH2 + sigOffset
        minGapID = i
      else if((gapH3 - sigOffset) <= minGap) then
        minGap   = gapH3 - sigOffset
        minGapID = i
      else if((gapH4 - sigOffset) <= minGap) then
        minGap   = gapH4 - sigOffset
        minGapID = i
      end if
    end if
  end do

  write(coll_sigmaSetUnit, "(a)")      ""
  write(coll_sigmaSetUnit, "(a)")      "# SUMMARY"
  write(coll_sigmaSetUnit, "(a)")      "# MinGap Collimator:  '"//trim(cdb_cName(minGapID))//"'"
  write(coll_sigmaSetUnit, "(a,i0)")   "# MinGap Coll ID:     ",minGapID
  write(coll_sigmaSetUnit, "(a,f0.6)") "# Min Gap Sigma:      ",minGap
  write(coll_sigmaSetUnit, "(a,i0)")   "# Pencil Initial:     ",ipencil
  flush(coll_sigmaSetUnit)
#ifdef CR
  coll_sigmaSetFilePos = coll_sigmaSetFilePos + 6
#endif

  if(ipencil > 0 .and. do_mingap) then
    write(coll_sigmaSetUnit, "(a,i0)") "# Pencil (do_mingap): ",minGapID
    flush(coll_sigmaSetUnit)
#ifdef CR
    coll_sigmaSetFilePos = coll_sigmaSetFilePos + 1
#endif
  end if
  write(coll_sigmaSetUnit, "(a,i0)") "# Seed before reinit: ",rnd_seed
  flush(coll_sigmaSetUnit)
#ifdef CR
  coll_sigmaSetFilePos = coll_sigmaSetFilePos + 1
  end if !if(coll_sigmaSetFilePos == -1) then
#endif

  call f_close(coll_sigmaSetUnit)

end subroutine coll_getMinGapID

! ================================================================================================ !
!  Do collimation analysis at last element.
!  If selecting, look at number of scattered particles at selected collimator.
! ================================================================================================ !
subroutine coll_doEfficiency

  use mod_common
  use mod_common_main
  use mod_common_track
  use mathlib_bouncer
  use coll_common

  integer j,ieff,ieffdpop
  real(kind=fPrec) xnorm,xpnorm,ynorm,ypnorm,xangle,yangle,driftx,drifty
  real(kind=fPrec) dpopmin,dpopmax,dnormx,dnormy,c_dpop
  real(kind=fPrec) sigX,sigY,sigX2,sigY2,cSigEffX,cSigEffY,cSigEff

  sigX2 = tbetax(ie)*c_emitx0_collgap
  sigY2 = tbetay(ie)*c_emity0_collgap
  sigX  = sqrt(sigX2)
  sigY  = sqrt(sigY2)

  do j=1,napx

    ! Do the binning in amplitude, only considering particles that were not absorbed before.
    if(part_abs_pos(j) == 0 .and. part_abs_turn(j) == 0 .and. part_select(j) == 1) then

      ! Populate the efficiency arrays at the end of each turn.
      if(ie == iu .and. dowrite_efficiency) then

        cSigEffX = ((xineff(j)*c1m3)**2 + (talphax(ie)*xineff(j)*c1m3 + tbetax(ie)*xpineff(j)*c1m3)**2)/sigX2
        cSigEffY = ((yineff(j)*c1m3)**2 + (talphay(ie)*yineff(j)*c1m3 + tbetay(ie)*ypineff(j)*c1m3)**2)/sigY2

        cSigEff  = sqrt(cSigEffX + cSigEffY)
        cSigEffX = sqrt(cSigEffX)
        cSigEffY = sqrt(cSigEffY)

        do ieff=1,numeff
          if(counted_r(j,ieff) == 0 .and. cSigEff >= rsig(ieff)) then
            neff(ieff)        = neff(ieff) + one
            counted_r(j,ieff) = 1
          end if

          if(counted_x(j,ieff) == 0 .and. cSigEffX >= rsig(ieff)) then
            neffx(ieff)       = neffx(ieff) + one
            counted_x(j,ieff) = 1
          end if

          if(counted_y(j,ieff) == 0 .and. cSigEffY >= rsig(ieff)) then
            neffy(ieff)       = neffy(ieff) + one
            counted_y(j,ieff) = 1
          end if

          ! 2D eff
          do ieffdpop=1,numeffdpop
            if(counted2d(j,ieff,ieffdpop) == 0 .and. abs((ejv(j)-c_enom)/c_enom) >= dpopbins(ieffdpop)) then
              neff2d(ieff,ieffdpop)      = neff2d(ieff,ieffdpop) + one
              counted2d(j,ieff,ieffdpop) = 1
            end if
          end do
        end do

        do ieffdpop=1,numeffdpop
          if(counteddpop(j,ieffdpop) == 0) then
            dpopmin = zero
            c_dpop  = abs((ejv(j)-c_enom)/c_enom)
            if(ieffdpop > 1) then
              dpopmin = dpopbins(ieffdpop-1)
            end if

            dpopmax = dpopbins(ieffdpop)
            if(c_dpop >= dpopmin .and. c_dpop < dpopmax) then
              npartdpop(ieffdpop) = npartdpop(ieffdpop) + 1
            end if
          end if

          if(counteddpop(j,ieffdpop) == 0 .and. abs((ejv(j)-c_enom)/c_enom) >= dpopbins(ieffdpop)) then
            neffdpop(ieffdpop)      = neffdpop(ieffdpop) + one
            counteddpop(j,ieffdpop) = 1
          end if
        end do
      end if

      ! Do an emittance drift
      driftx = driftsx*sigX
      drifty = driftsy*sigY

      if(ie == iu) then
        dnormx = driftx/sigX
        dnormy = drifty/sigY

        xnorm  = (xv1(j)*c1m3)/sigX
        xpnorm = (talphax(ie)*(xv1(j)*c1m3) + tbetax(ie)*(yv1(j)*c1m3))/sigX
        xangle = atan2_mb(xnorm,xpnorm)
        xnorm  = xnorm  + dnormx*sin_mb(xangle)
        xpnorm = xpnorm + dnormx*cos_mb(xangle)
        xv1(j) = c1e3*(xnorm*sigX)
        yv1(j) = c1e3*((xpnorm*sigX - (talphax(ie)*xv1(j))*c1m3)/tbetax(ie))

        ynorm  = (xv2(j)*c1m3)/sigY
        ypnorm = (talphay(ie)*(xv2(j)*c1m3) + tbetay(ie)*(yv2(j)*c1m3))/sigY
        yangle = atan2_mb(ynorm,ypnorm)
        ynorm  = ynorm  + dnormy*sin_mb(yangle)
        ypnorm = ypnorm + dnormy*cos_mb(yangle)
        xv2(j) = c1e3*(ynorm*sigY)
        yv2(j) = c1e3*((ypnorm*sigY - (talphay(ie)*xv2(j))*c1m3)/tbetay(ie))
      end if
    end if
  end do

end subroutine coll_doEfficiency

! ================================================================================================ !
!  WRITE TO FILES
! ================================================================================================ !

subroutine coll_writeImpactAbsorb

  use mod_common
  use mod_common_main
  use coll_common
#ifdef HDF5
  use hdf5_output
#endif

  integer j

  do j=1,napx
    ! First check for particle interaction at this collimator and this turn
    if(part_hit_pos(j) == ie .and. part_hit_turn(j) == iturn) then
      ! Fill the change in particle angle into histogram
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
#ifdef CR
        flush(coll_allImpactUnit)
        coll_allImpactFilePos = coll_allImpactFilePos + 1
#endif
#ifdef HDF5
      end if
#endif
      if(part_abs_pos(j) /= 0 .and. part_abs_turn(j) /= 0) then
#ifdef HDF5
        if(h5_useForCOLL) then
          call h5_prepareWrite(coll_hdf5_allAbsorb, 1)
          call h5_writeData(coll_hdf5_allAbsorb, 1, 1, partID(j))
          call h5_writeData(coll_hdf5_allAbsorb, 2, 1, iturn)
          call h5_writeData(coll_hdf5_allAbsorb, 3, 1, dcum(ie))
          call h5_finaliseWrite(coll_hdf5_allAbsorb)
        else
#endif
          write(coll_allAbsorbUnit,"(i8,1x,i8,1x,f10.3)") partID(j),iturn,dcum(ie)
#ifdef CR
        flush(coll_allAbsorbUnit)
        coll_allAbsorbFilePos = coll_allAbsorbFilePos + 1
#endif
#ifdef HDF5
        end if
#endif
      end if
    end if
  end do

end subroutine coll_writeImpactAbsorb

! ================================================================================================ !
!  For a SELECTED collimator only consider particles that were scattered on this selected
!  collimator at the first turn. All other particles are discarded.
!  - This is switched on with the DO_SELECT flag in the input file.
!  - Note that the part_select(j) flag defaults to 1 for all particles.
! ================================================================================================ !
subroutine coll_writeSelectedCollimator

  use crcoall
  use mod_units
  use mod_common
  use coll_common

  integer j, n_impact, num_selabs, num_surhit
  real(kind=fPrec) average, sigma, sum, sqsum

  n_impact   = 0
  num_selhit = 0
  num_surhit = 0
  num_selabs = 0
  sum        = zero
  sqsum      = zero

  if(dowrite_impact) then
#ifdef CR
    if(coll_impactFilePos == -1) then
#endif
      call f_requestUnit(coll_impactFile,coll_impactUnit)
      call f_open(unit=coll_impactUnit,file=coll_impactFile,formatted=.true.,mode="w",status="replace")

      write(coll_impactUnit,"(a1,1x,a14,1x,a16)") "#","impact","divergence"
      flush(coll_impactUnit)
#ifdef CR
      coll_impactFilePos = 1
    end if
#endif
  end if

  do j=1,napx
    if(part_hit_pos (j) == ie .and. part_hit_turn(j) == iturn) then
      num_selhit = num_selhit + 1
      if(part_abs_pos(j) == 0 .and. part_abs_turn(j) == 0) then
        num_surhit = num_surhit + 1
      else
        num_selabs = num_selabs + 1
      end if

      if(part_impact(j) < -half) then
        write(lerr,"(a,i0)") "COLL> ERROR Found invalid impact parameter ",part_impact(j)
        write(outlun,*)      "ERR>  Found invalid impact parameter ",part_impact(j)
        flush(outlun)
#ifdef CR
        coll_trackoutPos = coll_trackoutPos + 1
#endif
        call prror
      end if

      n_impact = n_impact + 1
      sum      = sum + part_impact(j)
      sqsum    = sqsum + part_impact(j)**2

      if(dowrite_impact) then
        write(coll_impactUnit,"(1pe16.9,1x,1pe16.9)") part_impact(j), part_indiv(j)
        flush(coll_impactUnit)
#ifdef CR
        coll_impactFilePos = coll_impactFilePos + 1
#endif
      end if

    else if(do_select .and. iturn == 1) then
      ! If we want to select only partciles interacting at the specified
      ! collimator then remove all other particles and reset the number
      ! of the absorbed particles to the selected collimator.
      part_select(j) = 0
      n_tot_absorbed = num_selabs
    end if
  end do

  if(n_impact > 0) then
    average = sum/n_impact
    if(sqsum/n_impact >= average**2) then
      sigma = sqrt(sqsum/n_impact - average**2)
    else
      sigma = zero
    end if
  end if

  write(lout,"(a,i8)")    "COLL> Selected collimator hits     : ", num_selhit
  write(lout,"(a,i8)")    "COLL> Number of impacts            : ", n_impact
  write(lout,"(a,i8)")    "COLL> Number of escaped protons    : ", num_surhit
  write(lout,"(a,f18.9)") "COLL> Average impact parameter [m] : ", average
  write(lout,"(a,f18.9)") "COLL> Sigma impact parameter [m]   : ", sigma

  if(dowrite_impact) then
    call f_close(coll_impactUnit)
  end if

end subroutine coll_writeSelectedCollimator

! ================================================================================================ !
!  Writing of tracks2.dat
! ================================================================================================ !
subroutine coll_writeTracks2(iMode)

  use coll_db
  use coll_common
  use mod_common
  use mod_common_main
  use mod_common_track
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif

  integer, intent(in) :: iMode

  integer j
  real(kind=fPrec) sigX2, sigY2, halfLen
  real(kind=fPrec) xj, xpj, yj, ypj, sj, pj
  real(kind=fPrec) xk, xpk, yk, ypk, sk

  sigX2 = tbetax(ie)*c_emitx0_collgap
  sigY2 = tbetay(ie)*c_emity0_collgap

  if(iMode == 1) then

    halfLen = half*cdb_cLength(icoll)

    sj = dcum(ie) - halfLen
    sk = dcum(ie) + halfLen

    do j=1,napx
      if(                                                                                 &
        part_hit_pos(j) == ie .and. part_hit_turn(j) == iturn .and.                       &
        part_abs_pos(j) == 0  .and. part_abs_turn(j) == 0     .and.                       &
        nhit_stage(j)    > 0  .and. xv1(j) < 99.0_fPrec .and. xv2(j) < 99.0_fPrec .and. ( &
          (xv1(j)*c1m3)**2 / sigX2 >= sigsecut2 .or.                                      &
          (xv2(j)*c1m3)**2 / sigY2 >= sigsecut2 .or.                                      &
          (xv1(j)*c1m3)**2 / sigX2 + (xv2(j)*c1m3)**2 / sigY2 >= sigsecut3                &
        )                                                                                 &
      ) then

        xpj = rcxp0(j)*c1e3 + torbxp(ie)
        xj  =  rcx0(j)*c1e3 +  torbx(ie) - halfLen*xpj
        ypj = rcyp0(j)*c1e3 + torbyp(ie)
        yj  =  rcy0(j)*c1e3 +  torby(ie) - halfLen*ypj
        xk  = xv1(j) + halfLen*yv1(j)
        xpk = yv1(j)
        yk  = xv2(j) + halfLen*yv2(j)
        ypk = yv2(j)
        pj  = (ejv(j) - c_enom)/c_enom

#ifdef HDF5
        if(h5_writeTracks2) then
          call h5tr2_writeLine(partID(j),iturn,sj,xj,xpj,yj,ypj,pj,nhit_stage(j))
          call h5tr2_writeLine(partID(j),iturn,sk,xk,xpk,yk,ypk,pj,nhit_stage(j))
        else
#endif
          write(coll_tracksUnit,"(i8,1x,i8,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)") &
            partID(j),iturn,sj,xj,xpj,yj,ypj,pj,nhit_stage(j)
          write(coll_tracksUnit,"(i8,1x,i8,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)") &
            partID(j),iturn,sk,xk,xpk,yk,ypk,pj,nhit_stage(j)
          flush(coll_tracksUnit)
#ifdef CR
          coll_tracksFilePos = coll_tracksFilePos + 2
#endif
#ifdef HDF5
        end if
#endif
      end if
    end do

  elseif(iMode == 2) then ! End of element record

    do j=1,napx
      if(                                                                                &
        part_abs_pos(j) == 0 .and. part_abs_turn(j) == 0 .and.                           &
        nhit_stage(j)    > 0 .and. xv1(j) < 99.0_fPrec .and. xv2(j) < 99.0_fPrec .and. ( &
          (xv1(j)*c1m3)**2 / sigX2 >= sigsecut2 .or.                                     &
          (xv2(j)*c1m3)**2 / sigY2 >= sigsecut2 .or.                                     &
          (xv1(j)*c1m3)**2 / sigX2 + (xv2(j)*c1m3)**2 / sigY2 >= sigsecut3               &
        )                                                                                &
      ) then

#ifdef HDF5
        if(h5_writeTracks2) then
          call h5tr2_writeLine(partID(j),iturn,dcum(ie),xv1(j),yv1(j),xv2(j),yv2(j),&
            (ejv(j)-c_enom)/c_enom,nhit_stage(j))
        else
#endif
          write(coll_tracksUnit,"(i8,1x,i8,1x,f10.2,4(1x,e12.5),1x,e11.3,1x,i4)") partID(j),iturn,&
            dcum(ie),xv1(j),yv1(j),xv2(j),yv2(j),(ejv(j)-c_enom)/c_enom,nhit_stage(j)
          flush(coll_tracksUnit)
#ifdef CR
          coll_tracksFilePos = coll_tracksFilePos + 1
#endif
#ifdef HDF5
        end if
#endif
      end if
    end do

  end if

end subroutine coll_writeTracks2

! ================================================================================================ !
!  Additional Routines under Compiler Flags
! ================================================================================================ !

#ifdef G4COLLIMATION
subroutine coll_doCollimator_Geant4(c_aperture,c_rotation,c_length,onesided)

  use, intrinsic :: iso_c_binding
  use crcoall
  use mod_common
  use mod_common_main
  use mod_common_track
  use coll_db
  use coll_common
  use geant4
  use string_tools
  use mathlib_bouncer
  use physical_constants
#ifdef ROOT
  use root_output
#endif

  real(kind=fPrec), intent(in) :: c_aperture
  real(kind=fPrec), intent(in) :: c_rotation
  real(kind=fPrec), intent(in) :: c_length
  logical,          intent(in) :: onesided

  integer j

  integer :: g4_lostc
  integer :: g4_npart
  integer :: part_hit_flag = 0
  integer :: part_abs_flag = 0
  real(kind=fPrec) x_tmp,y_tmp,xp_tmp,yp_tmp

  real(kind=fPrec) g4_time

  ! ien0,ien1: ion energy entering/leaving the collimator
  ! energy in MeV
  real(kind=fPrec)    :: ien0, ien1
  integer(kind=int16) :: nnuc0,nnuc1

  !! Add the geant4 geometry
  if(firstrun .and. iturn == 1) then
    call g4_add_collimator(cdb_cName(icoll), cdb_cMaterial(icoll), c_length, c_aperture, c_rotation, torbx(ie), torby(ie), &
         logical(onesided,kind=C_BOOL))
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

  call g4_set_maximum_particle_id(MaximumPartID)
  if(g4_debug .eqv. .true.) then
    write(lout,"(a,I11)") 'GEANT4>: Setting Maximum ParticleID: ', MaximumPartID
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
      flush(lout)
    end if

    x_tmp = rcx(j)
    y_tmp = rcy(j)
    xp_tmp = rcxp(j)
    yp_tmp = rcyp(j)
    rcx(j) =  x_tmp *cos_mb(c_rotation) + sin_mb(c_rotation)*y_tmp
    rcy(j) =  y_tmp *cos_mb(c_rotation) - sin_mb(c_rotation)*x_tmp
    rcxp(j) = xp_tmp*cos_mb(c_rotation) + sin_mb(c_rotation)*yp_tmp
    rcyp(j) = yp_tmp*cos_mb(c_rotation) - sin_mb(c_rotation)*xp_tmp

!! Translate particle position into time.
    g4_time = -(sigmv(j) * c1m3) / (beta0*clight)

!! Add all particles
    call g4_add_particle(rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j), pdgid(j), nzz(j), naa(j), nqq(j), nucm(j), &
      g4_time, partID(j), parentID(j), partWeight(j), spin_x(j), spin_y(j), spin_z(j))

! Log input energy + nucleons as per the FLUKA coupling
! rcp is in GeV
    nnuc0   = nnuc0 + naa(j)
    ien0    = ien0 + (rcp(j) * c1e3)
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

  call g4_get_maximum_particle_id(MaximumPartID)
  if(g4_debug .eqv. .true.) then
    write(lout,"(a,I11)") 'GEANT4>: Got Maximum ParticleID: ', MaximumPartID
    flush(lout)
  end if

  do j = 1, napx
!! Get the particle back + information
!! Remember C arrays start at 0, fortran at 1 here.
    call g4_collimate_return(j-1, rcx(j), rcy(j), rcxp(j), rcyp(j), rcp(j), pdgid(j), nucm(j), nzz(j), naa(j), nqq(j), &
      g4_time, partID(j), parentID(j), partWeight(j), &
      part_hit_flag, part_abs_flag, part_impact(j), part_indiv(j), part_linteract(j), spin_x(j), spin_y(j), spin_z(j))

    pstop (j) = .false.

    if(nqq(j) .eq. 0) then
      mtc (j) = zero
    else
      mtc (j) = (nqq(j)*nucm0)/(qq0*nucm(j))  ! hisix: mass to charge
    endif

!! Rotate back into the accelerator frame
    x_tmp   = rcx(j)
    y_tmp   = rcy(j)
    xp_tmp  = rcxp(j)
    yp_tmp  = rcyp(j)
    rcx(j)  = x_tmp *cos_mb(-one*c_rotation) + sin_mb(-one*c_rotation)*y_tmp
    rcy(j)  = y_tmp *cos_mb(-one*c_rotation) - sin_mb(-one*c_rotation)*x_tmp
    rcxp(j) = xp_tmp*cos_mb(-one*c_rotation) + sin_mb(-one*c_rotation)*yp_tmp
    rcyp(j) = yp_tmp*cos_mb(-one*c_rotation) - sin_mb(-one*c_rotation)*xp_tmp

!   L - (beta0 c t)
    sigmv(j) = (c_length - (beta0 * (clight * g4_time))) * c1e3

    part_impact(j) = 0
    part_indiv(j) = 0
    part_linteract(j) = 0

! Log output energy + nucleons as per the FLUKA coupling
    nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
    ien1        = ien1  + (rcp(j) * c1e3)                 ! outcoming energy

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
      flush(lout)
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
    write(unit208,"(2(i6,1x),e24.16)") icoll, (nnuc0-nnuc1), c1m3*(ien0-ien1)
    flush(unit208)
#ifdef CR
    fort208Pos = fort208Pos + 1
#endif
  end if

end subroutine coll_doCollimator_Geant4
#endif

! ================================================================================================ !
!  Echo the simulation settings to lout
! ================================================================================================ !
subroutine coll_echoSettings

  use parpro
  use crcoall
  use mod_common
  use mod_common_track
  use coll_db
  use coll_dist
  use coll_common
  use mod_units
  use string_tools

  integer i

#ifdef CR
  if(outlun == -1) then
#endif
  call f_requestUnit("colltrack.out", outlun)
  call f_open(unit=outlun,file="colltrack.out",formatted=.true.,mode="w",status="replace")
#ifdef CR
  end if
#endif

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
  flush(outlun)
#ifdef CR
  coll_trackoutPos = 19
#endif

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
  write(lout,"(a,e15.8)") 'COLL> Info: Emitx0_dist [um]    = ', emitnx0_dist*gammar
  write(lout,"(a,e15.8)") 'COLL> Info: Emity0_dist [um]    = ', emitny0_dist*gammar
  write(lout,"(a,e15.8)") 'COLL> Info: Emitx0_collgap [um] = ', emitnx0_collgap*gammar
  write(lout,"(a,e15.8)") 'COLL> Info: Emity0_collgap [um] = ', emitny0_collgap*gammar
  write(lout,"(a,e15.8)") 'COLL> Info: E0 [MeV]            = ', e0
  write(lout,"(a)")

  write(lout,"(a,i0)")    'COLL> Info: DIST_TYPES          = ', do_thisdis
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_NEX            = ', cdist_ampX
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_DEX            = ', cdist_smearX
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_NEY            = ', cdist_ampY
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_DEY            = ', cdist_smearY
  write(lout,"(a,a)")     'COLL> Info: DIST_FILE           = ', trim(cdist_fileName)
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_EN_ERROR       = ', cdist_spreadE
  write(lout,"(a,e15.8)") 'COLL> Info: DIST_BUNCHLENGTH    = ', cdist_bunchLen
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
  write(lout,"(a,l1)")    'COLL> Info: DOWRITE_EFFICIENCY  = ', dowrite_efficiency
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
  write(lout,"(a,l1)")    'COLL: Info: SYSTILT_ANTISYMM    = ', systilt_antisymm
  write(lout,"(a)")
  write(lout,"(a,i0)")    'COLL> Info: IPENCIL             = ', ipencil
  write(lout,"(a,e15.8)") 'COLL> Info: PENCIL_OFFSET       = ', pencil_offset
  write(lout,"(a,e15.8)") 'COLL> Info: PENCIL_RMSX         = ', pencil_rmsx
  write(lout,"(a,e15.8)") 'COLL> Info: PENCIL_RMSY         = ', pencil_rmsy
  write(lout,"(a,i0)")    'COLL> Info: PENCIL_DISTR        = ', pencil_distr
  write(lout,"(a)")
  write(lout,"(a,a)")     'COLL> Info: COLL_DB             = ', cdb_fileName
  write(lout,"(a)")
  write(lout,"(a,l1)")    'COLL> Info: dowrite_tracks       = ', dowrite_tracks
  write(lout,"(a)")
  write(lout,"(a,e15.8)") 'COLL> Info: SIGSECUT2           = ', sigsecut2
  write(lout,"(a,e15.8)") 'COLL> Info: SIGSECUT3           = ', sigsecut3
  write(lout,"(a)")
  write(lout,"(a,i0)")    'COLL> Info: NAPX                = ', napx
  write(lout,"(a,e15.8)") 'COLL> Info: Sigma_x0            = ', sqrt(tbetax(1)*c_emitx0_dist)
  write(lout,"(a,e15.8)") 'COLL> Info: Sigma_y0            = ', sqrt(tbetay(1)*c_emity0_dist)
  write(lout,"(a)")

  flush(lout)

end subroutine coll_echoSettings

! ================================================================================================================================ !
!  Begin Checkpoint Restart COPIED FROM APERTURE MODULE
! ================================================================================================================================ !
#ifdef CR

! ================================================================================================================================ !
subroutine coll_crcheck_readdata(fileUnit,readerr)

  use crcoall
  use coll_common
  use mod_alloc

  use coll_db, only : cdb_nColl_cr

  implicit none

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: readerr

  integer j,k,i

  call alloc(part_hit_pos_cr,   npart, 0,       "part_hit_pos_cr")
  call alloc(part_hit_turn_cr,  npart, 0,       "part_hit_turn_cr")
  call alloc(part_abs_pos_cr,   npart, 0,       "part_abs_pos_cr")
  call alloc(part_abs_turn_cr,  npart, 0,       "part_abs_turn_cr")
  call alloc(part_select_cr,    npart, 0,       "part_select_cr")
  call alloc(nabs_type_cr,      npart, 0,       "nabs_type_cr")
  call alloc(nhit_stage_cr,     npart, 0,       "nhit_stage_cr")
  call alloc(part_linteract_cr, npart, zero,    "part_linteract_cr")
  call alloc(part_indiv_cr,     npart, -c1m6,   "part_indiv_cr")
  call alloc(part_impact_cr,    npart, zero,    "part_impact_cr")

  call alloc(cry_proc_cr,       npart, -1,      "cry_proc_cr")
  call alloc(cry_proc_prev_cr,  npart, -1,      "cry_proc_prev_cr")
  call alloc(cry_proc_tmp_cr,   npart, -1,      "cry_proc_tmp_cr")

  call alloc(cn_impact_cr,     cdb_nColl_cr, 0,    "cn_impact")
  call alloc(cn_absorbed_cr,   cdb_nColl_cr, 0,    "cn_absorbed")
  call alloc(caverage_cr,      cdb_nColl_cr, zero, "caverage")
  call alloc(csigma_cr,        cdb_nColl_cr, zero, "csigma")
  call alloc(gap_rms_error_cr, cdb_nColl_cr, zero, "gap_rms_error")
  call alloc(csum_cr,          cdb_nColl_cr, zero, "csum")
  call alloc(csqsum_cr,        cdb_nColl_cr, zero, "csqsum")

  if(dowrite_amplitude) then
    call alloc(nampl_cr,    nblz, 0,    "nampl_cr")
    call alloc(sum_ax_cr,   nblz, zero, "sum_ax_cr")
    call alloc(sqsum_ax_cr, nblz, zero, "sqsum_ax_cr")
    call alloc(sum_ay_cr,   nblz, zero, "sum_ay_cr")
    call alloc(sqsum_ay_cr, nblz, zero, "sqsum_ay_cr")
    call alloc(xbob_cr,     nblz, zero, "xbob_cr")
    call alloc(ybob_cr,     nblz, zero, "ybob_cr")
    call alloc(xpbob_cr,    nblz, zero, "xpbob_cr")
    call alloc(ypbob_cr,    nblz, zero, "ypbob_cr")
  end if

  if(dowrite_efficiency) then
    call alloc(xineff_cr,      npart,                     zero, "xineff_cr")      ! (npart)
    call alloc(yineff_cr,      npart,                     zero, "yineff_cr")      ! (npart)
    call alloc(xpineff_cr,     npart,                     zero, "xpineff_cr")     ! (npart)
    call alloc(ypineff_cr,     npart,                     zero, "ypineff_cr")     ! (npart)
    call alloc(counted_r_cr,   npart, numeff,             0,    "counted_r_cr")   ! (npart,numeff)
    call alloc(counted_x_cr,   npart, numeff,             0,    "counted_x_cr")   ! (npart,numeff)
    call alloc(counted_y_cr,   npart, numeff,             0,    "counted_y_cr")   ! (npart,numeff)
    call alloc(counteddpop_cr, npart, numeffdpop,         0,    "counteddpop_cr") ! (npart,numeffdpop)
    call alloc(counted2d_cr,   npart, numeff, numeffdpop, 0,    "counted2d_cr")   ! (npart,numeff,numeffdpop)
    call alloc(npartdpop_cr,   numeffdpop,                0,    "npartdpop_cr")   ! (numeffdpop)
    call alloc(neff_cr,        numeff,                    zero, "neff_cr")        ! (numeff)
    call alloc(rsig_cr,        numeff,                    zero, "rsig_cr")        ! (numeff)
    call alloc(neffdpop_cr,    numeffdpop,                zero, "neffdpop_cr")    ! (numeffdpop)
    call alloc(dpopbins_cr,    numeffdpop,                zero, "dpopbins_cr")    ! (numeffdpop)
    call alloc(neff2d_cr,      numeff, numeffdpop,        zero, "neff2d_cr")      ! (numeff,numeffdpop)
    call alloc(neffx_cr,       numeff,                    zero, "neffx_cr")       ! (numeff)
    call alloc(neffy_cr,       numeff,                    zero, "neffy_cr")       ! (numeff)
  end if

  read(fileunit,err=100,end=100) &
  fort208Pos_CR,                 &
  coll_survivalFilePos_CR,       &
  coll_gapsFilePos_CR,           &
  coll_settingsFilePos_CR,       &
  coll_positionsFilePos_CR,      &
  coll_tracksFilePos_CR,         &
  coll_pencilFilePos_CR,         &
  coll_cryEntFilePos_CR,         &
  coll_cryExitFilePos_CR,        &
  coll_cryInterFilePos_CR,       &
  coll_ellipseFilePos_CR,        &
  coll_allImpactFilePos_CR,      &
  coll_allAbsorbFilePos_CR,      &
  coll_scatterFilePos_CR,        &
  coll_fstImpactFilePos_CR,      &
  coll_flukImpFilePos_CR,        &
  coll_flukImpAllFilePos_CR,     &
  coll_jawProfileFilePos_CR,     &
  coll_efficFilePos_CR,          &
  coll_efficDPFilePos_CR,        &
  coll_effic2DFilePos_CR,        &
  coll_summaryFilePos_CR,        &
  coll_ampFilePos_CR,            &
  coll_orbitCheckFilePos_CR,     &
  coll_sigmaSetFilePos_CR,       &
  coll_impactFilePos_CR,         &
  coll_trackoutPos_CR,           &
  lux_CR,                        &
  seed_CR,                       &
  k1_CR,                         &
  k2_CR,                         &
  n_tot_absorbed_cr,             &
  n_absorbed_cr,                 &
  nabs_total_cr,                 &
  nsurvive_cr,                   &
  nsurvive_end_cr,               &
  num_selhit_cr,                 &
  (part_hit_pos_cr(j)   ,j=1,npart), &
  (part_hit_turn_cr(j)  ,j=1,npart), &
  (part_abs_pos_cr(j)   ,j=1,npart), &
  (part_abs_turn_cr(j)  ,j=1,npart), &
  (part_select_cr(j)    ,j=1,npart), &
  (nabs_type_cr(j)      ,j=1,npart), &
  (nhit_stage_cr(j)     ,j=1,npart), &
  (part_linteract_cr(j) ,j=1,npart), &
  (part_indiv_cr(j)     ,j=1,npart), &
  (part_impact_cr(j)    ,j=1,npart), &

  (cry_proc_cr(j)       ,j=1,npart), &
  (cry_proc_prev_cr(j)  ,j=1,npart), &
  (cry_proc_tmp_cr(j)   ,j=1,npart), &

  (cn_impact_cr(j)      ,j=1,cdb_nColl_cr), &
  (cn_absorbed_cr(j)    ,j=1,cdb_nColl_cr), &
  (caverage_cr(j)       ,j=1,cdb_nColl_cr), &
  (csigma_cr(j)         ,j=1,cdb_nColl_cr), &
  (gap_rms_error_cr(j)  ,j=1,cdb_nColl_cr), &
  (csum_cr(j)           ,j=1,cdb_nColl_cr), &
  (csqsum_cr(j)         ,j=1,cdb_nColl_cr)

  if(dowrite_efficiency) then
    read(fileunit,err=100,end=100)                                  &
    (xineff_cr(j)          ,j=1,npart),                             &
    (yineff_cr(j)          ,j=1,npart),                             &
    (xpineff_cr(j)         ,j=1,npart),                             &
    (ypineff_cr(j)         ,j=1,npart),                             &
    ((counted_r_cr(i,j)    ,j=1,numeff),i=1,npart),                 & ! (npart,numeff)
    ((counted_x_cr(i,j)    ,j=1,numeff),i=1,npart),                 & ! (npart,numeff)
    ((counted_y_cr(i,j)    ,j=1,numeff),i=1,npart),                 & ! (npart,numeff)
    ((counteddpop_cr(i,j)  ,j=1,numeffdpop),i=1,npart),             & ! (npart,numeffdpop)
    (((counted2d_cr(i,j,k) ,k=1,numeffdpop),j=1,numeff),i=1,npart), & ! (npart,numeff,numeffdpop)
    (npartdpop_cr(j)       ,j=1,numeffdpop),                        & ! (numeffdpop)
    (neff_cr(j)            ,j=1,numeff),                            & ! (numeff)
    (rsig_cr(j)            ,j=1,numeff),                            & ! (numeff)
    (neffdpop_cr(j)        ,j=1,numeffdpop),                        & ! (numeffdpop)
    (dpopbins_cr(j)        ,j=1,numeffdpop),                        & ! (numeffdpop)
    ((neff2d_cr(i,j)       ,j=1,numeffdpop),i=1,numeff),            & ! (numeff,numeffdpop)
    (neffx_cr(j)           ,j=1,numeff),                            & ! (numeff)
    (neffy_cr(j)           ,j=1,numeff)                               ! (numeff)
  end if

  if(dowrite_amplitude) then
    read(fileunit,err=100,end=100)    &
    (nampl_cr(j)          ,j=1,nblz), &
    (sum_ax_cr(j)         ,j=1,nblz), &
    (sqsum_ax_cr(j)       ,j=1,nblz), &
    (sum_ay_cr(j)         ,j=1,nblz), &
    (sqsum_ay_cr(j)       ,j=1,nblz), &
    (xbob_cr(j)           ,j=1,nblz), &
    (ybob_cr(j)           ,j=1,nblz), &
    (xpbob_cr(j)          ,j=1,nblz), &
    (ypbob_cr(j)          ,j=1,nblz)
  end if
!end read()

  readerr = .false.
  return

100 continue
  readerr = .true.
  write(lout, "(a,i0,a)") "CR_CHECK> ERROR Reading C/R file fort.",fileUnit," in COLLIMATION"
  write(crlog,"(a,i0,a)") "CR_CHECK> ERROR Reading C/R file fort.",fileUnit," in COLLIMATION"
  flush(crlog)

end subroutine coll_crcheck_readdata

! ================================================================================================================================ !
subroutine coll_crcheck_positionFiles

  use parpro
  use crcoall
  use mod_units
  use coll_common
  use mod_common, only : fort208, unit208

  implicit none

  logical isOpen, fErr
  integer iError
  integer j
  character(len=mInputLn) aRecord

  call f_positionFile(fort208, unit208, fort208Pos, fort208Pos_CR)

  call f_positionFile(coll_survivalFile, coll_survivalUnit, coll_survivalFilePos, coll_survivalFilePos_CR)
  call f_positionFile(coll_gapsFile, coll_gapsUnit, coll_gapsFilePos, coll_gapsFilePos_CR)
  call f_positionFile(coll_settingsFile, coll_settingsUnit, coll_settingsFilePos, coll_settingsFilePos_CR)
  call f_positionFile(coll_positionsFile, coll_positionsUnit, coll_positionsFilePos, coll_positionsFilePos_CR)
  call f_positionFile(coll_tracksFile, coll_tracksUnit, coll_tracksFilePos, coll_tracksFilePos_CR)
  call f_positionFile(coll_pencilFile, coll_pencilUnit, coll_pencilFilePos, coll_pencilFilePos_CR)
  call f_positionFile(coll_cryEntFile, coll_cryEntUnit, coll_cryEntFilePos, coll_cryEntFilePos_CR)
  call f_positionFile(coll_cryExitFile, coll_cryExitUnit, coll_cryExitFilePos, coll_cryExitFilePos_CR)
  call f_positionFile(coll_cryInterFile, coll_cryInterUnit, coll_cryInterFilePos, coll_cryInterFilePos_CR)
  call f_positionFile(coll_ellipseFile, coll_ellipseUnit, coll_ellipseFilePos, coll_ellipseFilePos_CR)
  call f_positionFile(coll_allImpactFile, coll_allImpactUnit, coll_allImpactFilePos, coll_allImpactFilePos_CR)
  call f_positionFile(coll_allAbsorbFile, coll_allAbsorbUnit, coll_allAbsorbFilePos, coll_allAbsorbFilePos_CR)
  call f_positionFile(coll_scatterFile, coll_scatterUnit, coll_scatterFilePos, coll_scatterFilePos_CR)
  call f_positionFile(coll_fstImpactFile, coll_fstImpactUnit, coll_fstImpactFilePos, coll_fstImpactFilePos_CR)
  call f_positionFile(coll_flukImpFile, coll_flukImpUnit, coll_flukImpFilePos, coll_flukImpFilePos_CR)
  call f_positionFile(coll_flukImpAllFile, coll_flukImpAllUnit, coll_flukImpAllFilePos, coll_flukImpAllFilePos_CR)
  call f_positionFile(coll_jawProfileFile, coll_jawProfileUnit, coll_jawProfileFilePos, coll_jawProfileFilePos_CR)
  call f_positionFile(coll_sigmaSetFile, coll_sigmaSetUnit, coll_sigmaSetFilePos,coll_sigmaSetFilePos_CR)
  call f_positionFile(coll_impactFile, coll_impactUnit, coll_impactFilePos, coll_impactFilePos_CR)

  call f_positionFile("colltrack.out", outlun, coll_trackoutPos, coll_trackoutPos_CR)

!Files in the final processing
  call f_positionFile(coll_efficFile, coll_efficUnit, coll_efficFilePos, coll_efficFilePos_CR)
  call f_positionFile(coll_efficDPFile, coll_efficDPUnit, coll_efficDPFilePos, coll_efficDPFilePos_CR)
  call f_positionFile(coll_effic2DFile, coll_effic2DUnit, coll_effic2DFilePos, coll_effic2DFilePos_CR)
  call f_positionFile(coll_summaryFile, coll_summaryUnit, coll_summaryFilePos, coll_summaryFilePos_CR)
  call f_positionFile(coll_ampFile, coll_ampUnit, coll_ampFilePos, coll_ampFilePos_CR)
  call f_positionFile(coll_orbitCheckFile, coll_orbitCheckUnit, coll_orbitCheckFilePos, coll_orbitCheckFilePos_CR)

end subroutine coll_crcheck_positionFiles

! ================================================================================================================================ !
subroutine coll_crpoint(fileUnit,lerror)
  use crcoall
  use coll_common
  use mod_ranlux
  use parpro

  use coll_db, only : cdb_nColl
  use mod_units
  implicit none

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: lerror

  integer lux,seed,k1,k2
  integer i,j,k

  call rluxat(lux,seed,k1,k2)

  write(fileUnit,err=100)   &
    fort208Pos,             &
    coll_survivalFilePos,   &
    coll_gapsFilePos,       &
    coll_settingsFilePos,   &
    coll_positionsFilePos,  &
    coll_tracksFilePos,     &
    coll_pencilFilePos,     &
    coll_cryEntFilePos,     &
    coll_cryExitFilePos,    &
    coll_cryInterFilePos,   &
    coll_ellipseFilePos,    &
    coll_allImpactFilePos,  &
    coll_allAbsorbFilePos,  &
    coll_scatterFilePos,    &
    coll_fstImpactFilePos,  &
    coll_flukImpFilePos,    &
    coll_flukImpAllFilePos, &
    coll_jawProfileFilePos, &
    coll_efficFilePos,      &
    coll_efficDPFilePos,    &
    coll_effic2DFilePos,    &
    coll_summaryFilePos,    &
    coll_ampFilePos,        &
    coll_orbitCheckFilePos, &
    coll_sigmaSetFilePos,   &
    coll_impactFilePos,     &
    coll_trackoutPos,       &
    lux,                    &
    seed,                   &
    k1,                     &
    k2,                     &
    n_tot_absorbed,         &
    n_absorbed,             &
    nabs_total,             &
    nsurvive,               &
    nsurvive_end,           &
    num_selhit,             &
    (part_hit_pos(j)   ,j=1,npart), &
    (part_hit_turn(j)  ,j=1,npart), &
    (part_abs_pos(j)   ,j=1,npart), &
    (part_abs_turn(j)  ,j=1,npart), &
    (part_select(j)    ,j=1,npart), &
    (nabs_type(j)      ,j=1,npart), &
    (nhit_stage(j)     ,j=1,npart), &
    (part_linteract(j) ,j=1,npart), &
    (part_indiv(j)     ,j=1,npart), &
    (part_impact(j)    ,j=1,npart), &

    (cn_impact(j)      ,j=1,cdb_nColl),  &
    (cn_absorbed(j)    ,j=1,cdb_nColl),  &
    (caverage(j)       ,j=1,cdb_nColl),  &
    (csigma(j)         ,j=1,cdb_nColl),  &
    (gap_rms_error(j)  ,j=1,cdb_nColl),  &
    (csum(j)           ,j=1,cdb_nColl),  &
    (csqsum(j)         ,j=1,cdb_nColl)

  if(dowrite_efficiency) then
    write(fileUnit,err=100)                                      &
    (xineff(j)          ,j=1,npart),                             &
    (yineff(j)          ,j=1,npart),                             &
    (xpineff(j)         ,j=1,npart),                             &
    (ypineff(j)         ,j=1,npart),                             &
    ((counted_r(i,j)    ,j=1,numeff),i=1,npart),                 & ! (npart,numeff)
    ((counted_x(i,j)    ,j=1,numeff),i=1,npart),                 & ! (npart,numeff)
    ((counted_y(i,j)    ,j=1,numeff),i=1,npart),                 & ! (npart,numeff)
    ((counteddpop(i,j)  ,j=1,numeffdpop),i=1,npart),             & ! (npart,numeffdpop)
    (((counted2d(i,j,k) ,k=1,numeffdpop),j=1,numeff),i=1,npart), & ! (npart,numeff,numeffdpop)
    (npartdpop(j)       ,j=1,numeffdpop),                        & ! (numeffdpop)
    (neff(j)            ,j=1,numeff),                            & ! (numeff)
    (rsig(j)            ,j=1,numeff),                            & ! (numeff)
    (neffdpop(j)        ,j=1,numeffdpop),                        & ! (numeffdpop)
    (dpopbins(j)        ,j=1,numeffdpop),                        & ! (numeffdpop)
    ((neff2d(i,j)       ,j=1,numeffdpop),i=1,numeff),            & ! (numeff,numeffdpop)
    (neffx(j)           ,j=1,numeff),                            & ! (numeff)
    (neffy(j)           ,j=1,numeff)                               ! (numeff)
  end if

  if(dowrite_amplitude) then
    write(fileUnit,err=100)    &
    (nampl(j)      ,j=1,nblz), &
    (sum_ax(j)     ,j=1,nblz), &
    (sqsum_ax(j)   ,j=1,nblz), &
    (sum_ay(j)     ,j=1,nblz), &
    (sqsum_ay(j)   ,j=1,nblz), &
    (xbob(j)       ,j=1,nblz), &
    (ybob(j)       ,j=1,nblz), &
    (xpbob(j)      ,j=1,nblz), &
    (ypbob(j)      ,j=1,nblz)
  end if
!end write()

  flush(fileunit)
  return

100 continue
  lerror = .true.
  write(lout, "(a,i0,a)") "CR_POINT> ERROR Writing C/R file fort.",fileUnit," in COLLIMATION"
  write(crlog,"(a,i0,a)") "CR_POINT> ERROR Writing C/R file fort.",fileUnit," in COLLIMATION"
  flush(crlog)

end subroutine coll_crpoint

! ================================================================================================================================ !
subroutine coll_crstart()

  use coll_common
  use mod_ranlux
  use mod_alloc
  use coll_db, only : cdb_nColl

  implicit none

!copy file positions
  fort208Pos             = fort208Pos_CR
  coll_survivalFilePos   = coll_survivalFilePos_CR
  coll_gapsFilePos       = coll_gapsFilePos_CR
  coll_settingsFilePos   = coll_settingsFilePos_CR
  coll_positionsFilePos  = coll_positionsFilePos_CR
  coll_tracksFilePos     = coll_tracksFilePos_CR
  coll_pencilFilePos     = coll_pencilFilePos_CR
  coll_cryEntFilePos     = coll_cryEntFilePos_CR
  coll_cryExitFilePos    = coll_cryExitFilePos_CR
  coll_cryInterFilePos   = coll_cryInterFilePos_CR
  coll_ellipseFilePos    = coll_ellipseFilePos_CR
  coll_allImpactFilePos  = coll_allImpactFilePos_CR
  coll_allAbsorbFilePos  = coll_allAbsorbFilePos_CR
  coll_scatterFilePos    = coll_scatterFilePos_CR
  coll_fstImpactFilePos  = coll_fstImpactFilePos_CR
  coll_flukImpFilePos    = coll_flukImpFilePos_CR
  coll_flukImpAllFilePos = coll_flukImpAllFilePos_CR
  coll_jawProfileFilePos = coll_jawProfileFilePos_CR
  coll_efficFilePos      = coll_efficFilePos_CR
  coll_efficDPFilePos    = coll_efficDPFilePos_CR
  coll_effic2DFilePos    = coll_effic2DFilePos_CR
  coll_summaryFilePos    = coll_summaryFilePos_CR
  coll_ampFilePos        = coll_ampFilePos_CR
  coll_orbitCheckFilePos = coll_orbitCheckFilePos_CR
  coll_sigmaSetFilePos   = coll_sigmaSetFilePos_CR
  coll_impactFilePos     = coll_impactFilePos_CR
  coll_trackoutPos       = coll_trackoutPos_CR

!reset ranlux RNG state
  call rluxgo(lux_CR, seed_CR, k1_CR, k2_CR)

  n_tot_absorbed = n_tot_absorbed_cr
  n_absorbed     = n_absorbed_cr
  nabs_total     = nabs_total_cr
  nsurvive       = nsurvive_cr
  nsurvive_end   = nsurvive_end_cr
  num_selhit     = num_selhit_cr

  part_hit_pos(1:npart)   = part_hit_pos_cr(1:npart)
  part_hit_turn(1:npart)  = part_hit_turn_cr(1:npart)
  part_abs_pos(1:npart)   = part_abs_pos_cr(1:npart)
  part_abs_turn(1:npart)  = part_abs_turn_cr(1:npart)

  part_select(1:npart)    = part_select_cr(1:npart)
  nabs_type(1:npart)      = nabs_type_cr(1:npart)
  nhit_stage(1:npart)     = nhit_stage_cr(1:npart)
  part_linteract(1:npart) = part_linteract_cr(1:npart)
  part_indiv(1:npart)     = part_indiv_cr(1:npart)
  part_impact(1:npart)    = part_impact_cr(1:npart)

  cry_proc(1:npart)       = cry_proc_cr(1:npart)
  cry_proc_prev(1:npart)  = cry_proc_prev_cr(1:npart)
  cry_proc_tmp(1:npart)   = cry_proc_tmp_cr(1:npart)

  cn_impact(1:cdb_nColl)     = cn_impact_cr(1:cdb_nColl)
  cn_absorbed(1:cdb_nColl)   = cn_absorbed_cr(1:cdb_nColl)
  caverage(1:cdb_nColl)      = caverage_cr(1:cdb_nColl)
  csigma(1:cdb_nColl)        = csigma_cr(1:cdb_nColl)
  gap_rms_error(1:cdb_nColl) = gap_rms_error_cr(1:cdb_nColl)
  csum(1:cdb_nColl)          = csum_cr(1:cdb_nColl)
  csqsum(1:cdb_nColl)        = csqsum_cr(1:cdb_nColl)

  if(dowrite_amplitude) then
    nampl(1:nblz)    = nampl_cr(1:nblz)
    sum_ax(1:nblz)   = sum_ax_cr(1:nblz)
    sqsum_ax(1:nblz) = sqsum_ax_cr(1:nblz)
    sum_ay(1:nblz)   = sum_ay_cr(1:nblz)
    sqsum_ay(1:nblz) = sqsum_ay_cr(1:nblz)
    xbob(1:nblz)     = xbob_cr(1:nblz)
    ybob(1:nblz)     = ybob_cr(1:nblz)
    xpbob(1:nblz)    = xpbob_cr(1:nblz)
    ypbob(1:nblz)    = ypbob_cr(1:nblz)
  end if

  if(dowrite_efficiency) then
   xineff(1:npart)                          = xineff_cr(1:npart)
   yineff(1:npart)                          = yineff_cr(1:npart)
   xpineff(1:npart)                         = xpineff_cr(1:npart)
   ypineff(1:npart)                         = ypineff_cr(1:npart)
   counted_r(1:npart,1:numeff)              = counted_r_cr(1:npart,1:numeff)
   counted_x(1:npart,1:numeff)              = counted_x_cr(1:npart,1:numeff)
   counted_y(1:npart,1:numeff)              = counted_y_cr(1:npart,1:numeff)
   counteddpop(1:npart,1:numeffdpop)        = counteddpop_cr(1:npart,1:numeffdpop)
   counted2d(1:npart,1:numeff,1:numeffdpop) = counted2d_cr(1:npart,1:numeff,1:numeffdpop)
   npartdpop(1:numeffdpop)                  = npartdpop_cr(1:numeffdpop)
   neff(1:numeff)                           = neff_cr(1:numeff)
   rsig(1:numeff)                           = rsig_cr(1:numeff)
   neffdpop(1:numeffdpop)                   = neffdpop_cr(1:numeffdpop)
   dpopbins(1:numeffdpop)                   = dpopbins_cr(1:numeffdpop)
   neff2d(1:numeff,1:numeffdpop)            = neff2d_cr(1:numeff,1:numeffdpop)
   neffx(1:numeff)                          = neffx_cr(1:numeff)
   neffy(1:numeff)                          = neffy_cr(1:numeff)
  end if

  call dealloc(part_hit_pos_cr,   "part_hit_pos_cr")
  call dealloc(part_hit_turn_cr,  "part_hit_turn_cr")
  call dealloc(part_abs_pos_cr ,  "part_abs_pos_cr")
  call dealloc(part_abs_turn_cr,  "part_abs_turn_cr")
  call dealloc(part_select_cr,    "part_select_cr")
  call dealloc(nabs_type_cr,      "nabs_type_cr")
  call dealloc(nhit_stage_cr,     "nhit_stage_cr")
  call dealloc(part_linteract_cr, "part_linteract_cr")
  call dealloc(part_indiv_cr,     "part_indiv_cr")
  call dealloc(part_impact_cr,    "part_impact_cr")

  call dealloc(cry_proc_cr,       "cry_proc_cr")
  call dealloc(cry_proc_prev_cr,  "cry_proc_prev_cr")
  call dealloc(cry_proc_tmp_cr,   "cry_proc_tmp_cr")

  call dealloc(cn_impact_cr,      "cn_impact")
  call dealloc(cn_absorbed_cr,    "cn_absorbed")
  call dealloc(caverage_cr,       "caverage")
  call dealloc(csigma_cr,         "csigma")
  call dealloc(gap_rms_error_cr,  "gap_rms_error")
  call dealloc(csum_cr,           "csum")
  call dealloc(csqsum_cr,         "csqsum")

  if(dowrite_amplitude) then
    call dealloc(nampl_cr,    "nampl_cr")
    call dealloc(sum_ax_cr,   "sum_ax_cr")
    call dealloc(sqsum_ax_cr, "sqsum_ax_cr")
    call dealloc(sum_ay_cr,   "sum_ay_cr")
    call dealloc(sqsum_ay_cr, "sqsum_ay_cr")
    call dealloc(xbob_cr,     "xbob_cr")
    call dealloc(ybob_cr,     "ybob_cr")
    call dealloc(xpbob_cr,    "xpbob_cr")
    call dealloc(ypbob_cr,    "ypbob_cr")
  end if

  if(dowrite_efficiency) then
    call dealloc(xineff_cr,      "xineff_cr")
    call dealloc(yineff_cr,      "yineff_cr")
    call dealloc(xpineff_cr,     "xpineff_cr")
    call dealloc(ypineff_cr,     "ypineff_cr")
    call dealloc(counted_r_cr,   "counted_r_cr")
    call dealloc(counted_x_cr,   "counted_x_cr")
    call dealloc(counted_y_cr,   "counted_y_cr")
    call dealloc(counteddpop_cr, "counteddpop_cr")
    call dealloc(counted2d_cr,   "counted2d_cr")
    call dealloc(npartdpop_cr,   "npartdpop_cr")
    call dealloc(neff_cr,        "neff_cr")
    call dealloc(rsig_cr,        "rsig_cr")
    call dealloc(neffdpop_cr,    "neffdpop_cr")
    call dealloc(dpopbins_cr,    "dpopbins_cr")
    call dealloc(neff2d_cr,      "neff2d_cr")
    call dealloc(neffx_cr,       "neffx_cr")
    call dealloc(neffy_cr,       "neff_cr")
  end if

end subroutine coll_crstart

#endif
! ================================================================================================================================ !
!  End Checkpoint Restart
! ================================================================================================================================ !

end module collimation
