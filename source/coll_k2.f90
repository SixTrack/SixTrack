! ============================================================================ !
!  Collimation K2 Physics Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ============================================================================ !
module coll_k2

  use coll_common
  use floatPrecision
  use numerical_constants

  implicit none

  real(kind=fPrec), parameter :: tlcut = 0.0009982_fPrec

  integer,          private, save :: mat
  integer,          private, save :: nev
  integer,          private, save :: mcurr
  real(kind=fPrec), private, save :: csect(0:5,nmat)
  real(kind=fPrec), private, save :: xintl(nmat)
  real(kind=fPrec), private, save :: bn(nmat)
  real(kind=fPrec), private, save :: freep(nmat)
  real(kind=fPrec), private, save :: cgen(200,nmat)
  real(kind=fPrec), private, save :: zlm
  real(kind=fPrec), private, save :: zlm1
  real(kind=fPrec), private, save :: p0
  real(kind=fPrec), private, save :: length
  real(kind=fPrec), private, save :: xln15s
  real(kind=fPrec), private, save :: pptot
  real(kind=fPrec), private, save :: ppsd
  real(kind=fPrec), private, save :: ppel
  real(kind=fPrec), private, save :: ecmsq
  real(kind=fPrec), private, save :: bpp
  real(kind=fPrec), private, save :: xpsd
  real(kind=fPrec), private, save :: zpsd
  real(kind=fPrec), private, save :: psd

  ! RB DM 2014 added variables for FLUKA output
  real(kind=fPrec), private, save :: xInt,xpInt,yInt,ypInt,sInt
  real(kind=fPrec), private, save :: x,xp,z,zp,dpop

  ! Cross section inputs and material property database
  ! GRD CHANGED ON 2/2003 TO INCLUDE CODE FOR C, C2 from JBJ (rwa)
  ! Total number of materials are defined in nmat
  ! Number of real materials are defined in nrmat
  ! The last materials in nmat are 'vacuum' and 'black',see in sub. SCATIN
  ! Reference data at pRef=450Gev

  ! pp cross-sections and parameters for energy dependence
  real(kind=fPrec), parameter :: pptref = 0.04_fPrec
  real(kind=fPrec), parameter :: freeco = 1.618_fPrec

  ! Mean excitation energy (GeV) values added by Claudia for Bethe-Bloch implementation:
  real(kind=fPrec), parameter :: exenergy(nmat) = [ &
    63.7e-9_fPrec, 166.0e-9_fPrec, 322.0e-9_fPrec, 727.0e-9_fPrec, 823.0e-9_fPrec, 78.0e-9_fPrec, 78.0e-9_fPrec, &
    87.1e-9_fPrec, 152.9e-9_fPrec, 424.0e-9_fPrec, 320.8e-9_fPrec, 682.2e-9_fPrec, zero, c1e10 ]

  ! GRD IMPLEMENT CHANGES FROM JBJ, 2/2003 RWA
  real(kind=fPrec), public,  save :: anuc(nmat)  = &
    [ 9.01_fPrec,  26.98_fPrec,  63.55_fPrec, 183.85_fPrec, 207.19_fPrec,    12.01_fPrec,  12.01_fPrec,  &
     13.53_fPrec,  25.24_fPrec,  95.96_fPrec,  63.15_fPrec, 166.70_fPrec,     zero,         zero         ]
  real(kind=fPrec), public , save :: zatom(nmat) = &
    [ 4.00_fPrec,  13.00_fPrec,  29.00_fPrec,  74.00_fPrec,  82.00_fPrec,     6.00_fPrec,   6.00_fPrec,  &
      6.65_fPrec,  11.90_fPrec,  42.00_fPrec,  28.80_fPrec,  67.70_fPrec,     zero,         zero         ]
  real(kind=fPrec), public,  save :: rho(nmat)   = &
    [ 1.848_fPrec,  2.70_fPrec,   8.96_fPrec,  19.30_fPrec,   11.35_fPrec,    1.67_fPrec,   4.52_fPrec,  &
      2.500_fPrec,  5.40_fPrec,  10.22_fPrec,   8.93_fPrec,   18.00_fPrec,    zero,         zero         ]
  real(kind=fPrec), private, save :: emr(nmat)   = &
    [ 0.22_fPrec,   0.302_fPrec,  0.366_fPrec,  0.520_fPrec,   0.542_fPrec,   0.25_fPrec,   0.25_fPrec,  &
      0.25_fPrec,   0.308_fPrec,  0.481_fPrec,  0.418_fPrec,   0.578_fPrec,   zero,         zero         ]
  real(kind=fPrec), private, save :: hcut(nmat)  = &
    [ 0.02_fPrec,   0.02_fPrec,   0.01_fPrec,   0.01_fPrec,    0.01_fPrec,    0.02_fPrec,    0.02_fPrec, &
      0.02_fPrec,   0.02_fPrec,   0.02_fPrec,   0.02_fPrec,    0.02_fPrec,    zero,          zero        ]
  real(kind=fPrec), private, save :: radl(nmat)  = &
    [ 0.353_fPrec,  0.089_fPrec,  0.0143_fPrec, 0.0035_fPrec,  0.0056_fPrec,  0.2557_fPrec, 0.094_fPrec, &
      0.1193_fPrec, 0.0316_fPrec, 0.0096_fPrec, 0.0144_fPrec,  0.00385_fPrec, 1.0e12_fPrec, 1.0e12_fPrec ]

  ! Nuclear elastic slope from Schiz et al.,PRD 21(3010)1980
  ! MAY06-GRD value for Tungsten (W) not stated. Last 2 ones interpolated
  real(kind=fPrec), private, save :: bnref(nmat) = &
    [74.7_fPrec, 120.3_fPrec, 217.8_fPrec, 440.3_fPrec, 455.3_fPrec, 70.0_fPrec, 70.0_fPrec, &
     76.7_fPrec, 115.0_fPrec, 273.9_fPrec, 208.7_fPrec, 392.1_fPrec, zero,       zero        ]

  ! All cross-sections are in barns,nuclear values from RPP at 20geV
  ! Coulomb is integerated above t=tLcut[Gev2] (+-1% out Gauss mcs)

  ! in Cs and CsRef,1st index: Cross-sections for processes
  ! 0:Total, 1:absorption, 2:nuclear elastic, 3:pp or pn elastic
  ! 4:Single Diffractive pp or pn, 5:Coulomb for t above mcs

  ! Claudia 2013: updated cross section values. Unit: Barn. New 2013:
  real(kind=fPrec), private, save :: csref(0:5,nmat)
  data csref(0,1), csref(1,1), csref(5,1) /0.271_fPrec, 0.192_fPrec, 0.0035e-2_fPrec/
  data csref(0,2), csref(1,2), csref(5,2) /0.643_fPrec, 0.418_fPrec, 0.0340e-2_fPrec/
  data csref(0,3), csref(1,3), csref(5,3) /1.253_fPrec, 0.769_fPrec, 0.1530e-2_fPrec/
  data csref(0,4), csref(1,4), csref(5,4) /2.765_fPrec, 1.591_fPrec, 0.7680e-2_fPrec/
  data csref(0,5), csref(1,5), csref(5,5) /3.016_fPrec, 1.724_fPrec, 0.9070e-2_fPrec/
  data csref(0,6), csref(1,6), csref(5,6) /0.337_fPrec, 0.232_fPrec, 0.0076e-2_fPrec/
  data csref(0,7), csref(1,7), csref(5,7) /0.337_fPrec, 0.232_fPrec, 0.0076e-2_fPrec/
  data csref(0,8), csref(1,8), csref(5,8) /0.362_fPrec, 0.247_fPrec, 0.0094e-2_fPrec/
  data csref(0,9), csref(1,9), csref(5,9) /0.572_fPrec, 0.370_fPrec, 0.0279e-2_fPrec/
  data csref(0,10),csref(1,10),csref(5,10)/1.713_fPrec, 1.023_fPrec, 0.2650e-2_fPrec/
  data csref(0,11),csref(1,11),csref(5,11)/1.246_fPrec, 0.765_fPrec, 0.1390e-2_fPrec/
  data csref(0,12),csref(1,12),csref(5,12)/2.548_fPrec, 1.473_fPrec, 0.5740e-2_fPrec/

  ! Cprob to choose an interaction in iChoix
  real(kind=fPrec), private, save :: cprob(0:5,nmat)
  data cprob(0,1:nmat)/nmat*zero/
  data cprob(5,1:nmat)/nmat*one/

contains

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
subroutine collimate2(icoll, iturn, ie, c_material, c_length, c_rotation, c_aperture, c_offset, c_tilt,      &
  x_in, xp_in, y_in, yp_in, p_in, s_in, np, enom, lhit_pos, lhit_turn, part_abs_pos_local,   &
  part_abs_turn_local, impact, indiv, lint, onesided, flagsec, j_slices, nabs_type, linside)

  use crcoall
  use parpro
  use coll_common
  use mod_common, only : iexact, napx
  use mod_common_main, only : partID
  use mathlib_bouncer
  use mod_ranlux
#ifdef HDF5
  use hdf5_output
#endif

  implicit none

  integer, intent(in) :: icoll
  integer, intent(in) :: iturn
  integer, intent(in) :: ie

  character(len=4), intent(in)    :: c_material  ! material
  real(kind=fPrec), intent(in)    :: c_length    ! length in m
  real(kind=fPrec), intent(in)    :: c_rotation  ! rotation angle vs vertical in radian
  real(kind=fPrec), intent(in)    :: c_aperture  ! aperture in m
  real(kind=fPrec), intent(in)    :: c_offset    ! offset in m
  real(kind=fPrec), intent(inout) :: c_tilt(2)   ! tilt in radians

  logical onesided,hit
  integer nprim,j,nabs,nhit,np

  integer, allocatable :: lhit_pos(:) !(npart)
  integer, allocatable :: lhit_turn(:) !(npart)
  integer, allocatable :: part_abs_pos_local(:) !(npart)
  integer, allocatable :: part_abs_turn_local(:) !(npart)
  integer, allocatable :: nabs_type(:) !(npart)
!MAY2005

  logical linside(napx)
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

  real(kind=fPrec) x00,z00,p,sp,s,enom

!AUGUST2006 Added ran_gauss for generation of pencil/     ------- TW
!           sheet beam distribution  (smear in x and y)
!

  real(kind=fPrec) x_flk,xp_flk,y_flk,yp_flk,zpj
  real(kind=fPrec) x_Dump,xpDump,y_Dump,ypDump,s_Dump

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
    call prror
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
    if(x.lt.zero) then
      mirror = -one
      tiltangle = -one*c_tilt(2)
    end if

    if(x.ge.zero) then
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
    x  = (x - c_aperture/two) - mirror*c_offset
!
!++  Include collimator tilt
!
    if(tiltangle.gt.zero) then
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
      if(pencil_rmsx.eq.zero .and. pencil_rmsy.eq.zero) then
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

      if(pencil_distr.eq.2 .and.(pencil_rmsx.ne.zero.or.pencil_rmsy.ne.zero )) then
        x  = pencil_dx(icoll) + pencil_rmsx*(real(rndm4(),fPrec)-half)
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
      if(mirror.lt.zero) then
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
      write(coll_pencilUnit,'(f10.8,(2x,f10.8),(2x,f10.8),(2x,f10.8),(2x,f10.8))') x, xp, z, zp, tiltangle

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
    s     = zero
    keeps = zero
    zlm   = -one * length

    if(x.ge.zero) then

!++  Particle hits collimator and we assume interaction length ZLM equal
!++  to collimator length (what if it would leave collimator after
!++  small length due to angle???)
      zlm = length
      impact(j) = x
      indiv(j) = xp
    else if(xp.le.zero) then

!++  Particle does not hit collimator. Interaction length ZLM is zero.
      zlm = zero
    else
!
!++  Calculate s-coordinate of interaction point
      s = (-one*x) / xp
      if(s.le.0) then
        write(lout,*) 'S.LE.0 -> This should not happen'
        call prror
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
    if(drift_length.gt.zero) then
      if(iexact) then
        zpj = sqrt(one-xp**2-zp**2)
        x = x + drift_length*(xp/zpj)
        z = z + drift_length*(zp/zpj)
        sp = sp + drift_length
      else
        x  = x + xp* drift_length
        z  = z + zp * drift_length
        sp = sp + drift_length
      end if
    end if
!
!++  Now do the scattering part
!
    if (zlm.gt.zero) then
      if(.not.linside(j)) then
        ! first time particle hits collimator: entering jaw
        linside(j)=.true.
        if(dowrite_impact) then
          if ( tiltangle.gt.zero ) then
            x_Dump=(x+c_aperture/two+tiltangle*sp)*mirror+c_offset
          else
            x_Dump=(x+c_aperture/two+tiltangle*(sp-c_length))*mirror+c_offset
          end if
          xpDump=(xp+tiltangle)*mirror
          y_Dump=z
          ypDump=zp
          s_Dump=sp+real(j_slices-1,fPrec)*c_length
          write(coll_jawProfileUnit,"(3(1x,i7),5(1x,e17.9),1x,i1)") &
            icoll,iturn,partID(j),x_Dump,xpDump,y_Dump,ypDump,s_Dump,1
        end if
      end if
!JUNE2005
      s_impact = sp
!JUNE2005
      nhit = nhit + 1
!            WRITE(*,*) J,X,XP,Z,ZP,SP,DPOP
!     RB: add new input arguments to jaw icoll,iturn,partID for writeout
      call jaw(s,nabs,icoll,iturn,partID(j))

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
            call h5_writeData(coll_hdf5_fstImpacts, 1,  1, partID(j))
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
            write(coll_fstImpactUnit,'(i5,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))') &
              partID(j),iturn,icoll,nabs,                             &
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
        write(coll_flukImpAllUnit,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))') &
     &              icoll,c_rotation,                                   &
     &              sInt + sp + (real(j_slices,fPrec)-one) * c_length,  &
     &              x_flk*c1e3, xp_flk*c1e3, y_flk*c1e3, yp_flk*c1e3,   &
     &              nabs,partID(j),iturn
      end if

! standard FLUKA_impacts writeout of inelastic and single diffractive
      if((nabs.eq.1).OR.(nabs.eq.4)) then

!     SR, 29-08-2005: Include the slice numer!
        if(dowrite_impact) then
          write(coll_flukImpUnit,'(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))') &
     &icoll,c_rotation,                                                 &
     &sInt + sp + (real(j_slices,fPrec)-one) * c_length,                &
     &x_flk*c1e3, xp_flk*c1e3, y_flk*c1e3, yp_flk*c1e3,                 &
     &nabs,partID(j),iturn
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
        linside(j)=.false.
        if(dowrite_impact) then
          if ( tiltangle.gt.zero ) then
            x_Dump=(x+c_aperture/two+tiltangle*(s+sp))*mirror+c_offset
          else
            x_Dump=(x+c_aperture/two+tiltangle*(s+sp-c_length))*mirror+c_offset
          end if
          xpDump=(xp+tiltangle)*mirror
          y_Dump=z
          ypDump=zp
          s_Dump=s+sp+real(j_slices-1,fPrec)*c_length
          write(coll_jawProfileUnit,"(3(1x,i7),5(1x,e17.9),1x,i1)") &
            icoll,iturn,partID(j),x_Dump,xpDump,y_Dump,ypDump,s_Dump,2
        end if
        if(iexact) then
          zpj = sqrt(one-xp**2-zp**2)
          x   = x + drift_length*(xp/zpj)
          z   = z + drift_length*(zp/zpj)
          sp  = sp + drift_length
        else
          x  = x + xp * drift_length
          z  = z + zp * drift_length
          sp = sp + drift_length
        end if
      end if
      lint(j) = zlm - drift_length
    end if

!++  Transform back to particle coordinates with opening and offset
    if(x.lt.99.0d-3) then

!++  Include collimator tilt
      if(tiltangle.gt.zero) then
        x  = x  + tiltangle*c_length
        xp = xp + tiltangle
      else if(tiltangle.lt.zero) then
        x  = x + tiltangle*c_length
        xp = xp + tiltangle
        x  = x - sin_mb(tiltangle) * c_length
      end if

!++  Transform back to particle coordinates with opening and offset
      z00 = z
      x00 = x + mirror*c_offset
      x = (x + c_aperture/two) + mirror*c_offset

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
      s_in(j) = sp + (real(j_slices,fPrec)-one) * c_length

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
!c$$$     &            s+sp + (dble(j_slices)-1d0) * c_length,               &
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
!! scatin(plab)
!! Configure the K2 scattering routine cross sections
!!
!<
subroutine scatin(plab)

  use physical_constants
  use mathlib_bouncer
  use mod_funlux
#ifdef HDF5
  use hdf5_output
  use hdf5_tracks2
#endif

  implicit none

  integer ma,i
  real(kind=fPrec) plab
  real(kind=fPrec) tlow,thigh

  ecmsq = (two * pmap) * plab
#ifndef MERLINSCATTER
  xln15s=log_mb(0.15_fPrec*ecmsq)

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

  tlow=tlcut
  do ma=1,nrmat
    mcurr=ma
! prepare for Rutherford differential distribution
    thigh=hcut(ma)
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
    bn(ma) = (bnref(ma) * csect(0,ma)) / csref(0,ma)
! also correct inel-CS
    csect(1,ma) = (csref(1,ma) * csect(0,ma)) / csref(0,ma)
!
! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect(2,ma) = ((csect(0,ma) - csect(1,ma)) - csect(3,ma)) - csect(4,ma)
    csect(5,ma) = csref(5,ma)
! Now add Coulomb
    csect(0,ma) = csect(0,ma) + csect(5,ma)
! Interaction length in meter
  xintl(ma) = (c1m2*anuc(ma))/(((fnavo * rho(ma))*csect(0,ma))*1d-24)

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
subroutine jaw(s,nabs,icoll,iturn,ipart)

  use mathlib_bouncer
  use mod_ranlux
  use coll_common
#ifdef HDF5
  use hdf5_output
#endif

  implicit none

  integer nabs,inter,iturn,icoll,ipart,nabs_tmp ! RB: added variables icoll,iturn,ipart for writeout
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
      write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))') icoll, iturn, ipart, nabs_tmp, -one, zero, zero
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
      write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))') icoll, iturn, ipart, nabs_tmp, -one, zero, zero
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

10  zlm1=(-one*xintl(mat))*log_mb(real(rndm4(),fPrec))
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
    s=(zlm-rlen)+s
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
      write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))') icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
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
    s=(zlm-rlen)+s

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
      write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))') icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
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
  sInt=(zlm-rlen)+zlm1

  if(inter.eq.1) then
    s=(zlm-rlen)+zlm1

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
    write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))') icoll,iturn,ipart,nabs_tmp,-one,zero,zero
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
    s=(zlm-rlen)+zlm1
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
    write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))') icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
#ifdef HDF5
    endif
#endif
  end if

!++  Calculate the remaining interaction length and close the iteration
!++  loop.
  rlen=rlen-zlm1
  goto 10

end subroutine jaw

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

  x=(x/theta)/radl(mat)
  xp=xp/theta
  z=(z/theta)/radl(mat)
  zp=zp/theta
  rlen0=zlm1/radl(mat)
  rlen=rlen0
10    ae=bn0*x
  be=bn0*xp
  call soln3(ae,be,dh,rlen,s)
  if(s.lt.h) s=h
  call scamcs(x,xp,s,radl_mat)
  if(x.le.0.d0) then
   s=(rlen0-rlen)+s
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
  x=(x*theta)*radl(mat)
  xp=xp*theta
  z=(z*theta)*radl(mat)
  zp=zp*theta
end subroutine mcs

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
  use mathlib_bouncer
  use mod_ranlux

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

!>
!! tetat(t,p,tx,tz)
!! ???
!!
!<
subroutine tetat(t,p,tx,tz)

  use mod_ranlux

  implicit none

  real(kind=fPrec) t,p,tx,tz,va,vb,va2,vb2,r2,teta
  teta = sqrt(t)/p

! Generate sine and cosine of an angle uniform in [0,2pi](see RPP)
10 va  =(two*real(rndm4(),fPrec))-one
  vb = real(rndm4(),fPrec)
  va2 = va**2
  vb2 = vb**2
  r2 = va2 + vb2
  if ( r2.gt.one) go to 10
  tx = teta * ((two*va)*vb) / r2
  tz = teta * (va2 - vb2) / r2
  return
end subroutine tetat

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

!>
!! scamcs(xx,xxp,s,radl_mat)
!! ???
!<
subroutine scamcs(xx,xxp,s,radl_mat)

  use mathlib_bouncer
  use mod_ranlux

  implicit none

  real(kind=fPrec) v1,v2,r2,a,z1,z2,ss,s,xx,xxp,x0,xp0
  real(kind=fPrec) radl_mat

  x0=xx
  xp0=xxp

5 v1=2d0*real(rndm4(),fPrec)-1d0
  v2=2d0*real(rndm4(),fPrec)-1d0
  r2=v1**2+v2**2
  if(r2.ge.1.d0) goto 5

  a=sqrt((-2.d0*log_mb(r2))/r2)
  z1=v1*a
  z2=v2*a
  ss=sqrt(s)
  xx=x0+s*(xp0+(half*ss)*(one+0.038_fPrec*log_mb(s))*(z2+z1*0.577350269_fPrec)) !Claudia: added logarithmic part in mcs formula
  xxp=xp0+ss*z2*(one+0.038_fPrec*log_mb(s))
end subroutine scamcs

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
!! ruth(t)
!! Calculate the rutherford scattering cross section
!<
function ruth(t)

  use mathlib_bouncer

  implicit none

  real(kind=fPrec) ruth,t
  real(kind=fPrec) cnorm,cnform
  parameter(cnorm=2.607e-5_fPrec,cnform=0.8561e3_fPrec) ! DM: changed 2.607d-4 to 2.607d-5 to fix Rutherford bug

  ruth=(cnorm*exp_mb(((-one*real(t,fPrec))*cnform)*emr(mcurr)**2))*(zatom(mcurr)/real(t,fPrec))**2
end function ruth

!>
!! gettran(inter,xmat,p)
!! This function determines: GETTRAN - rms transverse momentum transfer
!! Note: For single-diffractive scattering the vector p of momentum
!! is modified (energy loss is applied)
!<
real(kind=fPrec) function gettran(inter,xmat,p)

  use mathlib_bouncer
  use mod_ranlux
  use mod_funlux

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
    gettran = (-one*log_mb(real(rndm4(),fPrec)))/bn(xmat)

  else if( inter .eq. 3 ) then
    gettran = (-one*log_mb(real(rndm4(),fPrec)))/bpp

  else if( inter .eq. 4 ) then
    xm2 = exp_mb( real(rndm4(),fPrec) * xln15s )
    p = p  * (one - xm2/ecmsq)
    if( xm2 .lt. two ) then
      bsd = two * bpp
    else if (( xm2 .ge. two ).and. ( xm2 .le. five )) then
      bsd = ((106.0_fPrec-17.0_fPrec*xm2) *  bpp )/ 36.0_fPrec
!    else if ( xm2 .gt. five ) then
    else !makes the compiler more happy
      bsd = (seven * bpp) / 12.0_fPrec
    end if
      gettran = (-one*log_mb(real(rndm4(),fPrec)))/bsd

  else if( inter.eq.5 ) then
    length=1
    call funlux( cgen(1,mat), xran, length)
    truth=xran(1)
    t=real(truth,fPrec)
    gettran = t
  end if
#else

  if( inter.eq.2 ) then
    gettran = (-one*log_mb(real(rndm4(),fPrec)))/bn(xmat)

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
    t=real(truth,fPrec)
    gettran = t
  end if

#endif
  return
end function gettran

!>
!! ichoix(ma)
!! Select a scattering type (elastic, sd, inelastic, ...)
!<
function ichoix(ma)
  use mod_ranlux
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

end module coll_k2
