! ============================================================================ !
!  Collimation K2 Physics Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ============================================================================ !
module coll_k2

  use coll_common
  use floatPrecision
  use numerical_constants

  implicit none

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
subroutine collimate2(c_material, c_length, c_rotation, c_aperture, c_offset, c_tilt,      &
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

end module coll_k2
