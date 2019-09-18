! ============================================================================ !
!  Collimation K2 Physics Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! ============================================================================ !
module coll_k2

  use floatPrecision
  use numerical_constants

  implicit none

  real(kind=fPrec), parameter :: tlcut = 0.0009982_fPrec

  integer,          private, save :: mcurr
  integer,          private, save :: mat
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

contains

!>
!! "Merlin" scattering collimation configuration
!! This routine pre-calcuates some varibles for
!! the nuclear properties
!<
subroutine k2coll_merlinInit

  use coll_materials

  integer i

! compute the electron densnity and plasma energy for each material
  do i=1, nmat
    edens(i) = k2coll_calcElectronDensity(zatom(i),rho(i),anuc(i))
    pleng(i) = k2coll_calcPlasmaEnergy(edens(i))
  end do

end subroutine k2coll_merlinInit

!>
!! K2 scattering collimation configuration
!<
subroutine k2coll_init
!nothing currently
end subroutine k2coll_init

! ================================================================================================ !
!  Collimation K2 Routine
! ~~~~~~~~~~~~~~~~~~~~~~~~
!  G. ROBERT-DEMOLAIZE, November 1st, 2004
!  Based on routines by JBJ. Changed by RA 2001
! ================================================================================================ !
subroutine k2coll_collimate(icoll, iturn, ie, c_length, c_rotation, c_aperture, c_offset, c_tilt,  &
  x_in, xp_in, y_in, yp_in, p_in, s_in, enom, lhit_pos, lhit_turn, part_abs_pos_local,             &
  part_abs_turn_local, impact, indiv, lint, onesided, flagsec, j_slices, nabs_type, linside)

  use parpro
  use crcoall
  use coll_db
  use coll_common
  use coll_materials
  use mod_common, only : iexact, napx
  use mod_common_main, only : partID
  use mathlib_bouncer
  use mod_ranlux
#ifdef HDF5
  use hdf5_output
#endif

  integer,          intent(in)    :: icoll        ! Collimator ID
  integer,          intent(in)    :: iturn        ! Turn number
  integer,          intent(in)    :: ie           ! Structure element index

  real(kind=fPrec), intent(in)    :: c_length     ! Collimator length in m
  real(kind=fPrec), intent(in)    :: c_rotation   ! Collimator rotation angle vs vertical in radians
  real(kind=fPrec), intent(in)    :: c_aperture   ! Collimator aperture in m
  real(kind=fPrec), intent(in)    :: c_offset     ! Collimator offset in m
  real(kind=fPrec), intent(inout) :: c_tilt(2)    ! Collimator tilt in radians

  real(kind=fPrec), intent(inout) :: x_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: xp_in(npart) ! Particle coordinate
  real(kind=fPrec), intent(inout) :: y_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: yp_in(npart) ! Particle coordinate
  real(kind=fPrec), intent(inout) :: p_in(npart)  ! Particle coordinate
  real(kind=fPrec), intent(inout) :: s_in(npart)  ! Particle coordinate

  real(kind=fPrec), intent(in)    :: enom         ! Reference momentum in GeV
  logical,          intent(in)    :: onesided

  integer,          intent(inout) :: lhit_pos(npart)
  integer,          intent(inout) :: lhit_turn(npart)
  integer,          intent(inout) :: part_abs_pos_local(npart)
  integer,          intent(inout) :: part_abs_turn_local(npart)
  integer,          intent(inout) :: nabs_type(npart)
  integer,          intent(inout) :: flagsec(npart)
  real(kind=fPrec), intent(inout) :: indiv(npart)
  real(kind=fPrec), intent(inout) :: lint(npart)
  real(kind=fPrec), intent(inout) :: impact(npart)
  logical,          intent(inout) :: linside(napx)

  ! Internal Variables

  logical hit
  integer nprim,j,nabs,nhit,j_slices

  real(kind=fPrec) keeps,fracab,drift_length,mirror,tiltangle
  real(kind=fPrec) x00,z00,p,sp,s
  real(kind=fPrec) x_flk,xp_flk,y_flk,yp_flk,zpj
  real(kind=fPrec) x_Dump,xpDump,y_Dump,ypDump,s_Dump
  real(kind=fPrec) cRot,sRot,cRRot,sRRot
  real(kind=fPrec) s_impact,xinn,xpinn,yinn,ypinn

  ! Initilaisation

  mat    = cdb_cMaterialID(icoll)
  length = c_length
  p0     = enom

  ! Initialise scattering processes
  call k2coll_scatin(p0)

  nhit   = 0
  fracab = zero
  mirror = one

  ! Compute rotation factors for collimator rotation
  cRot   = cos_mb(c_rotation)
  sRot   = sin_mb(c_rotation)
  cRRot  = cos_mb(-c_rotation)
  sRRot  = sin_mb(-c_rotation)

  do j=1,napx

    if(part_abs_pos_local(j) /= 0 .and. part_abs_turn_local(j) /= 0) then
      ! Don't do scattering process for particles already absorbed
      cycle
    end if

    impact(j) = -one
    lint(j)   = -one
    indiv(j)  = -one

    x      = x_in(j)
    xp     = xp_in(j)
    z      = y_in(j)
    zp     = yp_in(j)
    p      = p_in(j)
    sp     = zero
    dpop   = (p - p0)/p0
    x_flk  = zero
    y_flk  = zero
    xp_flk = zero
    yp_flk = zero

    ! Transform particle coordinates to get into collimator coordinate  system
    ! First do rotation into collimator frame
    x  =  x_in(j)*cRot + sRot*y_in(j)
    z  =  y_in(j)*cRot - sRot*x_in(j)
    xp = xp_in(j)*cRot + sRot*yp_in(j)
    zp = yp_in(j)*cRot - sRot*xp_in(j)

    ! For one-sided collimators consider only positive X. For negative X jump to the next particle
    if((onesided .and. x < zero) .and. ((icoll /= ipencil) .or. (iturn /= 1))) then
      cycle
    end if

    ! Now mirror at the horizontal axis for negative X offset
    if(x < zero) then
      mirror    = -one
      tiltangle = -one*c_tilt(2)
    end if
    if(x >= zero) then
      mirror    = one
      tiltangle = c_tilt(1)
    end if
    x  = mirror * x
    xp = mirror * xp

    ! Shift with opening and offset
    x = (x - c_aperture/two) - mirror*c_offset

    ! Include collimator tilt
    if(tiltangle > zero) then
      xp = xp - tiltangle
    end if
    if(tiltangle < zero) then
      x  = x + sin_mb(tiltangle) * c_length
      xp = xp - tiltangle
    end if

    ! For selected collimator, first turn reset particle distribution to simple pencil beam
    nprim = 3
    if(((icoll == ipencil .and. iturn == 1) .or. (iturn == 1 .and. ipencil == 999 .and. icoll <= nprim .and. &
       (j >= (icoll-1)*napx/nprim) .and. (j <= (icoll)*napx/nprim))) .and. (pencil_distr /= 3)) then
      ! RB addition : don't go in this if-statement if pencil_distr=3. This distribution is generated in main loop instead

      ! TW why did I set this to 0, seems to be needed for getting
      !    right amplitude => no "tilt" of jaw for the first turn !!!!
      c_tilt(1) = zero
      c_tilt(2) = zero

      ! Standard pencil beam as implemented by GRD ------- TW
      if(pencil_rmsx == zero .and. pencil_rmsy == zero) then
        x  = pencil_dx(icoll)
        xp = zero
        z  = zero
        zp = zero
      end if

      ! Rectangular (pencil-beam) sheet-beam with  ------ TW
      ! pencil_offset is the rectangulars center
      ! pencil_rmsx defines spread of impact parameter
      ! pencil_rmsy defines spread parallel to jaw surface
      if(pencil_distr == 0 .and.(pencil_rmsx /= zero .or. pencil_rmsy /= zero)) then
        ! how to assure that all generated particles are on the jaw ?!
        x  = pencil_dx(icoll)+pencil_rmsx*(real(rndm4(),fPrec)-half)
        xp = zero
        z  = pencil_rmsy*(real(rndm4(),fPrec)-half)
        zp = zero
      end if

      ! Gaussian (pencil-beam) sheet-beam with ------- TW
      ! pencil_offset is the mean  gaussian distribution
      ! pencil_rmsx defines spread of impact parameter
      ! pencil_rmsy defines spread parallel to jaw surface
      if(pencil_distr == 1 .and. (pencil_rmsx /= zero .or. pencil_rmsy /= zero)) then
        x  = pencil_dx(icoll) + pencil_rmsx*ran_gauss(two)
        ! all generated particles are on the jaw now
        x  = sqrt(x**2)
        xp = zero
        z  = pencil_rmsy*ran_gauss(two)
        zp = zero
      end if

      ! Gaussian (pencil-beam) sheet-beam with ------- TW
      ! pencil_offset is the mean  gaussian distribution
      ! pencil_rmsx defines spread of impact parameter
      !             here pencil_rmsx is not gaussian!!!
      ! pencil_rmsy defines spread parallel to jaw surface
      if(pencil_distr == 2 .and. (pencil_rmsx /= zero .or. pencil_rmsy /= zero)) then
        x  = pencil_dx(icoll) + pencil_rmsx*(real(rndm4(),fPrec)-half)
        ! all generated particles are on the jaw now
        x  = sqrt(x**2)
        xp = zero
        z  = pencil_rmsy*ran_gauss(two)
        zp = zero
      end if

      ! Selection of pos./neg. jaw  implemented by GRD ---- TW
      ! ensure that for onesided only particles on pos. jaw are created
      if(onesided) then
        mirror = one
      else
        if(rndm4() < half) then
          mirror = -one
        else
          mirror = one
        end if
      end if

      ! TW if c_tilt is set to zero before entering pencil beam
      !    section the assigning of the tilt will result in assigning zeros
      if(mirror < zero) then
        tiltangle = c_tilt(2)
      else
        tiltangle = c_tilt(1)
      end if
!
      write(coll_pencilUnit,"(f10.8,4(2x,f10.8))") x, xp, z, zp, tiltangle
    end if !if(( (icoll.eq.ipencil .and. iturn.eq.1) .or. (itu

    ! SR, 18-08-2005: after finishing the coordinate transformation,
    ! or the coordinate manipulations in case of pencil beams,
    ! write down the initial coordinates of the impacting particles
    xinn  = x
    xpinn = xp
    yinn  = z
    ypinn = zp

    ! particle passing above the jaw are discarded => take new event
    ! entering by the face, shorten the length (zlm) and keep track of
    ! entrance longitudinal coordinate (keeps) for histograms

    ! The definition is that the collimator jaw is at x>=0.

    ! 1) Check whether particle hits the collimator
    hit   = .false.
    s     = zero
    keeps = zero
    zlm   = -one * length

    if(x >= zero) then
      ! Particle hits collimator and we assume interaction length ZLM equal
      ! to collimator length (what if it would leave collimator after
      ! small length due to angle???)
      zlm       = length
      impact(j) = x
      indiv(j)  = xp
    else if(xp <= zero) then
      ! Particle does not hit collimator. Interaction length ZLM is zero.
      zlm = zero
    else
      ! Calculate s-coordinate of interaction point
      s = (-one*x)/xp
      if(s <= zero) then
        write(lerr,"(a)") "COLLK2> ERROR S <= zero. This should not happen!"
        call prror
      end if
      if(s < length) then
        zlm       = length - s
        impact(j) = zero
        indiv(j)  = xp
      else
        zlm = zero
      end if
    end if

    ! First do the drift part
    ! DRIFT PART
    drift_length = length - zlm
    if(drift_length > zero) then
      if(iexact) then
        zpj = sqrt(one-xp**2-zp**2)
        x   = x + drift_length*(xp/zpj)
        z   = z + drift_length*(zp/zpj)
        sp  = sp + drift_length
      else
        x  = x + xp* drift_length
        z  = z + zp * drift_length
        sp = sp + drift_length
      end if
    end if

    ! Now do the scattering part
    if(zlm > zero) then
      if(.not.linside(j)) then
        ! first time particle hits collimator: entering jaw
        linside(j) = .true.
        if(dowrite_impact) then
          if(tiltangle > zero) then
            x_Dump = (x+c_aperture/two+tiltangle*sp)*mirror+c_offset
          else
            x_Dump = (x+c_aperture/two+tiltangle*(sp-c_length))*mirror+c_offset
          end if
          xpDump = (xp+tiltangle)*mirror
          y_Dump = z
          ypDump = zp
          s_Dump = sp+real(j_slices-1,fPrec)*c_length
          write(coll_jawProfileUnit,"(3(1x,i7),5(1x,e17.9),1x,i1)") &
            icoll,iturn,partID(j),x_Dump,xpDump,y_Dump,ypDump,s_Dump,1
        end if
      end if

      s_impact = sp
      nhit = nhit + 1
      call k2coll_jaw(s,nabs,icoll,iturn,partID(j))

      nabs_type(j) = nabs

      ! SR+GRD: CREATE A FILE TO CHECK THE VALUES OF IMPACT PARAMETERS
      ! SR, 29-08-2005: Add to the longitudinal coordinates the position of the slice beginning
      if(dowrite_impact) then
        if(flagsec(j) == 0) then
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
            write(coll_fstImpactUnit,"(i5,1x,i7,1x,i2,1x,i1,2(1x,f5.3),8(1x,e17.9))") &
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

      ! SR, 18-08-2005: add also the initial coordinates of the impacting particles!
      lhit_pos(j)  = ie
      lhit_turn(j) = iturn

      ! September 2006
      ! If particle is absorbed then set x and y to 99.99 mm
      ! SR: before assigning new (x,y) for nabs=1, write the inelastic impact file .

      ! RB: writeout should be done for both inelastic and single diffractive. doing all transformations
      ! in x_flk and making the set to 99.99 mm conditional for nabs=1

      if(dowrite_impact .or. nabs == 1 .or. nabs == 4) then
        ! transform back to lab system for writeout.
        ! keep x,y,xp,yp unchanged for continued tracking, store lab system variables in x_flk etc

        x_flk  = xInt
        xp_flk = xpInt

        if(tiltangle > zero) then
          x_flk  = x_flk  + tiltangle*(sInt+sp)
          xp_flk = xp_flk + tiltangle
        else if(tiltangle < zero) then
          xp_flk = xp_flk + tiltangle
          x_flk  = x_flk - sin_mb(tiltangle) * (length -(sInt+sp))
        end if

        x_flk  = (x_flk + c_aperture/two) + mirror*c_offset
        x_flk  = mirror * x_flk
        xp_flk = mirror * xp_flk
        y_flk  = yInt   * cRRot - x_flk  * sRRot
        yp_flk = ypInt  * cRRot - xp_flk * sRRot
        x_flk  = x_flk  * cRRot + yInt   * sRRot
        xp_flk = xp_flk * cRRot + ypInt  * sRRot

        if(dowrite_impact) then
          ! write out all impacts to all_impacts.dat
          write(coll_flukImpAllUnit,"(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))") &
          icoll,c_rotation, (sInt+sp)+(real(j_slices,fPrec)-one)*c_length,             &
          x_flk*c1e3, xp_flk*c1e3, y_flk*c1e3, yp_flk*c1e3, nabs, partID(j), iturn
        end if

        ! Standard FLUKA_impacts writeout of inelastic and single diffractive
        if(nabs == 1 .or. nabs == 4) then
          if(dowrite_impact) then
            write(coll_flukImpUnit,"(i4,(1x,f6.3),(1x,f8.6),4(1x,e19.10),i2,2(1x,i7))") &
              icoll,c_rotation,                                                 &
              sInt + sp + (real(j_slices,fPrec)-one) * c_length,                &
              x_flk*c1e3, xp_flk*c1e3, y_flk*c1e3, yp_flk*c1e3,                 &
              nabs,partID(j),iturn
          end if

          ! Finally, the actual coordinate change to 99 mm
          if(nabs == 1) then
            fracab = fracab + 1
            x      = 99.99e-3_fPrec
            z      = 99.99e-3_fPrec
            part_abs_pos_local(j)  = ie
            part_abs_turn_local(j) = iturn
            lint(j) = zlm
          end if
        end if
      end if
    end if !if (zlm.gt.0.) then

    ! Do the rest drift, if particle left collimator early
    ! DRIFT PART
    if(nabs /= 1 .and. zlm > zero) then
      drift_length = (length-(s+sp))
      if(drift_length > c1m15) then
        linside(j) = .false.
        if(dowrite_impact) then
          if(tiltangle > zero) then
            x_Dump = (x+c_aperture/two+tiltangle*(s+sp))*mirror+c_offset
          else
            x_Dump = (x+c_aperture/two+tiltangle*(s+sp-c_length))*mirror+c_offset
          end if
          xpDump = (xp+tiltangle)*mirror
          y_Dump = z
          ypDump = zp
          s_Dump = s+sp+real(j_slices-1,fPrec)*c_length
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
    if(x < 99.0d-3) then
      ! Include collimator tilt
      if(tiltangle > zero) then
        x  = x  + tiltangle*c_length
        xp = xp + tiltangle
      else if(tiltangle < zero) then
        x  = x + tiltangle*c_length
        xp = xp + tiltangle
        x  = x - sin_mb(tiltangle) * c_length
      end if

      ! Transform back to particle coordinates with opening and offset
      z00 = z
      x00 = x + mirror*c_offset
      x   = (x + c_aperture/two) + mirror*c_offset

      ! Now mirror at the horizontal axis for negative X offset
      x  = mirror * x
      xp = mirror * xp

      ! Last do rotation into collimator frame
      x_in(j)  =  x*cRRot +  z*sRRot
      y_in(j)  =  z*cRRot -  x*sRRot
      xp_in(j) = xp*cRRot + zp*sRRot
      yp_in(j) = zp*cRRot - xp*sRRot

      if(((icoll == ipencil .and. iturn == 1) .or. &
         (iturn == 1 .and. ipencil == 999 .and. icoll <= nprim .and. (j >= (icoll-1)*napx/nprim) .and. &
         (j <= (icoll)*napx/nprim))) .and.(pencil_distr /= 3)) then
        ! RB: adding condition that this shouldn't be done if pencil_distr=3
        x00      = mirror * x00
        x_in(j)  = x00*cRRot + z00*sRRot
        y_in(j)  = z00*cRRot - x00*sRRot
        xp_in(j) = xp_in(j) + mirror*xp_pencil0
        yp_in(j) = yp_in(j) + mirror*yp_pencil0
        x_in(j)  = x_in(j)  + mirror*x_pencil(icoll)
        y_in(j)  = y_in(j)  + mirror*y_pencil(icoll)
      end if

      p_in(j) = (one + dpop) * p0
      s_in(j) = sp + (real(j_slices,fPrec)-one) * c_length
    else
      x_in(j) = x
      y_in(j) = z
    end if
  end do ! End of loop over all particles

end subroutine k2coll_collimate

!>
!! k2coll_scatin(plab)
!! Configure the K2 scattering routine cross sections
!<
subroutine k2coll_scatin(plab)

  use mod_funlux
  use coll_materials
  use mathlib_bouncer
  use physical_constants

  real(kind=fPrec), intent(in) :: plab

  real(kind=fPrec) tlow,thigh
  integer ma,i

  ecmsq = (two * pmap) * plab
#ifndef MERLINSCATTER
  xln15s = log_mb(0.15_fPrec*ecmsq)
  ! Claudia Fit from COMPETE collaboration points "arXiv:hep-ph/0206172v1 19Jun2002"
  pptot = 0.041084_fPrec-0.0023302_fPrec*log_mb(ecmsq)+0.00031514_fPrec*log_mb(ecmsq)**2
  ! Claudia used the fit from TOTEM for ppel (in barn)
  ppel = (11.7_fPrec-1.59_fPrec*log_mb(ecmsq)+0.134_fPrec*log_mb(ecmsq)**2)/c1e3
  ! Claudia updated SD cross that cointains renormalized pomeron flux (in barn)
  ppsd = (4.3_fPrec+0.3_fPrec*log_mb(ecmsq))/c1e3
#endif

#ifdef MERLINSCATTER
  ! No crlibm...
  call merlinscatter_setup(plab,rnd_seed)
  call merlinscatter_setdata(pptot,ppel,ppsd)
#endif

  ! Claudia new fit for the slope parameter with new data at sqrt(s)=7 TeV from TOTEM
  bpp = 7.156_fPrec+1.439_fPrec*log_mb(sqrt(ecmsq))

  ! unmeasured tungsten data,computed with lead data and power laws
  bnref(4) = bnref(5)*(anuc(4) / anuc(5))**(two/three)
  emr(4) = emr(5) * (anuc(4)/anuc(5))**(one/three)

  ! Compute cross-sections (CS) and probabilities + Interaction length
  ! Last two material treated below statement number 100

  tlow = tlcut
  do ma=1,nrmat

    mcurr = ma
    ! prepare for Rutherford differential distribution
    thigh = hcut(ma)
    call funlxp(k2coll_ruth, cgen(1,ma), tlow, thigh)

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

    ! Nuclear Elastic is TOT-inel-qel ( see definition in RPP)
    csect(2,ma) = ((csect(0,ma) - csect(1,ma)) - csect(3,ma)) - csect(4,ma)
    csect(5,ma) = csref(5,ma)

    ! Now add Coulomb
    csect(0,ma) = csect(0,ma) + csect(5,ma)

    ! Interaction length in meter
    xintl(ma) = (c1m2*anuc(ma))/(((fnavo * rho(ma))*csect(0,ma))*1e-24_fPrec)

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

end subroutine k2coll_scatin

!>
!! jaw(s,nabs,icoll,iturn,ipart)
!!     RB: adding as input arguments to jaw variables icoll,iturn,ipart
!!         these are only used for the writeout of particle histories
!!
!! Input:   ZLM is interaction length
!!          MAT is choice of material
!!
!! Output:  nabs = 1   Particle is absorped
!!          nabs = 4   Single-diffractive scattering
!!          dpop       Adjusted for momentum loss (dE/dx)
!!          s          Exit longitudinal position
!!
!! Physics:  If monte carlo interaction length greater than input
!!           interaction length, then use input interaction length
!!           Is that justified???
!<
subroutine k2coll_jaw(s, nabs, icoll, iturn, ipart)

  use mod_ranlux
  use coll_common
  use coll_materials
  use mathlib_bouncer
#ifdef HDF5
  use hdf5_output
#endif

  real(kind=fPrec), intent(inout) :: s
  integer,          intent(inout) :: nabs
  integer,          intent(in)    :: icoll
  integer,          intent(in)    :: iturn
  integer,          intent(in)    :: ipart

  real(kind=fPrec) m_dpodx,p,rlen,t,dxp,dzp,p1,zpBef,xpBef,pBef
  integer inter,nabs_tmp

  ! Note that the input parameter is dpop. Here the momentum p is constructed out of this input.
  p    = p0*(one+dpop)
  nabs = 0
  nabs_tmp = nabs

  if(mat == nmat) then
    ! Collimator treated as black absorber
    nabs = 1
    nabs_tmp = nabs
    s = zero

    if(dowrite_impact) then
      ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
      if(h5_useForCOLL) then
        call coll_hdf5_writeCollScatter(icoll, iturn, ipart, nabs_tmp, -one, zero, zero)
      else
#endif
      write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))') icoll, iturn, ipart, nabs_tmp, -one, zero, zero
#ifdef HDF5
      end if
#endif
    end if
    return
  else if(mat == nmat-1) then
    ! Collimator treated as drift
    s = zlm
    x = x+s*xp
    z = z+s*zp

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

  ! Initialize the interaction length to input interaction length
  rlen = zlm

  ! Do a step for a point-like interaction.
  ! Get monte-carlo interaction length.
10 continue
  zlm1     = (-one*xintl(mat))*log_mb(real(rndm4(),fPrec))
  nabs_tmp = 0  ! type of interaction reset before following scattering process
  xpBef    = xp ! save angles and momentum before scattering
  zpBef    = zp
  pBef     = p

  ! If the monte-carlo interaction length is longer than the
  ! remaining collimator length, then put it to the remaining
  ! length, do multiple coulomb scattering and return.
  ! LAST STEP IN ITERATION LOOP
  if(zlm1 > rlen) then
    zlm1 = rlen
    call k2coll_mcs(s)
    s = (zlm-rlen)+s
#ifdef MERLINSCATTER
    call merlinscatter_calc_ion_loss(p,edens(mat), pleng(mat),exenergy(mat),s,m_dpodx)
    p = p-m_dpodx
#else
    call k2coll_calcIonLoss(mat,p,rlen,m_dpodx)  ! DM routine to include tail
    p = p-m_dpodx*s
#endif

    dpop = (p-p0)/p0
    if(dowrite_impact) then
      ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
      if(h5_useForCOLL) then
        call coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef)
      else
#endif
      write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))') icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
#ifdef HDF5
      end if
#endif
    end if
    return
  end if

  ! Otherwise do multi-coulomb scattering.
  ! REGULAR STEP IN ITERATION LOOP
  call k2coll_mcs(s)

  ! Check if particle is outside of collimator (X.LT.0) after
  ! MCS. If yes, calculate output longitudinal position (s),
  ! reduce momentum (output as dpop) and return.
  ! PARTICLE LEFT COLLIMATOR BEFORE ITS END.
  if(x <= zero) then
    s = (zlm-rlen)+s

#ifdef MERLINSCATTER
    call merlinscatter_calc_ion_loss(p,edens(mat),pleng(mat),exenergy(mat),s,m_dpodx)
    p = p-m_dpodx
#else
    call k2coll_calcIonLoss(mat,p,rlen,m_dpodx)
    p = p-m_dpodx*s
#endif
    dpop = (p-p0)/p0

    if(dowrite_impact) then
      ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
      if(h5_useForCOLL) then
        call coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef)
      else
#endif
      write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e18.10))') icoll,iturn,ipart,nabs_tmp,(p-pBef)/pBef,xp-xpBef,zp-zpBef
#ifdef HDF5
      end if
#endif
    end if

    return
  end if

  ! Check whether particle is absorbed. If yes, calculate output
  ! longitudinal position (s), reduce momentum (output as dpop)
  ! and return.
  ! PARTICLE WAS ABSORBED INSIDE COLLIMATOR DURING MCS.

  inter = k2coll_ichoix(mat)
  nabs  = inter
  nabs_tmp = nabs

  ! RB, DM: save coordinates before interaction for writeout to FLUKA_impacts.dat
  xInt  = x
  xpInt = xp
  yInt  = z
  ypInt = zp
  sInt  = (zlm-rlen)+zlm1

  if(inter == 1) then
    s = (zlm-rlen)+zlm1

#ifdef MERLINSCATTER
    call merlinscatter_calc_ion_loss(p,edens(mat),pleng(mat),exenergy(mat),s,m_dpodx)
    p = p-m_dpodx
#else
    call k2coll_calcIonLoss(mat,p,rlen,m_dpodx)
    p = p-m_dpodx*s
#endif

    dpop = (p-p0)/p0

    ! write coll_scatter.dat for complete scattering histories
#ifdef HDF5
    if(h5_useForCOLL) then
      call coll_hdf5_writeCollScatter(icoll,iturn,ipart,nabs_tmp,-one,zero,zero)
    else
#endif
    write(coll_scatterUnit,'(1x,i2,2x,i4,2x,i5,2x,i1,3(2x,e14.6))') icoll,iturn,ipart,nabs_tmp,-one,zero,zero
#ifdef HDF5
    end if
#endif
    return
  end if

  ! Now treat the other types of interaction, as determined by ICHOIX:

  ! Nuclear-Elastic:          inter = 2
  ! pp Elastic:               inter = 3
  ! Single-Diffractive:       inter = 4    (changes momentum p)
  ! Coulomb:                  inter = 5

  ! As the single-diffractive interaction changes the momentum, save input momentum in p1.
  p1 = p

  ! Gettran returns some monte carlo number, that, as I believe, gives the rms transverse momentum transfer.
  t = k2coll_gettran(inter,mat,p)

  ! Tetat calculates from the rms transverse momentum transfer in
  ! monte-carlo fashion the angle changes for x and z planes. The
  ! angle change is proportional to SQRT(t) and 1/p, as expected.
  call k2coll_tetat(t,p,dxp,dzp)

  ! Apply angle changes
  xp = xp+dxp
  zp = zp+dzp

  ! Treat single-diffractive scattering.
  if(inter == 4) then

    ! added update for s
    s    = (zlm-rlen)+zlm1
    xpsd = dxp
    zpsd = dzp
    psd  = p1

    ! Add this code to get the momentum transfer also in the calling routine
    dpop = (p-p0)/p0
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
    end if
#endif
  end if

  ! Calculate the remaining interaction length and close the iteration loop.
  rlen = rlen-zlm1
  goto 10

end subroutine k2coll_jaw

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
subroutine k2coll_mcs(s)

  use coll_materials

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
  call k2coll_soln3(ae,be,dh,rlen,s)
  if(s.lt.h) s=h
  call k2coll_scamcs(x,xp,s,radl_mat)
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
20    call k2coll_scamcs(z,zp,s,radl_mat)
  s=s*radl(mat)
  x=(x*theta)*radl(mat)
  xp=xp*theta
  z=(z*theta)*radl(mat)
  zp=zp*theta
end subroutine k2coll_mcs

!>
!! k2coll_calcIonLoss(IS,PC,DZ,EnLo)
!! subroutine for the calculazion of the energy loss by ionization
!! Either mean energy loss from Bethe-Bloch, or higher energy loss, according to finite probability from cross section
!! written by DM for crystals, introduced in main code by RB
!<
subroutine k2coll_calcIonLoss(IS, PC, DZ, EnLo)

! IS material ID
! PC momentum in GeV
! DZ length traversed in material (meters)
! EnLo energy loss in GeV/meter

  use physical_constants
  use mathlib_bouncer
  use mod_ranlux
  use coll_materials

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

end subroutine k2coll_calcIonLoss

!>
!! k2coll_tetat(t,p,tx,tz)
!! ???
!!
!<
subroutine k2coll_tetat(t,p,tx,tz)

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
end subroutine k2coll_tetat

!>
!! k2coll_soln3(a,b,dh,smax,s)
!! ???
!<
subroutine k2coll_soln3(a,b,dh,smax,s)

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
      call k2coll_iterat(a,b,dh,s)
    end if
  else
    c=(-one*a)/b
    if(smax.lt.c) then
      if(smax**3.le.(a+b*smax)**2) then
        s=smax
        return
      else
        s=smax*half
        call k2coll_iterat(a,b,dh,s)
      end if
    else
      s=c*half
      call k2coll_iterat(a,b,dh,s)
    end if
  end if

end subroutine k2coll_soln3

!>
!! k2coll_scamcs(xx,xxp,s,radl_mat)
!! ???
!<
subroutine k2coll_scamcs(xx,xxp,s,radl_mat)

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
end subroutine k2coll_scamcs

subroutine k2coll_iterat(a,b,dh,s)

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

end subroutine k2coll_iterat

!>
!! k2coll_ruth(t)
!! Calculate the rutherford scattering cross section
!<
real(kind=fPrec) function k2coll_ruth(t)

  use mathlib_bouncer
  use coll_materials

  implicit none

  real(kind=fPrec) t,cnorm,cnform
  parameter(cnorm=2.607e-5_fPrec,cnform=0.8561e3_fPrec) ! DM: changed 2.607d-4 to 2.607d-5 to fix Rutherford bug

  k2coll_ruth=(cnorm*exp_mb(((-one*real(t,fPrec))*cnform)*emr(mcurr)**2))*(zatom(mcurr)/real(t,fPrec))**2
end function k2coll_ruth

!>
!! k2coll_gettran(inter,xmat,p)
!! This function determines: GETTRAN - rms transverse momentum transfer
!! Note: For single-diffractive scattering the vector p of momentum
!! is modified (energy loss is applied)
!<
real(kind=fPrec) function k2coll_gettran(inter,xmat,p)

  use mathlib_bouncer
  use mod_ranlux
  use mod_funlux
  use coll_materials

  implicit none

  integer, intent(in) :: inter,xmat
  real(kind=fPrec) :: p

  integer :: length
  real(kind=fPrec) :: t,xm2,bsd
  real(kind=fPrec) :: truth,xran(1)

  ! Neither if-statements below have an else, so defaultingfuction return to zero.
  k2coll_gettran = zero ! -Wmaybe-uninitialized

! inter=2: Nuclear Elastic, 3: pp Elastic, 4: Single Diffractive, 5:Coulomb
#ifndef MERLINSCATTER
  if( inter.eq.2 ) then
    k2coll_gettran = (-one*log_mb(real(rndm4(),fPrec)))/bn(xmat)

  else if( inter .eq. 3 ) then
    k2coll_gettran = (-one*log_mb(real(rndm4(),fPrec)))/bpp

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
      k2coll_gettran = (-one*log_mb(real(rndm4(),fPrec)))/bsd

  else if( inter.eq.5 ) then
    length=1
    call funlux( cgen(1,mat), xran, length)
    truth=xran(1)
    t=real(truth,fPrec)
    k2coll_gettran = t
  end if
#else

  if( inter.eq.2 ) then
    k2coll_gettran = (-one*log_mb(real(rndm4(),fPrec)))/bn(xmat)

  else if( inter .eq. 3 ) then
    call merlinscatter_get_elastic_t(k2coll_gettran)

  else if( inter .eq. 4 ) then
    call merlinscatter_get_sd_xi(xm2)
    call merlinscatter_get_sd_t(k2coll_gettran)
    p = p  * (one - (xm2/ecmsq))

  else if ( inter.eq.5 ) then
    length=1
    call funlux( cgen(1,mat) , xran, length)
    truth=xran(1)
    t=real(truth,fPrec)
    k2coll_gettran = t
  end if

#endif
  return
end function k2coll_gettran

!>
!! k2coll_ichoix(ma)
!! Select a scattering type (elastic, sd, inelastic, ...)
!<
integer function k2coll_ichoix(ma)
  use mod_ranlux
  use coll_materials
  integer ma,i
  real(kind=fPrec) aran
  aran=real(rndm4(),fPrec)
  i=1
10 if( aran.gt.cprob(i,ma) ) then
    i=i+1
    goto 10
  end if
  k2coll_ichoix=i
end function k2coll_ichoix

!>
!! k2coll_calcElectronDensity(AtomicNumber, Density, AtomicMass)
!! Function to calculate the electron density in a material
!! Should give the number per cubic meter
!<
real(kind=fPrec) function k2coll_calcElectronDensity(AtomicNumber, Density, AtomicMass)
  real(kind=fPrec) AtomicNumber, Density, AtomicMass
  real(kind=fPrec) Avogadro
  real(kind=fPrec) PartA, PartB
  parameter (Avogadro = 6.022140857e23_fPrec)
  PartA = AtomicNumber * Avogadro * Density
  !1e-6 factor converts to n/m^-3
  PartB = AtomicMass * c1m6
  k2coll_calcElectronDensity = PartA/PartB
  return
end function k2coll_calcElectronDensity

!>
!! k2coll_calcPlasmaEnergy(ElectronDensity)
!! Function to calculate the plasma energy in a material
!! CalculatePlasmaEnergy = (PlanckConstantBar * sqrt((ElectronDensity *(ElectronCharge**2)) / &
!!& (ElectronMass * FreeSpacePermittivity)))/ElectronCharge*eV;
!<
real(kind=fPrec) function k2coll_calcPlasmaEnergy(ElectronDensity)

  real(kind=fPrec) ElectronDensity
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
  k2coll_calcPlasmaEnergy=PlanckConstantBar*sqrtAB/ElectronCharge*c1m9

end function k2coll_calcPlasmaEnergy

end module coll_k2
