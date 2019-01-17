module cheby
  use parpro
  use floatPrecision
  implicit none

  ! A.Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 28-02-2018
  ! module for handling lenses with kicks expressed by Chebyshev polynomials

  integer,allocatable, save  :: icheby(:)              ! index of chebyshev lens
  integer, parameter     :: ncheby=20                 ! max number of chebyshev lenses treated (in LATTICE structure)
  integer, save          :: mcheby                    ! last chebyshev lens read
  integer, parameter     :: ncheby_tables=5           ! max number of chebyshev tables in memory (in SINGLE ELEMENT array)
  integer, save          :: mcheby_tables             ! last chebyshev table read

  ! variables to save parameters for tracking etc.
  real(kind=fPrec), save :: cheby_offset_x(ncheby), cheby_offset_y(ncheby)  ! hor./vert. offset [mm]
  real(kind=fPrec), save :: cheby_angle(ncheby)       ! rotation angle about the longitudinal axis [rad]
  integer, save          :: cheby_itable(ncheby)      ! index of chebyshev table
  real(kind=fPrec), save :: cheby_scalingFact(ncheby) ! scaling factor []

  ! tables with chebyshev coefficients
  integer, parameter     :: cheby_unit=107            ! unit for reading the chebyshev coefficients
  integer, parameter     :: cheby_max_order=60        ! max order of chebyshev polynomials
  integer, parameter     :: cheby_lFileName=16        ! length of filenames
  character(len=cheby_lFileName), save:: cheby_filename(ncheby_tables)  ! file names
  real(kind=fPrec), save :: cheby_coeffs(0:cheby_max_order,0:cheby_max_order,ncheby_tables) ! coefficients
  integer, save          :: cheby_maxOrder(ncheby_tables)  ! max order of the current map
  real(kind=fPrec), save :: cheby_refRadius(ncheby_tables) ! reference radius [mm]


contains

subroutine cheby_allocate_arrays
  use crcoall
  implicit none
  integer stat
    call alloc(icheby,nele,0,'icheby')
end subroutine cheby_allocate_arrays

subroutine cheby_expand_arrays(nele_new)
  implicit none
  integer, intent(in) :: nele_new
  call alloc(icheby,nele_new,0,'icheby')
end subroutine cheby_expand_arrays


  subroutine cheby_kick(jcheby)

    ! A. Mereghetti (CERN, BE-ABP-HSS)
    ! last modified: 20-02-2018
    ! apply kick of electron lens

    use crcoall
    use mod_common
    use mod_common_main
    use mathlib_bouncer
    use numerical_constants
    use physical_constants

    real(kind=fPrec) :: xx, yy, rr, frr, dxp, dyp
    real(kind=fPrec) :: theta, radio, angle_rad, brho_b, beta_b
    integer          :: j, jcheby
    logical          :: lrotate

    angle_rad = zero ! -Wmaybe-uninitialized

    ! rotation angle
    lrotate=cheby_angle(jcheby).ne.zero

    ! Brho of beam
    beta_b=e0f*gammar/nucm0
    brho_b=e0f/(clight*c1m6)

    do j=1,napx

       ! apply offset
       xx=xv1(j)-cheby_offset_x(jcheby)
       yy=xv2(j)-cheby_offset_y(jcheby)

       ! check that particle is within the domain of chebyshev polynomials
       rr=sqrt(xx**2+yy**2)
       if (rr.gt.cheby_refRadius(cheby_itable(jcheby))) then
          write(lout,*) 'ERROR in cheby_kick: particle at position (x,y,r): ',     &
               xv1(j), xv2(j), rr,' is outside radial domain of Chebyshev polinomials: ', &
               cheby_refRadius(cheby_itable(jcheby))
          call prror(-1)
       end if

       ! in case of non-zero tilt angle, rotate coordinates
       if (lrotate) then
          theta = atan2_mb(yy, xx)-angle_rad
          xx = rr * cos_mb(theta)
          yy = rr * sin_mb(theta)
       end if

       ! apply kick
       call cheby_getKick( xx, yy, dxp, dyp, cheby_itable(jcheby), brho_b, beta_b )
       ! take into account scaling factor
       dxp=dxp *cheby_scalingFact(jcheby)
       dyp=dyp *cheby_scalingFact(jcheby)

       ! in case cheby has a non-zero angle, rotate kicks
       if (lrotate) then
          ! NB: cheby_angle(jcheby) is the rotation angle of the cheby
          theta = atan2_mb(dyp, dxp)+angle_rad
          radio = sqrt(dxp**2 + dyp**2)
          dxp = radio * cos_mb(theta)
          dyp = radio * sin_mb(theta)
       end if

       ! apply kicks, taking into account magnetic rigidity of particle being tracked;
       yv1(j)=yv1(j)+dxp *oidpsv(j)
       yv2(j)=yv2(j)+dyp *oidpsv(j)
    end do
    return

  end subroutine cheby_kick


  subroutine cheby_comnul

    ! A. Mereghetti (CERN, BE-ABP-HSS)
    ! last modified: 28-02-2018
    ! always in main code

    use mod_common
    use mod_common_main
    use numerical_constants

    integer          :: i, j, i1, i2

    do i=1,nele
       icheby(i) = 0
    end do
    mcheby=0
    do i=1,ncheby
       cheby_itable(i)      = 0
       cheby_offset_x(i)    = zero
       cheby_offset_y(i)    = zero
       cheby_angle(i)       = zero
       cheby_scalingFact(i) = one
    end do

    ! table with coefficients of chebyshev polynominals
    mcheby_tables=0
    do i=1,ncheby_tables
       do j=1,cheby_lFileName
          cheby_filename(i)(j:j)=' '
       end do
       do i1=0,cheby_max_order
          do i2=0,cheby_max_order
             cheby_coeffs(i1,i2,i)=zero
          end do
       end do
       cheby_refRadius(i)=zero
       cheby_maxOrder(i)=0
    end do
    return
  end subroutine cheby_comnul


  subroutine cheby_getKick( xx, yy, dxp, dyp, iTable, brho_b, beta_b )

    ! A. Mereghetti (CERN, BE-ABP-HSS)
    ! last modified: 28-02-2018
    ! compute kicks from Chebyshev polinomials - see FermiLAB-FN-0972-APC

    use mathlib_bouncer
    use physical_constants
    use numerical_constants

    ! interface vars
    real(kind=fPrec) :: xx, yy, dxp, dyp, brho_b, beta_b
    integer          :: iTable

    ! temp vars
    real(kind=fPrec) :: uu, vv, Tx (0:cheby_maxOrder(iTable)), Ty (0:cheby_maxOrder(iTable)), &
                        kx, ky, Tpx(0:cheby_maxOrder(iTable)), Tpy(0:cheby_maxOrder(iTable)), &
                        fu, fv
    integer          :: ii, jj

    ! normalised variables
    uu=xx/cheby_refRadius(iTable)
    vv=yy/cheby_refRadius(iTable)
    ! normalisation factors of derivatives
    fu=(one-uu)*(one+uu)
    fv=(one-vv)*(one+vv)

    ! polynomials:
    Tx(0)=one
    Ty(0)=one
    Tx(1)=uu
    Ty(1)=vv
    do ii=2,cheby_max_order
       Tx(ii)=two*uu*Tx(ii-1)-Tx(ii-2)
       Ty(ii)=two*vv*Ty(ii-1)-Ty(ii-2)
    end do
    ! derivatives:
    Tpx(0)=zero
    Tpy(0)=zero
    do ii=1,cheby_max_order
       Tpx(ii)=real(ii,fPrec)*(Tx(ii-1)-uu*Tx(ii))/fu
       Tpy(ii)=real(ii,fPrec)*(Ty(ii-1)-vv*Ty(ii))/fv
    end do

    ! get kicks
    dxp=zero
    dyp=zero
    do ii=0,cheby_maxOrder(iTable)
       do jj=0,ii
          dxp=dxp+cheby_coeffs(jj,ii,iTable)*Tpx(jj)*Ty (ii-jj)
          dyp=dyp+cheby_coeffs(jj,ii,iTable)*Tx (jj)*Tpy(ii-jj)
       end do
    end do
    dxp=-dxp/cheby_refRadius(iTable)
    dyp=-dyp/cheby_refRadius(iTable)

    ! take into account Brho and beta
    dxp=dxp/(brho_b*clight*beta_b)
    dyp=dyp/(brho_b*clight*beta_b)

   end subroutine cheby_getKick


end module cheby
