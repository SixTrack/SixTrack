module cheby
  use parpro
  use floatPrecision
  use crcoall
  use numerical_constants, only : zero
  implicit none

  ! A.Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 25-02-2019
  ! module for handling maps expressed by Chebyshev polynomials

  integer, allocatable, save  :: icheby(:)            ! index of chebyshev lens
  integer, parameter          :: ncheby=20            ! max number of chebyshev lenses treated (in LATTICE structure)
  integer, save               :: mcheby=0             ! last chebyshev lens read
  integer, parameter          :: ncheby_tables=20     ! max number of chebyshev tables in memory (in SINGLE ELEMENT array)
  integer, save               :: mcheby_tables=0      ! last chebyshev table read

  ! variables to save parameters for tracking etc.
  integer, save          :: cheby_itable(ncheby)=0         ! index of chebyshev table
  real(kind=fPrec), save :: cheby_offset_x(ncheby)=zero    ! hor. offset [mm]
  real(kind=fPrec), save :: cheby_offset_y(ncheby)=zero    ! ver. offset [mm]
  real(kind=fPrec), save :: cheby_angle(ncheby)=zero       ! rotation angle about the longitudinal axis [rad]
  real(kind=fPrec), save :: cheby_scalingFact(ncheby)=zero ! scaling factor []

  ! tables with chebyshev coefficients
  integer, parameter     :: cheby_max_order=30        ! max order of chebyshev polynomials
  character(len=mFNameLen), save:: cheby_filename(ncheby_tables) = " "! file names
  real(kind=fPrec), save :: cheby_coeffs(0:cheby_max_order,0:cheby_max_order,ncheby_tables) = zero ! coefficients
  integer, save          :: cheby_maxOrder(ncheby_tables)  ! max order of the current map
  real(kind=fPrec), save :: cheby_refCurr(ncheby_tables)   = zero ! reference current [A]
  real(kind=fPrec), save :: cheby_refRadius(ncheby_tables) = zero ! reference radius [mm]
  real(kind=fPrec), save :: cheby_refBeta(ncheby_tables)   = zero ! reference e-beta []

contains

subroutine cheby_allocate_arrays
  implicit none
  integer stat
  call alloc(icheby,nele,0,'icheby')
end subroutine cheby_allocate_arrays

subroutine cheby_expand_arrays(nele_new)
  implicit none
  integer, intent(in) :: nele_new
  call alloc(icheby,nele_new,0,'icheby')
end subroutine cheby_expand_arrays


subroutine cheby_postInput
  use mathlib_bouncer
  use utils

  integer j
  logical exist
  
  ! Parse files with coefficients for Chebyshev polynomials:
   do j=1,mcheby
    inquire(file=cheby_filename(j), exist=exist)
    if(.not. exist) then
      write(lout,"(a)") "CHEBY> ERROR Problems with file with coefficients for Chebyshev polynominals: ", &
            trim(cheby_filename(j))
      call prror(-1)
    end if
    call parseChebyFile(j)
  end do

end subroutine cheby_postInput


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
  beta_b=e0f*gammar/pma
  brho_b=e0f/(clight*c1m6)

  do j=1,napx

    ! apply offset
    xx=xv1(j)-cheby_offset_x(jcheby)
    yy=xv2(j)-cheby_offset_y(jcheby)

    ! check that particle is within the domain of chebyshev polynomials
    rr=sqrt(xx**2+yy**2)
    if (rr.gt.cheby_refRadius(cheby_itable(jcheby))) then
      write(lout,"(a,3(e12.6,1x),a,e12.6)") "CHEBY> ERROR in cheby_kick: particle at position (x,y,r): ",&
        xv1(j), xv2(j), rr,' is outside radial domain of Chebyshev polinomials: ',cheby_refRadius(cheby_itable(jcheby))
      call prror
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

end subroutine cheby_kick


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

! ================================================================================================ !
!  Last modified: 2019-02-25
!  Rewritten by VKBO, June 2018
!  Read file with coefficients for chebyshev polynomials
!  ifile is index of file in table of chebyshev files
!  file is structured as:
!    keyword : value
!  keywords:
!  - I: reference current intensity of e-beam [A]
!  - Ek: reference kinetic energy of electrons [keV]
!  - rad: reference radius [mm]
!  comment line is headed by '#'
!  coefficients are give with the following syntax:
!  i j : value
!  where i->x and j->y
! ================================================================================================ !
subroutine parseChebyFile(ifile)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use physical_constants
  use crcoall
  use mod_common
  use mod_settings
  use string_tools
  use mod_units

  implicit none

  integer, intent(in) :: ifile

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mInputLn) inLine
  integer nSplit
  logical spErr,err

  integer iErr, ii, jj, fUnit
  real(kind=fPrec) tmpflt, beta, gamma

  ierr = 0
  write(lout,"(a)") "CHEBY> Parsing file with coefficients for Chebyshev polynomials "//trim(cheby_filename(ifile))
  call f_requestUnit(cheby_filename(ifile),fUnit)
  call f_open(unit=fUnit,file=cheby_filename(ifile),mode='r',err=err,formatted=.true.,status="old")

10 continue
  read(fUnit,"(a)",end=20,err=30) inLine
  if(inLine(1:1) == "#") goto 10

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "CHEBY> ERROR Failed to Chebyshev input line."
    goto 30
  end if

  if(inLine(1:1) == "I") then
    ! Read reference current of e-beam in e-lens
    if(nSplit < 3) then
      iErr = 2
      goto 30
    end if
    call chr_cast(lnSplit(3),tmpflt,spErr)
    cheby_refCurr(ifile) = tmpflt

  else if(inLine(1:2) == "Ek") then
    ! Read reference kinetic energy of e-beam in e-lens
    if(nSplit < 3) then
      iErr = 3
      goto 30
    end if
    call chr_cast(lnSplit(3),tmpflt,spErr)
    gamma = (tmpflt*c1m3)/pmae+one ! from kinetic energy
    cheby_refBeta(ifile) = sqrt((gamma+one)*(gamma-one))/gamma

  else if(inLine(1:3) == "rad") then
    ! Read reference radius e-beam in e-lens
    if(nSplit < 3) then
      iErr = 4
      goto 30
    end if
    call chr_cast(lnSplit(3),tmpflt,spErr)
    cheby_refRadius(ifile) = tmpflt

  else
    ! Read chebyshev coefficients
    if(nSplit /= 4) then
      iErr = 5
      goto 30
    end if
    call chr_cast(lnSplit(1),ii,spErr)
    if(ii > cheby_max_order) then
      iErr = 6
      goto 30
    end if
    call chr_cast(lnSplit(2),jj,spErr)
    if(jj > cheby_max_order) then
      iErr = 7
      goto 30
    end if
    call chr_cast(lnSplit(4),tmpflt,spErr)
    cheby_coeffs(ii,jj,ifile) = tmpflt

  end if ! close if for keyword identification
  goto 10

20 continue

  call f_close(fUnit)

  if(st_quiet < 2) then
    ! Echo parsed data (unless told to be quiet!)
    write(lout,"(a,i0)") "CHEBY> Coefficients for Chebyshev polynomials as from file "//&
      trim(cheby_filename(ifile))//" - #",ifile
    write(lout,"(a,e22.15)") "ELENS> * Reference current [A] : ",cheby_refCurr(ifile)
    write(lout,"(a,e22.15)") "ELENS> * reference beta     [] : ",cheby_refBeta(ifile)
    write(lout,"(a,e22.15)") "ELENS> * reference radius [mm] : ",cheby_refRadius(ifile)
    do ii=0,cheby_max_order
      do jj=0,cheby_max_order
        if(cheby_coeffs(ii,jj,ifile)/= zero) then
          write(lout,"(2(a,i4),a,e22.15)") "CHEBY> Order ",ii,",",jj," : ",cheby_coeffs(ii,jj,ifile)
        end if
      end do
    end do
  end if
  return

30 continue
  write(lout,"(a,i0,a)") "CHEBY> ERROR ",iErr," while parsing file "//trim(cheby_filename(ifile))
  call prror(-1)

end subroutine parseChebyFile


end module cheby
