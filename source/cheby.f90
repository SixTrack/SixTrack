module cheby
  use parpro
  use floatPrecision
  use crcoall
  use mod_alloc
  use numerical_constants, only : zero, one
  implicit none

  ! A.Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 25-02-2019
  ! module for handling maps expressed by Chebyshev polynomials

  integer, allocatable, save  :: icheby(:)            ! index of chebyshev lens
  integer, parameter          :: ncheby=20            ! max number of chebyshev lenses treated (in LATTICE structure)
  integer, save               :: mcheby=0             ! last chebyshev lens read
  integer, parameter          :: ncheby_tables=20     ! max number of chebyshev tables in memory (in SINGLE ELEMENT array)
  integer, save               :: mcheby_tables=0      ! last chebyshev table read
  integer, parameter          :: kzcheby=111          ! kz of chebyshev lenses

  ! variables to save parameters for tracking etc.
  integer, save          :: cheby_itable(ncheby)=0         ! index of chebyshev table
  real(kind=fPrec), save :: cheby_r2(ncheby)     = zero    ! outer radius R2 [mm] (optional)
  real(kind=fPrec), save :: cheby_r1(ncheby)     = zero    ! inner radius R1 [mm] (optional)
  real(kind=fPrec), save :: cheby_angle(ncheby)=zero       ! rotation angle about the longitudinal axis [rad] (optional)
  real(kind=fPrec), save :: cheby_offset_x(ncheby)=zero    ! hor. offset [mm] (optional)
  real(kind=fPrec), save :: cheby_offset_y(ncheby)=zero    ! ver. offset [mm] (optional)
  real(kind=fPrec), save :: cheby_L(ncheby)=-one           ! actual length of (thick) lens [m] (optional)
  real(kind=fPrec), save :: cheby_I(ncheby)=-one           ! actual powering of lens [A] (optional)
  real(kind=fPrec), save :: cheby_scalingFact(ncheby)=one  ! scaling factor [] (computed internally)

  ! tables with chebyshev coefficients
  integer, parameter     :: cheby_max_order=30             ! max order of chebyshev polynomials currently supported
  character(len=mFNameLen), save:: cheby_filename(ncheby_tables) = " "! file names
  real(kind=fPrec), save :: cheby_coeffs(0:cheby_max_order,0:cheby_max_order,ncheby_tables) = zero ! coefficients
  integer, save          :: cheby_maxOrder(ncheby_tables)  = 0    ! max order of the current map
  real(kind=fPrec), save :: cheby_refI(ncheby_tables) = one  ! reference current [A] (optional)
  real(kind=fPrec), save :: cheby_refR(ncheby_tables) = zero ! reference radius [mm] (mandatory)
  real(kind=fPrec), save :: cheby_refL(ncheby_tables) = one  ! reference length [m]  (optional)

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


! ================================================================================================ !
!  Parse Line for Chebyshev lens
!  Last modified: 2019-02-25
! ================================================================================================ !
subroutine cheby_parseInputLine(inLine, iLine, iErr)

  use mod_settings
  use sixtrack_input
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=mStrLen) tmpch
  integer nSplit, iElem, j, chIdx, tmpi1
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "CHEBY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return
  
  if(nSplit < 2) then
    write(lout,"(a,i0)") "CHEBY> ERROR Expected at least 2 input parameters, got ",nSplit
    iErr = .true.
    return
  end if

  iElem = -1
  tmpch = trim(lnSplit(1))
  do j=1,nele
    if(bez(j) == tmpch) then
      iElem = j
      exit
    end if
  end do
  if(iElem == -1) then
    write(lout,"(a)") "CHEBY> ERROR Element '"//trim(lnSplit(1))//"' not found in single element list."
    iErr = .true.
    return
  end if

  if(kz(iElem) /= kzcheby) then
    write(lout,"(2(a,i0),a)") "CHEBY> ERROR Element type is kz(",iElem,") = ",kz(iElem)," != ",kzcheby
    iErr = .true.
    return
  end if
  if(el(iElem) /= zero .or. ek(iElem) /= zero .or. ed(iElem) /= zero) then
    write(lout,"(a)")       "CHEBY> ERROR Length el(iElem) (Chebyshev lens is treated as thin element), "//&
      "and first and second field have to be zero:"
    write(lout,"(2(a,i0))") "CHEBY>       el(",iElem,") = ",el(iElem)," != 0"
    write(lout,"(2(a,i0))") "CHEBY>       ed(",iElem,") = ",ed(iElem)," != 0"
    write(lout,"(2(a,i0))") "CHEBY>       ek(",iElem,") = ",ek(iElem)," != 0"
    iErr = .true.
    return
  end if

  mcheby = mcheby+1
  if(mcheby > ncheby) then
    write(lout,"(2(a,i0))") "CHEBY> ERROR Too many Chebyshev lenses: ",mcheby,". Max is ",ncheby
    iErr = .true.
    return
  end if

  if(icheby(iElem) /= 0) then
    write(lout,"(a)") "CHEBY> ERROR The element '"//trim(bez(iElem))//"' was defined twice."
    iErr = .true.
    return
  end if
  icheby(iElem) = mcheby

  ! File with Chebyshev polynomials
  tmpch = trim(lnSplit(2))
  ! Check if profile has already been requested:
  chIdx = -1
  do tmpi1=1,mcheby_tables
    if(tmpch == cheby_filename(tmpi1)) then
      cheby_itable(icheby(iElem)) = tmpi1
      chIdx = tmpi1
      exit
    end if
  end do
  if(chIdx == -1) then
    ! Unsuccessful search
    mcheby_tables = mcheby_tables+1
    if(mcheby_tables > ncheby_tables) then
      write(lout,"(2(a,i0))") "CHEBY> ERROR Too many files with Chebyshev polynomials: ",mcheby_tables,&
        ". Max is ",ncheby_tables
      iErr = .true.
      return
    end if
    cheby_itable(icheby(iElem)) = mcheby_tables
    cheby_filename(tmpi1) = tmpch
  end if

  ! Additional geometrical infos:
  if(nSplit >= 3) call chr_cast(lnSplit(3),cheby_r2(icheby(iElem)),iErr)
  if(nSplit >= 4) call chr_cast(lnSplit(4),cheby_r1(icheby(iElem)),iErr)
  if(nSplit >= 5) call chr_cast(lnSplit(5),cheby_angle(icheby(iElem)),iErr)
  if(nSplit >= 6) call chr_cast(lnSplit(6),cheby_offset_x(icheby(iElem)),iErr)
  if(nSplit >= 7) call chr_cast(lnSplit(7),cheby_offset_y(icheby(iElem)),iErr)
  if(nSplit >= 8) call chr_cast(lnSplit(8),cheby_L(icheby(iElem)),iErr)
  if(nSplit >= 9) call chr_cast(lnSplit(9),cheby_I(icheby(iElem)),iErr)

  if(st_debug) then
    call sixin_echoVal("name",    bez(iElem),                                       "CHEBY",iLine)
    call sixin_echoVal("filename",trim(cheby_filename(cheby_itable(icheby(iElem)))),"CHEBY",iLine)
    if(nSplit >= 3) call sixin_echoVal("r2",      cheby_r2(icheby(iElem)),          "CHEBY",iLine)
    if(nSplit >= 4) call sixin_echoVal("r1",      cheby_r1(icheby(iElem)),          "CHEBY",iLine)
    if(nSplit >= 5) call sixin_echoVal("tilt",    cheby_angle(icheby(iElem)),       "CHEBY",iLine)
    if(nSplit >= 6) call sixin_echoVal("offset_x",cheby_offset_x(icheby(iElem)),    "CHEBY",iLine)
    if(nSplit >= 7) call sixin_echoVal("offset_y",cheby_offset_y(icheby(iElem)),    "CHEBY",iLine)
    if(nSplit >= 8) call sixin_echoVal("L"       ,cheby_L(icheby(iElem)),           "CHEBY",iLine)
    if(nSplit >= 9) call sixin_echoVal("I"       ,cheby_I(icheby(iElem)),           "CHEBY",iLine)
  end if

end subroutine cheby_parseInputLine


subroutine cheby_parseInputDone(iErr)

  use mod_common, only : kz,bez

  implicit none

  logical, intent(inout) :: iErr

  integer j

  ! Loop over single elements to check that they have been defined in the fort.3 block
  if(mcheby /= 0) then
    do j=1,nele
      if(kz(j) == kzcheby) then
        if( icheby(j) == 0) then
          write(lout,"(a)") "CHEBY> ERROR Chebyshev lens element '"//trim(bez(j))//"'not defined in fort.3."
          write(lout,"(a)") "CHEBY>       You must define every Chebyshev lens in the CHEB block."
          iErr = .true.
          return
        end if
      end if
    end do
  end if

end subroutine cheby_parseInputDone


subroutine cheby_postInput
  use utils
  use mod_common, only : kz,bez
  use mod_settings, only : st_quiet

  integer jj, kk
  logical exist
  
  ! Parse files with coefficients for Chebyshev polynomials:
   do jj=1,mcheby_tables
    inquire(file=cheby_filename(jj), exist=exist)
    if(.not. exist) then
      write(lout,"(a)") "CHEBY> ERROR Problems with file with coefficients for Chebyshev polynominals: ", &
            trim(cheby_filename(jj))
      call prror(-1)
    end if
    call parseChebyFile(jj)
  end do

  ! finalise setting-up of chebyshev lenses
  do jj=1,mcheby
    ! some checks and further post-processing of declared lines
    if (cheby_r2(jj)<=zero) cheby_r2(jj)=cheby_refR(cheby_itable(jj))
    if (cheby_r1(jj)<=zero) cheby_r1(jj)=zero
    if (cheby_L (jj)<=zero) then
      cheby_L (jj)=cheby_refL(cheby_itable(jj))
    else
      cheby_scalingFact(jj)=cheby_scalingFact(jj)*((cheby_L(jj))/(cheby_refL(cheby_itable(jj))))
    end if
    if (cheby_I (jj)<=zero) then
      cheby_I (jj)=cheby_refI(cheby_itable(jj))
    else
      cheby_scalingFact(jj)=cheby_scalingFact(jj)*((cheby_I(jj))/(cheby_refI(cheby_itable(jj))))
    end if
    if(st_quiet < 2) then
      write(lout,"(a)") ''
      write(lout,"(a,i0)") ' status of chebyshev lens #',jj
      do kk=1,nele
        if(kz(kk) == kzcheby .and. icheby(kk) == mcheby ) then
          write(lout,"(a,i0)") ' name:',trim(bez(kk))
          exit
        end if
      end do
      write(lout,"(a)")        " - filename         :",trim(cheby_filename(cheby_itable(jj)))
      write(lout,"(a,e22.15)") " - R2           [mm]:",cheby_r2(jj)
      write(lout,"(a,e22.15)") " - R1           [mm]:",cheby_r1(jj)
      write(lout,"(a,e22.15)") " - tilt angle  [rad]:",cheby_angle(jj)
      write(lout,"(a,e22.15)") " - hor offset   [mm]:",cheby_offset_x(jj)
      write(lout,"(a,e22.15)") " - ver offset   [mm]:",cheby_offset_y(jj)
      write(lout,"(a,e22.15)") " - lens length   [m]:",cheby_L(jj)
      write(lout,"(a,e22.15)") " - lens powering [A]:",cheby_I(jj)
    end if
  end do

end subroutine cheby_postInput


! ================================================================================================ !
!  Last modified: 2019-02-25
!  Rewritten by VKBO, June 2018
!  Read file with coefficients for chebyshev polynomials
!  ifile is index of file in table of chebyshev files
!  file is structured as:
!    keyword : value
!  keywords:
!  - I: reference current intensity [A] (optional)
!  - L: reference length of thick lens [m] (optional)
!  - R: reference radius [mm] (mandatory)
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
  use mod_common
  use mod_settings
  use string_tools
  use mod_units

  implicit none

  integer, intent(in) :: ifile

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mInputLn) inLine
  integer nSplit
  logical spErr, err, lDefI, lDefL

  integer ii, jj, fUnit

  lDefI=.true.
  lDefL=.true.

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
    ! Read reference current of lens
    if(nSplit < 3) then
      write(lout,"(a)") "CHEBY> ERROR Not enough arguments for expressing ref lens current [A]."
      write(lout,"(a)") "CHEBY>       Correct format:"
      write(lout,"(a)") "I : value"
      goto 30
    end if
    call chr_cast(lnSplit(3),cheby_refI(ifile),spErr)
    lDefI=.false.

  else if(inLine(1:1) == "L") then
    ! Read reference length of lens
    if(nSplit < 3) then
      write(lout,"(a)") "CHEBY> ERROR Not enough arguments for expressing ref lens length [m]."
      write(lout,"(a)") "CHEBY>       Correct format:"
      write(lout,"(a)") "L : value"
      goto 30
    end if
    call chr_cast(lnSplit(3),cheby_refL(ifile),spErr)
    lDefL=.false.

  else if(inLine(1:1) == "R") then
    ! Read reference radius e-beam in e-lens
    if(nSplit < 3) then
      write(lout,"(a)") "CHEBY> ERROR Not enough arguments for expressing ref lens radius [mm]."
      write(lout,"(a)") "CHEBY>       Correct format:"
      write(lout,"(a)") "R : value"
      goto 30
    end if
    call chr_cast(lnSplit(3),cheby_refR(ifile),spErr)

  else
    ! Read chebyshev coefficients
    if(nSplit /= 4) then
      goto 30
    end if
    call chr_cast(lnSplit(1),ii,spErr)
    if(ii > cheby_max_order) then
      goto 30
    end if
    call chr_cast(lnSplit(2),jj,spErr)
    if(jj > cheby_max_order) then
      goto 30
    end if
    call chr_cast(lnSplit(4),cheby_coeffs(ii,jj,ifile),spErr)

  end if ! close if for keyword identification
  goto 10

20 continue

  call f_close(fUnit)
  if (.not.lDefL.and.cheby_refL(ifile)<=zero) then
    write(lout,"(a)") "CHEBY> ERROR ref lens length [m] must be positive."
    goto 30
  end if
  if (cheby_refR(ifile)<=zero) then
    write(lout,"(a)") "CHEBY> ERROR ref lens radius [mm] must be positive."
    goto 30
  end if

  if(st_quiet < 2) then
    ! Echo parsed data (unless told to be quiet!)
    write(lout,"(a)" ) ""
    write(lout,"(a,i0)") "CHEBY> Coefficients for Chebyshev polynomials as from file "//&
      trim(cheby_filename(ifile))//" - #",ifile
    write(lout,"(a,e22.15)") "CHEBY> * Reference current [A] : ",cheby_refI(ifile)
    if (lDefI) write(lout,"(a)") "         --> default value!"
    write(lout,"(a,e22.15)") "CHEBY> * reference length  [m] : ",cheby_refL(ifile)
    if (lDefL) write(lout,"(a)") "         --> default value!"
    write(lout,"(a,e22.15)") "CHEBY> * reference radius [mm] : ",cheby_refR(ifile)
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
  write(lout,"(a)") "CHEBY> ERROR while parsing file "//trim(cheby_filename(ifile))
  call prror(-1)

end subroutine parseChebyFile


subroutine cheby_kick(jcheby)

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 20-02-2018
  ! apply kick of electron lens

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
    if (rr.gt.cheby_refR(cheby_itable(jcheby))) then
      write(lout,"(a,3(e12.6,1x),a,e12.6)") "CHEBY> ERROR in cheby_kick: particle at position (x,y,r): ",&
        xv1(j), xv2(j), rr,' is outside radial domain of Chebyshev polinomials: ',cheby_refR(cheby_itable(jcheby))
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
  uu=xx/cheby_refR(iTable)
  vv=yy/cheby_refR(iTable)
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
  dxp=-dxp/cheby_refR(iTable)
  dyp=-dyp/cheby_refR(iTable)

  ! take into account Brho and beta
  dxp=dxp/(brho_b*clight*beta_b)
  dyp=dyp/(brho_b*clight*beta_b)

 end subroutine cheby_getKick

end module cheby
