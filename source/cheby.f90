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
  integer, parameter          :: cheby_kz=42          ! kz of chebyshev lenses
  integer, parameter          :: cheby_ktrack=67      ! ktrack of chebyshev lenses

  ! variables to save parameters for tracking etc.
  integer, save          :: cheby_itable(ncheby)=0         ! index of chebyshev table
  real(kind=fPrec), save :: cheby_r2(ncheby)     = zero    ! outer radius R2 [mm] (optional)
  real(kind=fPrec), save :: cheby_r1(ncheby)     = zero    ! inner radius R1 [mm] (optional)
  real(kind=fPrec), save :: cheby_angle(ncheby)=zero       ! rotation angle about the longitudinal axis [deg] (optional)
  real(kind=fPrec), save :: cheby_offset_x(ncheby)=zero    ! hor. offset [mm] (optional)
  real(kind=fPrec), save :: cheby_offset_y(ncheby)=zero    ! ver. offset [mm] (optional)
  real(kind=fPrec), save :: cheby_I(ncheby)=-one           ! actual powering of lens [A] (optional)
  real(kind=fPrec), save :: cheby_scalingFact(ncheby)=one  ! scaling factor [] (computed internally)

  ! tables with chebyshev coefficients
  integer, parameter     :: cheby_max_order=30             ! max order of chebyshev polynomials currently supported
  character(len=mFNameLen), save:: cheby_filename(ncheby_tables) = " "! file names
  real(kind=fPrec), save :: cheby_coeffs(0:cheby_max_order,0:cheby_max_order,ncheby_tables) = zero ! coefficients
  integer, save          :: cheby_maxOrder(ncheby_tables)  = 0    ! max order of the current map
  real(kind=fPrec), save :: cheby_refI(ncheby_tables) = one  ! reference current [A] (optional)
  real(kind=fPrec), save :: cheby_refR(ncheby_tables) = zero ! reference radius [mm] (mandatory)

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

  if(kz(iElem) /= cheby_kz) then
    write(lout,"(3(a,i0))") "CHEBY> ERROR Element type is kz(",iElem,") = ",kz(iElem)," != ",cheby_kz
    iErr = .true.
    return
  end if
  if(el(iElem) /= zero .or. ek(iElem) /= zero .or. ed(iElem) /= zero) then
    write(lout,"(a)")       "CHEBY> ERROR Length el(iElem) (Chebyshev lens is treated as thin element), "//&
      "and first and second field have to be zero:"
    write(lout,"(2(a,i0),a)") "CHEBY>       el(",iElem,") = ",el(iElem)," != 0"
    write(lout,"(2(a,i0),a)") "CHEBY>       ed(",iElem,") = ",ed(iElem)," != 0"
    write(lout,"(2(a,i0),a)") "CHEBY>       ek(",iElem,") = ",ek(iElem)," != 0"
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
  if(nSplit >= 8) call chr_cast(lnSplit(8),cheby_I(icheby(iElem)),iErr)

  if(st_debug) then
    call sixin_echoVal("name",    bez(iElem),                                         "CHEBY",iLine)
    call sixin_echoVal("filename",trim(cheby_filename(cheby_itable(icheby(iElem)))),  "CHEBY",iLine)
    if(nSplit >= 3) call sixin_echoVal("r2 [mm]"      , cheby_r2(icheby(iElem)),      "CHEBY",iLine)
    if(nSplit >= 4) call sixin_echoVal("r1 [mm]"      , cheby_r1(icheby(iElem)),      "CHEBY",iLine)
    if(nSplit >= 5) call sixin_echoVal("tilt [deg]"   , cheby_angle(icheby(iElem)),   "CHEBY",iLine)
    if(nSplit >= 6) call sixin_echoVal("offset_x [mm]", cheby_offset_x(icheby(iElem)),"CHEBY",iLine)
    if(nSplit >= 7) call sixin_echoVal("offset_y [mm]", cheby_offset_y(icheby(iElem)),"CHEBY",iLine)
    if(nSplit >= 8) call sixin_echoVal("I [A]"        , cheby_I(icheby(iElem)),       "CHEBY",iLine)
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
      if(kz(j) == cheby_kz) then
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
!  use utils
  use mod_common, only : kz,bez
  use mod_settings, only : st_quiet

  integer jj, kk
  logical exist
  real(kind=fPrec) tmpFlt
  
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
    if (cheby_r2(jj)<cheby_r1(jj)) then
      ! swap
      tmpFlt=cheby_r2(jj)
      cheby_r2(jj)=cheby_r1(jj)
      cheby_r1(jj)=tmpFlt
    end if
    if (cheby_r2(jj)>cheby_refR(cheby_itable(jj))) then
      write(lout,"(a)") "CHEBY> ERROR R2 cannot be larger than domain of Chebyshev polynomials!"
      write(lout,"(a,1pe22.15,a,1pe22.15)") "CHEBY>       R2 [mm]: ",cheby_r2(jj), &
           " - reference radius [mm]:",cheby_refR(cheby_itable(jj))
      call prror(-1)
    end if
    if (cheby_I (jj)<=zero) then
      cheby_I (jj)=cheby_refI(cheby_itable(jj))
    else
      cheby_scalingFact(jj)=cheby_scalingFact(jj)*((cheby_I(jj))/(cheby_refI(cheby_itable(jj))))
    end if
    if(st_quiet < 2) then
      ! find name
      do kk=1,nele
        if(kz(kk) == cheby_kz .and. icheby(kk) == jj ) then
          exit
        end if
      end do
      write(lout,"(a)") ''
      write(lout,"(a,i0,a)") ' status of chebyshev lens #',jj," - name: '"//trim(bez(kk))//"'"
      write(lout,"(a)")          " - filename         : '"//trim(cheby_filename(cheby_itable(jj)))//"'"
      write(lout,"(a,1pe22.15)") " - R2           [mm]: ",cheby_r2(jj)
      write(lout,"(a,1pe22.15)") " - R1           [mm]: ",cheby_r1(jj)
      write(lout,"(a,1pe22.15)") " - tilt angle  [deg]: ",cheby_angle(jj)
      write(lout,"(a,1pe22.15)") " - hor offset   [mm]: ",cheby_offset_x(jj)
      write(lout,"(a,1pe22.15)") " - ver offset   [mm]: ",cheby_offset_y(jj)
      write(lout,"(a,1pe22.15)") " - lens powering [A]: ",cheby_I(jj)
      write(lout,"(a,1pe22.15)") " - scaling factor []: ",cheby_scalingFact(jj)
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

  use mod_common
  use mod_settings
  use string_tools
  use mod_units

  implicit none

  integer, intent(in) :: ifile

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mInputLn) inLine
  integer nSplit, ii, jj, fUnit
  logical spErr, err, lDefI

  lDefI=.true.

  write(lout,"(a)") "CHEBY> Parsing file with coefficients for Chebyshev polynomials "//trim(cheby_filename(ifile))
  call f_requestUnit(cheby_filename(ifile),fUnit)
  call f_open(unit=fUnit,file=cheby_filename(ifile),mode='r',err=err,formatted=.true.,status="old")
  if(err) then
    write(lout,"(a)") "CHEBY> ERROR Failed to open file."
    goto 40
  end if

10 continue
  read(fUnit,"(a)",end=20,err=30) inLine
  if(inLine(1:1) == "#") goto 10

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "CHEBY> ERROR Failed to split Chebyshev input line:"
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
    if(spErr) then
      write(lout,"(a)") "CHEBY> ERROR in casting ref lens current: "//trim(lnSplit(3))
      goto 30
    end if
    lDefI=.false.

  else if(inLine(1:1) == "R") then
    ! Read reference radius e-beam in e-lens
    if(nSplit < 3) then
      write(lout,"(a)") "CHEBY> ERROR Not enough arguments for expressing ref lens radius [mm]."
      write(lout,"(a)") "CHEBY>       Correct format:"
      write(lout,"(a)") "R : value"
      goto 30
    end if
    call chr_cast(lnSplit(3),cheby_refR(ifile),spErr)
    if(spErr) then
      write(lout,"(a)") "CHEBY> ERROR in casting ref lens radius: "//trim(lnSplit(3))
      goto 30
    end if

  else
    ! Read chebyshev coefficients
    if(nSplit /= 4) then
      write(lout,"(a)") "CHEBY> ERROR Not enough arguments for expressing Chebyshev coefficients [Vm]."
      write(lout,"(a)") "CHEBY>       Correct format:"
      write(lout,"(a)") "ii jj : value (ii->x,jj->y)"
      goto 30
    end if
    call chr_cast(lnSplit(1),ii,spErr)
    if(ii > cheby_max_order) then
      write(lout,"(2(a,i0))") "CHEBY> ERROR Too high order in Chebyshev polynomial - requested (hor):", &
            ii," - available:",cheby_max_order
      goto 30
    end if
    call chr_cast(lnSplit(2),jj,spErr)
    if(jj > cheby_max_order) then
      write(lout,"(2(a,i0))") "CHEBY> ERROR Too high order in Chebyshev polynomial - requested (ver):", &
            jj," - available:",cheby_max_order
      goto 30
    end if
    call chr_cast(lnSplit(4),cheby_coeffs(ii,jj,ifile),spErr)
    if(spErr) then
      write(lout,"(a)") "CHEBY> ERROR in casting Chebyshev coefficient: "//trim(lnSplit(4))
      goto 30
    end if
    cheby_maxOrder(ifile)=max(ii,jj,cheby_maxOrder(ifile))

  end if ! close if for keyword identification
  goto 10

20 continue

  call f_close(fUnit)
  if (cheby_refR(ifile)<=zero) then
    write(lout,"(a)") "CHEBY> ERROR ref lens radius [mm] must be positive."
    goto 30
  end if
  if (cheby_maxOrder(ifile)<2) then
    write(lout,"(a,i0,a)") "CHEBY> ERROR max order too low:",cheby_maxOrder(ifile)," - it should be at least 2."
    goto 30
  end if

  if(st_quiet < 2) then
    ! Echo parsed data (unless told to be quiet!)
    write(lout,"(a,i0)") "CHEBY> Coefficients for Chebyshev polynomials as from file '"//&
      trim(cheby_filename(ifile))//"' - #",ifile
    write(lout,"(a,1pe22.15)") "CHEBY> * Reference current [A] : ",cheby_refI(ifile)
    if (lDefI) write(lout,"(a)") "         --> default value!"
    write(lout,"(a,1pe22.15)") "CHEBY> * reference radius [mm] : ",cheby_refR(ifile)
    do ii=0,cheby_max_order
      do jj=0,cheby_max_order
        if(cheby_coeffs(ii,jj,ifile)/= zero) then
          write(lout,"(2(a,i4),a,1pe22.15)") "CHEBY> Order ",ii,",",jj," : ",cheby_coeffs(ii,jj,ifile)
        end if
      end do
    end do
    write(lout,"(a)" ) ""
  end if
  return

30 continue
  write(lout,"(a)") "CHEBY>       last line read:"
  write(lout,"(a)") trim(inLine)
  write(lout,"(a)") "CHEBY>       split fields:"
  do ii=1,nSplit
    write(lout,"('CHEBY>       - ',i2,': ',a)") ii,trim(lnSplit(ii))
  end do
40 continue
  write(lout,"(a)") "CHEBY> ERROR while parsing file "//trim(cheby_filename(ifile))
  call prror(-1)

end subroutine parseChebyFile


subroutine cheby_kick(i,ix,n)

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 01-03-2019
  ! apply kick of Chebyshev lenses

  use mod_common, only : betrel, napx, brho
  use mod_hions, only : moidpsv
  use mod_common_main
  use mathlib_bouncer
  use numerical_constants, only : zero, c180e0, pi
  use physical_constants, only: clight

  integer, intent(in) :: i
  integer, intent(in) :: ix
  integer, intent(in) :: n
  
  real(kind=fPrec) xx, yy, rr, frr, dxp, dyp
  real(kind=fPrec) theta, radio, angle_rad
  integer          j
  logical          lrotate

  ! rotation angle
  lrotate = cheby_angle(icheby(ix)).ne.zero
  angle_rad = (cheby_angle(icheby(ix))/c180e0)*pi

  do j=1,napx

    ! apply offset
    xx=xv1(j)-cheby_offset_x(icheby(ix))
    yy=xv2(j)-cheby_offset_y(icheby(ix))

    ! check that particle is within the domain of chebyshev polynomials
    rr=sqrt(xx**2+yy**2)
    if (rr.ge.cheby_r1(icheby(ix)).and.rr.lt.cheby_r2(icheby(ix))) then ! rr<r1 || rr>=r2 -> no kick from lens
      
      ! in case of non-zero tilt angle, rotate coordinates
      if (lrotate) then
        theta = atan2_mb(yy, xx)-angle_rad
        xx = rr * cos_mb(theta)
        yy = rr * sin_mb(theta)
      end if
      ! compute kick from cheby map
      call cheby_getKick( xx, yy, dxp, dyp, cheby_itable(icheby(ix)) )
      ! in case cheby has a non-zero angle, rotate kicks
      if (lrotate) then
        ! NB: cheby_angle(icheby(ix)) is the rotation angle of the cheby
        theta = atan2_mb(dyp, dxp)+angle_rad
        radio = sqrt(dxp**2 + dyp**2)
        dxp = radio * cos_mb(theta)
        dyp = radio * sin_mb(theta)
      end if
     
      ! take into account scaling factor, Brho of beam and its relativistic beta,
      !    and magnetic rigidity and relativistic beta of particle being tracked
      dxp=((dxp*cheby_scalingFact(icheby(ix)))/(brho*clight*betrel)*moidpsv(j))*rvv(j)
      dyp=((dyp*cheby_scalingFact(icheby(ix)))/(brho*clight*betrel)*moidpsv(j))*rvv(j)
      
      ! apply kicks
      yv1(j)=yv1(j)+dxp
      yv2(j)=yv2(j)+dyp
    end if
  end do

end subroutine cheby_kick


subroutine cheby_getKick( xx, yy, dxp, dyp, iTable )

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 01-03-2019
  ! compute kicks from Chebyshev polinomials - see FermiLAB-FN-0972-APC

  use numerical_constants, only : zero, one, two, c1m3, c1e3

  ! interface vars
  real(kind=fPrec), intent(in ) :: xx
  real(kind=fPrec), intent(in ) :: yy
  real(kind=fPrec), intent(out) :: dxp
  real(kind=fPrec), intent(out) :: dyp
  integer,          intent(in ) :: iTable

  ! temp vars
  real(kind=fPrec) :: uu, vv, Tx (0:cheby_maxOrder(iTable)), Ty (0:cheby_maxOrder(iTable)), &
                      kx, ky, Tpx(0:cheby_maxOrder(iTable)), Tpy(0:cheby_maxOrder(iTable)), &
                      fu, fv
  integer          :: nn, jj

  ! normalised variables
  uu=xx/cheby_refR(iTable)
  vv=yy/cheby_refR(iTable)
  ! normalisation factors of derivatives
  fu=(one-uu)*(one+uu)
  fv=(one-vv)*(one+vv)

  ! - polynomials:
  Tx(0)=one
  Ty(0)=one
  Tx(1)=uu
  Ty(1)=vv
  ! - derivatives:
  Tpx(0)=zero
  Tpy(0)=zero
  Tpx(1)=one
  Tpy(1)=one
  do nn=2,cheby_maxOrder(iTable)
     Tx(nn)=two*(uu*Tx(nn-1))-Tx(nn-2)
     Ty(nn)=two*(vv*Ty(nn-1))-Ty(nn-2)
     Tpx(nn)=(real(nn,fPrec)*(Tx(nn-1)-uu*Tx(nn)))/fu
     Tpy(nn)=(real(nn,fPrec)*(Ty(nn-1)-vv*Ty(nn)))/fv
  end do

  ! get kicks
  dxp=zero
  dyp=zero
  do nn=0,cheby_maxOrder(iTable)
     do jj=0,nn
        dxp=dxp+(cheby_coeffs(jj,nn-jj,iTable)*Tpx(jj))*Ty (nn-jj)
        dyp=dyp+(cheby_coeffs(jj,nn-jj,iTable)*Tx (jj))*Tpy(nn-jj)
     end do
  end do
  dxp=-(dxp*c1e3)/(cheby_refR(iTable)*c1m3) ! ref radius in [mm], kick in [mrad]
  dyp=-(dyp*c1e3)/(cheby_refR(iTable)*c1m3) ! ref radius in [mm], kick in [mrad]

 end subroutine cheby_getKick

end module cheby
