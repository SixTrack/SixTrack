! M. Fitterer, FNAL, A. Mereghtti, CERN
! last modified: 09-02-2018
! Common block for electron lens definition
module elens

  use parpro
  use floatPrecision
  use crcoall
  use mod_alloc

  implicit none

  ! size of table with elens data
  integer, parameter     :: nelens=125
  ! last elens read
  integer, save          :: melens

  ! index of elens:
  integer,allocatable, save          :: ielens(:) !(nele)

  ! variables to save elens parameters for tracking etc.
  integer, save          :: elens_type(nelens)      ! integer for elens type
                                                    ! 0 : Un-initialized.
                                                    ! 1 : Hollow annular elens, uniform profile
  real(kind=fPrec), save :: elens_theta_r2(nelens)    ! kick strength at R2 [mrad]
  real(kind=fPrec), save :: elens_r2(nelens)          ! outer radius R2 [mm]
  real(kind=fPrec), save :: elens_r1(nelens)          ! inner radius R1 [mm]
  real(kind=fPrec), save :: elens_offset_x(nelens), elens_offset_y(nelens)  ! hor./vert. offset of elens [mm]
  real(kind=fPrec), save :: elens_sig(nelens)         ! sig (Gaussian profile) [mm]
  real(kind=fPrec), save :: elens_geo_norm(nelens)    ! normalisation of f(r)
  real(kind=fPrec), save :: elens_len(nelens)         ! length of eLens (e-beam region) [m]
  real(kind=fPrec), save :: elens_I(nelens)           ! current of e-beam [A]
                                                      ! <0: e-beam opposite to beam
  real(kind=fPrec), save :: elens_Ek(nelens)          ! kinetic energy of e-beam [keV]
  logical, save          :: elens_lThetaR2(nelens)    ! flag for computing theta@R2
  integer, save          :: elens_iCheby(nelens)      ! mapping to the table with chebyshev coeffs
  real(kind=fPrec), save :: elens_cheby_angle(nelens) ! angle for getting the real bends [deg]
  ! file with chebyshev coefficients
  integer, parameter     :: nelens_cheby_tables=20    ! number of tables with chebyshev coefficients
  integer, parameter     :: elens_cheby_unit=107      ! unit for reading the chebyshev coefficients
  integer, parameter     :: elens_cheby_order=18      ! max order of chebyshev polynomials
  integer, save          :: melens_cheby_tables       ! tables available in memory
  character(len=16), save:: elens_cheby_filename(nelens_cheby_tables) ! names
  real(kind=fPrec), save :: elens_cheby_coeffs(0:elens_cheby_order,0:elens_cheby_order,nelens_cheby_tables)
  real(kind=fPrec), save :: elens_cheby_refCurr(nelens_cheby_tables) ! reference current [A]
  real(kind=fPrec), save :: elens_cheby_refRadius(nelens_cheby_tables) ! reference radius [mm]
  real(kind=fPrec), save :: elens_cheby_refBeta(nelens_cheby_tables) ! reference e-beta []

contains

subroutine elens_allocate_arrays
  use crcoall
  implicit none
  integer stat
    call alloc(ielens,nele,0,'ielens')
end subroutine elens_allocate_arrays

subroutine elens_expand_arrays(nele_new)
  implicit none
  integer, intent(in) :: nele_new
  call alloc(ielens,nele_new,0,'ielens')
end subroutine elens_expand_arrays

! ================================================================================================ !
!  Parse Elens Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-25
! ================================================================================================ !
subroutine elens_parseInputLine(inLine, iLine, iErr)

  use mod_settings
  use sixtrack_input
  use mathlib_bouncer
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mStrLen) tmpch
  real(kind=fPrec) tmpflt
  integer nSplit, iElem, j, chIdx, tmpi1, tmpi2, tmpi3
  logical spErr, tmpl

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "ELENS> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit < 7) then
    write(lout,"(a,i0)") "ELENS> ERROR Expected at least 7 input parameters, got ",nSplit
    iErr = .true.
    return
  end if

  iElem = -1
  do j=1,nele
    if(bez(j) == trim(lnSplit(1))) then
      iElem = j
      exit
    end if
  end do
  if(iElem == -1) then
    write(lout,"(a)") "ELENS> ERROR Element '"//trim(lnSplit(1))//"' not found in single element list."
    iErr = .true.
    return
  end if

  if(kz(iElem) /= 29) then
    write(lout,"(2(a,i0),a)") "ELENS> ERROR Element type is kz(",iElem,") = ",kz(iElem)," != 29"
    iErr = .true.
    return
  end if
  if(el(iElem) /= 0 .or. ek(iElem) /= 0 .or. ed(iElem) /= 0) then
    write(lout,"(a)")       "ELENS> ERROR Length el(iElem) (elens is treated as thin element), "//&
      "and first and second field have to be zero:"
    write(lout,"(2(a,i0))") "ELENS>       el(",iElem,") = ",el(iElem)," != 0"
    write(lout,"(2(a,i0))") "ELENS>       ed(",iElem,") = ",ed(iElem)," != 0"
    write(lout,"(2(a,i0))") "ELENS>       ek(",iElem,") = ",ek(iElem)," != 0"
    iErr = .true.
    return
  end if

  melens = melens+1
  if(melens > nelens) then
    write(lout,"(2(a,i0))") "ELENS> ERROR Too many elenses: ",melens,". Max is ",nelens
    iErr = .true.
    return
  end if

  ielens(iElem) = melens
  if(elens_type(ielens(iElem)) /= 0) then
    write(lout,"(a)") "ELENS> ERROR The element '"//trim(bez(iElem))//"' was defined twice."
    iErr = .true.
    return
  end if

  ! Parse the element
  select case (lnSplit(2))
  case("UNIFORM")
    elens_type(ielens(iElem)) = 1
  case("GAUSSIAN")
    elens_type(ielens(iElem)) = 2
    if(nSplit < 8) then
      write(lout,"(a,i0)") "ELENS> ERROR Expected at least 8 input parameters for GAUSSIAN, got ",nSplit
      iErr = .true.
      return
    end if
  case("CHEBYSHEV")
    write(lout,"(a)") "ELENS> ERROR CHEBYSHEV type not fully supported yet - elens name: '"//trim(bez(iElem))
    iErr = .true.
    return
!     elens_type(ielens(iElem)) = 3
!     if(nSplit < 8) then
!       write(lout,"(a,i0)") "ELENS> ERROR Expected at least 8 input parameters for CHEBYSHEV, got ",nSplit
!       iErr = .true.
!       return
!     end if
  case default
    write(lout,"(a)") "ELENS> ERROR Elens type '"//trim(lnSplit(2))//"' not recognized. Remember to use all UPPER CASE."
    iErr = .true.
    return
  end select

  call chr_cast(lnSplit(3),elens_theta_r2(ielens(iElem)),iErr)
  call chr_cast(lnSplit(4),elens_r2(ielens(iElem)),      iErr)
  call chr_cast(lnSplit(5),elens_r1(ielens(iElem)),      iErr)
  call chr_cast(lnSplit(6),elens_offset_x(ielens(iElem)),iErr)
  call chr_cast(lnSplit(7),elens_offset_y(ielens(iElem)),iErr)

  if(elens_type(ielens(iElem)) == 2) then
    ! GAUSSIAN profile of electrons: need also sigma of e-beam
    call chr_cast(lnSplit(8),elens_sig(ielens(iElem)),iErr)

!   elseif(elens_type(ielens(iElem)) == 3 )then
!     ! Profile of electrons given by Chebyshev polynomials: need also
!     !   name of file where coefficients are stored and angle
!     tmpch = trim(lnSplit(8))
!     call chr_cast(lnSplit(9),elens_cheby_angle(ielens(iElem)),iErr)
! 
!     ! Check if table of coefficients has already been requested:
!     chIdx = -1
!     do tmpi1=1,melens_cheby_tables
!       if(tmpch == elens_cheby_filename(tmpi1)) then
!         elens_iCheby(ielens(iElem)) = tmpi1
!         chIdx = tmpi1
!         exit
!       end if
!     end do
!     if(chIdx == -1) then
!       ! Unsuccessful search
!       melens_cheby_tables = melens_cheby_tables+1
!       if(melens_cheby_tables > nelens_cheby_tables) then
!         write(lout,"(2(a,i0))") "ELENS> ERROR Too many tables for Chebyshev coefficients: ",melens_cheby_tables,&
!           ". Max is ",nelens_cheby_tables
!         iErr = .true.
!         return
!       end if
!       elens_iCheby(ielens(iElem)) = melens_cheby_tables
!       elens_cheby_filename(tmpi1) = tmpch
!     end if
  end if

  ! Additional geometrical infos:
  ! Depending on profile, the position of these parameters change
  tmpi1 = 0
  tmpi2 = 0
  tmpi3 = 0
  if(elens_type(ielens(iElem)) == 1 .and. nSplit >= 10) then
    tmpi1 = 8
    tmpi2 = 9
    tmpi3 = 10
    elens_lThetaR2(ielens(iElem)) = .true.
  else if(elens_type(ielens(iElem)) == 2 .and. nSplit >= 11) then
    tmpi1 = 9
    tmpi2 = 10
    tmpi3 = 11
    elens_lThetaR2(ielens(iElem)) = .true.
!   else if(elens_type(ielens(iElem)) == 3 .and. nSplit >= 12) then
!     tmpi1 = 10
!     tmpi2 = 11
!     tmpi3 = 12
!     elens_lThetaR2(ielens(iElem)) = .true.
  end if

  if(elens_lThetaR2(ielens(iElem))) then
    call chr_cast(lnSplit(tmpi1),elens_len(ielens(iElem)),iErr)
    call chr_cast(lnSplit(tmpi2),elens_I(ielens(iElem)),  iErr)
    call chr_cast(lnSplit(tmpi3),elens_Ek(ielens(iElem)), iErr)
  end if

  ! sanity checks
  if(elens_r2(ielens(iElem)) < elens_r1(ielens(iElem))) then
    write(lout,"(a)") "ELENS> WARNING ELEN R2<R1. Inverting."
    tmpflt=elens_r2(ielens(iElem))
    elens_r2(ielens(iElem)) = elens_r1(ielens(iElem))
    elens_r1(ielens(iElem)) = tmpflt
  else if(elens_r2(ielens(iElem)) == elens_r1(ielens(iElem))) then
    write(lout,"(a)") "ELENS> ERROR ELEN R2=R1. Elens does not exist."
    iErr = .true.
    return
  end if
  if(elens_r2(ielens(iElem)) <= zero) then
    write(lout,"(a)") "ELENS> ERROR R2<=0!"
    iErr = .true.
    return
  end if
  if(elens_r1(ielens(iElem)) <= zero) then
    write(lout,"(a)") "ELENS> ERROR R1<=0!"
    iErr = .true.
    return
  end if
  if(elens_lThetaR2(ielens(iElem))) then
    if(elens_len(ielens(iElem)) <= zero) then
      write(lout,"(a)") "ELENS> ERROR L<0!"
      iErr = .true.
      return
    end if
    if(elens_I(ielens(iElem)) == zero) then
      write(lout,"(a)") "ELENS> ERROR I=0!"
      iErr = .true.
      return
    end if
    if(elens_Ek(ielens(iElem)) <= zero) then
      write(lout,"(a)") "ELENS> ERROR Ek<0! (e-beam)"
      iErr = .true.
      return
    end if
  end if
  if( elens_type(ielens(iElem)) == 2 ) then
    if ( elens_sig(ielens(iElem)).le.zero ) then
       write(lout,"(a)") "ELENS> ERROR sigma of electron beam <=0 in Elens '"//trim(bez(iElem))//"'."
       iErr = .true.
       return
    end if
  end if

  ! Proper normalisation
  if(elens_type(ielens(iElem)) == 1) then
    ! Uniform distribution
    elens_geo_norm(ielens(iElem)) = (elens_r2(ielens(iElem))+elens_r1(ielens(iElem)))*&
                                    (elens_r2(ielens(iElem))-elens_r1(ielens(iElem)))
  else if(elens_type(ielens(iElem)) == 2) then
    ! Gaussian distribution
    elens_geo_norm(ielens(iElem)) = exp_mb(-0.5*(elens_r1(ielens(iElem))/elens_sig(ielens(iElem)))**2)&
                                   -exp_mb(-0.5*(elens_r2(ielens(iElem))/elens_sig(ielens(iElem)))**2)
  end if

  if(st_debug) then
    call sixin_echoVal("name",    bez(iElem),                   "ELENS",iLine)
    call sixin_echoVal("type",    lnSplit(2),                   "ELENS",iLine)
    call sixin_echoVal("type",    elens_type(ielens(iElem)),    "ELENS",iLine)
    call sixin_echoVal("theta_r2",elens_theta_r2(ielens(iElem)),"ELENS",iLine)
    call sixin_echoVal("r1",      elens_r1(ielens(iElem)),      "ELENS",iLine)
    call sixin_echoVal("r2",      elens_r2(ielens(iElem)),      "ELENS",iLine)
    call sixin_echoVal("offset_x",elens_offset_x(ielens(iElem)),"ELENS",iLine)
    call sixin_echoVal("offset_y",elens_offset_y(ielens(iElem)),"ELENS",iLine)
    if(elens_type(ielens(iElem)) == 2) then
      call sixin_echoVal("sig",elens_sig(ielens(iElem)),"ELENS",iLine)
    end if
    if(elens_lThetaR2(ielens(iElem))) then
      call sixin_echoVal("L", elens_len(ielens(iElem)),"ELENS",iLine)
      call sixin_echoVal("I", elens_I(ielens(iElem)),  "ELENS",iLine)
      call sixin_echoVal("Ek",elens_Ek(ielens(iElem)), "ELENS",iLine)
    end if
  end if

end subroutine elens_parseInputLine

subroutine elens_parseInputDone(iErr)

  use mod_common, only : kz,bez

  implicit none

  logical, intent(inout) :: iErr

  integer j

  ! Loop over single elements to check that they have been defined in the fort.3 block
  if(melens /= 0) then
    do j=1,nele
      if(kz(j) == 29) then
        if(elens_type(ielens(j)) == 0) then
          write(lout,"(a)") "ELENS> ERROR Elens element '"//trim(bez(j))//"'not defined in fort.3."
          write(lout,"(a)") "ELENS>       You must define every wire in the WIRE block."
          iErr = .true.
          return
        end if
      end if
    end do
  end if

end subroutine elens_parseInputDone

subroutine elens_postInput
  
  use mod_common, only : e0,bez,kz
  use mod_hions, only : aa0, zz0
  
  integer j,jj
  logical exist

  ! Compute elens theta at R2, if requested by user
  do j=1,melens
    if(elens_lThetaR2(j)) then
      do jj=1,nele
        if(kz(jj)==29) then
          if (ielens(jj).eq.j) then
            exit
          end if
        end if
      end do
      call eLensTheta( j, e0 )
      write(lout,"(a,i0,a,e22.15)") "ELENS> New theta at r2 for elens #",j," named "//trim(bez(jj))//": ",elens_theta_r2(j)
    end if
  end do

  ! Parse files with coefficients for Chebyshev polynomials:
   do j=1,melens_cheby_tables
    inquire(file=elens_cheby_filename(j), exist=exist)
    if(.not. exist) then
      write(lout,"(a)") "ELENS> Problems with file with coefficients for Chebyshev polynominals: ",trim(elens_cheby_filename(j))
      call prror(-1)
    end if
    call parseChebyFile(j)
  end do

end subroutine elens_postInput

! ================================================================================================ !
!  Compute eLens theta at r2
!  input variables:
!  - length of eLens [m];
!  - current intensity of e-beam [A]
!  - kinetic energy of electrons [keV]
!  - total beam energy [MeV]
!  - outer radius [mm]
! ================================================================================================ !
subroutine eLensTheta( j, Etot )

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use physical_constants
  use crcoall
  use mod_common
  use mod_hions, only : nucm0, zz0

  implicit none

  integer j
  real(kind=fPrec) gamma, beta_e, beta_b, brho, Etot

  gamma  = ((elens_Ek(j)*c1m3)/pmae)+one ! from kinetic energy
  beta_e = sqrt((gamma+one)*(gamma-one))/(gamma)
  gamma  = Etot/nucm0                ! from total energy
  beta_b = sqrt((gamma+one)*(gamma-one))/(gamma)
  brho   = (Etot/(clight*c1m6))/zz0

  ! r2: from mm to m (c1m3)
  ! theta: from rad to mrad (c1e3)
  elens_theta_r2(j) = ((elens_len(j)*abs(elens_I(j)))/((((two*pi)*((eps0*clight)*clight))*brho)*(elens_r2(j)*c1m3)))*c1e3
  if(elens_I(j) < zero) then
    elens_theta_r2(j) = elens_theta_r2(j)*(one/(beta_e*beta_b)+one)
  else
    elens_theta_r2(j) = elens_theta_r2(j)*(one/(beta_e*beta_b)-one)
  end if

  if ( elens_type(j) == 2 ) then
     elens_theta_r2(j) = elens_theta_r2(j) * elens_geo_norm(j)
  end if

end subroutine eLensTheta

! ================================================================================================ !
!  Last modified: 2018-06-25
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

  implicit none

  integer, intent(in) :: ifile

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mInputLn) inLine
  integer nSplit
  logical spErr

  integer iErr, ii, jj
  real(kind=fPrec) tmpflt, beta, gamma

  ierr = 0
  write(lout,"(a)") "ELENS> Parsing file with coefficients for Chebyshev polynomials "//trim(elens_cheby_filename(ifile))
  open(elens_cheby_unit,file=elens_cheby_filename(ifile),status="old")

10 continue
  read(elens_cheby_unit,"(a)",end=20,err=30) inLine
  if(inLine(1:1) == "#") goto 10

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "ELENS> ERROR Failed to Chebyshev input line."
    goto 30
  end if

  if(inLine(1:1) == "I") then
    ! Read reference current of e-beam in e-lens
    if(nSplit < 3) then
      iErr = 2
      goto 30
    end if
    call chr_cast(lnSplit(3),tmpflt,spErr)
    elens_cheby_refCurr(ifile) = tmpflt

  else if(inLine(1:2) == "Ek") then
    ! Read reference kinetic energy of e-beam in e-lens
    if(nSplit < 3) then
      iErr = 3
      goto 30
    end if
    call chr_cast(lnSplit(3),tmpflt,spErr)
    gamma = (tmpflt*c1m3)/pmae+one ! from kinetic energy
    elens_cheby_refBeta(ifile) = sqrt((gamma+one)*(gamma-one))/(gamma)

  else if(inLine(1:3) == "rad") then
    ! Read reference radius e-beam in e-lens
    if(nSplit < 3) then
      iErr = 4
      goto 30
    end if
    call chr_cast(lnSplit(3),tmpflt,spErr)
    elens_cheby_refRadius(ifile) = tmpflt

  else
    ! Read chebyshev coefficients
    if(nSplit /= 4) then
      iErr = 5
      goto 30
    end if
    call chr_cast(lnSplit(1),ii,spErr)
    if(ii > elens_cheby_order) then
      iErr = 6
      goto 30
    end if
    call chr_cast(lnSplit(2),jj,spErr)
    if(jj > elens_cheby_order) then
      iErr = 7
      goto 30
    end if
    call chr_cast(lnSplit(4),tmpflt,spErr)
    elens_cheby_coeffs(ii,jj,ifile) = tmpflt

  end if ! close if for keyword identification
  goto 10

20 continue

  close(elens_cheby_unit)

  if(st_quiet < 2) then
    ! Echo parsed data (unless told to be quiet!)
    write(lout,"(a,i0)") "ELENS> Coefficients for Chebyshev polynomials as from file "//&
      trim(elens_cheby_filename(ifile))//" - #",ifile
    write(lout,"(a,e22.15)") "ELENS> * Reference current [A] : ",elens_cheby_refCurr(ifile)
    write(lout,"(a,e22.15)") "ELENS> * reference beta     [] : ",elens_cheby_refBeta(ifile)
    write(lout,"(a,e22.15)") "ELENS> * reference radius [mm] : ",elens_cheby_refRadius(ifile)
    do ii=0,elens_cheby_order
      do jj=0,elens_cheby_order
        if(elens_cheby_coeffs(ii,jj,ifile)/= zero) then
          write(lout,"(2(a,i4),a,e22.15)") "ELENS> Order ",ii,",",jj," : ",elens_cheby_coeffs(ii,jj,ifile)
        end if
      end do
    end do
  end if
  return

30 continue
  write(lout,"(a,i0,a)") "ELENS> ERROR ",iErr," while parsing file "//trim(elens_cheby_filename(ifile))
  call prror(-1)

end subroutine parseChebyFile

end module elens
