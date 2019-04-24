! M. Fitterer, FNAL, A. Mereghtti, CERN
! last modified: 11-04-2019
! Common block for electron lens definition
module elens

  use parpro
  use floatPrecision
  use crcoall
  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  ! size of table with elens data
  integer, parameter :: nelens = 125
  ! last elens read
  integer, save      :: melens = 0

  integer, parameter          :: elens_kz=29          ! kz of electron lenses
  integer, parameter          :: elens_ktrack=63      ! ktrack of electron lenses
  
  ! index of elens:
  integer,allocatable, save :: ielens(:) !(nele)

  ! variables to save elens parameters for tracking etc.
  integer, save          :: elens_type(nelens)         = 0       ! integer for elens type
                                                                 ! 0 : Un-initialized.
                                                                 ! 1 : uniform profile
                                                                 ! 2 : Gaussian profile
                                                                 ! 3 : radial profile from file
  real(kind=fPrec), save :: elens_theta_r2(nelens)     = zero    ! kick strength at R2 [mrad]
  real(kind=fPrec), save :: elens_r2(nelens)           = zero    ! outer radius R2 [mm]
  real(kind=fPrec), save :: elens_r1(nelens)           = zero    ! inner radius R1 [mm]
  real(kind=fPrec), save :: elens_offset_x(nelens)     = zero    ! hor offset of elens [mm]
  real(kind=fPrec), save :: elens_offset_y(nelens)     = zero    ! vert. offset of elens [mm]
  real(kind=fPrec), save :: elens_sig(nelens)          = zero    ! sig (Gaussian profile) [mm]
  real(kind=fPrec), save :: elens_geo_norm(nelens)     = zero    ! normalisation of f(r)
  real(kind=fPrec), save :: elens_len(nelens)          = zero    ! length of eLens (e-beam region) [m]
  real(kind=fPrec), save :: elens_I(nelens)            = zero    ! current of e-beam [A]
                                                                 ! <0: e-beam opposite to beam
  real(kind=fPrec), save :: elens_Ek(nelens)           = zero    ! kinetic energy of e-beam [keV]
  real(kind=fPrec), save :: elens_beta_e(nelens)                 ! relativistic beta of electrons
  logical, save          :: elens_lThetaR2(nelens)     = .false. ! flag for computing theta@R2
  logical, save          :: elens_lAllowUpdate(nelens) = .true.  ! Flag for disabling updating of kick,
                                                                 !  i.e. after DYNK has touched thetaR2,
                                                                 !  the energy update is disabled.
  logical, save          :: elens_lFox(nelens)         = .false. ! the kick from the elens should be taken into account
                                                                 !  when searching for the closed orbit with DA algebra
#ifdef CR
  logical, save          :: elens_lAllowUpdate_CR(nelens)
#endif
  
  ! radial profile
  integer, save          :: elens_iRadial(nelens)                ! mapping to the radial profile
  real(kind=fPrec), save :: elens_radial_fr1(nelens)             ! value of f(R1) in case of radial profiles from file [0:1]
  real(kind=fPrec), save :: elens_radial_fr2(nelens)             ! value of f(R2) in case of radial profiles from file [0:1]
  integer, save          :: elens_radial_mpoints(nelens)=2       ! how many points for polynomial interpolation
  integer, save          :: elens_radial_jguess(nelens)=-1       ! bin for guessed search
  ! - file handling and data storage:
  integer, parameter     :: nelens_radial_profiles=20 ! max number of radial profiles
  integer, save          :: melens_radial_profiles    ! radial profiles available in memory
  integer, parameter     :: elens_radial_dim=500      ! max number of points in radial profiles
  character(len=mFileName), save:: elens_radial_filename(nelens_radial_profiles) ! names
  real(kind=fPrec), save :: elens_radial_profile_R(0:elens_radial_dim,nelens_radial_profiles) ! [mm]
  real(kind=fPrec), save :: elens_radial_profile_J(0:elens_radial_dim,nelens_radial_profiles) ! [A]
  integer, save          :: elens_radial_profile_nPoints(nelens_radial_profiles)

contains

subroutine elens_allocate_arrays
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
  use string_tools
  use mod_common

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
    write(lerr,"(a)") "ELENS> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  if(nSplit < 7) then
    write(lerr,"(a,i0)") "ELENS> ERROR Expected at least 7 input parameters, got ",nSplit
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
    write(lerr,"(a)") "ELENS> ERROR Element '"//trim(lnSplit(1))//"' not found in single element list."
    iErr = .true.
    return
  end if

  if(kz(iElem) /= elens_kz) then
    write(lerr,"(3(a,i0))") "ELENS> ERROR Element type is kz(",iElem,") = ",kz(iElem)," != ",elens_kz
    iErr = .true.
    return
  end if
  if(el(iElem) /= zero .or. ek(iElem) /= zero .or. ed(iElem) /= zero) then
    write(lerr,"(a)")       "ELENS> ERROR Length el(iElem) (elens is treated as thin element), "//&
      "and first and second field have to be zero:"
    write(lerr,"(2(a,i0),a)") "ELENS>       el(",iElem,") = ",el(iElem)," != 0"
    write(lerr,"(2(a,i0),a)") "ELENS>       ed(",iElem,") = ",ed(iElem)," != 0"
    write(lerr,"(2(a,i0),a)") "ELENS>       ek(",iElem,") = ",ek(iElem)," != 0"
    iErr = .true.
    return
  end if

  melens = melens+1
  if(melens > nelens) then
    write(lerr,"(2(a,i0))") "ELENS> ERROR Too many elenses: ",melens,". Max is ",nelens
    iErr = .true.
    return
  end if

  ielens(iElem) = melens
  if(elens_type(ielens(iElem)) /= 0) then
    write(lerr,"(a)") "ELENS> ERROR The element '"//trim(bez(iElem))//"' was defined twice."
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
      write(lerr,"(a,i0)") "ELENS> ERROR Expected at least 8 input parameters for GAUSSIAN, got ",nSplit
      iErr = .true.
      return
    end if
  case("RADIAL")
    elens_type(ielens(iElem)) = 3
    if(nSplit < 8) then
      write(lerr,"(a,i0)") "ELENS> ERROR Expected at least 8 input parameters for RADIAL, got ",nSplit
      iErr = .true.
      return
    end if
  case("CHEBYSHEV")
    write(lerr,"(a)") "ELENS> ERROR CHEBYSHEV type not fully supported yet - elens name: '"//trim(bez(iElem))
    iErr = .true.
    return
!     elens_type(ielens(iElem)) = 4
!     if(nSplit < 8) then
!       write(lerr,"(a,i0)") "ELENS> ERROR Expected at least 8 input parameters for CHEBYSHEV, got ",nSplit
!       iErr = .true.
!       return
!     end if
  case default
    write(lerr,"(a)") "ELENS> ERROR Elens type '"//trim(lnSplit(2))//"' not recognized. Remember to use all UPPER CASE."
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

  elseif(elens_type(ielens(iElem)) == 3 )then
    ! Radial profile of electrons given by file: need also
    !   name of file where profile is stored
    tmpch = trim(lnSplit(8))

    ! Check if profile has already been requested:
    chIdx = -1
    do tmpi1=1,melens_radial_profiles
      if(tmpch == elens_radial_filename(tmpi1)) then
        elens_iRadial(ielens(iElem)) = tmpi1
        chIdx = tmpi1
        exit
      end if
    end do
    if(chIdx == -1) then
      ! Unsuccessful search
      melens_radial_profiles = melens_radial_profiles+1
      if(melens_radial_profiles > nelens_radial_profiles) then
        write(lerr,"(2(a,i0))") "ELENS> ERROR Too many radial profiles: ",melens_radial_profiles,&
          ". Max is ",nelens_radial_profiles
        iErr = .true.
        return
      end if
      elens_iRadial(ielens(iElem)) = melens_radial_profiles
      elens_radial_filename(tmpi1) = tmpch
    end if
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
  else if(elens_type(ielens(iElem)) == 3 .and. nSplit >= 11) then
    tmpi1 = 9
    tmpi2 = 10
    tmpi3 = 11
    elens_lThetaR2(ielens(iElem)) = .true.
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
    write(lerr,"(a)") "ELENS> ERROR ELEN R2=R1. Elens does not exist."
    iErr = .true.
    return
  end if
  if(elens_r2(ielens(iElem)) <= zero) then
    write(lerr,"(a)") "ELENS> ERROR R2<=0!"
    iErr = .true.
    return
  end if
  if(elens_r1(ielens(iElem)) <= zero) then
    write(lerr,"(a)") "ELENS> ERROR R1<=0!"
    iErr = .true.
    return
  end if
  if(elens_lThetaR2(ielens(iElem))) then
    if(elens_len(ielens(iElem)) <= zero) then
      write(lerr,"(a)") "ELENS> ERROR L<0!"
      iErr = .true.
      return
    end if
    if(elens_I(ielens(iElem)) == zero) then
      write(lerr,"(a)") "ELENS> ERROR I=0!"
      iErr = .true.
      return
    end if
    if(elens_Ek(ielens(iElem)) <= zero) then
      write(lerr,"(a)") "ELENS> ERROR Ek<0! (e-beam)"
      iErr = .true.
      return
    end if
  end if
  if( elens_type(ielens(iElem)) == 2 ) then
    if ( elens_sig(ielens(iElem)).le.zero ) then
       write(lerr,"(a)") "ELENS> ERROR sigma of electron beam <=0 in Elens '"//trim(bez(iElem))//"'."
       iErr = .true.
       return
    end if
  end if

  if(st_debug) then
    call sixin_echoVal("name",    bez(iElem),                   "ELENS",iLine)
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

  use mod_common, only : bez, kz

  implicit none

  logical, intent(inout) :: iErr

  integer jj, kk

  ! check Loop over single elements to check that they have been defined in the fort.3 block
  do jj=1,melens
    if(elens_type(jj)==0) then
      ! find name of elens (for printout purposes)
      do kk=1,nele
        if(kz(kk)==elens_kz) then
          if (ielens(kk).eq.jj) then
            exit
          end if
        end if
      end do
      ! report error
      write(lerr,"(a)") "ELENS> ERROR Type of elens not declared in fort.3 for element '"//trim(bez(kk))//"'"
      iErr = .true.
      return
    end if
  end do

end subroutine elens_parseInputDone

subroutine elens_postInput

  use mathlib_bouncer
  use utils, only : polinterp
  use mod_common, only : bez,kz

  integer j, jj, nlens, jguess
  logical exist

  ! Check that all elenses in fort.2 have a corresponding declaration in fort.3
  nlens=0
  do jj=1,nele
    if(kz(jj)==elens_kz) then
      if (ielens(jj).eq.0) then
        write(lerr,"(a,i0,a)") "ELENS> ERROR single element ",jj," named '"//trim(bez(jj))//"'"
        write(lerr,"(a)")      "ELENS>       does not have a corresponding line in ELEN block in fort.3"
        call prror
      elseif ( elens_type(ielens(jj))==0 ) then
        write(lerr,"(a,i0,a)") "ELENS> ERROR single element ",jj," named '"//trim(bez(jj))//"'"
        write(lerr,"(a)")      "ELENS>       had not been assigned a type"
        call prror
      else
        nlens=nlens+1
      end if
    end if
  end do
  if ( nlens.ne.melens ) then
    write(lerr,"(a,i0)") "ELENS> ERROR number of elenses declared in ELEN block in fort.3 ",melens
    write(lerr,"(a,i0)") "ELENS>       is not the same as the total number of elenses in lattice ",nlens
    call prror
  end if

  ! Parse files with radial profiles
   do j=1,melens_radial_profiles
    inquire(file=elens_radial_filename(j), exist=exist)
    if(.not. exist) then
      write(lerr,"(a)") "ELENS> ERROR Problems with file with radial profile: "//trim(elens_radial_filename(j))
      call prror
    end if
    call parseRadialProfile(j)
    call integrateRadialProfile(j)
    call normaliseRadialProfile(j)
  end do

  ! Proper normalisation
  do j=1,melens
    if(elens_type(j) == 1) then
      ! Uniform distribution
      elens_geo_norm(j) = (elens_r2(j)+elens_r1(j))*(elens_r2(j)-elens_r1(j))
    else if(elens_type(j) == 2) then
      ! Gaussian distribution
      elens_geo_norm(j) = exp_mb(-0.5*(elens_r1(j)/elens_sig(j))**2)&
                         -exp_mb(-0.5*(elens_r2(j)/elens_sig(j))**2)
    else if(elens_type(j) == 3) then
      ! Radial profile
      elens_radial_fr1(j) = polinterp( elens_r1(j), &
            elens_radial_profile_R(0:elens_radial_profile_nPoints(elens_iRadial(j)),elens_iRadial(j)), &
            elens_radial_profile_J(0:elens_radial_profile_nPoints(elens_iRadial(j)),elens_iRadial(j)), &
            elens_radial_profile_nPoints(elens_iRadial(j))+1, &
            elens_radial_mpoints(elens_iRadial(j)), jguess )
      elens_radial_fr2(j) = polinterp( elens_r2(j), &
            elens_radial_profile_R(0:elens_radial_profile_nPoints(elens_iRadial(j)),elens_iRadial(j)), &
            elens_radial_profile_J(0:elens_radial_profile_nPoints(elens_iRadial(j)),elens_iRadial(j)), &
            elens_radial_profile_nPoints(elens_iRadial(j))+1, &
            elens_radial_mpoints(elens_iRadial(j)), jguess  )
      elens_geo_norm(j) = elens_radial_fr2(j) -elens_radial_fr1(j)
    end if

    ! printout:
    ! - find name of elens
    do jj=1,nele
      if(kz(jj)==elens_kz) then
        if (ielens(jj).eq.j) then
          exit
        end if
      end if
    end do
    ! - report geometrical factor
    write(lout,"(a,i0,a,e22.15)") "ELENS> Geom. norm. fact. for elens #",j, &
         " named "//trim(bez(jj))//": ",elens_geo_norm(j)
    
    ! Compute elens theta at R2, if requested by user
    call eLensTheta(j)
    
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
subroutine eLensTheta(j)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants, only : zero, one, two, pi, c1e3, c1m3, c1m6
  use physical_constants, only: clight, pmae, eps0
  use mod_hions, only : zz0
  use mod_common, only : e0, betrel, brho, bez, kz
  use mod_settings, only : st_quiet

  implicit none

  integer j,jj
  real(kind=fPrec) gamma

  if(elens_lThetaR2(j) .and. elens_lAllowUpdate(j)) then
    ! the update of elens_radial_beta_e is not strictly needed here,
    !   apart from the case of elens_Ek is DYNK-ed
    gamma  = ((elens_Ek(j)*c1m3)/pmae)+one ! from kinetic energy
    elens_beta_e(j) = sqrt((one+one/gamma)*(one-one/gamma))
    
    ! r2: from mm to m (c1m3)
    ! theta: from rad to mrad (c1e3)
    elens_theta_r2(j) = ((elens_len(j)*abs(elens_I(j)))/ &
         ((((two*pi)*((eps0*clight)*clight))*brho)*(elens_r2(j)*c1m3)))*c1e3
    if(elens_I(j) < zero) then
      elens_theta_r2(j) = elens_theta_r2(j)*(one/(elens_beta_e(j)*betrel)+one)
    else
      elens_theta_r2(j) = elens_theta_r2(j)*(one/(elens_beta_e(j)*betrel)-one)
    end if
    
    if ( elens_type(j)>=2 ) elens_theta_r2(j) = elens_theta_r2(j) * elens_geo_norm(j)
    
    if(st_quiet < 2) then
      ! find name of elens
      do jj=1,nele
        if(kz(jj)==elens_kz) then
          if (ielens(jj).eq.j) then
            exit
          end if
        end if
      end do
      write(lout,"(a,i0,a,e22.15)") "ELENS> New theta at r2 for elens #",j, &
           " named "//trim(bez(jj))//": ",elens_theta_r2(j)
      if ( elens_type(j)>=2 ) write(lout,"(a,e22.15)") "ELENS>   ...considering also geom. norm. fact.: ", &
           elens_geo_norm(j)
    end if
  end if
 
end subroutine eLensTheta

! ================================================================================================ !
!  Last modified: 2018-09-10
!  Read file with radial profile of electron beam
!  ifile is index of file in table of radial profiles
!  file is structured as:
!     r[mm] j[A/cm2]
!  comment line is headed by '#'
! ================================================================================================ !
subroutine parseRadialProfile(ifile)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants, only: zero, c1e2
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
  logical spErr, err

  integer iErr, ii, fUnit
  real(kind=fPrec) tmpR, tmpJ

  ierr = 0
  ii = 0
  elens_radial_profile_R(ii,ifile) = zero
  elens_radial_profile_J(ii,ifile) = zero
  write(lout,"(a)") "ELENS> Parsing file with radial profile "//trim(elens_radial_filename(ifile))
  call f_requestUnit(elens_radial_filename(ifile),fUnit)
  call f_open(unit=fUnit,file=elens_radial_filename(ifile),mode='r',err=err,formatted=.true.,status="old")

10 continue
  read(fUnit,"(a)",end=20,err=30) inLine
  if(inLine(1:1) == "#") goto 10

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "ELENS> ERROR Failed to parse input line from radial profile."
    goto 30
  end if

  ! Read radial profile
  if(nSplit<2) then
    iErr = 1
    goto 30
  end if
  call chr_cast(lnSplit(1),tmpR,spErr)
  call chr_cast(lnSplit(2),tmpJ,spErr)
  if(tmpR<=elens_radial_profile_R(ii,ifile)) then
    iErr = 1
    write(lerr,"(a,i0)") "ELENS> ERROR radius not in increasing order at ii=",ii
    goto 30
  end if
  if(tmpJ>=0.0) then
    ii=ii+1
    if(ii>elens_radial_dim) then
      iErr = 2
      write(lerr,"(a,i0,a,i0)") "ELENS> ERROR too many points in radial profile: ",ii, &
           ". Max is ",elens_radial_dim
      goto 30
    end if
    elens_radial_profile_nPoints(ifile) = ii
    elens_radial_profile_R(ii,ifile) = tmpR
    elens_radial_profile_J(ii,ifile) = tmpJ/c1e2
  end if

  goto 10

20 continue

  call f_close(fUnit)
  write(lout,"(a,i0,a)") "ELENS> ...acquired ",elens_radial_profile_nPoints(ifile)," points."

  if(st_quiet < 2) then
    ! Echo parsed data (unless told to be quiet!)
    write(lout,"(a,i0)") "ELENS> Radial profile as from file "//&
      trim(elens_radial_filename(ifile))//" - #",ifile
    do ii=0,elens_radial_profile_nPoints(ifile)
      if(elens_radial_profile_J(ii,ifile)/= zero) then
        write(lout,"((a,i4),2(a,1pe22.15))") "ELENS> ",ii,",",elens_radial_profile_R(ii,ifile),",",elens_radial_profile_J(ii,ifile)
      end if
    end do
  end if
  return

30 continue
  write(lerr,"(a,i0,a)") "ELENS> ERROR ",iErr," while parsing file "//trim(elens_radial_filename(ifile))
  call prror(-1)

end subroutine parseRadialProfile

! ================================================================================================ !
!  Last modified: 2018-10-04
!  integrate radial profile of electron beam
!  ifile is index of file in table of radial profiles
!  original formula:
!     cdf(ii)=2pi*pdf(ii)*Dr*r_ave
!  becomes:
!     cdf(ii)=pi*pdf(ii)*(r(ii)-r(ii-1))*(r(ii)+r(ii-1))
! ================================================================================================ !
subroutine integrateRadialProfile(ifile)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use physical_constants
  use crcoall
  use utils, only: polintegrate
  use mod_alloc, only: alloc, dealloc
  use mod_settings, only: st_quiet

  implicit none

  integer, intent(in) :: ifile

  integer ii, nn
  real(kind=fPrec) tmpTot
  real(kind=fPrec), allocatable :: cumul(:)

  write(lout,"(a)") "ELENS> Normalising radial profile described in "//trim(elens_radial_filename(ifile))
  nn=elens_radial_profile_nPoints(ifile)
  call alloc(cumul,nn,zero,'cumul')
  tmpTot=polintegrate(elens_radial_profile_R(0:nn,ifile), elens_radial_profile_J(0:nn,ifile), &
                      nn+1, elens_radial_mpoints(ifile), 2, cumul)
  elens_radial_profile_J(0:nn,ifile)=cumul(1:nn+1)
  call dealloc(cumul,'cumul')
  
  if(st_quiet < 2) then
    write(lout,"(a,i0)") "ELENS> Integrated radial profile read from file "//&
      trim(elens_radial_filename(ifile))//" - #",ifile
    do ii=0,elens_radial_profile_nPoints(ifile)
      if(elens_radial_profile_J(ii,ifile)/= zero) then
        write(lout,"((a,i4),2(a,1pe22.15))") "ELENS> ",ii,",",elens_radial_profile_R(ii,ifile),",",elens_radial_profile_J(ii,ifile)
      end if
    end do
  end if
  write(lout,"(a,1pe22.15)") "ELENS> Total current in radial profile [A]: ", &
         elens_radial_profile_J(elens_radial_profile_nPoints(ifile),ifile)

end subroutine integrateRadialProfile

! ================================================================================================ !
!  Last modified: 2018-10-04
!  normalise integrated radial profiles of electron beam
!  ifile is index of file in table of radial profiles
! ================================================================================================ !
subroutine normaliseRadialProfile(ifile)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use physical_constants
  use crcoall

  implicit none

  integer, intent(in) :: ifile

  integer ii

  do ii=0,elens_radial_profile_nPoints(ifile)
    elens_radial_profile_J(ii,ifile)=elens_radial_profile_J(ii,ifile)/&
                                     elens_radial_profile_J(elens_radial_profile_nPoints(ifile),ifile)
  end do

end subroutine normaliseRadialProfile

! ================================================================================================ !
!  Last modified: 2019-04-11
!  compute kick from electron lens
! ================================================================================================ !
subroutine elens_kick(i,ix,n)

  use mod_common, only : betrel, napx
  use mod_hions, only : moidpsv
  use mod_common_main
  use mathlib_bouncer
  use numerical_constants, only : zero, one
  use utils, only : polinterp

  implicit none
  
  integer, intent(in) :: i
  integer, intent(in) :: ix
  integer, intent(in) :: n
  
  real(kind=fPrec) xx, yy, rr, frr
  integer          jj
  
  do jj=1,napx
    ! 1) apply offset of e-lens
    !    xx = x(proton) - elens_offset_x
    !    yy = y(proton) - elens_offset_y
    xx=xv1(jj)-elens_offset_x(ielens(ix))
    yy=xv2(jj)-elens_offset_y(ielens(ix))
    ! 2) calculate radius
    !    radial position of main beam relative to center of elens beam
    !    rr = sqrt(xx**2+yy**2)
    rr=sqrt(xx*xx+yy*yy)
    ! 3) calculate kick
    !    shape function: spatial charge density depends on type:
    !    0    if r < R1
    !    frr  if R1 < r < R2
    !    1    if r > R2
    if (rr.gt.elens_r1(ielens(ix))) then ! rr <= r1 -> no kick from elens
      if (rr.lt.elens_r2(ielens(ix))) then ! r1 < rr < r2
        select case (elens_type(ielens(ix)))
        case (1)
          ! UNIFORM: eLens with uniform profile
          ! formula: (r^2-r1^2)/(r2^2-r1^2)
          frr=( (rr+elens_r1(ielens(ix)))*(rr-elens_r1(ielens(ix))) )/elens_geo_norm(ielens(ix))
        case (2)
          ! GAUSSIAN: eLens with Gaussian profile
          ! formula: (exp(-r1^2/2sig^2)-exp(-r^2/2sig^2))/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))
          frr=( exp_mb(-0.5*(elens_r1(ielens(ix))/elens_sig(ielens(ix)))**2)    &
               -exp_mb(-0.5*(rr                  /elens_sig(ielens(ix)))**2) )/ &
                                  elens_geo_norm(ielens(ix))
        case (3)
          ! RADIAL PROFILE: eLens with radial profile as from file
          ! formula: (cumul_J(r)-cumul_J(r1))/(cumul_J(r2)-cumul_J(r1))
          frr=(polinterp( rr, &
                elens_radial_profile_R(0:elens_radial_profile_nPoints(elens_iRadial(ielens(ix))),elens_iRadial(ielens(ix))), &
                elens_radial_profile_J(0:elens_radial_profile_nPoints(elens_iRadial(ielens(ix))),elens_iRadial(ielens(ix))), &
                elens_radial_profile_nPoints(elens_iRadial(ielens(ix)))+1, &
                elens_radial_mpoints(ielens(ix)), elens_radial_jguess(ielens(ix)) )-elens_radial_fr1(ielens(ix))) &
                /elens_geo_norm(ielens(ix))
        case default
          write(lerr,"(a,i0,a)") "ELENS> ERROR elens_kick: elens_type=",elens_type(ielens(ix))," not recognized. "
          write(lerr,"(a)")      "ELENS>       Possible values for type are: 1, 2 and 3"
          call prror
        end select
      else ! r1 < r2 <= rr
        frr = one
      endif
      ! 'radial kick'
      frr = (((elens_theta_r2(ielens(ix))*elens_r2(ielens(ix)))/rr)*frr)*moidpsv(jj)
      if(elens_lThetaR2(ielens(ix))) then
        if(elens_I(ielens(ix)) < zero) then
          frr = frr*((rvv(jj)+elens_beta_e(ielens(ix))*betrel)/(one+elens_beta_e(ielens(ix))*betrel))
        else
          frr = frr*((rvv(jj)-elens_beta_e(ielens(ix))*betrel)/(one-elens_beta_e(ielens(ix))*betrel))
        end if
      endif
      yv1(jj)=yv1(jj)-(frr*xx)/rr
      yv2(jj)=yv2(jj)-(frr*yy)/rr
    endif
  end do
end subroutine elens_kick

! ================================================================================================ !
!  Last modified: 2019-04-11
!  compute kick from electron lens
! ================================================================================================ !
subroutine elens_kick_fox(i,ix)

  use mod_common, only : betrel, napx, mtcda
  use mod_hions, only : moidpsv
  use mod_common_main
  use mathlib_bouncer
  use numerical_constants, only : zero, one
  use utils, only : huntBin, polcof
  use mod_lie_dab, only : lnv, idao, rscrri, iscrda
  use mod_common_track, only : xxtr, yytr, issss, comt_daStart, comt_daEnd
  use mod_common_da

  implicit none
  
  integer, intent(in) :: i
  integer, intent(in) :: ix
  
  integer          :: jj, idaa, iLens, iRadial, nBin, nPoints, mPoints, kMin, kMax, kk
  integer          :: hh(lnv)=0
  real(kind=fPrec) :: rra, xclo, yclo, ele_r1, ele_r2, elenor, elesig, elefr1, elebet, rkk, tmpcof
  real(kind=fPrec), allocatable :: cof(:)

  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  save
!-----------------------------------------------------------------------
!FOX  B D ;
#include "include/dainicom.f90"
!FOX  D V DA INT XI NORD NVAR ; D V DA INT YI  NORD NVAR ;
!FOX  D V DA INT RR NORD NVAR ; D V DA INT FRR NORD NVAR ;
!FOX  D V RE INT XCLO ; D V RE INT YCLO ;
!FOX  D V RE INT BETREL ;
!FOX  D V RE INT ELE_R1 ; D V RE INT ELE_R2 ;
!FOX  D V RE INT ELEFR1 ;
!FOX  D V RE INT ELENOR ;
!FOX  D V RE INT ELESIG ;
!FOX  D V RE INT ELEBET ;
!FOX  D V DA INT LINLIN NORD NVAR ;
!FOX  D V DA INT TMPTMP NORD NVAR ;
!FOX  D V RE INT TMPCOF ;
!FOX  D V RE INT RKK ;
!FOX  D V RE INT ONE ; D V RE INT TWO ; D V RE INT HALF ;
!FOX  D V RE INT RRA ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
  call comt_daStart
!-----------------------------------------------------------------------

  iLens=ielens(ix)
  XCLO=elens_offset_x(iLens)
  YCLO=elens_offset_y(iLens)
  ELE_R1=elens_r1(iLens)
  ELE_R2=elens_r2(iLens)
  ELENOR=elens_geo_norm(iLens)
  ELESIG=elens_sig(iLens)
  ELEFR1=elens_radial_fr1(iLens)
  ELEBET=elens_beta_e(iLens)
  
  ! 1) apply offset of e-lens
  !    xx = x(proton) - elens_offset_x
  !    yy = y(proton) - elens_offset_y
!FOX    XI=(XX(1)-XCLO) ;
!FOX    YI=(XX(2)-YCLO) ;
  
  ! 2) calculate radius
  !    radial position of main beam relative to center of elens beam
  !    rr = sqrt(xx**2+yy**2)
!FOX    RR=SQRT(XI*XI+YI*YI) ;
  
  ! 3) calculate kick
  !    shape function: spatial charge density depends on type:
  !    0    if r < R1
  !    frr  if R1 < r < R2
  !    1    if r > R2
  call dapek(RR,hh,RRA)
  if (RRA.gt.elens_r1(iLens)) then ! rr <= r1 -> no kick from elens
    if (RRA.lt.elens_r2(iLens)) then ! r1 < rr < r2
      select case (elens_type(iLens))
      case (1)
        ! UNIFORM: eLens with uniform profile
        ! formula: (r^2-r1^2)/(r2^2-r1^2)
!FOX    FRR=((RR+ELE_R1)*(RR-ELE_R1))/ELENOR ;
         
      case (2)
        ! GAUSSIAN: eLens with Gaussian profile
        ! formula: (exp(-r1^2/2sig^2)-exp(-r^2/2sig^2))/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))
!FOX    FRR=(EXP(-HALF*((ELE_R1/ELESIG)*(ELE_R1/ELESIG)))-EXP(-HALF*((RR/ELESIG)*(RR/ELESIG))))/ELENOR ;

      case (3)
        ! RADIAL PROFILE: eLens with radial profile as from file
        ! formula: (cumul_J(r)-cumul_J(r1))/(cumul_J(r2)-cumul_J(r1))
        iRadial=elens_iRadial(iLens)
        nPoints=elens_radial_profile_nPoints(iRadial)
        mPoints=elens_radial_mpoints(iLens)
        nBin=huntBin(RRA,elens_radial_profile_R(0:nPoints,iRadial),nPoints+1,-1)
        kMin=min(max(nBin-(mPoints-1)/2,1),nPoints+2-mPoints)
        kMax=min(kMin+mPoints-1,nPoints+1)
        call alloc(cof,mPoints,zero,'cof')
        call polcof(elens_radial_profile_R(kMin:kMax,iRadial),elens_radial_profile_J(kMin:kMax,iRadial),mPoints,cof)
        TMPCOF=COF(1)
!FOX    LINLIN=TMPCOF ;
!FOX    TMPTMP=RR;
        do kk=2,mPoints
           TMPCOF=COF(kk)
           RKK=real(kk-1,fPrec)
!FOX       LINLIN=LINLIN+(TMPTMP*TMPCOF)/RKK ;
!FOX       TMPTMP=TMPTMP*RR ;
        end do
!FOX    FRR=(LINLIN-ELEFR1)/ELENOR ;
        call dealloc(cof,'cof')

      case default
        write(lerr,"(a,i0,a)") "ELENS> ERROR in elens_kick_fox: elens_type=",elens_type(ielens(ix))," not recognized. "
        write(lerr,"(a)")      "ELENS>       Possible values for type are: 1, 2 and 3"
        call prror
      end select
    else ! r1 < r2 <= rr
!FOX    FRR=ONE ;
    endif
    ! 'radial kick'
!FOX    FRR=(((ELE_R2*ELE_R2)/RR)*FRR)*(MTCDA/(ONE+DPDA)) ;
    if(elens_lThetaR2(ielens(ix))) then
      if(elens_I(ielens(ix)) < zero) then
!FOX    FRR=FRR*((RV+ELEBET*BETREL)/(ONE+ELEBET*BETREL)) ;
      else
!FOX    FRR=FRR*((RV-ELEBET*BETREL)/(ONE-ELEBET*BETREL)) ;
      end if
    endif
!FOX  YY(1)=YY(1)-(FRR*XI)/RR ;
!FOX  YY(2)=YY(2)-(FRR*YI)/RR ;
  endif

  call comt_daEnd
end subroutine elens_kick_fox

#ifdef CR
subroutine elens_crcheck(fileUnit,readErr)
  implicit none
  integer, intent(in)  :: fileUnit
  logical, intent(out) :: readErr

  integer j

  read(fileUnit,err=10,end=10) (elens_lAllowUpdate_CR(j), j=1, nelens)

  readErr = .false.
  return

10 continue

  write(lerr,"(a,i0)") "SIXTRACR> ERROR in elens_crcheck; fileUnit = ",fileUnit
  write(93,  "(a,i0)") "SIXTRACR> ERROR in elens_crcheck; fileUnit = ",fileUnit
  readErr = .true.

end subroutine elens_crcheck

subroutine elens_crpoint(fileUnit, writeErr,iErro)
  implicit none

  integer, intent(in)    :: fileUnit
  logical, intent(inout) :: writeErr
  integer, intent(inout) :: iErro

  integer j

  write(fileunit,err=10,iostat=iErro) (elens_lAllowUpdate(j), j=1, nelens)
  endfile(fileunit,iostat=iErro)
  backspace(fileunit,iostat=iErro)

  writeErr = .false.
  return

10 continue

  writeErr = .true.
  return

end subroutine elens_crpoint

subroutine elens_crstart
  implicit none
  elens_lAllowUpdate(1:nelens) = elens_lAllowUpdate_CR(1:nelens)
end subroutine elens_crstart

#endif

end module elens
