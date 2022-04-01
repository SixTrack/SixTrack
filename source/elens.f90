! M. Fitterer, FNAL, A. Mereghtti, CERN
! last modified: 02-03-2020
! Common block for electron lens definition
module elens

  use parpro
  use floatPrecision
  use numerical_constants, only : zero, one
  use physical_constants, only: pmae
  implicit none
  private
  public :: &
       elens_allocate_arrays, elens_expand_arrays, &
#ifdef CR
       elens_crcheck, elens_crpoint, elens_crstart, &
#endif
       elens_kick, elens_kick_fox, eLensTheta, &
       elens_parseInputDone, elens_parseInputLine, &
       elens_postInput, elens_postLinopt, &
       elens_setOptics

  integer, allocatable, public, save  :: ielens(:)          ! index of elens (nele)
  integer, public, save               :: nelens=0           ! number of elenses lenses actually in memory
  integer, save               :: nelens_radial_profiles=0   ! number of radial profiles available in memory
  integer, parameter, public  :: elens_kz=29                ! kz of electron lenses
  integer, parameter, public  :: elens_ktrack=63            ! ktrack of electron lenses
  logical, save               :: lRequestOptics=.false.     ! request optics calculations, since there is at least
                                                            !    a lens where either R1 or R2 was given in normalised units
  ! R1/R2/Rref expressed in normalised units (either betatron or dispersion)
  real(kind=fPrec), save      :: elens_emin_def=-one        ! default normalised emittance [m rad]
  real(kind=fPrec), save      :: elens_sigdpp_def=-one      ! default rms of distribution in delta_p / p []
  integer, save               :: elens_iSet_def=0           ! default iSet
  logical, save               :: elens_lFox_def=.true.      ! default lFox
  integer, save               :: elens_radial_mpoints_def=2 ! default number of points for interpolating radial profiles from ASCI files
  integer, save               :: elens_radial_mpoints_ori=3 ! default minimum number of points for interpolating radial profiles from ASCI files in case we are very close to the origin

  ! beam of the lens
  real(kind=fPrec), save      :: elens_beam_mass_def=pmae   ! default mass of lens beam [MeV/c2]
  real(kind=fPrec), save      :: elens_beam_chrg_def=-one   ! default charge of lens beam [e]
  
  ! variables to save elens parameters for tracking etc.
  integer, allocatable, save          :: elens_type(:)              ! integer for elens type (nelens)
                                                                    ! 0 : Un-initialized.
                                                                    ! 1 : uniform profile
                                                                    ! 2 : Gaussian profile
                                                                    ! 3 : radial profile from file
                                                                    ! 4 : degenerate case: wire
  real(kind=fPrec), allocatable, public, save :: elens_theta_ref(:) ! kick strength at Rref [mrad] (nelens)
  real(kind=fPrec), allocatable, save :: elens_r2(:)                ! outer radius R2 [mm] (nelens)
  real(kind=fPrec), allocatable, save :: elens_r1(:)                ! inner radius R1 [mm] (nelens)
  real(kind=fPrec), allocatable, save :: elens_rref(:)              ! radius at which theta_ref is given [mm] (nelens)
  real(kind=fPrec), allocatable, save :: elens_offset_x(:)          ! hor offset of elens [mm] (nelens)
  real(kind=fPrec), allocatable, save :: elens_offset_y(:)          ! vert. offset of elens [mm] (nelens)
  real(kind=fPrec), allocatable, save :: elens_sig(:)               ! sig (Gaussian profile) [mm] (nelens)
  real(kind=fPrec), allocatable, save :: elens_geo_norm(:)          ! normalisation of f(r) (nelens)
  real(kind=fPrec), allocatable, save :: elens_len(:)               ! length of eLens (e-beam region) [m] (nelens)
  real(kind=fPrec), allocatable, public, save :: elens_I(:)         ! current of e-beam [A] (nelens)
                                                                    ! <0: e-beam opposite to beam
  real(kind=fPrec), allocatable, public, save :: elens_Ek(:)        ! kinetic energy of e-beam [keV] (nelens)
  real(kind=fPrec), allocatable, save :: elens_beta_lens_beam(:)    ! relativistic beta of lens beam (nelens)
  logical, allocatable, save          :: elens_lThetaRref(:)        ! flag for computing theta@Rref (nelens)
  logical, allocatable, public, save  :: elens_lAllowUpdate(:)      ! Flag for disabling updating of kick, (nelens)
                                                                    !  i.e. after DYNK has touched thetaRref,
                                                                    !  the energy update is disabled.
  logical, allocatable, public, save  :: elens_lFox(:)              ! the kick from the elens should be taken into account
                                                                    !  when searching for the closed orbit with DA algebra (nelens)
  logical, allocatable, save          :: elens_lFull(:)             ! if .true., elens is full, i.e. not hollow and R1=0.0 (nelens)
  logical, allocatable, save          :: elens_lZeroThick(:)        ! if .true., elens has no thickness, i.e. R2=R1 (nelens)
#ifdef CR
  logical, allocatable, save          :: elens_lAllowUpdate_CR(:)   ! (nelens)
#endif
  ! R1/R2/Rref expressed in normalised units (either betatron or dispersion)
  real(kind=fPrec), allocatable, save :: elens_emin(:)              ! normalised emittance [m rad] (nelens)
  real(kind=fPrec), allocatable, save :: elens_sigdpp(:)            ! default rms of distribution in delta_p / p [] (nelens)
  integer, allocatable, save          :: elens_iSet(:)              ! which value of beta or disp should be taken (nelens):
                                                                    ! 0: none
                                                                    ! 1: min
                                                                    ! 2: max
                                                                    ! 3: average
                                                                    ! 4: quadratic average
                                                                    ! 5: geometric average
  real(kind=fPrec), allocatable, save :: elens_optVal(:)            ! value of optics function for normalised settings (nelens)
  integer, allocatable, save          :: elens_nUpdates(:)          ! how many slices are being used (nelens)
                                                                  
  ! beam of the lens
  real(kind=fPrec), allocatable, save :: elens_beam_mass(:)         ! mass of lens beam [MeV/c2] (nelens)
  real(kind=fPrec), allocatable, save :: elens_beam_chrg(:)         ! charge of lens beam [e] (nelens)
  
  ! radial profile
  integer, allocatable, save          :: elens_iRadial(:)           ! mapping to the radial profile (nelens)
  real(kind=fPrec), allocatable, save :: elens_radial_fr1(:)        ! value of f(R1) in case of radial profiles from file [0:1] (nelens)
  real(kind=fPrec), allocatable, save :: elens_radial_fr2(:)        ! value of f(R2) in case of radial profiles from file [0:1] (nelens)
  integer, allocatable, save          :: elens_radial_mpoints(:)    ! how many points for polynomial interpolation (nelens) (default: 2,
                                                                    !    i.e. linear interpolation
  integer, allocatable, save          :: elens_radial_jguess(:)     ! bin for guessed search (nelens)
  ! - file handling and data storage:
  character(len=:), allocatable, save :: elens_radial_filename(:)        ! names (nelens_radial_profiles)
  real(kind=fPrec), allocatable, save :: elens_radial_profile_R(:,:)     ! [mm] (0:elens_radial_dim,nelens_radial_profiles)
  real(kind=fPrec), allocatable, save :: elens_radial_profile_J(:,:)     ! [A]  (0:elens_radial_dim,nelens_radial_profiles)
  integer, allocatable, save          :: elens_radial_profile_nPoints(:) ! number of points in current radial profile (nelens_radial_profiles)

contains

subroutine elens_allocate_arrays
  use mod_alloc, only : alloc
  implicit none
  integer stat
  call alloc(ielens,nele,0,'ielens')
end subroutine elens_allocate_arrays

subroutine elens_expand_arrays(nele_new)
  use mod_alloc, only : alloc
  implicit none
  integer, intent(in) :: nele_new
  call alloc(ielens,nele_new,0,'ielens')
end subroutine elens_expand_arrays

subroutine elens_expand_arrays_lenses(nelens_new)
  use mod_alloc, only : alloc
  implicit none
  integer, intent(in) :: nelens_new
  ! elens charachteristics
  call alloc(elens_type           , nelens_new,                        0, 'elens_type'           )
  call alloc(elens_theta_ref      , nelens_new,                     zero, 'elens_theta_ref'      )
  call alloc(elens_rref           , nelens_new,                     zero, 'elens_rref'           )
  call alloc(elens_r2             , nelens_new,                     zero, 'elens_r2'             )
  call alloc(elens_r1             , nelens_new,                     zero, 'elens_r1'             )
  call alloc(elens_offset_x       , nelens_new,                     zero, 'elens_offset_x'       )
  call alloc(elens_offset_y       , nelens_new,                     zero, 'elens_offset_y'       )
  call alloc(elens_sig            , nelens_new,                     zero, 'elens_sig'            )
  call alloc(elens_geo_norm       , nelens_new,                     zero, 'elens_geo_norm'       )
  call alloc(elens_len            , nelens_new,                     zero, 'elens_len'            )
  call alloc(elens_I              , nelens_new,                     zero, 'elens_I'              )
  call alloc(elens_Ek             , nelens_new,                     zero, 'elens_Ek'             )
  call alloc(elens_beta_lens_beam , nelens_new,                     zero, 'elens_beta_lens_beam' )
  call alloc(elens_lThetaRref     , nelens_new,                  .false., 'elens_lThetaRref'     )
  call alloc(elens_lAllowUpdate   , nelens_new,                   .true., 'elens_lAllowUpdate'   )
  call alloc(elens_lFox           , nelens_new,           elens_lFox_def, 'elens_lFox'           )
  call alloc(elens_lFull          , nelens_new,                  .false., 'elens_lFull'          )
  call alloc(elens_lZeroThick     , nelens_new,                  .false., 'elens_lZeroThick'     )
#ifdef CR
  call alloc(elens_lAllowUpdate_CR, nelens_new,                  .false., 'elens_lAllowUpdate_CR')
#endif                                                            
  call alloc(elens_emin           , nelens_new,           elens_emin_def, 'elens_emin'           )
  call alloc(elens_sigdpp         , nelens_new,         elens_sigdpp_def, 'elens_sigdpp'         )
  call alloc(elens_iSet           , nelens_new,           elens_iSet_def, 'elens_iSet'           )
  call alloc(elens_optVal         , nelens_new,                     zero, 'elens_optVal'         )
  call alloc(elens_nUpdates       , nelens_new,                        0, 'elens_nUpdates'       )
  call alloc(elens_beam_mass      , nelens_new,      elens_beam_mass_def, 'elens_beam_mass'      )
  call alloc(elens_beam_chrg      , nelens_new,      elens_beam_chrg_def, 'elens_beam_chrg'      )
  call alloc(elens_iRadial        , nelens_new,                        0, 'elens_iRadial'        )
  call alloc(elens_radial_fr1     , nelens_new,                     zero, 'elens_radial_fr1'     )
  call alloc(elens_radial_fr2     , nelens_new,                     zero, 'elens_radial_fr2'     )
  call alloc(elens_radial_mpoints , nelens_new, elens_radial_mpoints_def, 'elens_radial_mpoints' )
  call alloc(elens_radial_jguess  , nelens_new,                       -1, 'elens_radial_jguess'  )
end subroutine elens_expand_arrays_lenses

subroutine elens_expand_arrays_rad_profiles(nelens_profiles_new)
  use mod_alloc, only : alloc
  implicit none
  integer, intent(in) :: nelens_profiles_new
  call alloc(elens_radial_filename,  mFileName,     nelens_profiles_new,  " ", 'elens_radial_filename' )
  call alloc(elens_radial_profile_R,         0,     nelens_profiles_new, zero, 'elens_radial_profile_R', &
                                             0,     nelens_radial_profiles )
  call alloc(elens_radial_profile_J,         0,     nelens_profiles_new, zero, 'elens_radial_profile_J', &
                                             0,     nelens_radial_profiles )
  call alloc(elens_radial_profile_nPoints,          nelens_profiles_new,    0, 'elens_radial_profile_nPoints' )
end subroutine elens_expand_arrays_rad_profiles

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
  use crcoall, only : lout, lerr

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mStrLen) tmpch
  real(kind=fPrec) tmpflt, tmpflt2
  integer nSplit, iElem, j, tmpi1, tmpi2, tmpi3
  logical spErr, tmpl, lfound

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "ELENS> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case (lnSplit(1))
     
  case("EMIN")
    if(nSplit<2 .or. nSplit>4 ) then
      write(lerr,"(a,i0)") "ELENS> ERROR Expected 1, 2 or 3 input parameters for EMIttance (Normalised) line, got ",nSplit-1
      write(lerr,"(a)")    "ELENS>       example:     EMIN  <emin> (min|max|ave|qve|geo) (ALL|BEF(ORE)|AFT(ER))"
      iErr = .true.
      return
    end if
   
    call chr_cast(lnSPlit(2),tmpflt,iErr)
    if ( tmpflt < zero ) then
      write(lerr,"(a,1pe22.15)") "ELENS> ERROR Negative values of normalised emittance are unacceptable, got: ", tmpflt
      iErr = .true.
      return
    end if
    if (nelens>0) then
      elens_emin(nelens)=tmpflt
      elens_sigdpp(nelens)=-one
    end if
    
    if (nSplit>=3) then
      select case (chr_toLower(trim(lnSplit(3))))
      case('min')
        tmpi2=1
        if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance with min beta"
      case('max')
        tmpi2=2
        if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance with max beta"
      case('ave')
        tmpi2=3
        if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance with average beta"
      case('qve')
        tmpi2=4
        if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance with quad-average beta"
      case('gve')
        tmpi2=5
        if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance with geometric average beta"
      case default
        write(lerr,"(a)") "ELENS> ERROR Unidentified third parameter of EMIN line, got: '"//trim(lnSplit(3))//"'"
        write(lerr,"(a)") "ELENS>       example:     EMIN  <emin> (min|max|ave|qve|geo) (ALL|BEF(ORE)|AFT(ER))"
        iErr = .true.
        return
      end select ! case (lnSplit(3))
    else
      tmpi2=1
      if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance with min beta (default)"
    end if
    if (nelens>0) elens_iSet(nelens)=tmpi2
       
    if (nSplit>=4) then
      select case (chr_toLower(trim(lnSplit(4))))
      case('all')
        do tmpi1=1,nelens-1
          elens_emin(tmpi1) = tmpflt
          elens_iSet(tmpi1) = tmpi2
          elens_sigdpp(tmpi1) = -one
        end do
        elens_emin_def=tmpflt
        elens_iSet_def=tmpi2
        elens_sigdpp_def=-one
        if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance to all e-lenses"
      case('bef','before')
        do tmpi1=1,nelens-1
          elens_emin(tmpi1) = tmpflt
          elens_iSet(tmpi1) = tmpi2
          elens_sigdpp(tmpi1) = -one
        end do
        if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance to all e-lenses "// &
             "declared before the current EMIN line"
      case('aft','after')
        elens_emin_def=tmpflt
        elens_iSet_def=tmpi2
        elens_sigdpp_def=-one
        if(st_debug) write(lout,"(a)") "ELENS> Applying read normalised emittance to all e-lenses "// &
             "declared after the current EMIN line"
      case default
        write(lerr,"(a)") "ELENS> ERROR Unidentified third parameter of EMIN line, got: '"//trim(lnSplit(4))//"'"
        write(lerr,"(a)") "ELENS>       example:     EMIN  <emin> (min|max|ave|qve|geo) (ALL|BEF(ORE)|AFT(ER))"
        iErr = .true.
        return
      end select ! case (lnSplit(4))
    end if

    if(st_debug) then
      call sixin_echoVal("normalised emittance [m rad]",tmpflt,"ELENS",iLine)
    end if
    
  case("FOX")
    if(nSplit<2 .or. nSplit>3) then
      write(lerr,"(a,i0)") "ELENS> ERROR Expected 1 or 2 input parameters for FOX line, got ",nSplit-1
      write(lerr,"(a)")    "ELENS>       example:     FOX  on|off|true|false (ALL|BEF(ORE)|AFT(ER))"
      iErr = .true.
      return
    end if
   
    call chr_cast(lnSPlit(2), tmpl,iErr)
    if (nelens>0) elens_lFox(nelens)=tmpl

    if (nSplit>=3) then
      select case (chr_toLower(trim(lnSplit(3))))
      case('all')
        do tmpi1=1,nelens-1
          elens_lFox(tmpi1) = tmpl
        end do
        elens_lFox_def=tmpl
        if(st_debug) write(lout,"(a)") "ELENS> Setting lFox as read to all e-lenses"
      case('bef','before')
        do tmpi1=1,nelens-1
          elens_lFox(tmpi1) = tmpl
        end do
        if(st_debug) write(lout,"(a)") "ELENS> Setting lFox as read to all e-lenses "// &
             "declared before the current FOX line"
      case('aft','after')
        elens_lFox_def=tmpl
        if(st_debug) write(lout,"(a)") "ELENS> Setting lFox as read to all e-lenses "// &
             "declared after the current FOX line"
      case default
        write(lerr,"(a)") "ELENS> ERROR Unidentified third parameter of FOX line, got: '"//trim(lnSplit(3))//"'"
        write(lerr,"(a)") "ELENS>       example:     FOX  on|off|true|false (ALL|BEF(ORE)|AFT(ER))"
        iErr = .true.
        return
      end select ! case (lnSplit(3))
    end if
    
    if(st_debug) then
      call sixin_echoVal("fox",tmpl,"ELENS",iLine)
    end if
     
  case("INTER")
    if(nSplit<2 .or. nSplit>3) then
      write(lerr,"(a,i0)") "ELENS> ERROR Expected 1 or 2 input parameters for INTERpolation line, got ",nSplit-1
      write(lerr,"(a)")    "ELENS>       example:     INTER  <elens_radial_mpoints> (ALL|BEF(ORE)|AFT(ER))"
      iErr = .true.
      return
    end if
    if ( elens_type(nelens).ne.3 ) then
      write(lout,"(a,i0)") "ELENS> WARNING INTERpolation setting for an ELENS type without radial profile - ignoring setting..."
      return
    end if
    call chr_cast(lnSPlit(2), tmpi2,iErr)
    if ( tmpi2 <=0 .or. tmpi2>20 ) then
      write(lerr,"(a,i0)") "ELENS> ERROR Unreasonable number of points for radial interpolation, got ",tmpi2
      write(lerr,"(a)")    "ELENS>       Please choose a value beterrn 1 and 20 (included)"
      iErr = .true.
      return
    end if
    if (nelens>0) elens_radial_mpoints(nelens)=tmpi2
    
    if (nSplit>=3) then
      select case (chr_toLower(trim(lnSplit(3))))
      case('all')
        do tmpi1=1,nelens-1
          elens_radial_mpoints(tmpi1) = tmpi2
        end do
        elens_radial_mpoints_def=tmpi2
        if(st_debug) write(lout,"(a)") "ELENS> Setting elens_radial_mpoints as read to all e-lenses"
      case('bef','before')
        do tmpi1=1,nelens-1
          elens_radial_mpoints(tmpi1) = tmpi2
        end do
        if(st_debug) write(lout,"(a)") "ELENS> Setting elens_radial_mpoints as read to all e-lenses "// &
             "declared before the current FOX line"
      case('aft','after')
        elens_radial_mpoints_def=tmpi2
        if(st_debug) write(lout,"(a)") "ELENS> Setting elens_radial_mpoints as read to all e-lenses "// &
             "declared after the current FOX line"
      case default
        write(lerr,"(a)") "ELENS> ERROR Unidentified third parameter of INTER line, got: '"//trim(lnSplit(3))//"'"
        write(lerr,"(a)") "ELENS>       example:     INTER  <elens_radial_mpoints> (ALL|BEF(ORE)|AFT(ER))"
        iErr = .true.
        return
      end select ! case (lnSplit(3))
    end if
    
    if(st_debug) then
      call sixin_echoVal("number of interpolation points",tmpi2,"ELENS",iLine)
    end if

  case("SIGDPP")
    if(nSplit<2 .or. nSplit>4 ) then
      write(lerr,"(a,i0)") "ELENS> ERROR Expected 1, 2 or 3 input parameters for SIGma Delta_P over P line, got ",nSplit-1
      write(lerr,"(a)")    "ELENS>       example:     SIDPP  <sigdpp> (min|max|ave|qve|geo) (ALL|BEF(ORE)|AFT(ER))"
      iErr = .true.
      return
    end if
   
    call chr_cast(lnSPlit(2), tmpflt,iErr)
    if ( tmpflt < zero ) then
      write(lerr,"(a,1pe22.15)") "ELENS> ERROR Negative values of rms of delta distribution are unacceptable, got: ", &
            tmpflt
      iErr = .true.
      return
    end if
    if (nelens>0) then
      elens_emin(nelens)=-one
      elens_sigdpp(nelens)=tmpflt
    end if

    if (nSplit>=3) then
      select case (chr_toLower(trim(lnSplit(3))))
      case('min')
        tmpi2=1
        if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution with min dispersion"
      case('max')
        tmpi2=2
        if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution with max dispersion"
      case('ave')
        tmpi2=3
        if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution with average dispersion"
      case('qve')
        tmpi2=4
        if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution with quad-average dispersion"
      case('geo')
        tmpi2=5
        if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution with geometric average beta"
      case default
        write(lerr,"(a)") "ELENS> ERROR Unidentified third parameter of SIGDPP line, got: '"//trim(lnSplit(3))//"'"
        write(lerr,"(a)") "ELENS>       example:     SIGDPP  <sigdpp> (min|max|ave|qve|geo) (ALL|BEF(ORE)|AFT(ER))"
        iErr = .true.
        return
      end select ! case (lnSplit(3))
    else
      tmpi2=1
      if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution with min dispersion (default)"
    end if
    if (nelens>0) elens_iSet(nelens)=tmpi2
       
    if (nSplit>=4) then
      select case (chr_toLower(trim(lnSplit(4))))
      case('all')
        do tmpi1=1,nelens-1
          elens_emin(tmpi1) = -one
          elens_iSet(tmpi1) = tmpi2
          elens_sigdpp(tmpi1) = tmpflt
        end do
        elens_emin_def=-one
        elens_sigdpp_def=tmpflt
        elens_iSet_def=tmpi2
        if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution to all e-lenses"
      case('before')
        do tmpi1=1,nelens-1
          elens_emin(tmpi1) = -one
          elens_iSet(tmpi1) = tmpi2
          elens_sigdpp(tmpi1) = tmpflt
        end do
        if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution to all e-lenses"// &
             "declared before the current SIGDPP line"
      case('after')
        elens_emin_def=-one
        elens_iSet_def=tmpi2
        elens_sigdpp_def=tmpflt
        if(st_debug) write(lout,"(a)") "ELENS> Applying read rms of delta distribution to all e-lenses"// &
             "declared after the current SIGDPP line"
      case default
        write(lerr,"(a)") "ELENS> ERROR Unidentified third parameter of SIGDPP line, got: '"//trim(lnSplit(4))//"'"
        write(lerr,"(a)") "ELENS>       example:     SIGDPP  <sigdpp> (min|max|ave|qve|geo) (ALL|BEF(ORE)|AFT(ER))"
        iErr = .true.
        return
      end select ! case (lnSplit(4))
    end if

    if(st_debug) then
      call sixin_echoVal("rms of delta distribution []",tmpflt,"ELENS",iLine)
    end if
    
  case("SPEC")
    if(nSplit<3.or.nSplit>4) then
      write(lerr,"(a,i0)") "ELENS> ERROR Expected 2 or 3 input parameters for SPEC line, got ",nSplit-1
      write(lerr,"(a)")    "ELENS>       example:     SPEC  m[MeV/c2]  Q[e] (ALL|BEF(ORE)|AFT(ER))"
      iErr = .true.
      return
    end if
   
    call chr_cast(lnSPlit(2),tmpflt,iErr)
    if (tmpflt<zero) then
      write(lerr,"(a,1pe22.15)") "ELENS> ERROR negative mass specified in SPEC line, got: ",tmpflt
      iErr = .true.
      return
    end if
    call chr_cast(lnSPlit(3),tmpflt2,iErr)
    if (nelens>0) then
      elens_beam_mass(nelens)=tmpflt
      elens_beam_chrg(nelens)=tmpflt2
    end if
  
    if (nSplit>=4) then
      select case (chr_toLower(trim(lnSplit(4))))
      case('all')
        do tmpi1=1,nelens-1
          elens_beam_mass(tmpi1) = tmpflt
          elens_beam_chrg(tmpi1) = tmpflt2
        end do
        elens_beam_mass_def=tmpflt
        elens_beam_chrg_def=tmpflt2
        if(st_debug) write(lout,"(a)") "ELENS> Applying read mass and charge to all e-lenses"
      case('bef','before')
        do tmpi1=1,nelens-1
          elens_beam_mass(tmpi1) = tmpflt
          elens_beam_chrg(tmpi1) = tmpflt2
        end do
        if(st_debug) write(lout,"(a)") "ELENS> Applying read mass and charge to all e-lenses "// &
             "declared before the current SPEC line"
      case('aft','after')
        elens_beam_mass_def=tmpflt
        elens_beam_chrg_def=tmpflt2
        if(st_debug) write(lout,"(a)") "ELENS> Applying read mass and charge to all e-lenses "// &
             "declared after the current SPEC line"
      case default
        write(lerr,"(a)") "ELENS> ERROR Unidentified fourth parameter of SPEC line, got: '"//trim(lnSplit(4))//"'"
        write(lerr,"(a)") "ELENS>       example:     SPEC  m[MeV/c2]  Q[e] (ALL|BEF(ORE)|AFT(ER))"
        iErr = .true.
        return
      end select ! case (lnSplit(4))
    end if
    
    if(st_debug) then
      call sixin_echoVal("mass of beam in lens [MeV/c2]",tmpflt,"ELENS",iLine)
      call sixin_echoVal("charge of beam in lens [e]",tmpflt2,"ELENS",iLine)
    end if
     
  case default
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
    
    if(ielens(iElem) /= 0) then
      write(lerr,"(a)") "ELENS> ERROR The element '"//trim(bez(iElem))//"' was defined twice."
      iErr = .true.
      return
    end if
    nelens = nelens+1
    call elens_expand_arrays_lenses(nelens)
    ielens(iElem) = nelens
    
    ! Parse the element
    select case (lnSplit(2))
    case("UNIFORM")
      elens_type(ielens(iElem)) = 1
    case("GAUSSIAN")
      elens_type(ielens(iElem)) = 2
      if(nSplit < 8) then
        write(lerr,"(a,i0)") "ELENS> ERROR Expected at least 6 input parameters for GAUSSIAN, got ",nSplit
        iErr = .true.
        return
      end if
    case("RADIAL")
      elens_type(ielens(iElem)) = 3
      if(nSplit < 8) then
        write(lerr,"(a,i0)") "ELENS> ERROR Expected at least 6 input parameters for RADIAL, got ",nSplit
        iErr = .true.
        return
      end if
    case("WIRE")
      elens_type(ielens(iElem)) = 4
    case default
      write(lerr,"(a)") "ELENS> ERROR Elens type '"//trim(lnSplit(2))//"' not recognized. Remember to use all UPPER CASE."
      iErr = .true.
      return
    end select ! case (lnSplit(2))
    
    call chr_cast(lnSplit(3),elens_theta_ref(ielens(iElem)),iErr)
    call chr_cast(lnSplit(4),elens_r2(ielens(iElem)),       iErr)
    call chr_cast(lnSplit(5),elens_r1(ielens(iElem)),       iErr)
    call chr_cast(lnSplit(6),elens_offset_x(ielens(iElem)), iErr)
    call chr_cast(lnSplit(7),elens_offset_y(ielens(iElem)), iErr)
    
    if(elens_type(ielens(iElem)) == 2) then
      ! GAUSSIAN profile of electrons: need also sigma of e-beam
      call chr_cast(lnSplit(8),elens_sig(ielens(iElem)),iErr)
    
    elseif(elens_type(ielens(iElem)) == 3 )then
      ! Radial profile of electrons given by file: need also
      !   name of file where profile is stored
      tmpch = trim(lnSplit(8))
    
      ! Check if profile has already been requested:
      lfound = .false.
      do tmpi1=1,nelens_radial_profiles
        if(tmpch == elens_radial_filename(tmpi1)) then
          elens_iRadial(ielens(iElem)) = tmpi1
          lfound = .true.
          exit
        end if
      end do
      if(.not.lfound) then
        ! Unsuccessful search
        nelens_radial_profiles = nelens_radial_profiles+1
        call elens_expand_arrays_rad_profiles(nelens_radial_profiles)
        elens_iRadial(ielens(iElem)) = nelens_radial_profiles
        elens_radial_filename(nelens_radial_profiles) = tmpch
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
      elens_lThetaRref(ielens(iElem)) = .true.
    else if(elens_type(ielens(iElem)) == 2 .and. nSplit >= 11) then
      tmpi1 = 9
      tmpi2 = 10
      tmpi3 = 11
      elens_lThetaRref(ielens(iElem)) = .true.
    else if(elens_type(ielens(iElem)) == 3 .and. nSplit >= 11) then
      tmpi1 = 9
      tmpi2 = 10
      tmpi3 = 11
      elens_lThetaRref(ielens(iElem)) = .true.
    end if
    
    if(elens_lThetaRref(ielens(iElem))) then
      call chr_cast(lnSplit(tmpi1),elens_len(ielens(iElem)),iErr)
      call chr_cast(lnSplit(tmpi2),elens_I(ielens(iElem)),  iErr)
      call chr_cast(lnSplit(tmpi3),elens_Ek(ielens(iElem)), iErr)
    end if
    
    ! sanity checks
    if(elens_r2(ielens(iElem)) < zero) then
      lRequestOptics=.true.
    else if(abs(elens_r2(ielens(iElem))) < pieni) then
      elens_r2(ielens(iElem))=zero
    end if
    if(elens_r1(ielens(iElem)) < zero) then
      lRequestOptics=.true.
    else if(abs(elens_r1(ielens(iElem))) < pieni) then
      elens_r1(ielens(iElem))=zero
    end if
    if(abs(elens_r2(ielens(iElem))) < abs(elens_r1(ielens(iElem)))) then
      write(lout,"(a)") "ELENS> WARNING ELEN R2<R1. Inverting."
      tmpflt=elens_r2(ielens(iElem))
      elens_r2(ielens(iElem)) = elens_r1(ielens(iElem))
      elens_r1(ielens(iElem)) = tmpflt
    end if
    if(elens_lThetaRref(ielens(iElem))) then
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
      if(elens_sig(ielens(iElem)) < zero) then
        lRequestOptics=.true.
      else if(abs(elens_sig(ielens(iElem))) < pieni) then
        write(lerr,"(a)") "ELENS> ERROR sigma of electron beam ==0 in Elens '"//trim(bez(iElem))//"'."
        iErr = .true.
        return
      end if
    end if
    if( elens_type(ielens(iElem)) == 4 ) then
      ! the user declares the e-lens as a degenerate one to a WIRE
      ! R1 must be 0.0, whereas R2 must be != 0.0 
      elens_lFull(ielens(iElem))=.true.
      elens_lZeroThick(ielens(iElem))=.true.
      if( elens_r1(ielens(iElem)) /= zero ) then
        elens_r1(ielens(iElem))=zero
        write(lout,"(a)") "ELENS> WARNING: Elens '"//trim(bez(iElem))//"' is a wire."
        write(lout,"(a)") "ELENS>          Forcing R1 to 0 mm."
      end if
      if( elens_r2(ielens(iElem)) == zero ) then
        write(lerr,"(a)") "ELENS> ERROR R2==0 for a WIRE e-lens!"
        iErr = .true.
        return
      end if
      elens_rref(ielens(iElem))=elens_r2(ielens(iElem))
      elens_r2(ielens(iElem))=zero
    else
      elens_lFull(ielens(iElem))=elens_r1(ielens(iElem))==zero
      elens_lZeroThick(ielens(iElem))=elens_r1(ielens(iElem))==elens_r2(ielens(iElem))
      if ( elens_lFull(ielens(iElem)) .and. elens_lZeroThick(ielens(iElem)) ) then
        ! degenerate e-lens to a wire, where user decided not to set Rref
        ! set it to 1mm
        elens_rref(ielens(iElem))=one
        write(lout,"(a)") "ELENS> WARNING: Elens '"//trim(bez(iElem))//"' is actually a wire (i.e. R1=R2=0.0)."
        write(lout,"(a)") "ELENS>          Forcing Rref to 1 mm, such that theta_ref keeps a meaning."
        write(lout,"(a)") "ELENS>          If this is not what you want to do, please use the WIRE type or the wire element."
      else
        elens_rref(ielens(iElem))=elens_r2(ielens(iElem))
      end if
    end if
    
    if(st_debug) then
      call sixin_echoVal("name",            bez(iElem),                    "ELENS",iLine)
      call sixin_echoVal("type",            elens_type(ielens(iElem)),     "ELENS",iLine)
      call sixin_echoVal("theta_ref [mrad]",elens_theta_ref(ielens(iElem)),"ELENS",iLine)
      call sixin_echoVal("r1 [mm]",         elens_r1(ielens(iElem)),       "ELENS",iLine)
      call sixin_echoVal("r2 [mm]",         elens_r2(ielens(iElem)),       "ELENS",iLine)
      call sixin_echoVal("rRef [mm]",       elens_rref(ielens(iElem)),     "ELENS",iLine)
      call sixin_echoVal("offset_x [mm]",   elens_offset_x(ielens(iElem)), "ELENS",iLine)
      call sixin_echoVal("offset_y [mm]",   elens_offset_y(ielens(iElem)), "ELENS",iLine)
      if(elens_type(ielens(iElem)) == 2) then
        call sixin_echoVal("sig [mm]",elens_sig(ielens(iElem)),"ELENS",iLine)
      end if
      if(elens_lThetaRref(ielens(iElem))) then
        call sixin_echoVal("L [m]",   elens_len(ielens(iElem)),"ELENS",iLine)
        call sixin_echoVal("I [A]",   elens_I(ielens(iElem)),  "ELENS",iLine)
        call sixin_echoVal("Ek [keV]",elens_Ek(ielens(iElem)), "ELENS",iLine)
      end if
      if ( elens_lFull(ielens(iElem)) ) then
        if ( elens_lZeroThick(ielens(iElem)) ) then
          write(lout,"(a)") "ELENS> Elens '"//trim(bez(iElem))//"' is actually a wire (i.e. R1=R2=0.0)."
        else
          write(lout,"(a)") "ELENS> Elens '"//trim(bez(iElem))//"' is full (i.e. R1=0.0, R2>0.0)."
        end if
      else 
        if ( elens_lZeroThick(ielens(iElem)) ) then
          write(lout,"(a)") "ELENS> Elens '"//trim(bez(iElem))//"' is a zero-thickness hollow electron lens (i.e. R1=R2>0.0)."
        else
          write(lout,"(a)") "ELENS> Elens '"//trim(bez(iElem))//"' is a regular hollow electron lens (i.e. R1>0.0, R2>0.0, R1<R2)."
        end if
      end if
    end if
 
  end select ! case (lnSplit(1))

end subroutine elens_parseInputLine

subroutine elens_parseInputDone(iErr)

  use mod_common, only : bez, kz, fort3
  use crcoall, only : lerr

  implicit none

  logical, intent(inout) :: iErr

  integer jj, kk

  ! check Loop over single elements to check that they have been defined in the fort.3 block
  do jj=1,nelens
    if(elens_type(jj)==0) then
      ! find name of elens (for printout purposes)
      do kk=1,nele
        if(kz(kk)==elens_kz) then
          if (ielens(kk)==jj) then
            exit
          end if
        end if
      end do
      ! report error
      write(lerr,"(a)") "ELENS> ERROR Type of elens not declared in "//trim(fort3)//" for element '"//trim(bez(kk))//"'"
      iErr = .true.
      return
    end if
  end do

end subroutine elens_parseInputDone

subroutine elens_postInput

  use mod_common, only : bez, kz, fort3, ilin
  use crcoall, only : lout, lerr

  integer j, jj, nlens
  logical lexist, lFoxOne

  ! Check that all elenses in fort.2 have a corresponding declaration in fort.3
  nlens=0
  do jj=1,nele
    if(kz(jj)==elens_kz) then
      if (ielens(jj)==0) then
        write(lerr,"(a,i0,a)") "ELENS> ERROR single element ",jj," named '"//trim(bez(jj))//"'"
        write(lerr,"(a)")      "ELENS>       does not have a corresponding line in ELEN block in "//trim(fort3)
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
  if ( nlens.ne.nelens ) then
    write(lerr,"(a,i0)") "ELENS> ERROR number of elenses declared in ELEN block in "//trim(fort3)//" ",nelens
    write(lerr,"(a,i0)") "ELENS>       is not the same as the total number of elenses in lattice ",nlens
    call prror
  end if

  ! forcing calculation of linear optics, if necessary
  if ( lRequestOptics ) then
    write(lout,"(a)") "ELENS> One of the lenses has R1<0 or R2<0 or Rref<0 hence expressed in normalised units"
    write(lout,"(a)") "ELENS>     (i.e. units of beam sigma). Calculation of linear optics functions is necessary"
    select case(ilin)
    case (0) 
      write(lout,"(a)") "ELENS> Forcing calculation of linear optics functions with ilin=1"
      ilin=1
    case (1)
      write(lout,"(a)") "ELENS> Calculation of linear optics functions already requested by user with ilin=1"
    case(2)
      lFoxOne=.false.
      do j=1,nelens
        if (elens_lFox(j)) then
          lFoxOne=.true.
          exit 
        end if
      end do
      if (lFoxOne) then
        write(lerr,"(a)") "ELENS> ERROR cannot compute optics with ilin=2 and at least an e-lens"
        write(lerr,"(a)") "ELENS>       with optics-dependent parameters"
        call prror
      end if
    case default
      write(lerr,"(a,i0)") "ELENS> Unknown value of ilin=",ilin
      call prror
   end select ! case(ilin)
  end if

  ! check that, in case R1/R2/Rref are declared by the user in terms of n-sigma, the provided info is consistent
  do j=1,nelens
    if (elens_r1(j)<zero .or. elens_r2(j)<zero .or. elens_rref(j)<zero .or. elens_sig(j)<zero) then
      ! printout:
      ! - find name of elens
      do jj=1,nele
        if(kz(jj)==elens_kz) then
          if (ielens(jj)==j) then
            exit
          end if
        end if
      end do
      write(lout,"(a,i0,a)") "ELENS> checking consistency of user input data for e-lens #",j," named "//trim(bez(jj))//"..."
      if (elens_r1(j)<zero) write(lout,"(a)") "ELENS> ...e-lens has R1<0"
      if (elens_r2(j)<zero) write(lout,"(a)") "ELENS> ...e-lens has R2<0"
      if (elens_rref(j)<zero) write(lout,"(a)") "ELENS> ...e-lens has Rref<0"
      if (elens_sig(j)<zero) write(lout,"(a)") "ELENS> ...e-lens has sig_el<0"
      if (elens_emin(j)>zero.and.elens_sigdpp(j)>zero) then
        write(lerr,"(a)") "ELENS> ERROR cannot express R1/R2/Rref/sig cuts in terms of both normalised emittance "// &
              "and rms of delta distribution."
        write(lerr,"(a)") "ELENS>       please choose one of the two!"
        write(lerr,"(a,1pe22.15)") "ELENS>       got emin [m rad]=",elens_emin(j)
        write(lerr,"(a,1pe22.15)") "ELENS>       got sigdpp []=",elens_sigdpp(j)
        call prror
      elseif (elens_emin(j)<zero.and.elens_sigdpp(j)<zero) then
        write(lerr,"(a)") "ELENS> ERROR cannot express R1/R2/Rref/sig cuts in terms of both normalised emittance "// &
              "or rms of delta distribution without specifying one of the two."
        write(lerr,"(a)") "ELENS>       please choose one of the two!"
        if (elens_emin(j).ne.-one) then
          ! echo value only if non-default one has been given
          write(lerr,"(a,1pe22.15)") "ELENS>       got emin [m rad]=",elens_emin(j)
        end if
        if (elens_sigdpp(j).ne.-one) then
          ! echo value only if non-default one has been given
          write(lerr,"(a,1pe22.15)") "ELENS>       got sigdpp []=",elens_sigdpp(j)
        end if
        call prror
      end if
      if (elens_iSet(j)==0) then
        write(lerr,"(a)") "ELENS> ERROR don't know which value of beta/dispersion to use"
        write(lerr,"(a)") "ELENS>       please choose one among min|max|ave|qve"
        call prror
      end if
      write(lout,"(a)") "ELENS> ...all checks passed!"
    end if
  end do
  
  ! Parse files with radial profiles
   do j=1,nelens_radial_profiles
    inquire(file=elens_radial_filename(j), exist=lexist)
    if(.not. lexist) then
      write(lerr,"(a)") "ELENS> ERROR Problems with file with radial profile: "//trim(elens_radial_filename(j))
      call prror
    end if
    call parseRadialProfile(j)
    call integrateRadialProfile(j)
    call normaliseRadialProfile(j)
  end do

end subroutine elens_postInput

! ================================================================================================ !
!  operations after linopt has been called
!  - compute r1, r2, rRef and sig_el out of optics, if needed;
!  - compute geometrical normalisation factors;
! ================================================================================================ !
subroutine elens_postLinopt

  use mathlib_bouncer
  use mod_utils, only : polinterp
  use mod_common, only : bez, kz, e0f, nucm0, fort3
  use mod_settings, only : st_debug
  use numerical_constants, only : c1e3, two
  use crcoall, only : lout, lerr

  integer j, jj, jguess
  real(kind=fPrec) oldVal
  logical lPrint

  do j=1,nelens
     
    ! find name of elens (for printout purposes)
    do jj=1,nele
      if(kz(jj)==elens_kz) then
        if (ielens(jj)==j) then
          exit
        end if
      end if
    end do

    if (elens_r1(j)<zero.or.elens_r2(j)<zero.or.elens_rref(j)<zero.or.elens_sig(j)<zero) then
      if (elens_iSet(j)==3 ) then
        elens_optVal(j)=elens_optVal(j)/real(elens_nUpdates(j),fPrec)
      elseif (elens_iSet(j)==4 ) then
        elens_optVal(j)=elens_optVal(j)/sqrt(real(elens_nUpdates(j),fPrec))
      end if
    end if

    lPrint=.false.
    ! compute R1 out of normalised settings
    if (elens_r1(j)<zero) then
      oldVal=elens_r1(j)
      if (elens_emin(j)>zero) then ! betatron cut
        elens_r1(j)=abs(elens_r1(j))*sqrt(elens_optVal(j)*elens_emin(j)/(e0f/nucm0))
      else ! momentum cut
        elens_r1(j)=abs(elens_r1(j))*(elens_optVal(j)*elens_sigdpp(j))
      end if
      elens_r1(j)=elens_r1(j)*c1e3
      ! ...and printout:
      write(lout,"(a,i0,a)") "ELENS> Recomputing R1 of e-lens #",j," named "//trim(bez(jj))//": "
      write(lout,"(a,1pe22.15)") "ELENS> - original value [sig]=",oldVal
      write(lout,"(a,1pe22.15)") "ELENS> - new value [mm]=",elens_r1(j)
      lPrint=.true.
    end if
    ! compute R2 out of normalised settings
    if (elens_r2(j)<zero) then
      oldVal=elens_r2(j)
      if (elens_emin(j)>zero) then ! betatron cut
        elens_r2(j)=abs(elens_r2(j))*sqrt(elens_optVal(j)*elens_emin(j)/(e0f/nucm0))
      else ! momentum cut
        elens_r2(j)=abs(elens_r2(j))*(elens_optVal(j)*elens_sigdpp(j))
      end if
      elens_r2(j)=elens_r2(j)*c1e3
      ! ...and printout:
      write(lout,"(a,i0,a)") "ELENS> Recomputing R2 of e-lens #",j," named "//trim(bez(jj))//": "
      write(lout,"(a,1pe22.15)") "ELENS> - original value [sig]=",oldVal
      write(lout,"(a,1pe22.15)") "ELENS> - new value [mm]=",elens_r2(j)
      lPrint=.true.
    end if
    ! compute Rref out of normalised settings
    if (elens_rref(j)<zero) then
      oldVal=elens_rref(j)
      if (elens_emin(j)>zero) then ! betatron cut
        elens_rref(j)=abs(elens_rref(j))*sqrt(elens_optVal(j)*elens_emin(j)/(e0f/nucm0))
      else ! momentum cut
        elens_rref(j)=abs(elens_rref(j))*(elens_optVal(j)*elens_sigdpp(j))
      end if
      elens_rref(j)=elens_rref(j)*c1e3
      ! ...and printout:
      write(lout,"(a,i0,a)") "ELENS> Recomputing Rref of e-lens #",j," named "//trim(bez(jj))//": "
      write(lout,"(a,1pe22.15)") "ELENS> - original value [sig]=",oldVal
      write(lout,"(a,1pe22.15)") "ELENS> - new value [mm]=",elens_rref(j)
      lPrint=.true.
    end if
    ! compute electron sigma out of normalised settings
    if (elens_sig(j)<zero) then
      oldVal=elens_sig(j)
      if (elens_emin(j)>zero) then ! betatron cut
        elens_sig(j)=abs(elens_sig(j))*sqrt(elens_optVal(j)*elens_emin(j)/(e0f/nucm0))
      else ! momentum cut
        elens_sig(j)=abs(elens_sig(j))*(elens_optVal(j)*elens_sigdpp(j))
      end if
      elens_sig(j)=elens_sig(j)*c1e3
      ! ...and printout:
      write(lout,"(a,i0,a)") "ELENS> Recomputing sigma of e-beam of e-lens #",j," named "//trim(bez(jj))//": "
      write(lout,"(a,1pe22.15)") "ELENS> - original value [sig]=",oldVal
      write(lout,"(a,1pe22.15)") "ELENS> - new value [mm]=",elens_sig(j)
      lPrint=.true.
    end if
    if(st_debug.and.lPrint) then
      write(lout,"(a)"         ) "ELENS> ...using the following parameters:"
      write(lout,"(a,1pe22.15)") "ELENS> - normalised emittance [m rad]=",elens_emin(j)
      write(lout,"(a,1pe22.15)") "ELENS> - geometrical emittance [m rad]=",elens_emin(j)/(e0f/nucm0)
      write(lout,"(a,1pe22.15)") "ELENS> - momentum of reference particle [MeV/c]=",e0f
      write(lout,"(a,1pe22.15)") "ELENS> - mass of reference particle [MeV/c2]=",nucm0
      write(lout,"(a,1pe22.15)") "ELENS> - rms of delta distributione []=",elens_sigdpp(j)
      write(lout,"(a,1pe22.15)") "ELENS> - beta/disp function [m]=",elens_optVal(j)
    end if

    ! check R1 and R2 against map
    if (elens_type(j)==3 ) then
      if ( elens_r1(j)< elens_radial_profile_R(0,elens_iRadial(j)) ) then
        write(lerr,"(a,i0,a)") "ELENS> ERROR on elens #",j, " named "//trim(bez(jj))//": "
        write(lerr,"(a)")      "ELENS>       R1 declared in "//trim(fort3)//" falls outside range of map" // &
             " contained in "//trim(elens_radial_filename(elens_iRadial(j)))
        write(lerr,"(2(a,1pe22.15),a)") "ELENS>       Rmin=",elens_radial_profile_R(0,elens_iRadial(j)), &
             " mm, Rmax=",elens_radial_profile_R(elens_radial_profile_nPoints(elens_iRadial(j)),elens_iRadial(j))," mm"
        call prror
      end if
      if ( elens_r2(j)> elens_radial_profile_R(elens_radial_profile_nPoints(elens_iRadial(j)),elens_iRadial(j)) ) then
        write(lerr,"(a,i0,a)") "ELENS> ERROR on elens #",j, " named "//trim(bez(jj))//": "
        write(lerr,"(a)")      "ELENS>       R2 declared in "//trim(fort3)//" falls outside range of map" // &
             " contained in "//trim(elens_radial_filename(elens_iRadial(j)))
        write(lerr,"(2(a,1pe22.15),a)") "ELENS>       Rmin=",elens_radial_profile_R(0,elens_iRadial(j)), &
             " mm, Rmax=",elens_radial_profile_R(elens_radial_profile_nPoints(elens_iRadial(j)),elens_iRadial(j))," mm"
        call prror
      end if
    end if

    ! compute geometrical normalisation factor (only if lens is not wire!)
    if ( .not.elens_lFull(j) .or. .not.elens_lZeroThick(j) ) then
      select case (elens_type(j))
      case(1) ! Uniform distribution
        elens_geo_norm(j) = (elens_r2(j)+elens_r1(j))*(elens_r2(j)-elens_r1(j))
      case(2) ! Gaussian distribution
        elens_geo_norm(j) = exp_mb(-(((elens_r1(j)/elens_sig(j))*(elens_r1(j)/elens_sig(j)))/two)) &
                           -exp_mb(-(((elens_r2(j)/elens_sig(j))*(elens_r2(j)/elens_sig(j)))/two))
      case(3) ! Radial profile
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
      end select ! case (elens_type(j))
      ! ...and printout:
      write(lout,"(a,i0,a,1pe22.15)") "ELENS> Geometrical normalisation factor for elens #",j, &
           " named "//trim(bez(jj))//": ",elens_geo_norm(j)
    end if

    ! Compute elens theta at R2, if requested by user
    call eLensTheta(j)

  end do

end subroutine elens_postLinopt

! ================================================================================================ !
!  Compute geometrical normalisation factor
! ================================================================================================ !
subroutine elens_setOptics(iElem, bAlpha, bBeta, bOrbit, bOrbitP, bDisp, bDispP)

  use parpro,     only : nele
  use mod_common, only : bez
  use numerical_constants, only : zero, two
  use mod_settings, only : st_debug
  use crcoall, only : lout

  integer,          intent(in) :: iElem
  real(kind=fPrec), intent(in) :: bAlpha(2)
  real(kind=fPrec), intent(in) :: bBeta(2)
  real(kind=fPrec), intent(in) :: bOrbit(2)
  real(kind=fPrec), intent(in) :: bOrbitP(2)
  real(kind=fPrec), intent(in) :: bDisp(2)
  real(kind=fPrec), intent(in) :: bDispP(2)

  if(iElem<1 .or. iElem>nele .or. &
       ( elens_r1  (ielens(iElem))>zero .and. elens_r2(ielens(iElem))>zero .and. &
         elens_rref(ielens(iElem))>zero .and. elens_sig(ielens(iElem))>zero ) ) return

  elens_nUpdates(ielens(iElem))=elens_nUpdates(ielens(iElem))+1

  if (elens_emin(ielens(iElem))>zero) then ! betatron cut
    select case (elens_iSet(ielens(iElem)))
    case(1) ! min
       if (elens_optVal(ielens(iElem))==zero) then ! first time
         elens_optVal(ielens(iElem))=min(bBeta(1),bBeta(2))
       else   
         elens_optVal(ielens(iElem))=min(bBeta(1),bBeta(2),elens_optVal(ielens(iElem)))
       end if
    case(2) ! max
       elens_optVal(ielens(iElem))=max(bBeta(1),bBeta(2),elens_optVal(ielens(iElem)))
    case(3) ! average
       elens_optVal(ielens(iElem))=elens_optVal(ielens(iElem))+(bBeta(1)+bBeta(2))/two
    case(4) ! quadratic average
       elens_optVal(ielens(iElem))=sqrt((elens_optVal(ielens(iElem))*elens_optVal(ielens(iElem)))+(bBeta(1)**2+bBeta(2)**2)/two)
    case(5) ! geometric average
       elens_optVal(ielens(iElem))=sqrt((elens_optVal(ielens(iElem))*elens_optVal(ielens(iElem)))*(bBeta(1)*bBeta(2)))
    end select ! case (elens_iSet(ielens(iElem)))
  else ! momentum cut
    select case (elens_iSet(ielens(iElem)))
    case(1) ! min
       if (elens_optVal(ielens(iElem))==zero) then ! first time
         elens_optVal(ielens(iElem))=min(abs(bDisp(1)),abs(bDisp(2)))
       else   
         elens_optVal(ielens(iElem))=min(abs(bDisp(1)),abs(bDisp(2)),elens_optVal(ielens(iElem)))
       end if
    case(2) ! max
       elens_optVal(ielens(iElem))=max(abs(bDisp(1)),abs(bDisp(2)),elens_optVal(ielens(iElem)))
    case(3) ! average
       elens_optVal(ielens(iElem))=elens_optVal(ielens(iElem))+(abs(bDisp(1))+abs(bDisp(2)))/two
    case(4) ! quadratic average
       elens_optVal(ielens(iElem))=sqrt((elens_optVal(ielens(iElem))*elens_optVal(ielens(iElem)))+(bDisp(1)**2+bDisp(2)**2)/two)
    case(5) ! geometric average
       elens_optVal(ielens(iElem))=sqrt((elens_optVal(ielens(iElem))*elens_optVal(ielens(iElem)))*(abs(bDisp(1))*abs(bDisp(2))))
    end select ! case (elens_iSet(ielens(iElem)))
  end if
 
  write(lout,"(a)") "ELENS> LinOpt for element: '"//trim(bez(iElem))//"'"
  write(lout,"(a,1pe22.15)") "ELENS> recorded value of beta/disp:",elens_optVal(ielens(iElem))
  
  if(st_debug) then
    write(lout,"(a,2(1x,f16.6))") "ELENS>        * Alpha X/Y:        ",bAlpha
    write(lout,"(a,2(1x,f16.6))") "ELENS>        * Beta X/Y:         ",bBeta
    write(lout,"(a,2(1x,f16.6))") "ELENS>        * Dispersion X/Y:   ",bDisp
    write(lout,"(a,2(1x,f16.6))") "ELENS>        * Dispersion XP/YP: ",bDispP
    write(lout,"(a,2(1x,f16.6))") "ELENS>        * Orbit X/Y:        ",bOrbit
    write(lout,"(a,2(1x,f16.6))") "ELENS>        * Orbit XP/YP:      ",bOrbitP
  end if
 
end subroutine elens_setOptics
  
! ================================================================================================ !
!  Compute eLens theta at rRef
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
  use physical_constants, only: clight, eps0
  use mod_common, only : e0, beta0, brho, bez, kz, zz0
  use mod_settings, only : st_quiet
  use crcoall, only : lout

  implicit none

  integer j,jj
  real(kind=fPrec) gamma_lens_beam

  if(elens_lThetaRref(j) .and. elens_lAllowUpdate(j)) then
    ! the update of elens_beta_lens_beam is not strictly needed here,
    !   apart from the case of elens_Ek is DYNK-ed
    gamma_lens_beam  = ((elens_Ek(j)*c1m3)/elens_beam_mass(j))+one ! from kinetic energy
    elens_beta_lens_beam(j) = sqrt((one+one/gamma_lens_beam)*(one-one/gamma_lens_beam))
    
    ! rRef: from mm to m (c1m3)
    ! theta: from rad to mrad (c1e3)
    elens_theta_ref(j) = one*((elens_len(j)*abs(elens_I(j)))/ &
         ((((two*pi)*((eps0*clight)*clight))*brho)*(elens_rref(j)*c1m3)))*c1e3
    elens_theta_ref(j) = sign(elens_theta_ref(j),elens_beam_chrg(j))
    elens_theta_ref(j) = elens_theta_ref(j)*(one/(elens_beta_lens_beam(j)*beta0) &
         -sign(one,elens_I(j)/elens_beam_chrg(j)))

    if(st_quiet < 2) then
      ! find name of elens
      do jj=1,nele
        if(kz(jj)==elens_kz) then
          if (ielens(jj)==j) then
            exit
          end if
        end if
      end do
      write(lout,"(a,i0,a,1pe22.15)") "ELENS> New theta at Rref for elens #",j, &
           " named "//trim(bez(jj))//": ",elens_theta_ref(j)
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
  use crcoall, only : lout, lerr
  use mod_common
  use mod_settings
  use string_tools
  use mod_units
  use mod_utils, only: checkArray
  use mod_alloc, only : alloc

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
  elens_radial_profile_R(:,ifile) = zero
  elens_radial_profile_J(:,ifile) = zero

  write(lout,"(a)") "ELENS> Parsing file with radial profile "//trim(elens_radial_filename(ifile))
  call f_requestUnit(elens_radial_filename(ifile),fUnit)
  call f_open(unit=fUnit,file=elens_radial_filename(ifile),mode='r',err=err,formatted=.true.,status="old")

10 continue
  read(fUnit,"(a)",end=20,err=30) inLine
  if(inLine(1:1) == "#") goto 10 ! comment line
  if(len(trim(inLine)) == 0 ) goto 10 ! empty line

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
  ii=ii+1
  if(ii>=size(elens_radial_profile_R,1)-1) then
    call alloc(elens_radial_profile_R, ii, ifile, zero, 'elens_radial_profile_R', 0, 1 )
    call alloc(elens_radial_profile_J, ii, ifile, zero, 'elens_radial_profile_J', 0, 1 )
  end if
  elens_radial_profile_nPoints(ifile) = ii
  elens_radial_profile_R(ii,ifile) = tmpR
  elens_radial_profile_J(ii,ifile) = tmpJ/c1e2 ! from [A/cm2] to [A/mm2]

  goto 10

20 continue

  call f_freeUnit(fUnit)
  if ( .not. checkArray(elens_radial_profile_R(:,ifile),elens_radial_profile_nPoints(ifile)) ) then
    write(lerr,"(a,i0,a,i0)") "ELENS> ERROR in radial profile"
    go to 30
  end if

  ! set current density at r(0) as at r(1), to properly compute the cumulative curve
  elens_radial_profile_J(0,ifile)=elens_radial_profile_J(1,ifile)

  write(lout,"(a,i0,a)") "ELENS> ...acquired ",elens_radial_profile_nPoints(ifile)," points."
  if(st_quiet < 2) then
    ! Echo parsed data (unless told to be quiet!)
    write(lout,"(a,i0,a)") "ELENS> Radial profile as from file "//&
      trim(elens_radial_filename(ifile))//" - #",ifile," - showing: ii, R[mm], J[A/mm2]"
    write(lout,"(a)") "ELENS> NB: point at ii=0 (i.e. R=0) added automatically by SixTrack"
    do ii=0,elens_radial_profile_nPoints(ifile)
      write(lout,"((a,i4),2(a,1pe22.15))") "ELENS> ",ii,",",elens_radial_profile_R(ii,ifile),",",elens_radial_profile_J(ii,ifile)
    end do
  end if
  return

30 continue
  write(lerr,"(a,i0,a)") "ELENS> ERROR ",iErr," while parsing file "//trim(elens_radial_filename(ifile))
  call prror

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
  use crcoall, only : lout
  use mod_utils, only: polintegrate
  use mod_alloc, only: alloc, dealloc
  use mod_settings, only: st_quiet

  implicit none

  integer, intent(in) :: ifile

  integer ii, nn
  real(kind=fPrec) tmpTot
  real(kind=fPrec), allocatable :: cumul(:)

  write(lout,"(a)") "ELENS> Normalising radial profile described in "//trim(elens_radial_filename(ifile))
  if(st_quiet < 2) flush(lout)
  nn=elens_radial_profile_nPoints(ifile)
  call alloc(cumul,nn+1,zero,'cumul')
  if(st_quiet < 2) then
    write(lout,"(a,i0,a)") "ELENS> Allocating cumul array to ",nn+1," elements"
    flush(lout)
  end if
  tmpTot=polintegrate(elens_radial_profile_R(0:nn,ifile), elens_radial_profile_J(0:nn,ifile), &
                      nn+1, elens_radial_mpoints(ifile), 2, cumul)
  elens_radial_profile_J(0:nn,ifile)=cumul(1:nn+1)
  call dealloc(cumul,'cumul')
  if(st_quiet < 2) then
    write(lout,"(a)") "ELENS> De-allocating cumul array"
    flush(lout)
  end if
  
  if(st_quiet < 2) then
    write(lout,"(a,i0)") "ELENS> Integrated radial profile read from file "//&
      trim(elens_radial_filename(ifile))//" - #",ifile
    write(lout,"(a)") "ELENS> NB: point at ii=0 (i.e. R=0) added automatically by SixTrack"
    do ii=0,elens_radial_profile_nPoints(ifile)
      write(lout,"((a,i4),2(a,1pe22.15))") "ELENS> ",ii,",",elens_radial_profile_R(ii,ifile),",",elens_radial_profile_J(ii,ifile)
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

  implicit none

  integer, intent(in) :: ifile

  integer ii

  do ii=0,elens_radial_profile_nPoints(ifile)
    elens_radial_profile_J(ii,ifile)=elens_radial_profile_J(ii,ifile)/&
                                     elens_radial_profile_J(elens_radial_profile_nPoints(ifile),ifile)
  end do

end subroutine normaliseRadialProfile

! ================================================================================================ !
!  Last modified: 2020-03-02
!  compute kick from electron lens
! ================================================================================================ !
subroutine elens_kick(i,ix,n)

  use mod_common, only : beta0, napx, bez
  use mod_common_main, only : xv1, xv2, yv1, yv2, moidpsv, rvv
  use mathlib_bouncer
  use numerical_constants, only : zero, one, two, c1m15, c1m7
  use mod_utils, only : polinterp
  use crcoall, only : lerr, lout

  implicit none
  
  integer, intent(in) :: i
  integer, intent(in) :: ix
  integer, intent(in) :: n
  
  real(kind=fPrec) xx, yy, rr, frr, tmpBB, epsilon, gteps, lteps, r1oSigSq, rroSigSq, cl2ori
  integer          jj, mPoints

  epsilon=c1m15
  gteps=one+epsilon
  lteps=one-epsilon
  cl2ori=c1m7

  ! for Gaussian distribution, prepare some useful numbers
  r1oSigSq=zero
  if ( elens_type(ielens(ix)) == 2 ) then
    r1oSigSq=(elens_r1(ielens(ix))/elens_sig(ielens(ix)))*(elens_r1(ielens(ix))/elens_sig(ielens(ix)))
  end if
  
  do jj=1,napx
     
    ! 1) apply offset of e-lens
    !    xx = x(proton) - elens_offset_x
    !    yy = y(proton) - elens_offset_y
    xx=xv1(jj)-elens_offset_x(ielens(ix))
    yy=xv2(jj)-elens_offset_y(ielens(ix))
    if (abs(xx)<epsilon) xx=zero
    if (abs(yy)<epsilon) yy=zero
    
    ! 2) calculate radius
    !    radial position of main beam relative to center of elens beam
    !    rr = sqrt(xx**2+yy**2)
    rr=sqrt((xx+yy)*(xx+yy)-two*(xx*yy))
    if (abs(rr)<epsilon) rr=zero
    
    ! 3) calculate kick
    !    shape function: spatial charge density depends on type:
    !    0    if r < R1
    !    frr  if R1 < r < R2
    !    1    if r >= R2
    if ( rr > elens_r1(ielens(ix))*lteps .or. elens_lFull(ielens(ix)) ) then ! rr<R1: no kick from elens
      if ( rr < elens_r2(ielens(ix))*lteps .and. .not.elens_lZeroThick(ielens(ix)) ) then ! R1<=rr<R2: type-dependent kick
          
        select case (elens_type(ielens(ix)))
            
        case (1) ! UNIFORM: eLens with uniform profile
          if ( elens_lFull(ielens(ix)) ) then
            ! formula: r/Rref
            frr=rr/elens_Rref(ielens(ix))
          else
            ! formula: (r^2-r1^2)/(r2^2-r1^2)*Rref/r
            frr=(((rr+elens_r1(ielens(ix)))*(rr-elens_r1(ielens(ix))))/elens_geo_norm(ielens(ix)))*(elens_rref(ielens(ix))/rr)
          end if
         
        case (2)
          ! GAUSSIAN: eLens with Gaussian profile
          if ( elens_lFull(ielens(ix)) ) then
            if ( rr <= cl2ori ) then
              ! formula: r*Rref/2sig^2/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))
              frr=(((rr/elens_sig(ielens(ix)))*(elens_rref(ielens(ix))/elens_sig(ielens(ix))))/two)/elens_geo_norm(ielens(ix))
            else
              ! formula: (1-exp(-r^2/2sig^2))/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))*Rref/r
              rroSigSq=(rr/elens_sig(ielens(ix)))*(rr/elens_sig(ielens(ix)))
              frr=((one-exp_mb(-(rroSigSq/two)))/elens_geo_norm(ielens(ix)))*(elens_rref(ielens(ix))/rr)
            end if     
          else
            ! formula: (exp(-r1^2/2sig^2)-exp(-r^2/2sig^2))/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))*Rref/r
            rroSigSq=(rr/elens_sig(ielens(ix)))*(rr/elens_sig(ielens(ix)))
            frr=((exp_mb(-(r1oSigSq/two))-exp_mb(-(rroSigSq/two)))/elens_geo_norm(ielens(ix)))*(elens_rref(ielens(ix))/rr)
          end if
           
        case (3)
          ! RADIAL PROFILE: eLens with radial profile as from file
          ! formula: (cumul_J(r)-cumul_J(r1))/(cumul_J(r2)-cumul_J(r1))*Rref/r
          mPoints=elens_radial_mpoints(ielens(ix))
          if ( elens_lFull(ielens(ix)) .and. rr <= cl2ori ) mPoints=elens_radial_mpoints_ori
          frr=polinterp( rr, &
                elens_radial_profile_R(0:elens_radial_profile_nPoints(elens_iRadial(ielens(ix))),elens_iRadial(ielens(ix))), &
                elens_radial_profile_J(0:elens_radial_profile_nPoints(elens_iRadial(ielens(ix))),elens_iRadial(ielens(ix))), &
                elens_radial_profile_nPoints(elens_iRadial(ielens(ix)))+1, &
                mPoints, elens_radial_jguess(ielens(ix)) )-elens_radial_fr1(ielens(ix))
          frr=(frr/elens_geo_norm(ielens(ix)))*(elens_rref(ielens(ix))/rr)
           
        case default
          write(lerr,"(a,i0,a)") "ELENS> ERROR elens_kick: elens_type=",elens_type(ielens(ix))," not recognized. "
          write(lerr,"(a)")      "ELENS>       Possible values for type are: 1, 2 and 3"
          call prror
          
        end select
       
      else ! rr>=R2
        if ( rr > elens_r2(ielens(ix))*gteps ) then
          ! rr>R2: formula: Rref/r
          frr=elens_rref(ielens(ix))/rr
        else
          if ( elens_lFull(ielens(ix)) .and. elens_lZeroThick(ielens(ix)) ) then
            ! degenerate e-lens to a wire: rr=R2=0.0!!!
            write(lerr,"(a,i0,a)") "ELENS> ERROR E-lens # ",ielens(ix)," named '"//trim(bez(ix))// &
                  "' degerate to a wire (R1=R2=0.0)"
            write(lerr,"(a)")      "ELENS>       Kick cannot be computed at wire position!"
            call prror
          else 
            ! rr=R2: formula: 1
            frr=one
          end if
        end if
        
      endif
     
      ! 'radial kick'
      frr = ( elens_theta_ref(ielens(ix))*frr )*moidpsv(jj)
      if(elens_lThetaRref(ielens(ix))) then
        tmpBB=sign(elens_beta_lens_beam(ielens(ix))*beta0,elens_I(ielens(ix))/elens_beam_chrg(ielens(ix)))
        frr = frr*((rvv(jj)-tmpBB)/(one-tmpBB))
      endif
     
      ! apply 'radial kick'
      yv1(jj)=yv1(jj)+frr*(xx/rr)
      yv2(jj)=yv2(jj)+frr*(yy/rr)
      
    endif ! rr<R1: no kick from elens
   
  end do ! do jj=1,napx

end subroutine elens_kick

! ================================================================================================ !
!  Last modified: 2020-03-02
!  compute kick from electron lens
!  inspired by wireda
! ================================================================================================ !
subroutine elens_kick_fox(i,ix)

  use mod_common, only : beta0, mtcda, bez
  use mod_settings, only : st_debug
  use crcoall, only : lout, lerr
  use mod_common_main
  use numerical_constants, only : zero, one, two, c1m15, c1m7
  use mod_utils, only : huntBin, polcof, polinterp
  use mod_lie_dab, only : lnv, idao, rscrri, iscrda
  use mod_common_da
  use mod_alloc, only: alloc, dealloc

  implicit none
  
  integer, intent(in) :: i
  integer, intent(in) :: ix
  
  integer          :: iLens, iRadial, nBin, nPoints, mPoints, kMin, kMax, kk
  real(kind=fPrec) :: rra, frra, xa, ya, xffset, yffset, ele_r1, ele_r2, ele_rr, elenor, elesig, elebet, tmpcof
  real(kind=fPrec) :: eletrf, epsilon, gteps, lteps, cl2ori
  real(kind=fPrec), allocatable :: cof(:)

! for FOX
  integer          :: idaa
  integer          :: hh(lnv)=0
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  save

!-----------------------------------------------------------------------
!FOX  B D ;
#include "include/dainicom.f90"
!FOX  D V DA INT XI     NORD NVAR ;
!FOX  D V DA INT YI     NORD NVAR ;
!FOX  D V DA INT RR     NORD NVAR ;
!FOX  D V DA INT FRR    NORD NVAR ;
!FOX  D V DA INT RR_SQ  NORD NVAR ;
!FOX  D V RE INT XFFSET ;
!FOX  D V RE INT YFFSET ;
!FOX  D V RE INT BETA0  ;
!FOX  D V RE INT ELE_R1 ;
!FOX  D V RE INT ELE_R2 ;
!FOX  D V RE INT ELE_RR ;
!FOX  D V RE INT ELETRF ;
!FOX  D V RE INT ELENOR ;
!FOX  D V RE INT ELESIG ;
!FOX  D V RE INT ELEBET ;
!FOX  D V DA INT TMPRR  NORD NVAR ;
!FOX  D V RE INT TMPCOF ;
!FOX  D V RE INT ONE    ;
!FOX  D V RE INT ZERO   ;
!FOX  D V RE INT TWO    ;
!FOX  D V RE INT HALF   ;
!FOX  D V RE INT RRA    ;
!FOX  D V RE INT FRRA   ;
!FOX  D V RE INT XA     ;
!FOX  D V RE INT YA     ;
!FOX  E D ;
!FOX  1 if(1==1) then
!-----------------------------------------------------------------------

  epsilon=c1m15
  gteps=one+epsilon
  lteps=one-epsilon
  cl2ori=c1m7

  iLens=ielens(ix)
  XFFSET=elens_offset_x(iLens)
  YFFSET=elens_offset_y(iLens)
  ELE_R1=elens_r1(iLens)
  ELE_R2=elens_r2(iLens)
  ELE_RR=elens_rref(iLens)
  ELETRF=elens_theta_ref(iLens)
  ELENOR=elens_geo_norm(iLens)
  ELESIG=elens_sig(iLens)
  ELEBET=elens_beta_lens_beam(iLens)
  
  FRRA=zero
  RRA=zero
  XA=zero
  YA=zero
  
!FOX    FRR=ZERO ;

  if (st_debug) then
    write(lout,'(2(a,i0))')'ELENS> ELENS_KICK_FOX for i=',i,' - ix=',ix
    write(lout,'(a)')'ELENS> ELENS_KICK_FOX closed orbit BEFORE elens:'
    call dapri(XX(1),6)
    call dapri(XX(2),6)
    call dapri(YY(1),6)
    call dapri(YY(2),6)
  end if
  
  ! 1) apply offset of e-lens
  !    xx = x(proton) - elens_offset_x
  !    yy = y(proton) - elens_offset_y
!FOX  XI=(XX(1)-XFFSET) ;
!FOX  YI=(XX(2)-YFFSET) ;
  call dapek(XI,hh,XA)
  call dapek(YI,hh,YA)
  if (abs(XA)<epsilon) then
!FOX  XI=ZERO ;
    XA=zero
  end if
  if (abs(YA)<epsilon) then
!FOX  YI=ZERO ;
    YA=zero
  end if
  if (st_debug) then
    write(lout,'(2(a,1pe23.16))')'ELENS> ELENS_KICK_FOX computing at XA=',XA,' - YA=',YA
  end if
  
  ! 2) calculate radius
  !    radial position of main beam relative to center of elens beam
  !    rr = sqrt(xx**2+yy**2)
!FOX  RR_SQ=(XI+YI)*(XI+YI)-(TWO*XI)*YI ;
  call dapek(RR_SQ,hh,RRA)
  if ( abs(RRA)>epsilon**2 ) then
!FOX  RR=SQRT(RR_SQ) ;
  else
!FOX  RR_SQ=ZERO ;
!FOX  RR=ZERO ;
  end if
  call dapek(RR,hh,RRA)
  if (st_debug) then
    write(lout,'(a,1pe23.16)')   'ELENS> ELENS_KICK_FOX computing at RRA=',RRA
    write(lout,'(2(a,1pe23.16))')'ELENS>                when R1=',elens_r1(iLens),' and R2=',elens_r2(iLens)
    flush(lout)
  end if

  ! 3) calculate kick
  !    shape function: spatial charge density depends on type:
  !    0    if r <= R1
  !    frr  if R1 < r < R2
  !    1    if r >= R2
  if ( RRA > elens_r1(iLens)*lteps .or. elens_lFull(iLens) ) then ! rr<R1: no kick from elens
    if ( RRA < elens_r2(iLens)*lteps.and..not.elens_lZeroThick(iLens) ) then ! R1<=rr<R2: type-dependent kick
       
      select case (elens_type(iLens))
        
      case (1) ! UNIFORM: eLens with uniform profile
        if ( elens_lFull(iLens) ) then
          ! formula: r/Rref
!FOX      FRR=RR/ELE_RR ;
          if (st_debug) write(lout,'(a)') "ELENS> ELENS_KICK_FOX: elens type UNIFORM, formula for full lens"
        else
          ! formula: (r^2-r1^2)/(r2^2-r1^2)*Rref/r
!FOX      FRR=(((RR+ELE_R1)*(RR-ELE_R1))/ELENOR)*(ELE_RR/RR) ;
          if (st_debug) write(lout,'(a)') "ELENS> ELENS_KICK_FOX: elens type UNIFORM, formula for hollow lens"
        end if
         
      case (2)
        ! GAUSSIAN: eLens with Gaussian profile
        if ( elens_lFull(iLens) ) then
          if ( RRA <= cl2ori ) then
            ! formula: r*Rref/2sig^2/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))
!FOX        FRR=(((RR/ELESIG)*(ELE_RR/ELESIG))/TWO)/ELENOR ;
            if (st_debug) write(lout,'(a)') "ELENS> ELENS_KICK_FOX: elens type GAUSSIAN, linearised formula"
          else
            ! formula: (1-exp(-r^2/2sig^2))/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))*Rref/r
!FOX        FRR=((ONE-EXP(-HALF*((RR/ELESIG)*(RR/ELESIG))))/ELENOR)*(ELE_RR/RR) ;
            if (st_debug) write(lout,'(a)') "ELENS> ELENS_KICK_FOX: elens type GAUSSIAN, formula for full lens"
          end if     
        else
          ! formula: (exp(-r1^2/2sig^2)-exp(-r^2/2sig^2))/(exp(-r1^2/2sig^2)-exp(-r2^2/2sig^2))*Rref/r
!FOX      FRR=((EXP(-HALF*((ELE_R1/ELESIG)*(ELE_R1/ELESIG)))-EXP(-HALF*((RR/ELESIG)*(RR/ELESIG))))/ELENOR)*(ELE_RR/RR) ;
          if (st_debug) write(lout,'(a)') "ELENS> ELENS_KICK_FOX: elens type GAUSSIAN, formula for hollow lens"
        end if
         
      case (3)
        ! RADIAL PROFILE: eLens with radial profile as from file
        ! formula: (cumul_J(r)-cumul_J(r1))/(cumul_J(r2)-cumul_J(r1))
        if (st_debug) write(lout,'(a)') "ELENS> ELENS_KICK_FOX: elens type RADIAL"
        iRadial=elens_iRadial(iLens)
        nPoints=elens_radial_profile_nPoints(iRadial)
        mPoints=elens_radial_mpoints(iLens)
        if ( elens_lFull(iLens) .and. RRA <= cl2ori ) then
          if (st_debug) write(lout,'(a,6(1X,i5))') "ELENS> ELENS_KICK_FOX: forcing mPoints to:", elens_radial_mpoints_ori
          mPoints=elens_radial_mpoints_ori
        end if
        nBin=huntBin(RRA,elens_radial_profile_R(0:nPoints,iRadial),nPoints+1,-1)-1
        kMin=min(max(nBin-(mPoints-1)/2,1),nPoints+2-mPoints)
        kMax=min(kMin+mPoints-1,nPoints+1)
        if (st_debug) write(lout,'(a,6(1X,i5))') "ELENS> ELENS_KICK_FOX: iRadial, nPoints, mPoints, nBin, kMin, kMax:", &
             iRadial, nPoints, mPoints, nBin, kMin, kMax
        call alloc(cof,kMax-kMin+1,zero,'cof')
        if (st_debug) then
          do kk=kMin,kMax
            write(lout,'(a,1X,i0,2(1X,1pe23.16))') "ELENS> ELENS_KICK_FOX: kk, R [mm], JJ [A/mm2]:", &
                kk, elens_radial_profile_R(kk,iRadial), elens_radial_profile_J(kk,iRadial)
          end do
        end if 
        call polcof(elens_radial_profile_R(kMin:kMax,iRadial),elens_radial_profile_J(kMin:kMax,iRadial),kMax-kMin+1,cof)
        if (st_debug) then
          do kk=1,mPoints
            write(lout,'(a,1X,i5,1X,1pe23.16)') "ELENS> ELENS_KICK_FOX: order, coefficient:", kk-1, cof(kk)
          end do
        end if 
        TMPCOF=COF(1)
!FOX    FRR=TMPCOF ;
!FOX    TMPRR=RR;
        do kk=2,kMax-kMin+1
          TMPCOF=COF(kk)
!FOX      FRR=FRR+(TMPRR*TMPCOF) ;
!FOX      TMPRR=TMPRR*RR ;
        end do
        if (st_debug) then
          call dapek(FRR,hh,FRRA)
          write(lout,'(a,1pe23.16)') "ELENS> ELENS_KICK_FOX: FRRA 1:", FRRA
        end if
!FOX    FRR=(FRR/ELENOR)*(ELE_RR/RR) ;
        if (st_debug) then
          call dapek(FRR,hh,FRRA)
          write(lout,'(a,1pe23.16)') "ELENS> ELENS_KICK_FOX: FRRA 2:", FRRA
        end if
        call dealloc(cof,'cof')
        
      case default
        write(lerr,"(a,i0,a)") "ELENS> ERROR in elens_kick_fox: elens_type=",elens_type(iLens)," not recognized. "
        write(lerr,"(a)")      "ELENS>       Possible values for type are: 1, 2 and 3"
        call prror
        
      end select

    else ! rr>=R2
      if ( rr > elens_r2(iLens)*gteps ) then
        ! rr>R2: formula: Rref/r
!FOX    FRR=ELE_RR/RR ;
      else
        if ( elens_lFull(iLens) .and. elens_lZeroThick(iLens) ) then
          ! degenerate e-lens to a wire: rr=R2=0.0!!!
          write(lerr,"(a,i0,a)") "ELENS> ERROR E-lens # ",iLens," named '"//trim(bez(ix))// &
                "' degerate to a wire (R1=R2=0.0)"
          write(lerr,"(a)")      "ELENS>       Kick cannot be computed at wire position!"
          call prror
        else 
          ! rr=RREF: formula: 1
!FOX      FRR=ONE ;
        end if
      end if
    endif
   
    ! 'radial kick'
!FOX    FRR=(ELETRF*FRR)*(MTCDA/(ONE+DPDA)) ;
    if(elens_lThetaRref(iLens)) then
      if(elens_I(iLens)/elens_beam_chrg(iLens) < zero) then
!FOX    FRR=FRR*((RV+ELEBET*BETA0)/(ONE+ELEBET*BETA0)) ;
      else
!FOX    FRR=FRR*((RV-ELEBET*BETA0)/(ONE-ELEBET*BETA0)) ;
      end if
    end if
    if (st_debug) then
      call dapek(FRR,hh,FRRA)
      write(lout,'(a,1pe23.16)') "ELENS> ELENS_KICK_FOX: FRRA 3:", FRRA
    end if
    
!FOX  YY(1)=YY(1)+(FRR*XI)/RR ;
!FOX  YY(2)=YY(2)+(FRR*YI)/RR ;
  end if

  if (st_debug) then
    call dapek(FRR,hh,FRRA)
    write(lout,'(2(a,1pe23.16))')'ELENS> ELENS_KICK_FOX computed at RRA= ',RRA,' - FRRA= ',FRRA
    write(lout,'(a)')'ELENS> ELENS_KICK_FOX closed orbit AFTER elens:'
    call dapri(XX(1),6)
    call dapri(XX(2),6)
    call dapri(YY(1),6)
    call dapri(YY(2),6)
  end if
  
! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION

end subroutine elens_kick_fox

#ifdef CR
subroutine elens_crcheck(fileUnit,readErr)

  use crcoall

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: readErr

  integer j

  read(fileUnit,err=10,end=10) (elens_lAllowUpdate_CR(j), j=1, nelens)

  readErr = .false.
  return

10 continue

  write(lerr, "(a,i0,a)") "CR_CHECK> ERROR Reading C/R file fort.",fileUnit," in ELENS"
  write(crlog,"(a,i0,a)") "CR_CHECK> ERROR Reading C/R file fort.",fileUnit," in ELENS"
  flush(crlog)
  readErr = .true.

end subroutine elens_crcheck

subroutine elens_crpoint(fileUnit, writeErr)

  use crcoall

  integer, intent(in)  :: fileUnit
  logical, intent(out) :: writeErr

  integer j

  write(fileunit,err=10) (elens_lAllowUpdate(j), j=1, nelens)
  flush(fileunit)

  writeErr = .false.
  return

10 continue

  write(lerr, "(a,i0,a)") "CR_POINT> ERROR Writing C/R file fort.",fileUnit," in ELENS"
  write(crlog,"(a,i0,a)") "CR_POINT> ERROR Writing C/R file fort.",fileUnit," in ELENS"
  flush(crlog)
  writeErr = .true.

end subroutine elens_crpoint

subroutine elens_crstart
  elens_lAllowUpdate(1:nelens) = elens_lAllowUpdate_CR(1:nelens)
end subroutine elens_crstart

#endif

end module elens
