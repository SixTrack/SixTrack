module cheby
  use floatPrecision
  use numerical_constants, only : zero, one
  implicit none
  private
  public :: &
       ! specific to allocate arrays
       cheby_allocate_arrays, cheby_expand_arrays, &
       ! specific to FOX
       cheby_lFox, icheby, cheby_kz, cheby_ktrack, cheby_kick, cheby_kick_fox, &
       ! specific to input parsing
       cheby_parseInputLine, cheby_parseInputDone, cheby_postInput, &
       ! specific to DYNK
       cheby_resetI, cheby_setScaleKick, cheby_I

  ! A.Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 02-04-2020
  ! module for handling maps expressed by Chebyshev polynomials

  integer, allocatable, save  :: icheby(:)              ! index of chebyshev lens
  integer, save               :: ncheby=0               ! number of chebyshev lenses actually in memory
  integer, save               :: ncheby_mapEchoes=0     ! number of requested echoes of maps
  integer, save               :: ncheby_tables=0        ! number of chebyshev tables in memory
  integer, parameter          :: cheby_kz=42            ! kz of chebyshev lenses
  integer, parameter          :: cheby_ktrack=67        ! ktrack of chebyshev lenses
  integer, parameter          :: cheby_maxNterms=21     ! max number of terms in Chebyshev polynomials (required by FOX)
                                                        ! please keep this synched with TX, TY, TPX and TPY
                                                        !    in the FOX part of this module
  logical, save               :: cheby_lFox_def=.true.  ! default lFox

  ! variables to save parameters for tracking etc.
  integer,          allocatable, save :: cheby_itable(:)      ! index of chebyshev table
  real(kind=fPrec), allocatable, save :: cheby_r2(:)          ! outer radius R2 [mm] (optional)
  real(kind=fPrec), allocatable, save :: cheby_r1(:)          ! inner radius R1 [mm] (optional)
  real(kind=fPrec), allocatable, save :: cheby_angle(:)       ! rotation angle about the longitudinal axis [deg] (optional)
  real(kind=fPrec), allocatable, save :: cheby_offset_x(:)    ! hor. offset [mm] (optional)
  real(kind=fPrec), allocatable, save :: cheby_offset_y(:)    ! ver. offset [mm] (optional)
  real(kind=fPrec), allocatable, save :: cheby_I(:)           ! actual powering of lens [A] (optional)
  real(kind=fPrec), allocatable, save :: cheby_scalingFact(:) ! scaling factor [] (computed internally)
  logical,          allocatable, save :: cheby_lFox(:)        ! the kick from the chebyshev map should be taken into account
                                                              !   when searching for the closed orbit with DA algebra  
  ! show map
  integer,          allocatable, save :: cheby_iLens(:)       ! lens for which a map echo is requested
  character(len=:), allocatable, save :: cheby_mapFileName(:) ! file name
  real(kind=fPrec), allocatable, save :: cheby_mapXmin(:)     ! xMin [mm]
  real(kind=fPrec), allocatable, save :: cheby_mapXmax(:)     ! xMax [mm]
  integer,          allocatable, save :: cheby_mapNx(:)       ! number of intervals
  real(kind=fPrec), allocatable, save :: cheby_mapYmin(:)     ! yMin [mm]
  real(kind=fPrec), allocatable, save :: cheby_mapYmax(:)     ! yMax [mm]
  integer,          allocatable, save :: cheby_mapNy(:)       ! number of intervals

  ! tables with chebyshev coefficients
  character(len=:), allocatable, save :: cheby_filename(:)    ! file names
  real(kind=fPrec), allocatable, save :: cheby_coeffs(:,:,:)  ! coefficients
  integer,          allocatable, save :: cheby_maxOrder(:)    ! max order of the current map
  real(kind=fPrec), allocatable, save :: cheby_refI(:)        ! reference current [A] (optional)
  real(kind=fPrec), allocatable, save :: cheby_refR(:)        ! reference radius [mm] (mandatory)

contains

subroutine cheby_allocate_arrays
  use mod_alloc, only : alloc
  use parpro, only : nele
  implicit none
  call alloc(icheby,nele,0,'icheby')
end subroutine cheby_allocate_arrays

subroutine cheby_expand_arrays(nele_new)
  use mod_alloc, only : alloc
  implicit none
  integer, intent(in) :: nele_new
  call alloc(icheby,nele_new,0,'icheby')
end subroutine cheby_expand_arrays

subroutine cheby_expand_arrays_lenses(ncheby_new)
  use mod_alloc, only : alloc
  implicit none
  integer, intent(in) :: ncheby_new
  ! chebyshev lens charachteristics
  call alloc(cheby_itable     ,            ncheby_new,               0, 'cheby_itable'     )
  call alloc(cheby_r2         ,            ncheby_new,            zero, 'cheby_r2'         )
  call alloc(cheby_r1         ,            ncheby_new,            zero, 'cheby_r1'         )
  call alloc(cheby_angle      ,            ncheby_new,            zero, 'cheby_angle'      )
  call alloc(cheby_offset_x   ,            ncheby_new,            zero, 'cheby_offset_x'   )
  call alloc(cheby_offset_y   ,            ncheby_new,            zero, 'cheby_offset_y'   )
  call alloc(cheby_I          ,            ncheby_new,            -one, 'cheby_I'          )
  call alloc(cheby_scalingFact,            ncheby_new,             one, 'cheby_scalingFact')
  call alloc(cheby_lFox       ,            ncheby_new,  cheby_lFox_def, 'cheby_lFox'       )
end subroutine cheby_expand_arrays_lenses

subroutine cheby_expand_arrays_map_echo(ncheby_mapEchoes_new)
  use parpro, only : mFileName
  use mod_alloc, only : alloc
  implicit none
  integer, intent(in) :: ncheby_mapEchoes_new
  ! map
  call alloc(cheby_iLens      ,            ncheby_mapEchoes_new,              0, 'cheby_iLens'      )
  call alloc(cheby_mapFileName, mFileName, ncheby_mapEchoes_new,            " ", 'cheby_mapFileName')
  call alloc(cheby_mapXmin    ,            ncheby_mapEchoes_new, -5.0e+00_fPrec, 'cheby_mapXmin'    )
  call alloc(cheby_mapXmax    ,            ncheby_mapEchoes_new,  5.0e+00_fPrec, 'cheby_mapXmax'    )
  call alloc(cheby_mapNx      ,            ncheby_mapEchoes_new,            100, 'cheby_mapNx'      )
  call alloc(cheby_mapYmin    ,            ncheby_mapEchoes_new, -5.0e+00_fPrec, 'cheby_mapYmin'    )
  call alloc(cheby_mapYmax    ,            ncheby_mapEchoes_new,  5.0e+00_fPrec, 'cheby_mapYmax'    )
  call alloc(cheby_mapNy      ,            ncheby_mapEchoes_new,            100, 'cheby_mapNy'      )
end subroutine cheby_expand_arrays_map_echo

subroutine cheby_expand_arrays_tables(ncheby_tables_new)
  use parpro, only : mFileName
  use mod_alloc, only : alloc
  implicit none
  integer, intent(in) :: ncheby_tables_new
  call alloc(cheby_filename   , mFileName,       ncheby_tables_new,  " ", 'cheby_filename'          )
  call alloc(cheby_coeffs     ,            0, 0, ncheby_tables_new, zero, 'cheby_coeffs'  , 0, 0, 1 )
  call alloc(cheby_maxOrder   ,                  ncheby_tables_new,    0, 'cheby_maxOrder'          )
  call alloc(cheby_refI       ,                  ncheby_tables_new,  one, 'cheby_refI'              )
  call alloc(cheby_refR       ,                  ncheby_tables_new, zero, 'cheby_refR'              )
end subroutine cheby_expand_arrays_tables

! ================================================================================================ !
!  Parse Line for Chebyshev lens
!  Last modified: 2019-02-25
! ================================================================================================ !
subroutine cheby_parseInputLine(inLine, iLine, iErr)

  use crcoall, only : lerr, lout
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
  integer nSplit, iElem, j, chIdx, tmpi1
  logical spErr, tmpl

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "CHEBY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(lnSplit(1))

  case("SHOW")
    if(nSplit < 3) then
      write(lerr,"(a,i0)") "CHEBY> ERROR Expected at least 3 input parameters, got ",nSplit
      iErr = .true.
      return
    end if

    iElem = -1
    tmpch = trim(lnSplit(2))
    do j=1,nele
      if(bez(j)==tmpch.and.icheby(j)/=0) then
        iElem = j
        exit
      end if
    end do
    if(iElem == -1) then
      write(lerr,"(a)") "CHEBY> ERROR Element '"//trim(lnSplit(2))//"' not a Chebyshev lens."
      write(lerr,"(a)") "CHEBY>       Either you mis-typed the element name or"
      write(lerr,"(a)") "CHEBY>       SHOW line comes before the declaration of the lens"
      iErr = .true.
      return
    end if
    do tmpi1=1,ncheby_mapEchoes
      if (cheby_mapFileName(tmpi1)==lnSplit(3)) then
        write(lerr,"(a)") "CHEBY> ERROR File '"//trim(lnSplit(3))//"' already in use."
        iErr = .true.
        return
      end if
    end do
    ncheby_mapEchoes = ncheby_mapEchoes+1
    call cheby_expand_arrays_map_echo(ncheby_mapEchoes)
    cheby_mapFileName(ncheby_mapEchoes) = trim(lnSplit(3))
    cheby_iLens(ncheby_mapEchoes) = icheby(iElem)
    if ( nSplit > 3 ) then
      if ( nSplit < 9 ) then
        write(lerr,"(a,i0)") "CHEBY> ERROR Expected 8 input parameters, got ",nSplit
        write(lerr,"(a)")    "CHEBY> format of SHOW line:"
        write(lerr,"(a)")    "SHOW name filename xmin[mm] xmax[mm] nx ymin[mm] ymax[mm] ny"
        iErr = .true.
        return
      endif
      call chr_cast(lnSplit(4),cheby_mapXmin(ncheby_mapEchoes),iErr)
      call chr_cast(lnSplit(5),cheby_mapXmax(ncheby_mapEchoes),iErr)
      call chr_cast(lnSplit(6),cheby_mapNx  (ncheby_mapEchoes),iErr)
      call chr_cast(lnSplit(7),cheby_mapYmin(ncheby_mapEchoes),iErr)
      call chr_cast(lnSplit(8),cheby_mapYmax(ncheby_mapEchoes),iErr)
      call chr_cast(lnSplit(9),cheby_mapNy  (ncheby_mapEchoes),iErr)
    end if

    if(st_debug) then
      call sixin_echoVal("name",trim(bez(iElem)),                              "CHEBY",iLine)
      call sixin_echoVal("filename", trim(cheby_mapFileName(ncheby_mapEchoes)),"CHEBY",iLine)
      call sixin_echoVal("xmin [mm]",cheby_mapXmin(ncheby_mapEchoes),          "CHEBY",iLine)
      call sixin_echoVal("xmax [mm]",cheby_mapXmax(ncheby_mapEchoes),          "CHEBY",iLine)
      call sixin_echoVal("Nx     []",cheby_mapNx  (ncheby_mapEchoes),          "CHEBY",iLine)
      call sixin_echoVal("ymin [mm]",cheby_mapYmin(ncheby_mapEchoes),          "CHEBY",iLine)
      call sixin_echoVal("ymin [mm]",cheby_mapYmax(ncheby_mapEchoes),          "CHEBY",iLine)
      call sixin_echoVal("Ny     []",cheby_mapNy  (ncheby_mapEchoes),          "CHEBY",iLine)
    end if

  case("FOX")
    if(nSplit<2 .or. nSplit>3) then
      write(lerr,"(a,i0)") "CHEBY> ERROR Expected 1 or 2 input parameters for FOX line, got ",nSplit-1
      write(lerr,"(a)")    "CHEBY>       example:     FOX  on|off|true|false (ALL|BEF(ORE)|AFT(ER))"
      iErr = .true.
      return
    end if
   
    call chr_cast(lnSPlit(2), tmpl,iErr)
    if (ncheby>0) cheby_lFox(ncheby)=tmpl

    if (nSplit>=3) then
      select case (chr_toLower(trim(lnSplit(3))))
      case('all')
        do tmpi1=1,ncheby-1
          cheby_lFox(tmpi1) = tmpl
        end do
        cheby_lFox_def=tmpl
        if(st_debug) write(lout,"(a)") "CHEBY> Setting lFox as read to all chebyshev lenses"
      case('bef','before')
        do tmpi1=1,ncheby-1
          cheby_lFox(tmpi1) = tmpl
        end do
        if(st_debug) write(lout,"(a)") "CHEBY> Setting lFox as read to all chebyshev lenses "// &
             "declared before the current FOX line"
      case('aft','after')
        cheby_lFox_def=tmpl
        if(st_debug) write(lout,"(a)") "CHEBY> Setting lFox as read to all chebyshev lenses "// &
             "declared after the current FOX line"
      case default
        write(lerr,"(a)") "CHEBY> ERROR Unidentified third parameter of FOX line, got: '"//trim(lnSplit(3))//"'"
        write(lerr,"(a)") "CHEBY>       example:     FOX  on|off|true|false (ALL|BEF(ORE)|AFT(ER))"
        iErr = .true.
        return
      end select ! case (lnSplit(3))
    end if
    
    if(st_debug) then
      call sixin_echoVal("fox",tmpl,"CHEBY",iLine)
    end if
     
  case default

    if(nSplit < 2) then
      write(lerr,"(a,i0)") "CHEBY> ERROR Expected at least 2 input parameters, got ",nSplit
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
      write(lerr,"(a)") "CHEBY> ERROR Element '"//trim(lnSplit(1))//"' not found in single element list."
      iErr = .true.
      return
    end if

    if(kz(iElem) /= cheby_kz) then
      write(lerr,"(3(a,i0))") "CHEBY> ERROR Element type is kz(",iElem,") = ",kz(iElem)," != ",cheby_kz
      iErr = .true.
      return
    end if
    if(el(iElem) /= zero .or. ek(iElem) /= zero .or. ed(iElem) /= zero) then
      write(lerr,"(a)")       "CHEBY> ERROR Length el(iElem) (Chebyshev lens is treated as thin element), "//&
        "and first and second field have to be zero:"
      write(lerr,"(2(a,i0),a)") "CHEBY>       el(",iElem,") = ",el(iElem)," != 0"
      write(lerr,"(2(a,i0),a)") "CHEBY>       ed(",iElem,") = ",ed(iElem)," != 0"
      write(lerr,"(2(a,i0),a)") "CHEBY>       ek(",iElem,") = ",ek(iElem)," != 0"
      iErr = .true.
      return
    end if

    if(icheby(iElem) /= 0) then
      write(lerr,"(a)") "CHEBY> ERROR The element '"//trim(bez(iElem))//"' was defined twice."
      iErr = .true.
      return
    end if
    ncheby = ncheby+1
    call cheby_expand_arrays_lenses(ncheby)
    icheby(iElem) = ncheby

    ! File with Chebyshev polynomials
    tmpch = trim(lnSplit(2))
    ! Check if profile has already been requested:
    chIdx = -1
    do tmpi1=1,ncheby_tables
      if(trim(tmpch) == trim(cheby_filename(tmpi1))) then
        cheby_itable(icheby(iElem)) = tmpi1
        chIdx = tmpi1
        exit
      end if
    end do
    if(chIdx == -1) then
      ! Unsuccessful search
      ncheby_tables = ncheby_tables+1
      if(st_debug) write(lout,"(a,i0,a)") "CHEBY> adding table" ,ncheby_tables,"'"//trim(tmpch)//"'"
      call cheby_expand_arrays_tables(ncheby_tables)
      cheby_itable(icheby(iElem)) = ncheby_tables
      cheby_filename(tmpi1) = trim(tmpch)
    end if

    ! Additional geometrical infos:
    if(nSplit >= 3) call chr_cast(lnSplit(3),cheby_r2(icheby(iElem)),iErr)
    if(nSplit >= 4) call chr_cast(lnSplit(4),cheby_r1(icheby(iElem)),iErr)
    if(nSplit >= 5) call chr_cast(lnSplit(5),cheby_angle(icheby(iElem)),iErr)
    if(nSplit >= 6) call chr_cast(lnSplit(6),cheby_offset_x(icheby(iElem)),iErr)
    if(nSplit >= 7) call chr_cast(lnSplit(7),cheby_offset_y(icheby(iElem)),iErr)
    if(nSplit >= 8) call chr_cast(lnSplit(8),cheby_I(icheby(iElem)),iErr)

    if(st_debug) then
      call sixin_echoVal("name",    trim(bez(iElem)),                                   "CHEBY",iLine)
      call sixin_echoVal("filename",trim(cheby_filename(cheby_itable(icheby(iElem)))),  "CHEBY",iLine)
      if(nSplit >= 3) call sixin_echoVal("r2 [mm]"      , cheby_r2(icheby(iElem)),      "CHEBY",iLine)
      if(nSplit >= 4) call sixin_echoVal("r1 [mm]"      , cheby_r1(icheby(iElem)),      "CHEBY",iLine)
      if(nSplit >= 5) call sixin_echoVal("tilt [deg]"   , cheby_angle(icheby(iElem)),   "CHEBY",iLine)
      if(nSplit >= 6) call sixin_echoVal("offset_x [mm]", cheby_offset_x(icheby(iElem)),"CHEBY",iLine)
      if(nSplit >= 7) call sixin_echoVal("offset_y [mm]", cheby_offset_y(icheby(iElem)),"CHEBY",iLine)
      if(nSplit >= 8) call sixin_echoVal("I [A]"        , cheby_I(icheby(iElem)),       "CHEBY",iLine)
    end if

  end select

end subroutine cheby_parseInputLine


subroutine cheby_parseInputDone(iErr)

  use mod_common, only : kz,bez

  implicit none

  logical, intent(inout) :: iErr

end subroutine cheby_parseInputDone


subroutine cheby_postInput

  use crcoall, only : lout, lerr
  use parpro, only : nele
  use mod_common, only : kz,bez,fort3
  use mod_settings, only : st_quiet

  integer ii, jj, kk, ncheb
  logical exist
  real(kind=fPrec) tmpFlt

  ! Check that all chebyshev lenses in fort.2 have a corresponding declaration in fort.3
  ncheb=0
  do jj=1,nele
    if(kz(jj)==cheby_kz) then
      if (icheby(jj).eq.0) then
        write(lerr,"(a,i0,a)") "CHEBY> ERROR single element ",jj," named '"//trim(bez(jj))//"'"
        write(lerr,"(a)")      "CHEBY>       does not have a corresponding line in CHEB block in "//trim(fort3)
        call prror
      else
        ncheb=ncheb+1
      end if
    end if
  end do
  if ( ncheb.ne.ncheby ) then
    write(lerr,"(a,i0)") "CHEBY> ERROR number of chebyshev lenses declared in CHEB block in "//trim(fort3)//" ",ncheby
    write(lerr,"(a,i0)") "CHEBY>       is not the same as the total number of chebyshev lenses in lattice ",ncheb
    call prror
  end if

  ! Parse files with coefficients for Chebyshev polynomials:
  do jj=1,ncheby_tables
    inquire(file=cheby_filename(jj), exist=exist)
    if(.not. exist) then
      write(lerr,"(a)") "CHEBY> ERROR Problems with file with coefficients for Chebyshev polynominals: ", &
            trim(cheby_filename(jj))
      call prror
    end if
    call parseChebyFile(jj)
  end do

  ! finalise setting-up of chebyshev lenses
  do jj=1,ncheby

    ! find name, to get ready for error messages
    do kk=1,nele
      if(kz(kk) == cheby_kz .and. icheby(kk) == jj ) then
        exit
      end if
    end do

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
      write(lerr,"(a)")      "CHEBY> ERROR R2 cannot be larger than domain of Chebyshev polynomials!"
      write(lerr,"(a,1pe22.15,a,1pe22.15)") "CHEBY>       R2 [mm]: ",cheby_r2(jj), &
           " - reference radius [mm]:",cheby_refR(cheby_itable(jj))
      goto 10
    end if
    if (cheby_r1(jj)<zero) then
      write(lerr,"(a)")      "CHEBY> ERROR R1 cannot be lower than zero!"
      goto 10
    end if
    if (cheby_r2(jj)<zero) then
      ! a sanity check, most probably will never be triggered, but better stay safe
      write(lerr,"(a)")      "CHEBY> ERROR R2 cannot be lower than zero!"
      goto 10
    end if
    if (cheby_I(jj)<=zero) then
      call cheby_resetI( jj )
    else
      call cheby_setScaleKick(jj)
    end if

    ! checks on maps
    do ii=1,ncheby_mapEchoes
      if ( cheby_iLens(ii) .eq. jj) then
        if (cheby_mapXmax(ii)<cheby_mapXmin(ii)) then
          ! swap
          tmpFlt=cheby_mapXmax(ii)
          cheby_mapXmax(ii)=cheby_mapXmin(ii)
          cheby_mapXmin(ii)=tmpFlt
        end if
        if (cheby_mapXmax(ii)==cheby_mapXmin(ii)) then
          write(lerr,"(a)")                     "CHEBY> ERROR X-extremes for map coincide!"
          write(lerr,"(a,1pe22.15,a,1pe22.15)") "CHEBY>       xmin [mm]: ",cheby_mapXmin(ii), &
                                                          " - xmax [mm]: ",cheby_mapXmax(ii)
          goto 10
        end if
        if (cheby_mapNx(ii)<1) then
          write(lerr,"(a)")    "CHEBY> ERROR wrong X-stepping for map!"
          write(lerr,"(a,i0)") "CHEBY>       must be >0 - got: ",cheby_mapNx(ii)
          goto 10
        end if
        if (cheby_mapYmax(ii)<cheby_mapYmin(ii)) then
          ! swap
          tmpFlt=cheby_mapYmax(ii)
          cheby_mapYmax(ii)=cheby_mapYmin(ii)
          cheby_mapYmin(ii)=tmpFlt
        end if
        if (cheby_mapYmax(ii)==cheby_mapYmin(ii)) then
          write(lerr,"(a)")                     "CHEBY> ERROR Y-extremes for map coincide!"
          write(lerr,"(a,1pe22.15,a,1pe22.15)") "CHEBY>       ymin [mm]: ",cheby_mapYmin(ii), &
                                                          " - ymax [mm]: ",cheby_mapYmax(ii)
          goto 10
        end if
        if (cheby_mapNy(ii)<1) then
          write(lerr,"(a)")    "CHEBY> ERROR wrong Y-stepping for map!"
          write(lerr,"(a,i0)") "CHEBY>       must be >0 - got: ",cheby_mapNy(ii)
          goto 10
        end if
      end if
    end do

    if(st_quiet < 2) then
      write(lout,"(a)") ''
      write(lout,"(a,i0,a)")       "CHEBY> status of chebyshev lens #",jj," - name: '"//trim(bez(kk))//"'"
      write(lout,"(a)")            "CHEBY> - filename         : '"//trim(cheby_filename(cheby_itable(jj)))//"'"
      write(lout,"(a,1pe22.15)")   "CHEBY> - R2           [mm]: ",cheby_r2(jj)
      write(lout,"(a,1pe22.15)")   "CHEBY> - R1           [mm]: ",cheby_r1(jj)
      write(lout,"(a,1pe22.15)")   "CHEBY> - tilt angle  [deg]: ",cheby_angle(jj)
      write(lout,"(a,1pe22.15)")   "CHEBY> - hor offset   [mm]: ",cheby_offset_x(jj)
      write(lout,"(a,1pe22.15)")   "CHEBY> - ver offset   [mm]: ",cheby_offset_y(jj)
      write(lout,"(a,1pe22.15)")   "CHEBY> - lens powering [A]: ",cheby_I(jj)
      write(lout,"(a,1pe22.15)")   "CHEBY> - scaling factor []: ",cheby_scalingFact(jj)
      do ii=1,ncheby_mapEchoes
        if ( cheby_iLens(ii) .eq. jj) then
          write(lout,"(a)")          "CHEBY> requested dump of potential map:"
          write(lout,"(a)")          "CHEBY> - map     filename : '"//trim(cheby_mapFileName(ii))//"'"
          write(lout,"(a,1pe22.15)") "CHEBY> - xmin         [mm]: ",cheby_mapXmin(ii)
          write(lout,"(a,1pe22.15)") "CHEBY> - xmax         [mm]: ",cheby_mapXmax(ii)
          write(lout,"(a,i0)")       "CHEBY> - Nx             []: ",cheby_mapNx  (ii)
          write(lout,"(a,1pe22.15)") "CHEBY> - ymin         [mm]: ",cheby_mapYmin(ii)
          write(lout,"(a,1pe22.15)") "CHEBY> - ymin         [mm]: ",cheby_mapYmax(ii)
          write(lout,"(a,i0)")       "CHEBY> - Ny             []: ",cheby_mapNy  (ii)
          call cheby_potentialMap(jj,kk)
        end if
      end do
    end if
  end do

  return

10 continue
   write(lout,"(a,i0,a)") "CHEBY>       concerned lens #", jj," - name: '"//trim(bez(kk))//"'"
   call prror

end subroutine cheby_postInput


! ================================================================================================ !
!  Last modified: 2020-04-07
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
!  or in tabular format:
!  TABLE : y   |   0    1    2    3  ...
!  TABLE : x 0 | 0.0  1.2  2.4 -3.6  ...
! ================================================================================================ !
subroutine parseChebyFile(ifile)

  use crcoall, only : lout, lerr
  use mod_alloc, only : alloc, dealloc
  use mod_common
  use mod_settings
  use string_tools
  use mod_units

  implicit none

  integer, intent(in) :: ifile

  character(len=:), allocatable   :: lnSplit(:)
  integer, allocatable :: jOrder(:)
  character(len=mInputLn) inLine
  integer nSplit, ii, jj, fUnit, iDim, jDim, jMax
  logical spErr, err, lDefI

  lDefI=.true.

  write(lout,"(a)") "CHEBY> Parsing file with coefficients for Chebyshev polynomials '"//trim(cheby_filename(ifile))//"'..."
  call f_requestUnit(cheby_filename(ifile),fUnit)
  call f_open(unit=fUnit,file=cheby_filename(ifile),mode='r',err=err,formatted=.true.,status="old")
  if(err) then
    write(lerr,"(a)") "CHEBY> ERROR Failed to open file."
    goto 40
  end if

10 continue
  read(fUnit,"(a)",end=20,err=30) inLine
  if(inLine(1:1) == "#") goto 10

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "CHEBY> ERROR Failed to split Chebyshev input line"
    goto 30
  end if
 
  if(inLine(1:1) == "I") then
    ! Read reference current of lens
    if(nSplit < 3) then
      write(lerr,"(a)") "CHEBY> ERROR Not enough arguments for expressing ref lens current [A]."
      write(lerr,"(a)") "CHEBY>       Correct format:"
      write(lerr,"(a)") "CHEBY>       I : value"
      goto 30
    end if
    call chr_cast(lnSplit(3),cheby_refI(ifile),spErr)
    if(spErr) then
      write(lerr,"(a)") "CHEBY> ERROR in casting ref lens current: "//trim(lnSplit(3))
      goto 30
    end if
    lDefI=.false.

  else if(inLine(1:1) == "R") then
    ! Read reference radius e-beam in e-lens
    if(nSplit < 3) then
      write(lerr,"(a)") "CHEBY> ERROR Not enough arguments for expressing ref lens radius [mm]."
      write(lerr,"(a)") "CHEBY>       Correct format:"
      write(lerr,"(a)") "CHEBY>       R : value"
      goto 30
    end if
    call chr_cast(lnSplit(3),cheby_refR(ifile),spErr)
    if(spErr) then
      write(lerr,"(a)") "CHEBY> ERROR in casting ref lens radius: "//trim(lnSplit(3))
      goto 30
    end if

  else if(chr_toUpper(inLine(1:5)) == "TABLE" ) then
    ! Read chebyshev coefficients in tabular format
    if (chr_toLower(lnSplit(3)) == "y" ) then
      ! line with header
      if ( lnSplit(2) /= ":" .or. lnSplit(4) /= "|" ) then
        write(lerr,"(a)") "CHEBY> ERROR Incorret format of header line."
        goto 50
      end if
      ! - in case, de-allocate
      if ( allocated( jOrder ) ) call dealloc(jOrder,'jOrder')
      ! - allocate
      call alloc( jOrder, nSplit-4, 0, 'jOrder' )
      ! - parse header and store it in memory
      do jj=1,nSplit-4
        call chr_cast(lnSplit(jj+4),jOrder(jj),spErr)
        if ( jOrder(jj) < 0 ) then
          write(lerr,"(a)") "CHEBY> ERROR Header line contains a negative index."
          goto 50
        end if
      end do
      ! - find max order in header
      jMax=maxval( jOrder )
    elseif (chr_toLower(lnSplit(3)) == "x" ) then
      ! line with actual coefficients
      if ( .not. allocated( jOrder ) ) then
        write(lerr,"(a)") "CHEBY> ERROR Coefficient line comes before/without header line."
        goto 50
      end if
      if ( lnSplit(2) /= ":" .or. lnSplit(5) /= "|" ) then
        write(lerr,"(a)") "CHEBY> ERROR Incorret format of line of coefficients."
        goto 50
      end if
      if ( nSplit-5 > size(jOrder) ) then
        write(lerr,"(a)") "CHEBY> ERROR Incorret format of line of coefficients."
        goto 50
      end if
      call chr_cast(lnSplit(4),ii,spErr)
      ! add new coefficients
      iDim = size(cheby_coeffs,1)
      jDim = size(cheby_coeffs,2)
      if(ii>=iDim)   call alloc(cheby_coeffs,     ii, jDim-1, ncheby_tables, zero, 'cheby_coeffs', 0, 0, 1 )
      iDim = size(cheby_coeffs,1)
      if(jMax>=jDim) call alloc(cheby_coeffs, iDim-1,   jMax, ncheby_tables, zero, 'cheby_coeffs', 0, 0 ,1 )
      do jj=1,nSplit-5
        call chr_cast(lnSplit(jj+5),cheby_coeffs(ii,jOrder(jj),ifile),spErr)
      end do
      cheby_maxOrder(ifile)=max(ii,jMax,cheby_maxOrder(ifile))
      if ( cheby_maxOrder(ifile) > cheby_maxNterms-1 ) then
        write(lerr,"(2(a,i0))") "CHEBY> ERROR order of map too high: ",cheby_maxOrder(ifile), &
             " > ",cheby_maxNterms-1
        write(lerr,"(a)")       "CHEBY>       please increase cheby_maxNterms parameter"
        goto 30
      end if
    else
      write(lerr,"(a)") "CHEBY> ERROR Un-recognized line for tabular format of coefficients."
      goto 50
    end if
    
  else
    ! Read chebyshev coefficients in format: "ii jj : value"
    if(nSplit /= 4) then
      write(lerr,"(a)") "CHEBY> ERROR Not enough arguments for expressing Chebyshev coefficients [Vm]."
      write(lerr,"(a)") "CHEBY>       Correct format:"
      write(lerr,"(a)") "CHEBY>       ii jj : value (ii->x,jj->y)"
      goto 30
    end if
    call chr_cast(lnSplit(1),ii,spErr)
    if(spErr) then
      write(lerr,"(a)") "CHEBY> ERROR in casting first index of cheby coeff: "//trim(lnSplit(1))
      goto 30
    end if
    call chr_cast(lnSplit(2),jj,spErr)
    if(spErr) then
      write(lerr,"(a)") "CHEBY> ERROR in casting second index of cheby coeff: "//trim(lnSplit(2))
      goto 30
    end if
    iDim = size(cheby_coeffs,1)
    jDim = size(cheby_coeffs,2)
    if(ii>=iDim) call alloc(cheby_coeffs,     ii, jDim-1, ncheby_tables, zero, 'cheby_coeffs', 0, 0, 1 )
    iDim = size(cheby_coeffs,1)
    if(jj>=jDim) call alloc(cheby_coeffs, iDim-1,     jj, ncheby_tables, zero, 'cheby_coeffs', 0, 0 ,1 )
    call chr_cast(lnSplit(4),cheby_coeffs(ii,jj,ifile),spErr)
    if(spErr) then
      write(lerr,"(a)") "CHEBY> ERROR in casting Chebyshev coefficient: "//trim(lnSplit(4))
      goto 30
    end if
    cheby_maxOrder(ifile)=max(ii,jj,cheby_maxOrder(ifile))
    if ( cheby_maxOrder(ifile) > cheby_maxNterms-1 ) then
      write(lerr,"(2(a,i0))") "CHEBY> ERROR order of map too high: ",cheby_maxOrder(ifile), &
                                            " > ",cheby_maxNterms-1
      write(lerr,"(a)")       "CHEBY>       please increase cheby_maxNterms parameter"
      goto 30
    end if

  end if ! close if for keyword identification
  goto 10

20 continue

  call f_freeUnit(fUnit)
  if (cheby_refR(ifile)<=zero) then
    write(lerr,"(a)") "CHEBY> ERROR ref lens radius [mm] must be positive."
    goto 30
  end if
  if (cheby_maxOrder(ifile)<2) then
    write(lerr,"(a,i0,a)") "CHEBY> ERROR max order too low:",cheby_maxOrder(ifile)," - it should be at least 2."
    goto 30
  end if
  if ( allocated( jOrder ) ) call dealloc(jOrder,'jOrder')
  write(lout,"('CHEBY> ',a)") "...done."

  if(st_quiet < 2) then
    ! Echo parsed data (unless told to be quiet!)
    if (st_debug) then
      ! detailed echo
      call cheby_echoCoefficients_detailed(ifile,lDefI)
    else 
      ! tabbed echo (more human-readable)
      call cheby_echoCoefficients(ifile,lDefI)
    end if
  end if
  return

50 continue
  write(lerr,"(a)") "CHEBY>       Example of header line:"
  write(lerr,"(a)") "CHEBY>       TABLE : Y   |     0     1     2    etc..."
  write(lerr,"(a)") "CHEBY>       Example of coefficient line:"
  write(lerr,"(a)") "CHEBY>       TABLE : x 1 |   0.0  -0.1   2.2    etc..."
30 continue
  write(lout,"(a)") "CHEBY>       last line read:"
  write(lout,"(a)") trim(inLine)
  write(lout,"(a)") "CHEBY>       split fields:"
  do ii=1,nSplit
    write(lout,"('CHEBY>       - ',i2,': ',a)") ii,trim(lnSplit(ii))
  end do
40 continue
  write(lerr,"(a)") "CHEBY> ERROR while parsing file "//trim(cheby_filename(ifile))
  call prror

end subroutine parseChebyFile


subroutine cheby_echoCoefficients(ifile,lDefI)
  use crcoall, only : lout
  use mod_alloc, only : alloc, dealloc
  implicit none
  integer, intent(in) :: ifile
  logical, intent(in) :: lDefI
  integer, allocatable :: jOrder(:)
  integer ii, jj, kk, maxLen
  character(len=:), allocatable :: ruler, header, tmpLine
  character(len=10) :: tmpString

  ! write header
  write(lout,"(a,i0)") "CHEBY> Coefficients for Chebyshev polynomials as from file '"//&
    trim(cheby_filename(ifile))//"' - #",ifile
  write(lout,"(a,1pe9.2)")   "CHEBY> Reference current [A] : ",cheby_refI(ifile)
  if (lDefI) write(lout,"(a)") "       --> default value!"
  write(lout,"(a,1pe9.2)")   "CHEBY> reference radius [mm] : ",cheby_refR(ifile)

  ! initialise writing strings
  maxLen=7+10*cheby_maxOrder(ifile)
  allocate(character(len=maxLen) :: ruler)
  allocate(character(len=maxLen) :: header)
  allocate(character(len=maxLen) :: tmpLine)
  ruler=" "
  header=" "
  tmpLine=" "

  ! initialise header and ruler
  ruler =repeat("-",7)
  header=repeat(" ",5)//"|"

  ! get non-NULL columns
  do jj=0,cheby_maxOrder(ifile)
    do ii=0,cheby_maxOrder(ifile)
      if ( cheby_coeffs(ii,jj,ifile) /= zero ) then
        if ( allocated( jOrder ) ) then
          call alloc( jOrder, size(jOrder)+1, jj, 'jOrder' )
        else
          call alloc( jOrder,              1, jj, 'jOrder' )
        end if
        ruler=trim(ruler)//repeat("-",10)
        write(tmpString,'(1X,i9)') jj
        header=trim(header)//tmpString
        exit
      end if
    end do
  end do

  ! write table
  write(lout,"('CHEBY> ',a)") trim(ruler)
  write(lout,"('CHEBY> ',a)") trim(header)
  write(lout,"('CHEBY> ',a)") trim(ruler)
  do ii=0,cheby_maxOrder(ifile)
    do kk=1,size(jOrder)
      if ( cheby_coeffs(ii,jOrder(kk),ifile) /= zero ) then
        ! the line should be written
        write(tmpString,"(i4,1X,'|')") ii
        tmpLine=tmpString
        do jj=1,size(jOrder)
          if ( cheby_coeffs(ii,jOrder(jj),ifile) == zero ) then
            tmpLine=trim(tmpLine)//repeat(" ",10)
          else   
            write(tmpString,'(1X,1pe9.2)') cheby_coeffs(ii,jOrder(jj),ifile)
            tmpLine=trim(tmpLine)//tmpString
          end if
        end do
        write(lout,"('CHEBY> ',a)") trim(tmpLine)
        exit
      end if
    end do
  end do
  write(lout,"('CHEBY> ',a)") trim(ruler)
  write(lout,"('CHEBY> ',a)") ""
  
  call dealloc(jOrder,'jOrder')
  deallocate(ruler)
  deallocate(header)
  deallocate(tmpLine)
  return
end subroutine cheby_echoCoefficients


subroutine cheby_echoCoefficients_detailed(ifile,lDefI)
  use crcoall, only : lout
  implicit none
  integer, intent(in) :: ifile
  logical, intent(in) :: lDefI
  integer ii, jj
  
  write(lout,"(a,i0)") "CHEBY> Coefficients for Chebyshev polynomials as from file '"//&
    trim(cheby_filename(ifile))//"' - #",ifile
  write(lout,"(a,1pe22.15)")   "CHEBY> Reference current [A] : ",cheby_refI(ifile)
  if (lDefI) write(lout,"(a)") "       --> default value!"
  write(lout,"(a,1pe22.15)")   "CHEBY> reference radius [mm] : ",cheby_refR(ifile)
  do ii=0,cheby_maxOrder(ifile)
    do jj=0,cheby_maxOrder(ifile)
      if(cheby_coeffs(ii,jj,ifile)/= zero) then
        write(lout,"(2(a,i4),a,1pe22.15)") "CHEBY> Order ",ii,",",jj," : ",cheby_coeffs(ii,jj,ifile)
      end if
    end do
  end do
  write(lout,"(a)" ) ""
  return
end subroutine cheby_echoCoefficients_detailed


subroutine cheby_kick(i,ix,n)

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 02-04-2020
  ! apply kick of Chebyshev lenses

  use mod_common, only : beta0, napx, brho
  use mod_common_main
  use mathlib_bouncer
  use numerical_constants, only : zero, c180e0, pi, two, c1m15
  use physical_constants, only: clight

  implicit none

  integer, intent(in) :: i
  integer, intent(in) :: ix
  integer, intent(in) :: n

  real(kind=fPrec) xx, yy, ax, ay, rr, dxp, dyp
  real(kind=fPrec) theta, angle_rad, epsilon, gteps, lteps, r1lteps, r2gteps
  integer          jj
  logical          lrotate

  epsilon=c1m15
  gteps=one+epsilon
  lteps=one-epsilon
  r1lteps=cheby_r1(icheby(ix))*lteps
  r2gteps=cheby_r2(icheby(ix))*gteps
  
  ! rotation angle
  lrotate = cheby_angle(icheby(ix)).ne.zero
  angle_rad = (cheby_angle(icheby(ix))/c180e0)*pi

  do jj=1,napx

    ! apply offset
    xx=xv1(jj)-cheby_offset_x(icheby(ix))
    yy=xv2(jj)-cheby_offset_y(icheby(ix))
    if (abs(xx)<epsilon) xx=zero
    if (abs(yy)<epsilon) yy=zero

    ! in case of non-zero tilt angle, rotate coordinates
    if (lrotate) then
      ! rr = sqrt(xx**2+yy**2)
      rr=sqrt((xx+yy)*(xx+yy)-two*(xx*yy))
      if (abs(rr)<epsilon) then
        xx=zero
        yy=zero
      else
        theta = atan2_mb(yy, xx)-angle_rad
        xx = rr * cos_mb(theta)
        yy = rr * sin_mb(theta)
      end if
    end if
   
    ! check that particle is within the domain of chebyshev polynomials
    ax=abs(xx)
    ay=abs(yy)
    ! (x,y)<r1 or ( (xx>r2) || (yy>r2) ): no kick from lens
    if ( (ax < r1lteps .and. ay < r1lteps) .or. (ax > r2gteps .or. ay > r2gteps) ) cycle

    ! compute kick from cheby map
    call cheby_getKick( xx, yy, dxp, dyp, cheby_itable(icheby(ix)) )
    if (abs(dxp)<epsilon) dxp=zero
    if (abs(dyp)<epsilon) dyp=zero
    if ( dxp == zero .and. dyp == zero ) cycle
    
    ! in case cheby has a non-zero angle, rotate kicks
    if (lrotate) then
      ! NB: cheby_angle(icheby(ix)) is the rotation angle of the cheby
      theta = atan2_mb(dyp, dxp)+angle_rad
      rr = sqrt((dxp+dyp)*(dxp+dyp)-two*(dxp*dyp))
      dxp = rr * cos_mb(theta)
      dyp = rr * sin_mb(theta)
    end if

    ! take into account scaling factor, Brho of beam and its relativistic beta,
    !    and magnetic rigidity and relativistic beta of particle being tracked
    dxp=(((dxp*cheby_scalingFact(icheby(ix)))/(brho*(clight*beta0)))*moidpsv(jj))*rvv(jj)
    dyp=(((dyp*cheby_scalingFact(icheby(ix)))/(brho*(clight*beta0)))*moidpsv(jj))*rvv(jj)

    ! apply kicks
    yv1(jj)=yv1(jj)+dxp
    yv2(jj)=yv2(jj)+dyp

  end do

end subroutine cheby_kick


subroutine cheby_kick_fox(i,ix)

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 06-04-2020
  ! apply kick of Chebyshev lenses (FOX)
  ! NB: for derivatives of polynomials, do not follow what is suggested by
  !     G.Stancari, but get d(Tn(u))/du from its definition

  use mod_common, only : beta0, mtcda, brho
  use mod_settings, only : st_debug
  use crcoall, only : lout
  use numerical_constants, only : zero, one, c180e0, pi, c1m15, two, pieni
  use physical_constants, only: clight
  use mod_lie_dab, only : lnv, idao, rscrri, iscrda
  use mod_common_da
  
  implicit none

  integer, intent(in) :: i
  integer, intent(in) :: ix

  integer          :: jj, morder, nn, mm, ll, oo
  real(kind=fPrec) :: rra, tta, xa, ya, ax, ay, xclo, yclo, angrad, cscal, dxpa, dypa, refr, coeff
  real(kind=fPrec) :: epsilon, gteps, lteps, r1lteps, r2gteps
  logical          lrotate

! for FOX
  integer          :: idaa
  integer          :: hh(lnv)=0
  common/daele/alda,asda,aldaq,asdaq,smida,xx,yy,dpda,dpda1,sigmda,ej1,ejf1,rv
  save

!-----------------------------------------------------------------------
!FOX  B D ;
#include "include/dainicom.f90"
!FOX  D V DA INT XI    NORD NVAR ;
!FOX  D V DA INT YI    NORD NVAR ;
!FOX  D V DA INT DXP   NORD NVAR ;
!FOX  D V DA INT DYP   NORD NVAR ;
!FOX  D V DA INT RR    NORD NVAR ;
!FOX  D V DA INT RR_SQ NORD NVAR ;  
!FOX  D V DA INT RADIO NORD NVAR ;
!FOX  D V DA INT UU    NORD NVAR ;
!FOX  D V DA INT VV    NORD NVAR ;
!FOX  D V DA INT THETA NORD NVAR ;
!FOX  D V DA INT TX    NORD NVAR 21 ;
!FOX  D V DA INT TY    NORD NVAR 21 ;
!FOX  D V DA INT TPX   NORD NVAR 21 ;
!FOX  D V DA INT TPY   NORD NVAR 21 ;
!FOX  D V RE INT XCLO ;
!FOX  D V RE INT YCLO ;
!FOX  D V RE INT CSCAL ;
!FOX  D V RE INT BETA0 ;
!FOX  D V RE INT BRHO ;
!FOX  D V RE INT CLIGHT ;
!FOX  D V RE INT PI ;
!FOX  D V RE INT XA ;
!FOX  D V RE INT YA ;
!FOX  D V RE INT REFR ;
!FOX  D V RE INT DXPA ;
!FOX  D V RE INT DYPA ;
!FOX  D V RE INT COEFF ;
!FOX  D V RE INT ANGRAD ;
!FOX  D V RE INT ONE ;
!FOX  D V RE INT ZERO ;
!FOX  D V RE INT C1E3 ;
!FOX  D V RE INT C1M3 ;
!FOX  D V RE INT TWO ;
!FOX  D V RE INT NN ;
!FOX  D V RE INT JJ ;
!FOX  D V RE INT MM ;
!FOX  D V RE INT LL ;
!FOX  D V RE INT OO ;
!FOX  E D ;
!FOX  1 if(1.eq.1) then
!-----------------------------------------------------------------------

  epsilon=c1m15
  gteps=one+epsilon
  lteps=one-epsilon
  r1lteps=cheby_r1(icheby(ix))*lteps
  r2gteps=cheby_r2(icheby(ix))*gteps
  
  XCLO=cheby_offset_x(icheby(ix))
  YCLO=cheby_offset_y(icheby(ix))
  CSCAL=cheby_scalingFact(icheby(ix))
  REFR=cheby_refR(cheby_itable(icheby(ix)))
  MORDER=cheby_maxOrder(cheby_itable(icheby(ix)))
  
  RRA=zero
  XA=zero
  YA=zero

  ! rotation angle
  lrotate = cheby_angle(icheby(ix)).ne.zero
  angrad = (cheby_angle(icheby(ix))/c180e0)*pi

  if (st_debug) then
    write(lout,'(3(a,i0))')'CHEBY> CHEBY_KICK_FOX for i= ',i,' - ix= ',ix, ' - icheby(ix)= ',icheby(ix)
    write(lout,'(a)')      'CHEBY> CHEBY_KICK_FOX closed orbit BEFORE cheby map:'
    call dapri(XX(1),6)
    call dapri(XX(2),6)
    call dapri(YY(1),6)
    call dapri(YY(2),6)
  end if
 
  ! apply offset
!FOX  XI=XX(1)-XCLO ;
!FOX  YI=XX(2)-YCLO ;
  call dapek(XI,hh,XA)
  if (abs(XA)<epsilon) then
!FOX  XI=ZERO ;
  end if
  call dapek(YI,hh,YA)
  if (abs(YA)<epsilon) then
!FOX  YI=ZERO ;
  end if
    
  ! in case of non-zero tilt angle, rotate coordinates
  if (lrotate) then
    if (st_debug) then
      write(lout,'(1(a,1pe23.16))')'CHEBY> CHEBY_KICK_FOX taking into account angle [deg]= ',cheby_angle(icheby(ix))
    end if
    ! radial position of main beam relative to center of cheby beam
    ! rr = sqrt(xx**2+yy**2)
!FOX  RR_SQ=(XI+YI)*(XI+YI)-(TWO*XI)*YI ;
    call dapek(RR_SQ,hh,RRA)
    if ( abs(RRA)>epsilon**2 ) then
!FOX  RR=SQRT(RR_SQ) ;
      ! no ATAN2 in FOX!
      call dapek(XI,hh,XA)
      call dapek(YI,hh,YA)
      if ( XA == zero ) then
        if ( YA > zero ) then
!FOX      THETA= PI/TWO ;
        else
!FOX      THETA=-PI/TWO ;
        end if
      else
!FOX    THETA=ATAN(YI/XI) ;
        if ( XA < zero ) then
          if ( YA >= zero ) then
!FOX        THETA=THETA+PI ;
          else
!FOX        THETA=THETA-PI ;
          end if
        end if
      end if      
!FOX  THETA=THETA-ANGRAD ;
!FOX  XI=RR*COS(THETA) ;
!FOX  YI=RR*SIN(THETA) ;
    else
!FOX  RR=ZERO ;
!FOX  XI=ZERO ;
!FOX  YI=ZERO ;
    end if
  end if
  if (st_debug) then
    call dapek(XI,hh,XA)
    call dapek(YI,hh,YA)
    call dapek(RR,hh,RRA)
    call dapek(THETA,hh,TTA)
    write(lout,'(2(a,1pe23.16))')'CHEBY> CHEBY_KICK_FOX computing at XA[mm]= ',XA,' - YA[mm]= ',YA
    write(lout,'(2(a,1pe23.16))')'CHEBY>                corresponding to RRA= ',RRA,' - THETA[deg]= ',TTA/pi*c180e0
  end if
    
  ! check that particle is within the domain of chebyshev polynomials
  call dapek(XI,hh,XA)
  call dapek(YI,hh,YA)
  ax=abs(XA) ;
  ay=abs(YA) ;
  ! (x,y)<r1 or ( (xx>r2) || (yy>r2) ): no kick from lens
  if (.not. ( (ax < r1lteps .and. ay < r1lteps) .or. (ax > r2gteps .or. ay > r2gteps) ) ) then
     
    if (st_debug) write(lout,'(2(a,1pe23.16))')'CHEBY> CHEBY_KICK_FOX computing kick when R1[mm]= ',&
         cheby_r1(icheby(ix)),' and R2[mm]= ',cheby_r2(icheby(ix))
    
    ! compute kick
    
    ! normalised variables
!FOX  UU=XI/REFR ;
!FOX  VV=YI/REFR ;
    
    ! fox: T/TP arrays go from 1 to NN
    ! - polynomials:
!FOX  TX(1)=ONE ;
!FOX  TY(1)=ONE ;
!FOX  TX(2)=UU ;
!FOX  TY(2)=VV ;
    ! - derivatives:
!FOX  TPX(1)=ZERO ;
!FOX  TPY(1)=ZERO ;
!FOX  TPX(2)=ONE ;
!FOX  TPY(2)=ONE ;
    do NN=3,morder+1
      MM=NN-1
      OO=NN-2
!FOX  TX(NN)=(TWO*(UU*TX(MM)))-TX(OO) ;
!FOX  TY(NN)=(TWO*(VV*TY(MM)))-TY(OO) ;
      ! use definition of derivative applied to definition of polynomial:
      ! d/du Tn(u)= d/du[ 2u*T_(n-1)(u)-T_(n-2)(u) ] = 2T_(n-1)(u) +2u*T'_(n-1)(u) -T'_(n-2)(u)
      !           = 2(T_(n-1)(u) +u*T'_(n-1)(u))-T'_(n-2)(u)
!FOX  TPX(NN)=TWO*(TX(MM)+UU*TPX(MM))-TPX(OO) ;
!FOX  TPY(NN)=TWO*(TY(MM)+VV*TPY(MM))-TPY(OO) ;
    end do
    
    ! get kicks
!FOX  DXP=ZERO ;
!FOX  DYP=ZERO ;
    do NN=0,morder
      do JJ=0,NN
        COEFF=cheby_coeffs(JJ,NN-JJ,cheby_itable(icheby(ix)))
        if ( COEFF /= zero ) then
          MM=JJ+1
          LL=NN-JJ+1
!FOX      DXP=DXP+(COEFF*TPX(MM))*TY (LL) ;
!FOX      DYP=DYP+(COEFF*TX (MM))*TPY(LL) ;
        end if
      end do
    end do
    ! ref radius in [mm], kick in [mrad]   
!FOX  DXP=-(DXP*C1E3)/(REFR*C1M3) ;
!FOX  DYP=-(DYP*C1E3)/(REFR*C1M3) ;
  
    call dapek(DXP,hh,DXPA)
    if (abs(DXPA)<epsilon) then
!FOX  DXP=ZERO ;
    end if       
    call dapek(DYP,hh,DYPA)
    if (abs(DYPA)<epsilon) then
!FOX  DYP=ZERO ;
    end if
    
    call dapek(DXP,hh,DXPA)
    call dapek(DYP,hh,DYPA)
    if ( DXPA /= zero .and. DYPA /= zero ) then
      if (lrotate) then
!FOX    RR=SQRT(DXP*DXP+DYP*DYP) ;
        ! no ATAN2 in FOX!
        if ( DXPA == zero ) then
          if ( DYPA > zero ) then
!FOX        THETA= PI/TWO ;
          else
!FOX        THETA=-PI/TWO ;
          end if
        else
!FOX      THETA=ATAN(DYP/DXP) ;
          if ( DXPA < zero ) then
            if ( DYPA >= zero ) then
!FOX          THETA=THETA+PI ;
            else
!FOX          THETA=THETA-PI ;
            end if
          end if
        end if      
!FOX    THETA=THETA+ANGRAD ;
!FOX    DXP=RR*COS(THETA) ;
!FOX    DYP=RR*SIN(THETA) ;
      end if
  
      ! take into account scaling factor, Brho of beam and its relativistic beta,
      !    and magnetic rigidity and relativistic beta of particle being tracked
!FOX  DXP=(((DXP*CSCAL)/(BRHO*(CLIGHT*BETA0)))*(MTCDA/(ONE+DPDA)))*RV ;
!FOX  DYP=(((DYP*CSCAL)/(BRHO*(CLIGHT*BETA0)))*(MTCDA/(ONE+DPDA)))*RV ;

      ! apply kicks
!FOX  YY(1)=YY(1)+DXP ;
!FOX  YY(2)=YY(2)+DYP ;
    end if
  end if

  if (st_debug) then
    call dapek(DXP,hh,DXPA)
    call dapek(DYP,hh,DYPA)
    write(lout,'(2(a,1pe23.16))')'CHEBY> CHEBY_KICK_FOX - computed kick: DXPA[mrad]= ', &
      DXPA,' - DYPA[mrad]= ',DYPA
    write(lout,'(a)')'CHEBY> CHEBY_KICK_FOX closed orbit AFTER cheby lens:'
    call dapri(XX(1),6)
    call dapri(XX(2),6)
    call dapri(YY(1),6)
    call dapri(YY(2),6)
  end if

! Do not remove or modify the comment below.
!     DADAL AUTOMATIC INCLUSION

end subroutine cheby_kick_fox


subroutine cheby_potentialMap(iLens,ix)

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 02-04-2020
  ! dump map of potential

  use crcoall, only : lout, lerr
  use mod_common, only : bez
  use mod_common_main
  use mathlib_bouncer
  use numerical_constants, only : zero, c180e0, pi, c1m15
  use mod_units

  integer, intent(in) :: iLens
  integer, intent(in) :: ix

  real(kind=fPrec) xx, yy, rr, zz, dx, dy, xxr, yyr, xxn, yyn, ax, ay
  real(kind=fPrec) theta, radio, angle_rad, epsilon, gteps, lteps, r1lteps, r2gteps
  integer          ii, jj, inside, fUnit
  logical          lrotate, err

  epsilon=c1m15
  gteps=one+epsilon
  lteps=one-epsilon
  r1lteps=cheby_r1(icheby(ix))*lteps
  r2gteps=cheby_r2(icheby(ix))*gteps
  
  write(lout,"(a)") "CHEBY> Dumping potential map..."
  call f_requestUnit(cheby_mapFileName(iLens),fUnit)
  call f_open(unit=fUnit,file=cheby_mapFileName(iLens),mode='w',err=err,formatted=.true.,status="replace")
  if(err) then
    write(lerr,"(a)") "CHEBY> ERROR Failed to open file."
    call prror
  end if

  ! rotation angle
  lrotate = cheby_angle(iLens).ne.zero
  angle_rad = (cheby_angle(iLens)/c180e0)*pi
  dx=(cheby_mapXmax(iLens)-cheby_mapXmin(iLens))/real(cheby_mapNx(iLens),fPrec)
  dy=(cheby_mapYmax(iLens)-cheby_mapYmin(iLens))/real(cheby_mapNy(iLens),fPrec)

  ! write header
  write(fUnit,'("# ",a,i0)')       "name of chebyshev lens: '"//trim(bez(ix))//"' - id: ",iLens
  write(fUnit,"('# ',a)")          "filename              : '"//trim(cheby_filename(cheby_itable(iLens)))//"'"
  write(fUnit,"('# ',a,1pe22.15)") "R2                [mm]: ",cheby_r2(iLens)
  write(fUnit,"('# ',a,1pe22.15)") "R1                [mm]: ",cheby_r1(iLens)
  write(fUnit,"('# ',a,1pe22.15)") "tilt angle       [deg]: ",cheby_angle(iLens)
  write(fUnit,"('# ',a,1pe22.15)") "hor offset        [mm]: ",cheby_offset_x(iLens)
  write(fUnit,"('# ',a,1pe22.15)") "ver offset        [mm]: ",cheby_offset_y(iLens)
  write(fUnit,"('# ',a,1pe22.15)") "lens powering      [A]: ",cheby_I(iLens)
  write(fUnit,"('# ',a,1pe22.15)") "scaling factor      []: ",cheby_scalingFact(iLens)
  write(fUnit,"('# ',a)")          "map filename          : '"//trim(cheby_mapFileName(iLens))//"'"
  write(fUnit,"('# ',a,1pe22.15)") "xmin              [mm]: ",cheby_mapXmin(iLens)
  write(fUnit,"('# ',a,1pe22.15)") "xmax              [mm]: ",cheby_mapXmax(iLens)
  write(fUnit,"('# ',a,i0)")       "Nx                  []: ",cheby_mapNx  (iLens)
  write(fUnit,"('# ',a,1pe22.15)") "ymin              [mm]: ",cheby_mapYmin(iLens)
  write(fUnit,"('# ',a,1pe22.15)") "ymax              [mm]: ",cheby_mapYmax(iLens)
  write(fUnit,"('# ',a,i0)")       "Ny                  []: ",cheby_mapNy  (iLens)
  write(fUnit,"('# ',a)") "x [mm], y [mm], x_map [mm], y_map [mm], V [V m], inside [0:False,1:True]"

  ! get map
  do jj=0,cheby_mapNy(iLens)
    yy=cheby_mapYmin(iLens)+(real(jj,fPrec)*dy) ! mesh y-coordinate
    if (jj==cheby_mapNy(iLens)) yy=cheby_mapYmax(iLens)
    yyr=yy-cheby_offset_y(iLens)  ! take into account y-offset of Cheby map

    do ii=0,cheby_mapNx(iLens)
      xx=cheby_mapXmin(iLens)+(real(ii,fPrec)*dx) ! mesh x-coordinate
      if (ii==cheby_mapNx(iLens)) xx=cheby_mapXmax(iLens)
      xxr=xx-cheby_offset_x(iLens)  ! take into account x-offset of Cheby map

      ! in case of non-zero tilt angle, rotate coordinates
      if (lrotate) then
        rr=sqrt(xxr**2+yyr**2)
        theta = atan2_mb(yyr, xxr)-angle_rad
        xxn = rr * cos_mb(theta)
        yyn = rr * sin_mb(theta)
      else
        theta=zero
        xxn=xxr
        yyn=yyr
      end if
     
      ! check that particle is within the domain of chebyshev polynomials
      inside=1
      ax=abs(xxn)
      ay=abs(yyn)
      ! (x,y)<r1 or ( (xx>r2) || (yy>r2) ): no kick from lens
      if ( (ax < r1lteps .and. ay < r1lteps) .or. (ax > r2gteps .or. ay > r2gteps) ) inside=0
      ! compute kick from cheby map
      call cheby_getPotential( xxn, yyn, zz, cheby_itable(iLens) )
      write(fUnit,'(5(1X,1pe22.15),1X,i0)') xx, yy, xxn, yyn, zz, inside

    end do
  end do

  call f_freeUnit(fUnit)

end subroutine cheby_potentialMap


subroutine cheby_getKick( xx, yy, dxp, dyp, iTable )

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 06-04-2020
  ! compute kicks from Chebyshev polinomials - see FermiLAB-FN-0972-APC
  ! coordinates and kicks are in the map reference frame!
  ! NB: for derivatives of polynomials, do not follow what is suggested by
  !     G.Stancari, but get d(Tn(u))/du from its definition

  use numerical_constants, only : zero, one, two, c1m3, c1e3

  ! interface vars
  real(kind=fPrec), intent(in ) :: xx
  real(kind=fPrec), intent(in ) :: yy
  real(kind=fPrec), intent(out) :: dxp
  real(kind=fPrec), intent(out) :: dyp
  integer,          intent(in ) :: iTable

  ! temp vars
  real(kind=fPrec) :: uu, vv, Tx (0:cheby_maxOrder(iTable)), Ty (0:cheby_maxOrder(iTable)), &
                              Tpx(0:cheby_maxOrder(iTable)), Tpy(0:cheby_maxOrder(iTable))
  integer          :: nn, jj

  ! normalised variables
  uu=xx/cheby_refR(iTable)
  vv=yy/cheby_refR(iTable)

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
    ! use definition of derivative applied to definition of polynomial:
    ! d/du Tn(u)= d/du[ 2u*T_(n-1)(u)-T_(n-2)(u) ] = 2T_(n-1)(u) +2u*T'_(n-1)(u) -T'_(n-2)(u)
    !           = 2(T_(n-1)(u) +u*T'_(n-1)(u))-T'_(n-2)(u)
    Tpx(nn)=two*(Tx(nn-1)+uu*Tpx(nn-1))-Tpx(nn-2)
    Tpy(nn)=two*(Ty(nn-1)+vv*Tpy(nn-1))-Tpy(nn-2)
  end do

  ! get kicks
  dxp=zero
  dyp=zero
  do nn=0,cheby_maxOrder(iTable)
    do jj=0,nn
      if (cheby_coeffs(jj,nn-jj,iTable) /= zero ) then
        dxp=dxp+(cheby_coeffs(jj,nn-jj,iTable)*Tpx(jj))*Ty (nn-jj)
        dyp=dyp+(cheby_coeffs(jj,nn-jj,iTable)*Tx (jj))*Tpy(nn-jj)
      end if
    end do
  end do
  dxp=-(dxp*c1e3)/(cheby_refR(iTable)*c1m3) ! ref radius in [mm], kick in [mrad]
  dyp=-(dyp*c1e3)/(cheby_refR(iTable)*c1m3) ! ref radius in [mm], kick in [mrad]

 end subroutine cheby_getKick


subroutine cheby_getPotential( xx, yy, zz, iTable )

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 01-03-2019
  ! compute potential from Chebyshev polinomials - see FermiLAB-FN-0972-APC
  ! coordinates and potential are in the map reference frame!

  use numerical_constants, only : zero, one, two

  ! interface vars
  real(kind=fPrec), intent(in ) :: xx
  real(kind=fPrec), intent(in ) :: yy
  real(kind=fPrec), intent(out) :: zz
  integer,          intent(in ) :: iTable

  ! temp vars
  real(kind=fPrec) :: uu, vv, Tx(0:cheby_maxOrder(iTable)), Ty(0:cheby_maxOrder(iTable))
  integer          :: nn, jj

  ! normalised variables
  uu=xx/cheby_refR(iTable)
  vv=yy/cheby_refR(iTable)

  ! - polynomials:
  Tx(0)=one
  Ty(0)=one
  Tx(1)=uu
  Ty(1)=vv
  do nn=2,cheby_maxOrder(iTable)
    Tx(nn)=two*(uu*Tx(nn-1))-Tx(nn-2)
    Ty(nn)=two*(vv*Ty(nn-1))-Ty(nn-2)
  end do

  ! get potential
  zz=zero
  do nn=0,cheby_maxOrder(iTable)
    do jj=0,nn
      if (cheby_coeffs(jj,nn-jj,iTable) /= zero ) zz=zz+(cheby_coeffs(jj,nn-jj,iTable)*Tx(jj))*Ty(nn-jj)
    end do
  end do

 end subroutine cheby_getPotential

subroutine cheby_setScaleKick( iCheby )

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 05-03-2019
  ! compute scaling factor for chebyshev lenses

  ! interface vars
  integer, intent(in) :: iCheby
  cheby_scalingFact(iCheby)=cheby_I(iCheby)/cheby_refI(cheby_itable(iCheby))

end subroutine cheby_setScaleKick

subroutine cheby_resetI( iCheby, Inew )

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 02-04-2020
  ! reset current of chebyshev lenses

  ! interface vars
  integer,                    intent(in) :: iCheby
  real(kind=fPrec), optional, intent(in) :: Inew

  real(kind=fPrec) Iupdate

  Iupdate = cheby_refI(cheby_itable(iCheby))
  if(present(Inew)) Iupdate = Inew
  cheby_I(iCheby) = Iupdate

end subroutine cheby_resetI

end module cheby
