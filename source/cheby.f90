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
  integer, save               :: ncheby=0             ! number of chebyshev lenses actually in memory
  integer, save               :: ncheby_mapEchoes=0   ! number of requested echoes of maps
  integer, save               :: ncheby_tables=0      ! number of chebyshev tables in memory
  integer, parameter          :: cheby_kz=42          ! kz of chebyshev lenses
  integer, parameter          :: cheby_ktrack=67      ! ktrack of chebyshev lenses

  ! variables to save parameters for tracking etc.
  integer,          allocatable, save :: cheby_itable(:)      ! index of chebyshev table
  real(kind=fPrec), allocatable, save :: cheby_r2(:)          ! outer radius R2 [mm] (optional)
  real(kind=fPrec), allocatable, save :: cheby_r1(:)          ! inner radius R1 [mm] (optional)
  real(kind=fPrec), allocatable, save :: cheby_angle(:)       ! rotation angle about the longitudinal axis [deg] (optional)
  real(kind=fPrec), allocatable, save :: cheby_offset_x(:)    ! hor. offset [mm] (optional)
  real(kind=fPrec), allocatable, save :: cheby_offset_y(:)    ! ver. offset [mm] (optional)
  real(kind=fPrec), allocatable, save :: cheby_I(:)           ! actual powering of lens [A] (optional)
  real(kind=fPrec), allocatable, save :: cheby_scalingFact(:) ! scaling factor [] (computed internally)
  
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
  implicit none
  integer stat
  call alloc(icheby,nele,0,'icheby')
end subroutine cheby_allocate_arrays

subroutine cheby_expand_arrays(nele_new)
  implicit none
  integer, intent(in) :: nele_new
  call alloc(icheby,nele_new,0,'icheby')
end subroutine cheby_expand_arrays

subroutine cheby_expand_arrays_lenses(ncheby_new)
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
end subroutine cheby_expand_arrays_lenses

subroutine cheby_expand_arrays_map_echo(ncheby_mapEchoes_new)
  implicit none
  integer, intent(in) :: ncheby_mapEchoes_new
  ! map
  call alloc(cheby_iLens      ,            ncheby_mapEchoes_new,               0, 'cheby_iLens'      )
  call alloc(cheby_mapFileName, mFileName, ncheby_mapEchoes_new,             " ", 'cheby_mapFileName')
  call alloc(cheby_mapXmin    ,            ncheby_mapEchoes_new, -10.0e+00_fPrec, 'cheby_mapXmin'    )
  call alloc(cheby_mapXmax    ,            ncheby_mapEchoes_new,  10.0e+00_fPrec, 'cheby_mapXmax'    )
  call alloc(cheby_mapNx      ,            ncheby_mapEchoes_new,             100, 'cheby_mapNx'      )
  call alloc(cheby_mapYmin    ,            ncheby_mapEchoes_new, -10.0e+00_fPrec, 'cheby_mapYmin'    )
  call alloc(cheby_mapYmax    ,            ncheby_mapEchoes_new,  10.0e+00_fPrec, 'cheby_mapYmax'    )
  call alloc(cheby_mapNy      ,            ncheby_mapEchoes_new,             100, 'cheby_mapNy'      )
end subroutine cheby_expand_arrays_map_echo

subroutine cheby_expand_arrays_tables(ncheby_tables_new)
  implicit none
  integer, intent(in) :: ncheby_tables_new
  call alloc(cheby_filename   , mFileName,       ncheby_tables_new,  " ", 'cheby_filename'                      )
  call alloc(cheby_coeffs     ,            0, 0, ncheby_tables_new, zero, 'cheby_coeffs'  , 0, 0, ncheby_tables )
  call alloc(cheby_maxOrder   ,                  ncheby_tables_new,    0, 'cheby_maxOrder'                      )
  call alloc(cheby_refI       ,                  ncheby_tables_new,  one, 'cheby_refI'                          )
  call alloc(cheby_refR       ,                  ncheby_tables_new, zero, 'cheby_refR'                          )
end subroutine cheby_expand_arrays_tables

! ================================================================================================ !
!  Parse Line for Chebyshev lens
!  Last modified: 2019-02-25
! ================================================================================================ !
subroutine cheby_parseInputLine(inLine, iLine, iErr)

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
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "CHEBY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(lnSplit(1))

  case("SHOW")
    if(nSplit < 3) then
      write(lout,"(a,i0)") "CHEBY> ERROR Expected at least 3 input parameters, got ",nSplit
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
      write(lout,"(a)") "CHEBY> ERROR Element '"//trim(lnSplit(2))//"' not a Chebyshev lens."
      write(lout,"(a)") "CHEBY>       Either you mis-typed the element name or"
      write(lout,"(a)") "CHEBY>       SHOW line comes before the declaration of the lens"
      iErr = .true.
      return
    end if
    do tmpi1=1,ncheby_mapEchoes
      if (cheby_mapFileName(tmpi1)==lnSplit(3)) then
        write(lout,"(a)") "CHEBY> ERROR File '"//trim(lnSplit(3))//"' already in use."
        iErr = .true.
        return
      end if
    end do
    ncheby_mapEchoes = ncheby_mapEchoes+1
    call cheby_expand_arrays_map_echo(ncheby_mapEchoes)
    cheby_mapFileName(ncheby_mapEchoes) = trim(lnSplit(3))
    cheby_iLens(ncheby_mapEchoes) = icheby(iElem)
    if (nSplit > 3 ) then
      if (nSplit<9) then
        write(lout,"(a,i0)") "CHEBY> ERROR Expected 8 input parameters, got ",nSplit
        write(lout,"(a)")    "CHEBY> format of SHOW line:"
        write(lout,"(a)")    "SHOW name filename xmin[mm] xmax[mm] nx ymin[mm] ymax[mm] ny"
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
     
  case default
  
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
  
    if(icheby(iElem) /= 0) then
      write(lout,"(a)") "CHEBY> ERROR The element '"//trim(bez(iElem))//"' was defined twice."
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
      if(tmpch == cheby_filename(tmpi1)) then
        cheby_itable(icheby(iElem)) = tmpi1
        chIdx = tmpi1
        exit
      end if
    end do
    if(chIdx == -1) then
      ! Unsuccessful search
      ncheby_tables = ncheby_tables+1
      call cheby_expand_arrays_tables(ncheby_tables)
      cheby_itable(icheby(iElem)) = ncheby_tables
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
!  use utils
  use mod_common, only : kz,bez
  use mod_settings, only : st_quiet

  integer ii, jj, kk, ncheb
  logical exist
  real(kind=fPrec) tmpFlt
  
  ! Check that all chebyshev lenses in fort.2 have a corresponding declaration in fort.3
  ncheb=0
  do jj=1,nele
    if(kz(jj)==cheby_kz) then
      if (icheby(jj).eq.0) then
        write(lout,"(a,i0,a)") "CHEBY> ERROR single element ",jj," named '"//trim(bez(jj))//"'"
        write(lout,"(a)")      "CHEBY>       does not have a corresponding line in CHEB block in fort.3"
        call prror
      else
        ncheb=ncheb+1
      end if
    end if
  end do
  if ( ncheb.ne.ncheby ) then
    write(lout,"(a,i0)") "CHEBY> ERROR number of chebyshev lenses declared in CHEB block in fort.3 ",ncheby
    write(lout,"(a,i0)") "CHEBY>       is not the same as the total number of chebyshev lenses in lattice ",ncheb
    call prror
  end if

  ! Parse files with coefficients for Chebyshev polynomials:
  do jj=1,ncheby_tables
    inquire(file=cheby_filename(jj), exist=exist)
    if(.not. exist) then
      write(lout,"(a)") "CHEBY> ERROR Problems with file with coefficients for Chebyshev polynominals: ", &
            trim(cheby_filename(jj))
      call prror(-1)
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
      write(lout,"(a)")      "CHEBY> ERROR R2 cannot be larger than domain of Chebyshev polynomials!"
      write(lout,"(a,1pe22.15,a,1pe22.15)") "CHEBY>       R2 [mm]: ",cheby_r2(jj), &
           " - reference radius [mm]:",cheby_refR(cheby_itable(jj))
      goto 10 
    end if
    if (cheby_r1(jj)==zero) then
      write(lout,"(a)")      "CHEBY> ERROR R1 cannot be zero for the time being!"
      goto 10 
    end if
    if (cheby_I (jj)<=zero) then
      cheby_I (jj)=cheby_refI(cheby_itable(jj))
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
          write(lout,"(a)")                     "CHEBY> ERROR X-extremes for map coincide!"
          write(lout,"(a,1pe22.15,a,1pe22.15)") "CHEBY>       xmin [mm]: ",cheby_mapXmin(ii), &
                                                          " - xmax [mm]: ",cheby_mapXmax(ii)
          goto 10 
        end if
        if (cheby_mapNx(ii)<1) then
          write(lout,"(a)")    "CHEBY> ERROR wrong X-stepping for map!"
          write(lout,"(a,i0)") "CHEBY>       must be >0 - got: ",cheby_mapNx(ii)
          goto 10 
        end if
        if (cheby_mapYmax(ii)<cheby_mapYmin(ii)) then
          ! swap
          tmpFlt=cheby_mapYmax(ii)
          cheby_mapYmax(ii)=cheby_mapYmin(ii)
          cheby_mapYmin(ii)=tmpFlt
        end if
        if (cheby_mapYmax(ii)==cheby_mapYmin(ii)) then
          write(lout,"(a)")                     "CHEBY> ERROR Y-extremes for map coincide!"
          write(lout,"(a,1pe22.15,a,1pe22.15)") "CHEBY>       ymin [mm]: ",cheby_mapYmin(ii), &
                                                          " - ymax [mm]: ",cheby_mapYmax(ii)
          goto 10 
        end if
        if (cheby_mapNy(ii)<1) then
          write(lout,"(a)")    "CHEBY> ERROR wrong Y-stepping for map!"
          write(lout,"(a,i0)") "CHEBY>       must be >0 - got: ",cheby_mapNy(ii)
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
   call prror(-1)
   
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
  integer nSplit, ii, jj, fUnit, iDim, jDim
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
    write(lout,"(a)") "CHEBY> ERROR Failed to split Chebyshev input line"
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
    if(spErr) then
      write(lout,"(a)") "CHEBY> ERROR in casting first index of cheby coeff: "//trim(lnSplit(1))
      goto 30
    end if
    call chr_cast(lnSplit(2),jj,spErr)
    if(spErr) then
      write(lout,"(a)") "CHEBY> ERROR in casting second index of cheby coeff: "//trim(lnSplit(2))
      goto 30
    end if
    iDim = size(cheby_coeffs,1)
    jDim = size(cheby_coeffs,2)
    if(ii>=iDim) call alloc(cheby_coeffs,     ii, jDim-1, ifile, zero, 'cheby_coeffs', 0, 0, 1 )
    iDim = size(cheby_coeffs,1)
    if(jj>=jDim) call alloc(cheby_coeffs, iDim-1,     jj, ifile, zero, 'cheby_coeffs', 0, 0 ,1 )
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
    do ii=0,cheby_maxOrder(ifile)
      do jj=0,cheby_maxOrder(ifile)
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
  
  real(kind=fPrec) xx, yy, rr, dxp, dyp
  real(kind=fPrec) theta, radio, angle_rad
  integer          jj
  logical          lrotate

  ! rotation angle
  lrotate = cheby_angle(icheby(ix)).ne.zero
  angle_rad = (cheby_angle(icheby(ix))/c180e0)*pi

  do jj=1,napx

    ! apply offset
    xx=xv1(jj)-cheby_offset_x(icheby(ix))
    yy=xv2(jj)-cheby_offset_y(icheby(ix))

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
      dxp=(((dxp*cheby_scalingFact(icheby(ix)))/(brho*(clight*betrel)))*moidpsv(jj))*rvv(jj)
      dyp=(((dyp*cheby_scalingFact(icheby(ix)))/(brho*(clight*betrel)))*moidpsv(jj))*rvv(jj)
      
      ! apply kicks
      yv1(jj)=yv1(jj)+dxp
      yv2(jj)=yv2(jj)+dyp
    end if
  end do

end subroutine cheby_kick


subroutine cheby_potentialMap(iLens,ix)

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 01-03-2019
  ! dump map of potential

  use mod_common, only : bez
  use mod_common_main
  use mathlib_bouncer
  use numerical_constants, only : zero, c180e0, pi
  use mod_units

  integer, intent(in) :: iLens
  integer, intent(in) :: ix
  
  real(kind=fPrec) xx, yy, rr, zz, dx, dy, xxr, yyr, xxn, yyn
  real(kind=fPrec) theta, radio, angle_rad
  integer          ii, jj, inside, fUnit
  logical          lrotate, err

  write(lout,"(a)") "CHEBY> Dumping potential map..."
  call f_requestUnit(cheby_mapFileName(iLens),fUnit)
  call f_open(unit=fUnit,file=cheby_mapFileName(iLens),mode='w',err=err,formatted=.true.,status="replace")
  if(err) then
    write(lout,"(a)") "CHEBY> ERROR Failed to open file."
    call prror(-1)
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
    yy=cheby_mapYmin(iLens)+(real(jj,fPrec)*dy) ! mesh coordinate
    if (jj==cheby_mapNy(iLens)) yy=cheby_mapYmax(iLens)
    yyr=yy-cheby_offset_y(iLens)  ! point in ref sys of Cheby map
     
    do ii=0,cheby_mapNx(iLens)
      xx=cheby_mapXmin(iLens)+(real(ii,fPrec)*dx) ! mesh coordinate
      if (ii==cheby_mapNx(iLens)) xx=cheby_mapXmax(iLens)
      xxr=xx-cheby_offset_x(iLens)  ! point in ref sys of Cheby map

      ! check that particle is within the domain of chebyshev polynomials
      rr=sqrt(xxr**2+yyr**2)
      inside=0
      if (rr.ge.cheby_r1(iLens).and.rr.lt.cheby_r2(iLens)) inside=1 ! kick only if rr>=r1 && rr<r2
      ! in case of non-zero tilt angle, rotate coordinates
      if (lrotate) then
        theta = atan2_mb(yyr, xxr)-angle_rad
        xxn = rr * cos_mb(theta)
        yyn = rr * sin_mb(theta)
      else
        theta=zero
        xxn=xxr
        yyn=yyr
      end if
      ! compute kick from cheby map
      call cheby_getPotential( xxn, yyn, zz, cheby_itable(iLens) )
      write(fUnit,'(5(1X,1pe22.15),1X,i0)') xx, yy, xxn, yyn, zz, inside
        
    end do
  end do

  call f_close(fUnit)
  
end subroutine cheby_potentialMap


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
                      fu, fv, Tpx(0:cheby_maxOrder(iTable)), Tpy(0:cheby_maxOrder(iTable))
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
    Tx(nn)=(two*(uu*Tx(nn-1)))-Tx(nn-2)
    Ty(nn)=(two*(vv*Ty(nn-1)))-Ty(nn-2)
    Tpx(nn)=(real(nn,fPrec)*(Tx(nn-1)-(uu*Tx(nn))))/fu
    Tpy(nn)=(real(nn,fPrec)*(Ty(nn-1)-(vv*Ty(nn))))/fv
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


subroutine cheby_getPotential( xx, yy, zz, iTable )

  ! A. Mereghetti (CERN, BE-ABP-HSS)
  ! last modified: 01-03-2019
  ! compute potential from Chebyshev polinomials - see FermiLAB-FN-0972-APC

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

  ! get kicks
  zz=zero
  do nn=0,cheby_maxOrder(iTable)
    do jj=0,nn
      zz=zz+(cheby_coeffs(jj,nn-jj,iTable)*Tx(jj))*Ty(nn-jj)
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
  
end module cheby
