module fma

  use string_tools, only : getfields_l_max_string

  implicit none

  integer, parameter     :: fma_max           = 200                ! Max. number of FMAs
  integer, parameter     :: fma_nturn_max     = 10000              ! Max. number of turns used for fft
  integer, private, save :: fma_numfiles      = 0                  ! Number of FMAs
  logical, public,  save :: fma_flag          = .false.            ! FMA input block exists
  logical, private, save :: fma_writeNormDUMP = .true.             ! Writing out the normalized DUMP files
  character(len=:), allocatable, private, save :: fma_fname(:)     ! Name of input file from dump
  character(len=:), allocatable, private, save :: fma_method(:)    ! Method used to find the tunes
  integer,          allocatable, private, save :: fma_first(:)     ! First turn used for FMA
  integer,          allocatable, private, save :: fma_last(:)      ! Last turn used for FMA
  integer,          allocatable, private, save :: fma_norm_flag(:) ! Normalize phase space before FFT

contains

subroutine fma_allocate
  use mod_alloc
  use parpro
  call alloc(fma_fname,  mStrLen, fma_max, str_dSpace, "fma_fname")
  call alloc(fma_method, mStrLen, fma_max, str_dSpace, "fma_method")
  call alloc(fma_first,           fma_max, 0,          "fma_first")
  call alloc(fma_last,            fma_max, 0,          "fma_last")
  call alloc(fma_norm_flag,       fma_max, 1,          "fma_norm_flag")
end subroutine fma_allocate

  subroutine fma_error(ierro,str,subroutine_name)
    !-----------------------------------------------------------------------*
    !  FMA                                                                  *
    !  M.Fitterer & R. De Maria & K.Sjobak, BE-ABP/HSS                      *
    !  last modified: 04-01-2016                                            *
    !  purpose: error messages for fma analysis                             *
    !-----------------------------------------------------------------------*
      use crcoall
    implicit none
    integer,       intent(in)  :: ierro
    character (*), intent (in) :: subroutine_name
    character (*), intent (in) :: str             !error message
    if(ierro.ne.0) then
       write(lout,*) 'ERROR in ',subroutine_name,': ', str, ', iostat=',ierro
       call prror(-1)
    endif
  end subroutine fma_error

subroutine fma_parseInputline(inLine,iErr)

  use crcoall
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr, cErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "FMA> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(lnSplit(1) == "NoNormDUMP") then
    fma_writeNormDUMP = .false.
    return
  endif

  if(fma_numfiles >= fma_max) then
    write(lout,"(a,i0,a)") "FMA> ERROR You can only do ",fma_max," number of FMAs"
    iErr = .true.
    return
  end if

  fma_numfiles=fma_numfiles+1

  if(nSplit == 1 .or. nSplit == 4 .or. nSplit >= 6) then
    write(lout,"(a,i0)") "FMA> ERROR Wrong number of input parameters. Expected 2, 3 or 5, got ",nSplit
    iErr = .true.
    return
  end if

  fma_fname(fma_numfiles)  = trim(lnSplit(1))
  fma_method(fma_numfiles) = trim(lnSplit(2))
  if(nSplit == 2) then
    fma_norm_flag(fma_numfiles) = 1 ! default: normalize phase space
  else if(nSplit == 3) then
    call chr_cast(lnSplit(3),fma_norm_flag(fma_numfiles),cErr)
  else if(nSplit == 5) then
    call chr_cast(lnSplit(3),fma_norm_flag(fma_numfiles),cErr)
    call chr_cast(lnSplit(4),fma_first(fma_numfiles),    cErr)
    call chr_cast(lnSplit(5),fma_last(fma_numfiles),     cErr)
  end if

  ! Input sanity checks
  if(.not.(&
          fma_method(fma_numfiles) == "TUNELASK"  &
     .or. fma_method(fma_numfiles) == "TUNEFFTI"  &
     .or. fma_method(fma_numfiles) == "TUNEFFT"   &
     .or. fma_method(fma_numfiles) == "TUNEAPA"   &
     .or. fma_method(fma_numfiles) == "TUNEFIT"   &
     .or. fma_method(fma_numfiles) == "TUNENEWT"  &
     .or. fma_method(fma_numfiles) == "TUNEABT2"  &
     .or. fma_method(fma_numfiles) == "TUNEABT"   &
     .or. fma_method(fma_numfiles) == "TUNENEWT1" &
#ifdef NAFF
     .or. fma_method(fma_numfiles) == "NAFF"      &
#endif
    )) then
    write(lout,"(a,i0)") "FMA> ERROR The method '"//trim(fma_method(fma_numfiles))//"' is unknown. FMA index = ",fma_numfiles
    write(lout,"(a)")    "FMA>       Please use one of TUNELASK, TUNEFFTI, TUNEFFT, "// &
      "TUNEAPA, TUNEFIT, TUNENEWT, TUNEABT, TUNEABT2, TUNENEWT1"// &
#ifdef NAFF
      ", NAFF"// &
#endif
      "."
    write(lout,"(a)")    "FMA>       Note that it is case-sensitive, so use uppercase only."
    iErr = .true.
    return
  end if

  if(.not.(fma_norm_flag(fma_numfiles).eq.0 .or. fma_norm_flag(fma_numfiles).eq.1)) then
    write(lout,"(2(a,i0))") "FMA> ERROR Expected fma_norm_flag = 1 or 0, Got: ",fma_norm_flag(fma_numfiles),&
      ", FMA index =",fma_numfiles
    iErr = .true.
    return
  end if

  fma_flag = .true.

end subroutine fma_parseInputline

  subroutine fma_postpr
    !-----------------------------------------------------------------------*
    !  FMA                                                                  *
    !  M.Fitterer & R. De Maria & K.Sjobak, BE-ABP/HSS                      *
    !  last modified: 04-01-2016                                            *
    !  purpose: return files used for fma analysis                          *
    !           -> calculate particle amplitudes and tunes using the        *
    !              normalized coordinates for input files                   *
    !              fma_fname(fma_numfiles)                                  *
    !  output format: q1,q2,q3,eps1_min,eps2_min,eps3_min,eps1_max,         *
    !                 eps2_max,eps3_max,eps1_avg, eps2_avg,eps3_avg,        *
    !                 eps1_0,eps2_0,eps3_0,phi1_0,phi2_0,phi3_0             *
    !-----------------------------------------------------------------------*
    use floatPrecision
    use string_tools
    use mathlib_bouncer
    use platofma

    use dump, only : dump_fname, dumpfmt, dumpunit, dumpfirst, dumplast, dumptas, dumpclo, dumptasinv

    !numbers (zero,one,two etc.)
    use numerical_constants, only : zero, one, c1e3

    use crcoall
    use parpro
    use mod_common
    use mod_commont

    implicit none

    integer :: i,j,k,l,m,n                    ! for do loops
    integer :: num_modes                      ! 3 for 6D tracking, 2 for 4D tracking.
    integer :: fma_npart,fma_tfirst,fma_tlast ! local variables to check input files
    logical :: lopen                          ! flag to check if file is already open
    logical :: lexist                         ! flag to check if file fma_fname exists
    logical :: lread                          ! flag for file reading
    character(len=getfields_l_max_string) :: ch,ch1
    character getfields_fields(getfields_n_max_fields)*(getfields_l_max_string) ! Array of fields
    integer   getfields_nfields                                                 ! Number of identified fields
    integer   getfields_lfields(getfields_n_max_fields)                         ! Length of each what:
    logical   getfields_lerr                                                    ! An error flag
    character filefields_fields ( getfields_n_max_fields )*( getfields_l_max_string )
    integer filefields_nfields
    integer filefields_lfields( getfields_n_max_fields )
    logical filefields_lerr
    real(kind=fPrec) round_near

    integer, dimension(:,:),allocatable :: turn ! Current turn no (particle, rel. turn no)
    integer, dimension(:),allocatable :: nturns ! Number of turns to analyze for this particle
    logical hasNormDumped(-1:nele)              ! Have we written a normDump file for this element before?
    integer fma_nturn (fma_max)                 ! Number of turns used for fft for this FMA
    real(kind=fPrec), dimension(:,:,:),allocatable :: xyzv,nxyzv ! phase space variables (x,x',y,y',z,dE/E)
                                                                 ! [mm,mrad,mm,mrad,mm,1.e-3],
                                                                 ! normalized phase space variables [sqrt(m) 1.e-3]
    real(kind=fPrec), dimension(:,:,:),allocatable :: epsnxyzv   ! normalized emittances
    integer :: dump_last_turn                                    ! auxiliary variable for loop over turns
#ifdef NAFF
    interface
       REAL(C_DOUBLE) function tunenaff (x,xp,maxn,plane_idx,norm_flag) BIND(C)
         use, intrinsic :: ISO_C_BINDING
         IMPLICIT NONE
         REAL(C_DOUBLE), INTENT(in), DIMENSION(1) :: x,xp
         INTEGER(C_INT), INTENT(in), VALUE :: maxn, plane_idx, norm_flag
       end function tunenaff
    end interface

    ! Need to pass a single dimension array to NAFF,
    !  since the stride in the xyzv/nxyzv arrays are difficult to pass correctly to C++.
    ! (We can't interpret the struct that Fortran is passing us;
    !  see the naff_interface.cpp for more info                 )
    real(kind=fPrec), dimension(:), allocatable :: naff_xyzv1,naff_xyzv2
#endif
    ! dummy variables for readin + normalisation + loops
    integer :: id,kt,counter,thisturn
    real(kind=fPrec) :: pos
    real(kind=fPrec), dimension(6) :: xyzvdummy,xyzvdummy2,nxyzvdummy !phase space variables x,x',y,y',sig,delta
    real(kind=fPrec), dimension(3) :: q123 !tune q1,q2,q3
    real(kind=fPrec), dimension(3) :: eps123_0,eps123_min,eps123_max,eps123_avg !initial,minimum,maximum,average emittance
    real(kind=fPrec), dimension(3) :: phi123_0  !initial phase

#ifdef BOINC
    character(len=256) filename
#endif

#ifdef FIO
    ! Do not support FIO, it is not supported by any compilers.
    write (lout,*) "FIO not supported in FMA!"
    call prror(-1)
#endif

    allocate(turn(napx,fma_nturn_max),   &
         xyzv(napx,fma_nturn_max,6),     &
         nxyzv(napx,fma_nturn_max,6),    &
         epsnxyzv(napx,fma_nturn_max,3), &
         STAT=i)
    if (i.ne.0) then
       write(lout,*) "Error in fma_postpr: Cannot ALLOCATE arrays 'turn,xyzv,nxyzv,epsnxyzv'"//&
            " of size proportional to napx*fma_nturn_max."
       call prror(-1)
    endif

    allocate(nturns(napx),STAT=i)
    if (i.ne.0) then
       write(lout,*) "Error in fma_postpr: Cannot ALLOCATE array 'nturns' of size proportional to napx."
       call prror(-1)
    endif

#ifdef NAFF
    allocate(naff_xyzv1(fma_nturn_max), naff_xyzv2(fma_nturn_max), STAT=i )
#endif
    if (i.ne.0) then
       write(lout,*) "Error in fma_postpr: Cannot ALLOCATE arrays 'naff_xyzv1,naff_xyzv2' of size fma_nturn_max."
       call prror(-1)
    endif

    ! Initialize the hasNormDumped array
    do i=-1,nele
       hasNormDumped(i)=.false.
    end do

    ! fma_six = data file for storing the results of the FMA analysis
    inquire(unit=2001001,opened=lopen)
    if(lopen) then
       write(lout,*) "ERROR in FMA: Tried to open unit 2001001 for file 'fma_sixtrack', but it was already taken?"
       call prror(-1)
    endif
#ifdef BOINC
    call boincrf("fma_sixtrack",filename)
    open(2001001,file=filename,      status='replace',iostat=ierro,action='write',form='formatted')
#else
    open(2001001,file='fma_sixtrack',status='replace',iostat=ierro,action='write',form='formatted')
#endif
    call fma_error(ierro,'cannot open file fma_sixtrack for writing!','fma_postpr')

    if (idp.eq.0 .or. ition.eq.0) then
       num_modes = 2          !4D tracking
    else
       num_modes = 3          !6D tracking
    endif

    ! write the header of fma_sixtrack
    write(2001001,'(a)') '# eps1*,eps2*,eps3* all in 1.e-6*m, phi* [rad]'
    write(2001001,'(a)') '# inputfile method id q1 q2 q3 eps1_min '       //&
         'eps2_min eps3_min eps1_max eps2_max eps3_max eps1_avg eps2_avg' //&
         ' eps3_avg eps1_0 eps2_0 eps3_0 phi1_0 phi2_0 phi3_0 norm_flag'  //&
         ' first_turn last_turn'

    ! start FMA analysis: loop over all files, calculate tunes, write output file
    do i=1,fma_numfiles
       lexist=.false.
       do j=-1,nele !START: loop over dump files = loop over single elements
          if(trim(stringzerotrim(fma_fname(i))).eq. &
             trim(stringzerotrim(dump_fname(j)))) then
             lexist=.true.     !set lexist = true if the file fma_fname(j) exists
             write(lout,*) 'start FMA analysis using file ', trim(stringzerotrim(fma_fname(i))), &
                  ': number of particles=',napx, &
                  ', first turn=',fma_first(i),', last turn=',fma_last(i)

             ! check the format, if dumpfmt != 2,3 (physical) or 7,8 (normalized) then abort
             if(.not. (dumpfmt(j).eq.2 .or. dumpfmt(j).eq.3 .or. &
                       dumpfmt(j).eq.7 .or. dumpfmt(j).eq.8)) then
                call fma_error(-1,'input file has wrong format! Choose format=2,3,7 or 8 in DUMP block.','fma_postpr')
             endif
             ! open dump file for reading, resume to original position before exiting the subroutine
             inquire(unit=dumpunit(j),opened=lopen)
             if(lopen) then
                close(dumpunit(j))
             else ! file has to be open if nothing went wrong
                call fma_error(-1,'Expected file '//trim(stringzerotrim(dump_fname(j)))//' to be open','fma_postpr')
             endif

             if (dumpfmt(j).eq.2 .or. dumpfmt(j).eq.7) then
#ifdef BOINC
                call boincrf(dump_fname(j),filename)
                open(dumpunit(j),file=filename,status='old',iostat=ierro,action='read')
#else
                open(dumpunit(j),file=dump_fname(j),status='old',iostat=ierro,action='read')
#endif
                !Create error message to be used in case ierro.ne.0
                write(ch,'(a,1x,I5,1x,a)') "cannot open file '"//trim(stringzerotrim(dump_fname(j)))//"' (dumpfmt=",dumpfmt(j),')'

                call fma_error(ierro, ch, 'fma_postpr')
             else if (dumpfmt(j).eq.3 .or. dumpfmt(j).eq.8) then
#ifdef BOINC
                call boincrf(dump_fname(j),filename)
                open(dumpunit(j),file=filename,status='old',iostat=ierro,action='read',form='unformatted')
#else
                open(dumpunit(j),file=dump_fname(j),status='old',iostat=ierro,action='read',form='unformatted')
#endif
             else
                write(lout,*) 'Error in fma_postpr, got dumpfmt=',dumpfmt(j),'expected 2,3,7 or 8.'
                call prror(-1)
             endif
             ! define first/last turn for FMA
             ! if first and last turn are not defined in FMA block,
             ! take all turns saved in DUMP file
             if (fma_first(i) .eq. 0 .and. fma_last(i) .eq. 0) then
                fma_first(i) = dumpfirst(j)
                fma_last(i)  = dumplast(j)
             endif
             ! if -1 take the last turn of the dump file
             ! or the maximum number of turns if dumplast = -1
             if (fma_last(i) .eq. -1) then
                if (dumplast(j) .eq. -1) then
                   fma_last(i) = numl
                else
                   fma_last(i) = dumplast(j)
                endif
             endif
             ! now check that first turn are compatible with
             ! turns saved in dump file
             if (fma_first(i) .lt. dumpfirst(j)) then
                write(lout,*) 'ERROR in fma_postpr: First turn in FMA block is smaller than first turn in DUMP block '//&
                     'fma_first=',fma_first(i),'< dumpfirst=',dumpfirst(j),'fma_post_pr! This cannot work!'
                call prror(-1)
             endif
             ! now check last turn
             ! if fma_last = -1, we already have fma_last = numl
             ! check if fma_last < 0 and !=-1
             if (fma_last(i) .le. 0) then
                write(lout,*) &
                     "ERROR in fma_postpr: Last turn in FMA block must be -1 or a positive integer, "// &
                     'but fma_last=',fma_last(i),'!'
                call prror(-1)
             endif
             ! if fma_last >0 check that fma_last < dump_last
             if (dumplast(j) .eq. -1) then
                if (fma_last(i) .gt. numl) then
                   write(lout,*) &
                        'ERROR in fma_postpr: Last turn in FMA block is larger than number of turns tracked fma_last=',&
                        fma_last(i),'> turns tracked=',numl,'!'
                endif
             else
                if (fma_last(i) .gt. dumplast(j)) then
                   write(lout,*) &
                        'ERROR in fma_postpr: Last turn in FMA block is larger than number of turns tracked'//&
                        ' in DUMP block fma_last=',fma_last(i),'> dumplast=',dumplast(j),'!'
                endif
             endif
             ! now we can set the number of turns used for the FMA required for the PLATO routines
             fma_nturn(i) = fma_last(i)-fma_first(i)+1
             do m=1,napx
                nturns(m)=fma_nturn(i)
             enddo
             if(fma_nturn(i).gt.fma_nturn_max) then
                write(lout,*) 'ERROR in fma_postpr: only ', fma_nturn_max,' turns allowed for fma and ',fma_nturn(i),' used!'
                write(lout,*) '->reset fma_nturn_max > ', fma_nturn_max
                call prror(-1)
             endif
             ! now we can start reading in the file
             if ( dumpfmt(j).eq.2 .or. dumpfmt(j) .eq. 7) then ! ASCII -> skip the header
                counter=1
                do
                   read(dumpunit(j),'(A)',iostat=ierro) ch
                   call fma_error(ierro,'while reading file '//dump_fname(j),'fma_postpr')
                   ch1=adjustl(trim(ch))
                   if(ch1(1:1).ne.'#')  exit
                   if(counter>500) then
                      call fma_error(ierro, &
                           "Something is wrong with your dumpfile '"//trim(stringzerotrim(dump_fname(j)))//&
                           "'; found more than 500 header lines",&
                           'fma_postpr')
                   endif
                   counter=counter+1
                enddo
                backspace(dumpunit(j),iostat=ierro)
             endif
             ! format 7 and 8 use normalized coordinates -> set fma_norm_flag =1
             if(dumpfmt(j) .eq. 7 .or. dumpfmt(j) .eq. 8) then
                if ( fma_norm_flag(i) .ne. 1 ) then
                   ! For format 7 and 8, the particles are already normalized by the DUMP block
                   write(lout,*) "ERROR in fma_postpr() for FMA #",i
                   write(lout,*) "Cannot do FMA on physical coordinates if normalized DUMP is used (dump format = 7 or 8)!"
                   call prror(-1)
                endif
             else ! Reading physical coordinates
                if ( fma_norm_flag(i) .eq. 1 ) then
                   ! Have a matrix that's not zero (i.e. did we put a 6d LINE block?)
                   if ( dumptas(j,1,1).eq.zero .and. &
                        dumptas(j,1,2).eq.zero .and. &
                        dumptas(j,1,3).eq.zero .and. &
                        dumptas(j,1,4).eq.zero         ) then
                      write(lout,*) "ERROR in FMA with normalized coordinates:"
                      write(lout,*) "The normalization matrix appears to not be set?"
                      write(lout,*) "Did you forget to put a 6D LINE block?"
                      call prror(-1)
                   endif
                   if(idp.eq.0 .or. ition.eq.0) then ! We're in the 4D case
                      if(imc.ne.1) then !Energy scan
                         write(lout,*) "ERROR in FMA with normalized coordinates:"
                         write(lout,*) "Energy scan (imc != 1) not supported!"
                         call prror(-1)
                      endif
                      if(j.ne.-1) then !Not at StartDUMP
                         write(lout,*) "ERROR in FMA with normalized coordinates:"
                         write(lout,*) "4D only supported for StartDUMP!"
                         call prror(-1)
                      endif
                   endif
                endif
             endif
             ! - now we have done all checks
             if ( fma_writeNormDUMP .and. .not.(dumpfmt(j).eq.7 .or. dumpfmt(j).eq.8) .and. .not.hasNormDumped(j) ) then
                write(lout,*) "FMA: Writing normalized DUMP for '"//trim(stringzerotrim(dump_fname(j)))// "'..."
                ! Dump normalized particle amplitudes for debugging (200101+i*10)
                inquire(unit=200101+i*10,opened=lopen)
                if(lopen) then
                   write(lout,*) "ERROR in FMA: Tried to open unit",200101+i*10, &
                        "for file 'NORM_"//dump_fname(j)//"', but it was already taken?!?"
                   call prror(-1)
                endif
#ifdef BOINC
                call boincrf('NORM_'//dump_fname(j),filename)
                open(200101+i*10,file=filename,              status='replace',iostat=ierro,action='write') ! nx,nx',ny,ny'
#else
                open(200101+i*10,file='NORM_'//dump_fname(j),status='replace',iostat=ierro,action='write') ! nx,nx',ny,ny'
#endif
                !  units: dumptas, dumptasinv, dumpclo [mm,mrad,mm,mrad,1]
                !  note: closed orbit dumpclo already converted in linopt part to [mm,mrad,mm,mrad,1]
                !        tas matrix in linopt part in [mm,mrad,mm,mrad,1.e-3]

                ! - write closed orbit in header of file with normalized phase space coordinates (200101+i*10)
                !   units: x,xp,y,yp,sig,dp/p = [mm,mrad,mm,mrad,1] (note: units are already changed in linopt part)
                write(200101+i*10,'(a,1x,6(1X,1PE16.9))') &
                     '# closorb', &
                     dumpclo(j,1),dumpclo(j,2),dumpclo(j,3), &
                     dumpclo(j,4),dumpclo(j,5),dumpclo(j,6)
                ! - write tas-matrix and its inverse in header of file with normalized phase space coordinates (200101+i*10)
                !   units: x,px,y,py,sig,dp/p [mm,mrad,mm,mrad,1]
                write(200101+i*10,'(a,1x,36(1X,1PE16.9))') &
                     '# tamatrix [mm,mrad,mm,mrad,1]', &
                     dumptas(j,1,1),dumptas(j,1,2),dumptas(j,1,3),dumptas(j,1,4), &
                     dumptas(j,1,5),dumptas(j,1,6),dumptas(j,2,1),dumptas(j,2,2), &
                     dumptas(j,2,3),dumptas(j,2,4),dumptas(j,2,5),dumptas(j,2,6), &
                     dumptas(j,3,1),dumptas(j,3,2),dumptas(j,3,3),dumptas(j,3,4), &
                     dumptas(j,3,5),dumptas(j,3,6),dumptas(j,4,1),dumptas(j,4,2), &
                     dumptas(j,4,3),dumptas(j,4,4),dumptas(j,4,5),dumptas(j,4,6), &
                     dumptas(j,5,1),dumptas(j,5,2),dumptas(j,5,3),dumptas(j,5,4), &
                     dumptas(j,5,5),dumptas(j,5,6),dumptas(j,6,1),dumptas(j,6,2), &
                     dumptas(j,6,3),dumptas(j,6,4),dumptas(j,6,5),dumptas(j,6,6)
                write(200101+i*10,'(a,1x,36(1X,1PE16.9))') &
                     '# inv(tamatrix)', &
                     dumptasinv(j,1,1),dumptasinv(j,1,2),dumptasinv(j,1,3), &
                     dumptasinv(j,1,4),dumptasinv(j,1,5),dumptasinv(j,1,6), &
                     dumptasinv(j,2,1),dumptasinv(j,2,2),dumptasinv(j,2,3), &
                     dumptasinv(j,2,4),dumptasinv(j,2,5),dumptasinv(j,2,6), &
                     dumptasinv(j,3,1),dumptasinv(j,3,2),dumptasinv(j,3,3), &
                     dumptasinv(j,3,4),dumptasinv(j,3,5),dumptasinv(j,3,6), &
                     dumptasinv(j,4,1),dumptasinv(j,4,2),dumptasinv(j,4,3), &
                     dumptasinv(j,4,4),dumptasinv(j,4,5),dumptasinv(j,4,6), &
                     dumptasinv(j,5,1),dumptasinv(j,5,2),dumptasinv(j,5,3), &
                     dumptasinv(j,5,4),dumptasinv(j,5,5),dumptasinv(j,5,6), &
                     dumptasinv(j,6,1),dumptasinv(j,6,2),dumptasinv(j,6,3), &
                     dumptasinv(j,6,4),dumptasinv(j,6,5),dumptasinv(j,6,6)

                write(200101+i*10,'(a)') &
                     '# id turn pos[m] nx[1.e-3 sqrt(m)]' //&
                     ' npx[1.e-3 sqrt(m)] ny[1.e-3 sqrt(m)]'//&
                     ' npy[1.e-3 sqrt(m)] nsig[1.e-3 sqrt(m)]'//&
                     ' ndp/p[1.e-3 sqrt(m)] kt'
             endif !END IF fma_writeNormDUMP

             ! Read in particle amplitudes a(part,turn), x,xp,y,yp,sigma,dE/E [mm,mrad,mm,mrad,mm,1]
             ! TODO: This logic breaks apart if there are particle losses;
             !  it is checked for, but it only triggers a "call prror(-1)".

             ! If normalization within FMA, we now have to always write the full NORM_* file
             ! Otherwise  one would overwrite the NORM_* file constantly if different FMAs are done
             ! on the same DUMP file

             if (dumplast(j) .eq. -1) then
                dump_last_turn = numl
             else
                dump_last_turn = dumplast(j)
             endif

             !Loop over all turns in the DUMP file;
             ! this is neccessary since we're writing normalized DUMP files.
             do k=dumpfirst(j),dump_last_turn !loop over turns, use the dump files

                !loop over particles
                do l=1,napx
                   if (dumpfmt(j).eq.2 .or. dumpfmt(j).eq.7) then  ! Read an ASCII dump
#ifndef CRLIBM
                      read(dumpunit(j),*,iostat=ierro) &
                           id,thisturn,pos,xyzvdummy(1),xyzvdummy(2),xyzvdummy(3),xyzvdummy(4),xyzvdummy(5),xyzvdummy(6),kt
                      if(ierro.gt.0) then
                         write(ch,'(a,1x,I5,1x,a)') &
                              "while reading  particles from file '"//trim(stringzerotrim(dump_fname(j)))// &
                              "' (dumpfmt=",dumpfmt(j),')'
                         call fma_error(ierro,ch,'fma_postpr') !read error
                      endif
#else
                      read(dumpunit(j),'(a)', iostat=ierro) ch
                      if(ierro.gt.0) then
                         ! read error
                         call fma_error(ierro,'while reading  particles from file '//trim(stringzerotrim(dump_fname(j))) // &
                              '. Check that tracked turns is larger than the number of turns used for FFT!', 'fma_postpr')
                      endif
                      call getfields_split(ch,filefields_fields,filefields_lfields,filefields_nfields,filefields_lerr)
                      ! error in getfields_split while reading
                      if( filefields_lerr ) then
                         call fma_error(-1,'while reading  particles from file '//trim(stringzerotrim(dump_fname(j))) // &
                              'in function getfields_split','fma_postpr')
                      endif
                      ! check if number of fields is correct
                      if( filefields_nfields  .ne. 10 ) then
                         write(lout,*) &
                              "ERROR in fma_postpr while reading particles from file '"//&
                              trim(stringzerotrim(dump_fname(j)))//&
                              "'. 10 fields expected from getfields_split, got ", filefields_nfields, ' and ch=',ch
                         call prror(-1)
                      endif

                      read(filefields_fields(1) (1:filefields_lfields(1)),*) id
                      read(filefields_fields(2) (1:filefields_lfields(2)),*) thisturn
                      pos = round_near(ierro, filefields_lfields(3)+1, filefields_fields(3) )
                      if (ierro.ne.0) call rounderr(ierro, filefields_fields, 3, pos)
                      xyzvdummy(1) = round_near(ierro, filefields_lfields(4)+1, filefields_fields(4) )
                      if (ierro.ne.0) call rounderr(ierro, filefields_fields,4, xyzvdummy(1))
                      xyzvdummy(2) = round_near(ierro, filefields_lfields(5)+1, filefields_fields(5) )
                      if (ierro.ne.0) call rounderr(ierro, filefields_fields,5, xyzvdummy(2))
                      xyzvdummy(3) = round_near(ierro, filefields_lfields(6)+1, filefields_fields(6) )
                      if (ierro.ne.0) call rounderr(ierro,filefields_fields,6, xyzvdummy(3))
                      xyzvdummy(4) = round_near(ierro, filefields_lfields(7)+1, filefields_fields(7) )
                      if (ierro.ne.0) call rounderr(ierro,filefields_fields,7, xyzvdummy(4))
                      xyzvdummy(5) = round_near(ierro, filefields_lfields(8)+1, filefields_fields(8) )
                      if (ierro.ne.0) call rounderr(ierro,filefields_fields,8, xyzvdummy(5))
                      xyzvdummy(6) = round_near(ierro, filefields_lfields(9)+1, filefields_fields(9) )
                      if (ierro.ne.0) call rounderr(ierro,filefields_fields,9, xyzvdummy(6))
                      read(filefields_fields(10) (1:filefields_lfields(10)),*) kt
#endif

                   else if (dumpfmt(j).eq.3 .or. dumpfmt(j).eq.8) then ! Read a binary dump
                      read(dumpunit(j),iostat=ierro) &
                           id,thisturn,pos,xyzvdummy(1),xyzvdummy(2),xyzvdummy(3),xyzvdummy(4),xyzvdummy(5),xyzvdummy(6),kt
                      if(ierro.gt.0) then
                         write(ch,'(a,1x,I5,1x,a)') &
                              "while reading particles from file '"// trim(stringzerotrim(dump_fname(j))) // &
                              "' (dumpfmt=",dumpfmt(j),')'
                         call fma_error(ierro,ch,'fma_postpr') !read error
                      endif
                   endif

                   !Check for losses
                   if (l.ne.id .or. k.ne.thisturn) then
                      if (k .lt. nturns(l)+fma_first(l)-1) then
                         nturns(l) = k-fma_first(l)
                      endif

                      !TODO: Actually handle those losses.
                      write(lout,*) "ERROR when reading DUMP file #",j, "for FMA #",i
                      write(lout,*) "Expected turn and particle ID =", k,l
                      write(lout,*) "Got turn and particle ID =", thisturn,id
                      write(lout,*) "Reading probably got unsynchronized because of particle losses,"//&
                           " which is currently not handled in FMA."
                      call prror(-1)
                   endif

                   !Normalization
                   if (dumpfmt(j).eq.2 .or.dumpfmt(j).eq.3) then
                      ! Case: The file isn't pre-normalized -> Compute normalization
                      !
                      ! At this point fma_norm_flag doesn't matter;
                      ! we anyway compute the normalized coordinates.
                      !
                      ! units: dumptas, dumptasinv, dumpclo [mm,mrad,mm,mrad,1]

                      ! remove closed orbit -> check units used in dumpclo (is x' or px used?)
                      do m=1,6
                         xyzvdummy2(m)=xyzvdummy(m)-dumpclo(j,m)
                      enddo

                      !For use in with normalized coordinates:
                      ! convert to canonical variables
                      xyzvdummy2(2)=xyzvdummy2(2) * ((one+xyzvdummy2(6))+dumpclo(j,6))
                      xyzvdummy2(4)=xyzvdummy2(4) * ((one+xyzvdummy2(6))+dumpclo(j,6))

                      ! normalize nxyz=dumptasinv*xyz2
                      do m=1,6
                         nxyzvdummy(m)=zero
                         do n=1,6
                            nxyzvdummy(m)=nxyzvdummy(m) + dumptasinv(j,m,n)*xyzvdummy2(n)
                         enddo
                         ! convert nxyzvdummy(6) to 1.e-3 sqrt(m)
                         ! unit: nx,npx,ny,npy,nsig,ndelta all in [1.e-3 sqrt(m)]
                         if(m.eq.6) then
                            nxyzvdummy(m)=nxyzvdummy(m)*c1e3
                         endif
                      enddo

                      ! Write normalized particle amplitudes
                      ! (only when reading physical coordinates)
                      if (fma_writeNormDUMP .and. .not.hasNormDumped(j) ) then
                         write(200101+i*10, '(2(1x,I8),1X,F12.5,6(1X,1PE16.9),1X,I8)') &
                              id,thisturn,pos,nxyzvdummy(1),nxyzvdummy(2),nxyzvdummy(3),nxyzvdummy(4),nxyzvdummy(5),nxyzvdummy(6),kt
                      endif

                   else if (dumpfmt(j).eq.7 .or. dumpfmt(j).eq.8) then
                      ! Case: we are already normalized;
                      ! just copy the data into the relevant array
                      do m=1,6
                         nxyzvdummy(m) = xyzvdummy(m)
                      end do
                   endif ! END IF already normalized or not

                   ! Copy the data into the final arrays
                   if (thisturn.ge.fma_first(i) .and. thisturn.le.fma_last(i) ) then

                      turn(l,k-fma_first(i)+1) = thisturn

                      do m=1,6
                         ! for FMA in physical coordinates, convert units to [mm,mrad,mm,mrad,mm,1.e-3]
                         if(m.eq.6) then
                            xyzv(l,k-fma_first(i)+1,m)=xyzvdummy(m)*c1e3
                         else
                            xyzv(l,k-fma_first(i)+1,m)=xyzvdummy(m)
                         endif

                         nxyzv(l,k-fma_first(i)+1,m) = nxyzvdummy(m)

                         ! calculate emittance of mode 1,2,3
                         if(mod(m,2).eq.0) then
                            epsnxyzv(l,k-fma_first(i)+1,m/2) = nxyzvdummy((m-1))**2+nxyzvdummy(m)**2
                         endif
                      enddo

                   endif !END if fma_first <= thisturn <= fma_last

                enddo ! END loop over particles l
             enddo ! END loop over turns k

             ! Calculate tunes of particles using the methods in plato_seq.f and NAFF
             !  for fma_norm_flag == 0: use physical coordinates x,x',y,y',sig,dp/p
             !  for fma_norm_flag == 1: use normalized coordinates
             do l=1,napx ! loop over particles
                !TODO particle losses - detect if nturns(l) is too small & skip that particle.
                ! (probably just write a line of mostly zeros to the file)

                do m=1,num_modes ! loop over modes (hor.,vert.,long.)
                   select case( trim(stringzerotrim(fma_method(i))) )
                   case('TUNELASK')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) =  tunelask( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) =  tunelask(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

                   case('TUNEFFTI')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) =  tuneffti( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) =  tuneffti(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

                   case('TUNEFFT')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) =   tunefft( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) =   tunefft(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

                   case('TUNEAPA')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) =   tuneapa( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) =   tuneapa(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

                   case('TUNEFIT')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) =   tunefit( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) =   tunefit(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

                   case('TUNENEWT')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) =  tunenewt( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) =  tunenewt(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

                   case('TUNEABT2')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) =  tuneabt2( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) =  tuneabt2(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

                   case('TUNEABT')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) =   tuneabt( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) =   tuneabt(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

                   case('TUNENEWT1')
                      if(fma_norm_flag(i) .eq. 0) then
                         q123(m) = tunenewt1( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
                      else
                         q123(m) = tunenewt1(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
                      endif

#ifdef NAFF
                   case("NAFF")
                      ! write(lout,*) "DBG", fma_nturn(i),l
                      ! write(lout,*) "DBG", nxyzv(l,1,2*(m-1)+1), nxyzv(l,1,2*m)
                      !
                      ! write(lout,*) size(xyzv(l,fma_first(i):fma_last(i),2*(m-1)+1))
                      ! write(lout,*) size(xyzv(l,fma_first(i):fma_last(i),2*m))

                      flush(lout)  ! F2003 does specify a FLUSH statement.
                      ! However NAFF should NOT be chatty...

                      ! do n=1,fma_nturn(i)
                      !    write(*,*) n, nxyzv(l,n,2*(m-1)+1), nxyzv(l,n,2*m)
                      ! enddo
                      ! write(*,*) ""

                      ! Copy the relevant contents of the arrays
                      ! into a new temporary array with stride=1
                      ! for passing to C++.
                      if(fma_norm_flag(i) .eq. 0) then
                         naff_xyzv1 = xyzv (l, 1:nturns(l), 2*(m-1)+1)
                         naff_xyzv2 = xyzv (l, 1:nturns(l), 2*m)
                      else
                         naff_xyzv1 = nxyzv(l,1:nturns(l),  2*(m-1)+1)
                         naff_xyzv2 = nxyzv(l,1:nturns(l),  2*m)
                      endif

                      q123(m) = tunenaff(naff_xyzv1, naff_xyzv2, nturns(l), m, fma_norm_flag(i) )

                      flush(lout)
                      ! stop
#endif

                   case default
                      call fma_error(-1,'FMA method '//trim(stringzerotrim(fma_method(i)))// &
                           ' not known! Note that the method name must be in capital letters!','fma_postpr')
                   end select

                   ! mode 3 rotates anticlockwise, mode 1 and 2 rotate clockwise -> synchroton tune is negative,
                   ! but define it as convention positive
                   if(m.eq.3) q123(m)=one-q123(m)

                   !Some general calculations
                   eps123_0(m)=epsnxyzv(l,1,m) ! initial amplitude
                   phi123_0(m)=atan_mb(nxyzv(l,1,2*m)/nxyzv(l,1,2*(m-1)+1))! inital phase
                   eps123_min(m)=minval( epsnxyzv(l,1:nturns(l),m) )       ! minimum emittance
                   eps123_max(m)=maxval( epsnxyzv(l,1:nturns(l),m) )       ! maximum emittance
                   eps123_avg(m)=sum(epsnxyzv(l,1:nturns(l),m))/nturns(l)  ! average emittance
                enddo
                if ( num_modes .eq. 2 ) then
                   q123(3)=zero
                   eps123_min(3)=zero
                   eps123_max(3)=zero
                   eps123_avg(3)=zero
                   eps123_0(3)=zero
                   phi123_0(3)=zero
                endif

                ! Write the FMA output file "fma_sixtrack"
                ! TODO losses: fma_first and fma_last may not be the right start/stop variables...
                write(2001001,'(2(1x,A20),1x,I8,18(1X,1PE16.9),3(1X,I8))') &
                     trim(stringzerotrim(fma_fname(i))), &
                     trim(stringzerotrim(fma_method(i))),l,q123(1),q123(2),q123(3), &
                     eps123_min(1),eps123_min(2),eps123_min(3),eps123_max(1), &
                     eps123_max(2),eps123_max(3),eps123_avg(1),eps123_avg(2), &
                     eps123_avg(3),eps123_0(1),eps123_0(2),eps123_0(3), &
                     phi123_0(1),phi123_0(2),phi123_0(3),fma_norm_flag(i), &
                     fma_first(i),fma_last(i)

             enddo ! END loop over particles l

             if (fma_writeNormDUMP .and. .not.hasNormDumped(j)) then
                ! filename NORM_* (normalized particle amplitudes)
                close(200101+i*10)
                hasNormDumped(j) = .true.
             endif

             ! resume initial position of dumpfile = end of file
             close(dumpunit(j))
             if (dumpfmt(j).eq.2 .or. dumpfmt(j).eq. 7) then !ASCII
#ifdef BOINC
                call boincrf(dump_fname(j),filename)
                open(dumpunit(j),file=filename,      status='old', form='formatted', &
                     action='readwrite', position='append', iostat=ierro)
#else
                open(dumpunit(j),file=dump_fname(j), status='old', form='formatted', &
                     action='readwrite', position='append', iostat=ierro)
#endif
                call fma_error(ierro, "while resuming file '"// &
                     trim(stringzerotrim(dump_fname(j)))//"' (dumpfmt=2)", 'fma_postpr')
             elseif (dumpfmt(j).eq.3 .or. dumpfmt(j).eq.8) then !BINARY
#ifdef BOINC
                call boincrf(dump_fname(j),filename)
                open(dumpunit(j),file=filename,      status='old', form='unformatted', &
                     action='readwrite', position='append', iostat=ierro)
#else
                open(dumpunit(j),file=dump_fname(j), status='old', form='unformatted', &
                     action='readwrite', position='append', iostat=ierro)
#endif
                call fma_error(ierro, "while resuming file '"// &
                     trim(stringzerotrim(dump_fname(j)))//"' (dumpfmt=3)", 'fma_postpr')
             endif
          endif !END: if fma_fname(i) matches dump_fname(j)

          ! if file has been already found, jump to next file fma_fname(i)
          if( lexist ) then
             exit
          endif
       enddo !END: loop over dump files
       if(.not. lexist) then !if no dumpfile has been found, raise error and abort
          call fma_error(-1,'dump file '//trim(stringzerotrim(fma_fname(i)))//&
               ' does not exist! Please check that filenames in FMA block agree with the ones in the DUMP block!', 'fma_postpr')
       endif
    enddo !END: loop over fma files
    close(2001001) !filename: fma_sixtrack

    deallocate(turn, nturns, xyzv, nxyzv, epsnxyzv)
#ifdef NAFF
    deallocate(naff_xyzv1, naff_xyzv2)
#endif
  end subroutine fma_postpr

end module fma
