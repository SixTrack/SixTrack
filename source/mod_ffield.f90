! ================================================================================================ !
! INTERFACE BETWEEN SIXTRACT AND FRINGE-FIELD ROUTINE
! B.Dalena, T.Pugnat and A.Simona
! Last modified: 2019-01-18
! ================================================================================================ !
module mod_ffield
  ! ------------------------------------------------------------------------------------------------ !
  ! Mod from SixTrack
  ! ------------------------------------------------------------------------------------------------ !
!  use floatPrecision

  use ffTable_n_Tracks
  ! ------------------------------------------------------------------------------------------------ !

  implicit none

  
  ! ------------------------------------------------------------------------------------------------ !
  ! FFIELD Function
  ! ------------------------------------------------------------------------------------------------ !
  public  :: ffield_mod_init
  public  :: ffield_mod_expand_arrays
  public  :: ffield_mod_expand_arrays_study
  public  :: ffield_mod_end
  public  :: ffield_parseInputLine
  public  :: ffield_parsingDone
  public  :: ffield_mod_link
  public  :: ffield_genAntiQuad
  public  :: ffield_enterQuad
  public  :: ffield_exitQuad

  private :: ffield_mod_ChckQuad


  ! ------------------------------------------------------------------------------------------------ !
  ! FFIELD parameters
  ! ------------------------------------------------------------------------------------------------ !
  integer, public, save :: ff_status              ! Status routine: 0=off, 1=init, 2=loaded, 3=AQ ready
  integer, public, save :: ffNLn, ffNLn_max       ! Nb. and nb. max of Quadrupoles with ffNLn < ffNLn_max
  integer, public, save :: ffMSn, ffMSn_max       ! Nb. and nb. max of Multipole to skip with ffMSn < ffMSn_max
  integer, public, save :: ffNLFile, ffNLFile_max ! Nb. and nb. max of files with ffNLFile < ffNLFile_max
!  real(kind=fPrec), public, save :: ffdelta         ! Max ffdelta posible 

  real(kind=fPrec) :: r0_2 
  parameter (r0_2=6.4e-3_fPrec)                   ! Maximum radius in quad
  ! ------------------------------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------------------------------ !
  ! FFIELD table
  ! ------------------------------------------------------------------------------------------------ !
  logical         , allocatable, private, save :: ffInQuad(:)  ! (1:npart)        Check if particle enter the Quad
  integer         , allocatable, public, save :: ffindex(:)    ! (0:nele)         Table with the index of the Quad in our study (0 = not studied)
  integer         , allocatable, public, save :: ffQ2File(:,:) ! (1:ffNLn, 1:2)   Link Quad/Files
!  real(kind=fPrec), allocatable, public, save :: ffParam(:,:)  ! (1:ffNLFile,1:6) Kin, Lin, Corin, Kex, Lex, Corex
  real(kind=fPrec), allocatable, public, save :: ffParam(:,:)  ! (1:ffNLFile,1:2) Lgth. in Quadrupoles and total lgth.
  character(len=:), allocatable, public, save :: ffQNames(:)   ! (1:ffNLn)        Name of Quadrupoles
  character(len=:), allocatable, public, save :: ffMSNames(:)  ! (1:ffMSn)        Name of Multipoles skip
  character(len=:), allocatable, public, save :: ffFNames(:)   ! (1:ffNLFile)     Name of Files

  type(ffTable_n_Track), allocatable, public, save :: ffTable(:)  ! (1:ffNLFile)
  ! ------------------------------------------------------------------------------------------------ !
  
  contains


  ! ================================================================================================ !
  !  Init the module
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2018-10-18
  ! ================================================================================================ !
  subroutine ffield_mod_init(npart, nele)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use mod_alloc,           only : alloc
    use numerical_constants, only : zero
    use parpro,              only : mNameLen, str_nmSpace, mStrLen, str_dSpace

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    integer, intent(in) :: npart, nele
    
    
    ! initialization of the variable
    ! ---------------------------------------------------------------------------------------------- !
    ff_status   =0
    ffNLn       =1
    ffNLn_max   =20
    ffMSn       =1
    ffMSn_max   =20
    ffNLFile    =1
    ffNLFile_max=20

    
    ! allocation of the memory
    ! ---------------------------------------------------------------------------------------------- !
    call alloc(ffInQuad, npart,           .false., 'ffInQuad')
    call alloc(ffindex,  nele,            0,       'ffindex')
    call alloc(ffQ2File, ffNLFile_max, 2, 0,       'ffQ2File')
    call alloc(ffParam,  ffNLn_max,    2, zero,    'ffParam')

    call alloc(ffQNames,  mNameLen, ffNLn_max,    str_nmSpace, 'ffQNames')
    call alloc(ffMSNames, mNameLen, ffMSn_max,    str_nmSpace, 'ffMSNames')
    call alloc(ffFNames,  mStrLen,  ffNLFile_max, str_dSpace,  'ffFNames')

  end subroutine ffield_mod_init



  ! ================================================================================================ !
  !  Expand array in module common between SixTrack and our Study
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2018-10-18
  ! ================================================================================================ !
  subroutine ffield_mod_expand_arrays(npart_new, nele_new)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use mod_alloc,           only : resize

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    integer, intent(in) :: npart_new, nele_new

    ! resizing of the memory
    ! ---------------------------------------------------------------------------------------------- !
    call resize(ffInQuad, npart_new, .false., 'ffInQuad')
    call resize(ffindex,  nele_new,  0,       'ffindex')

  end subroutine ffield_mod_expand_arrays



  ! ================================================================================================ !
  !  Expand array in module only for the study
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2018-10-18
  ! ================================================================================================ !
  subroutine ffield_mod_expand_arrays_study(ffNLn_max_new, ffMSn_max_new, ffNLFile_max_new)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use mod_alloc,           only : resize
    use numerical_constants, only : zero
    use parpro,              only : mNameLen, str_nmSpace, mStrLen, str_dSpace

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    integer, intent(in) :: ffNLn_max_new, ffMSn_max_new, ffNLFile_max_new

    ! resizing of the memory
    ! ---------------------------------------------------------------------------------------------- !
    if (ffNLn_max /= ffNLn_max_new) then
      call resize(ffParam,            ffNLn_max_new, 2, zero,        'ffParam')
      call resize(ffQNames, mNameLen, ffNLn_max_new,    str_nmSpace, 'ffQNames')

      ffNLn_max = ffNLn_max_new

    endif

    if (ffMSn_max /= ffMSn_max_new) then
      call resize(ffMSNames, mNameLen, ffMSn_max_new,    str_nmSpace, 'ffMSNames')

      ffMSn_max = ffMSn_max_new

    endif

    if (ffNLFile_max /= ffNLFile_max_new) then
      call resize(ffQ2File,           ffNLFile_max_new, 2, 0,           'ffQ2File')
      call resize(ffFNames, mStrLen , ffNLFile_max_new,    str_dSpace , 'ffFNames')

      ffNLFile_max = ffNLFile_max_new

    endif

  end subroutine ffield_mod_expand_arrays_study



  ! ================================================================================================ !
  !  Close and clean memory in module
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2018-10-18
  ! ================================================================================================ !
  subroutine ffield_mod_end()
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use mod_alloc,           only : dealloc

    implicit none


    ! routine variables
    ! ---------------------------------------------------------------------------------------------- !
    integer :: i  ! iterator

    ! free the memory
    ! ---------------------------------------------------------------------------------------------- !
    if (allocated(ffInQuad))  call dealloc(ffindex,  'ffInQuad')
    if (allocated(ffindex))   call dealloc(ffindex,  'ffindex')
    if (allocated(ffQ2File))  call dealloc(ffQ2File, 'ffQ2File')
    if (allocated(ffParam))   call dealloc(ffParam,  'ffParam')


    if (allocated(ffTable)) then
      do i=1,ffNLFile
        call ffTable(i)%free()
      end do
      deallocate(ffTable)
    end if

    if (allocated(ffQNames))  deallocate(ffQNames) !,  'ffQNames')
    if (allocated(ffMSNames)) deallocate(ffMSNames) !, 'ffMSNames')
    if (allocated(ffFNames))  deallocate(ffFNames) !,  'ffFNames')
  end subroutine ffield_mod_end



  ! ================================================================================================ !
  !  Parse FField Input Line
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2019-01-17
  ! ================================================================================================ !
  subroutine ffield_parseInputLine(inLine, iLine, iErr)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use string_tools, only : chr_rpad, chr_split, chr_cast
    use mod_common,   only : bez
    use parpro,       only : mNameLen, mStrLen

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    character(len=*), intent(in)    :: inLine
    integer,          intent(in)    :: iLine
    logical,          intent(inout) :: iErr

    ! routine variables
    ! ---------------------------------------------------------------------------------------------- !
    character(len=:), allocatable   :: lnSplit(:)
    character(len=:), allocatable   :: elemName
    integer :: nSplit, i
    logical :: spErr

    ! Routine start
    ! ---------------------------------------------------------------------------------------------- !
    if (ff_status==0) then
      write(lout,"(a)") "FFIELD> Start to read the input."
      ff_status=1
    endif

    ! cut the string
    ! ---------------------------------------------------------------------------------------------- !
    call chr_split(inLine, lnSplit, nSplit, spErr)
    if(spErr) then
      write(lout,"(a)") "FFIELD> ERROR Failed to parse input line."
      iErr = .true.
      return
    end if

    ! select case
    ! ---------------------------------------------------------------------------------------------- !
    select case(lnSplit(1)(1:4))
  
    case("FFQN")   ! FFQN,    Quad name,     no file start,     no file end
      !     * Name a Quadrupole for the study
      ! -------------------------------------------------------------------------------------------- !
      if(nSplit < 4) then
        write(lout,"(a,i0)") "FFIELD> ERROR FFQN line must have at less 4 values, got ",nSplit
        iErr = .true.
        return
      end if

      elemName = chr_rpad(lnSplit(2), mNameLen)

      ! Check that the name is unique
      do i=1,ffNLn
        if(ffQNames(i) == elemName) then
          write(lout,"(a,i0)") "FFIELD> ERROR Quadrupole element '"//trim(elemName)//"' is not unique."
          iErr = .true.
          return
        end if
      end do
      
      ffQNames(ffNLn) = elemName

      call chr_cast(lnSplit(3),ffQ2File(ffNLn,1),  iErr)
      call chr_cast(lnSplit(4),ffQ2File(ffNLn,2),  iErr)
      if(iErr) return

      ffNLn = ffNLn+1
      
      if (ffNLn>=ffNLn_max-1) then
        call ffield_mod_expand_arrays_study(ffNLn_max+10, ffMSn_max, ffNLFile_max)
      end if
      
            
    case("FFMS") 
      !     * Name a Multipole that will be skip in the study
      ! -------------------------------------------------------------------------------------------- !
      if(nSplit < 2) then
        write(lout,"(a,i0)") "FFIELD> ERROR FFMS line must have at less 2 values, got ",nSplit
        iErr = .true.
        return
      end if

      elemName = chr_rpad(lnSplit(2), mNameLen)

      ! Check that the name is unique
      do i=1,ffMSn
        if(ffMSNames(i) == elemName) then
          write(lout,"(a,i0)") "FFIELD> ERROR Multipole to skip element '"//trim(elemName)//"' is not unique."
          iErr = .true.
          return
        end if
      end do
      
      ffMSNames(ffMSn) = elemName

      ffMSn = ffMSn+1
      
      if (ffMSn>=ffMSn_max-1) then
        call ffield_mod_expand_arrays_study(ffNLn_max, ffMSn_max+10, ffNLFile_max)
      end if
      
            
    case("FFFI")   ! FFFI,    File name,     Length equi.,      Total length file
      !     * Name a File for the study
      ! -------------------------------------------------------------------------------------------- !
      if(nSplit < 4) then
        write(lout,"(a,i0)") "FFIELD> ERROR FFFI line must have at less 4 values, got ",nSplit
        iErr = .true.
        return
      end if
      
      elemName = chr_rpad(lnSplit(2), mStrLen)
      ! Check that the name is unique
      do i=1,ffNLFile
        if(ffFNames(i) == elemName) then
          write(lout,"(a,i0)") "FFIELD> ERROR File '"//trim(elemName)//"' is not unique."
          iErr = .true.
          return
        end if
      end do
      
      ffFNames(ffNLFile) = elemName
      call chr_cast(lnSplit(3),ffParam(ffNLFile,1),  iErr) ! Lin
      call chr_cast(lnSplit(4),ffParam(ffNLFile,2),  iErr) ! Lgth
      if(iErr) return

      ffNLFile = ffNLFile+1
      
      if (ffNLFile>=ffNLFile_max) then
        call ffield_mod_expand_arrays_study(ffNLn_max, ffMSn_max, ffNLFile_max+10)
      end if
      
    end select
  
  end subroutine ffield_parseInputLine



  ! ================================================================================================ !
  !  Parse FField Done
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2019-01-16
  ! ================================================================================================ !
  subroutine ffield_parsingDone()
    write(lout,"(a,i0)") "FFIELD> All input ready!"
  end subroutine ffield_parsingDone





  ! ================================================================================================ !
  !  Init the module
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2019-01-18
  ! ================================================================================================ !
  subroutine ffield_mod_link(iErr)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use mod_common,          only : e0
    use physical_constants,  only : clight,pmap
    use numerical_constants, only : one, c1e6

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    logical, intent(inout) :: iErr

    ! routine variables
    ! ---------------------------------------------------------------------------------------------- !
    integer :: i
    real(kind=fPrec) :: norm=1
    real(kind=fPrec) :: beta0,gamma0r
    
    if (ff_status==1) then
      ffNLn = ffNLn-1
      ffMSn = ffMSn-1
      ffNLFile = ffNLFile-1
      gamma0r=pmap/e0 !      gamma0=e0/pmap
      beta0=sqrt(one-gamma0r*gamma0r) ! = sqrt(one-one/(gamma0*gamma0))
!      p0=beta0*e0*c1e6/clight
      norm=clight/(beta0*e0*c1e6)   ! [c/eV] = 1/p0

#ifdef DEBUG
      ! Debug
      ! -------------------------------------------------------------------------------------------- !
      do i=1,ffNLn
        write(lout,"(a,i0)") "FFIELD> DEBUG Quadrupole element '"//trim(ffQNames(i))//"' will be studied."
      end do

      do i=1,ffNLFile
        write(lout,"(a,i0)") "FFIELD> DEBUG File '"//trim(ffFNames(i))//"' will be loaded."
      end do
#endif
      
      ! Generate the array of type(ffTable_n_Track)
      ! -------------------------------------------------------------------------------------------- !
      allocate(ffTable(1:ffNLFile))
      do i=1,ffNLFile
        call ffTable(i)%set(trim(ffFNames(i)),ffParam(i,1),ffParam(i,2))
      end do
      
      ! Check if the Quadrupole ask for the study is in the lattice
      ! -------------------------------------------------------------------------------------------- !
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      write(lout,"(a,i0)") "FFIELD> |      Summary of the quadrupole for the Fringe Field study      |"
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      do i=1,ffNLn
        if ((ffQ2File(i,1)>ffNLn).or.(ffQ2File(i,1)<1).or.(ffQ2File(i,2)>ffNLn).or.(ffQ2File(i,2)<1)) then
          write(lout,"(a,i0)") "FFIELD> ERROR FFQN: Wrong choise of file for the Quadupole's head. Check '"//trim(ffQNames(i))//"'."
          iErr = .true.
          return
        end if
        
        call ffield_mod_ChckQuad(i,norm,iErr)
        
      end do
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      write(lout,"(a,i0)") "FFIELD>"
      
      ! Check if the multipole skip for the study is in the lattice
      ! -------------------------------------------------------------------------------------------- !
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      write(lout,"(a,i0)") "FFIELD> |      Summary of the multipole for the Fringe Field skipped     |"
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      do i=1,ffMSn
        call ffield_mod_ChckMulti(i,iErr)
      end do
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      write(lout,"(a,i0)") "FFIELD>"

      ! Check if the Quadrupole ask for the study is in the lattice
      ! -------------------------------------------------------------------------------------------- !
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      write(lout,"(a,i0)") "FFIELD> |         Summary of the file for the Fringe Field study         |"
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      do i=1,ffNLFile
        if (ffTable(i)%chk_Status == 2) then 
          write(lout,"(a)") "FFIELD>   * '"//ffTable(i)%ffFNames//"' loaded!"
        else if (ffTable(i)%chk_Status == 1) then
          write(lout,"(a)") "FFIELD>   * '"//ffTable(i)%ffFNames//"' waiting!"
        end if
      end do
      write(lout,"(a,i0)") "FFIELD> +----------------------------------------------------------------+"
      
      ff_status=2
    end if

  end subroutine ffield_mod_link





  ! ================================================================================================ !
  !  Load file if Quadrupole exist in the lattice
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2019-01-18
  ! ================================================================================================ !
  subroutine ffield_mod_ChckQuad(i,norm,iErr)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use parpro    ,   only : nele, npart
    use mod_common,   only : bez

    implicit none

    ! Interface variables
    ! ---------------------------------------------------------------------------------------------- !
    logical,          intent(inout) :: iErr
    integer,          intent(in)    :: i
    real(kind=fPrec), intent(in)    :: norm 

    ! Subroutine variables
    ! ---------------------------------------------------------------------------------------------- !
    integer :: j           ! Iterator
    
    ! Check quadrupole in the lattice
    ! ---------------------------------------------------------------------------------------------- !
    do j=1,nele
      if (bez(j) == ffQNames(i)) then
        ffindex(j)=i
        call ffTable(ffQ2File(i,1))%load(norm,iErr)
        call ffTable(ffQ2File(i,2))%load(norm,iErr)
        write(lout,"(a,i0)") "FFIELD>   * '"//trim(ffQNames(i))//"' ready!"
        return
      end if
    end do

    write(lout,"(a,i0)") "FFIELD>   * '"//trim(ffQNames(i))//"' ignored!"
  end subroutine ffield_mod_ChckQuad





  ! ================================================================================================ !
  !  Load file if Multipole exist in the lattice
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2019-01-18
  ! ================================================================================================ !
  subroutine ffield_mod_ChckMulti(i,iErr)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use parpro    ,   only : nele, npart
    use mod_common,   only : bez

    implicit none

    ! Interface variables
    ! ---------------------------------------------------------------------------------------------- !
    logical,          intent(inout) :: iErr
    integer,          intent(in)    :: i

    ! Subroutine variables
    ! ---------------------------------------------------------------------------------------------- !
    integer :: j           ! Iterator
    
    ! Check Multipole not in FFQN
    ! ---------------------------------------------------------------------------------------------- !
    do j=1,ffNLn
      if (ffQNames(j)==ffMSNames(i)) then
        write(lout,"(a,i0)") "FFIELD> ERROR: Multipole cannot be in 'FFQN' AND 'FFMS'. Check '"//trim(ffMSNames(i))//"'."
        iErr=.true.
        return
      end if
    end do
    
    ! Check multipole in the lattice
    ! ---------------------------------------------------------------------------------------------- !
    do j=1,nele
      if (bez(j) == ffMSNames(i)) then
        ffindex(j)=-i
        write(lout,"(a,i0)") "FFIELD>   * '"//trim(ffMSNames(i))//"' will be skipped!"
        return
      end if
    end do

    write(lout,"(a,i0)") "FFIELD>   * '"//trim(ffMSNames(i))//"' not found!"
!    end if
  end subroutine ffield_mod_ChckMulti





  ! ================================================================================================ !
  !  Generate the AntiQuadrupole matrix 
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2019-01-21
  ! ================================================================================================ !
  subroutine ffield_genAntiQuad()
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use mod_common,          only : napx
    use mod_commonmn,        only : dpsv
    use numerical_constants, only : one, two, c1e7

    implicit none

    ! Subroutine variables
    ! ---------------------------------------------------------------------------------------------- !
    real(kind=fPrec) :: ffdelta 
    integer :: iFile, j        ! Iterator
    integer :: nbDlt           ! Nb of point in the interval [-ffdelta,ffdelta]
    
    if (ff_status==2) then
      ffdelta=-one
      do j=1,napx
        if (ffdelta<abs(dpsv(j))) ffdelta=abs(dpsv(j))
      enddo

      nbDlt=max(1,ceiling(abs( two*(ffdelta)*c1e7 )))
      do iFile=1,ffNLFile
        call ffTable(iFile)%AQgen(ffdelta,nbDlt)
      end do

      ff_status=3
    end if
    
  end subroutine ffield_genAntiQuad





  ! ================================================================================================ !
  !  Lie 2 Tracking if particles enter the Quadupole
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2019-01-21
  ! ================================================================================================ !
  subroutine ffield_enterQuad(ffi)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use parpro,              only : nblo
    use mod_common,          only : napx, ic, tiltc, tilts
    use mod_commont,         only : strack
    use mod_commonmn,        only : xv, yv, oidpsv, dpsv, rvv, ejv, zsiv, xsiv, nlostp
    use numerical_constants, only : half, one, c1e3, c1m3, c1m6

    implicit none

    ! Interface variables
    ! ---------------------------------------------------------------------------------------------- !
    integer, intent(in):: ffi

    ! Subroutine variables
    ! ---------------------------------------------------------------------------------------------- !
    integer          :: iFile
    integer          :: ffj,  k                    ! iterator
    integer          :: itDlt                      ! 
    real(kind=fPrec) :: x, px, y, py               ! Transverse canonical parameter in new referencial
    real(kind=fPrec) :: x_tp, px_tp, y_tp, py_tp   ! Temp transverse canonical parameter

    real(kind=fPrec) :: ffdelta!, gam0, betabeta0     ! 
    real(kind=fPrec) :: zb!, sigma_s                ! 
    real(kind=fPrec) :: LoutQ!, Ldpsv1, Ldpsv2      ! 


    iFile=ffQ2File(ffindex(ic(ffi)-nblo),1)

! <<<<<<<<<<<<<<<<< Debug
!write(lout,*)"FFIELD> Debug ************ IN ************"
!write(lout,*)"FFIELD> Debug [etp ffj]  ->  x, px, y, py"
! <<<<<<<<<<<<<<<<< Debug
    do ffj=1,napx
      ! Save data
      ! -------------------------------------------------------------------------------------------- !
      x = xv(1,ffj);  y = xv(2,ffj);  px= yv(1,ffj);  py= yv(2,ffj)

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [0, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug

!  <<<<<<< IN
      ! Return to beginning of the Quad
      ! -------------------------------------------------------------------------------------------- !
      !             * HE
      x_tp = x - strack(ffi+2)*half*px
      y_tp = y - strack(ffi+2)*half*py
      

      !             * Riccardo's Quad head
!      x_tp = x - strack(ffi-1)*px
!      y_tp = y - strack(ffi-1)*py

      x=x_tp;   y=y_tp;
!  <<<<<<< IN

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [1, ",ffj,"]  ->  ",x, px, y, py,"   -> strack(ffi+2)*half*px",strack(ffi+2)*half*px
!endif
! <<<<<<<<<<<<<<<<< Debug


      ! Rotation error for the quad
      ! -------------------------------------------------------------------------------------------- !
      x_tp = ((x-xsiv(1,ffi))*tiltc(ffi) + (y-zsiv(1,ffi))*tilts(ffi))*c1m3! mm -> m
      y_tp = ((y-zsiv(1,ffi))*tiltc(ffi) - (x-xsiv(1,ffi))*tilts(ffi))*c1m3! mm -> m
      px_tp= ((px           )*tiltc(ffi) + (py           )*tilts(ffi))*c1m3*(one+dpsv(ffj))
      py_tp= ((py           )*tiltc(ffi) - (px           )*tilts(ffi))*c1m3*(one+dpsv(ffj))
      x=x_tp;   y=y_tp;   px=px_tp;   py=py_tp;   

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [2, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug
     
      ! Selection of the particle that are only in the radius (r = 0.08m)
      ! -------------------------------------------------------------------------------------------- !
      if (x*x+y*y<=r0_2) then
        ! Check particle enter the Quad
        ffInQuad(nlostp(ffj))=.true.

        !   - 
        ffdelta    = dpsv(ffj);
!        gam0     = gamma0;
        zb       = 0;
!        sigma_s  = 0;
!        betabeta0= ejfv(ffj)/ejv(ffj);
!        betabeta0= betabeta0*betabeta0*rvv(ffj);
        LoutQ    = (ffTable(iFile)%Lgth-ffTable(iFile)%Lin)
!        Ldpsv1   = LoutQ*oidpsv(ffj);
!        Ldpsv2   = oidpsv(ffj)*oidpsv(ffj);

!  <<<<<<< IN
  	!   - Initial repositionning (AntiDrift)
        x_tp=x-LoutQ*oidpsv(ffj)*px;
        y_tp=y-LoutQ*oidpsv(ffj)*py;
        x=x_tp;   y=y_tp;
!  <<<<<<< IN

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [3, ",ffj,"]  ->  ",x, px, y, py,'  > oidpsv=',oidpsv(ffj)
!endif
! <<<<<<<<<<<<<<<<< Debug

  	!   - Compute Fringe Field using asymplectic Map (Lie2)
        call ffTable(iFile)%Lie2(x,px,y,py,zb,oidpsv(ffj))

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [4, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug

!  <<<<<<< IN
        !   - Check AQ matrix are computed for a ffdelta in [Tdpsv(1),Tdpsv(nbDlt)]
        if (abs(ffdelta)>ffTable(iFile)%Tdpsv(ffTable(iFile)%nbDlt)+c1m6) then
          call ffield_genAntiQuad()
        end if

        !   - Find the AQ matrix for the right dpsv
        itDlt=0
        do k=1,ffTable(iFile)%nbDlt
          if (abs(dpsv(ffj)-ffTable(iFile)%Tdpsv(k))<c1m6) then
            itDlt=k
          endif
        enddo
        if (itDlt==0) itDlt=1

        !   - Final repositionning (AntiQuad)
        x_tp =ffTable(iFile)%TAQx(1,1,itDlt)*x + ffTable(iFile)%TAQx(1,2,itDlt)*px
        px_tp=ffTable(iFile)%TAQx(2,1,itDlt)*x + ffTable(iFile)%TAQx(2,2,itDlt)*px
        y_tp =ffTable(iFile)%TAQy(1,1,itDlt)*y + ffTable(iFile)%TAQy(1,2,itDlt)*py
        py_tp=ffTable(iFile)%TAQy(2,1,itDlt)*y + ffTable(iFile)%TAQy(2,2,itDlt)*py
        x=x_tp;   px=px_tp;   y=y_tp;   py=py_tp;

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [5, ",ffj,"]  ->  ",x, px, y, py
!write(lout,*)"FFIELD> [5, ",ffj,"]  ->  ",itDlt," / ",ffTable(iFile)%nbDlt
!endif
! <<<<<<<<<<<<<<<<< Debug

        x_tp=x+LoutQ*oidpsv(ffj)*px;
        y_tp=y+LoutQ*oidpsv(ffj)*py;
        x=x_tp;   y=y_tp;
!  <<<<<<< IN

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [6, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug

        ! Change to SixTrack referenciale
        ! ------------------------------------------------------------------------------------------ !
        px_tp= (tiltc(ffi)*px - tilts(ffi)*py)*c1e3*oidpsv(ffj)
        py_tp= (tilts(ffi)*px + tiltc(ffi)*py)*c1e3*oidpsv(ffj)
        x_tp = (tiltc(ffi)*x  - tilts(ffi)*y )*c1e3 + xsiv(1,ffi)   ! m -> mm
        y_tp = (tilts(ffi)*x  + tiltc(ffi)*y )*c1e3 + zsiv(1,ffi)   ! m -> mm
        x=x_tp;   px=px_tp;   y=y_tp;   py=py_tp;

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [7, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug

!  <<<<<<< IN
        ! Return to beginning of the Quad
        ! ------------------------------------------------------------------------------------------ !
        !             * HE
        x_tp = x + strack(ffi+2)*half*px
        y_tp = y + strack(ffi+2)*half*py

        !             * Riccardo's Quad head
!        x_tp = x + strack(ffi-1)*px
!        y_tp = y + strack(ffi-1)*py

        x=x_tp;   y=y_tp;
!  <<<<<<< IN

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [8, ",ffj,"]  ->  ",x, px, y, py
!write(lout,*)"FFIELD> Debug ****************************"
!endif
! <<<<<<<<<<<<<<<<< Debug

        ! Save data
        ! ------------------------------------------------------------------------------------------ !
        xv(1,ffj) = x;  xv(2,ffj) = y;  yv(1,ffj) = px;  yv(2,ffj) = py;
      else
        ffInQuad(nlostp(ffj))=.false.
      end if
    end do

! <<<<<<<<<<<<<<<<< Debug
!stop 1
! <<<<<<<<<<<<<<<<< Debug!

  end subroutine ffield_enterQuad





  ! ================================================================================================ !
  !  Lie 2 Tracking if particles left the Quadupole
  !  B. Dalena, T. Pugnat and A. Simona from CEA
  !  V.K. Berglyd Olsen, BE-ABP-HSS
  !  Last modified: 2019-01-21
  ! ================================================================================================ !
  subroutine ffield_exitQuad(ffi)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use parpro,              only : nblo
    use mod_common,          only : napx, ic, tiltc, tilts
    use mod_commont,         only : strack
    use mod_commonmn,        only : xv, yv, oidpsv, dpsv, rvv, ejv, zsiv, xsiv, nlostp
    use numerical_constants, only : half, one, c1e3, c1m3, c1m6

    implicit none

    ! Interface variables
    ! ---------------------------------------------------------------------------------------------- !
    integer, intent(in):: ffi

    ! Subroutine variables
    ! ---------------------------------------------------------------------------------------------- !
    integer          :: iFile
    integer          :: ffj,  k                    ! iterator
    integer          :: itDlt                      ! 
    real(kind=fPrec) :: x, px, y, py               ! Transverse canonical parameter in new referencial
    real(kind=fPrec) :: x_tp, px_tp, y_tp, py_tp   ! Temp transverse canonical parameter

    real(kind=fPrec) :: ffdelta!, gam0, betabeta0     ! 
    real(kind=fPrec) :: zb!, sigma_s                ! 
    real(kind=fPrec) :: LoutQ!, Ldpsv1, Ldpsv2      ! 

    
!  <<<<<<< OUT
    iFile=ffQ2File(ffindex(ic(ffi)-nblo),2)
!  <<<<<<< OUT

! <<<<<<<<<<<<<<<<< Debug
!write(lout,*)"FFIELD> Debug ************ EX ************"
!write(lout,*)"FFIELD> Debug [etp ffj]  ->  x, px, y, py"
! <<<<<<<<<<<<<<<<< Debug
    do ffj=1,napx
      ! Save data
      ! -------------------------------------------------------------------------------------------- !
      x = xv(1,ffj);  y = xv(2,ffj);  px= yv(1,ffj);  py= yv(2,ffj);

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [0, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug

!  <<<<<<< OUT
      ! Return to beginning of the Quad
      ! -------------------------------------------------------------------------------------------- !
      !             * HE
      x_tp = x + strack(ffi+2)*half*px
      y_tp = y + strack(ffi+2)*half*py


      !             * Riccardo's Quad head
!      x_tp = x + strack(ffi-1)*px
!      y_tp = y + strack(ffi-1)*py

      x=x_tp;   y=y_tp;
!  <<<<<<< OUT

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [1, ",ffj,"]  ->  ",x, px, y, py,"   -> strack(ffi+2)*half*px",strack(ffi+2)*half*px
!endif
! <<<<<<<<<<<<<<<<< Debug

      ! Rotation error for the quad
      ! -------------------------------------------------------------------------------------------- !
      x_tp = ((x-xsiv(1,ffi))*tiltc(ffi) + (y-zsiv(1,ffi))*tilts(ffi))*c1m3! mm -> m
      y_tp = ((y-zsiv(1,ffi))*tiltc(ffi) - (x-xsiv(1,ffi))*tilts(ffi))*c1m3! mm -> m
      px_tp= ((px           )*tiltc(ffi) + (py           )*tilts(ffi))*c1m3*(one+dpsv(ffj))
      py_tp= ((py           )*tiltc(ffi) - (px           )*tilts(ffi))*c1m3*(one+dpsv(ffj))
      x=x_tp;   px=px_tp;   y=y_tp;   py=py_tp;

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [2, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug

      ! Selection of the particle that are only in the radius (r = 0.08m)
      ! -------------------------------------------------------------------------------------------- !
      if (ffInQuad(nlostp(ffj))) then
        ! Check particle enter the Quad
        ffInQuad(nlostp(ffj))=.false.

        !   - 
        ffdelta    = dpsv(ffj);
!        gam0     = gamma0;
        zb       = 0;
!        sigma_s  = 0;
!        betabeta0= ejfv(ffj)/ejv(ffj);
!        betabeta0= betabeta0*betabeta0*rvv(ffj);
        LoutQ    = (ffTable(iFile)%Lgth-ffTable(iFile)%Lin)
!        Ldpsv1   = LoutQ*oidpsv(ffj);
!        Ldpsv2   = oidpsv(ffj)*oidpsv(ffj);


!  <<<<<<< OUT
        !   - Check AQ matrix are computed for a ffdelta in [Tdpsv(1),Tdpsv(nbDlt)]
        if (abs(ffdelta)>ffTable(iFile)%Tdpsv(ffTable(iFile)%nbDlt)+c1m6) then
          call ffield_genAntiQuad()
        end if


        !   - Find the AQ matrix for the right dpsv
        itDlt=0
        do k=1,ffTable(iFile)%nbDlt
          if (abs(dpsv(ffj)-ffTable(iFile)%Tdpsv(k))<c1m6) then
            itDlt=k
          endif
        enddo
        if (itDlt==0) itDlt=1


        !   - Final repositionning (AntiQuad)
        x_tp=x+LoutQ*oidpsv(ffj)*px;
        y_tp=y+LoutQ*oidpsv(ffj)*py;
        x=x_tp;   y=y_tp;

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [3, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug


        x_tp =ffTable(iFile)%TAQx(1,1,itDlt)*x + ffTable(iFile)%TAQx(1,2,itDlt)*px
        px_tp=ffTable(iFile)%TAQx(2,1,itDlt)*x + ffTable(iFile)%TAQx(2,2,itDlt)*px
        y_tp =ffTable(iFile)%TAQy(1,1,itDlt)*y + ffTable(iFile)%TAQy(1,2,itDlt)*py
        py_tp=ffTable(iFile)%TAQy(2,1,itDlt)*y + ffTable(iFile)%TAQy(2,2,itDlt)*py
        x=x_tp;   px=px_tp;   y=y_tp;   py=py_tp;
!  <<<<<<< OUT

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [4, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug


  	!   - Compute Fringe Field using asymplectic Map (Lie2)
        call ffTable(iFile)%Lie2(x,px,y,py,zb,oidpsv(ffj))

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [5, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug


!  <<<<<<< OUT
  	!   - Initial repositionning (AntiDrift)
        x_tp=x-LoutQ*oidpsv(ffj)*px;
        y_tp=y-LoutQ*oidpsv(ffj)*py;
        x=x_tp;   y=y_tp;
!  <<<<<<< OUT

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [6, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug


        ! Change to SixTrack referenciale
        ! ------------------------------------------------------------------------------------------ !
        px_tp= (tiltc(ffi)*px - tilts(ffi)*py)*c1e3*oidpsv(ffj)
        py_tp= (tilts(ffi)*px + tiltc(ffi)*py)*c1e3*oidpsv(ffj)
        x_tp = (tiltc(ffi)*x  - tilts(ffi)*y )*c1e3 + xsiv(1,ffi)   ! m -> mm
        y_tp = (tilts(ffi)*x  + tiltc(ffi)*y )*c1e3 + zsiv(1,ffi)   ! m -> mm
        x=x_tp;   px=px_tp;   y=y_tp;   py=py_tp;

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [7, ",ffj,"]  ->  ",x, px, y, py
!endif
! <<<<<<<<<<<<<<<<< Debug


!  <<<<<<< OUT
        ! Return to beginning of the Quad
        ! ------------------------------------------------------------------------------------------ !
        !             * HE
        x_tp = x - strack(ffi+2)*half*px
        y_tp = y - strack(ffi+2)*half*py


      !             * Riccardo's Quad head
!        x_tp = x - strack(ffi-1)*px
!        y_tp = y - strack(ffi-1)*py

        x=x_tp;   y=y_tp;
!  <<<<<<< OUT

! <<<<<<<<<<<<<<<<< Debug
!if (ffj==1) then
!write(lout,*)"FFIELD> [8, ",ffj,"]  ->  ",x, px, y, py
!write(lout,*)"FFIELD> Debug ****************************"
!endif
! <<<<<<<<<<<<<<<<< Debug


        ! Save data
        ! ------------------------------------------------------------------------------------------ !
        xv(1,ffj) = x;  xv(2,ffj) = y;  yv(1,ffj) = px;  yv(2,ffj) = py;
      end if
    end do

  end subroutine ffield_exitQuad

end module mod_ffield



