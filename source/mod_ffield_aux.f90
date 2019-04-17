! ================================================================================================ !
!  Fringe Field Tracking Module
! ~~~~~~~~~~~~~~~~~~~~~~~~
!  B. Dalena, T. Pugnat and A. Simona from CEA
!  Last Modified: 2019-02-01
!
!  Usage
! ~~~~~~~
!  Declaration:   type(ffTable_n_Track) aFile
!  Set:           aFile%set("ffFNames",Kin,Lgth) => Set File's name and info
!  Load:          aFile%load()                   => Load File
!  Free:          aFile%free()                   => Free memory
!
!  Note: None yet
! ================================================================================================ !
module ffTable_n_Tracks
  ! ------------------------------------------------------------------------------------------------ !
  ! Mod from SixTrack
  ! ------------------------------------------------------------------------------------------------ !
  use floatPrecision
  use, intrinsic :: iso_fortran_env, only : int32

  implicit none

  type, public :: ffTable_n_Track
    character(len=:), allocatable, public :: ffFNames      ! Filename of the Vec. Pot. coefficient
    integer(kind=2),               public :: chk_Status    ! Check statut of the file (0=empty, 1=ready,2=loaded)

    integer,                       public :: n,max_i       ! Max exposant for x
    integer,                       public :: m,max_j       ! Max exposant for y
    integer,                       public :: s             ! Number of point in z

    real(kind=fPrec),              public :: dz            ! Step size in z
    real(kind=fPrec),              public :: norm          ! normalization
    real(kind=fPrec),              public :: Lgth          ! Total file length
    real(kind=fPrec),              public :: Lin           ! Length in Quad
    real(kind=fPrec),              public :: r0_2          ! Maximum radius in quad

    ! Table for the vector potential Ax
    integer,                       public :: lx            ! Number of coef per step for Ax
    integer,          allocatable, public :: ij_TAx(:,:,:) ! (1:2,1:lx,1:s) Table of indices ij
    real(kind=fPrec), allocatable, public :: TAx(:,:)      ! (    1:lx,1:s) Table of coefficients

    ! Table for the vector potential Ay
    integer,                       public :: ly            ! Number of coef per step for Ay
    integer,          allocatable, public :: ij_TAy(:,:,:) ! (1:2,1:ly,1:s) Table of indices ij
    real(kind=fPrec), allocatable, public :: TAy(:,:)      ! (    1:ly,1:s) Table of coefficients

    ! Table for the vector potential Az
    integer,                       public :: lz            ! Number of coef per step for Az
    integer,          allocatable, public :: ij_TAz(:,:,:) ! (1:2,1:lz,1:s) Table of indices ij
    real(kind=fPrec), allocatable, public :: TAz(:,:)      ! (    1:lz,1:s) Table of coefficients

    ! Matrix for the antiquadripole
    integer :: nbDlt
    real(kind=fPrec), allocatable, public :: Tdpsv(:)      ! (1:nbDlt) dpsv value for AQ-case
!    real(kind=fPrec), allocatable, public :: AQx(:,:,:)    ! (1:2,1:2,1:npart) AntiQuad matrix
!    real(kind=fPrec), allocatable, public :: AQy(:,:,:)    ! (1:2,1:2,1:npart) AntiQuad matrix
    real(kind=fPrec), allocatable, public :: TAQx(:,:,:)   ! (1:2,1:2,1:nbDlt) Table of AntiQuad matrix
    real(kind=fPrec), allocatable, public :: TAQy(:,:,:)   ! (1:2,1:2,1:nbDlt) Table of AntiQuad matrix
  contains

    procedure, public, pass(this)  :: set    => Tset
    procedure, public, pass(this)  :: load   => Tload
    procedure, public, pass(this)  :: free   => Tfree
    procedure, public, pass(this)  :: AQgen  => TAQgen
    procedure, public, pass(this)  :: Lie2   => TLie2

    procedure, private, pass(this)  :: Tset
    procedure, private, pass(this)  :: Tload
    procedure, private, pass(this)  :: Tfree
    procedure, private, pass(this)  :: TAQgen
    procedure, private, pass(this)  :: TLie2

    procedure, private, pass(this)  :: HornerDX_Az
    procedure, private, pass(this)  :: HornerDY_Az
    procedure, private, pass(this)  :: Horner2D_Ax
    procedure, private, pass(this)  :: Horner2D_Ay
    procedure, private, pass(this)  :: HornerDYIntX_Ax
    procedure, private, pass(this)  :: HornerDXIntY_Ay

    procedure, private, nopass      :: ReadExpMax
    procedure, private, nopass      :: ReadVectPotCoeff


!    procedure, private, pass(left)  :: assignStrStr
!    procedure, private, pass(left)  :: assignStrChr
!    procedure, private, pass(right) :: assignChrStr

  end type ffTable_n_Track

  interface ffTable_n_Track
    module procedure constructT
  end interface ffTable_n_Track

!  interface assignment(=)
!    module procedure assignStrStr
!    module procedure assignStrChr
!    module procedure assignChrStr
!  end interface

!  interface operator(==)
!    module procedure compStrStr
!    module procedure compStrChr
!    module procedure compChrStr
!  end interface

!  interface operator(/=)
!    module procedure compNStrStr
!    module procedure compNStrChr
!    module procedure compNChrStr
!  end interface

  private :: Tset
  private :: Tload
  private :: Tfree
  private :: TLie2
  private :: ReadExpMax
  private :: ReadVectPotCoeff
  private :: TAQgen
  private :: HornerDX_Az
  private :: HornerDY_Az
  private :: Horner2D_Ax
  private :: Horner2D_Ay
  private :: HornerDYIntX_Ax
  private :: HornerDXIntY_Ay

contains


  ! ================================================================ !
  !  Constructor
  ! ================================================================ !
  type(ffTable_n_Track) function constructT()
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use parpro,              only : mPathName
    use numerical_constants, only : zero

    implicit none

    if (allocated(constructT%ffFNames)) deallocate(constructT%ffFNames)
    allocate(character(len=mPathName) :: constructT%ffFNames)
    constructT%ffFNames=" "
    constructT%chk_Status=0
    constructT%n=0
    constructT%m=0
    constructT%s=0
    constructT%max_i=0
    constructT%max_j=0

    constructT%dz  =zero
    constructT%norm=zero
    constructT%Lgth=zero
    constructT%Lin =zero
    constructT%r0_2=6.4e-3_fPrec

    constructT%lx=0
    constructT%ly=0
    constructT%lz=0

    constructT%nbDlt=0
  end function constructT


  ! ================================================================ !
  !  Set
  ! ================================================================ !
  subroutine Tset(this, nFile, LQin, Length, r0)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, c1m12

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(inout) :: this
    character(len=*),       intent(in)    :: nFile
    real(kind=fPrec),       intent(in)    :: LQin
    real(kind=fPrec),       intent(in)    :: Length
    real(kind=fPrec),       intent(in)    :: r0     ! Physical aperture

    this%r0_2=6.4e-3_fPrec    ! Default value

    this%ffFNames   = TRIM(ADJUSTL(nFile))
    this%Lin        = LQin
    this%Lgth       = Length
    this%chk_Status = 1
    if(r0>c1m12) this%r0_2= r0*r0
  end subroutine Tset


  ! ================================================================ !
  !  Free
  ! ================================================================ !
  subroutine Tfree(this)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use mod_alloc,           only : dealloc
    use numerical_constants, only : zero

    implicit none


    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(inout) :: this

    if (this%chk_Status/=0) then
      if(allocated(this%ffFNames)) deallocate(this%ffFNames)

      if(allocated(this%ij_TAx))   call dealloc(this%ij_TAx, 'this%ij_TAx')
      if(allocated(this%TAx))      call dealloc(this%TAx,    'this%TAx')
      this%lx=0

      if(allocated(this%ij_TAy))   call dealloc(this%ij_TAy, 'this%ij_TAy')
      if(allocated(this%TAy))      call dealloc(this%TAy,    'this%TAy')
      this%ly=0

      if(allocated(this%ij_TAz))   call dealloc(this%ij_TAz, 'this%ij_TAz')
      if(allocated(this%TAz))      call dealloc(this%TAz,    'this%TAz')
      this%lz=0

      if(allocated(this%Tdpsv))    call dealloc(this%Tdpsv,  'this%Tdpsv')

      if(allocated(this%TAQx))     call dealloc(this%TAQx,   'this%TAQx')
      if(allocated(this%TAQy))     call dealloc(this%TAQy,   'this%TAQy')

      this%n=0
      this%m=0
      this%s=0

      this%dz  =zero
      this%norm=zero
      this%Lgth=zero
      this%Lin =zero
      this%chk_Status=0
    end if
  end subroutine Tfree


  ! ================================================================ !
  !  Load
  ! ================================================================ !
  subroutine Tload(this,norm,iErr)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use mod_units, only : f_requestUnit
    use crcoall,   only : lout, lerr

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(inout) :: this
    logical,                intent(inout) :: iErr
    real(kind=fPrec),       intent(in)    :: norm      ! Normalization

    integer :: lun


    if (this%chk_Status<1) then
      write(lerr,"(a,i1)") "FFIELD> ERROR while loading file '"//trim(this%ffFNames)//"' -> status ",this%chk_Status
      iErr = .true.
      return

    else if (this%chk_Status == 1) then
      ! Request to find unit
      ! -------------------------------------------------------------------------------------------- !
      call f_requestUnit(this%ffFNames,lun)

      ! Detect the mamixum exposant for x and y, and the number of point in z
      ! -------------------------------------------------------------------------------------------- !
      call ReadExpMax(this,lun,iErr)
      if (iErr) return

      ! Generate the vecto potential table
      ! -------------------------------------------------------------------------------------------- !
      call ReadVectPotCoeff(this,lun,norm,iErr)
      if (iErr) return

      this%s=this%s-1

      this%chk_Status=2

    end if
  end subroutine Tload


  subroutine ReadExpMax(this,lun,iErr)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, c1e12, c1m12
    use crcoall,             only : lout, lerr
    use parpro,              only : mInputLn
    use mod_units,           only : f_open, f_close
    use string_tools,        only : chr_split, chr_cast

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(inout) :: this
    integer,                intent(inout) :: lun  ! Unit for file to be open
    logical,                intent(inout) :: iErr ! Error return

    character(len=mInputLn) :: inLine             ! Line of the file
    integer :: nSplit                             ! Nb element in line
    integer :: istat                              ! Check. accecibility of the file
    integer :: line                               ! Nb of line in file
    integer :: n,m,s                              ! Exposant max for x and y, Nb of point in z
    integer :: expx,expy,expz                     ! Expo. for x, y and z
    real(kind=fPrec) :: st, sm1                   ! Parametter for the detection of new step in z
    character(len=:), allocatable :: lnSplit(:)   ! (1:nSplit) Splited line

    ! Open file
    ! ---------------------------------------------------------------------------------------------- !
    call f_open(unit=lun,file=trim(this%ffFNames),formatted=.true.,mode='r',err=iErr,status="old")
    if (iErr) then
!      write(lerr,"(a)")"FFIELD> ERROR in ReadExpMax(): Error opening file '"//trim(this%ffFNames)//"'"
      return
    end if

    ! Initialization
    ! ---------------------------------------------------------------------------------------------- !
    istat=0
    n=0
    m=0
    s=1
    line=1
    sm1=c1e12

    ! Read file
    ! ---------------------------------------------------------------------------------------------- !
    do while(istat==0)
      st=zero; expx=0; expy=0;! expz=0;

      !     - Read line
      read(lun,"(a)",iostat=istat) inLine

      !     - Split line
      call chr_split(inLine, lnSplit, nSplit, iErr)

      !     - Check error in line
      if(iErr) then
        write(lerr,"(a)") "FFIELD> ERROR in ReadExpMax(): Failed to read '"//trim(this%ffFNames)//"'"
        return
      end if
      if(nSplit < 7)then
        write(lerr,"(a)") "FFIELD> ERROR in ReadExpMax(): Wrong number of element in '"//trim(this%ffFNames)//"'"
        write(lerr,"(a)") "FFIELD> ERROR in ReadExpMax(): Line '"//trim(inLine)//"'"
        iErr = .true.
        return
      end if

      !     - Save data
      call chr_cast(lnSplit(1),st,iErr)
      call chr_cast(lnSplit(2),expx,iErr)
      call chr_cast(lnSplit(3),expy,iErr)
!      call chr_cast(lnSplit(4),expz,iErr)
      if (iErr) return

      if (st>sm1+c1m12) then                  ! Detect new step in z
         s = s + 1                                  ! Compte step in z
      endif
      if (expx>n) n=expx                            ! Max. expo. in x
      if (expy>m) m=expy                            ! Max. expo. in y
      sm1=st
      line=line+1                                   ! File size
    end do

    ! Save data
    ! ---------------------------------------------------------------------------------------------- !
    this%n=n
    this%m=m
    this%s=s

    ! Close file
    ! ---------------------------------------------------------------------------------------------- !
    call f_close(lun)

  end subroutine ReadExpMax


  subroutine ReadVectPotCoeff(this,lun,norm,iErr)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, c1e12, c1m12
    use crcoall,             only : lout, lerr
    use parpro,              only : mInputLn
    use mod_alloc,           only : alloc, dealloc
    use mod_units,           only : f_open, f_close
    use string_tools,        only : chr_split, chr_cast

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(inout) :: this
    integer,                intent(inout) :: lun       ! Unit for file to be open
    logical,                intent(inout) :: iErr      ! Error return
    real(kind=fPrec),       intent(in)    :: norm      ! Normalization

    ! Subroutine variables
    ! ---------------------------------------------------------------------------------------------- !
    character(len=mInputLn) :: inLine              ! Line of the file
    integer :: nSplit                              ! Nb element in line
    logical :: CoefSave                            ! Check if coef. already exist
    integer :: istat                               ! Check accecibility of the file
    integer :: ind                                 ! Nb. max of different coefficient
    integer :: l,k                                 ! Iterator
    integer :: linex,liney,linez                   ! Nb. max of coeff/step for Ax,Ay and Az
    integer :: tlinex,tliney,tlinez                ! Temp. nb. max of coeff/step for Ax,Ay and Az
    integer :: sline                               ! Nb of line in file
    integer :: expx, expy, expz                    ! Expo. for x, y and z
    real(kind=fPrec) :: ax, ay, az                 ! Coeff. for x, y and z
    real(kind=fPrec) :: st, sm1                    ! Parametter for the detection of new step in z
    real(kind=fPrec) :: zin                        ! Position of the first step
    real(kind=fPrec) :: dz                         ! Step size in z
    character(len=:)   , allocatable :: lnSplit(:) ! (1:nSplit)       Splited line
    integer, allocatable :: tpij_Ax(:,:,:)         ! (1:2,1:ind,0:s)) Temp. table of indices ij for Ax
    integer, allocatable :: tpij_Ay(:,:,:)         ! (1:2,1:ind,0:s)) Temp. table of indices ij for Ay
    integer, allocatable :: tpij_Az(:,:,:)         ! (1:2,1:ind,0:s)) Temp. table of indices ij for Az

    real(kind=fPrec), allocatable :: tpTAx(:,:)    ! (    1:ind,1:s)) Temp. table for coeff. of Ax
    real(kind=fPrec), allocatable :: tpTAy(:,:)    ! (    1:ind,1:s)) Temp. table for coeff. of Ay
    real(kind=fPrec), allocatable :: tpTAz(:,:)    ! (    1:ind,1:s)) Temp. table for coeff. of Az

    ! Allocation of memory for tables
    ! ---------------------------------------------------------------------------------------------- !
    ind= this%n*this%m/2+1
    tlinex=1; tliney=1; tlinez=1
    linex=1;  liney=1;  linez=1
    call alloc(tpTAx,      ind, this%s, zero, 'tpTAx')   ! (    1:ind,1:s))
    call alloc(tpTAy,      ind, this%s, zero, 'tpTAy')   ! (    1:ind,1:s))
    call alloc(tpTAz,      ind, this%s, zero, 'tpTAz')   ! (    1:ind,1:s))
    call alloc(tpij_Ax, 2, ind, this%s, 0,    'tpij_Ax') ! (1:2,1:ind,0:s))
    call alloc(tpij_Ay, 2, ind, this%s, 0,    'tpij_Ay') ! (1:2,1:ind,0:s))
    call alloc(tpij_Az, 2, ind, this%s, 0,    'tpij_Az') ! (1:2,1:ind,0:s))

    ! Open file
    ! ---------------------------------------------------------------------------------------------- !
    call f_open(unit=lun,file=trim(this%ffFNames),formatted=.true.,mode='r',err=iErr,status="old")
    if (iErr) then
!      write(lerr,"(a)")"FFIELD> ERROR in ReadVectPotCoeff(): Error opening file '"//trim(this%ffFNames)//"'"
      return
    end if

    ! Initialize parameters
    ! ---------------------------------------------------------------------------------------------- !
    sm1=c1e12
    istat=0
    sline=1
    dz=zero

    ! Read coef file
    ! ---------------------------------------------------------------------------------------------- !
    do while(istat==0)
      st=zero; expx=0;  expy=0;! expz=0;
      ax=zero; ay=zero; az=zero;
      !     - Read line
      read(lun,"(a)",iostat=istat) inLine

      !     - Split line
      call chr_split(inLine, lnSplit, nSplit, iErr)

      !     - Check error in line
      if(iErr) then
        write(lerr,"(a)") "FFIELD> ERROR in ReadVectPotCoeff(): Fail to read '"//trim(this%ffFNames)//"'"
        return
      end if
      if(nSplit < 7)then
        write(lerr,"(a)") "FFIELD> ERROR in ReadVectPotCoeff(): Wrong number of element in '"//trim(this%ffFNames)//"'"
        write(lerr,"(a)") "FFIELD> ERROR in ReadVectPotCoeff(): Line '"//trim(inLine)//"'"
        return
      end if

      !     - Save data
      call chr_cast(lnSplit(1),st,  iErr)
      call chr_cast(lnSplit(2),expx,iErr)
      call chr_cast(lnSplit(3),expy,iErr)
!      call chr_cast(lnSplit(4),expz,iErr)
      call chr_cast(lnSplit(5),ax,  iErr)
      call chr_cast(lnSplit(6),ay,  iErr)
      call chr_cast(lnSplit(7),az,  iErr)
      if (iErr) return


      !     - Detect step in z
      if (st>sm1+c1m12) then
        dz=dz+st-sm1
        sline = sline + 1
        tlinex=1; tliney=1; tlinez=1
      endif
      if (sline == 1) zin = st


      !     - Detect coef non trivial for Ax
      if (ax/=zero) then
        CoefSave=.False.

        if (tlinex/=1) then
          do l=1,tlinex                ! Check coef already saved
            if ((tpij_Ax(1,l,sline)==expx).AND.(tpij_Ax(2,l,sline)==expy)) then
              tpTAx(l,sline)=tpTAx(l,sline) + ax
              CoefSave=.True.
            endif
          enddo
        endif

        if (CoefSave.eqv..False.) then ! Create new coef
          tpij_Ax(1,tlinex,sline)=expx
          tpij_Ax(2,tlinex,sline)=expy
          tpTAx(tlinex,sline)=tpTAx(tlinex,sline) + ax
          tlinex=tlinex+1
          CoefSave=.True.
        endif
      endif

      !     - Detect coef non trivial for Ay
      if (ay/=zero) then
        CoefSave=.False.

        if (tliney/=1) then
          do l=1,tliney                ! Check coef already saved
            if ((tpij_Ay(1,l,sline)==expx).AND.(tpij_Ay(2,l,sline)==expy)) then
              tpTAy(l,sline)=tpTAy(l,sline) + ay
              CoefSave=.True.
            endif
          enddo
        endif

        if (CoefSave.eqv..False.) then ! Create new coef
          tpij_Ay(1,tliney,sline)=expx
          tpij_Ay(2,tliney,sline)=expy
          tpTAy(tliney,sline)=tpTAy(tliney,sline) + ay
          tliney=tliney+1
          CoefSave=.True.
        endif
      endif

      !     - Detect coef non trivial for Az
      if (az/=zero) then
        CoefSave=.False.

        if (tlinez/=1) then
          do l=1,tlinez                ! Check coef already saved
            if ((tpij_Az(1,l,sline)==expx).AND.(tpij_Az(2,l,sline)==expy)) then
              tpTAz(l,sline)=tpTAz(l,sline) + az
              CoefSave=.True.
            endif
          enddo
        endif

        if (CoefSave.eqv..False.) then ! Create new coef
          tpij_Az(1,tlinez,sline)=expx
          tpij_Az(2,tlinez,sline)=expy
          tpTAz(tlinez,sline)=tpTAz(tlinez,sline) + az
          tlinez=tlinez+1
          CoefSave=.True.
        endif
      endif

      !     - Update tables size
      if (linex<tlinex-1) then
        linex=tlinex-1
      endif
      if (liney<tliney-1) then
        liney=tliney-1
      endif
      if (linez<tlinez-1) then
        linez=tlinez-1
      endif

      if (istat==0) then
        sm1=st
      endif
      if (sline>this%s) then           ! Detect error in z
        write(lerr,"(a)")"FFIELD> ERROR in ReadExpMax(): Wrong nb. of step for '"//this%ffFNames//"'"
        iErr=.False.
        return
      endif
    enddo
    this%Lgth = sm1 - zin

    ! Normalize coefficients
    ! -------------------------------------------------------------------------
    this%lx=linex; this%ly=liney; this%lz=linez
    do k=1,this%s
      do l=1,linex
        tpTAx(l,k)=tpTAx(l,k)*norm
      enddo

      do l=1,liney
        tpTAy(l,k)=tpTAy(l,k)*norm
      enddo

      do l=1,linez
        tpTAz(l,k)=tpTAz(l,k)*norm
      enddo
    enddo
    this%dz=dz/dble(sline)

    ! Reduce table size
    ! -------------------------------------------------------------------------
    call alloc(this%TAx,       linex, this%s, zero, 'this%TAx')    ! (line,z)
    call alloc(this%TAy,       liney, this%s, zero, 'this%TAy')    ! (line,z)
    call alloc(this%TAz,       linez, this%s, zero, 'this%TAz')    ! (line,z)
    call alloc(this%ij_TAx, 2, linex, this%s, 0,    'this%ij_TAx') ! (1 pour i et 2 pour j,line,z)
    call alloc(this%ij_TAy, 2, liney, this%s, 0,    'this%ij_TAy') ! (1 pour i et 2 pour j,line,z)
    call alloc(this%ij_TAz, 2, linez, this%s, 0,    'this%ij_TAz') ! (1 pour i et 2 pour j,line,z)

    this%TAx(1:linex,1:this%s)     =tpTAx(1:linex,1:this%s)
    this%TAy(1:liney,1:this%s)     =tpTAy(1:liney,1:this%s)
    this%TAz(1:linez,1:this%s)     =tpTAz(1:linez,1:this%s)
    this%ij_TAx(:,1:linex,1:this%s)=tpij_Ax(:,1:linex,1:this%s)
    this%ij_TAy(:,1:liney,1:this%s)=tpij_Ay(:,1:liney,1:this%s)
    this%ij_TAz(:,1:linez,1:this%s)=tpij_Az(:,1:linez,1:this%s)

    call dealloc(tpTAx,'tpTAx')
    call dealloc(tpTAy,'tpTAy')
    call dealloc(tpTAz,'tpTAz')
    call dealloc(tpij_Ax,'tpij_Ax')
    call dealloc(tpij_Ay,'tpij_Ay')
    call dealloc(tpij_Az,'tpij_Az')

    ! Close file
    ! ---------------------------------------------------------------------------------------------- !
    call f_close(lun)

  end subroutine ReadVectPotCoeff


  ! ================================================================ !
  !  AQgen
  ! ================================================================ !
  subroutine TAQgen(this,delta,nbDlt)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, half, one, two
    use mod_alloc,           only : alloc, dealloc

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(inout) :: this
    real(kind=fPrec),       intent(in)    :: delta     ! Delta
    integer,                intent(in)    :: nbDlt     ! Nb of point in the interval [-delta,delta]


    ! Subroutine variables
    ! ---------------------------------------------------------------------------------------------- !
    integer :: i,l,k                                   ! iterator
    real(kind=fPrec) :: p0sp                           ! 1/(delta+1) = p0/p
    real(kind=fPrec) :: a1,b1,c1,d1                    ! Matrix coef. for Qx normal
    real(kind=fPrec) :: a2,b2,c2,d2                    ! Matrix coef. for Qy normal
    real(kind=fPrec) :: atp,btp,ctp,dtp,bcmad          ! Temp. value for matrix coef.
    real(kind=fPrec) :: tp,tp1,tp2,tp3
    real(kind=fPrec) :: C0x,C0y                        ! Vect. Pot. normal coef for Ax and Ay


    if (this%chk_Status==2) then

      ! Initialisation of size array
      ! -------------------------------------------------------------------------------------------- !
      this%nbDlt=nbDlt !max(1,ceiling(abs( two*(delta)*c1e6 )))

      ! Clean memory
      ! -------------------------------------------------------------------------------------------- !
      if (allocated(this%Tdpsv)) call dealloc(this%Tdpsv,'this%Tdpsv')
      if (allocated(this%TAQx))  call dealloc(this%TAQx, 'this%TAQx')
      if (allocated(this%TAQy))  call dealloc(this%TAQy, 'this%TAQy')
      call alloc(this%Tdpsv,       this%nbDlt, zero, 'this%Tdpsv')
      call alloc(this%TAQx,  2, 2, this%nbDlt, zero, 'this%TAQx')
      call alloc(this%TAQy,  2, 2, this%nbDlt, zero, 'this%TAQy')

      ! Generate the Matrix of the antiquad for delta /= 0
      ! -------------------------------------------------------------------------------------------- !
      if (this%nbDlt==1) then  ! * for delta == 0
        k=1
        this%Tdpsv(k)=delta
      else                     ! * for delta /= 0
        do k=1,this%nbDlt
          this%Tdpsv(k)=(two*(k-1)/(this%nbDlt-1)- one)*delta
        end do
      end if

      do k=1,this%nbDlt
        ! Initialisation of the parameter
        ! ------------------------------------------------------------------------------------------ !
        a1=one;   b1=zero;   c1=zero;   d1=one;
        a2=one;   b2=zero;   c2=zero;   d2=one;
        p0sp=one/(one+this%Tdpsv(k))

        do i=1,this%s !-1
          C0x=zero; C0y=zero

          do l=1,this%lz
            if ((this%ij_TAz(1,l,i)==2).and.(this%ij_TAz(2,l,i)==0)) then
              C0x=two*this%TAz(l,i)
            endif

            if ((this%ij_TAz(1,l,i)==0).and.(this%ij_TAz(2,l,i)==2)) then
              C0y=two*this%TAz(l,i)
            endif
          enddo


          ! Approx
!          atp= (one+p0sp*C0x*this%dz*this%dz*half)*a1 + ((one+p0sp*C0x*this%dz*this%dz*sixth)*p0sp*this%dz)*c1
!          btp= (one+p0sp*C0x*this%dz*this%dz*half)*b1 + ((one+p0sp*C0x*this%dz*this%dz*sixth)*p0sp*this%dz)*d1
!          ctp= (C0x*this%dz*(one+p0sp*C0x*this%dz*this%dz*sixth))*a1 + (one+p0sp*C0x*this%dz*this%dz*half)*c1
!          dtp= (C0x*this%dz*(one+p0sp*C0x*this%dz*this%dz*sixth))*b1 + (one+p0sp*C0x*this%dz*this%dz*half)*d1

          ! Thin DKD
!          atp= (one+p0sp*C0x*this%dz*this%dz*half)*a1 + ((one+p0sp*C0x*this%dz*this%dz*0.25_fPrec)*p0sp*this%dz)*c1
!          btp= (one+p0sp*C0x*this%dz*this%dz*half)*b1 + ((one+p0sp*C0x*this%dz*this%dz*0.25_fPrec)*p0sp*this%dz)*d1
!          ctp= (C0x*this%dz)*a1 + (one+p0sp*C0x*this%dz*this%dz*half)*c1
!          dtp= (C0x*this%dz)*b1 + (one+p0sp*C0x*this%dz*this%dz*half)*d1
!          a1=atp; b1=btp; c1=ctp; d1=dtp;

          ! Thin KDK
!          atp= (one+((p0sp*C0x)*(this%dz*this%dz))*half)*a1 + (p0sp*this%dz)*c1
!          btp= (one+((p0sp*C0x)*(this%dz*this%dz))*half)*b1 + (p0sp*this%dz)*d1
!          ctp= (one+((p0sp*C0x)*(this%dz*this%dz))*0.25_fPrec)*(C0x*this%dz)*a1 + (one+((p0sp*C0x)*(this%dz*this%dz)))*c1
!          dtp= (one+((p0sp*C0x)*(this%dz*this%dz))*0.25_fPrec)*(C0x*this%dz)*b1 + (one+((p0sp*C0x)*(this%dz*this%dz)))*d1
!          a1=atp; b1=btp; c1=ctp; d1=dtp;

          tp = (p0sp*C0x)*(this%dz*this%dz)
          tp1= one+tp*half
          tp2= p0sp*this%dz
          tp3=(one+tp*0.25_fPrec)*(C0x*this%dz)
          atp= tp1*a1 + tp2*c1; btp= tp1*b1 + tp2*d1; ctp= tp3*a1 + tp1*c1; dtp= tp3*b1 + tp1*d1
          a1=atp; b1=btp; c1=ctp; d1=dtp;


          ! Approx
!          atp= (one+p0sp*C0y*this%dz*this%dz*half)*a2 + ((one+p0sp*C0y*this%dz*this%dz*sixth)*p0sp*this%dz)*c2
!          btp= (one+p0sp*C0y*this%dz*this%dz*half)*b2 + ((one+p0sp*C0y*this%dz*this%dz*sixth)*p0sp*this%dz)*d2
!          ctp= (C0y*this%dz*(one+p0sp*C0y*this%dz*this%dz*sixth))*a2 + (one+p0sp*C0y*this%dz*this%dz*half)*c2
!           dtp= (C0y*this%dz*(one+p0sp*C0y*this%dz*this%dz*sixth))*b2 + (one+p0sp*C0y*this%dz*this%dz*half)*d2

          ! Thin DKD
!          atp= (one+p0sp*C0y*this%dz*this%dz*half)*a2 + ((one+p0sp*C0y*this%dz*this%dz*0.25_fPrec)*p0sp*this%dz)*c2
!          btp= (one+p0sp*C0y*this%dz*this%dz*half)*b2 + ((one+p0sp*C0y*this%dz*this%dz*0.25_fPrec)*p0sp*this%dz)*d2
!          ctp= (C0y*this%dz)*a2 + (one+p0sp*C0y*this%dz*this%dz*half)*c2
!          dtp= (C0y*this%dz)*b2 + (one+p0sp*C0y*this%dz*this%dz*half)*d2
!          a2=atp; b2=btp; c2=ctp; d2=dtp;

          ! Thin KDK
!          atp= (one+p0sp*C0y*this%dz*this%dz*half)*a2 + (p0sp*this%dz)*c2
!          btp= (one+p0sp*C0y*this%dz*this%dz*half)*b2 + (p0sp*this%dz)*d2
!          ctp= (one+p0sp*C0y*this%dz*this%dz*0.25_fPrec)*(C0y*this%dz)*a2 + (one+p0sp*C0y*this%dz*this%dz*half)*c2
!          dtp= (one+p0sp*C0y*this%dz*this%dz*0.25_fPrec)*(C0y*this%dz)*b2 + (one+p0sp*C0y*this%dz*this%dz*half)*d2
!          a2=atp; b2=btp; c2=ctp; d2=dtp;


          tp = (p0sp*C0y)*(this%dz*this%dz)
          tp1= one+tp*half
          tp2= p0sp*this%dz
          tp3=(one+tp*0.25_fPrec)*(C0y*this%dz)
          atp= tp1*a2 + tp2*c2; btp= tp1*b2 + tp2*d2;  ctp= tp3*a2 + tp1*c2; dtp= tp3*b2 + tp1*d2
          a2=atp; b2=btp; c2=ctp; d2=dtp;

        end do

        ! Inversion of the matrix
        ! ---------------------------------------------------------------------------------------- !
!        bcmad=one!/(a1*d1-b1*c1)
!        atp=d1*bcmad;  btp=-b1*bcmad;  ctp=-c1*bcmad;  dtp=a1*bcmad;
!        a1=atp; b1=btp; c1=ctp; d1=dtp
        atp=d1; btp=-b1; ctp=-c1; dtp=a1;
        a1=atp; b1=btp;  c1=ctp;  d1=dtp;

!        bcmad=one!/(a2*d2-b2*c2)
!        atp=d2*bcmad;  btp=-b2*bcmad;  ctp=-c2*bcmad;  dtp=a2*bcmad;
!        a2=atp; b2=btp; c2=ctp; d2=dtp
        atp=d2; btp=-b2; ctp=-c2; dtp=a2;
        a2=atp; b2=btp;  c2=ctp;  d2=dtp;

        this%TAQx(1,1,k)= a1; this%TAQx(1,2,k)= b1; this%TAQx(2,1,k)= c1; this%TAQx(2,2,k)= d1
        this%TAQy(1,1,k)= a2; this%TAQy(1,2,k)= b2; this%TAQy(2,1,k)= c2; this%TAQy(2,2,k)= d2

      end do
    endif
  end subroutine TAQgen


  ! ================================================================ !
  !  Lie2
  ! ================================================================ !
  subroutine TLie2(this,x,px,y,py,zb,deltap1)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : half
    use mathlib_bouncer,     only : log10_mb

    implicit none

    ! interface variables
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(inout) :: this
    real(kind=fPrec),       intent(in)    :: deltap1
    real(kind=fPrec),       intent(inout) :: zb
    real(kind=fPrec),       intent(inout) :: x,px,y,py ! Transverse canonical parameter in new referencial

    ! Subroutine variables
    ! ---------------------------------------------------------------------------------------------- !
    integer          :: i       ! Iterator
    real(kind=fPrec) :: valA    ! Coefficient given by Horner subroutine
    real(kind=fPrec) :: dzover2 ! Dsigma/2
    real(kind=fPrec) :: g2d2inv ! 1/(delta+1)
    real(kind=fPrec) :: log_tmp

    ! Initialize the step size in z
    ! ---------------------------------------------------------------------------------------------- !
    dzover2=half*this%dz !loc

    ! Check size vectors xpow and ypow   (Prevent SIGFPE)
    ! ---------------------------------------------------------------------------------------------- !
    log_tmp=abs(log10_mb(x))
    if (log_tmp*(this%n)>230) then
      this%max_i=230/log_tmp-1
    else
      this%max_i=this%n
    endif
    log_tmp=abs(log10_mb(y))
    if (log_tmp*(this%m)>230) then
      this%max_j=230/log_tmp-1
    else
      this%max_j=this%m
    endif

    ! Initialize the step size in z
    ! ---------------------------------------------------------------------------------------------- !
    do i=1,this%s !-1
      !             * h1
!      zb=zb-dzover2                                    ! ???????????

      !             * h2
      call this%HornerDX_Az(x, y, i, valA)
      px=px+dzover2*valA
      call this%HornerDY_Az(x, y, i, valA)
      py=py+dzover2*valA

      !             * h3
      !                 ** Change of variable
      call this%Horner2D_Ax(x, y, i, valA)
      px=px-valA
      call this%HornerDYIntX_Ax(x, y, i, valA)
      py=py-valA

      x=x+dzover2*(px*deltap1)
!      zb=zb-dzover2*px*px*g2d2inv


      !                 ** Change of variable
      call this%Horner2D_Ax(x, y, i, valA)
      px=px+valA
      call this%HornerDYIntX_Ax(x, y, i, valA)
      py=py+valA


      !             * h4
      !                 ** Change of variable
      call this%HornerDXIntY_Ay(x, y, i, valA)
      px=px-valA
      call this%Horner2D_Ay(x, y, i, valA)
      py=py-valA

      y=y+this%dz*(py*deltap1)
!      zb=zb-dz*py*py*g2d2inv


      !                 ** Change of variable
      call this%HornerDXIntY_Ay(x, y, i, valA)
      px=px+valA
      call this%Horner2D_Ay(x, y, i, valA)
      py=py+valA


      !             * h3
      !                 ** Change of variable
      call this%Horner2D_Ax(x, y, i, valA)
      px=px-valA
      call this%HornerDYIntX_Ax(x, y, i, valA)
      py=py-valA

      x=x+dzover2*(px*deltap1)
!      zb=zb-dzover2*px*px*g2d2inv


      !                 ** Change of variable
      call this%Horner2D_Ax(x, y, i, valA)
      px=px+valA
      call this%HornerDYIntX_Ax(x, y, i, valA)
      py=py+valA


      !             * h2
      call this%HornerDX_Az(x, y, i, valA)
      px=px+dzover2*valA
      call this%HornerDY_Az(x, y, i, valA)
      py=py+dzover2*valA

      !             * h1
!      zb=zb-dzover2
    end do
  end subroutine TLie2


  subroutine HornerDX_Az(this,x,y,z,resultat)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, one

    implicit none

    ! Input parameter
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(in)  :: this
    integer,                intent(in)  :: z        ! Position in z
    real(kind=fPrec),       intent(in)  :: x,y      ! Transverse canonical parameter in new referencial
    real(kind=fPrec),       intent(out) :: resultat ! Return value

    ! Subroutine parameter
    ! ---------------------------------------------------------------------------------------------- !
    integer :: l                                        ! Indice for x and y in loop
    integer :: max_i_l,max_j_l                          ! To prevent the SIGFPE
    real(kind=fPrec) :: xpow(0:this%n-1),ypow(0:this%m) ! Power of x and y
    real(kind=fPrec) :: r0                              ! Variable for the vector potential computation
    real(kind=fPrec) :: log_tmp,di

    ! Initialize vectors xpow and ypow
    ! ---------------------------------------------------------------------------------------------- !
    max_i_l=this%max_i-1;  max_j_l=this%max_j

    xpow(0)=one
    do l=1,max_i_l
      xpow(l)=xpow(l-1)*x
    enddo
    ypow(0)=one
    do l=1,max_j_l
      ypow(l)=ypow(l-1)*y
    enddo

    ! Compute return value
    ! ---------------------------------------------------------------------------------------------- !
    r0=zero
    do l=1,this%lz
      di=real(this%ij_TAz(1,l,z), fPrec)
      if ((di>zero).and.(this%ij_TAz(1,l,z)-1<=max_i_l).and.(this%ij_TAz(2,l,z)<=max_j_l)) then
        r0 = r0 + di*(xpow(this%ij_TAz(1,l,z)-1)*(ypow(this%ij_TAz(2,l,z))*this%TAz(l,z)))
      endif
    enddo
    resultat=r0
  end subroutine HornerDX_Az


  subroutine HornerDY_Az(this,x,y,z,resultat)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, one

    implicit none

    ! Input parameter
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(in)  :: this
    integer,                intent(in)  :: z        ! Position in z
    real(kind=fPrec),       intent(in)  :: x,y      ! Transverse canonical parameter in new referencial
    real(kind=fPrec),       intent(out) :: resultat ! Return value

    ! Subroutine parameter
    ! ---------------------------------------------------------------------------------------------- !
    integer :: l                                        ! Indice for x and y in loop
    integer :: max_i_l,max_j_l                          ! To prevent the SIGFPE
    real(kind=fPrec) :: xpow(0:this%n),ypow(0:this%m-1) ! Power of x and y
    real(kind=fPrec) :: r0                              ! Variable for the vector potential computation
    real(kind=fPrec) :: log_tmp,dj

    ! Initialize vectors xpow and ypow
    ! ---------------------------------------------------------------------------------------------- !
    max_i_l=this%max_i;    max_j_l=this%max_j-1

    xpow(0)=one
    do l=1,max_i_l
      xpow(l)=xpow(l-1)*x
    enddo
    ypow(0)=one
    do l=1,max_j_l
      ypow(l)=ypow(l-1)*y
    enddo

    ! Compute return value
    ! ---------------------------------------------------------------------------------------------- !
    r0=zero
    do l=1,this%lz
      dj=real(this%ij_TAz(2,l,z), fPrec)
      if ((dj>zero).and.(this%ij_TAz(1,l,z)<=max_i_l).and.(this%ij_TAz(2,l,z)-1<=max_j_l)) then
        r0 = r0 + dj*(xpow(this%ij_TAz(1,l,z))*(ypow(this%ij_TAz(2,l,z)-1)*this%TAz(l,z)))
      endif
    enddo
    resultat=r0
  end subroutine HornerDY_Az


  subroutine Horner2D_Ax(this,x,y,z,resultat)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, one

    implicit none

    ! Input parameter
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(in)  :: this
    integer,                intent(in)  :: z        ! Position in z
    real(kind=fPrec),       intent(in)  :: x,y      ! Transverse canonical parameter in new referencial
    real(kind=fPrec),       intent(out) :: resultat ! Return value

    ! Subroutine parameter
    ! ---------------------------------------------------------------------------------------------- !
    integer :: l                                      ! Indice for x and y in loop
    integer :: max_i_l,max_j_l                        ! To prevent the SIGFPE
    real(kind=fPrec) :: xpow(0:this%n),ypow(0:this%m) ! Power of x and y
    real(kind=fPrec) :: r0                            ! Variable for the vector potential computation
    real(kind=fPrec) :: log_tmp

    ! Initialize vectors xpow and ypow
    ! ---------------------------------------------------------------------------------------------- !
    max_i_l=this%max_i;    max_j_l=this%max_j

    xpow(0)=one
    do l=1,max_i_l
      xpow(l)=xpow(l-1)*x
    enddo
    ypow(0)=one
    do l=1,max_j_l
      ypow(l)=ypow(l-1)*y
    enddo

    ! Compute return value
    ! ---------------------------------------------------------------------------------------------- !
    r0=zero
    do l=1,this%lx
      if ((this%ij_TAx(1,l,z)<=max_i_l).and.(this%ij_TAx(2,l,z)<=max_j_l)) then
        r0 = r0 + xpow(this%ij_TAx(1,l,z))*(ypow(this%ij_TAx(2,l,z))*this%TAx(l,z))
      endif
    enddo
    resultat=r0

  end subroutine Horner2D_Ax


  subroutine Horner2D_Ay(this,x,y,z,resultat)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, one

    implicit none

    ! Input parameter
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(in)  :: this
    integer,                intent(in)  :: z        ! Position in z
    real(kind=fPrec),       intent(in)  :: x,y      ! Transverse canonical parameter in new referencial
    real(kind=fPrec),       intent(out) :: resultat ! Return value

    ! Subroutine parameter
    ! ---------------------------------------------------------------------------------------------- !
    integer :: l                                      ! Indice for x and y in loop
    integer :: max_i_l,max_j_l                        ! To prevent the SIGFPE
    real(kind=fPrec) :: xpow(0:this%n),ypow(0:this%m) ! Power of x and y
    real(kind=fPrec) :: r0                            ! Variable for the vector potential computation
    real(kind=fPrec) :: log_tmp

    ! Initialize vectors xpow and ypow
    ! ---------------------------------------------------------------------------------------------- !
    max_i_l=this%max_i;    max_j_l=this%max_j

    xpow(0)=one
    do l=1,max_i_l
      xpow(l)=xpow(l-1)*x
    enddo
    ypow(0)=one
    do l=1,max_j_l
      ypow(l)=ypow(l-1)*y
    enddo

    ! Compute return value
    ! ---------------------------------------------------------------------------------------------- !
    r0=zero
    do l=1,this%ly
      if ((this%ij_TAy(1,l,z)<=max_i_l).and.(this%ij_TAy(2,l,z)<=max_j_l)) then
        r0 = r0 + xpow(this%ij_TAy(1,l,z))*(ypow(this%ij_TAy(2,l,z))*this%TAy(l,z))
      endif
    enddo
    resultat=r0

  end subroutine Horner2D_AY


  subroutine HornerDYIntX_Ax(this,x,y,z,resultat)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, one

    implicit none

    ! Input parameter
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(in)  :: this
    integer,                intent(in)  :: z        ! Position in z
    real(kind=fPrec),       intent(in)  :: x,y      ! Transverse canonical parameter in new referencial
    real(kind=fPrec),       intent(out) :: resultat ! Return value

    ! Subroutine parameter
    ! ---------------------------------------------------------------------------------------------- !
    integer :: l                                          ! Indice for x and y in loop
    integer :: max_i_l,max_j_l                            ! To prevent the SIGFPE
    real(kind=fPrec) :: xpow(0:this%n+1),ypow(0:this%m-1) ! Power of x and y
    real(kind=fPrec) :: r0                                ! Variable for vector potential computation
    real(kind=fPrec) :: log_tmp,di,dj

    ! Initialize vectors xpow and ypow
    ! ---------------------------------------------------------------------------------------------- !
    max_i_l=this%max_i+1;  max_j_l=this%max_j-1

    xpow(0)=one
    do l=1,max_i_l
      xpow(l)=xpow(l-1)*x
    enddo
    ypow(0)=one
    do l=1,max_j_l
      ypow(l)=ypow(l-1)*y
    enddo

    ! Compute return value
    ! ---------------------------------------------------------------------------------------------- !
    r0=zero
    do l=1,this%lx
      di=real(this%ij_TAx(1,l,z) + 1, fPrec)
      dj=real(this%ij_TAx(2,l,z)    , fPrec)
      if ((dj>0).and.(this%ij_TAx(1,l,z)+1<=max_i_l).and.(this%ij_TAx(2,l,z)-1<=max_j_l)) then
        r0 = r0 + (dj*(xpow(this%ij_TAx(1,l,z) + 1)*(ypow(this%ij_TAx(2,l,z) - 1)*this%TAx(l,z))))/di
      endif
    enddo
    resultat=r0
  end subroutine HornerDYIntX_Ax


  subroutine HornerDXIntY_Ay(this,x,y,z,resultat)
    ! Mod from SixTrack
    ! ---------------------------------------------------------------------------------------------- !
    use numerical_constants, only : zero, one

    implicit none

    ! Input parameter
    ! ---------------------------------------------------------------------------------------------- !
    class(ffTable_n_Track), intent(in)  :: this
    integer,                intent(in)  :: z        ! Position in z
    real(kind=fPrec),       intent(in)  :: x,y      ! Transverse canonical parameter in new referencial
    real(kind=fPrec),       intent(out) :: resultat ! Return value

    ! Subroutine parameter
    ! ---------------------------------------------------------------------------------------------- !
    integer :: l                                          ! Indice for x and y in loop
    integer :: max_i_l,max_j_l                            ! To prevent the SIGFPE
    real(kind=fPrec) :: xpow(0:this%n-1),ypow(0:this%m+1) ! Power of x and y
    real(kind=fPrec) :: r0                                ! Variable for vector potential computation
    real(kind=fPrec) :: log_tmp,di,dj

    ! Initialize vectors xpow and ypow
    ! ---------------------------------------------------------------------------------------------- !
    max_i_l=this%max_i-1;  max_j_l=this%max_j+1

    xpow(0)=one
    do l=1,max_i_l
      xpow(l)=xpow(l-1)*x
    enddo
    ypow(0)=one
    do l=1,max_j_l
      ypow(l)=ypow(l-1)*y
    enddo

    ! Compute return value
    ! ---------------------------------------------------------------------------------------------- !
    r0=zero
    do l=1,this%ly
      di=real(this%ij_TAy(1,l,z)    , fPrec)
      dj=real(this%ij_TAy(2,l,z) + 1, fPrec)
      if ((di>0).and.(this%ij_TAy(1,l,z)-1<=max_i_l).and.(this%ij_TAy(2,l,z)+1<=max_j_l)) then
        r0 = r0 + (di*(xpow(this%ij_TAy(1,l,z) - 1)*(ypow(this%ij_TAy(2,l,z) + 1)*this%TAy(l,z))))/dj
      endif
    enddo
    resultat=r0

  end subroutine HornerDXIntY_Ay

end module ffTable_n_Tracks

