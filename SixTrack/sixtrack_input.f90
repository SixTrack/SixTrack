! ================================================================================================ !
!  SixTrack Input
! ~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-18
! ================================================================================================ !
module sixtrack_input
  
  use crcoall
  use numerical_constants
  use mod_alloc
  use string_tools
  use mod_common
  
  implicit none
  
  ! Record of encountered blocks
  character(len=:), allocatable, public, save :: sixin_cBlock(:) ! Name of block
  integer,          allocatable, public, save :: sixin_iBlock(:) ! Line of block
  logical,          allocatable, public, save :: sixin_lBlock(:) ! Block closed
  integer,                       public, save :: sixin_nBlock    ! Number of blocks
  
  ! Single Element Variables
  integer,                       public, save :: sixin_nSing
  integer,                       public, save :: sixin_ncy2
  character(len=:), allocatable, public, save :: sixin_bez0(:) ! (str_maxName)(nele)
  character(len=3), parameter,   public       :: sixin_cavity = "CAV"
  
  ! Block Definition Variables
  integer,                       public, save :: sixin_nBloc
  character(len=:), allocatable, public, save :: sixin_beze(:,:)
  character(len=:), allocatable, public, save :: sixin_ilm(:)
  integer,                       public, save :: sixin_k0
  
  ! Structure Input Variables
  integer,                       public, save :: sixin_nStru
  integer,                       public, save :: sixin_icy
  character(len=2), parameter,   public       :: sixin_go = "GO"
  
contains

subroutine sixin_checkBlock(blockName, blockOpened, blockClosed, blockLine)
  
  character(len=*), intent(in)  :: blockName
  logical,          intent(out) :: blockOpened
  logical,          intent(out) :: blockClosed
  integer,          intent(out) :: blockLine
  
  integer i
  
  blockLine   = 0
  blockOpened = .false.
  blockClosed = .false.
  
  if(len(blockName) < 4) then
    write(lout,"(a)") "INPUT> WARNING Unknown blockname '"//blockName//"'"
    return
  end if
  
  ! We should of course not try to open a block "NEXT"
  if(blockName == "NEXT") return
  
  if(allocated(sixin_cBlock)) then ! Assuming the others are allocated too
    ! Check status of block, and increment line number if it exists
    do i=1,sixin_nBlock
      if(sixin_cBlock(i) == blockName) then
        ! Block already opened, so don't add it to the list.
        ! Just increment line number and return.
        blockOpened     = .true.
        blockClosed     = sixin_lBlock(i)
        blockLine       = sixin_iBlock(i) + 1
        sixin_iBlock(i) = blockLine
        return
      end if
    end do
    sixin_nBlock = sixin_nBlock + 1
  else
    ! If the array isn't allocated, it obviously doesn't contain anything
    sixin_nBlock = 1
  end if
  
  ! New block. Expand the arrays.
  call alloc(sixin_cBlock,4,sixin_nBlock,"    ", "sixin_cBlock")
  call alloc(sixin_iBlock,  sixin_nBlock,0,      "sixin_iBlock")
  call alloc(sixin_lBlock,  sixin_nBlock,.false.,"sixin_lBlock")
  
  sixin_cBlock(sixin_nBlock)(1:4) = blockName(1:4)
  sixin_iBlock(sixin_nBlock)      = 0
  sixin_lBlock(sixin_nBlock)      = .false.
  
  blockOpened = .true.
  
  write(lout,"(a)") "INPUT> Opened block '"//blockName//"'"
  
end subroutine sixin_checkBlock

subroutine sixin_closeBlock(blockName)
  
  character(len=*), intent(in) :: blockName
  
  integer i
  
  do i=1,sixin_nBlock
    if(sixin_cBlock(i) == blockName) then
      sixin_lBlock(i) = .true.
      return
    end if
  end do
  
  write(lout,"(a)") "INPUT> Closed block '"//blockName//"'"
  
end subroutine sixin_closeBlock

subroutine sixin_blockReport
  
  integer i
  
  write(lout,"(a)") "INPUT> Finished parsing input file(s)."
  write(lout,"(a)") "INPUT> Parsed the following blocks:"
  do i=1,sixin_nBlock
    write(lout,"(a,i0,a)") "INPUT> * "//sixin_cBlock(i)//" block with ",sixin_iBlock(i)," lines"
  end do
  
end subroutine sixin_blockReport

subroutine sixin_parseInputLineSING(inLine, iLine, iErr)
  
  use parpro_scale
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: elemName
  integer nSplit
  
  integer i
  
  call chr_split(inLine, lnSplit, nSplit)
  
  if(nSplit <= 2) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element line must have more than 2 values, got ",nSplit
    iErr = .true.
    return
  end if
  
  elemName = chr_trimZero(lnSplit(1))
  if(len(elemName) > str_maxName) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element name too long. Max length is ",str_maxName
    iErr = .true.
    return
  end if
  
  ! Check that the name is unique
  do i=1,sixin_nSing-1
    if(bez(i) == elemName) then
      write(lout,"(a,i0)") "GEOMETRY> ERROR Single element '"//elemName//"' is not unique."
      iErr = .true.
      return
    end if
  end do
  
  ! Save Values
  if(nSplit > 1) kz(sixin_nSing)   = chr_toInt(lnSplit(2))
  if(nSplit > 2) ed(sixin_nSing)   = chr_toReal(lnSplit(3))
  if(nSplit > 3) ek(sixin_nSing)   = chr_toReal(lnSplit(4))
  if(nSplit > 4) el(sixin_nSing)   = chr_toReal(lnSplit(5))
  if(nSplit > 5) bbbx(sixin_nSing) = chr_toReal(lnSplit(6))
  if(nSplit > 6) bbby(sixin_nSing) = chr_toReal(lnSplit(7))
  if(nSplit > 7) bbbs(sixin_nSing) = chr_toReal(lnSplit(8))
  
  if(kz(sixin_nSing) == 25) then
    ed(sixin_nSing) = ed(sixin_nSing)/two
    ek(sixin_nSing) = ek(sixin_nSing)/two
  endif
  
  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  if((kz(sixin_nSing) == 4 .or. kz(sixin_nSing) == 5) .and. abs(el(sixin_nSing)) > pieni) then
    ed(sixin_nSing) = -one*ed(sixin_nSing)
  end if
  
  ! CAVITIES
  if(abs(kz(sixin_nSing)) == 12) then
    if(abs(ed(sixin_nSing)) > pieni .and. abs(ek(sixin_nSing)) > pieni) then
      sixin_ncy2    = sixin_ncy2 + 1
      itionc(sixin_nSing) = kz(sixin_nSing)/abs(kz(sixin_nSing))
      kp(sixin_nSing)     = 6
    end if
  end if
  
  !--------------------------------------------
  ! Handled by initialize_element subroutine:
  !--------------------------------------------
  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  ! THIN LENS (+/- 1-10)
  ! MULTIPOLES (11)
  ! CAVITY (+/- 12)
  ! CRABCAVITY (23/-23) / CC multipoles order 2/3/4 (+/- 23/26/27/28)
  ! ELECTRON LENSE (29)
  call initialize_element(sixin_nSing,.true.)
  
  ! ACDIPOLE
  if(abs(kz(sixin_nSing)) == 16) then
    if(abs(ed(sixin_nSing)) <= pieni) then
      kz(sixin_nSing) = 0
      ed(sixin_nSing) = zero
      ek(sixin_nSing) = zero
      el(sixin_nSing) = zero
    else
      acdipph(sixin_nSing) = el(sixin_nSing)
      el(sixin_nSing)      = zero
    end if
  end if
  
  ! General
  if(abs(el(sixin_nSing)) > pieni .and. kz(sixin_nSing) /= 0) then
    ithick = 1
  end if
  
  ! Expand Arrays
  if(sixin_nSing > nele-2) then
    call expand_arrays(nele+50, npart, nblz, nblo)
    ! The value of nele will have been updated here
    call resize(sixin_bez0, str_maxName, nele, str_nmZeros, "sixin_bez0")
  end if
  
  if(abs(kz(sixin_nSing)) /= 12 .or. (abs(kz(sixin_nSing)) == 12 .and. sixin_ncy2 == 0)) then
    kp(sixin_nSing) = 0
  end if
  
  bez(sixin_nSing)        = elemName
  sixin_bez0(sixin_nSing) = elemName
  
  !If no active RF cavities are seen so far in the single element list,
  ! add a CAV element to the end of the list.
  ! This is then overwritten when reading the next element, so that if
  ! and only if no active RF cavities are found, a CAV element can be
  ! used in the structure to enable 6D tracking using the parameters
  ! from the SYNC block.
  if(sixin_ncy2 == 0) then
    sixin_nSing             = sixin_nSing + 1
    il                      = sixin_nSing
    bez(sixin_nSing)        = sixin_cavity
    sixin_bez0(sixin_nSing) = sixin_cavity
    kp(sixin_nSing)         = 6
  else
    il          = sixin_nSing
    sixin_nSing = sixin_nSing + 1
  end if
  
end subroutine sixin_parseInputLineSING

subroutine sixin_parseInputLineBLOC(inLine, iLine, iErr)
  
  use parpro_scale
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: blocName
  integer nSplit
  
  integer i, j, ka, ke
  logical eFound, isCont
  character(len=str_maxName) ilm0(40)
  
  call chr_split(inLine, lnSplit, nSplit, isCont)
  
  if(nSplit < 2 .and. .not. isCont) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Block definition line must be at least 2 values, got ",nSplit
    iErr = .true.
    return
  end if
  
  if(nSplit > 40) then
    write(lout,"(a)") "GEOMETRY> ERROR Block definition cannot have more then 40 elements."
    iErr = .true.
    return
  end if
  
  ! If first line, read super period information
  if(iLine == 1) then
    mper = chr_toInt(lnSplit(1))
    if(mper > nper) then
      write(lout,"(a,i0)") "GEOMETRY> ERROR Block definition number of super periods is too large. Max value is ",nper
      iErr = .true.
      return
    end if
    if(mper > nSplit+1) then
      write(lout,"(a)") "GEOMETRY> ERROR Block definition number of super periods does not match the number of values."
      iErr = .true.
      return
    end if
    do i=1,mper
      msym(i) = chr_toInt(lnSplit(i+1))
    end do
    
    ! Init variables
    call alloc(sixin_ilm,  str_maxName,       nelb, str_nmSpace, "sixin_ilm")
    call alloc(sixin_beze, str_maxName, nblo, nelb, str_nmSpace, "sixin_beze")
    
    ! No need to parse anything more for this line
    return
  end if
  
  ! Parse normal line, iLine > 1
  do i=1,40
    ilm0(i) = str_nmSpace
  end do
  
  if(isCont) then
    ! This line continues the previous block
    blocName = str_nmSpace
    do i=1,nSplit
      ilm0(i) = chr_trimZero(lnSplit(i))
    end do
  else
    blocName = chr_trimZero(lnSplit(1))
    do i=1,nSplit-1
      ilm0(i) = chr_trimZero(lnSplit(i+1))
    end do
  end if
  
  if(blocName /= str_nmSpace) then
    sixin_nBloc = sixin_nBloc + 1 ! Current BLOC number
    if(sixin_nBloc > nblo-1) then
      call expand_arrays(nele, npart, nblz, nblo+50)
      call alloc(sixin_beze, str_maxName, nblo, nelb, str_nmSpace, "sixin_beze")
    end if
    bezb(sixin_nBloc) = blocName
    sixin_k0          = 0
    mblo              = sixin_nBloc ! Update total number of BLOCs
  end if
  
  ka = sixin_k0 + 1
  ke = sixin_k0 + 40
  
  do i=ka, ke
    if(i > nelb) then
      write(lout,"(a,2(i0,a))") "GEOMETRY> ERROR Block definitions can only have ",nelb," elements. ",i," given."
      iErr = .true.
      return
    end if
    sixin_ilm(i) = ilm0(i-sixin_k0)
    if(sixin_ilm(i) == str_nmSpace) exit
    
    mel(sixin_nBloc)          = i            ! Number of single elements in this block
    sixin_beze(sixin_nBloc,i) = sixin_ilm(i) ! Name of the current single element
    
    ! Search for the single element idx j
    eFound = .false.
    do j=1,il 
      if(sixin_bez0(j) == sixin_ilm(i)) then
        eFound = .true.
        exit
      end if
    end do
    if(eFound) then
      ! Block sixin_nBloc / sub-element i has single element index j
      mtyp(sixin_nBloc,i) = j
      if(kz(j) /= 8) then
        ! Count block length (kz=8 -> edge focusing->skip!)
        elbe(sixin_nBloc) = elbe(sixin_nBloc) + el(j)
      end if
    else
      write(lout,"(a)") "GEOMETRY> ERROR Unknown element '"//sixin_ilm(i)//"' in block definitions."
      iErr = .true.
      return
    end if
  end do
  
  sixin_k0 = i-1
  
end subroutine sixin_parseInputLineBLOC

subroutine sixin_parseInputLineSTRU(inLine, iLine, iErr)

  use parpro_scale
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: expLine
  integer nSplit
  
  integer i, j
  character(len=str_maxName) ilm0(40)
  
  do i=1,40
    ilm0(i) = str_nmSpace
  end do
  
  expLine = chr_expandBrackets(inLine)
  call chr_split(expLine, lnSplit, nSplit)
  
  if(nSplit > 40) then
    write(lout,"(a)") "GEOMETRY> ERROR Structure input line cannot have more then 40 elements."
    iErr = .true.
    return
  end if
  
  do i=1,nSplit
    ilm0(i) = chr_trimZero(lnSplit(i))
  end do
  
  do i=1,40
    
    if(ilm0(i) == str_nmSpace) cycle
    if(ilm0(i) == sixin_go) then
      kanf = sixin_nStru + 1
      cycle
    end if
    
    sixin_nStru = sixin_nStru + 1
    
    do j=1,mblo ! is it a BLOC?
      if(bezb(j) == ilm0(i)) then
        ic(sixin_nStru) = j
        cycle
      end if
    end do
    
    do j=1,il ! is it a SINGLE ELEMENT?
      if(sixin_bez0(j) == ilm0(i)) then
        ic(sixin_nStru) = j+nblo
        if(sixin_bez0(j) == sixin_cavity) sixin_icy = sixin_icy+1
        cycle
      end if
    end do
  end do
  
  mbloz = sixin_nStru
!     if(mbloz.gt.nblz-2) call prror(21)
  
end subroutine sixin_parseInputLineSTRU

subroutine sixin_parseInputLineDISP(inLine, iErr)
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: elemName
  integer nSplit
  
  integer i
  real(kind=fPrec) xpl0, xrms0, zpl0, zrms0
  
  call chr_split(inLine, lnSplit, nSplit)
  
  xpl0  = zero
  xrms0 = zero
  zpl0  = zero
  zrms0 = zero
  
  if(nSplit < 2) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Displacement of element line must have more than 1 values, got ",nSplit
    iErr = .true.
    return
  end if
  
  elemName = chr_trimZero(lnSplit(1))
  if(len(elemName) > str_maxName) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Displacement of element name too long. Max length is ",str_maxName
    iErr = .true.
    return
  end if
  
  ! Save Values
  if(nSplit > 1) xpl0  = chr_toReal(lnSplit(2))
  if(nSplit > 2) xrms0 = chr_toReal(lnSplit(3))
  if(nSplit > 3) zpl0  = chr_toReal(lnSplit(4))
  if(nSplit > 4) zrms0 = chr_toReal(lnSplit(5))
  
  do i=1,il
    if(elemName /= bez(i)) cycle
    
    xpl(i)  = xpl0
    xrms(i) = xrms0
    zpl(i)  = zpl0
    zrms(i) = zrms0
    
    ! Insertion for AC dipole
    if(abs(kz(i)) == 16) then
      nturn1(i) = int(xpl0)
      nturn2(i) = int(xrms0)
      nturn3(i) = int(zpl0)
      nturn4(i) = int(zrms0)
      xpl(i)    = zero
      xrms(i)   = zero
      zpl(i)    = zero
      zrms(i)   = zero
      if(xrms0 == zero .and. zpl0 == zero .and. zrms0 == zero) then
        write(lout,"(a)") "INPUT> INFO DISP Block: AC dipole disregarded, 0-length."
        kz(i) = 0
        ed(i) = zero
        ek(i) = zero
      end if
    end if
    
  end do
  
end subroutine sixin_parseInputLineDISP

end module sixtrack_input
