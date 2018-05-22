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
  integer,                       public, save :: sixin_ncy2
  character(len=:), allocatable, public, save :: sixin_bez0(:) ! (max_name_len)(nele)
  ! character(len=3), parameter,   public,      :: sixin_cavity = "CAV"
  
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

subroutine sixin_parseInputLineSING(inLine, iElem, iErr)
  
  use parpro_scale
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iElem
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: elemName
  integer nSplit
  
  integer i
  
  call chr_split(inLine, lnSplit, nSplit)
  
  ! write(lout,"(a,i0,a)") "SING> Line(",nSplit,"): '"//chr_trimZero(inLine)//"'"
  ! do i=1,nSplit
  !   write(lout,"(a,i3,a)") "SING> ",i,": '"//chr_trimZero(lnSplit(i))//"'"
  ! end do
  
  if(nSplit <= 2) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element line must have more than 2 values, got ",nSplit
    iErr = .true.
    return
  end if
  
  elemName = chr_trimZero(lnSplit(1))
  if(len(elemName) > max_name_len) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element name too long. Max length is ",max_name_len
    iErr = .true.
    return
  end if
  
  ! Check that the name is unique
  do i=1,iElem-1
    if(bez(i) == elemName) then
      write(lout,"(a,i0)") "GEOMETRY> ERROR Single element '"//elemName//"' is not unique."
      iErr = .true.
      return
    end if
  end do
  
  ! Save Values
  if(nSplit > 1) kz(iElem)   = chr_toInt(lnSplit(2))
  if(nSplit > 2) ed(iElem)   = chr_toReal(lnSplit(3))
  if(nSplit > 3) ek(iElem)   = chr_toReal(lnSplit(4))
  if(nSplit > 4) el(iElem)   = chr_toReal(lnSplit(5))
  if(nSplit > 5) bbbx(iElem) = chr_toReal(lnSplit(6))
  if(nSplit > 6) bbby(iElem) = chr_toReal(lnSplit(7))
  if(nSplit > 7) bbbs(iElem) = chr_toReal(lnSplit(8))
  
  if(kz(iElem) == 25) then
    ed(iElem) = ed(iElem)/two
    ek(iElem) = ek(iElem)/two
  endif
  
  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  if((kz(iElem) == 4 .or. kz(iElem) == 5) .and. abs(el(iElem)) > pieni) then
    ed(iElem) = -one*ed(iElem)
  end if
  
  ! CAVITIES
  if(abs(kz(iElem)) == 12) then
    if(abs(ed(iElem)) > pieni .and. abs(ek(iElem)) > pieni) then
      sixin_ncy2    = sixin_ncy2 + 1
      itionc(iElem) = kz(iElem)/abs(kz(iElem))
      kp(iElem)     = 6
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
  call initialize_element(iElem,.true.)
  
  ! ACDIPOLE
  if(abs(kz(iElem)) == 16) then
    if(abs(ed(iElem)) <= pieni) then
      kz(iElem) = 0
      ed(iElem) = zero
      ek(iElem) = zero
      el(iElem) = zero
    else
      acdipph(iElem) = el(iElem)
      el(iElem)      = zero
    end if
  end if
  
  ! General
  if(abs(el(iElem)) > pieni .and. kz(iElem) /= 0) then
    ithick = 1
  end if
  
  ! Expand Arrays
  if(iElem > nele-2) then
    call expand_arrays(nele+50, npart, nblz, nblo)
    ! The value of nele will have been updated here
    call resize(sixin_bez0, max_name_len, nele, repeat(char(0),max_name_len), "sixin_bez0")
  end if
  
  if(abs(kz(iElem)) /= 12 .or. (abs(kz(iElem)) == 12 .and. sixin_ncy2 == 0)) then
    kp(iElem) = 0
  end if
  
  bez(iElem)        = elemName
  sixin_bez0(iElem) = elemName
  
  !If no active RF cavities are seen so far in the single element list,
  ! add a CAV element to the end of the list.
  ! This is then overwritten when reading the next element, so that if
  ! and only if no active RF cavities are found, a CAV element can be
  ! used in the structure to enable 6D tracking using the parameters
  ! from the SYNC block.
  if(sixin_ncy2 == 0) then
    iElem = iElem + 1
    il    = iElem
    bez(iElem)        = "CAV"
    sixin_bez0(iElem) = "CAV"
    kp(iElem)         = 6
  else
    il    = iElem
    iElem = iElem + 1
  end if
  
end subroutine sixin_parseInputLineSING

subroutine sixin_parseInputLineBLOC(inLine, iElem, iErr)
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iElem
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: elemName
  integer nSplit
  
  integer i
  
  call chr_split(inLine, lnSplit, nSplit)
  
  if(nSplit < 2) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Block definition line must be at least 2 values, got ",nSplit
    iErr = .true.
    return
  end if
  
  ! If first line, read super period information
  if(iElem == 0) then
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
    return ! No need to parse anything more for this line
  end if
  
  ! Parse normal line, iElem > 0
  
  
end subroutine sixin_parseInputLineBLOC

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
  if(len(elemName) > max_name_len) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Displacement of element name too long. Max length is ",max_name_len
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
