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
! character(len=3), parameter,   public,      :: sixin_cavity = "CAV"
  
  ! Block Definition Variables
  integer,                       public, save :: sixin_nBloc
  character(len=:), allocatable, public, save :: sixin_beze(:,:)
  character(len=:), allocatable, public, save :: sixin_ilm(:)
  integer,                       public, save :: sixin_k0
  
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
  if(len(elemName) > str_maxName) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element name too long. Max length is ",str_maxName
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
    call resize(sixin_bez0, str_maxName, nele, str_nmZeros, "sixin_bez0")
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
    sixin_nBloc = 0
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
  character(len=:), allocatable   :: blocName
  integer nSplit
  
  integer i
  character(len=str_maxName) ilm0(40)
  
  call chr_split(inLine, lnSplit, nSplit)
  
  do i=1,40
    ilm0(i) = str_nmSpace
  end do
  
!     i2=1
!     ! Look for repetition with syntax N( ... )
!     do 420 ii=1,80
!       if(ch(ii:ii).eq.kl) then !kl='('
!         if(ii.gt.1) then

!           do jj=1,ii-1
!             if(ch(jj:jj).ne.' ') goto 380
!           end do

!         endif
!         iw=1
!         goto 390
! 380     read(ch(:ii-1),*) iw
! 390     ia=i
!         iw0=iw-1
!         i2=ii+1
!         goto 430
!       endif
!       if(ch(ii:ii).eq.kr) then !kr=')'
!         if(iw0.le.0) goto 330
!         idi=i-ia
!         do 410 k=1,iw0
!           do j=1,idi
!             ic(i+j)=ic(i+j-idi)
!           end do
!           i=i+idi
! 410     continue
!         mbloz=i
!         goto 330
!       endif
! 420 continue
!     ! Create the structure
! 430 call intepr(3,i2,ch,ch1)
! ! reading character strings so OK
!     read(ch1,*) (ilm0(k),k=1,40)
!     do 490 k=1,40
!       if(ilm0(k).eq.idum) goto 490
!       if(ilm0(k).eq.go) goto 480
!       i=i+1
!       do 440 j=1,mblo !is it a BLOC?
!         if(bezb(j).eq.ilm0(k)) goto 470
! 440   continue
!       do 450 l=1,il   !is it a SINGLE ELEMENT?
!         if(sixin_bez0(l).eq.ilm0(k)) goto 460
! 450   continue
!       ! It was neither BLOC or SINGLE ELEMENT! ERROR!
!       erbez=ilm0(k)
!       call prror(20)
      
!       ! Handle SINGLE ELEMENT
! 460   continue
!       ic(i)=l+nblo
!       if(sixin_bez0(l).eq.cavi) icy=icy+1
!       goto 490

!       !Handle BLOC
! 470   ic(i)=j
!       goto 490
!       !Handle GO
! 480   kanf=i+1
! 490 continue
!     mbloz=i
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
