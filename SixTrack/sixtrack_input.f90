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
  use mod_commons
  
  implicit none
  
  ! Global Flags
  logical,                       public, save :: sixin_debug     ! Global debug flag
  
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
  
  interface sixin_echoVal
    module procedure sixin_echoVal_int
    module procedure sixin_echoVal_real64
  end interface sixin_echoVal
  
contains

! ================================================================================================ !
!  BLOCK PARSING RECORD
! ================================================================================================ !
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

subroutine sixin_debugInfo
  
  if(.not. sixin_debug) return
  
  write(lout,"(a)") "INPUT> SixTrack Input Debugging ENABLED"
#ifdef CRLIBM
  write(lout,"(a)") "INPUT> DEBUG CRLIBM is ON"
#else
  write(lout,"(a)") "INPUT> DEBUG CRLIBM is OFF"
#endif
#ifdef FIO
  write(lout,"(a)") "INPUT> DEBUG FIO is ON"
#else
  write(lout,"(a)") "INPUT> DEBUG FIO is OFF"
#endif
  
end subroutine sixin_debugInfo

subroutine sixin_echoVal_int(varName, varVal, blockName, lineNo)
  character(len=*), intent(in) :: varName
  integer,          intent(in) :: varVal
  character(len=*), intent(in) :: blockName
  integer,          intent(in) :: lineNo
  if(lineNo < 10) then
    write(lout,"(2(a,i0))") "INPUT> DEBUG "//blockName//":",lineNo,"  "//&
      chr_padSpace(varName,10)//" =  ",varVal
  else
    write(lout,"(2(a,i0))") "INPUT> DEBUG "//blockName//":",lineNo," "//&
      chr_padSpace(varName,10)//" =  ",varVal
  end if
end subroutine sixin_echoVal_int

subroutine sixin_echoVal_real64(varName, varVal, blockName, lineNo)
  character(len=*),  intent(in) :: varName
  real(kind=real64), intent(in) :: varVal
  character(len=*),  intent(in) :: blockName
  integer,           intent(in) :: lineNo
  if(lineNo < 10) then
    write(lout,"(a,i0,a,e22.15)") "INPUT> DEBUG "//blockName//":",lineNo,"  "//&
      chr_padSpace(varName,10)//" = ",varVal
  else
    write(lout,"(a,i0,a,e22.15)") "INPUT> DEBUG "//blockName//":",lineNo," "//&
      chr_padSpace(varName,10)//" = ",varVal
  end if
end subroutine sixin_echoVal_real64

! ================================================================================================ !
!  LINE PARSIN G ROUTINES
! ================================================================================================ !

! ================================================================================================ !
!  Parse Single Element Line
!  Rewritten from code from DATEN
! ================================================================================================ !
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
  
  elemName = chr_padSpace(chr_trimZero(lnSplit(1)),str_maxName)
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
  if(nSplit > 1) call chr_toInt(lnSplit(2),kz(sixin_nSing),iErr)
  if(nSplit > 2) call chr_toReal(lnSplit(3),ed(sixin_nSing),iErr)
  if(nSplit > 3) call chr_toReal(lnSplit(4),ek(sixin_nSing),iErr)
  if(nSplit > 4) call chr_toReal(lnSplit(5),el(sixin_nSing),iErr)
  if(nSplit > 5) call chr_toReal(lnSplit(6),bbbx(sixin_nSing),iErr)
  if(nSplit > 6) call chr_toReal(lnSplit(7),bbby(sixin_nSing),iErr)
  if(nSplit > 7) call chr_toReal(lnSplit(8),bbbs(sixin_nSing),iErr)
  
  if(iErr) return
  
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

! ================================================================================================ !
!  Parse Block Definitions Line
!  Rewritten from code from DATEN
! ================================================================================================ !
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
    call chr_toInt(lnSplit(1),mper,iErr)
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
      call chr_toInt(lnSplit(i+1),msym(i),iErr)
    end do
    if(iErr) return
    
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
  
  if(isCont) then                             ! This line continues the previous BLOC
    blocName = str_nmSpace                    ! No name returned, set an empty BLOC name
    do i=1,nSplit                             ! All elements are sub-elements. Save to buffer.
      ilm0(i) = chr_padSpace(chr_trimZero(lnSplit(i)),str_maxName)
    end do
  else                                        ! This is a new BLOC
    blocName = chr_padSpace(chr_trimZero(lnSplit(1)),str_maxName)
    do i=1,nSplit-1                           ! Save the rest to buffer
      ilm0(i) = chr_padSpace(chr_trimZero(lnSplit(i+1)),str_maxName)
    end do
  end if
  
  if(blocName /= str_nmSpace) then            ! We have a new BLOC
    sixin_nBloc = sixin_nBloc + 1             ! Increment the BLOC number
    if(sixin_nBloc > nblo-1) then             ! Expand arrays if needed
      ! For now, don't expand. This is incompatible with the offset in ic() array
      ! call expand_arrays(nele, npart, nblz, nblo+50)
      ! call alloc(sixin_beze, str_maxName, nblo, nelb, str_nmSpace, "sixin_beze")
      write(lout,"(a,i0)") "GEOMETRY> ERROR Too many block definitions. Max is ",nblo
      iErr = .true.
      return
    end if
    bezb(sixin_nBloc) = blocName              ! Set the BLOC name in bezb
    sixin_k0          = 0                     ! Reset the single element counter
    mblo              = sixin_nBloc           ! Update total number of BLOCs
  end if
  
  ka = sixin_k0 + 1
  ke = sixin_k0 + 40
  
  do i=ka, ke
    if(i > nelb) then
      write(lout,"(a,2(i0,a))") "GEOMETRY> ERROR Block definitions can only have ",&
        nelb," elements. ",i," given."
      iErr = .true.
      return
    end if
    sixin_ilm(i) = ilm0(i-sixin_k0)           ! Append to sub-element buffer
    if(sixin_ilm(i) == str_nmSpace) exit      ! No more sub-elements to append
    mel(sixin_nBloc)          = i             ! Update number of single elements in this block
    sixin_beze(sixin_nBloc,i) = sixin_ilm(i)  ! Name of the current single element
    
    eFound = .false.
    do j=1,il                                 ! Search for the single element index
      if(sixin_bez0(j) == sixin_ilm(i)) then
        eFound = .true.
        exit
      end if
    end do
    if(eFound) then                           ! Handle element found
      mtyp(sixin_nBloc,i) = j                 ! Save single element index
      if(kz(j) /= 8) then                     ! Count block length (kz=8 - edge focusing: skip)
        elbe(sixin_nBloc) = elbe(sixin_nBloc) + el(j)
      end if
    else                                      ! If not, the input is invalid
      write(lout,"(a)") "GEOMETRY> ERROR Unknown single element '"//&
        chr_strip(sixin_ilm(i))//"' in block definitions."
      iErr = .true.
      return
    end if
  end do
  
  sixin_k0 = i-1
  
end subroutine sixin_parseInputLineBLOC

! ================================================================================================ !
!  Parse Structure Input Line
!  Rewritten from code from DATEN
! ================================================================================================ !
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
    ilm0(i) = chr_padSpace(chr_trimZero(lnSplit(i)),str_maxName)
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
  if(mbloz > nblz-2) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Structure input block has too many elements. Max is ",(nblz-2)
    iErr = .true.
    return
    ! call prror(21)
  end if
  
end subroutine sixin_parseInputLineSTRU

! ================================================================================================ !
!  Parse Displacement Block Line
!  Rewritten from code from DATEN
! ================================================================================================ !
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
  if(nSplit > 1) call chr_toReal(lnSplit(2),xpl0,iErr)
  if(nSplit > 2) call chr_toReal(lnSplit(3),xrms0,iErr)
  if(nSplit > 3) call chr_toReal(lnSplit(4),zpl0,iErr)
  if(nSplit > 4) call chr_toReal(lnSplit(5),zrms0,iErr)
  if(iErr) return
  
  if(sixin_debug) then
    call sixin_echoVal("xpl0",xpl0,"DISP",0)
    call sixin_echoVal("xrms0",xrms0,"DISP",0)
    call sixin_echoVal("zpl0",zpl0,"DISP",0)
    call sixin_echoVal("zrms0",zrms0,"DISP",0)
  end if
  
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
        write(lout,"(a)") "INPUT> INFO DISP Block: AC dipole disregarded, 0 length."
        kz(i) = 0
        ed(i) = zero
        ek(i) = zero
      end if
    end if
    
  end do
  
end subroutine sixin_parseInputLineDISP

! ================================================================================================ !
!  Parse Initial Coordinates Line
!  Rewritten from code from DATEN
! ================================================================================================ !
subroutine sixin_parseInputLineINIT(inLine, iLine, iErr)

  use parpro_scale
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: expLine
  integer nSplit
  
  integer i
  
  call chr_split(inLine, lnSplit, nSplit)
  
  if(nSplit < 1) then
    write(lout,"(a,i0,a)") "PARAM> ERROR INIT block line ",iLine," did not receive any values."
    iErr = .true.
    return
  end if
  
  select case(iLine)
  case(1) ! Line One
    if(nSplit > 0) call chr_toInt(lnSplit(1),itra,iErr)  ! Number of particles
    if(nSplit > 1) call chr_toReal(lnSplit(2),chi0,iErr) ! Starting phase of the initial coordinate
    if(nSplit > 2) call chr_toReal(lnSplit(3),chid,iErr) ! Phase difference between particles
    if(nSplit > 3) call chr_toReal(lnSplit(4),rat,iErr)  ! Emittance ratio
    if(nSplit > 4) call chr_toInt(lnSplit(5),iver,iErr)  ! Vertical coordinates switch
    if(sixin_debug) then
      call sixin_echoVal("itra",itra,"INIT",iLine)
      call sixin_echoVal("chi0",chi0,"INIT",iLine)
      call sixin_echoVal("chid",chid,"INIT",iLine)
      call sixin_echoVal("rat", rat, "INIT",iLine)
      call sixin_echoVal("iver",iver,"INIT",iLine)
    end if
    if(iErr) return
    if(itra < 0 .or. itra > 2) then
      write(lout,"(a,i0,a)") "PARAM> ERROR INIT First value (itra) can only be 0, 1 or 2, but ",itra," given."
      iErr = .true.
      return
    end if
    if(iver < 0 .or. iver > 1) then
      write(lout,"(a,i0,a)") "PARAM> ERROR INIT Fifth value (iver) can only be 0 or 1, but ",iver," given."
      iErr = .true.
      return
    end if
    
  case(2)  ! x [mm] coordinate of particle 1
    call chr_toReal(lnSplit(1),exz(1,1),iErr)
    if(sixin_debug) call sixin_echoVal("exz(1,1)",exz(1,1),"INIT",iLine)
    if(iErr) return
    
  case(3)  ! xp [mrad] coordinate of particle 1
    call chr_toReal(lnSplit(1),exz(1,2),iErr)
    if(sixin_debug) call sixin_echoVal("exz(1,2)",exz(1,2),"INIT",iLine)
    if(iErr) return
    
  case(4)  ! y [mm] coordinate of particle 1
    call chr_toReal(lnSplit(1),exz(1,3),iErr)
    if(sixin_debug) call sixin_echoVal("exz(1,3)",exz(1,3),"INIT",iLine)
    if(iErr) return
    
  case(5)  ! yp [mrad] coordinate of particle 1
    call chr_toReal(lnSplit(1),exz(1,4),iErr)
    if(sixin_debug) call sixin_echoVal("exz(1,4)",exz(1,4),"INIT",iLine)
    if(iErr) return
    
  case(6)  ! Path length difference 1(sigma = s−v0*t) [mm] of particle 1
    call chr_toReal(lnSplit(1),exz(1,5),iErr)
    if(sixin_debug) call sixin_echoVal("exz(1,5)",exz(1,5),"INIT",iLine)
    if(iErr) return
    
  case(7)  ! dp/p0 of particle 1
    call chr_toReal(lnSplit(1),exz(1,6),iErr)
    if(sixin_debug) call sixin_echoVal("exz(1,6)",exz(1,6),"INIT",iLine)
    if(iErr) return
    
  case(8)  ! x [mm] coordinate of particle 2
    call chr_toReal(lnSplit(1),exz(2,1),iErr)
    if(sixin_debug) call sixin_echoVal("exz(2,1)",exz(2,1),"INIT",iLine)
    if(iErr) return
    
  case(9)  ! xp [mrad] coordinate of particle 2
    call chr_toReal(lnSplit(1),exz(2,2),iErr)
    if(sixin_debug) call sixin_echoVal("exz(2,2)",exz(2,2),"INIT",iLine)
    if(iErr) return
    
  case(10) ! y [mm] coordinate of particle 2
    call chr_toReal(lnSplit(1),exz(2,3),iErr)
    if(sixin_debug) call sixin_echoVal("exz(2,3)",exz(2,3),"INIT",iLine)
    if(iErr) return
    
  case(11) ! yp [mrad] coordinate of particle 2
    call chr_toReal(lnSplit(1),exz(2,4),iErr)
    if(sixin_debug) call sixin_echoVal("exz(2,4)",exz(2,4),"INIT",iLine)
    if(iErr) return
    
  case(12) ! Path length difference 1(sigma = s−v0*t) [mm] of particle 2
    call chr_toReal(lnSplit(1),exz(2,5),iErr)
    if(sixin_debug) call sixin_echoVal("exz(2,5)",exz(2,5),"INIT",iLine)
    if(iErr) return
    
  case(13) ! dp/p0 of particle 2
    call chr_toReal(lnSplit(1),exz(2,6),iErr)
    if(sixin_debug) call sixin_echoVal("exz(2,6)",exz(2,6),"INIT",iLine)
    if(iErr) return
    
  case(14) ! energy [MeV] of the reference particle
    call chr_toReal(lnSplit(1),e0,iErr)
    if(sixin_debug) call sixin_echoVal("e0",e0,"INIT",iLine)
    if(iErr) return
    
  case(15) ! energy [MeV] of particle 1
    call chr_toReal(lnSplit(1),ej(1),iErr)
    if(sixin_debug) call sixin_echoVal("ej(1)",ej(1),"INIT",iLine)
    if(iErr) return
    
  case(16) ! energy [MeV] of particle 2
    call chr_toReal(lnSplit(1),ej(2),iErr)
    if(sixin_debug) call sixin_echoVal("ej(2)",ej(2),"INIT",iLine)
    if(iErr) return
    
  case default
    write(lout,"(a,i0,a)") "PARAM> ERROR Unexpected line number ",iLine," in INIT block."
    iErr = .true.
    return
    
  end select
  
end subroutine sixin_parseInputLineINIT

! ================================================================================================ !
!  Parse Tracking Parameters Line
!  Rewritten from code from DATEN
! ================================================================================================ !
subroutine sixin_parseInputLineTRAC(inLine, iLine, iErr)

  use parpro_scale
  use mod_commont
  use mod_commond
  
  implicit none
  
  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr
  
  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: expLine
  integer nSplit
  
  call chr_split(inLine, lnSplit, nSplit)
  
  select case(iLine)
  case(1)
    if(nSplit < 7) then
      write(lout,"(a,i0,a)") "PARAM> ERROR TRAC block line 1 should be at least 7 values, but ",nSplit," given."
      iErr = .true.
      return
    end if
    if(nSplit > 0)  call chr_toInt(lnSplit(1),numl,iErr)     ! Number of turns in the forward direction
    if(nSplit > 1)  call chr_toInt(lnSplit(2),numlr,iErr)    ! Number of turns in the backward direction
    if(nSplit > 2)  call chr_toInt(lnSplit(3),napx,iErr)     ! Number of amplitude variations (i.e. particle pairs)
    if(nSplit > 3)  call chr_toReal(lnSplit(4),amp(1),iErr)  ! End amplitude
    if(nSplit > 4)  call chr_toReal(lnSplit(5),amp0,iErr)    ! Start amplitude
    if(nSplit > 5)  call chr_toInt(lnSplit(6),ird,iErr)      ! Ignored
    if(nSplit > 6)  call chr_toInt(lnSplit(7),imc,iErr)      ! Number of variations of the relative momentum deviation dp/p
    if(nSplit > 7)  call chr_toInt(lnSplit(8),niu(1),iErr)   ! Unknown
    if(nSplit > 8)  call chr_toInt(lnSplit(9),niu(2),iErr)   ! Unknown
    if(nSplit > 9)  call chr_toInt(lnSplit(10),numlcp,iErr)  ! CR: How often to write checkpointing files
    if(nSplit > 10) call chr_toInt(lnSplit(11),numlmax,iErr) ! CR: Maximum amount of turns; default is 1e6
    
    ! Default nnmul to numl
    nnuml = numl
    ! Defualt numlcp to 1000
    if(numlcp == 0) numlcp = 1000
    
    if(sixin_debug) then
      call sixin_echoVal("numl",numl,"TRAC",iLine)
      call sixin_echoVal("numlr",numlr,"TRAC",iLine)
      call sixin_echoVal("napx",napx,"TRAC",iLine)
      call sixin_echoVal("amp(1)",amp(1),"TRAC",iLine)
      call sixin_echoVal("amp0",amp0,"TRAC",iLine)
      call sixin_echoVal("ird",ird,"TRAC",iLine)
      call sixin_echoVal("imc",imc,"TRAC",iLine)
      call sixin_echoVal("niu(1)",niu(1),"TRAC",iLine)
      call sixin_echoVal("niu(2)",niu(2),"TRAC",iLine)
      call sixin_echoVal("numlcp",numlcp,"TRAC",iLine)
      call sixin_echoVal("numlmax",numlmax,"TRAC",iLine)
    end if
    if(iErr) return
    
  case(2)
    if(nSplit < 4) then
      write(lout,"(a,i0,a)") "PARAM> ERROR TRAC block line 2 should be at least 4 values, but ",nSplit," given."
      iErr = .true.
      return
    end if
    if(nSplit > 0) call chr_toInt(lnSplit(1),idz(1),iErr) ! Coupling on/off
    if(nSplit > 1) call chr_toInt(lnSplit(2),idz(2),iErr) ! Coupling on/off
    if(nSplit > 2) call chr_toInt(lnSplit(3),idfor,iErr)  ! Closed orbit and initial coordinates
    if(nSplit > 3) call chr_toInt(lnSplit(4),irew,iErr)   ! Disable rewind
    if(nSplit > 4) call chr_toInt(lnSplit(5),iclo6,iErr)  ! Calculate the 6D closed orbit
    if(idz(1) < 0 .or. idz(1) > 1) then
      write(lout,"(a,i0,a)") "PARAM> ERROR TRAC first value (idz(1)) can only be 0 or 1, but ",idz(1)," given."
      iErr = .true.
      return
    end if
    if(idz(2) < 0 .or. idz(2) > 1) then
      write(lout,"(a,i0,a)") "PARAM> ERROR TRAC second value (idz(2)) can only be 0 or 1, but ",idz(2)," given."
      iErr = .true.
      return
    end if
    if(idfor < 0 .or. idfor > 2) then
      write(lout,"(a,i0,a)") "PARAM> ERROR TRAC third value (idfor) can only be 0, 1 or 2, but ",idfor," given."
      iErr = .true.
      return
    end if
    if(irew < 0 .or. irew > 1) then
      write(lout,"(a,i0,a)") "PARAM> ERROR TRAC fourth value (irew) can only be 0 or 1, but ",irew," given."
      iErr = .true.
      return
    end if
    if(iclo6 /= 0 .and. iclo6 /= 1 .and. iclo6 /= 2 .and. iclo6 /= 5 .and. iclo6 /= 6) then
      write(lout,"(a,i0,a)") "PARAM> ERROR TRAC fifth value (iclo6) can only be 0, 1, 2, 5 or 6, but ",iclo6," given."
      iErr = .true.
      return
    end if
    
    if(sixin_debug) then
      call sixin_echoVal("idz(1)",idz(1),"TRAC",iLine)
      call sixin_echoVal("idz(2)",idz(2),"TRAC",iLine)
      call sixin_echoVal("idfor",idfor,"TRAC",iLine)
      call sixin_echoVal("irew",irew,"TRAC",iLine)
      call sixin_echoVal("iclo6",iclo6,"TRAC",iLine)
    end if
    if(iErr) return
    
    if(iclo6 == 5 .or. iclo6 == 6) then
      iclo6  = iclo6-4
      iclo6r = 1
    end if
    if(iclo6 == 2 .and. idfor == 0) idfor = 1
    if(iclo6 == 1 .or.  iclo6 == 2) nsix  = 0
    
  case(3)
    if(nSplit < 7) then
      write(lout,"(a,i0,a)") "PARAM> ERROR TRAC block line 3 should be at least 7 values, but ",nSplit," given."
      iErr = .true.
      return
    end if
    if(nSplit > 0) call chr_toInt(lnSplit(1),nde(1),iErr) ! Number of turns at flat bottom
    if(nSplit > 1) call chr_toInt(lnSplit(2),nde(2),iErr) ! Number of turns for the energy ramping
    if(nSplit > 2) call chr_toInt(lnSplit(3),nwr(1),iErr) ! Every nth turn coordinates will be written
    if(nSplit > 3) call chr_toInt(lnSplit(4),nwr(2),iErr) ! Every nth turn coordinates in the ramping region will be written
    if(nSplit > 4) call chr_toInt(lnSplit(5),nwr(3),iErr) ! Every nth turn at the flat top a write out of the coordinates
    if(nSplit > 5) call chr_toInt(lnSplit(6),nwr(4),iErr) ! Every nth turn coordinates are written to unit 6.
    if(nSplit > 6) call chr_toInt(lnSplit(7),ntwin,iErr)  ! Flag for calculated distance of phase space
    if(nSplit > 7) call chr_toInt(lnSplit(8),ibidu,iErr)  ! Switch to create or read binary dump
    if(nSplit > 8) call chr_toInt(lnSplit(9),iexact,iErr) ! Switch to enable exact solution of the equation of motion
    
    if(sixin_debug) then
      call sixin_echoVal("nde(1)",nde(1),"TRAC",iLine)
      call sixin_echoVal("nde(2)",nde(2),"TRAC",iLine)
      call sixin_echoVal("nwr(1)",nwr(1),"TRAC",iLine)
      call sixin_echoVal("nwr(2)",nwr(2),"TRAC",iLine)
      call sixin_echoVal("nwr(3)",nwr(3),"TRAC",iLine)
      call sixin_echoVal("nwr(4)",nwr(4),"TRAC",iLine)
      call sixin_echoVal("ntwin",ntwin,"TRAC",iLine)
      call sixin_echoVal("ibidu",ibidu,"TRAC",iLine)
      call sixin_echoVal("iexact",iexact,"TRAC",iLine)
    end if
    if(iErr) return
    
  case default
    write(lout,"(a,i0,a)") "PARAM> ERROR Unexpected line number ",iLine," in TRAC block."
    iErr = .true.
    return
  end select
  
end subroutine sixin_parseInputLineTRAC

end module sixtrack_input