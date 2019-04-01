! ================================================================================================ !
!  SixTrack Input Module
! ~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2018-05-18
!  Updated: 2019-04-01
! ================================================================================================ !
module sixtrack_input

  use floatPrecision
  use numerical_constants, only : zero, one, c1m3

  implicit none

  ! Record of encountered blocks
  character(len=:), allocatable, private, save :: sixin_cBlock(:) ! Name of block
  integer,          allocatable, private, save :: sixin_uBlock(:) ! Unit of block
  integer,          allocatable, private, save :: sixin_iBlock(:) ! Line of block
  logical,          allocatable, private, save :: sixin_lBlock(:) ! Block closed
  integer,                       private, save :: sixin_nBlock    ! Number of blocks

  integer,                       public,  save :: sixin_ncy2  = 0
  integer,                       public,  save :: sixin_icy   = 0

  ! Linear Optics Variables
  integer,                       private, save :: sixin_ilin0 = 1

  ! Synchrotron Oscillations
  real(kind=fPrec),              public,  save :: sixin_alc  = c1m3
  real(kind=fPrec),              public,  save :: sixin_harm = one
  real(kind=fPrec),              public,  save :: sixin_phag = zero
  real(kind=fPrec),              public,  save :: sixin_u0   = zero

  ! Multipole Coefficients
  integer,                       private, save :: sixin_im = 0

  ! RF-multipoles
  integer,                       private, save :: sixin_rfm = 0

  ! Beam-Beam Elements
  real(kind=fPrec),              public,  save :: sixin_emitNX = zero
  real(kind=fPrec),              public,  save :: sixin_emitNY = zero

  ! "Phase Trombone" Element
  integer,                       private, save :: sixin_imtr0 = 0

  ! Settings
  logical,                       private, save :: sixin_forcePartSummary = .false.
  logical,                       private, save :: sixin_forceWriteFort12 = .false.

  interface sixin_echoVal
    module procedure sixin_echoVal_int
    module procedure sixin_echoVal_real32
    module procedure sixin_echoVal_real64
    module procedure sixin_echoVal_real128
    module procedure sixin_echoVal_char
    module procedure sixin_echoVal_logical
  end interface sixin_echoVal

  private :: sixin_echoVal_int
  private :: sixin_echoVal_real32
  private :: sixin_echoVal_real64
  private :: sixin_echoVal_real128
  private :: sixin_echoVal_char
  private :: sixin_echoVal_logical

contains

! ================================================================================================ !
!  BLOCK PARSING RECORD
! ================================================================================================ !
subroutine sixin_checkBlock(blockName, blockUnit, blockOpened, blockClosed, blockLine, blockCount)

  use crcoall
  use mod_alloc
  use mod_settings

  character(len=*), intent(in)  :: blockName
  integer,          intent(in)  :: blockUnit
  logical,          intent(out) :: blockOpened
  logical,          intent(out) :: blockClosed
  integer,          intent(out) :: blockLine
  integer,          intent(out) :: blockCount

  integer i

  blockLine   = 0
  blockOpened = .false.
  blockClosed = .false.
  blockCount  = 1

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
        if(sixin_lBlock(i)) then
          ! Block already exists, but is closed.
          ! Count it, and continue.
          blockCount = blockCount + 1
        else
          ! Block exists and is opened,
          ! Just increment line number and return.
          blockOpened     = .true.
          blockClosed     = sixin_lBlock(i)
          blockLine       = sixin_iBlock(i) + 1
          sixin_iBlock(i) = blockLine
          return
        end if
      end if
    end do
    sixin_nBlock = sixin_nBlock + 1
  else
    ! If the array isn't allocated, it obviously doesn't contain anything
    sixin_nBlock = 1
  end if

  ! New block. Expand the arrays.
  call alloc(sixin_cBlock,4,sixin_nBlock,"    ", "sixin_cBlock")
  call alloc(sixin_uBlock,  sixin_nBlock,0,      "sixin_uBlock")
  call alloc(sixin_iBlock,  sixin_nBlock,0,      "sixin_iBlock")
  call alloc(sixin_lBlock,  sixin_nBlock,.false.,"sixin_lBlock")

  sixin_cBlock(sixin_nBlock)(1:4) = blockName(1:4)
  sixin_uBlock(sixin_nBlock)      = blockUnit
  sixin_iBlock(sixin_nBlock)      = 0
  sixin_lBlock(sixin_nBlock)      = .false.

  blockOpened = .true.

  if(st_debug) then
    write(lout,"(a)") "INPUT> Opened block '"//blockName//"'"
  end if

end subroutine sixin_checkBlock

subroutine sixin_closeBlock(blockName)

  use crcoall
  use mod_settings

  character(len=*), intent(in) :: blockName

  integer i

  do i=1,sixin_nBlock
    if(sixin_cBlock(i) == blockName .and. .not. sixin_lBlock(i)) then
      sixin_lBlock(i) = .true.
      if(st_debug) then
        write(lout,"(a)") "INPUT> Closed block '"//blockName//"'"
      end if
      return
    end if
  end do

  write(lout,"(a)") "INPUT> WARNING Could not close block '"//blockName//"'"

end subroutine sixin_closeBlock

subroutine sixin_blockReport

  use parpro
  use crcoall

  integer i

  write(lout,"(a)") str_divLine
  write(lout,"(a)") ""
  write(lout,"(a)") "    Finished parsing input file(s)."
  write(lout,"(a)") "    Parsed the following blocks:"
  do i=1,sixin_nBlock
    write(lout,"(a,i5,a,i0)") "     * "//sixin_cBlock(i)//" block "//&
      "with ",(sixin_iBlock(i)-1)," line(s) from fort.",sixin_uBlock(i)
  end do
  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine

end subroutine sixin_blockReport

subroutine sixin_echoVal_int(varName, varVal, blockName, lineNo)

  use crcoall
  use string_tools

  character(len=*), intent(in) :: varName
  integer,          intent(in) :: varVal
  character(len=*), intent(in) :: blockName
  integer,          intent(in) :: lineNo
  character(len=2) :: lineNm

  if(lineNo == -1) then
    lineNm = "PP"
  else if(lineNo < 10) then
    write(lineNm,"(i1,1x)") lineNo
  else
    write(lineNm,"(i2)") lineNo
  end if
  write(lout,"(a,i0)") "INPUT> DEBUG "//blockName//":"//lineNm//" "//chr_rpad(varName,10)//" =  ",varVal

end subroutine sixin_echoVal_int

subroutine sixin_echoVal_real32(varName, varVal, blockName, lineNo)

  use crcoall
  use string_tools

  character(len=*),  intent(in) :: varName
  real(kind=real32), intent(in) :: varVal
  character(len=*),  intent(in) :: blockName
  integer,           intent(in) :: lineNo
  character(len=2) :: lineNm

  if(lineNo == -1) then
    lineNm = "PP"
  else if(lineNo < 10) then
    write(lineNm,"(i1,1x)") lineNo
  else
    write(lineNm,"(i2)") lineNo
  end if
  write(lout,"(a,e13.6)") "INPUT> DEBUG "//blockName//":"//lineNm//" "//chr_rpad(varName,10)//" = ",varVal

end subroutine sixin_echoVal_real32

subroutine sixin_echoVal_real64(varName, varVal, blockName, lineNo)

  use crcoall
  use string_tools

  character(len=*),  intent(in) :: varName
  real(kind=real64), intent(in) :: varVal
  character(len=*),  intent(in) :: blockName
  integer,           intent(in) :: lineNo
  character(len=2) :: lineNm

  if(lineNo == -1) then
    lineNm = "PP"
  else if(lineNo < 10) then
    write(lineNm,"(i1,1x)") lineNo
  else
    write(lineNm,"(i2)") lineNo
  end if
  write(lout,"(a,e22.15)") "INPUT> DEBUG "//blockName//":"//lineNm//" "//chr_rpad(varName,10)//" = ",varVal

end subroutine sixin_echoVal_real64

subroutine sixin_echoVal_real128(varName, varVal, blockName, lineNo)

  use crcoall
  use string_tools

  character(len=*),   intent(in) :: varName
  real(kind=real128), intent(in) :: varVal
  character(len=*),   intent(in) :: blockName
  integer,            intent(in) :: lineNo
  character(len=2) :: lineNm

  if(lineNo == -1) then
    lineNm = "PP"
  else if(lineNo < 10) then
    write(lineNm,"(i1,1x)") lineNo
  else
    write(lineNm,"(i2)") lineNo
  end if
  write(lout,"(a,e41.34)") "INPUT> DEBUG "//blockName//":"//lineNm//" "//chr_rpad(varName,10)//" = ",varVal

end subroutine sixin_echoVal_real128

subroutine sixin_echoVal_char(varName, varVal, blockName, lineNo)

  use crcoall
  use string_tools

  character(len=*), intent(in) :: varName
  character(len=*), intent(in) :: varVal
  character(len=*), intent(in) :: blockName
  integer,          intent(in) :: lineNo
  character(len=2) :: lineNm

  if(lineNo == -1) then
    lineNm = "PP"
  else if(lineNo < 10) then
    write(lineNm,"(i1,1x)") lineNo
  else
    write(lineNm,"(i2)") lineNo
  end if
  write(lout,"(a)") "INPUT> DEBUG "//blockName//":"//lineNm//" "//chr_rpad(varName,10)//" = '"//varVal//"'"

end subroutine sixin_echoVal_char

subroutine sixin_echoVal_logical(varName, varVal, blockName, lineNo)

  use crcoall
  use string_tools

  character(len=*), intent(in) :: varName
  logical,          intent(in) :: varVal
  character(len=*), intent(in) :: blockName
  integer,          intent(in) :: lineNo
  character(len=2) :: lineNm

  if(lineNo == -1) then
    lineNm = "PP"
  else if(lineNo < 10) then
    write(lineNm,"(i1,1x)") lineNo
  else
    write(lineNm,"(i2)") lineNo
  end if
  if(varVal) then
    write(lout,"(a)") "INPUT> DEBUG "//blockName//":"//lineNm//" "//chr_rpad(varName,10)//" = True"
  else
    write(lout,"(a)") "INPUT> DEBUG "//blockName//":"//lineNm//" "//chr_rpad(varName,10)//" = False"
  end if

end subroutine sixin_echoVal_logical

! ================================================================================================ !
!  LINE PARSING ROUTINES
! ================================================================================================ !

! ================================================================================================ !
!  Parse Global Settings Block Line
! ================================================================================================ !
subroutine sixin_parseInputLineSETT(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, i
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(lnSplit(1))

  case("DEBUG")
    st_debug = .true.
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

  case("CRKILLSWITCH")
#ifdef CR
    if(nSplit < 2) then
      write(lout,"(a)") "INPUT> ERROR CRKILLSWITCH requires at least one turn number."
    end if
    st_killswitch = .true.
    write(lout,"(a)") "INPUT> C/R kill switch ENABLED"
    allocate(st_killturns(nSplit-1))
    do i=1,nSplit-1
      call chr_cast(lnSplit(i+1),st_killturns(i),iErr)
      write(lout,"(a,i0)") "INPUT>  * Will kill after turn ",st_killturns(i)
    end do
#else
    write(lout,"(a)") "INPUT> Ignoring CRKILLSWITCH flag. Not using CR version of SixTrack. "
#endif

  case("PRINT")
    st_print = .true.
    write(lout,"(a)") "INPUT> Printout of input parameters ENABLED"

  case("PRINT_DCUM")
    print_dcum = .true.
    write(lout,"(a)") "INPUT> Printout of dcum array ENABLED"

  case("PARTICLESUMMARY","PARTSUMMARY")
    if(nSplit > 1) then
      call chr_cast(lnSplit(2),st_partsum,iErr)
      sixin_forcePartSummary = .true.
    else
      st_partsum = .true.
      sixin_forcePartSummary = .true.
    end if
    if(st_partsum) then
      write(lout,"(a,i0)") "INPUT> Printing of particle summary is ENABLED"
    else
      write(lout,"(a,i0)") "INPUT> Printing of particle summary is DISABLED"
    end if

  case("WRITEFORT12")
    if(nSplit > 1) then
      call chr_cast(lnSplit(2),st_writefort12,iErr)
      sixin_forceWriteFort12 = .true.
    else
      st_writefort12 = .true.
      sixin_forceWriteFort12 = .true.
    end if
    if(st_writefort12) then
      write(lout,"(a,i0)") "INPUT> Writing of fort.12 after tracking is ENABLED"
    else
      write(lout,"(a,i0)") "INPUT> Writing of fort.12 after tracking is DISABLED"
    end if

  case("INITIALSTATE")
    if(nSplit /= 2 .and. nSplit /= 3) then
      write(lout,"(a,i0)") "INPUT> ERROR INITIALSTATE takes 1 or 2 values, got ",nSplit-1
      iErr = .true.
      return
    end if
    select case(lnSplit(2))
    case("binary")
      st_initialstate = 1
    case("text")
      st_initialstate = 2
    case default
      write(lout,"(a)") "INPUT> ERROR INITIALSTATE first value be either 'binary' or 'text', got '"//trim(lnSplit(2))//"'"
      iErr = .true.
      return
    end select
    if(nSplit == 3) then
      if(lnSplit(3) == "ions") then
        st_initialstate = st_initialstate + 2
      else
        write(lout,"(a)") "INPUT> ERROR INITIALSTATE second value must be 'ions', got '"//trim(lnSplit(3))//"'"
        iErr = .true.
        return
      end if
    end if
    select case(st_initialstate)
    case(1)
      write(lout,"(a,i0)") "INPUT> Particle initial state will be dumped as a binary file"
    case(2)
      write(lout,"(a,i0)") "INPUT> Particle initial state will be dumped as a text file"
    case(3)
      write(lout,"(a,i0)") "INPUT> Particle initial state will be dumped as a binary file with ion data included"
    case(4)
      write(lout,"(a,i0)") "INPUT> Particle initial state will be dumped as a text file with ion data included"
    end select

  case("FINALSTATE")
    if(nSplit /= 2 .and. nSplit /= 3) then
      write(lout,"(a,i0)") "INPUT> ERROR FINALSTATE takes 1 or 2 values, got ",nSplit-1
      iErr = .true.
      return
    end if
    select case(lnSplit(2))
    case("binary")
      st_finalstate = 1
    case("text")
      st_finalstate = 2
    case default
      write(lout,"(a)") "INPUT> ERROR FINALSTATE first value must be either 'binary' or 'text', got '"//trim(lnSplit(2))//"'"
      iErr = .true.
      return
    end select
    if(nSplit == 3) then
      if(lnSplit(3) == "ions") then
        st_finalstate = st_finalstate + 2
      else
        write(lout,"(a)") "INPUT> ERROR FINALSTATE second value must be 'ions', got '"//trim(lnSplit(3))//"'"
        iErr = .true.
        return
      end if
    end if
    select case(st_finalstate)
    case(1)
      write(lout,"(a,i0)") "INPUT> Particle final state will be dumped as a binary file"
    case(2)
      write(lout,"(a,i0)") "INPUT> Particle final state will be dumped as a text file"
    case(3)
      write(lout,"(a,i0)") "INPUT> Particle final state will be dumped as a binary file with ion data included"
    case(4)
      write(lout,"(a,i0)") "INPUT> Particle final state will be dumped as a text file with ion data included"
    end select

  case("QUIET")
    if(nSplit > 1) then
      call chr_cast(lnSplit(2),st_quiet,iErr)
    else
      st_quiet = 1
    end if
    write(lout,"(a,i0)") "INPUT> SixTrack Quiet level set to: ",st_quiet

  end select

end subroutine sixin_parseInputLineSETT

! ================================================================================================ !
!  Parse Displacement Block Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-05-xx
! ================================================================================================ !
subroutine sixin_parseInputLineDISP(inLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: elemName
  integer nSplit
  logical spErr

  integer i
  real(kind=fPrec) xpl0, xrms0, zpl0, zrms0

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "DISP> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  xpl0  = zero
  xrms0 = zero
  zpl0  = zero
  zrms0 = zero

  if(nSplit < 2) then
    write(lout,"(a,i0)") "DISP> ERROR Displacement of element line must have more than 1 values, got ",nSplit
    iErr = .true.
    return
  end if

  elemName = trim(lnSplit(1))
  if(len(elemName) > mNameLen) then
    write(lout,"(a,i0)") "DISP> ERROR Displacement of element name too long. Max length is ",mNameLen
    iErr = .true.
    return
  end if

  ! Save Values
  if(nSplit > 1) call chr_cast(lnSplit(2),xpl0, iErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),xrms0,iErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),zpl0, iErr)
  if(nSplit > 4) call chr_cast(lnSplit(5),zrms0,iErr)
  if(iErr) return

  if(st_debug) then
    call sixin_echoVal("xpl0", xpl0, "DISP",0)
    call sixin_echoVal("xrms0",xrms0,"DISP",0)
    call sixin_echoVal("zpl0", zpl0, "DISP",0)
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
        write(lout,"(a)") "DISP> INFO AC dipole disregarded, 0 length."
        kz(i) = 0
        ed(i) = zero
        ek(i) = zero
      end if
    end if

  end do

end subroutine sixin_parseInputLineDISP

! ================================================================================================ !
!  Parse Initial Coordinates Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-xx
! ================================================================================================ !
subroutine sixin_parseInputLineINIT(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common
  use mod_commons

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: expLine
  integer nSplit
  logical spErr

  integer i

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "INIT> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit < 1) then
    write(lout,"(a,i0,a)") "INIT> ERROR Block line ",iLine," did not receive any values."
    iErr = .true.
    return
  end if

  select case(iLine)
  case(1) ! Line One

    if(nSplit > 0) call chr_cast(lnSplit(1),itra,iErr) ! Number of particles
    if(nSplit > 1) call chr_cast(lnSplit(2),chi0,iErr) ! Starting phase of the initial coordinate
    if(nSplit > 2) call chr_cast(lnSplit(3),chid,iErr) ! Phase difference between particles
    if(nSplit > 3) call chr_cast(lnSplit(4),rat, iErr) ! Emittance ratio
    if(nSplit > 4) call chr_cast(lnSplit(5),iver,iErr) ! Vertical coordinates switch

    if(st_debug) then
      call sixin_echoVal("itra",itra,"INIT",iLine)
      call sixin_echoVal("chi0",chi0,"INIT",iLine)
      call sixin_echoVal("chid",chid,"INIT",iLine)
      call sixin_echoVal("rat", rat, "INIT",iLine)
      call sixin_echoVal("iver",iver,"INIT",iLine)
    end if

    if(iErr) return
    if(itra < 0 .or. itra > 2) then
      write(lout,"(a,i0,a)") "INIT> ERROR First value (itra) can only be 0, 1 or 2, but ",itra," given."
      iErr = .true.
      return
    end if

    if(iver < 0 .or. iver > 1) then
      write(lout,"(a,i0,a)") "INIT> ERROR Fifth value (iver) can only be 0 or 1, but ",iver," given."
      iErr = .true.
      return
    end if

  case(2)  ! x [mm] coordinate of particle 1
    call chr_cast(lnSplit(1),exz(1,1),iErr)
    if(st_debug) call sixin_echoVal("exz(1,1)",exz(1,1),"INIT",iLine)
    if(iErr) return

  case(3)  ! xp [mrad] coordinate of particle 1
    call chr_cast(lnSplit(1),exz(1,2),iErr)
    if(st_debug) call sixin_echoVal("exz(1,2)",exz(1,2),"INIT",iLine)
    if(iErr) return

  case(4)  ! y [mm] coordinate of particle 1
    call chr_cast(lnSplit(1),exz(1,3),iErr)
    if(st_debug) call sixin_echoVal("exz(1,3)",exz(1,3),"INIT",iLine)
    if(iErr) return

  case(5)  ! yp [mrad] coordinate of particle 1
    call chr_cast(lnSplit(1),exz(1,4),iErr)
    if(st_debug) call sixin_echoVal("exz(1,4)",exz(1,4),"INIT",iLine)
    if(iErr) return

  case(6)  ! Path length difference 1(sigma = s−v0*t) [mm] of particle 1
    call chr_cast(lnSplit(1),exz(1,5),iErr)
    if(st_debug) call sixin_echoVal("exz(1,5)",exz(1,5),"INIT",iLine)
    if(iErr) return

  case(7)  ! dp/p0 of particle 1
    call chr_cast(lnSplit(1),exz(1,6),iErr)
    if(st_debug) call sixin_echoVal("exz(1,6)",exz(1,6),"INIT",iLine)
    if(iErr) return

  case(8)  ! x [mm] coordinate of particle 2
    call chr_cast(lnSplit(1),exz(2,1),iErr)
    if(st_debug) call sixin_echoVal("exz(2,1)",exz(2,1),"INIT",iLine)
    if(iErr) return

  case(9)  ! xp [mrad] coordinate of particle 2
    call chr_cast(lnSplit(1),exz(2,2),iErr)
    if(st_debug) call sixin_echoVal("exz(2,2)",exz(2,2),"INIT",iLine)
    if(iErr) return

  case(10) ! y [mm] coordinate of particle 2
    call chr_cast(lnSplit(1),exz(2,3),iErr)
    if(st_debug) call sixin_echoVal("exz(2,3)",exz(2,3),"INIT",iLine)
    if(iErr) return

  case(11) ! yp [mrad] coordinate of particle 2
    call chr_cast(lnSplit(1),exz(2,4),iErr)
    if(st_debug) call sixin_echoVal("exz(2,4)",exz(2,4),"INIT",iLine)
    if(iErr) return

  case(12) ! Path length difference 1(sigma = s−v0*t) [mm] of particle 2
    call chr_cast(lnSplit(1),exz(2,5),iErr)
    if(st_debug) call sixin_echoVal("exz(2,5)",exz(2,5),"INIT",iLine)
    if(iErr) return

  case(13) ! dp/p0 of particle 2
    call chr_cast(lnSplit(1),exz(2,6),iErr)
    if(st_debug) call sixin_echoVal("exz(2,6)",exz(2,6),"INIT",iLine)
    if(iErr) return

  case(14) ! energy [MeV] of the reference particle
    call chr_cast(lnSplit(1),e0,iErr)
    if(st_debug) call sixin_echoVal("e0",e0,"INIT",iLine)
    if(iErr) return

  case(15) ! energy [MeV] of particle 1
    call chr_cast(lnSplit(1),ej(1),iErr)
    if(st_debug) call sixin_echoVal("ej(1)",ej(1),"INIT",iLine)
    if(iErr) return

  case(16) ! energy [MeV] of particle 2
    call chr_cast(lnSplit(1),ej(2),iErr)
    if(st_debug) call sixin_echoVal("ej(2)",ej(2),"INIT",iLine)
    if(iErr) return

  case default
    write(lout,"(a,i0)") "INIT> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end select

end subroutine sixin_parseInputLineINIT

! ================================================================================================ !
!  Parse Tracking Parameters Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-xx
! ================================================================================================ !
subroutine sixin_parseInputLineTRAC(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common
  use mod_commons
  use mod_common_da
  use mod_common_track

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: expLine
  integer nSplit, iDummy
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)
  case(1)

    if(nSplit < 7) then
      write(lout,"(a,i0,a)") "TRAC> ERROR Line 1 should be at least 7 values, but ",nSplit," given."
      iErr = .true.
      return
    end if

    if(nSplit > 0)  call chr_cast(lnSplit(1), numl,   iErr) ! Number of turns in the forward direction
    if(nSplit > 1)  call chr_cast(lnSplit(2), numlr,  iErr) ! Number of turns in the backward direction
    if(nSplit > 2)  call chr_cast(lnSplit(3), napx,   iErr) ! Number of amplitude variations (i.e. particle pairs)
    if(nSplit > 3)  call chr_cast(lnSplit(4), amp(1), iErr) ! End amplitude
    if(nSplit > 4)  call chr_cast(lnSplit(5), amp0,   iErr) ! Start amplitude
    if(nSplit > 5)  call chr_cast(lnSplit(6), ird,    iErr) ! Ignored
    if(nSplit > 6)  call chr_cast(lnSplit(7), imc,    iErr) ! Number of variations of the relative momentum deviation dp/p
    if(nSplit > 7)  call chr_cast(lnSplit(8), niu(1), iErr) ! Unknown
    if(nSplit > 8)  call chr_cast(lnSplit(9), niu(2), iErr) ! Unknown
    if(nSplit > 9)  call chr_cast(lnSplit(10),numlcp, iErr) ! CR: How often to write checkpointing files
    if(nSplit > 10) call chr_cast(lnSplit(11),numlmax,iErr) ! CR: Maximum amount of turns; default is 1e6

    ! Default nnuml to numl
    nnuml = numl

    if(st_debug) then
      call sixin_echoVal("numl",   numl,   "TRAC",iLine)
      call sixin_echoVal("numlr",  numlr,  "TRAC",iLine)
      call sixin_echoVal("napx",   napx,   "TRAC",iLine)
      call sixin_echoVal("amp(1)", amp(1), "TRAC",iLine)
      call sixin_echoVal("amp0",   amp0,   "TRAC",iLine)
      call sixin_echoVal("ird",    ird,    "TRAC",iLine)
      call sixin_echoVal("imc",    imc,    "TRAC",iLine)
      call sixin_echoVal("niu(1)", niu(1), "TRAC",iLine)
      call sixin_echoVal("niu(2)", niu(2), "TRAC",iLine)
      call sixin_echoVal("numlcp", numlcp, "TRAC",iLine)
      call sixin_echoVal("numlmax",numlmax,"TRAC",iLine)
    end if
    if(iErr) return

#ifndef STF
    if(napx > 32) then
      write(lout,"(a)") "TRAC> ERROR To run SixTrack with more than 32 particle pairs, it has to be compiled with the STF flag."
      iErr = .true.
      return
    end if
#endif

    if(napx*2 > npart) then
      call expand_arrays(nele, napx*2, nblz, nblo)
    end if

    if(napx > 32 .and. sixin_forcePartSummary .eqv. .false.) then
      write(lout,"(a)") "TRAC> NOTE More than 64 particles requested, switching off printing of particle summary."
      st_partsum = .false.
    end if

    if(napx > 32 .and. sixin_forceWriteFort12 .eqv. .false.) then
      write(lout,"(a)") "TRAC> NOTE More than 64 particles requested, switching off wriritng of fort.12."
      st_writefort12 = .false.
    end if

  case(2)

    if(nSplit < 4) then
      write(lout,"(a,i0,a)") "TRAC> ERROR Line 2 should be at least 4 values, but ",nSplit," given."
      iErr = .true.
      return
    end if

    if(nSplit > 0) call chr_cast(lnSplit(1),idz(1),iErr) ! Coupling on/off
    if(nSplit > 1) call chr_cast(lnSplit(2),idz(2),iErr) ! Coupling on/off
    if(nSplit > 2) call chr_cast(lnSplit(3),idfor, iErr) ! Closed orbit and initial coordinates
    if(nSplit > 3) call chr_cast(lnSplit(4),irew,  iErr) ! Disable rewind
    if(nSplit > 4) call chr_cast(lnSplit(5),iclo6, iErr) ! Calculate the 6D closed orbit

    if(idz(1) < 0 .or. idz(1) > 1) then
      write(lout,"(a,i0,a)") "TRAC> ERROR First value (idz(1)) can only be 0 or 1, but ",idz(1)," given."
      iErr = .true.
      return
    end if
    if(idz(2) < 0 .or. idz(2) > 1) then
      write(lout,"(a,i0,a)") "TRAC> ERROR Second value (idz(2)) can only be 0 or 1, but ",idz(2)," given."
      iErr = .true.
      return
    end if
    if(idfor < 0 .or. idfor > 2) then
      write(lout,"(a,i0,a)") "TRAC> ERROR Third value (idfor) can only be 0, 1, or 2, but ",idfor," given."
      iErr = .true.
      return
    end if
    if(irew < 0 .or. irew > 1) then
      write(lout,"(a,i0,a)") "TRAC> ERROR Fourth value (irew) can only be 0 or 1, but ",irew," given."
      iErr = .true.
      return
    end if
    if(iclo6 /= 0 .and. iclo6 /= 1 .and. iclo6 /= 2 .and. iclo6 /= 5 .and. iclo6 /= 6) then
      write(lout,"(a,i0,a)") "TRAC> ERROR Fifth value (iclo6) can only be 0, 1, 2, 5 or 6, but ",iclo6," given."
      iErr = .true.
      return
    end if

    if(st_debug) then
      call sixin_echoVal("idz(1)",idz(1),"TRAC",iLine)
      call sixin_echoVal("idz(2)",idz(2),"TRAC",iLine)
      call sixin_echoVal("idfor", idfor, "TRAC",iLine)
      call sixin_echoVal("irew",  irew,  "TRAC",iLine)
      call sixin_echoVal("iclo6", iclo6, "TRAC",iLine)
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
      write(lout,"(a,i0,a)") "TRAC> ERROR Line 3 should be at least 7 values, but ",nSplit," given."
      iErr = .true.
      return
    end if

    nwr(4) = 10000

    if(nSplit > 0) call chr_cast(lnSplit(1),nde(1),  iErr) ! Number of turns at flat bottom
    if(nSplit > 1) call chr_cast(lnSplit(2),nde(2),  iErr) ! Number of turns for the energy ramping
    if(nSplit > 2) call chr_cast(lnSplit(3),nwr(1),  iErr) ! Every nth turn coordinates will be written
    if(nSplit > 3) call chr_cast(lnSplit(4),nwr(2),  iErr) ! Every nth turn coordinates in the ramping region will be written
    if(nSplit > 4) call chr_cast(lnSplit(5),nwr(3),  iErr) ! Every nth turn at the flat top a write out of the coordinates
    if(nSplit > 5) call chr_cast(lnSplit(6),nwr(4),  iErr) ! Every nth turn coordinates are written to unit 6.
    if(nSplit > 6) call chr_cast(lnSplit(7),ntwin,   iErr) ! Flag for calculated distance of phase space
    if(nSplit > 7) call chr_cast(lnSplit(8),iDummy,  iErr) ! No longer in use. Formerly ibidu
    if(nSplit > 8) call chr_cast(lnSplit(9),iexact,  iErr) ! Switch to enable exact solution of the equation of motion
    if(nSplit > 9) call chr_cast(lnSplit(10),curveff,iErr) ! Switch to include curvatures effect on multipoles..

    if(st_debug) then
      call sixin_echoVal("nde(1)",nde(1),  "TRAC",iLine)
      call sixin_echoVal("nde(2)",nde(2),  "TRAC",iLine)
      call sixin_echoVal("nwr(1)",nwr(1),  "TRAC",iLine)
      call sixin_echoVal("nwr(2)",nwr(2),  "TRAC",iLine)
      call sixin_echoVal("nwr(3)",nwr(3),  "TRAC",iLine)
      call sixin_echoVal("nwr(4)",nwr(4),  "TRAC",iLine)
      call sixin_echoVal("ntwin", ntwin,   "TRAC",iLine)
      call sixin_echoVal("iexact",iexact,  "TRAC",iLine)
      call sixin_echoVal("curveff",curveff,"TRAC",iLine)
    end if
    if(iErr) return

  case default
    write(lout,"(a,i0)") "TRAC> ERROR Unexpected line number ",iLine
    iErr = .true.
    return
  end select

end subroutine sixin_parseInputLineTRAC

! ================================================================================================ !
!  Parse Differential Algebra Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-xx
! ================================================================================================ !
subroutine sixin_parseInputLineDIFF(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common
  use mod_common_da

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) ilm0(40)
  integer i, j1, j2, nSplit
  logical spErr

  do i=1,40
    ilm0(i) = " "
  end do

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "DIFF> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 1) then

    idial = 1
    numlr = 0
    napx  = 1
    imc   = 1

    if(nSplit > 0) call chr_cast(lnSplit(1),nord, iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),nvar, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),preda,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),nsix, iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),ncor, iErr)

    if(st_debug) then
      call sixin_echoVal("nord", nord, "DIFF",iLine)
      call sixin_echoVal("nvar", nvar, "DIFF",iLine)
      call sixin_echoVal("preda",preda,"DIFF",iLine)
      call sixin_echoVal("nsix", nsix, "DIFF",iLine)
      call sixin_echoVal("ncor", ncor, "DIFF",iLine)
    end if
    if(iErr) return

    if(nvar <= 4) ition = 0
    if(nord <= 0 .or. nvar <= 0) then
      write(lout,"(a)") "DIFF> ERROR Order and number of variables have to be larger than zero to "//&
        "calculate a differential algebra map."
      iErr = .true.
      return
    end if

  else

    if(nSplit /= ncor) then
      write(lout,"(2(a,i0))") "DIFF> ERROR Expected line > 1 to have ",ncor," elements, got ",nSplit
      iErr = .true.
      return
    end if
    do i=1,ncor
      ilm0(i) = chr_rpad(lnSplit(i),mNameLen)
    end do

  end if

  if(iclo6 == 1 .or. iclo6 == 2) nsix = 0
  if(nvar /= 6) then
    nsix  = 0
    iclo6 = 0
  end if
  if(nvar == 5) then
    idp    = 1
    ition  = 1
    hsy(1) = zero
  end if

  if(iLine == 1) then
    if(nsix /= 1) nsix = 0
    if(nord > nema) then
      write(lout,"(2(a,i0))") "DIFF> ERROR Maximum order of the one turn map is  ",nema,", got ",nord
      iErr = .true.
      return
    end if
    nvar2 = nvar
    return
  else
    if(ncor > mcor) then
      write(lout,"(2(a,i0))") "DIFF> ERROR Maximum number of extra parameters is  ",mcor,", got ",ncor
      iErr = .true.
      return
    end if
    if(ncor > 0) then
      OUTER: do j1=1,ncor
        INNER: do j2=1,il
          if(ilm0(j1) == bez(j2)) then
            if(el(j2) /= zero .or. kz(j2) > 10) then
              write(lout,"(a)") "DIFF> ERROR Only single kick elements allowed for map calculation"
              iErr = .true.
              return
            end if
            ipar(j1) = j2
            exit OUTER
          end if
        end do INNER
      end do OUTER
    else
      ncor = 0
      write(lout,"(a)") "DIFF> INFOR No extra parameters for the map specified"
    end if
    nvar = nvar2+ncor
  end if

end subroutine sixin_parseInputLineDIFF

! ================================================================================================ !
!  Parse Chromaticity Adjustment Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-xx
! ================================================================================================ !
subroutine sixin_parseInputLineCHRO(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common
  use mod_common_track

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen)      :: tmp_is(2)
  integer nSplit,i,ichrom0
  logical spErr

  save :: tmp_is,ichrom0

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "CHRO> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    ichrom0   = 0
    tmp_is(:) = " "

    if(nSplit > 0) tmp_is(1) = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),cro(1),   iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),ichrom0,  iErr)

    if(st_debug) then
      call sixin_echoVal("bez_is(1)",tmp_is(1),"CHRO",iLine)
      call sixin_echoVal("cro(1)",   cro(1),   "CHRO",iLine)
      call sixin_echoVal("ichrom0",  ichrom0,  "CHRO",iLine)
    end if
    if(iErr) return

  case(2)

    if(nSplit > 0) tmp_is(2) = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),cro(2),   iErr)

    if(st_debug) then
      call sixin_echoVal("bez_is(2)",tmp_is(2),"CHRO",iLine)
      call sixin_echoVal("cro(1)",   cro(2),   "CHRO",iLine)
    end if
    if(iErr) return

    do i=1,il
      if(tmp_is(1) == bez(i)) is(1) = i
      if(tmp_is(2) == bez(i)) is(2) = i
    end do
    if(ichrom0 >= 1 .and. ichrom0 <= 3) ichrom = ichrom0

    if(st_debug) then
      call sixin_echoVal("is(1)", is(1), "CHRO",iLine)
      call sixin_echoVal("is(2)", is(2), "CHRO",iLine)
      call sixin_echoVal("ichrom",ichrom,"CHRO",iLine)
    end if

  case default
    write(lout,"(a,i0)") "CHRO> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end select

end subroutine sixin_parseInputLineCHRO

! ================================================================================================ !
!  Parse Tune Adjustment Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-xx
! ================================================================================================ !
subroutine sixin_parseInputLineTUNE(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen)      :: tmp_iq(5)
  integer nSplit,i,nLines
  logical spErr

  save :: tmp_iq,nLines

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "TUNE> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    nLines    = 1
    tmp_iq(:) = " "

    if(nSplit > 0) tmp_iq(1) = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),qw0(1),iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),iqmod6,iErr)

    select case(iqmod6)
    case(1)
      iqmod  = 1
      iqmod6 = 0
    case(2)
      iqmod6 = 1
      iqmod  = 0
    case(3)
      iqmod  = 1
      iqmod6 = 1
    case default
      iqmod  = 1
      iqmod6 = 0
    end select

    if(st_debug) then
      call sixin_echoVal("tmp_iq(1)",tmp_iq(1),"TUNE",iLine)
      call sixin_echoVal("qw0(1)",   qw0(1),   "TUNE",iLine)
      call sixin_echoVal("iqmod",    iqmod,    "TUNE",iLine)
      call sixin_echoVal("iqmod6",   iqmod6,   "TUNE",iLine)
    end if
    if(iErr) return

  case(2)

    nLines = 2

    if(nSplit > 0) tmp_iq(2) = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),qw0(2),iErr)

    if(st_debug) then
      call sixin_echoVal("tmp_iq(2)",tmp_iq(2),"TUNE",iLine)
      call sixin_echoVal("qw0(2)",   qw0(2),   "TUNE",iLine)
    end if
    if(iErr) return

  case(3)

    nLines = 3

    if(nSplit > 0) tmp_iq(3) = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),qw0(3),iErr)

    if(st_debug) then
      call sixin_echoVal("tmp_iq(3)",tmp_iq(3),"TUNE",iLine)
      call sixin_echoVal("qw0(3)",   qw0(3),   "TUNE",iLine)
    end if
    if(iErr) return

  case(4)

    nLines = 4

    if(nSplit > 0) tmp_iq(4) = lnSplit(1)
    if(nSplit > 1) tmp_iq(5) = lnSplit(2)

    if(st_debug) then
      call sixin_echoVal("tmp_iq(4)",tmp_iq(4),"TUNE",iLine)
      call sixin_echoVal("tmp_iq(5)",tmp_iq(5),"TUNE",iLine)
    end if

  case(-1) ! Postprocessing

    if(nLines == 2 .or. nLines == 3) then
      if(abs(qw0(1)) > pieni .and. abs(qw0(2)) > pieni) then
        do i=1,il
          if(tmp_iq(1) == bez(i)) iq(1) = i
          if(tmp_iq(2) == bez(i)) iq(2) = i
        end do
        if(st_debug) then
          call sixin_echoVal("iq(1)",iq(1),"TUNE",-1)
          call sixin_echoVal("iq(2)",iq(2),"TUNE",-1)
        end if
      else
        write(lout,"(a)") "TUNE> Desired TUNE adjustment is zero. Block ignored."
        iqmod  = 0
        iqmod6 = 0
      end if
    else if(nLines == 4) then
      if(abs(qw0(1)) > pieni .and. abs(qw0(2)) > pieni .and. abs(qw0(3)) > pieni) then
        do i=1,il
          if(tmp_iq(1) == bez(1)) iq(1)  = i
          if(tmp_iq(2) == bez(1)) iq(2)  = i
          if(tmp_iq(3) == bez(1)) iq(3)  = i
          if(tmp_iq(4) == bez(1)) kpa(i) = 1
          if(tmp_iq(5) == bez(1)) kpa(i) = 2
        end do
        if(st_debug) then
          call sixin_echoVal("iq(1)",iq(1),"TUNE",-1)
          call sixin_echoVal("iq(2)",iq(2),"TUNE",-1)
          call sixin_echoVal("iq(3)",iq(3),"TUNE",-1)
        end if
      else
        write(lout,"(a)") "TUNE> Desired TUNE adjustment is zero. Block ignored."
        iqmod  = 0
        iqmod6 = 0
      endif
    else
      write(lout,"(a,i0)") "TUNE> ERROR Expected 2, 3 or 4 lines. Got ",nLines
      iErr = .true.
      return
    end if

  end select

end subroutine sixin_parseInputLineTUNE

! ================================================================================================ !
!  Parse Linear Optics Calculation Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-xx
! ================================================================================================ !
subroutine sixin_parseInputLineLINE(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) mode
  integer nSplit,i
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "LINE> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 1) then

    nlin = 0
    ilin = 1

    if(nSplit > 0) mode = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),nt,         iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),sixin_ilin0,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),ntco,       iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),eui,        iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),euii,       iErr)

    select case(mode)
    case("ELEMENT")
      iprint = 0
    case("BLOCK")
      iprint = 1
    case default
      write(lout,"(a)") "LINE> ERROR Valid modes are 'BLOCK' or 'ELEMENT'"
      iErr = .true.
    end select

    if(sixin_ilin0 == 1 .or. sixin_ilin0 == 2) then
      ilin = sixin_ilin0
    else
      write(lout,"(a)") "LINE> ERROR Only 1 (4D) and 2 (6D) are valid options for ilin."
      iErr = .true.
    end if

    if(st_debug) then
      call sixin_echoVal("mode",mode,"LINE",iLine)
      call sixin_echoVal("nt",  nt,  "LINE",iLine)
      call sixin_echoVal("ilin",ilin,"LINE",iLine)
      call sixin_echoVal("ntco",ntco,"LINE",iLine)
      call sixin_echoVal("eui", eui, "LINE",iLine)
      call sixin_echoVal("euii",euii,"LINE",iLine)
    end if
    if(iErr) return

  else

    do i=1,nSplit
      nlin = nlin + 1
      if(nlin > nele) then
        write(lout,"(2(a,i0))") "LINE> ERROR Too many elements for linear optics write out. Max is ",nele," got ",nlin
        iErr = .true.
        return
      end if
      bezl(nlin) = trim(lnSplit(i))
    end do

  end if

end subroutine sixin_parseInputLineLINE

! ================================================================================================ !
!  Parse Synchrotron Oscillations Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-xx
! ================================================================================================ !
subroutine sixin_parseInputLineSYNC(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common
  use mod_common_da,   only : nvar
  use mathlib_bouncer, only : cos_mb

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  real(kind=fPrec) cosy,halc,halc2,halc3,qigam,pmat,qbet
  integer          nSplit,i,ix
  logical          spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "SYNC> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    if(nSplit > 0) call chr_cast(lnSplit(1),sixin_harm,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),sixin_alc, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),sixin_u0,  iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),sixin_phag,iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),tlen,      iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),pma,       iErr)
    if(nSplit > 6) call chr_cast(lnSplit(7),ition,     iErr)
    if(nSplit > 7) call chr_cast(lnSplit(8),dppoff,    iErr)

    if(st_debug) then
      call sixin_echoVal("harm",  sixin_harm,"SYNC",iLine)
      call sixin_echoVal("alc",   sixin_alc, "SYNC",iLine)
      call sixin_echoVal("u0",    sixin_u0,  "SYNC",iLine)
      call sixin_echoVal("phag",  sixin_phag,"SYNC",iLine)
      call sixin_echoVal("tlen",  tlen,      "SYNC",iLine)
      call sixin_echoVal("pma",   pma,       "SYNC",iLine)
      call sixin_echoVal("ition", ition,     "SYNC",iLine)
      call sixin_echoVal("dppoff",dppoff,    "SYNC",iLine)
    end if
    if(iErr) return

  case(2)

    if(nSplit > 0) call chr_cast(lnSplit(1),dpscor,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),sigcor,iErr)

    if(st_debug) then
      call sixin_echoVal("dpscor",dpscor,"SYNC",iLine)
      call sixin_echoVal("sigcor",sigcor,"SYNC",iLine)
    end if
    if(iErr) return

    if(abs(pma-pmap) <= c1m1) pmat = pmap
    if(abs(pma-pmae) <= c1m1) pmat = pmae
    if(pmat /= pmap .and. pmat /= pmae) then
      write(lout,"(a)") "SYNC> WARNING Particle is neither proton nor electron"
    endif
    if(pma < pieni) then
      write(lout,"(a)") "SYNC> ERROR Kinetic energy of the particle is less than or equal to zero"
      iErr = .true.
      return
    end if
    crad = (crade*pmae)/pma
    if(abs(tlen) <= pieni) then
      write(lout,"(a)") "SYNC> ERROR Please include length of the machine."
      iErr = .true.
      return
    end if
    if(sixin_ncy2 == 0) then
      ncy = sixin_icy*mper
      idp = 1
      if(ncy == 0) then
        idp = 0
        write(lout,"(a)") "SYNC> No cavities specified."
      end if
      phas = sixin_phag*rad
      if(ncy /= 0) then
        hsy(1) = sixin_u0/real(ncy,fPrec)
      else
        hsy(1) = sixin_u0
      end if
      if(nvar == 5) then
        idp    = 1
        ition  = 1
        hsy(1) = zero
      end if
      halc   = sixin_harm*sixin_alc
      halc2  = sixin_harm/tlen
      hsy(3) = (two*pi)*halc2
      cosy   = cos_mb(phas)
      qigam  = (pma**2/e0)/e0
      qbet   = one-qigam
      halc3  = ((((((-one*(qigam-sixin_alc))*real(ition,fPrec))*sixin_harm)*sixin_u0)/e0)*cosy)/((two*pi)*qbet)
      qs     = sqrt(halc3)
      if(halc3 < zero) then
        write(lout,"(a)") "SYNC> ERROR Either your frequency is shifted by 180 degrees,"
        write(lout,"(a)") "SYNC>       then change the sign of ition in this block,"
        write(lout,"(a)") "SYNC>       or your alpha-p is wrongly introducd."
        iErr = .true.
        return
      end if
    else
      idp = 1
      ncy = 0
      do i=1,mper*mbloz
        ix = ic(i)
        if(ix > nblo) then
          ix = ix-nblo
          if(abs(kz(ix)) == 12) ncy = ncy+1
        end if
      end do
      do i=1,il
        if(abs(kz(i)) == 12) then
          hsyc(i) = ((two*pi)*ek(i))/tlen
          if(nvar == 5) then
            ition = 1
            ed(i) = zero
          end if
        end if
      end do
    endif

    if(st_debug) then
      call sixin_echoVal("phas",phas,"SYNC",iLine)
      call sixin_echoVal("idp", idp, "SYNC",iLine)
      call sixin_echoVal("ncy", ncy, "SYNC",iLine)
    end if
    if(iErr) return

  case default
    write(lout,"(a,i0)") "SYNC> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end select

end subroutine sixin_parseInputLineSYNC

! ================================================================================================ !
!  Parse Multipole Coefficient Line for KZ=11
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-18
! ================================================================================================ !
subroutine sixin_parseInputLineMULT(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) imn
  real(kind=fPrec) ak0d,akad,bk0d,bkad,r0,r0a
  integer          nSplit,i,nmul,iil
  logical          spErr

  real(kind=fPrec) :: benki = zero

  save nmul,iil,r0,r0a,benki

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "MULT> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 1) then

    if(nSplit > 0) imn = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),r0,   iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),benki,iErr)

    iil      = -1
    nmul     = 1
    r0a      = one
    sixin_im = sixin_im + 1

    benkc(sixin_im) = benki
    r00(sixin_im)   = r0

    do i=1,il
      if(imn == bez(i)) then
        irm(i) = sixin_im
        iil    = i
        exit
      end if
    end do

    if(iil == -1) then
      write(lout,"(a)") "MULT> ERROR Single element '"//trim(imn)//"' not found in single element list."
      iErr = .true.
      return
    end if

    if(st_debug) then
      call sixin_echoVal("imn",  imn,  "MULT",iLine)
      call sixin_echoVal("r0",   r0,   "MULT",iLine)
      call sixin_echoVal("benki",benki,"MULT",iLine)
    end if
    if(iErr) return

  else

    bk0d = zero
    bkad = zero
    ak0d = zero
    akad = zero

    if(nSplit > 0) call chr_cast(lnSplit(1),bk0d,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),bkad,iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),ak0d,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),akad,iErr)

    if(st_debug) then
      call sixin_echoVal("bk0d",bk0d,"MULT",iLine)
      call sixin_echoVal("bkad",bkad,"MULT",iLine)
      call sixin_echoVal("ak0d",ak0d,"MULT",iLine)
      call sixin_echoVal("akad",akad,"MULT",iLine)
      call sixin_echoVal("r0a", r0a, "MULT",iLine)
      call sixin_echoVal("nmul",nmul,"MULT",iLine)
    end if
    if(iErr) return

    ! Set nmu for the current single element (j)
    ! to the currently highest multipole seen (i)
    ! Changed so also 0 is considered to be a mutipole, since it might be changed later by dynk

    nmu(iil) = nmul

    bk0(sixin_im,nmul) = (benki*bk0d)/r0a
    ak0(sixin_im,nmul) = (benki*ak0d)/r0a
    bka(sixin_im,nmul) = (benki*bkad)/r0a
    aka(sixin_im,nmul) = (benki*akad)/r0a
    nmul = nmul + 1
    r0a  = r0a*r0
    if(nmul > mmul+1) then
      write(lout,"(a,i0)") "MULT> ERROR The order of multipoles is too large. Maximum is ",mmul
      iErr = .true.
      return
    end if

  end if

end subroutine sixin_parseInputLineMULT

! ================================================================================================ !
!  Parse RF Multipoles
!  Last modified: 2018-12-31
! ================================================================================================ !
subroutine sixin_parseInputLineRFMU(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) imn
  real(kind=fPrec) namp0,nphase0,samp0,sphase0, freq0
  integer          nSplit,i,nmul,iil
  logical          spErr

  save nmul,iil

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "RFMU> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 1) then

    if(nSplit > 0) imn = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),freq0,   iErr)


    iil      = -1
    nmul     = 1
    sixin_rfm = sixin_rfm + 1
    freq_rfm(sixin_rfm) = freq0

    do i=1,il
      if(imn == bez(i)) then
        irm_rf(i) = sixin_rfm
        iil    = i
        exit
      end if
    end do

    if(iil == -1) then
      write(lout,"(a)") "RFMU> ERROR Single element '"//trim(imn)//"' not found in single element list."
      iErr = .true.
      return
    end if

    if(st_debug) then
      call sixin_echoVal("imn",imn,"RFMU",iLine)
    end if

    if(iErr) return

  else

    namp0   = zero
    nphase0 = zero
    samp0   = zero
    sphase0 = zero
    if(nSplit > 0) call chr_cast(lnSplit(1),namp0,  iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),nphase0,iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),samp0,  iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),sphase0,iErr)
    if(st_debug) then
      call sixin_echoVal("namp0",  namp0,  "RFMU",iLine)
      call sixin_echoVal("nphase0",nphase0,"RFMU",iLine)
      call sixin_echoVal("samp0",  samp0,  "RFMU",iLine)
      call sixin_echoVal("sphase0",sphase0,"RFMU",iLine)
    end if
    if(iErr) return

    ! Set nmu for the current single element (j)
    ! to the currently highest multipole seen (i)
    ! Changed so also 0 is considered to be a mutipole, since it might be changed later by dynk

    nmu_rf(sixin_rfm)        = nmul
    norrfamp(sixin_rfm,nmul) = namp0
    norrfph(sixin_rfm,nmul)  = nphase0
    skrfamp(sixin_rfm,nmul)  = samp0
    skrfph(sixin_rfm,nmul)   = sphase0
    nmul = nmul + 1

    if(nmul > mmul+1) then
      write(lout,"(a,i0)") "RFMU> ERROR The order of multipoles is too large. Maximum is ",mmul
      iErr = .true.
      return
    end if

  end if

end subroutine sixin_parseInputLineRFMU

! ================================================================================================ !
!  Parse Sub-Resonance Calculation Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-20
!  Note: This block is not covered by any tests
! ================================================================================================ !
subroutine sixin_parseInputLineSUBR(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer          nSplit
  logical          spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "SUBR> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 1) then

    if(nSplit > 0) call chr_cast(lnSplit(1),nta, iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),nte, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),qxt, iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),qzt, iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),tam1,iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),tam2,iErr)
    if(nSplit > 6) call chr_cast(lnSplit(7),ipt, iErr)
    if(nSplit > 7) call chr_cast(lnSplit(8),totl,iErr)

    if(st_debug) then
      call sixin_echoVal("nta", nta, "SUBR",iLine)
      call sixin_echoVal("nte", nte, "SUBR",iLine)
      call sixin_echoVal("qxt", qxt, "SUBR",iLine)
      call sixin_echoVal("qzt", qzt, "SUBR",iLine)
      call sixin_echoVal("tam1",tam1,"SUBR",iLine)
      call sixin_echoVal("tam2",tam2,"SUBR",iLine)
      call sixin_echoVal("ipt", ipt, "SUBR",iLine)
      call sixin_echoVal("totl",totl,"SUBR",iLine)
    end if
    if(iErr) return

    if(nta < 2 .or. nte < nta .or. nte > 9) then
      write(lout,"(a)") "SUBR> ERROR Chosen orders of resonances can not be calculated."
      iErr = .true.
      return
    end if

  else

    write(lout,"(a,i0)") "SUBR> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end if

end subroutine sixin_parseInputLineSUBR

! ================================================================================================ !
!  Parse Organisation of Random Numbers Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-22
!  Note: This block is not covered by any tests
! ================================================================================================ !
subroutine sixin_parseInputLineORGA(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) elemOne
  integer nSplit, i, j0, j1, imo
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "ORGA> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  iorg = iorg + 1
  elemOne = " "
  if(nSplit > 0) elemOne      = trim(lnSplit(1))
  if(nSplit > 1) bezr(2,iorg) = trim(lnSplit(2))
  if(nSplit > 2) bezr(3,iorg) = trim(lnSplit(3))

  if(iorg == 1) then
    write(lout,"(a)") "ORGA>               |"//&
    " Own Random Num   |"                    //&
    " Same Random Numbers                 |" //&
    " Same Multipole Coefficients         |"
    write(lout,"(a)") "ORGA>               +"//&
    "------------------+"                    //&
    "------------------+------------------+" //&
    "------------------+------------------+"
  end if

  if(elemOne /= "MULT" .and. elemOne /= " ") then
    if(bezr(2,iorg) == " ") then
      write(lout,"(a,i4,a)") "ORGA> Elements ",iLine," |"//&
        " "//elemOne(1:16)//" |"                         //&
        "                  |"                            //&
        "                  |"                            //&
        "                  |"                            //&
        "                  |"
     else
      write(lout,"(a,i4,a)") "ORGA> Elements ",iLine," |"//&
        "                  |"                            //&
        " "//elemOne(1:16)//" |"                         //&
        " "//bezr(2,iorg)(1:16)//" |"              //&
        "                  |"                            //&
        "                  |"
     end if
  end if
  if(elemOne /= "MULT") bezr(1,iorg) = elemOne
  if(elemOne == "MULT" .and. bezr(2,iorg) /= " " .and. bezr(3,iorg) /= " ") then
    write(lout,"(a,i4,a)") "ORGA> Elements ",iLine," |"//&
      "                  |"                            //&
      "                  |"                            //&
      "                  |"                            //&
      " "//bezr(2,iorg)(1:16)//" |"              //&
      " "//bezr(3,iorg)(1:16)//" |"
    sixin_im = sixin_im + 1

    j0 = 0
    j1 = 0

    do i=1,il
      if(bez(i) == bezr(2,iorg)) j1 = i
      if(bez(i) == bezr(3,iorg)) j0 = i
    end do
    if(j0 == 0 .or. j1 == 0 .or. kz(j0) == 11 .or. kz(j1) == 11) then
      write(lout,"(a)") "ORGA> ERROR Multipole coefficients cannot be set equal."
      iErr = .true.
      return
    end if

    irm(j0)   = sixin_im
    benkc(j0) = benkc(j1)
    r00(j0)   = r00(j1)
    imo       = irm(j1)
    nmu(j0)   = nmu(j1)

    do i=1,nmu(j0)
      bk0(sixin_im,i)=bk0(imo,i)
      bka(sixin_im,i)=bka(imo,i)
      ak0(sixin_im,i)=ak0(imo,i)
      aka(sixin_im,i)=aka(imo,i)
    end do

  end if

end subroutine sixin_parseInputLineORGA

! ================================================================================================ !
!  Parse Iteration Errors Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-22
! ================================================================================================ !
subroutine sixin_parseInputLineITER(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "ITER> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    if(nSplit > 0) call chr_cast(lnSplit(1),itco,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),dma, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),dmap,iErr)

    if(st_debug) then
      call sixin_echoVal("itco",itco,"ITER",iLine)
      call sixin_echoVal("dma", dma, "ITER",iLine)
      call sixin_echoVal("dmap",dmap,"ITER",iLine)
    end if
    if(iErr) return

  case(2)

    if(nSplit > 0) call chr_cast(lnSplit(1),itqv,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),dkq, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),dqq, iErr)

    if(st_debug) then
      call sixin_echoVal("itqv",itqv,"ITER",iLine)
      call sixin_echoVal("dkq", dkq, "ITER",iLine)
      call sixin_echoVal("dqq", dqq,"ITER",iLine)
    end if
    if(iErr) return

  case(3)

    if(nSplit > 0) call chr_cast(lnSplit(1),itcro,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),dsm0, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),dech, iErr)

    if(st_debug) then
      call sixin_echoVal("itcro",itcro,"ITER",iLine)
      call sixin_echoVal("dsm0", dsm0, "ITER",iLine)
      call sixin_echoVal("dech", dech, "ITER",iLine)
    end if
    if(iErr) return

  case(4)

    if(nSplit > 0) call chr_cast(lnSplit(1),de0,    iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),ded,    iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),dsi,    iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),aper(1),iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),aper(2),iErr)

    if(st_debug) then
      call sixin_echoVal("de0",    de0,    "ITER",iLine)
      call sixin_echoVal("ded",    ded,    "ITER",iLine)
      call sixin_echoVal("dsi",    dsi,    "ITER",iLine)
      call sixin_echoVal("aper(1)",aper(1),"ITER",iLine)
      call sixin_echoVal("aper(2)",aper(2),"ITER",iLine)
    end if
    if(iErr) return

  case default
    write(lout,"(a,i0)") "ITER> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end select

end subroutine sixin_parseInputLineITER

! ================================================================================================ !
!  Parse Orbit Correction Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-23
! ================================================================================================ !
subroutine sixin_parseInputLineORBI(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, i, iElem
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "ORBI> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 1) then

    if(nSplit > 0) call chr_cast(lnSplit(1),sigma0(1), iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),sigma0(2), iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),ncorru,    iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),ncorrep,   iErr)

    if(st_debug) then
      call sixin_echoVal("sigmax", sigma0(1), "ORBI",iLine)
      call sixin_echoVal("sigmay", sigma0(2), "ORBI",iLine)
      call sixin_echoVal("ncorru", ncorru,    "ORBI",iLine)
      call sixin_echoVal("ncorrep",ncorrep,   "ORBI",iLine)
    end if
    if(iErr) return

  else

    if(nSplit /= 2) then
      write(lout,"(a,i0)") "ORBI> ERROR Expected 2 parameters for line > 2, got ",nSplit
      write(lout,"(a)")    "ORBI>       If your file has for instance HCOR=name, replace the = with a space."
      write(lout,"(a)")    "ORBI>       Name/value pairs with = is no longer supported in SixTrack for consistency between blocks."
      iErr = .true.
      return
    end if

    ! Look up element
    iElem = -1
    do i=1,il
      if(lnSplit(2) == bez(i)) then
        iElem = i
        exit
      end if
    end do

    if(iElem == -1) then
      write(lout,"(a)") "ORBI> ERROR Unknown element name '"//lnSplit(2)//"'"
      iErr = .true.
      return
    end if

    select case(lnSplit(1))

    case("HCOR")
      if(kp(iElem) == -4 .or. kp(iElem) == 3 .or. kp(iElem) == -3) then
        write(lout,"(a)") "ORBI> ERROR An element for closed orbit correction can be only either a horizontal monitor "
        write(lout,"(a)") "ORBI>       or a vertical monitor or a horizontal corrector or a vertical corrector."
        iErr = .true.
        return
      end if
      if(kz(iElem) /= 1 .and. kz(iElem) /= 11) then
        write(lout,"(a)") "ORBI> ERROR For closed orbit correctors only dipoles of legth zero or multipole lenses are allowed."
        iErr = .true.
        return
      end if
      kp(iElem) = 4

    case("VCOR")
      if(kp(iElem) == 4 .or. kp(iElem) == 3 .or. kp(iElem) == -3) then
        write(lout,"(a)") "ORBI> ERROR An element for closed orbit correction can be only either a horizontal monitor "
        write(lout,"(a)") "ORBI>       or a vertical monitor or a horizontal corrector or a vertical corrector."
        iErr = .true.
        return
      end if
      if(kz(iElem) /= -1 .and. kz(iElem) /= 11) then
        write(lout,"(a)") "ORBI> ERROR For closed orbit correctors only dipoles of legth zero or multipole lenses are allowed."
        iErr = .true.
        return
      end if
      kp(iElem) = -4

    case("HMON")
      if(kp(iElem) == 4 .or. kp(iElem) == -4 .or. kp(iElem) == -3) then
        write(lout,"(a)") "ORBI> ERROR An element for closed orbit correction can be only either a horizontal monitor "
        write(lout,"(a)") "ORBI>       or a vertical monitor or a horizontal corrector or a vertical corrector."
        iErr = .true.
        return
      end if
      kp(iElem) = 3

    case("VMON")
      if(kp(iElem) == 4 .or. kp(iElem) == -4 .or. kp(iElem) == 3) then
        write(lout,"(a)") "ORBI> ERROR An element for closed orbit correction can be only either a horizontal monitor "
        write(lout,"(a)") "ORBI>       or a vertical monitor or a horizontal corrector or a vertical corrector."
        iErr = .true.
        return
      end if
      kp(iElem) = -3

    case default
      write(lout,"(a)") "ORBI> ERROR Only correctors with the keywords HCOR/VCOR"
      write(lout,"(a)") "ORBI>       or monitors with the keywords HMON/VMON are allowed."
      iErr = .true.
      return

    end select

  end if

end subroutine sixin_parseInputLineORBI

! ================================================================================================ !
!  Parse Combination of Elements Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-22
! ================================================================================================ !
subroutine sixin_parseInputLineCOMB(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) elemName, elemComb(20)
  integer nSplit, nComb, i, j, ii, ico
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "COMB> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine > ncom) then
    write(lout,"(a,i0)") "COMB> ERROR Maximum number of combinations is ",ncom
    iErr = .true.
    return
  end if

  if(nSplit < 3 .or. nSplit > 41 .or. modulo(nSplit,2) /= 1) then
    write(lout,"(a,i0)") "COMB> ERROR Expected one element name + 1 to 20 pairs, got ",nSplit
    iErr = .true.
    return
  end if

  icoe        = iLine
  elemName    = trim(lnSplit(1))
  nComb       = (nSplit-1)/2
  elemComb(:) = " "
  do i=1,nComb
    call chr_cast(lnSplit(2*i),ratio(icoe,i),iErr)
    elemComb(i) = trim(lnSplit(2*i+1))
  end do
  if(st_debug) then
    call sixin_echoVal("element",elemName,"COMB",iLine)
  end if

  do i=1,il
    if(elemName == bez(i)) then
      kp(i)        = 5
      icomb0(icoe) = i
      ratioe(i)    = one
    end if
    do j=1,nComb
      if(elemComb(j) == bez(i)) then
        icomb(icoe,j) = i
        ratioe(i)     = ratio(icoe,j)
      end if
    end do
  end do

  ii = icomb0(icoe)
  if(ii == 0) then
    write(lout,"(a)") "COMB> ERROR Element '"//trim(elemName)//"' not found."
    iErr = .true.
    return
  end if

  do i=1,nComb
    ico = icomb(icoe,i)
    if(ico == ii) then
      write(lout,"(a)") "COMB> ERROR You cannot combine an element with itself."
      iErr = .true.
      return
    end if
    if(ico == 0) cycle
    write(lout,"(a,e13.6)") "COMB> "//bez(ii)(1:20)//" : "//bez(ico)(1:20)//" : ",ratio(icoe,i)
    iratioe(ico) = ii
    if(el(ii) <= pieni) then
      if(el(ico) <= pieni) then
        ed(ico) = ed(ii)*ratio(icoe,i)
      else
        ek(ico) = ed(ii)*ratio(icoe,i)
      end if
    else
      if(el(ico) <= pieni) then
        ed(ico) = ek(ii)*ratio(icoe,i)
      else
        ek(ico) = ek(ii)*ratio(icoe,i)
      end if
    end if
  end do

end subroutine sixin_parseInputLineCOMB

! ================================================================================================ !
!  Parse Resonance Compensation Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-23
! ================================================================================================ !
subroutine sixin_parseInputLineRESO(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) name(10)
  integer nSplit, j, k
  logical spErr

  save :: name

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "RESO> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    if(nSplit > 0) call chr_cast(lnSplit(1),nre,   iErr)
    if(nre /= 0) then
      if(nSplit > 1) call chr_cast(lnSplit(2),npp,   iErr)
      if(nSplit > 2) call chr_cast(lnSplit(3),nrr(1),iErr)
      if(nSplit > 3) call chr_cast(lnSplit(4),nrr(2),iErr)
      if(nSplit > 4) call chr_cast(lnSplit(5),nrr(3),iErr)
      if(nSplit > 5) call chr_cast(lnSplit(6),ipr(1),iErr)
      if(nSplit > 6) call chr_cast(lnSplit(7),ipr(2),iErr)
      if(nSplit > 7) call chr_cast(lnSplit(8),ipr(3),iErr)
    end if

    if(st_debug) then
      call sixin_echoVal("nre",   nre,   "RESO",iLine)
      call sixin_echoVal("npp",   npp,   "RESO",iLine)
      call sixin_echoVal("nrr(1)",nrr(1),"RESO",iLine)
      call sixin_echoVal("nrr(2)",nrr(2),"RESO",iLine)
      call sixin_echoVal("nrr(3)",nrr(3),"RESO",iLine)
      call sixin_echoVal("ipr(1)",ipr(1),"RESO",iLine)
      call sixin_echoVal("ipr(2)",ipr(2),"RESO",iLine)
      call sixin_echoVal("ipr(3)",ipr(3),"RESO",iLine)
    end if
    if(iErr) return

    if(nre /= 0 .and. (npp < 2 .or. npp > nrco)) then
      write(lout,"(a,i0)") "RESO> ERROR Order of compensation can not be larger than ",nrco
      iErr = .true.
      return
    end if
    if(nre < 0 .or. nre > 3) then
      write(lout,"(a)") "RESO> ERROR Only up to 3 resonances can be compensated."
      iErr = .true.
      return
    end if
    if(abs(nrr(1)) > npp .or. abs(nrr(2)) > npp .or. abs(nrr(3)) > npp) then
      write(lout,"(a)") "RESO> ERROR Resonance type is out of the range of the resonance order."
      iErr = .true.
      return
    end if

  case(2)

    if(nSplit > 0) call chr_cast(lnSplit(1),nur,  iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),nu(1),iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),nu(2),iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),nu(3),iErr)

    if(st_debug) then
      call sixin_echoVal("nur",  nur,  "RESO",iLine)
      call sixin_echoVal("nu(1)",nu(1),"RESO",iLine)
      call sixin_echoVal("nu(2)",nu(2),"RESO",iLine)
      call sixin_echoVal("nu(3)",nu(3),"RESO",iLine)
    end if
    if(iErr) return

    if(nur < 0 .or. nur > 3) then
      write(lout,"(a)") "RESO> ERROR Only up to 3 sub-resonances can be compensated."
      iErr = .true.
      return
    end if
    if(nu(1) > 9 .or. nu(2) > 9 .or. nu(3) > 9 .or. nu(1) < 0 .or. nu(2) < 0 .or. nu(3) < 0) then
      write(lout,"(a)") "RESO> ERROR The multipole order for the sub-resonance compensation should not exceed 9."
      iErr = .true.
      return
    end if

  case(3)

    if(nSplit > 0) call chr_cast(lnSplit(1),totl,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),qxt, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),qzt, iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),tam1,iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),tam2,iErr)

    if(st_debug) then
      call sixin_echoVal("totl",totl,"RESO",iLine)
      call sixin_echoVal("qxt", qxt, "RESO",iLine)
      call sixin_echoVal("qzt", qzt, "RESO",iLine)
      call sixin_echoVal("tam1",tam1,"RESO",iLine)
      call sixin_echoVal("tam2",tam2,"RESO",iLine)
    end if
    if(iErr) return

  case(4)

    name(:) = " "

    if(nSplit > 0) name(1) = trim(lnSplit(1))
    if(nSplit > 1) name(2) = trim(lnSplit(2))
    if(nSplit > 2) name(3) = trim(lnSplit(3))
    if(nSplit > 3) name(4) = trim(lnSplit(4))
    if(nSplit > 4) name(5) = trim(lnSplit(5))
    if(nSplit > 5) name(6) = trim(lnSplit(6))

    if(st_debug) then
      call sixin_echoVal("namel",name(1),"RESO",iLine)
      call sixin_echoVal("name2",name(2),"RESO",iLine)
      call sixin_echoVal("name3",name(3),"RESO",iLine)
      call sixin_echoVal("name4",name(4),"RESO",iLine)
      call sixin_echoVal("name5",name(5),"RESO",iLine)
      call sixin_echoVal("name6",name(6),"RESO",iLine)
    end if
    if(iErr) return

  case(5)

    if(nSplit > 0) call chr_cast(lnSplit(1),nch,iErr)
    if(nch /= 0) then
      if(nSplit > 1) name(7) = trim(lnSplit(2))
      if(nSplit > 2) name(8) = trim(lnSplit(3))
    end if

    if(st_debug) then
      call sixin_echoVal("nch",  nch,    "RESO",iLine)
      call sixin_echoVal("name7",name(7),"RESO",iLine)
      call sixin_echoVal("name8",name(8),"RESO",iLine)
    end if
    if(iErr) return

  case(6)

    if(nSplit > 0) call chr_cast(lnSplit(1),nqc,iErr)
    if(nqc /= 0) then
      if(nSplit > 1) name(9)  = trim(lnSplit(2))
      if(nSplit > 2) name(10) = trim(lnSplit(3))
      if(nSplit > 3) call chr_cast(lnSplit(4),qw0(1),iErr)
      if(nSplit > 4) call chr_cast(lnSplit(5),qw0(2),iErr)
    end if

    if(st_debug) then
      call sixin_echoVal("nqc",   nqc,     "RESO",iLine)
      call sixin_echoVal("name9", name(9), "RESO",iLine)
      call sixin_echoVal("name10",name(10),"RESO",iLine)
      call sixin_echoVal("qw0(1)",qw0(1),  "RESO",iLine)
      call sixin_echoVal("qw0(2)",qw0(2),  "RESO",iLine)
    end if
    if(iErr) return

    outer: do k=1,10
      inner: do j=1,il
        if(name(k) /= bez(j)) cycle inner
        ire(k) = j
        if(nre == 1 .and. k < 3 .and. abs(kz(j)) /= npp) then
          write(lout,"(a)") "RESO> ERROR With the specified elements the resonance cannot be compensated."
          write(lout,"(a)") "RESO> ERROR Resonance order and element type # must be the same."
          iErr = .true.
          return
        end if
        if(nre == 2 .and. k < 5 .and. abs(kz(j)) /= npp) then
          write(lout,"(a)") "RESO> ERROR With the specified elements the resonance cannot be compensated."
          write(lout,"(a)") "RESO> ERROR Resonance order and element type # must be the same."
          iErr = .true.
          return
        end if
        if(nre == 3 .and. k < 7 .and. abs(kz(j)) /= npp) then
          write(lout,"(a)") "RESO> ERROR With the specified elements the resonance cannot be compensated."
          write(lout,"(a)") "RESO> ERROR Resonance order and element type # must be the same."
          iErr = .true.
          return
        end if
        if(nch == 1 .and. (k == 7 .or. k == 8)  .and. kz(j) /= 3) then
          write(lout,"(a)") "RESO> ERROR Elements specified for resonance compensation are not sextupoles."
          iErr = .true.
          return
        end if
        if(nqc == 1 .and. (k == 9 .or. k == 10) .and. kz(j) /= 2) then
          write(lout,"(a)") "RESO> ERROR Elements specified for resonance compensation are not quadrupoles."
          iErr = .true.
          return
        end if
        cycle outer
      end do inner
      if((nre == 1 .and.  k < 3) .or. &
         (nre == 2 .and.  k < 5) .or. &
         (nre == 3 .and.  k < 7) .or. &
         (nch == 1 .and. (k == 7 .or. k == 8)) .or. &
         (nqc == 1 .and. (k == 9 .or. k == 10))) then
        write(lout,"(a)") "RESO> ERROR Element is not in the element list."
        iErr = .true.
        return
      end if
    end do outer

    irmod2 = 1

  case default
    write(lout,"(a,i0)") "RESO> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end select

end subroutine sixin_parseInputLineRESO

! ================================================================================================ !
!  Parse Search for Optimum Places to Compensate Resonances Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-23
! ================================================================================================ !
subroutine sixin_parseInputLineSEAR(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) name(40)
  integer nSplit, j, k, k0, ka, ke, ki
  logical spErr

  save :: name, k0

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "SEAR> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    if(nSplit > 0) call chr_cast(lnSplit(1),qxt, iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),qzt, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),tam1,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),tam2,iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),totl,iErr)

    if(st_debug) then
      call sixin_echoVal("qxt", qxt, "SEAR",iLine)
      call sixin_echoVal("qzt", qzt, "SEAR",iLine)
      call sixin_echoVal("tam1",tam1,"SEAR",iLine)
      call sixin_echoVal("tam2",tam2,"SEAR",iLine)
      call sixin_echoVal("totl",totl,"SEAR",iLine)
    end if
    if(iErr) return

  case(2)

    if(nSplit > 0) call chr_cast(lnSplit(1),mesa,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),mp,  iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),m21, iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),m22, iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),m23, iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),ise1,iErr)
    if(nSplit > 6) call chr_cast(lnSplit(7),ise2,iErr)
    if(nSplit > 7) call chr_cast(lnSplit(8),ise3,iErr)

    if(st_debug) then
      call sixin_echoVal("mesa",mesa,"SEAR",iLine)
      call sixin_echoVal("mp",  mp,  "SEAR",iLine)
      call sixin_echoVal("m21", m21, "SEAR",iLine)
      call sixin_echoVal("m22", m22, "SEAR",iLine)
      call sixin_echoVal("m23", m23, "SEAR",iLine)
      call sixin_echoVal("ise1",ise1,"SEAR",iLine)
      call sixin_echoVal("ise2",ise2,"SEAR",iLine)
      call sixin_echoVal("ise3",ise3,"SEAR",iLine)
    end if
    if(iErr) return

    k0 = 0

  case default

    name(:) = " "

    ka = k0 + 1
    ke = k0 + nSplit
    do k=ka,ke
      if(k > nele) then
        write(lout,"(a)") "SEAR> ERROR Cannot have more elements than in the single elements block."
        iErr = .true.
        return
      end if
      if(k > mesa) return
      ki = k-k0
      do j=1,il
        if(bez(j) == lnSplit(ki)) then
          isea(k) = j
          if(abs(kz(j)) /= mp) then
            write(lout,"(a)") "SEAR> ERROR With the specified elements the resonance cannot be compensated."
            write(lout,"(a)") "SEAR> ERROR Resonance order and element type # must be the same."
            iErr = .true.
            return
          else
            ise = 1
          end if
          exit
        end if
      end do
      if(isea(k) == 0) then
        write(lout,"(a)") "SEAR> ERROR Element is not in the element list."
        iErr = .true.
        return
      end if
    end do
    k0 = k-1

  end select

end subroutine sixin_parseInputLineSEAR

! ================================================================================================ !
!  Parse Post-Processing Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-23
! ================================================================================================ !
subroutine sixin_parseInputLinePOST(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "POST> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)
    toptit(1) = trim(inLine)

  case(2)

    if(nSplit > 0)  call chr_cast(lnSplit(1) ,iav,   iErr)
    if(nSplit > 1)  call chr_cast(lnSplit(2) ,nstart,iErr)
    if(nSplit > 2)  call chr_cast(lnSplit(3) ,nstop, iErr)
    if(nSplit > 3)  call chr_cast(lnSplit(4) ,iwg,   iErr)
    if(nSplit > 4)  call chr_cast(lnSplit(5) ,dphix, iErr)
    if(nSplit > 5)  call chr_cast(lnSplit(6) ,dphiz, iErr)
    if(nSplit > 6)  call chr_cast(lnSplit(7) ,iskip, iErr)
    if(nSplit > 7)  call chr_cast(lnSplit(8) ,iconv, iErr)
    if(nSplit > 8)  call chr_cast(lnSplit(9) ,imad,  iErr)
    if(nSplit > 9)  call chr_cast(lnSplit(10),cma1,  iErr)
    if(nSplit > 10) call chr_cast(lnSplit(11),cma2,  iErr)

    if(st_debug) then
      call sixin_echoVal("iav",   iav,   "POST",iLine)
      call sixin_echoVal("nstart",nstart,"POST",iLine)
      call sixin_echoVal("nstop", nstop, "POST",iLine)
      call sixin_echoVal("iwg",   iwg,   "POST",iLine)
      call sixin_echoVal("dphix", dphix, "POST",iLine)
      call sixin_echoVal("dphiz", dphiz, "POST",iLine)
      call sixin_echoVal("iskip", iskip, "POST",iLine)
      call sixin_echoVal("iconv", iconv, "POST",iLine)
      call sixin_echoVal("imad",  imad,  "POST",iLine)
      call sixin_echoVal("cma1",  cma1,  "POST",iLine)
      call sixin_echoVal("cma2",  cma2,  "POST",iLine)
    end if
    if(iErr) return

#ifdef STF
    if(imad == 1) then
      write(lout,"(a)") "POST> ERROR imad not supported when SixTrack is built with STF enabled."
      iErr = .true.
      return
    end if
#endif

  case(3)

    if(nSplit > 0) call chr_cast(lnSplit(1),qx0, iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),qz0, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),ivox,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),ivoz,iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),ires,iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),dres,iErr)
    if(nSplit > 6) call chr_cast(lnSplit(7),ifh, iErr)
    if(nSplit > 7) call chr_cast(lnSplit(8),dfft,iErr)

    if(st_debug) then
      call sixin_echoVal("qx0", qx0, "POST",iLine)
      call sixin_echoVal("qz0", qz0, "POST",iLine)
      call sixin_echoVal("ivox",ivox,"POST",iLine)
      call sixin_echoVal("ivoz",ivoz,"POST",iLine)
      call sixin_echoVal("ires",ires,"POST",iLine)
      call sixin_echoVal("dres",dres,"POST",iLine)
      call sixin_echoVal("ifh", ifh, "POST",iLine)
      call sixin_echoVal("dfft",dfft,"POST",iLine)
    end if
    if(iErr) return

  case(4)

    if(nSplit > 0) call chr_cast(lnSplit(1),kwtype,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),itf,   iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),icr,   iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),idis,  iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),icow,  iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),istw,  iErr)
    if(nSplit > 6) call chr_cast(lnSplit(7),iffw,  iErr)
    if(nSplit > 7) call chr_cast(lnSplit(8),nprint,iErr)
    if(nSplit > 8) call chr_cast(lnSplit(9),ndafi, iErr)

    if(st_debug) then
      call sixin_echoVal("kwtype",kwtype,"POST",iLine)
      call sixin_echoVal("itf",   itf,   "POST",iLine)
      call sixin_echoVal("icr",   icr,   "POST",iLine)
      call sixin_echoVal("idis",  idis,  "POST",iLine)
      call sixin_echoVal("icow",  icow,  "POST",iLine)
      call sixin_echoVal("istw",  istw,  "POST",iLine)
      call sixin_echoVal("iffw",  iffw,  "POST",iLine)
      call sixin_echoVal("nprint",nprint,"POST",iLine)
      call sixin_echoVal("ndafi", ndafi, "POST",iLine)
    end if
    if(iErr) return

    kwtype = 0
    icr    = 0
    if(iskip  <  0) iskip  = 1
    if(nprint /= 1) nprint = 0
    if(nstart <  0) nstart = 0
    if(nstop  <  0) nstop  = 0
    if(nstop < nstart) then
      nstart = 0
      nstop  = 0
    end if
    if(iconv /= 1) iconv = 0
    if(abs(cma1) <= pieni) cma1 = one
    cma1 = cma1*c1e3
    if(abs(cma2) <= pieni) cma2 = one
    ipos = 1 ! Turn postprocessing ON.

  case default
    write(lout,"(a,i0)") "POST> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end select

end subroutine sixin_parseInputLinePOST

! ================================================================================================ !
!  Parse Decoupling of Motion in the Transverse Planes Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-23
! ================================================================================================ !
subroutine sixin_parseInputLineDECO(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) name(6)
  integer nSplit, i, j, k
  logical spErr

  save :: name

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "DECO> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    name(:) = " "

    if(nSplit > 0) name(1) = lnSplit(1)
    if(nSplit > 1) name(2) = lnSplit(2)
    if(nSplit > 2) name(3) = lnSplit(3)
    if(nSplit > 3) name(4) = lnSplit(4)

    if(st_debug) then
      call sixin_echoVal("name1",name(1),"DECO",iLine)
      call sixin_echoVal("name2",name(2),"DECO",iLine)
      call sixin_echoVal("name3",name(3),"DECO",iLine)
      call sixin_echoVal("name4",name(4),"DECO",iLine)
    end if

    iskew = 1

  case(2)

    if(nSplit > 0) name(5) = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),qwsk(1),iErr)

    if(st_debug) then
      call sixin_echoVal("name5",  name(5),"DECO",iLine)
      call sixin_echoVal("qwsk(1)",qwsk(1),"DECO",iLine)
    end if
    if(iErr) return

  case(3)

    if(nSplit > 0) name(6) = lnSplit(1)
    if(nSplit > 1) call chr_cast(lnSplit(2),qwsk(2),iErr)

    if(st_debug) then
      call sixin_echoVal("name6",  name(6),"DECO",iLine)
      call sixin_echoVal("qwsk(2)",qwsk(2),"DECO",iLine)
    end if
    if(iErr) return

    iskew = 2

    do i=1,6
      do j=1,il
        if(iskew == 2 .and. i > 4) return
        if(bez(j) == name(i)) then
          if(i <= 4) then
            if(kz(j) /= -2) then
              write(lout,"(a)") "DECO> ERROR Elements specified is not a skew quadrupole."
              iErr = .true.
              return
            end if
          else
            if(kz(j) /= 2) then
              write(lout,"(a)") "DECO> ERROR Elements specified is not a quadrupole."
              iErr = .true.
              return
            end if
          end if
          nskew(i) = j
          do k=1,6
            if(nskew(k) /= 0 .and. (nskew(k) == nskew(i)) .and. (k /= i)) then
              write(lout,"(a)") "DECO> ERROR Same element specified twice."
              iErr = .true.
              return
            end if
          end do
        end if
      end do
    end do

  case default
    write(lout,"(a,i0)") "DECO> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end select

end subroutine sixin_parseInputLineDECO

! ================================================================================================ !
!  Parse Normal Forms Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-23
! ================================================================================================ !
subroutine sixin_parseInputLineNORM(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use mod_common
  use mod_common_da

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "NORM> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    ! FIXME: This should be moved to post ENDE checks
    if(idial == 0 .and. numl == 0) then
      write(lout,"(a)") "NORM> ERROR Normal forms analysis impossible. The transfer map does not exist!"
      iErr = .true.
      return
    end if
    inorm = 1

    if(nSplit > 0) call chr_cast(lnSplit(1),nordf,iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),nvarf,iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),nord1,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),idptr,iErr)

    if(st_debug) then
      call sixin_echoVal("nordf",nordf,"NORM",iLine)
      call sixin_echoVal("nvarf",nvarf,"NORM",iLine)
      call sixin_echoVal("nord1",nord1,"NORM",iLine)
      call sixin_echoVal("idptr",idptr,"NORM",iLine)
    end if

    if(nord /= 0 .and. nordf > nord+1) then
      imod1 = 1
    end if
    if(nvar /= 0 .and. nvarf> nvar) then
      nvarf = nvar
      imod2 = 1
    end if
    if(idptr < 0 .or. idptr > 6) then
      idptr = 0
    end if

  case default
    write(lout,"(a,i0)") "NORM> ERROR Unexpected line number ",iLine
    iErr = .true.
    return

  end select

end subroutine sixin_parseInputLineNORM

! ================================================================================================ !
!  Parse Beam–Beam Element Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-23
! ================================================================================================ !
subroutine sixin_parseInputLineBEAM(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use parbeam, only : beam_expflag
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) elemName
  real(kind=fPrec) xang,xstr,xplane
  integer nSplit, ibsix, j
  logical spErr, beamXStr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "BEAM> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 1) then

    if(lnSplit(1) == "EXPERT") then
      beam_expflag = 1
      write(lout,"(a)") "BEAM> EXPERT mode enabled."
      return
    end if

    write(lout,"(a)") "BEAM> Reading old style beam block."
    write(lout,"(a)") "BEAM>    To convert to the new format, copy-paste these lines into the BEAM block in fort.3,"
    write(lout,"(a)") "BEAM> replacing line 2 onwards. Then write EXPERT on the first line of the BEAM block, above"
    write(lout,"(a)") "BEAM> the current first line. Finally, in the SINGLE ELEMENTS list (normally in fort.2) set "
    write(lout,"(a)") "BEAM> the parameters of all beam-beam lenses (type 20) to 0.0."
    write(lout,"(a)") "BEAM>    This procedure produces a new set of input files that should have bit-for-bit iden-"
    write(lout,"(a)") "BEAM> tical results to this one."
    write(lout,"(a)") "BEAM>    The easiest way to check this is to run both simulations side-by-side and compare"
    write(lout,"(a)") "BEAM> the standard output in a text diff tool like meld. If the results are not identical,"
    write(lout,"(a)") "BEAM> this is a bug; please report it to the developers."
#ifndef CRLIBM
    write(lout,"(a)") "BEAM> WARNING This sixtrack binary was not compiled with crlibm, conversion will not be exact."
#endif

    if(nSplit > 0) call chr_cast(lnSplit(1),partnum,     iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),sixin_emitNX,iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),sixin_emitNY,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),sigz,        iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),sige,        iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),ibeco,       iErr)
    if(nSplit > 6) call chr_cast(lnSplit(7),ibtyp,       iErr)
    if(nSplit > 7) call chr_cast(lnSplit(8),lhc,         iErr)
    if(nSplit > 8) call chr_cast(lnSplit(9),ibbc,        iErr)

    if(st_debug) then
      call sixin_echoVal("partnum",partnum,     "BEAM",iLine)
      call sixin_echoVal("emitnx", sixin_emitNX,"BEAM",iLine)
      call sixin_echoVal("emitny", sixin_emitNY,"BEAM",iLine)
      call sixin_echoVal("sigz",   sigz,        "BEAM",iLine)
      call sixin_echoVal("sige",   sige,        "BEAM",iLine)
      call sixin_echoVal("ibeco",  ibeco,       "BEAM",iLine)
      call sixin_echoVal("ibtyp",  ibtyp,       "BEAM",iLine)
      call sixin_echoVal("lhc",    lhc,         "BEAM",iLine)
      call sixin_echoVal("ibbc",   ibbc,        "BEAM",iLine)
    end if
    if(iErr) return

    if(nSplit /= 9) then
      write(lout,"(a,i0)") "BEAM> WARNING (not EXPERT). First line should have 9 fields, got ",nSplit
    end if

    if(sixin_emitNX <= pieni .or. sixin_emitNY <= pieni) then
      write(lout,"(a)") "BEAM> ERROR Either normalised emittances or the resulting sigma values equal to zero."
      iErr = .true.
      return
    end if

    if(ibeco /= 0 .and. ibeco /= 1) ibeco = 1
    if(ibtyp /= 0 .and. ibtyp /= 1) ibtyp = 0
    if(ibbc  /= 0 .and. ibbc  /= 1) ibbc  = 0
    if(lhc    < 0 .or.  lhc    > 2) lhc   = 1

    nbeam = 1

    if(ibtyp == 1) call wzset

  else

    beamXStr = .false.
    if(nSplit == 5) then
      beamXStr = .true.
    elseif(nSplit == 4) then
      beamXStr = .false.
    else
      write(lout,"(a,i0)") "BEAM> ERROR Number of arguments in line 2 is expected to be 4 or 5, got ",nSplit
      iErr = .true.
      return
    end if

    if(nSplit > 0) elemName = trim(lnSplit(1))
    if(nSplit > 1) call chr_cast(lnSplit(2),ibsix, iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),xang,  iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),xplane,iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),xstr,  iErr)

    if(st_debug) then
      call sixin_echoVal("name",  elemName,"BEAM",iLine)
      call sixin_echoVal("ibsix", ibsix,   "BEAM",iLine)
      call sixin_echoVal("xang",  xang,    "BEAM",iLine)
      call sixin_echoVal("xplane",xplane,  "BEAM",iLine)
      call sixin_echoVal("xstr",  xstr,    "BEAM",iLine)
    end if
    if(iErr) return

    if(.not. beamXStr) then
      write(lout,"(a)") "BEAM> WARNING No xstr present, assuming xstr = xang."
      xstr = xang
    end if

    if(ibsix < 0) ibsix = 0
    do j=1,il
      if(bez(j) == elemName .and. kz(j) == 20) then
        ibb6d       = 1
        parbe(j,2)  = real(ibsix,fPrec)
        parbe(j,1)  = xang
        parbe(j,3)  = xplane
        parbe(j,18) = xstr
        exit
      end if
    end do

  end if

end subroutine sixin_parseInputLineBEAM

! ================================================================================================ !
!  Parse Beam–Beam Element Line EXPERT
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-23
! ================================================================================================ !
subroutine sixin_parseInputLineBEAM_EXP(inLine, iLine, iErr)

  use crcoall
  use string_tools
  use mod_settings
  use parbeam, only : beam_expflag
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) elemName
  real(kind=fPrec) sxx,syy,sxy,separx,separy,mm(11)
  integer nSplit, n6D, ibsix, j
  logical spErr, beamXStr

  save :: n6D,elemName,ibsix,sxx,syy,sxy,separx,separy,mm

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "BEAM> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine == 2) then

    if(nSplit > 0) call chr_cast(lnSplit(1),partnum,     iErr)
    if(nSplit > 1) call chr_cast(lnSplit(2),sixin_emitNX,iErr)
    if(nSplit > 2) call chr_cast(lnSplit(3),sixin_emitNY,iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),sigz,        iErr)
    if(nSplit > 4) call chr_cast(lnSplit(5),sige,        iErr)
    if(nSplit > 5) call chr_cast(lnSplit(6),ibeco,       iErr)
    if(nSplit > 6) call chr_cast(lnSplit(7),ibtyp,       iErr)
    if(nSplit > 7) call chr_cast(lnSplit(8),lhc,         iErr)
    if(nSplit > 8) call chr_cast(lnSplit(9),ibbc,        iErr)

    if(st_debug) then
      call sixin_echoVal("partnum",partnum,     "BEAM",iLine)
      call sixin_echoVal("emitnx", sixin_emitNX,"BEAM",iLine)
      call sixin_echoVal("emitny", sixin_emitNY,"BEAM",iLine)
      call sixin_echoVal("sigz",   sigz,        "BEAM",iLine)
      call sixin_echoVal("sige",   sige,        "BEAM",iLine)
      call sixin_echoVal("ibeco",  ibeco,       "BEAM",iLine)
      call sixin_echoVal("ibtyp",  ibtyp,       "BEAM",iLine)
      call sixin_echoVal("lhc",    lhc,         "BEAM",iLine)
      call sixin_echoVal("ibbc",   ibbc,        "BEAM",iLine)
    end if
    if(iErr) return

    if(nSplit /= 9) then
      write(lout,"(a,i0)") "BEAM> WARNING (not EXPERT). First line should have 9 fields, got ",nSplit
    end if

    if(sixin_emitNX <= pieni .or. sixin_emitNY <= pieni) then
      write(lout,"(a)") "BEAM> ERROR Either normalised emittances or the resulting sigma values equal to zero."
      iErr = .true.
      return
    end if

    if(ibeco /= 0 .and. ibeco /= 1) ibeco = 1
    if(ibtyp /= 0 .and. ibtyp /= 1) ibtyp = 0
    if(ibbc  /= 0 .and. ibbc  /= 1) ibbc  = 0
    if(lhc    < 0 .or.  lhc    > 2) lhc   = 1

    nbeam = 1

    if(ibtyp == 1) call wzset

    n6D = 0

  else

    if(n6D == 0) then

      mm(:)  = zero
      sxx    = zero
      syy    = zero
      separx = zero
      separy = zero
      sxy    = zero

      if(nSplit > 0) elemName = trim(lnSplit(1))
      if(nSplit > 1) call chr_cast(lnSplit(2),ibsix, iErr)
      if(nSplit > 2) call chr_cast(lnSplit(3),sxx,   iErr)
      if(nSplit > 3) call chr_cast(lnSplit(4),syy,   iErr)
      if(nSplit > 4) call chr_cast(lnSplit(5),separx,iErr)
      if(nSplit > 5) call chr_cast(lnSplit(6),separy,iErr)
      if(nSplit > 6) call chr_cast(lnSplit(7),mm(1), iErr)
      if(nSplit > 7) call chr_cast(lnSplit(8),sxy,   iErr)

      if(ibsix == 0) then
        if( nSplit < 7 .or. nSplit > 8 ) then
          write(lout,"(a,i0)") "BEAM> ERROR First line of a 4D element definition should have 7  or 8 fields, got ",nSplit
          iErr = .true.
          return
        end if
        if(st_debug) then
          write(lout,"(a)") "BEAM> DEBUG New 4D element encountered."
        end if
        n6D = 0
      else if(ibsix > 0) then
        n6D = 1
        if(nSplit /= 6) then
          write(lout,"(a,i0)") "BEAM> ERROR First line of a 6D element definition should have 6 fields, got ",nSplit
          iErr = .true.
          return
        end if
        if(st_debug) then
          write(lout,"(a)") "BEAM> DEBUG New 6D element encountered."
        end if
      else
        write(lout,"(a,i0)") "BEAM> ERROR Expected number of slices >= 0, but got ",ibsix
        iErr = .true.
        return
      end if

      if(st_debug) then
        call sixin_echoVal("name",  elemName,"BEAM",iLine)
        call sixin_echoVal("ibsix", ibsix,   "BEAM",iLine)
        call sixin_echoVal("Sxx",   sxx,     "BEAM",iLine)
        call sixin_echoVal("Syy",   syy,     "BEAM",iLine)
        call sixin_echoVal("Sxy",   sxy,     "BEAM",iLine)
        call sixin_echoVal("separx",separx,  "BEAM",iLine)
        call sixin_echoVal("separy",separy,  "BEAM",iLine)
        call sixin_echoVal("strrat",mm(1),   "BEAM",iLine)
      end if
      if(iErr) return

      if(n6D == 0) then
        ! Save 4D data
        do j=1,il
          if(bez(j) == elemName) then
            if(kz(j) /= 20) then
              write(lout,"(a,i0,a)") "BEAM> ERROR Found element '"//bez(j)//"' type ",kz(j), ", but expected type 20."
              iErr = .true.
              return
            else
              if(parbe(j,5) /= zero .or. parbe(j,6) /= zero .or. ptnfac(j)  /= zero .or. &
                 bbbx(j)    /= zero .or. bbby(j)    /= zero .or. bbbs(j)    /= zero) then
                write(lout,"(a)") "BEAM> ERROR Using EXPERT mode, but element '"//bez(j)//&
                  " does not have ed=ek=el=bbbx=bbby=bbbs=0.0 in the SINGLE ELEMENTS list."
                iErr = .true.
                return
              end if
              parbe(j,17) = 0
              parbe(j,2)  = real(ibsix,fPrec)
              parbe(j,1)  = sxx
              parbe(j,3)  = syy
              parbe(j,5)  = separx
              parbe(j,6)  = separy
              ptnfac(j)   = mm(1)
              parbe(j,13) = sxy
            end if
          end if
        end do
      end if

    else if(n6D == 1) then

      if(st_debug) then
        write(lout,"(a)") "BEAM> DEBUG Second 6D line encountered."
      end if
      if(nSplit /= 5) then
        write(lout,"(a,i0)") "BEAM> ERROR Second line of a 6D element definition should have 5 fields, got ",nSplit
        iErr = .true.
        return
      end if

      if(nSplit > 0) call chr_cast(lnSplit(1),mm(1),iErr)
      if(nSplit > 1) call chr_cast(lnSplit(2),mm(2),iErr)
      if(nSplit > 2) call chr_cast(lnSplit(3),mm(3),iErr)
      if(nSplit > 3) call chr_cast(lnSplit(4),mm(4),iErr)
      if(nSplit > 4) call chr_cast(lnSplit(5),mm(5),iErr)

      if(st_debug) then
        call sixin_echoVal("Sxx",  mm(1),"BEAM",iLine)
        call sixin_echoVal("Sxxp", mm(2),"BEAM",iLine)
        call sixin_echoVal("Sxpxp",mm(3),"BEAM",iLine)
        call sixin_echoVal("Syy",  mm(4),"BEAM",iLine)
        call sixin_echoVal("Syyp", mm(5),"BEAM",iLine)
      end if
      if(iErr) return

      n6D = 2

    else if(n6D == 2) then

      if(st_debug) then
        write(lout,"(a)") "BEAM> DEBUG Third 6D line encountered."
      end if
      if(nSplit /= 6) then
        write(lout,"(a,i0)") "BEAM> ERROR Tird line of a 6D element definition should have 6 fields, got ",nSplit
        iErr = .true.
        return
      end if

      if(nSplit > 0) call chr_cast(lnSplit(1),mm(6), iErr)
      if(nSplit > 1) call chr_cast(lnSplit(2),mm(7), iErr)
      if(nSplit > 2) call chr_cast(lnSplit(3),mm(8), iErr)
      if(nSplit > 3) call chr_cast(lnSplit(4),mm(9), iErr)
      if(nSplit > 4) call chr_cast(lnSplit(5),mm(10),iErr)
      if(nSplit > 5) call chr_cast(lnSplit(6),mm(11),iErr)

      if(st_debug) then
        call sixin_echoVal("Sypyp", mm(6), "BEAM",iLine)
        call sixin_echoVal("Sxy",   mm(7), "BEAM",iLine)
        call sixin_echoVal("Sxyp",  mm(8), "BEAM",iLine)
        call sixin_echoVal("Sxpy",  mm(9), "BEAM",iLine)
        call sixin_echoVal("Sxpyp", mm(10),"BEAM",iLine)
        call sixin_echoVal("strrat",mm(11),"BEAM",iLine)
      end if
      if(iErr) return

      ! Save 6D data
      do j=1,il
        if(bez(j) == elemName) then
          if(kz(j) /= 20) then
            write(lout,"(a,i0,a)") "BEAM> ERROR Found element '"//bez(j)//"' type ",kz(j), ", but expected type 20."
            iErr = .true.
            return
          else
            if(parbe(j,5) /= zero .or. parbe(j,6) /= zero .or. ptnfac(j)  /= zero .or. &
               bbbx(j)    /= zero .or. bbby(j)    /= zero .or. bbbs(j)    /= zero) then
              write(lout,"(a)") "BEAM> ERROR Using EXPERT mode, but element '"//bez(j)//&
                "' does not have ed=ek=el=bbbx=bbby=bbbs=0.0 in the SINGLE ELEMENTS list."
              iErr = .true.
              return
            end if
            parbe(j,17) = 1
            parbe(j,2)  = real(ibsix,fPrec)
            parbe(j,1)  = sxx
            parbe(j,3)  = syy
            parbe(j,5)  = separx
            parbe(j,6)  = separy
            parbe(j,7)  = mm(1)
            parbe(j,8)  = mm(2)
            parbe(j,9)  = mm(3)
            parbe(j,10) = mm(4)
            parbe(j,11) = mm(5)
            parbe(j,12) = mm(6)
            parbe(j,13) = mm(7)
            parbe(j,14) = mm(8)
            parbe(j,15) = mm(9)
            parbe(j,16) = mm(10)
            ptnfac(j)   = mm(11)
          end if
        end if
     end do

     n6D = 0

    end if

  end if

end subroutine sixin_parseInputLineBEAM_EXP

! ================================================================================================ !
!  Parse “Phase Trombone” Element Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-06-25
! ================================================================================================ !
subroutine sixin_parseInputLineTROM(inLine, iLine, iErr)

  use crcoall
  use mod_alloc
  use string_tools
  use mod_settings
  use mod_common

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) elemName
  real(kind=fPrec) cloOrb(6), matElems(6,6)
  integer nSplit, nLines, iElem, i, l, m, n
  logical spErr

  save :: elemName, iElem, cloOrb, matElems, l, m, n, nLines

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "TROM> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(iLine)

  case(1)

    if(nSplit /= 1) then
      write(lout,"(a,i0)") "TROM> ERROR Expected 1 element name for line 1, got ",nSplit
      iErr = .true.
      return
    end if

    elemName = trim(lnSplit(1))

    iElem = -1
    do i=1,il
      if(bez(i) == elemName) then
        iElem = i
        exit
      end if
    end do
    if(iElem == -1) then
      write(lout,"(a)") "TROM> ERROR Could not find element '"//elemName//"' in single element list."
      iErr = .true.
      return
    end if

    if(st_debug) then
      call sixin_echoVal("name",  elemName,"TROM",iLine)
      call sixin_echoVal("elemid",iElem,   "TROM",iLine)
    end if
    if(iErr) return

    cloOrb(:)     = zero
    matElems(:,:) = zero
    nLines        = 1

    l = 0
    m = 0
    n = 1

  case(2,3)

    if(nSplit /= 3) then
      write(lout,"(2(a,i0))") "TROM> ERROR Expected 3 values for line ",iLine,", got ",nSplit
      iErr = .true.
      return
    end if

    l = l + 3

    call chr_cast(lnSplit(1),cloOrb(l-2),iErr)
    call chr_cast(lnSplit(2),cloOrb(l-1),iErr)
    call chr_cast(lnSplit(3),cloOrb(l),  iErr)
    if(iErr) return

    nLines = nLines + 1

  case(4,5,6,7,8,9,10,11,12,13,14,15)

    if(nSplit /= 3) then
      write(lout,"(2(a,i0))") "TROM> ERROR Expected 3 values for line ",iLine,", got ",nSplit
      iErr = .true.
      return
    end if

    m = m + 3
    if(m > 6) then
      n = n + 1
      m = 3
    end if

    call chr_cast(lnSplit(1),matElems(n,m-2),iErr)
    call chr_cast(lnSplit(2),matElems(n,m-1),iErr)
    call chr_cast(lnSplit(3),matElems(n,m),  iErr)
    if(iErr) return

    nLines = nLines + 1

  case(-1)

    if(nLines /= 15) then
      write(lout,"(a,i0)") "TROM> ERROR Each trombone block takes exactly 15 lines, got ",nLines
      write(lout,"(a)")    "TROM>       If you neeed multiple TROM elements, add multiple TROM blocks."
      iErr = .true.
      return
    end if

    sixin_imtr0 = sixin_imtr0 + 1
    imtr(iElem) = sixin_imtr0
    ntr         = sixin_imtr0
    call alloc(cotr,sixin_imtr0,6,  zero,"cotr")
    call alloc(rrtr,sixin_imtr0,6,6,zero,"rrtr")
    cotr(sixin_imtr0,1:6)     = cloOrb(1:6)
    rrtr(sixin_imtr0,1:6,1:6) = matElems(1:6,1:6)

    if(st_debug) then
      call sixin_echoVal("cx", cloOrb(1),"TROM",iLine)
      call sixin_echoVal("cx'",cloOrb(2),"TROM",iLine)
      call sixin_echoVal("cy", cloOrb(3),"TROM",iLine)
      call sixin_echoVal("cy'",cloOrb(4),"TROM",iLine)
      call sixin_echoVal("cz", cloOrb(5),"TROM",iLine)
      call sixin_echoVal("cz'",cloOrb(6),"TROM",iLine)
      write(lout,"(a)") "INPUT> DEBUG TROM:PP Matrix Elements:"
      do i=1,6
        write(lout,"(a,6(1x,e15.8),a)") "INPUT> DEBUG TROM:PP [ ",rrtr(sixin_imtr0,i,1:6)," ]"
      end do
    end if

  end select

end subroutine sixin_parseInputLineTROM

end module sixtrack_input
