! ================================================================================================ !
!  Random Fluctuation Starting Number Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  If besides mean values for the multipole errors (Gaussian) random errors should be considered,
!  this module is used to set the start value for the random generator.
!
!  Moved from main code and updated by V.K. Berglyd Olsen, June 2018
!  Last modified: 2018-06-13
! ================================================================================================ !
module mod_fluc

  use crcoall
  use floatPrecision
  use parpro,         only : nmac
  use mod_common,     only : izu0,mmac,mcut
  use mod_settings,   only : st_debug
  use sixtrack_input, only : sixin_echoVal

  implicit none

  real(kind=fPrec), allocatable, public, save :: fluc_errExt(:,:) ! The errors from fort.16
  character(len=:), allocatable, public, save :: fluc_bezExt(:)   ! The name of the element
  integer,          allocatable, public, save :: fluc_ixExt(:)    ! The index of the element
  integer,                       public, save :: fluc_nExt        ! Number of multipoles

contains

subroutine fluc_parseInputLine(inLine, iLine, iErr, mout)

  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr
  integer,          intent(out)   :: mout

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "FLUC> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(iLine > 1) then
    write(lout,"(a)") "FLUC> ERROR This block only takes one line."
    iErr = .true.
    return
  end if

  if(nSplit > 0) call chr_cast(lnSplit(1),izu0,iErr)
  if(nSplit > 1) call chr_cast(lnSplit(2),mmac,iErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),mout,iErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),mcut,iErr)

  if(st_debug) then
    call sixin_echoVal("izu0",izu0,"FLUC",iLine)
    call sixin_echoVal("mmac",mmac,"FLUC",iLine)
    call sixin_echoVal("mout",mout,"FLUC",iLine)
    call sixin_echoVal("mcut",mcut,"FLUC",iLine)
  end if

  ! Process variables
  mcut = iabs(mcut)

  if(mmac > nmac) then
    write(lout,"(a,i0)") "FLUC> ERROR Maximum number of seeds for vectorisation is ",nmac
    iErr = .true.
    return
  end if

end subroutine fluc_parseInputLine

subroutine fluc_readFort8
end subroutine fluc_readFort8

subroutine fluc_readFort16

  use mod_alloc
  use string_tools
  use numerical_constants, only : zero
  use mod_common,          only : il,bez,icext,nblz,nblo,ic
  use mod_settings,        only : st_debug

  implicit none

  character(len=1024) :: inLine
  character(len=:), allocatable :: lnSplit(:)
  character(len=str_maxName) bezExt
  integer mType, lineNo, ioStat, nSplit, lMode, nVals, mVal, mExt, iStru, iSing, nExt
  logical iErr, isOpen, inSing
  integer i, j

  mType  = 0
  lineNo = 0
  lMode  = 0

  fluc_nExt = 0

  inquire(unit=16, opened=isOpen)
  if(isOpen) close(16)
  open(16,file="fort.16")
  rewind(16)

  write(lout,"(a)") "FLUC> Reading multipole errors from fort.16"

10 continue
  read(16,"(a)",end=20,iostat=ioStat) inLine
  lineNo = lineNo + 1
  if(ioStat /= 0) then
    write(lout,"(a,i0)") "FLUC> ERROR fort.16 iostat = ",ioStat
    goto 30
  end if
  
  call chr_split(inLine, lnSplit, nSplit, iErr)
  if(iErr) goto 30
  if(nSplit == 0) goto 10
  if(lMode == 0) then
    ! We're expecting an element name
    if(nSplit /= 1) then
      write(lout,"(a)") "FLUC> ERROR Expected a single element name."
      goto 30
    end if
    bezExt = trim(lnSplit(1))

    ! Check that the element exists
    inSing = .false.
    do i=1,il
      if(bez(i) == bezExt) then
        inSing = .true.
        exit
      end if
    end do
    if(inSing) then
      
      lMode = 1
      nVals = 0
      mVal  = 0
      
      if(allocated(fluc_ixExt)) then
        mExt = size(fluc_ixExt,1)
      else
        mExt = 0
      end if
      
      fluc_nExt = fluc_nExt + 1
      if(fluc_nExt > mExt) then
        call alloc(fluc_errExt,40,         mExt+10,zero,       "fluc_errExt")
        call alloc(fluc_bezExt,str_maxName,mExt+10,str_nmSpace,"fluc_bezExt")
        call alloc(fluc_ixExt,             mExt+10,0,          "fluc_ixExt")
      end if
      fluc_bezExt(fluc_nExt) = bezExt
      
    else
      write(lout,"(a)") "FLUC> ERROR Unknown element name '"//trim(bezExt)//"'."
      goto 30
    end if
  else
    ! We're expecting 40 float values
    nVals = nVals + nSplit
    if(nVals > 40) then
      write(lout,"(a,i0)") "FLUC> ERROR Each element takes exactly 40 values, got ",nVals
      goto 30
    end if
    do i=1,nSplit
      mVal = mVal + 1
      call chr_cast(lnSplit(i),fluc_errExt(mVal,fluc_nExt),iErr)
      if(iErr) goto 30
    end do
    if(mVal == 40) then
      lMode = 0
    end if
  end if

  goto 10

20 continue
  close(16)
  
  ! We require that the elements in fort.16 have the same order as in the STRUcture block.
  ! A repeated element name implies it is the next one of that name in the sequenc, so we will
  ! iterate over the index array ic() and check for a matchin order to build the index icext.
  nExt = 1
  do iStru=1,nblz
    iSing = ic(iStru)
    if(iSing <= nblo) cycle
    iSing = iSing - nblo
    if(bez(iSing) == fluc_bezExt(nExt)) then
      fluc_ixExt(nExt) = iStru
      icext(iStru)     = nExt
      nExt = nExt + 1
    end if
  end do

  write(lout,"(a,i0,a)") "FLUC> Read ",fluc_nExt," values from fort.16"
  if(st_debug) then
    do i=1,fluc_nExt
      write(lout,"(a,i5,a,i0)") "FLUC> ",i,": '"//trim(fluc_bezExt(i))//"' with ID ",fluc_ixExt(i)
      do j=1,40,4
        write(lout,"(a,i2,a,4(2x,e14.7))")  "FLUC>  * ",j,": ",&
          fluc_errExt(j,i),fluc_errExt(j+1,i),fluc_errExt(j+2,i),fluc_errExt(j+3,i)
      end do
    end do
  end if

  return

30 continue
  write(lout,"(a,i0,a)") "FLUC> ERROR Line ",lineNo," in fort.16"
  write(lout,"(a)")      "FLUC>       '"//trim(inLine)//"'"
  call prror(-1)

end subroutine fluc_readFort16

end module mod_fluc
