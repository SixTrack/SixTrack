! ================================================================================================ !
!  Random Fluctuation Starting Number Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  If besides mean values for the multipole errors (Gaussian) random errors should be considered,
!  this module is used to set the start value for the random generator.
!
!  Moved from main code and updated by V.K. Berglyd Olsen, June 2018
!  Last modified: 2018-06-14
! ================================================================================================ !
module mod_fluc

  use crcoall
  use floatPrecision

  implicit none

  real(kind=fPrec), allocatable, public, save :: fluc_errAlign(:,:) ! The alignemnt errors from fort.8
  character(len=:), allocatable, public, save :: fluc_bezAlign(:)   ! The name of the element
  integer,          allocatable, public, save :: fluc_ixAlign(:)    ! The index of the element
  integer,                       public, save :: fluc_nAlign        ! Number of multipoles in fort.8

  real(kind=fPrec), allocatable, public, save :: fluc_errExt(:,:)   ! The errors from fort.16
  character(len=:), allocatable, public, save :: fluc_bezExt(:)     ! The name of the element
  integer,          allocatable, public, save :: fluc_ixExt(:)      ! The index of the element
  integer,                       public, save :: fluc_nExt          ! Number of multipoles in fort.16

  integer,                       public, save :: fluc_mRead         ! Flag determining what files to read

  integer, public, save :: fluc_iSeed1
  integer, public, save :: fluc_iSeed2

contains

subroutine fluc_parseInputLine(inLine, iLine, iErr)

  use string_tools
  use parpro,         only : nmac
  use mod_common,     only : izu0,mmac,mcut
  use mod_settings,   only : st_debug
  use sixtrack_input, only : sixin_echoVal

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

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

  fluc_mRead = 0

  if(nSplit > 0) call chr_cast(lnSplit(1),izu0,      iErr)
  if(nSplit > 1) call chr_cast(lnSplit(2),mmac,      iErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),fluc_mRead,iErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),mcut,      iErr)

  if(st_debug) then
    call sixin_echoVal("izu0",izu0,      "FLUC",iLine)
    call sixin_echoVal("mmac",mmac,      "FLUC",iLine)
    call sixin_echoVal("mout",fluc_mRead,"FLUC",iLine)
    call sixin_echoVal("mcut",mcut,      "FLUC",iLine)
  end if

  ! Process variables
  mcut        = iabs(mcut)
  fluc_iSeed1 = izu0
  fluc_iSeed2 = 0

  if(mmac > nmac) then
    write(lout,"(a,i0)") "FLUC> ERROR Maximum number of seeds for vectorisation is ",nmac
    iErr = .true.
    return
  end if

  if(fluc_mRead < 0 .or. fluc_mRead > 7) then
    write(lout,"(a,i0)") "FLUC> ERROR I/O options (mout) must be a number between 0 and 7, got ",fluc_mRead
    iErr = .true.
    return
  end if

  if(mcut == 0) then
    write(lout,"(a)")      "FLUC> No cut on random distribution."
  else
    write(lout,"(a,i0,a)") "FLUC> Random distribution has been cut to ",mcut," sigma."
  end if

end subroutine fluc_parseInputLine

subroutine fluc_readInputs

  use mod_common, only : mout2

  implicit none

  select case(fluc_mRead)
  case(1)
    call fluc_readFort16
  case(2)
    mout2=1
  case(3)
    call fluc_readFort16
    mout2=1
  case(4)
    call fluc_readFort8
  case(5)
    call fluc_readFort8
    call fluc_readFort16
  case(6)
    mout2=1
    call fluc_readFort8
  case(7)
    mout2=1
    call fluc_readFort16
    call fluc_readFort8
  end select

  call fluc_moreRandomness

end subroutine fluc_readInputs

subroutine fluc_moreRandomness

  use mod_alloc
  use mod_ranecu
  use parpro,              only : nzfz
  use mod_common,          only : zfz,mcut
  use numerical_constants, only : zero

  implicit none

  integer, parameter :: newRnd = 5000
  real(kind=fPrec)   :: tmpRnd(newRnd)

  call recuin(fluc_iSeed1, fluc_iSeed2)
  call ranecu(tmpRnd, newRnd, mcut)
  call recuut(fluc_iSeed1, fluc_iSeed2)

  if(nzfz == -1) nzfz = 0
  call alloc(zfz, nzfz+newRnd, zero, "zfz")
  zfz(nzfz+1:nzfz+newRnd) = tmpRnd(1:newRnd)
  nzfz = nzfz + newRnd

end subroutine fluc_moreRandomness

subroutine fluc_randomReport

  use parpro,              only : str_divLine,nzfz
  use mod_common,          only : izu0,zfz
  use numerical_constants, only : zero

  implicit none

  real(kind=fPrec) rSum, rSqSum, rMean, rDev
  integer i

  rSum   = zero
  rSqSum = zero

  do i=1,nzfz
    rSum=rSum + zfz(i)
  end do

  rMean = rSum/real(nzfz,fPrec)
  do i=1,nzfz
    rSqSum = rSqSum + (zfz(i)-rMean)**2
  end do

  rDev = sqrt(rSqSum/real(nzfz,fPrec))

  write(lout,"(a)")       ""
  write(lout,"(a)")       "    FLUC Block Random Numbers (zfz):"
  write(lout,"(a,i0)")    "     * Generator Seed    = ",izu0
  write(lout,"(a,i0)")    "     * Numbers Generated = ",nzfz
  write(lout,"(a,f15.7)") "     * Mean Value        = ",rMean
  write(lout,"(a,f15.7)") "     * Deviation         = ",rDev
  write(lout,"(a)")       ""
  write(lout,"(a)")       str_divLine

end subroutine fluc_randomReport

subroutine fluc_readFort8

  use mod_alloc
  use mod_units
  use string_tools
  use parpro,              only : mNameLen,str_nmSpace
  use numerical_constants, only : zero
  use mod_common,          only : il,bez,icextal,nblz,nblo,ic
  use mod_settings,        only : st_debug

  implicit none

  character(len=1024) :: inLine
  character(len=:), allocatable :: lnSplit(:)
  real(kind=fPrec) alignx, alignz, tilt
  integer lineNo, ioStat, nSplit, mAlign, nAlign, iStru, iSing
  logical iErr, isOpen, inSing
  integer i, j

  lineNo      = 0
  fluc_nAlign = 0

  inquire(unit=8, opened=isOpen)
  if(isOpen) close(8)
  call units_openUnits(unit=8,fileName="fort.8",formatted=.true.,mode="r",err=iErr)
  if(iErr) then
    write(lout,"(a)") "FLUC> ERROR Failed to open fort.8"
    goto 30
  end if
  rewind(8)

  write(lout,"(a)") "FLUC> Reading alignment errors from fort.8"

10 continue
  read(8,"(a)",end=20,iostat=ioStat) inLine
  lineNo = lineNo + 1
  if(ioStat /= 0) then
    write(lout,"(a,i0)") "FLUC> ERROR fort.16 iostat = ",ioStat
    goto 30
  end if

  call chr_split(inLine, lnSplit, nSplit, iErr)
  if(iErr) goto 30
  if(nSplit == 0) goto 10

  if(allocated(fluc_ixAlign)) then
    mAlign = size(fluc_ixAlign,1)
  else
    mAlign = 0
  end if

  fluc_nAlign = fluc_nAlign + 1
  if(fluc_nAlign > mAlign) then
    call alloc(fluc_errAlign,3,       mAlign+50,zero,       "fluc_errAlign")
    call alloc(fluc_bezAlign,mNameLen,mAlign+50,str_nmSpace,"fluc_bezAlign")
    call alloc(fluc_ixAlign,          mAlign+50,0,          "fluc_ixAlign")
  end if

  if(nSplit > 0) fluc_bezAlign(fluc_nAlign) = trim(lnSplit(1))
  if(nSplit > 1) call chr_cast(lnSplit(2),fluc_errAlign(1,fluc_nAlign),iErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),fluc_errAlign(2,fluc_nAlign),iErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),fluc_errAlign(3,fluc_nAlign),iErr)

  goto 10

20 continue
  close(8)
  if(fluc_nAlign == 0) then
    write(lout,"(a)") "FLUC> Reading of fort.8 requested in FLUC block, but no elements read."
    return
  end if

  ! We require that the elements in fort.8 have the same order as in the STRUcture block.
  ! A repeated element name implies it is the next one of that name in the sequenc, so we will
  ! iterate over the index array ic() and check for a match in order to build the index icext.
  nAlign = 1
  do iStru=1,nblz
    iSing = ic(iStru)
    if(iSing <= nblo) cycle
    iSing = iSing - nblo
    if(bez(iSing) == fluc_bezAlign(nAlign)) then
      fluc_ixAlign(nAlign) = iStru
      icextal(iStru)       = nAlign
      nAlign = nAlign + 1
      if(nAlign > fluc_nAlign) exit
    end if
  end do
  if(nAlign /= fluc_nAlign+1) then
    write(lout,"(a)")       "FLUC> ERROR Did not find all the elements in fort.8 in the structure."
    write(lout,"(a)")       "FLUC>       You either have a non-existing element somewhere in the file,"
    write(lout,"(a)")       "FLUC>       or too many references to the same element name."
    write(lout,"(2(a,i0))") "FLUC>       Found ",(nAlign-1)," elements out of ",fluc_nAlign
    call prror(-1)
  end if

  write(lout,"(a,i0,a)") "FLUC> Read ",fluc_nAlign," values from fort.8"
  ! if(st_debug) then
  !   do i=1,fluc_nAlign
  !     write(lout,"(a,i5,2(a,i0))") "FLUC> ",i,": '"//trim(fluc_bezAlign(i))//"' "//&
  !       "with ID ",fluc_ixAlign(i)," = ",icextal(fluc_ixAlign(i))
  !     write(lout,"(a,3(2x,e14.7))")  "FLUC> Values: ",fluc_errAlign(1,i),fluc_errAlign(2,i),fluc_errAlign(3,i)
  !   end do
  ! end if

  return

30 continue
  write(lout,"(a,i0,a)") "FLUC> ERROR fort.8:",lineNo," '"//trim(inLine)//"'"
  call prror(-1)

end subroutine fluc_readFort8

subroutine fluc_readFort16

  use mod_alloc
  use mod_units
  use string_tools
  use parpro,              only : mNameLen,str_nmSpace
  use numerical_constants, only : zero
  use mod_common,          only : il,bez,icext,nblz,nblo,ic
  use mod_settings,        only : st_debug

  implicit none

  character(len=1024) :: inLine
  character(len=:), allocatable :: lnSplit(:)
  character(len=mNameLen) bezExt
  integer mType, lineNo, ioStat, nSplit, lMode, nVals, mVal, mExt, iStru, iSing, nExt
  logical iErr, isOpen, inSing
  integer i, j

  mType  = 0
  lineNo = 0
  lMode  = 0

  fluc_nExt = 0

  inquire(unit=16, opened=isOpen)
  if(isOpen) close(16)
  call units_openUnits(unit=16,fileName="fort.16",formatted=.true.,mode="r",err=iErr)
  if(iErr) then
    write(lout,"(a)") "FLUC> ERROR Failed to open fort.16"
    goto 30
  end if
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
        call alloc(fluc_errExt,40,      mExt+50,zero,       "fluc_errExt")
        call alloc(fluc_bezExt,mNameLen,mExt+50,str_nmSpace,"fluc_bezExt")
        call alloc(fluc_ixExt,          mExt+50,0,          "fluc_ixExt")
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
  if(fluc_nExt == 0) then
    write(lout,"(a)") "FLUC> Reading of fort.16 requested in FLUC block, but no elements read."
    return
  end if

  ! We require that the elements in fort.16 have the same order as in the STRUcture block.
  ! A repeated element name implies it is the next one of that name in the sequenc, so we will
  ! iterate over the index array ic() and check for a match in order to build the index icext.
  nExt = 1
  do iStru=1,nblz
    iSing = ic(iStru)
    if(iSing <= nblo) cycle
    iSing = iSing - nblo
    if(bez(iSing) == fluc_bezExt(nExt)) then
      fluc_ixExt(nExt) = iStru
      icext(iStru)     = nExt
      nExt = nExt + 1
      if(nExt > fluc_nExt) exit
    end if
  end do

  if(nExt /= fluc_nExt+1) then
    write(lout,"(a)")       "FLUC> ERROR Did not find all the elements in fort.16 in the structure."
    write(lout,"(a)")       "FLUC>       You either have a non-existing element somewhere in the file,"
    write(lout,"(a)")       "FLUC>       or too many references to the same element name."
    write(lout,"(2(a,i0))") "FLUC>       Found ",(nExt-1)," elements out of ",fluc_nExt
    call prror(-1)
  end if

  write(lout,"(a,i0,a)") "FLUC> Read ",fluc_nExt," values from fort.16"
  ! if(st_debug) then
  !   do i=1,fluc_nExt
  !     write(lout,"(a,i5,a,i0)") "FLUC> ",i,": '"//trim(fluc_bezExt(i))//"' with ID ",fluc_ixExt(i)
  !     do j=1,40,4
  !       write(lout,"(a,i2,a,4(2x,e14.7))")  "FLUC>  * ",j,": ",&
  !         fluc_errExt(j,i),fluc_errExt(j+1,i),fluc_errExt(j+2,i),fluc_errExt(j+3,i)
  !     end do
  !   end do
  ! end if

  return

30 continue
  write(lout,"(a,i0,a)") "FLUC> ERROR fort.16:",lineNo," '"//trim(inLine)//"'"
  call prror(-1)

end subroutine fluc_readFort16

! The code for reading fort.30 has not been implementet, as it seemed to be out of date anyway.
! Pasting it here for reference in case it needs to be put back in.
! -------------------------------------------------------------------------------------------------
! izu=0
! iexnum=0
! if(mout4.eq.1) then
!   read(30,10020,end=1591)
!   rewind 30
!   do 1590 i=1,mper*mbloz
!     ix=ic(i)
!     if(ix.gt.nblo) then
!       ix=ix-nblo
!       kpz=kp(ix)
!       kzz=kz(ix)
!       if(kpz.eq.6.or.kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 1590
!       if(kzz.eq.15) goto 1590
!       izu=izu+3
!       read(30,10020,end=1591,iostat=ierro) ch
!       if(ierro.gt.0) call prror(87)
!       lineno30=lineno30+1
!       call intepr(1,1,ch,ch1)
!       read(ch1,*) ilm0(1),zfz(izu-2)
!       iexnum=iexnum+1
!       if(kz(ix).eq.11) izu=izu+2*mmul
!     endif
! 1590   continue
!   if(iexnum.gt.0) then
!     write(lout,*)
!     write(lout,*)'          Single (random) kick errors read in from external file'
!     write(lout,*)
!     write(lout,*) '        From file fort.30 :',iexnum,' values read in.'
!     write(lout,*)
!   endif
!   iexread=0
!   ifiend8=0
!   iexnum=0
!   rewind 30
!   do 1593 i=1,mper*mbloz
!     ix=ic(i)
!     if(ix.gt.nblo) then
!       ix=ix-nblo
!       if(iexread.eq.0) then
! 1595         ilm0(1)=' '
! ! READ IN HORIZONTAL AND VERTICAL MISALIGNMENT AND TILT
!         if(ifiend8.eq.0) then
!           read(30,10020,end=1594,iostat=ierro) ch
!           if(ierro.gt.0) call prror(87)
!           lineno30=lineno30+1
!         else
!           goto 1594
!         endif
!         call intepr(1,1,ch,ch1)
!         read(ch1,*) ilm0(1),dummy,alignx,alignz,tilt
!         if(((abs(alignx)+abs(alignz))+abs(tilt)).le.pieni) goto 1595
!         iexnum=iexnum+1
!         bezext(iexnum)=ilm0(1)
!         iexread=1
!         goto 1596
! 1594         ifiend8=1
!         do 1597 j=1,iexnum
!           if(bez(ix).eq.bezext(j)) call prror(87)
! 1597         continue
! 1596         continue
!       endif
!       if(ilm0(1).eq.bez(ix)) then
!         icextal(i)=ix
!         extalign(i,1)=alignx
!         extalign(i,2)=alignz
!         extalign(i,3)=tilt
!         iexread=0
!         goto 1593
!       endif
!     endif
! 1593   continue
! 1591   continue
! endif

end module mod_fluc
