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

  real(kind=fPrec), allocatable, public, save :: fluc_errZFZ(:,:)   ! The alignemnt errors from fort.30
  character(len=:), allocatable, public, save :: fluc_bezZFZ(:)     ! The name of the element
  integer,          allocatable, public, save :: fluc_ixZFZ(:)      ! The index of the element
  integer,                       public, save :: fluc_nZFZ          ! Number of multipoles in fort.30

  integer,                       public, save :: fluc_mRead         ! Flag determining what files to read

  integer, public, save :: fluc_iSeed1
  integer, public, save :: fluc_iSeed2

contains

subroutine fluc_parseInputLine(inLine, iLine, iErr)

  use string_tools
  use mod_common,     only : izu0,mcut
  use mod_settings,   only : st_debug
  use sixtrack_input, only : sixin_echoVal

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, iDummy
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
  if(nSplit > 1) call chr_cast(lnSplit(2),iDummy,    iErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),fluc_mRead,iErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),mcut,      iErr)

  if(st_debug) then
    call sixin_echoVal("izu0",izu0,      "FLUC",iLine)
    call sixin_echoVal("mout",fluc_mRead,"FLUC",iLine)
    call sixin_echoVal("mcut",mcut,      "FLUC",iLine)
  end if

  ! Process variables
  mcut        = iabs(mcut)
  fluc_iSeed1 = izu0
  fluc_iSeed2 = 0

  if(iDummy > 1) then
    write(lout,"(a,i0)") "FLUC> ERROR Multiple seeds for vectorisation (mmac) is no longer supported. You requested ",iDummy
    iErr = .true.
    return
  end if

  if(fluc_mRead < 0 .or. fluc_mRead > 15) then
    write(lout,"(a,i0)") "FLUC> ERROR I/O options (mout) must be a number between 0 and 15, got ",fluc_mRead
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

  fluc_nExt   = 0 ! For fort.16
  fluc_nAlign = 0 ! For fort.8
  fluc_nZFZ   = 0 ! For fort.30

  if(iand(fluc_mRead, 1) == 1) call fluc_readFort16
  if(iand(fluc_mRead, 2) == 2) mout2 = 1
  if(iand(fluc_mRead, 4) == 4) call fluc_readFort8
  if(iand(fluc_mRead, 8) == 8) call fluc_readFort30

  call fluc_moreRandomness

end subroutine fluc_readInputs

subroutine fluc_moreRandomness

  use mod_alloc
  use mod_ranecu
  use parpro,              only : nzfz
  use mod_common,          only : zfz,mcut
  use numerical_constants, only : zero

  implicit none

  integer            :: i
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

! ================================================================================================ !
!  Reading alignment errors from fort.8
!  Last modified: 2018-06-15
! ================================================================================================ !
subroutine fluc_readFort8

  use mod_alloc
  use mod_units
  use string_tools
  use parpro,              only : mNameLen,str_nmSpace
  use numerical_constants, only : zero
  use mod_common,          only : il,bez,icextal,nblz,nblo,ic

  implicit none

  character(len=1024) :: inLine
  character(len=:), allocatable :: lnSplit(:)
  real(kind=fPrec) alignx, alignz, tilt
  integer lineNo, ioStat, nSplit, mAlign, nAlign, iStru, iSing
  logical iErr, isFile, inSing
  integer i

  lineNo = 0
  inLine = " "

  inquire(file="fort.8", exist=isFile)
  if(.not.isFile) return

  call f_close(8)

  call f_open(unit=8,file="fort.8",formatted=.true.,mode="r",err=iErr)
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
    write(lout,"(a,i0)") "FLUC> ERROR fort.8 iostat = ",ioStat
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
  call f_close(8)
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

  return

30 continue
  write(lout,"(a,i0,a)") "FLUC> ERROR fort.8:",lineNo," '"//trim(inLine)//"'"
  call prror(-1)

end subroutine fluc_readFort8

! ================================================================================================ !
!  Reading multipole errors from fort.16
!  Last modified: 2018-06-14
! ================================================================================================ !
subroutine fluc_readFort16

  use mod_alloc
  use mod_units
  use string_tools
  use parpro,              only : mNameLen,str_nmSpace
  use numerical_constants, only : zero
  use mod_common,          only : il,bez,icext,nblz,nblo,ic

  implicit none

  character(len=1024) :: inLine
  character(len=:), allocatable :: lnSplit(:)
  character(len=mNameLen) bezExt
  integer mType, lineNo, ioStat, nSplit, lMode, nVals, mVal, mExt, iStru, iSing, nExt
  logical iErr, isFile, inSing
  integer i

  mVal   = 0
  mType  = 0
  lineNo = 0
  lMode  = 0
  inLine = " "

  inquire(file="fort.16", exist=isFile)
  if(.not.isFile) return

  call f_close(16)

  call f_open(unit=16,file="fort.16",formatted=.true.,mode="r",err=iErr)
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
  call f_close(16)
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

  return

30 continue
  write(lout,"(a,i0,a)") "FLUC> ERROR fort.16:",lineNo," '"//trim(inLine)//"'"
  call prror(-1)

end subroutine fluc_readFort16

! ================================================================================================ !
!  Reading single random kick errors from fort.30
!  Last modified: 2018-06-18
!  Note: There were no tests in the test suite for reading these file when the code was moved here
!    from daten. The implementation is therefore not tested. The values read are appended to the
!    end of the arrays used for alignement errors from fort.8 and to an additional array for the
!    values to overwrite in zfz.
! ================================================================================================ !
subroutine fluc_readFort30

  use mod_alloc
  use mod_units
  use string_tools
  use parpro,              only : mNameLen,str_nmSpace,mmul
  use numerical_constants, only : zero
  use mod_common,          only : il,bez,icextal,nblz,nblo,ic,kp,kz

  implicit none

  character(len=1024) :: inLine
  character(len=:), allocatable :: lnSplit(:)
  real(kind=fPrec) alignx, alignz, tilt
  integer lineNo, ioStat, nSplit, mZFZ, nZFZ, iStru, iSing, iZ
  logical iErr, isFile, inSing
  integer i

  lineNo = 0
  inLine = " "

  inquire(file="fort.30", exist=isFile)
  if(.not.isFile) return

  call f_close(30)

  call f_open(unit=30,file="fort.30",formatted=.true.,mode="r",err=iErr)
  if(iErr) then
    write(lout,"(a)") "FLUC> ERROR Failed to open fort.30"
    goto 30
  end if
  rewind(30)

  write(lout,"(a)") "FLUC> Reading single (random) kick errors from fort.30"

10 continue
  read(30,"(a)",end=20,iostat=ioStat) inLine
  lineNo = lineNo + 1
  if(ioStat /= 0) then
    write(lout,"(a,i0)") "FLUC> ERROR fort.30 iostat = ",ioStat
    goto 30
  end if

  call chr_split(inLine, lnSplit, nSplit, iErr)
  if(iErr) goto 30
  if(nSplit == 0) goto 10

  if(allocated(fluc_ixZFZ)) then
    mZFZ = size(fluc_ixZFZ,1)
  else
    mZFZ = 0
  end if

  fluc_nZFZ = fluc_nZFZ + 1
  if(fluc_nZFZ > mZFZ) then
    call alloc(fluc_errZFZ,4,       mZFZ+50,zero,       "fluc_errZFZ")
    call alloc(fluc_bezZFZ,mNameLen,mZFZ+50,str_nmSpace,"fluc_bezZFZ")
    call alloc(fluc_ixZFZ,          mZFZ+50,0,          "fluc_ixZFZ")
  end if

  if(nSplit > 0) fluc_bezZFZ(fluc_nZFZ) = trim(lnSplit(1))
  if(nSplit > 1) call chr_cast(lnSplit(2),fluc_errZFZ(1,fluc_nZFZ),iErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),fluc_errZFZ(2,fluc_nZFZ),iErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),fluc_errZFZ(3,fluc_nZFZ),iErr)
  if(nSplit > 4) call chr_cast(lnSplit(5),fluc_errZFZ(4,fluc_nZFZ),iErr)

  goto 10

20 continue
  call f_close(30)
  if(fluc_nZFZ == 0) then
    write(lout,"(a)") "FLUC> Reading of fort.30 requested in FLUC block, but no elements read."
    return
  end if

  ! We require that the elements in fort.30 have the same order as in the STRUcture block.
  ! A repeated element name implies it is the next one of that name in the sequenc, so we will
  ! iterate over the index array ic() and check for a match in order to build the index icext.
  iZ   = 0
  nZFZ = 1
  do iStru=1,nblz
    iSing = ic(iStru)
    if(iSing <= nblo) cycle
    iSing = iSing - nblo
    if(kp(iSing) == 6)  cycle
    if(kz(iSing) == 0)  cycle
    if(kz(iSing) == 15) cycle
    if(kz(iSing) == 20) cycle
    if(kz(iSing) == 22) cycle
    if(bez(iSing) == fluc_bezZFZ(nZFZ)) then
      fluc_ixZFZ(nZFZ) = iStru
      icextal(iStru)   = -nZFZ ! Note: This overwrites the mapping index from the fort.8 array to the fort.30 negative index
      iZ = iZ + 3
      fluc_ixZFZ(nZFZ) = iZ - 2
      nZFZ = nZFZ + 1
      if(kz(iSing) == 11) iZ = iZ + 2*mmul
      if(nZFZ > fluc_nZFZ) exit
    end if
  end do
  if(nZFZ /= fluc_nZFZ+1) then
    write(lout,"(a)")       "FLUC> ERROR Did not find all the elements in fort.30 in the structure."
    write(lout,"(a)")       "FLUC>       You either have a non-existing element somewhere in the file,"
    write(lout,"(a)")       "FLUC>       or too many references to the same element name."
    write(lout,"(2(a,i0))") "FLUC>       Found ",(nZFZ-1)," elements out of ",fluc_nZFZ
    call prror(-1)
  end if

  write(lout,"(a,i0,a)") "FLUC> Read ",fluc_nZFZ," values from fort.30"

  return

30 continue
  write(lout,"(a,i0,a)") "FLUC> ERROR fort.30:",lineNo," '"//trim(inLine)//"'"
  call prror(-1)

end subroutine fluc_readFort30

! ================================================================================================ !
!  Write modified geometry file fort.4
!  Last modified: 2018-06-17
! ================================================================================================ !
subroutine fluc_writeFort4

  use floatPrecision
  use mod_common,          only : ncororb,sm,ek
  use numerical_constants, only : zero,pieni
  use string_tools
  use mod_units

  implicit none

  character(len=:), allocatable :: lnSplit(:)
  character(len=1024) inLine
  character(len=mNameLen) elemName
  real(kind=fPrec) inVal(6)
  integer          lineNo2, nSplit, ii, iKZ
  logical          iErr

  ii      = 0
  lineNo2 = 0

  call f_open(unit=2,file="fort.2",formatted=.true.,mode="r",err=fErr)
  call f_open(unit=4,file="fort.4",formatted=.true.,mode="w",err=iErr)

  rewind(2)
10 continue
  read(2,"(a)",end=90) inLine
  lineNo2 = lineNo2 + 1
  if(inLine(:1) == "/" .or. inLine(:1) == "!") then
    write(4,"(a)") trim(inLine)
    goto 10
  elseif(inLine(:4) == "SING") then
    write(4,"(a)") trim(inLine)
  else
    return
  end if

20 continue
  read(2,"(a)",end=90) inLine
  lineNo2 = lineNo2 + 1
  if(inLine(:1) == "/" .or. inLine(:1) == "!") then
    write(4,"(a)") trim(inLine)
  else
    ii = ii + 1
    if(inLine(:4) /= "NEXT") then
      iKZ      = 0
      inVal(:) = zero
      call chr_split(inLine,lnSplit,nSplit,iErr)
      if(nSplit == 0) goto 20
      if(nSplit > 0) elemName = trim(lnSplit(1))
      if(nSplit > 1) call chr_cast(lnSplit(2),iKZ,     iErr)
      if(nSplit > 2) call chr_cast(lnSplit(3),inVal(1),iErr)
      if(nSplit > 3) call chr_cast(lnSplit(4),inVal(2),iErr)
      if(nSplit > 4) call chr_cast(lnSplit(5),inVal(3),iErr)
      if(nSplit > 5) call chr_cast(lnSplit(6),inVal(4),iErr)
      if(nSplit > 6) call chr_cast(lnSplit(7),inVal(5),iErr)
      if(nSplit > 7) call chr_cast(lnSplit(8),inVal(6),iErr)
      if(ikz == 11) then
        write(4,"(a48,1x,i2,6(1x,e22.15))") elemName,iKZ,inVal(1),inVal(2),inVal(3),inVal(4),inVal(5),inVal(6)
      else
        if(abs(inVal(3)) <= pieni) then
          if(ncororb(ii) == 0) then
            write(4,"(a48,1x,i2,6(1x,e22.15))") elemName,iKZ,sm(ii),inVal(2),inVal(3),inVal(4),inVal(5),inVal(6)
          else
            write(4,"(a48,1x,i2,6(1x,e22.15))") elemName,iKZ,sm(ii),ek(ii),inVal(3),inVal(4),inVal(5),inVal(6)
          end if
        else
          write(4,"(a48,1x,i2,6(1x,e22.15))") elemName,iKZ,inVal(1),ek(ii),inVal(3),inVal(4),inVal(5),inVal(6)
        end if
      end if
    else
      write(4,"(a)") trim(inLine)
      goto 30
    end if
  end if
  goto 20

30 continue
  read(2,"(a)",end=90) inLine
  lineNo2 = lineNo2 + 1
  write(4,"(a)") trim(inLine)
  goto 30

90 continue
  call f_close(2)
  call f_close(4)
  return
end subroutine fluc_writeFort4

end module mod_fluc
