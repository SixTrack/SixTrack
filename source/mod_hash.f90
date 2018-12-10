! ================================================================================================ !
!  HASH MODULE
!  V.k. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-30
!
!  This module provides an interface to the MD5 implementation written by Ronald L. Rivest (MIT).
!  The source code is available under source/md5.
! ================================================================================================ !
module mod_hash

  use floatPrecision

  implicit none

  character(len=:), allocatable, private, save :: hash_listHashFiles(:)
  logical,          allocatable, private, save :: hash_isAscii(:)
  integer,                       private, save :: hash_nHashFiles  =  0
  logical,                       private, save :: hash_selfTestOK  = .false.
  character(len=8),              parameter     :: hash_sumFileName = "hash.md5"

  ! C Interface
  interface

    ! Direct interfaces to the MD5 main functions.
    subroutine hash_md5Init(nInst) bind(C, name="md5wrapper_md5Init")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value, intent(in) :: nInst
    end subroutine hash_md5Init

    subroutine hash_md5Update(ctxID, inStr, strLen) bind(C, name="md5wrapper_md5Update")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value,   intent(in) :: ctxID
      character(kind=C_CHAR,len=1), intent(in) :: inStr
      integer(kind=C_INT), value,   intent(in) :: strLen
    end subroutine hash_md5Update

    subroutine hash_md5FinalC(ctxID, md5Vals, md5Size) bind(C, name="md5wrapper_md5Final")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value, intent(in)    :: ctxID
      integer(kind=C_INT),        intent(inout) :: md5Vals(*)
      integer(kind=C_INT), value, intent(in)    :: md5Size
    end subroutine hash_md5FinalC

    ! Interfaces to complete digest functions
    subroutine hash_digestStringC(inStr, strLen, md5Vals, md5Size) bind(C, name="md5wrapper_digestString")
      use, intrinsic :: iso_c_binding
      character(kind=C_CHAR,len=1), intent(in)    :: inStr
      integer(kind=C_INT), value,   intent(in)    :: strLen
      integer(kind=C_INT),          intent(inout) :: md5Vals(*)
      integer(kind=C_INT), value,   intent(in)    :: md5Size
    end subroutine hash_digestStringC

    subroutine hash_digestFileC(fileName, strLen, md5Vals, md5Size, isAscii) bind(C, name="md5wrapper_digestFile")
      use, intrinsic :: iso_c_binding
      character(kind=C_CHAR,len=1), intent(in)    :: fileName
      integer(kind=C_INT), value,   intent(in)    :: strLen
      integer(kind=C_INT),          intent(inout) :: md5Vals(*)
      integer(kind=C_INT), value,   intent(in)    :: md5Size
      integer(kind=C_INT), value,   intent(in)    :: isAscii
    end subroutine hash_digestFileC

  end interface

contains

! ================================================================================================ !
!  INIT THE HASH MODULE
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-30
!  Performs a self test to check that it complies with the RFC1321 standard.
! ================================================================================================ !
subroutine hash_initialise

  use crcoall

  character(len=32) md5Digest, md5Valid
  character(len=6)  tStatus

  write(lout,"(a)") "HASH> Library Self Test:"
  hash_selfTestOK = .true.

  call hash_digestString("",md5Digest)
  md5Valid = "d41d8cd98f00b204e9800998ecf8427e"
  if(md5Valid == md5Digest) then
    tStatus = "Passed"
  else
    tStatus = "Failed"
    hash_selfTestOK = .false.
  end if
  write(lout,"(a)") "HASH>   #1 : "//md5Digest//" < ''"
  write(lout,"(a)") "HASH>    1 : "//md5Valid//" > "//tStatus

  call hash_digestString("message digest",md5Digest)
  md5Valid = "f96b697d7cb7938d525a2f31aaf161d0"
  if(md5Valid == md5Digest) then
    tStatus = "Passed"
  else
    tStatus = "Failed"
    hash_selfTestOK = .false.
  end if
  write(lout,"(a)") "HASH>   #2 : "//md5Digest//" < 'message digest'"
  write(lout,"(a)") "HASH>    2 : "//md5Valid//" > "//tStatus

  call hash_digestString("abcdefghijklmnopqrstuvwxyz",md5Digest)
  md5Valid = "c3fcd3d76192e4007dfb496cca67e13b"
  if(md5Valid == md5Digest) then
    tStatus = "Passed"
  else
    tStatus = "Failed"
    hash_selfTestOK = .false.
  end if
  write(lout,"(a)") "HASH>   #3 : "//md5Digest//" < 'abcdefghijklmnopqrstuvwxyz'"
  write(lout,"(a)") "HASH>    3 : "//md5Valid//" > "//tStatus

end subroutine hash_initialise

! ================================================================================================ !
!  INPUT LINE PARSING
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-30
! ================================================================================================ !
subroutine hash_parseInputLine(inLine, iErr)

  use crcoall
  use string_tools
  use mod_alloc

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable :: lnSplit(:)
  integer nSplit
  logical spErr, cErr
  logical tmpIsAscii

  cErr = .false.

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lout,"(a)") "HASH> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(trim(lnSplit(1)))

  case("MD5SUM")
    if(nSplit /= 2 .and. nSplit /= 3) then
      write(lout,"(a,i3)") "HASH> ERROR MD5SUM expected 1 or 2 arguments, got ",nSplit-1
      iErr = .true.
      return
    end if
    if (nSplit == 3) then
      call chr_cast(trim(lnSplit(3)), tmpIsAscii, cErr)
      if (cErr) then
        write(lout,"(a)") "HASH> ERROR MD5SUM optional argument expected a logical variable, but '"// &
             trim(lnSplit(3))//"' didn't match."
        call prror
      endif
    else
      tmpIsAscii = .false.
    endif
    if(len_trim(lnSplit(2)) > 255) then
      write(lout,"(a)") "HASH> ERROR MD5SUM filename is too long. Max is 255."
      iErr = .true.
      return
    end if

    hash_nHashFiles = hash_nHashFiles + 1
    call alloc(hash_listHashFiles, 255, hash_nHashFiles,     " ", "hash_listHashFiles")
    call alloc(hash_isAscii,            hash_nHashFiles, .false., "hash_isAscii")
    hash_listHashFiles(hash_nHashFiles) = trim(lnSplit(2))
    hash_isAscii(hash_nHashFiles)       = tmpIsAscii
  case default
    write(lout,"(a)") "HASH> ERROR Unknown keyword '"//trim(lnSplit(1))//"'"
    iErr = .true.
    return

  end select

end subroutine hash_parseInputLine

! ================================================================================================ !
!  COMPUTE MD5SUMS
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-12-04
!  Computes the md5 digest of files listed in the HASH block
! ================================================================================================ !
subroutine hash_fileSums

  use parpro
  use crcoall
  use file_units

  character(len=32) md5Digest
  integer           nFile

  integer hash_sumFileUnit
#ifdef WIN32
  integer hash_sumFileUnit_win32
#endif

  if(hash_nHashFiles == 0) return

  call funit_requestUnit(hash_sumFileName, hash_sumFileUnit)
  open(hash_sumFileUnit, file=hash_sumFileName)
#ifdef WIN32
  call funit_requestUnit(hash_sumFileName, hash_sumFileUnit_win32)
  open(hash_sumFileUnit_win32, file=hash_sumFileName//".win32")
#endif

  if(hash_selfTestOK .eqv. .false.) then
    write(hash_sumFileUnit,"(a)") "HASH library failed self test. No hashes written."
    flush(hash_sumFileUnit)
    close(hash_sumFileUnit)
    return
  end if

  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine
  write(lout,"(a)") ""
  write(lout,"(a)") "    Computing MD5 Hash of Files"
  write(lout,"(a)") "   ============================="
  do nFile=1,hash_nHashFiles
    call hash_digestFile(trim(hash_listHashFiles(nFile)), md5Digest, hash_isAscii(nFile))
#ifdef WIN32
    write(hash_sumFileUnit_win32,"(a32,2x,a)") md5Digest,trim(hash_listHashFiles(nFile))//".tmp"
    flush(hash_sumFileUnit_win32)
#endif

    write(hash_sumFileUnit,      "(a32,2x,a)") md5Digest,trim(hash_listHashFiles(nFile))
    write(lout,                  "(a36,2x,a)") md5Digest,trim(hash_listHashFiles(nFile))
    flush(hash_sumFileUnit)
  end do
  close(hash_sumFileUnit)
#ifdef WIN32
  close(hash_sumFileUnit_win32)
#endif

end subroutine hash_fileSums

! ================================================================================================ !
!  Wrapper Subroutines for the Interface
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-30
! ================================================================================================ !
subroutine hash_md5Final(ctxID, md5Digest)

  use, intrinsic :: iso_c_binding
  use string_tools

  integer,           intent(in)  :: ctxID
  character(len=32), intent(out) :: md5Digest

  integer(kind=C_INT) tmpVals(16)
  integer i

  tmpVals(:) = 0
  call hash_md5FinalC(ctxID, tmpVals, 16)
  write(md5Digest,"(16(z2.2))") tmpVals
  md5Digest = chr_toLower(md5Digest)

end subroutine hash_md5Final

subroutine hash_digestString(inStr, md5Digest)

  use, intrinsic :: iso_c_binding
  use string_tools

  character(len=*),  intent(in)  :: inStr
  character(len=32), intent(out) :: md5Digest

  integer(kind=C_INT) tmpVals(16)
  integer i

  tmpVals(:) = 0
  call hash_digestStringC(inStr//char(0), len(inStr), tmpVals, 16)
  write(md5Digest,"(16(z2.2))") tmpVals
  md5Digest = chr_toLower(md5Digest)

end subroutine hash_digestString

subroutine hash_digestFile(fileName, md5Digest, isAscii)

  use, intrinsic :: iso_c_binding
  use string_tools

  character(len=*),  intent(in)  :: fileName
  character(len=32), intent(out) :: md5Digest
  logical,           intent(in)  :: isAscii

  integer tmpIsAscii

  integer(kind=C_INT) tmpVals(16)
  integer i

  if (isAscii) then
    tmpIsAscii = 1
  else
    tmpIsAscii = 0
  end if

  tmpVals(:) = 0
  call hash_digestFileC(fileName//char(0), len(fileName)+1, tmpVals, 16, tmpIsAscii)
  if(tmpVals(1) == -1) then
    md5Digest = "***** ERROR File Not Found *****"
  else
    write(md5Digest,"(16(z2.2))") tmpVals
    md5Digest = chr_toLower(md5Digest)
  end if

end subroutine hash_digestFile

end module mod_hash
