! ================================================================================================ !
!  ZIPF Module
! ~~~~~~~~~~~~~
!  K.N. Sjobak, BE-ABP-HSS
!  Last Modified: 2018-04-27
!  Compress selected output files into a zip file at the end of the simulation
! ================================================================================================ !
module zipf

  use parpro

  implicit none

  integer, save :: zipf_numFiles = 0

  character(len=:), allocatable, private, save :: zipf_fileNames(:)           ! Name of files to pack into the zip file.
  character(len=mFileName),      private, save :: zipf_outFile = "Sixout.zip" ! Name of output file (Default: Sixout.zip)
  integer,                       private, save :: zipf_zipLevel = 3           ! Compression level, 0-9

  interface
    subroutine minizip_zip(zipFile, inFiles, nFiles, compLevel, iErr, lenOne, lenTwo) bind(C, name="minizip_zip")
      use, intrinsic :: iso_c_binding
      character(kind=C_CHAR,len=1), intent(in)  :: zipFile
      character(kind=C_CHAR,len=1), intent(in)  :: inFiles
      integer(kind=C_INT), value,   intent(in)  :: nFiles
      integer(kind=C_INT), value,   intent(in)  :: compLevel
      integer(kind=C_INT),          intent(out) :: iErr
      integer(kind=C_INT), value,   intent(in)  :: lenOne
      integer(kind=C_INT), value,   intent(in)  :: lenTwo
    end subroutine minizip_zip

    subroutine minizip_unzip(zipFile, toDir, iErr, lenOne, lenTwo) bind(C, name="minizip_unzip")
      use, intrinsic :: iso_c_binding
      character(kind=C_CHAR,len=1), intent(in)  :: zipFile
      character(kind=C_CHAR,len=1), intent(in)  :: toDir
      integer(kind=C_INT),          intent(out) :: iErr
      integer(kind=C_INT), value,   intent(in)  :: lenOne
      integer(kind=C_INT), value,   intent(in)  :: lenTwo
    end subroutine minizip_unzip
  end interface

contains

subroutine zipf_parseInputLine(inLine,iErr)

  use crcoall
  use mod_alloc
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable :: lnSplit(:)
  integer nSplit, i
  logical spErr

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lerr,"(a)") "ZIPF> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(lnSplit(1))

  case("OUTFILE")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "ZIPF> ERROR OUTFILE takes 1 argument, got ",nSplit-1
      write(lerr,"(a)")    "ZIPF>       OUTFILE zipFileName"
      iErr = .true.
      return
    end if
    zipf_outFile = lnSplit(2)

  case("ZIPLEVEL")
    if(nSplit /= 2) then
      write(lerr,"(a,i0)") "ZIPF> ERROR ZIPLEVEL level takes 1 input parameter, got ",nSplit-1
      write(lerr,"(a)")    "ZIPF>       ZIPLEVEL 0-9"
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2), zipf_zipLevel, iErr)
    if(zipf_zipLevel < 0 .or. zipf_zipLevel > 9) then
      write(lerr,"(a,i0)") "ZIPF> ERROR ZIPLEVEL level must be between 0 and 9, got ",zipf_zipLevel
      iErr = .true.
      return
    end if

  case default

    call alloc(zipf_fileNames, mFileName, zipf_numFiles + nSplit, " ", "zipf_fileNames")
    do i=1,nSplit
      zipf_fileNames(zipf_numFiles + i) = trim(lnSplit(i))
    end do
    zipf_numFiles = zipf_numFiles + nSplit

  end select

end subroutine zipf_parseInputLine

subroutine zipf_dozip

  use crcoall
  use mod_alloc
  use string_tools

  integer iErr
#ifdef BOINC
  integer i
  character(len=256)               boincZipFile
  character(len=:), allocatable :: boincNames(:)
#endif

#ifdef BOINC
  call alloc(boincNames, 256, zipf_numFiles, " ", "boincNames")
#endif

  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine
  write(lout,"(a)") ""

#ifdef BOINC
  ! For BOINC, we may need to translate the filenames.
  call boincrf(trim(zipf_outFile), boincZipFile)
  do i=1,zipf_numFiles
    call boincrf(trim(zipf_fileNames(i)), boincNames(i))
    boincNames(i) = trim(boincNames(i))
  end do
#endif

#ifdef ZLIB
#ifdef BOINC
  call minizip_zip(trim(boincZipFile),boincNames(1),zipf_numFiles,zipf_zipLevel,iErr,len_trim(boincZipFile),256)
#else
  call minizip_zip(trim(zipf_outFile),zipf_fileNames(1),zipf_numFiles,zipf_zipLevel,iErr,len_trim(zipf_outFile),mFileName)
#endif
  if(iErr /= 0) then
    write(lerr,"(a,i0)") "ZIPF> WARNING MiniZip returned error code ",iErr
  end if
#else
  write(lout,"(a)") "ZIPF> *** No ZLIB in this SixTrack *** "
#endif

  write(lout,"(a)") "ZIPF> Done!"

end subroutine zipf_dozip

end module zipf
