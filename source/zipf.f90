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

  character(len=:), allocatable, save :: zipf_fileNames(:)  ! Name of files to pack into the zip file.
  character(len=mFileName),      save :: zipf_outFile = " " ! Name of output file (Default: Sixout.zip)

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
  integer nSplit
  logical spErr

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lerr,"(a)") "ZIPF> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  if(nSplit /= 1) then
    write(lerr,"(a,i3,3a)") "ZIPF> ERROR Expected 1 filename per line, got ",nSplit
    iErr = .true.
    return
  end if

  zipf_numFiles = zipf_numFiles + 1
  call alloc(zipf_fileNames, mFileName, zipf_numFiles, " ", "zipf_fileNames")
  zipf_fileNames(zipf_numFiles) = trim(lnSplit(1))

end subroutine zipf_parseInputLine

subroutine zipf_parseInputDone(iErr)

  use crcoall
  use string_tools

  logical, intent(inout) :: iErr

  integer ii

  zipf_outFile = "Sixout.zip" ! Output name fixed for now
  write(lout,"(a)")    "ZIPF> Output file name = '"//trim(zipf_outFile)//"'"
  write(lout,"(a,i0)") "ZIPF> Number of files to pack = ",zipf_numFiles
  write(lout,"(a)")    "ZIPF> Files:"
  do ii=1,zipf_numFiles
    write(lout,"(a,i5,a)") "ZIPF>  * ",ii,": '"//trim(zipf_fileNames(ii))//"'"
  end do

  if(.not.(zipf_numFiles > 0)) then
    write(lerr,"(a)") "ZIPF> ERROR No files specified in input block"
    iErr = .true.
    return
  end if

end subroutine zipf_parseInputDone

subroutine zipf_dozip

  use crcoall
  use mod_alloc
  use string_tools

  integer iErr
#ifdef BOINC
  integer i
  character(len=256)               boincZipFile
  character(len=:), allocatable :: boincNames(:)

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

#ifdef LIBARCHIVE
  write(lout,"(a)") "LIBARCHIVE> Compressing file '"//trim(zipf_outFile)//"' ..."
  ! The f_write_archive function will handle the conversion from Fortran to C-style strings
  ! NOTE: Last two arguments of the C function are implicitly passed from FORTRAN, there is no need to do it explicitly.
#ifdef BOINC
  call f_write_archive(trim(boincZipFile),boincNames,zipf_numFiles)
#else
  call f_write_archive(trim(zipf_outFile),zipf_fileNames,zipf_numFiles)
#endif
#elif ZLIB
#ifdef BOINC
  call minizip_zip(trim(boincZipFile),boincNames(1),zipf_numFiles,9,iErr,len_trim(boincZipFile),256)
#else
  call minizip_zip(trim(zipf_outFile),zipf_fileNames(1),zipf_numFiles,9,iErr,len_trim(zipf_outFile),mFileName)
#endif
  if(iErr /= 0) then
    write(lerr,"(a,i0)") "ZIPF> ERROR MiniZip returned error code ",iErr
    call prror
  end if
#else
  write(lout,"(a)") "ZIPF> *** No Archive Libraries in this SixTrack *** "
#endif

  write(lout,"(a)") "ZIPF> Done!"

end subroutine zipf_dozip

end module zipf
