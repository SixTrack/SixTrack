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
  character(len=mFNameLen),      save :: zipf_outFile = " " ! Name of output file (Default: Sixout.zip)

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
    write(lout,"(a)") "ZIPF> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit /= 1) then
    write(lout,"(a,i3,3a)") "ZIPF> ERROR Expected 1 filename per line, got ",nSplit
    iErr = .true.
    return
  end if

  zipf_numFiles = zipf_numFiles + 1
  call alloc(zipf_fileNames, mFNameLen, zipf_numFiles, " ", "zipf_fileNames")
  zipf_fileNames(zipf_numFiles) = trim(lnSplit(1))

end subroutine zipf_parseInputLine

subroutine zipf_parseInputDone

  use crcoall
  use string_tools

  implicit none

  integer ii

  zipf_outFile = "Sixout.zip" ! Output name fixed for now
  write(lout,"(a)")     "ZIPF> Output file name = '"//trim(zipf_outFile)//"'"
  write(lout,"(a,i0)")  "ZIPF> Number of files to pack = ",zipf_numFiles
  write(lout,"(a)")     "ZIPF> Files:"
  do ii=1,zipf_numFiles
    write(lout,"(a,i5,a)") "ZIPF>  * ",ii,": '"//trim(zipf_fileNames(ii))//"'"
  end do

  if(.not.(zipf_numFiles > 0)) then
    write(lout,"(a)") "ZIPF> ERROR Block was empty; no files specified!"
    call prror(-1)
  endif

end subroutine zipf_parseInputDone

subroutine zipf_dozip

  use crcoall
  use mod_alloc
  use string_tools

  implicit none

#ifdef BOINC
  character(len=256)               zipf_outFile_boinc
  character(len=:), allocatable :: zipf_fileNames_boinc(:)
  integer ii

  call alloc(zipf_fileNames_boinc, 256, zipf_numFiles, " ", "zipf_fileNames_boinc")
#endif

!+if libarchive
! Having an actual explicit interface would be nice - this is 90% there,
! however some logic should probably be changed (here, in libArchive_Fwrapper.c,
! and in the test program). For now, if it passes CTEST, we think it works...
!    interface
!       subroutine f_write_archive(outname,filenames,numfiles,outname_len,filenames_len) &
!            bind(C,name="f_write_archive_")
!         use, intrinsic :: iso_c_binding, only : c_int
!         implicit none
!         character(len=1), intent(in) :: outname
!         character(len=1), intent(in) :: filenames(*)
!         integer(kind=c_int), intent(in) :: numfiles
!         integer(kind=c_int), intent(in), VALUE :: outname_len
!         integer(kind=c_int), intent(in), VALUE :: filenames_len
!       end subroutine f_write_archive
!    end interface
!+ei

  write(lout,"(a)") "ZIPF> Compressing file '"//trim(zipf_outFile)//"' ..."

#ifdef LIBARCHIVE
#ifdef BOINC
  ! For BOINC, we may need to translate the filenames.
  call boincrf(trim(zipf_outFile), zipf_outFile_boinc)

  do ii=1,zipf_numFiles
    call boincrf(trim(zipf_fileNames(ii)), zipf_fileNames_boinc(ii))
    zipf_fileNames_boinc(ii) = trim(zipf_fileNames_boinc(ii))
  end do

  ! The f_write_archive function will handle the conversion from Fortran to C-style strings
  ! NOTE: Last two arguments of the C function are implicitly passed from FORTRAN, there is no need to do it explicitly.
  call f_write_archive(trim(zipf_outFile_boinc),zipf_fileNames_boinc,zipf_numFiles)
#else
  call f_write_archive(trim(zipf_outFile),zipf_fileNames,zipf_numFiles)
#endif
#else
  ! If not libarchive, the zipf subroutine shall just be a stub.
  ! And anyway daten should not accept the block, so this is somewhat redundant.
  write(lout,"(a)") "ZIPF> *** No libArchive in this SixTrack *** "
#endif

  write(lout,"(a)") "ZIPF> Done!"

end subroutine zipf_dozip

end module zipf
