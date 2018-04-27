! ================================================================================================ !
!  ZIPF Module
! ~~~~~~~~~~~~~
!  K.N. Sjobak, BE-ABP-HSS
!  Last Modified: 2018-04-27
!  Compress selected output files into a zip file at the end of the simulation
! ================================================================================================ !
module zipf
  
  use end_sixtrack
  use string_tools
  
  implicit none
  
  integer, save :: zipf_numFiles
  
  character(len=:), allocatable, save :: zipf_fileNames(:) ! Name of files to pack into the zip file.
  character(len=str_maxLen),     save :: zipf_outFile      ! Name of output file (Default: Sixout.zip)
  
contains

subroutine zipf_parseInputLine(ch)
  
  use crcoall
  use mod_alloc
  
  implicit none
  
  character(len=*), intent(in) :: ch
  
  character gFields(str_maxFields)*(str_maxLen) ! Array of fields
  integer   nFields                             ! Number of identified fields
  integer   lFields(str_maxFields)              ! Length of each what:
  logical   errFields                           ! An error flag
  
  ! Read filenames
  call getfields_split(ch, gFields, lFields, nFields, errFields)
  if(errFields) call prror(-1)
  
  if (nFields /= 1) then
    write(lout,"(a)")       "ERROR in ZIPF:"
    write(lout,"(a,i3,3a)") "Expected 1 filename per line, got ",nFields,", line = '",ch,"'"
    call prror(-1)
  end if
  
  zipf_numFiles = zipf_numFiles + 1
  call resize(zipf_fileNames,str_maxLen, zipf_numFiles, str_dZeros, "zipf_fileNames")
  zipf_fileNames(zipf_numFiles)(1:lFields(1)) = gFields(1)(1:lFields(1))
  
end subroutine zipf_parseInputLine

subroutine zipf_parseInputDone
  
  use crcoall
  use string_tools
  
  implicit none
  
  integer ii
  
  zipf_outFile(1:10) = "Sixout.zip" ! Output name fixed for now
  write(lout,"(a)")     "**** ZIPF ****"
  write(lout,"(3a)")    " Output file name = '",trim(chr_trimZero(zipf_outFile)),"'"
  write(lout,"(a,i5)")  " Number of files to pack = ",zipf_numFiles
  write(lout,"(a)")     " Files:"
  write(lout,"(i5,2a)") (ii,": ",trim(chr_trimZero(zipf_fileNames(ii))),ii=1,zipf_numFiles)
  
  if(.not.(zipf_numFiles > 0)) then
    write(lout,"(a)") "ERROR in ZIPF:"
    write(lout,"(a)") " ZIPF block was empty;"
    write(lout,"(a)") " no files specified!"
    call prror(-1)
  endif
  
#ifndef LIBARCHIVE
  write(lout,"(a)") "ERROR in ZIPF:"
  write(lout,"(a)") " ZIPF needs LIBARCHIVE to work,"
  write(lout,"(a)") " but this SixTrack was compiled without it."
  call prror(-1)
#endif
  
end subroutine zipf_parseInputDone

subroutine zipf_comnul
  
  use mod_alloc
  
  implicit none
  
  zipf_numFiles = 0
  zipf_outFile  = str_dZeros
  
  call alloc(zipf_fileNames, str_maxLen, 1, str_dZeros, "zipf_fileNames")
  
end subroutine zipf_comnul

subroutine zipf_dozip
  
  use crcoall
  
  implicit none
  
#ifdef BOINC
    character(stringzerotrim_maxlen) zipf_outFile_boinc
    character(stringzerotrim_maxlen) zipf_fileNames_boinc(zipf_maxfiles)
    integer ii
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
  
  write(lout,"(3a)") "ZIPF: Compressing file '", trim(chr_trimZero(zipf_outFile)),"'..."
  
#ifdef LIBARCHIVE
#ifdef BOINC
  ! For BOINC, we may need to translate the filenames.
  call boincrf(trim(stringzerotrim(zipf_outFile)), zipf_outFile_boinc )
  
  do ii=1,zipf_numFiles
    call boincrf(trim(chr_trimZero(zipf_fileNames(ii))), zipf_fileNames_boinc(ii))
    zipf_fileNames_boinc(ii) = trim(zipf_fileNames_boinc(ii))
  end do
  
  ! The f_write_archive function will handle the conversion from Fortran to C-style strings
  call f_write_archive(trim(zipf_outFile_boinc),zipf_fileNames_boinc,zipf_numFiles)
#else
  call f_write_archive(zipf_outFile,zipf_fileNames,zipf_numFiles)
#endif
#else
  ! If not libarchive, the zipf subroutine shall just be a stub.
  ! And anyway daten should not accept the block, so this is somewhat redundant.
  write(lout,"(a)") " *** No libArchive in this SixTrack *** "
#endif
  
  write(lout,"(a)") "Done!"
  
end subroutine zipf_dozip

end module zipf
