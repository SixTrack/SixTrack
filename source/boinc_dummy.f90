subroutine boinc_time_to_checkpoint(timech)
  implicit none
  integer timech
  timech=1
end subroutine boinc_time_to_checkpoint

subroutine boincrf(myname,filename)
  implicit none
  character(*) myname
  character(256) filename
  filename=myname
end subroutine boincrf

subroutine boinc_finish(flag)
  implicit none
  integer flag
end subroutine boinc_finish

subroutine boinc_fraction_done(f)
  implicit none
  double precision f
end subroutine boinc_fraction_done

subroutine boinc_init
  implicit none
end subroutine boinc_init

subroutine boinc_checkpoint_completed
  implicit none
end subroutine boinc_checkpoint_completed

subroutine boinc_zip(mode,zipfile,path,strl1,strl2)
  implicit none
  integer mode,strl1,strl2
  character(*) zipfile,path
end subroutine boinc_zip

subroutine boinc_zipitall
  implicit none
end subroutine boinc_zipitall

subroutine boinc_unzip
  implicit none
end subroutine boinc_unzip
