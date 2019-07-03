module mod_dist

implicit none

    ! C Interface
  interface
  subroutine dist_readfile(fileName) bind(C, name="readfile")
 		use, intrinsic :: iso_c_binding
        character(kind=C_CHAR,len=1), intent(in) :: fileName
 end subroutine dist_readfile
    end interface 
end module