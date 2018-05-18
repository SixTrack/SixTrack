! ================================================================================================ !
!  SixTrack Input
! ~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-18
! ================================================================================================ !
module sixtrack_input
  
  use crcoall
  use end_sixtrack
  use numerical_constants
  use mod_alloc
  use string_tools
  use mod_common
  
  implicit none
  
  ! Single Element Variables
  integer,                       public, save :: stin_ncy2
  character(len=:), allocatable, public, save :: bez0(:)   ! (max_name_len)(nele)
  
contains

subroutine stin_parseInputLineSING(inLine, iElem)
  
  use parpro_scale
  
  implicit none
  
  character(len=*), intent(in)  :: inLine
  integer,          intent(in)  :: iElem
  
  character(len=:), allocatable :: lnSplit(:)
  character(len=:), allocatable :: elemName
  integer nSplit
  
  integer i
  
  call chr_split(inLine, lnSplit, nSplit)
  
  ! write(lout,"(a,i0,a)") "SING> Line(",nSplit,"): '"//chr_trimZero(inLine)//"'"
  ! do i=1,nSplit
  !   write(lout,"(a,i3,a)") "SING> ",i,": '"//chr_trimZero(lnSplit(i))//"'"
  ! end do
  
  if(nSplit < 2) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element line must be at least 2 values, got ",nSplit
    call prror(-1)
  end if
  
  elemName = chr_trimZero(lnSplit(1))
  if(len(elemName) > max_name_len) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element name too long. Max length is ",max_name_len
    call prror(-1)
  end if
  
  ! Check that the name is unique
  do i=1,iElem-1
    if(bez(i) == elemName) then
      write(lout,"(a,i0)") "GEOMETRY> ERROR Single element '"//elemName//"' is not unique."
      call prror(-1)
    end if
  end do
  
  ! Save Values
  if(nSplit > 1) kz(iElem)   = chr_toInt(chr_trimZero(lnSplit(2)))
  if(nSplit > 2) ed(iElem)   = chr_toReal(chr_trimZero(lnSplit(3)))
  if(nSplit > 3) ek(iElem)   = chr_toReal(chr_trimZero(lnSplit(4)))
  if(nSplit > 4) el(iElem)   = chr_toReal(chr_trimZero(lnSplit(5)))
  if(nSplit > 5) bbbx(iElem) = chr_toReal(chr_trimZero(lnSplit(6)))
  if(nSplit > 6) bbby(iElem) = chr_toReal(chr_trimZero(lnSplit(7)))
  if(nSplit > 7) bbbs(iElem) = chr_toReal(chr_trimZero(lnSplit(8)))
 
  if(kz(iElem) == 25) then
    ed(iElem) = ed(iElem)/two
    ek(iElem) = ek(iElem)/two
  endif
  
  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  if((kz(iElem) == 4 .or. kz(iElem) == 5) .and. abs(el(iElem)) > pieni) then
    ed(iElem) = -one*ed(iElem)
  end if
  
  ! CAVITIES
  if(abs(kz(iElem)) == 12) then
    if(abs(ed(iElem)) > pieni .and. abs(ek(iElem)) > pieni) then
      stin_ncy2     = stin_ncy2 + 1
      itionc(iElem) = kz(iElem)/abs(kz(iElem))
      kp(iElem)     = 6
    end if
  end if
  
  !--------------------------------------------
  ! Handled by initialize_element subroutine:
  !--------------------------------------------
  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  ! THIN LENS (+/- 1-10)
  ! MULTIPOLES (11)
  ! CAVITY (+/- 12)
  ! CRABCAVITY (23/-23) / CC multipoles order 2/3/4 (+/- 23/26/27/28)
  ! ELECTRON LENSE (29)
  call initialize_element(iElem,.true.)
  
  ! ACDIPOLE
  if(abs(kz(iElem)) == 16) then
    if(abs(ed(iElem)) <= pieni) then
      kz(iElem) = 0
      ed(iElem) = zero
      ek(iElem) = zero
      el(iElem) = zero
    else
      acdipph(iElem) = el(iElem)
      el(iElem)      = zero
    end if
  end if
  
  ! General
  if(abs(el(iElem)) > pieni .and. kz(iElem) /= 0) then
    ithick = 1
  end if
  
  ! Expand Arrays
  if(iElem  > nele-2) then
    call expand_arrays(nele+50, npart, nblz, nblo)
    ! The value of nele will have been updated here
    call resize(bez0, max_name_len, nele, repeat(char(0),max_name_len), "bez0")
  end if
  
  if(abs(kz(iElem)) /= 12 .or. (abs(kz(iElem)) == 12 .and. stin_ncy2 == 0)) then
    kp(iElem) = 0
  end if
  
  bez(iElem)  = elemName
  bez0(iElem) = elemName
  
end subroutine stin_parseInputLineSING

end module sixtrack_input
