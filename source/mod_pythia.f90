! ================================================================================================ !
!
!  SixTrack - Pythia8 Interface Module
!
!  V.K. Berglyd Olsen, BE-ABP-HSS, 2018
!
! ================================================================================================ !

#ifdef PYTHIA

module pythia

  use, intrinsic :: iso_c_binding
  use crcoall

  implicit none

contains

subroutine pythia_init(partType) bind(C, name="pythiaWrapper_init")
  integer(kind=C_INT), value, intent(in) :: partType
end subroutine pythia_init

end module pythia

#endif
