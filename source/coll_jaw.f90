module coll_jaw

  use floatPrecision
  use numerical_constants, only : zero, one

  implicit none

  type, private :: jaw_typeJaw
    integer                       :: nSlices      = 0       ! The number of slices
    real(kind=fPrec)              :: polyFit(2,6) = zero    ! The polynomical fit factors
    real(kind=fPrec)              :: fitScale(2)  = one     ! The scale for jaw 1 and 2
    logical                       :: reCentre(2)  = .false. ! Recentre the collimator to new minimum
    real(kind=fPrec)              :: fitLength    = zero    ! Slice length
    real(kind=fPrec), allocatable :: fitAper(:)             ! Slice aperture
    real(kind=fPrec), allocatable :: fitOffset(:)           ! Slice offset
    real(kind=fPrec), allocatable :: fitTilt(:,:)           ! Slice tilt
  end type jaw_typeJaw

  type(jaw_typeJaw), allocatable, private, save :: jaw_fitData(:) ! The list of jaw fit parameters
  integer,                        private, save :: jaw_nColl = 0  ! The number of collimaters with jaw fits

contains

subroutine jaw_addJawFit(nSlices, polyFit, fitScale, reCentre, jawID)

  integer,          intent(in)  :: nSlices
  real(kind=fPrec), intent(in)  :: polyFit(2,6)
  real(kind=fPrec), intent(in)  :: fitScale(2)
  logical,          intent(in)  :: reCentre(2)
  integer,          intent(out) :: jawID

  type(jaw_typeJaw), allocatable :: jawTmp(:)

  if(allocated(jaw_fitData)) then
    allocate(jawTmp(jaw_nColl+1))
    jawTmp(1:jaw_nColl) = jaw_fitData(1:jaw_nColl)
    call move_alloc(jawTmp, jaw_fitData)
    jaw_nColl = jaw_nColl + 1
  else
    allocate(jaw_fitData(1))
    jaw_nColl = 1
  end if

  jaw_fitData(jaw_nColl)%nSlices  = nSlices
  jaw_fitData(jaw_nColl)%polyFit  = polyFit
  jaw_fitData(jaw_nColl)%fitScale = fitScale
  jaw_fitData(jaw_nColl)%reCentre = reCentre

  jawID = jaw_nColl

end subroutine jaw_addJawFit

end module coll_jaw
