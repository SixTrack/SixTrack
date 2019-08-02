module coll_jaw

  use floatPrecision
  use numerical_constants, only : zero, one

  implicit none

  type, private :: jaw_fitType
    character(len=:), allocatable :: fitName                ! Fit name
    real(kind=fPrec)              :: fitParams(6) = zero    ! Fit parameters
    real(kind=fPrec)              :: fitScale     = one     ! Fit scale
    logical                       :: reCentre     = .false. ! Recentre the collimator to new minimum
  end type jaw_fitType

  type, private :: jaw_collType
    integer                       :: nSlices      = 0       ! The number of slices
    real(kind=fPrec)              :: fitLength    = zero    ! Slice length
    real(kind=fPrec), allocatable :: fitData(:,:)           ! Slice data
    real(kind=fPrec), allocatable :: fitTilt(:,:)           ! Slice tilt
  end type jaw_collType

  type(jaw_fitType),  allocatable, private, save :: jaw_fitData(:)  ! List of jaw fit parameters
  type(jaw_collType), allocatable, private, save :: jaw_collData(:) ! List of collimator fit slices
  integer,                         private, save :: jaw_nFits  = 0  ! Count of jaw_fitData
  integer,                         private, save :: jaw_nColls = 0  ! Count of jaw_collData

contains

integer function jaw_getSliceCount(sliceID)
  integer, intent(in) :: sliceID
  jaw_getSliceCount = jaw_collData(sliceID)%nSlices
end function jaw_getSliceCount

subroutine jaw_addJawFit(fitName, fitParams, fitScale, reCentre, fitID)

  character(len=*), intent(in)  :: fitName
  real(kind=fPrec), intent(in)  :: fitParams(6)
  real(kind=fPrec), intent(in)  :: fitScale
  logical,          intent(in)  :: reCentre
  integer,          intent(out) :: fitID

  type(jaw_fitType), allocatable :: jawTmp(:)

  if(allocated(jaw_fitData)) then
    allocate(jawTmp(jaw_nFits+1))
    jawTmp(1:jaw_nFits) = jaw_fitData(1:jaw_nFits)
    call move_alloc(jawTmp, jaw_fitData)
    jaw_nFits = jaw_nFits + 1
  else
    allocate(jaw_fitData(1))
    jaw_nFits = 1
  end if

  jaw_fitData(jaw_nFits)%fitName   = fitName
  jaw_fitData(jaw_nFits)%fitParams = fitParams
  jaw_fitData(jaw_nFits)%fitScale  = fitScale
  jaw_fitData(jaw_nFits)%reCentre  = reCentre

  fitID = jaw_nFits

end subroutine jaw_addJawFit

subroutine jaw_computeFit(collName, fitID, nSlices, cLength, cTilt, cOffset, sliceID)

  use crcoall
  use mathlib_bouncer
  use numerical_constants

  character(len=*), intent(in)  :: collName
  integer,          intent(in)  :: fitID(2)
  integer,          intent(in)  :: nSlices
  real(kind=fPrec), intent(in)  :: cLength
  real(kind=fPrec), intent(in)  :: cTilt(2)
  real(kind=fPrec), intent(in)  :: cOffset
  integer,          intent(out) :: sliceID

  type(jaw_collType), allocatable :: collTmp(:)
  real(kind=fPrec) sX(nSlices+1), sX1(nSlices+1), sX2(nSlices+1), sY1(nSlices+1), sY2(nSlices+1)
  real(kind=fPrec) angle1(nSlices), angle2(nSlices)
  real(kind=fPrec) fac1(6), fac2(6), scale1, scale2, maxY
  real(kind=fPrec) fitData(2,nSlices+1), fitTilt(2,nSlices+1)
  integer i

  if(allocated(jaw_collData)) then
    allocate(collTmp(jaw_nColls+1))
    collTmp(1:jaw_nColls) = jaw_collData(1:jaw_nColls)
    call move_alloc(collTmp, jaw_collData)
    jaw_nColls = jaw_nColls + 1
  else
    allocate(jaw_collData(1))
    jaw_nColls = 1
  end if

  fac1   = jaw_fitData(fitID(1))%fitParams
  fac2   = jaw_fitData(fitID(2))%fitParams
  scale1 = jaw_fitData(fitID(1))%fitScale
  scale2 = jaw_fitData(fitID(2))%fitScale

  ! Calculate longitudinal positions of slices and corresponding heights and angles from the fit parameters.
  ! Note: here, take (nSlices+1) points in order to calculate the tilt angle of the last slice!

  do i=1,nSlices+1
    ! Deformation of the jaws scaled with length
    sX(i)  = (real(i-1,fPrec)*cLength)/real(nSlices,fPrec)
    sY1(i) = fac1(1) + fac1(2)*sX(i) + (fac1(3)/cLength)*sX(i)**2 + fac1(4)*sX(i)**3 + fac1(5)*sX(i)**4 + fac1(6)*sX(i)**5
    sY2(i) = fac2(1) + fac2(2)*sX(i) + (fac2(3)/cLength)*sX(i)**2 + fac2(4)*sX(i)**3 + fac2(5)*sX(i)**4 + fac2(6)*sX(i)**5

    ! Apply the slicing scaling factors (ssf's):
    sY1(i) =       scale1 * sY1(i)
    sY2(i) = -one*(scale2 * sY2(i))

    ! Coordinates rotated of the tilt
    sX1(i) =  sX(i)*cos_mb(cTilt(1)) - sY1(i)*sin_mb(cTilt(1))
    sX2(i) =  sX(i)*cos_mb(cTilt(2)) - sY2(i)*sin_mb(cTilt(2))
    sY1(i) = sY1(i)*cos_mb(cTilt(1)) +  sX(i)*sin_mb(cTilt(1))
    sY2(i) = sY2(i)*cos_mb(cTilt(2)) +  sX(i)*sin_mb(cTilt(2))
  end do

  do i=1,nSlices
    ! Sign of the angle defined differently for the two jaws
    angle1(i) = ((sY1(i+1) - sY1(i)) / (sX1(i+1)-sX1(i)))
    angle2(i) = ((sY2(i+1) - sY2(i)) / (sX2(i+1)-sX2(i)))
  end do

  ! For both jaws, look for the 'deepest' point (closest point to beam)
  ! Then, shift the vectors such that this closest point defines the nominal aperture
  ! Index here must go up to (nSlices+1) in case the last point is the
  ! closest (and also for the later calculation of 'a_tmp1' and 'a_tmp2')
  ! The recentring flag, as given in 'fort.3' chooses whether we recentre the deepest point or not
  if(jaw_fitData(fitID(1))%reCentre) then
    maxY = c1e6
    do i=1,nSlices+1
      if(sY1(i) < maxY) then
        maxY = sY1(i)
      end if
    end do
  else
    maxY = zero
  end if
  do i=1,nSlices+1
    sY1(i) = sY1(i) - maxY ! + (half*c_aperture)
  end do

  if(jaw_fitData(fitID(2))%reCentre) then
    maxY = -c1e6
    do i=1,nSlices+1
      if(sY2(i) > maxY) then
        maxY = sY2(i)
      end if
    end do
  else
    maxY = zero
  end if
  do i=1,nSlices+1
    sY2(i) = sY2(i) - maxY ! - (half*c_aperture)
  end do

  fitData(1,:) = sY1(:)
  fitData(2,:) = sY2(:)
  fitTilt(1,:) = angle1(:)
  fitTilt(2,:) = angle2(:)

  jaw_collData(jaw_nColls)%nSlices   = nSlices
  jaw_collData(jaw_nColls)%fitLength = cLength/real(nSlices,fPrec)
  jaw_collData(jaw_nColls)%fitData   = fitData
  jaw_collData(jaw_nColls)%fitTilt   = fitTilt

  sliceID = jaw_nColls

  write(lout,"(a,i0,a)") "COLLJAW> Collimator "//trim(collName)//" sliced in ",nSlices," slices"

end subroutine jaw_computeFit

subroutine jaw_getFitData(sliceID, nSlice, cLength, cAperture, cOffset, cTilt)

  use numerical_constants

  integer,          intent(in)    :: sliceID
  integer,          intent(in)    :: nSlice
  real(kind=fPrec), intent(out)   :: cLength
  real(kind=fPrec), intent(inout) :: cAperture
  real(kind=fPrec), intent(inout) :: cOffset
  real(kind=fPrec), intent(out)   :: cTilt(2)

  real(kind=fPrec) val1, val2

  cLength  = jaw_collData(sliceID)%fitLength
  cTilt(:) = jaw_collData(sliceID)%fitTilt(:,nSlice)

  if(cTilt(1) > zero) then
    val1 = jaw_collData(sliceID)%fitData(1,nSlice)
  else
    val1 = jaw_collData(sliceID)%fitData(1,nSlice+1)
  end if

  if(cTilt(2) < zero) then
    val2 = jaw_collData(sliceID)%fitData(2,nSlice)
  else
    val2 = jaw_collData(sliceID)%fitData(2,nSlice+1)
  end if

  val1 = val1 + half*cAperture
  val2 = val2 - half*cAperture

  cAperture = val1 - val2
  cOffset   = half*(val1 + val2) + cOffset

end subroutine jaw_getFitData

end module coll_jaw
