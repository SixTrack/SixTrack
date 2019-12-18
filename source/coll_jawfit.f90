! ================================================================================================ !
!  Collimation Jaw Fit Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-01
!  Updated: 2019-09-10
! ================================================================================================ !
module coll_jawfit

  use floatPrecision
  use numerical_constants, only : zero, one

  implicit none

  integer, public, parameter :: jaw_fitNameLen = 16

  type, private :: type_fitParams
    character(len=jaw_fitNameLen) :: fitName                ! Fit name
    real(kind=fPrec)              :: fitParams(6) = zero    ! Fit parameters
  end type type_fitParams

  type, private :: type_sliceParams
    integer                       :: nSlices      = 0       ! The number of slices
    real(kind=fPrec)              :: fitLength    = zero    ! Slice length
    real(kind=fPrec), allocatable :: fitData(:,:)           ! Slice data
    real(kind=fPrec), allocatable :: fitTilt(:,:)           ! Slice tilt
  end type type_sliceParams

  type(type_fitParams),   allocatable, private, save :: jaw_fitData(:)     ! List of jaw fit parameters
  type(type_sliceParams), allocatable, private, save :: jaw_sliceData(:)   ! List of collimator fit slices
  integer,                             private, save :: jaw_nFitData   = 0 ! Count of jaw_fitData
  integer,                             private, save :: jaw_nSliceData = 0 ! Count of jaw_sliceData

contains

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-01
!  Updated: 2019-09-10
!  Add a jaw fit model to the storage
! ================================================================================================ !
subroutine jaw_addJawFit(fitName, fitParams, fitID, fErr)

  use crcoall

  character(len=*), intent(in)    :: fitName
  real(kind=fPrec), intent(in)    :: fitParams(6)
  integer,          intent(out)   :: fitID
  logical,          intent(inout) :: fErr

  type(type_fitParams), allocatable :: jawTmp(:)

  if(allocated(jaw_fitData)) then
    allocate(jawTmp(jaw_nFitData+1))
    jawTmp(1:jaw_nFitData) = jaw_fitData(1:jaw_nFitData)
    call move_alloc(jawTmp, jaw_fitData)
    jaw_nFitData = jaw_nFitData + 1
  else
    allocate(jaw_fitData(1))
    jaw_nFitData = 1
  end if

  fitID = jaw_getFitID(fitName)
  if(fitID /= -1) then
    write(lerr,"(a)") "COLLJAW> ERROR Duplicate jaw fit name '"//trim(fitName)//"'"
    fErr = .true.
    return
  end if

  jaw_fitData(jaw_nFitData)%fitName   = fitName
  jaw_fitData(jaw_nFitData)%fitParams = fitParams

  fitID = jaw_nFitData

  write(lout,"(a)") "COLLJAW> Added jaw fit profile '"//trim(fitName)//"'"

end subroutine jaw_addJawFit

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-01
!  Updated: 2019-09-10
!  Compute the fit parameters for a given collimator using a stored fit model.
!  Most of the fit code was extracted from the main collimation module.
!  Original comments from original authors have been preserved.
! ================================================================================================ !
subroutine jaw_computeFit(collName, fitID, nSlices, fitScale, reCentre, cLength, cTilt, cOffset, sliceID)

  use crcoall
  use mathlib_bouncer
  use numerical_constants

  character(len=*), intent(in)  :: collName
  integer,          intent(in)  :: fitID(2)
  integer,          intent(in)  :: nSlices
  real(kind=fPrec), intent(in)  :: fitScale(2)
  logical,          intent(in)  :: reCentre(2)
  real(kind=fPrec), intent(in)  :: cLength
  real(kind=fPrec), intent(in)  :: cTilt(2)
  real(kind=fPrec), intent(in)  :: cOffset
  integer,          intent(out) :: sliceID

  type(type_sliceParams), allocatable :: collTmp(:)
  real(kind=fPrec) sX(nSlices+1), sX1(nSlices+1), sX2(nSlices+1), sY1(nSlices+1), sY2(nSlices+1)
  real(kind=fPrec) angle1(nSlices), angle2(nSlices)
  real(kind=fPrec) fac1(6), fac2(6), scale1, scale2
  real(kind=fPrec) fitData(2,nSlices+1), fitTilt(2,nSlices)
  integer i

  if(allocated(jaw_sliceData)) then
    allocate(collTmp(jaw_nSliceData+1))
    collTmp(1:jaw_nSliceData) = jaw_sliceData(1:jaw_nSliceData)
    call move_alloc(collTmp, jaw_sliceData)
    jaw_nSliceData = jaw_nSliceData + 1
  else
    allocate(jaw_sliceData(1))
    jaw_nSliceData = 1
  end if

  fac1   = jaw_fitData(fitID(1))%fitParams
  fac2   = jaw_fitData(fitID(2))%fitParams
  scale1 = fitScale(1)
  scale2 = fitScale(2)

  sX(:)  = zero
  sX1(:) = zero
  sX2(:) = zero
  sY1(:) = zero
  sY2(:) = zero

  ! Calculate longitudinal positions of slices and corresponding heights and angles from the fit parameters.
  ! Note: here, take (nSlices+1) points in order to calculate the tilt angle of the last slice!

  do i=1,nSlices+1
    ! Deformation of the jaws scaled with length
    sX(i)  = (real(i-1,fPrec)*cLength)/real(nSlices,fPrec)
    sY1(i) = ((((fac1(1) + fac1(2)*sX(i)) + (fac1(3)/cLength)*sX(i)**2) + fac1(4)*sX(i)**3) + fac1(5)*sX(i)**4) + fac1(6)*sX(i)**5
    sY2(i) = ((((fac2(1) + fac2(2)*sX(i)) + (fac2(3)/cLength)*sX(i)**2) + fac2(4)*sX(i)**3) + fac2(5)*sX(i)**4) + fac2(6)*sX(i)**5

    ! Apply the slicing scaling factors
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
  if(reCentre(1)) then
    sY1 = sY1 - minval(sY1)
  end if
  if(reCentre(2)) then
    sY2 = sY2 - maxval(sY2)
  end if

  fitData(1,:) = sY1(:)
  fitData(2,:) = sY2(:)
  fitTilt(1,:) = angle1(:)
  fitTilt(2,:) = angle2(:)

  jaw_sliceData(jaw_nSliceData)%nSlices   = nSlices
  jaw_sliceData(jaw_nSliceData)%fitLength = cLength/real(nSlices,fPrec)
  jaw_sliceData(jaw_nSliceData)%fitData   = fitData
  jaw_sliceData(jaw_nSliceData)%fitTilt   = fitTilt

  sliceID = jaw_nSliceData

  write(lout,"(a,i0,a)") "COLLJAW> Collimator '"//trim(collName)//"' sliced in ",nSlices," slices"

end subroutine jaw_computeFit

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-01
!  Updated: 2019-08-02
!  Return the adjusted length, aperture and offset for a specific collimator slice
! ================================================================================================ !
subroutine jaw_getFitSliceValues(sliceID, iSlice, cLength, cAperture, cOffset, cTilt)

  use numerical_constants

  integer,          intent(in)    :: sliceID
  integer,          intent(in)    :: iSlice
  real(kind=fPrec), intent(out)   :: cLength
  real(kind=fPrec), intent(inout) :: cAperture
  real(kind=fPrec), intent(inout) :: cOffset
  real(kind=fPrec), intent(out)   :: cTilt(2)

  real(kind=fPrec) val1, val2

  cLength  = jaw_sliceData(sliceID)%fitLength
  cTilt(:) = jaw_sliceData(sliceID)%fitTilt(:,iSlice)

  if(cTilt(1) > zero) then
    val1 = jaw_sliceData(sliceID)%fitData(1,iSlice)
  else
    val1 = jaw_sliceData(sliceID)%fitData(1,iSlice+1)
  end if

  if(cTilt(2) < zero) then
    val2 = jaw_sliceData(sliceID)%fitData(2,iSlice)
  else
    val2 = jaw_sliceData(sliceID)%fitData(2,iSlice+1)
  end if

  val1 = val1 + half*cAperture
  val2 = val2 - half*cAperture

  cAperture = val1 - val2
  cOffset   = half*(val1 + val2) + cOffset

end subroutine jaw_getFitSliceValues

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-01
!  Updated: 2019-09-10
!  Return the number of slices for a given jaw fit
! ================================================================================================ !
integer function jaw_getSliceCount(sliceID)
  integer, intent(in) :: sliceID
  if(sliceID < 1 .or. sliceID > jaw_nSliceData) then ! Invalid index
    jaw_getSliceCount = -1
    return
  end if
  jaw_getSliceCount = jaw_sliceData(sliceID)%nSlices
end function jaw_getSliceCount

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-10
!  Updated: 2019-09-10
!  Return the index of a fit by name
! ================================================================================================ !
integer function jaw_getFitID(fitName)

  character(len=*), intent(in) :: fitName

  integer i, fitID

  fitID = -1
  do i=1,jaw_nFitData
    if(jaw_fitData(i)%fitName == fitName) then
      fitID = i
      exit
    end if
  end do
  jaw_getFitID = fitID

end function jaw_getFitID

end module coll_jawfit
