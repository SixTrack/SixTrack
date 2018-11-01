! ================================================================================================ !
!  SixTrack Particles Module
!  V.K. Berglyd Olsen, K.N. Sjobak, BE-ABP-HSS
!  Last modified: 2018-08-12
! ================================================================================================ !
module mod_particles

  use crcoall
  use floatPrecision

  implicit none

contains

subroutine part_allocate
end subroutine part_allocate

subroutine part_expand
end subroutine part_expand

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-10-31
!  Moved from maincr. Applies the closed orbit correction to the particle coordinates.
! ================================================================================================ !
subroutine part_applyClosedOrbit

  use mod_common
  use mod_commons
  use mod_commonmn

  implicit none

  integer j

  if(iclo6 == 2) then
    xv(1,1:napx)  = xv(1,1:napx)  +  clo6v(1,1:napx)
    yv(1,1:napx)  = yv(1,1:napx)  + clop6v(1,1:napx)
    xv(2,1:napx)  = xv(2,1:napx)  +  clo6v(2,1:napx)
    yv(2,1:napx)  = yv(2,1:napx)  + clop6v(2,1:napx)
    sigmv(1:napx) = sigmv(1:napx) +  clo6v(3,1:napx)
    dpsv(1:napx)  = dpsv(1:napx)  + clop6v(3,1:napx)
  else if(idfor == 0) then
    xv(1,1:napx)  = xv(1,1:napx)  +   clov(1,1:napx) * real(idz(1),fPrec)
    yv(1,1:napx)  = yv(1,1:napx)  +  clopv(1,1:napx) * real(idz(1),fPrec)
    xv(2,1:napx)  = xv(2,1:napx)  +   clov(2,1:napx) * real(idz(2),fPrec)
    yv(2,1:napx)  = yv(2,1:napx)  +  clopv(2,1:napx) * real(idz(2),fPrec)
  end if
  call part_updateEnergy(e0,3)

end subroutine part_applyClosedOrbit

! ================================================================================================ !
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-13
!  Updates the relevant particle arrays after the particle energies or the reference energy changed.
!  If only the ejv array has been changed, pass e0 to refEnergy
! ================================================================================================ !
subroutine part_updateEnergy(refEnergy, refArray)

  use mod_hions
  use mod_common
  use mod_commont
  use mod_commonmn
  use numerical_constants
  use physical_constants

  implicit none

  real(kind=fPrec), intent(in) :: refEnergy
  integer,          intent(in) :: refArray

  real(kind=fPrec) e0o, e0fo
  integer          j

  ! Modify the reference particle
  e0o  = e0
  e0fo = e0f
  if(e0 /= refEnergy) then
    e0     = refEnergy
    e0f    = sqrt(e0**2 - nucm0**2)
    gammar = nucm0/e0
  end if

  if(e0 <= pieni) then
    write(lout,"(a)") "PART> ERROR Reference energy ~= 0"
    call prror(-1)
  end if

  select case(refArray)
  case(1) ! Update from energy array
    do j=1, napx
      ejfv(j) = sqrt(ejv(j)**2 - nucm(j)**2)        ! Momentum [MeV/c]
      dpsv(j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f   ! Delta_p/p0 = delta
    end do
  case(2) ! Update from momentum array
    do j=1, napx
      ejv(j)  = sqrt(ejfv(j)**2 + nucm(j)**2)       ! Energy [MeV]
      dpsv(j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f   ! Delta_p/p0 = delta
    end do
  case(3) ! Update from delta array
    do j=1, napx
      ejfv(j) = ((nucm(j)/nucm0)*(dpsv(j)+one))*e0f ! Momentum [MeV/c]
      ejv(j)  = sqrt(ejfv(j)**2 + nucm(j)**2)       ! Energy [MeV]
    end do
  case default
    write(lout,"(a)") "PART> ERROR Internal error in part_updateEnergy"
    call prror(-1)
  end select

  ! Modify the Energy Dependent Arrays
  do j=1, napx
    dpsv1(j)    = (dpsv(j)*c1e3)/(one + dpsv(j))
    dpd(j)      = one + dpsv(j)                      ! For thick tracking
    dpsq(j)     = sqrt(dpd(j))                       ! For thick tracking
    oidpsv(j)   = one/(one + dpsv(j))
    moidpsv(j)  = mtc(j)/(one + dpsv(j))             ! Relative rigidity offset (mod_hions) [MV/c^2]
    omoidpsv(j) = ((one-mtc(j))*oidpsv(j))*c1e3
    rvv(j)      = (ejv(j)*e0f)/(e0*ejfv(j))          ! Beta_0 / beta(j)
  end do

  if(e0 /= e0o) then
    ! Also update sigmv with the new beta0 = e0f/e0
    sigmv = ((e0f*e0o)/(e0fo*e0))*sigmv
  end if

  if(ithick == 1) call synuthck

end subroutine part_updateEnergy

end module mod_particles
