! ================================================================================================ !
!  SixTrack BEAM Module
!  V.K. Berglyd Olsen, K.N. Sjobak, BE-ABP-HSS
!  Last modified: 2018-08-12
! ================================================================================================ !
module mod_beam

  use crcoall
  use floatPrecision
  use numerical_constants
  use physical_constants

  implicit none

contains

subroutine beam_allocate
end subroutine beam_allocate

subroutine beam_expand
end subroutine beam_expand

! ================================================================================================ !
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-08-12
!  Updates the relevant particle arrays after the particle energies or the reference energy changed.
!  If only the ejfv array has been changed, pass e0 to refEnergy
! ================================================================================================ !
subroutine beam_updateParticleEnergy(refEnergy)

  use mod_hions
  use mod_common
  use mod_commont
  use mod_commonmn

  implicit none

  real(kind=fPrec), intent(in) :: refEnergy

  real(kind=fPrec) e0o, e0fo
  integer          j

  ! Modify the reference particle
  e0o    = e0
  e0fo   = e0f
  e0     = refEnergy
  e0f    = sqrt(e0**2 - nucm0**2)
  gammar = nucm0/e0

  ! Modify the Energy
  do j=1, napx
    dpsv(j)    = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
    dpsv1(j)   = (dpsv(j)*c1e3)/(one + dpsv(j))
    dpd(j)     = one + dpsv(j)
    dpsq(j)    = sqrt(dpd(j))
    oidpsv(j)  = one/(one + dpsv(j))
    moidpsv(j) = mtc(j)/(one + dpsv(j))
    rvv(j)     = (ejv(j)*e0f)/(e0*ejfv(j))

    ! Also update sigmv with the new beta0 = e0f/e0
    if(e0 /= e0o) sigmv(j) = ((e0f*e0o)/(e0fo*e0))*sigmv(j)
  end do
  if(ithick == 1) call synuthck

end subroutine beam_updateParticleEnergy

end module mod_beam
