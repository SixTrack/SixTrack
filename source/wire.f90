module wire

  use floatPrecision
  use parpro, only : nele, nblz
  use mod_alloc

  implicit none

  ! A. Patapenka (NIU), M. Fitterer (FNAL)
  ! Common block for wire definition
  ! variables to save wire parameters for tracking etc.

  ! wire current [A]
  real(kind=fPrec), allocatable, save :: wire_current(:) !(nele)
  ! integrated length of the wire [m]
  real(kind=fPrec), allocatable, save :: wire_lint(:) !(nele)
  ! physical length of the wire [m]
  real(kind=fPrec), allocatable, save :: wire_lphys(:) !(nele)

  ! integer to include or not closed orbit in the separation between beam and wire
  ! 0  : Un-initialized if wire element not found
  ! +1 : dispx is the distance between x0=y0=0 and the wire
  ! -1 : dispx is the distance between the closed orbit and the wire
  !    x=y=0    <->   xco     <->    xwire
  !               closed orbit    wire position
  ! wire_flagco = +1: dispx = xwire -> rx = x + xsep
  ! wire_flagco = -1: dispx = xwire - xco -> rx = x - xco + xsep
  ! -> rx = x + xwire
  integer, allocatable, save          :: wire_flagco(:) !(nele)
  ! hor./vert. displacement of the wire [mm]
  real(kind=fPrec), allocatable, save :: wire_dispx(:),wire_dispy(:) !(nele)
  ! hor./vert. tilt of the wire [degrees] -90 < tilty < 90,
  !  uses the same definition as the DISP block
  real(kind=fPrec), allocatable, save :: wire_tiltx(:), wire_tilty(:) !(nele)


  ! wire parameters for closed orbit calculation (FOX part)
  ! for FOX length of variable names must be smaller 8
  integer, parameter :: wire_max = 350 ! max. number of wires (same as BB interactions)
  real(kind=fPrec), save :: wire_clo(6,wire_max) ! closed orbit at wire

  ! wire number for each structure element (default = 0 if no wire)
  integer, allocatable, save :: wire_num(:) ! (nblz)

contains

subroutine wire_expand_arrays(nele_new, nblz_new)

  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: nblz_new

  call alloc(wire_current,nele_new,zero,'wire_current')
  call alloc(wire_lint,nele_new,zero,'wire_lint')
  call alloc(wire_lphys,nele_new,zero,'wire_lphys')

  call alloc(wire_flagco,nele_new,0,'wire_flagco')
  call alloc(wire_num,nblz_new,0,"wire_num")

  call alloc(wire_dispx,nele_new,zero,'wire_dispx')
  call alloc(wire_dispy,nele_new,zero,'wire_dispy')
  call alloc(wire_tiltx,nele_new,zero,'wire_tiltx')
  call alloc(wire_tilty,nele_new,zero,'wire_tilty')

end subroutine wire_expand_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Here Tobias will put the code

end module wire
