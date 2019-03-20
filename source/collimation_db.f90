! ================================================================================================ !
!  Collimator Database Module
! ================================================================================================ !
module collimation_db

  use parpro
  use floatPrecision

  implicit none

  ! Database arrays
  integer,                       public, save :: cdb_nColl = 0    ! Number of collimators
  character(len=:), allocatable, public, save :: cdb_cName(:)     ! Collimator name
  character(len=:), allocatable, public, save :: cdb_cMaterial(:) ! Collimator material
  integer,          allocatable, public, save :: cdb_cFamily(:)   ! Collimator family
  real(kind=fPrec), allocatable, public, save :: cdb_cNSig(:)     ! Collimator sigma
  real(kind=fPrec), allocatable, public, save :: cdb_cLength(:)   ! Collimator length
  real(kind=fPrec), allocatable, public, save :: cdb_cOffset(:)   ! Collimator offset
  real(kind=fPrec), allocatable, public, save :: cdb_cRotation(:) ! Collimator rotation
  real(kind=fPrec), allocatable, public, save :: cdb_cBx(:)       ! Collimator beta x
  real(kind=fPrec), allocatable, public, save :: cdb_cBy(:)       ! Collimator beta y
  real(kind=fPrec), allocatable, public, save :: cdb_cTilt(:,:)   ! Collimator tilt
  logical,          allocatable, public, save :: cdb_cFound(:)    ! Found in lattice

  ! Family Arrays
  integer,                       private, save :: cdb_nfam  = 0   ! Number of collimator families
  character(len=:), allocatable, private, save :: cdb_famName(:)  ! Family name
  real(kind=fPrec), allocatable, private, save :: cdb_famNSig(:)  ! Family sigma

contains

subroutine cdb_allocDB

  use mod_alloc
  use numerical_constants

  call alloc(cdb_cName,     mNameLen, cdb_nColl, " ",     "cdb_cName")
  call alloc(cdb_cMaterial, mNameLen, cdb_nColl, " ",     "cdb_cMaterial")
  call alloc(cdb_cFamily,   4,        cdb_nColl, " ",     "cdb_cFamily")
  call alloc(cdb_cNSig,               cdb_nColl, -1,      "cdb_cNSig")
  call alloc(cdb_cLength,             cdb_nColl, zero,    "cdb_cLength")
  call alloc(cdb_cOffset,             cdb_nColl, zero,    "cdb_cOffset")
  call alloc(cdb_cRotation,           cdb_nColl, zero,    "cdb_cRotation")
  call alloc(cdb_cBx,                 cdb_nColl, zero,    "cdb_cBx")
  call alloc(cdb_cBy,                 cdb_nColl, zero,    "cdb_cBy")
  call alloc(cdb_cTilt,               cdb_nColl, 2, zero, "cdb_cTilt")
  call alloc(cdb_cFound,              cdb_nColl, .false., "cdb_cFound")

end subroutine cdb_allocDB

subroutine cdb_allocFam

  use mod_alloc
  use numerical_constants

  call alloc(cdb_famName, 16, cdb_nFam, " ",  "cdb_famName")
  call alloc(cdb_famNSig,     cdb_nFam, zero, "cdb_famNSig")

end subroutine cdb_allocFam

end module collimation_db
