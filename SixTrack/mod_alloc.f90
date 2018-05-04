! ================================================================================================ !
!  SixTrack Array Alloc Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  J. Molson, K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-04-26
! ================================================================================================ !

module mod_alloc
  
  use, intrinsic :: iso_fortran_env
 !use, intrinsic :: ieee_arithmetic
  use crcoall
  
  implicit none
  
  ! The total quantity of memory allocated (bits)
  integer(kind=int64) :: allocated_bits = 0
  
  ! The maximum quantity of memory allocated (bits)
  integer(kind=int64) :: maximum_bits = 0
  
 !integer(kind=int64), parameter :: gbyte = 1024*1024*1024*8
  integer(kind=int64), parameter :: mbyte = 1024*1024*8
  integer(kind=int64), parameter :: kbyte = 1024*8
  integer(kind=int64), parameter :: byte  = 8

! Function/subroutine overloading follows below:

!here we could probably use fPrec and save a lot of code
!some arrays are fixed to a specific variable type, and therefore we are prevented from doing this :(
interface alloc
  
  module procedure resize1di16  ! 1D int16
  
  module procedure resize1di32  ! 1D int32
  module procedure resize2di32  ! 2D int32
  module procedure resize3di32  ! 3D int32
  
  module procedure resize1di64  ! 1D int64
  module procedure resize2di64  ! 2D int64
  module procedure resize3di64  ! 3D int64
  
  module procedure resize1dr32  ! 1D real32
  module procedure resize2dr32  ! 2D real32
  module procedure resize3dr32  ! 3D real32
  module procedure resize4dr32  ! 4D real32
  
  module procedure resize1dr64  ! 1D real64
  module procedure resize2dr64  ! 2D real64
  module procedure resize3dr64  ! 3D real64
  module procedure resize4dr64  ! 4D real64
  
  module procedure resize1dr128 ! 1D real128
  module procedure resize2dr128 ! 2D real128
  module procedure resize3dr128 ! 3D real128
  module procedure resize4dr128 ! 4D real128
  
  module procedure resize1dc    ! 1D character
  module procedure resize2dc    ! 2D character
  
  module procedure resize1ds    ! 1D string
  
  module procedure resize1dl    ! 1D logical
  module procedure resize2dl    ! 2D logical
  
end interface alloc

interface resize
  
  module procedure resize1di16  ! 1D int16
  
  module procedure resize1di32  ! 1D int32
  module procedure resize2di32  ! 2D int32
  module procedure resize3di32  ! 3D int32
  
  module procedure resize1di64  ! 1D int64
  module procedure resize2di64  ! 2D int64
  module procedure resize3di64  ! 3D int64
  
  module procedure resize1dr32  ! 1D real32
  module procedure resize2dr32  ! 2D real32
  module procedure resize3dr32  ! 3D real32
  module procedure resize4dr32  ! 4D real32
  
  module procedure resize1dr64  ! 1D real64
  module procedure resize2dr64  ! 2D real64
  module procedure resize3dr64  ! 3D real64
  module procedure resize4dr64  ! 4D real64
  
  module procedure resize1dr128 ! 1D real128
  module procedure resize2dr128 ! 2D real128
  module procedure resize3dr128 ! 3D real128
  module procedure resize4dr128 ! 4D real128
  
  module procedure resize1dc    ! 1D character
  module procedure resize2dc    ! 2D character
  
  module procedure resize1ds    ! 1D string
  
  module procedure resize1dl    ! 1D logical
  module procedure resize2dl    ! 2D logical
  
end interface resize

interface dealloc
  module procedure dealloc1dr32
  module procedure dealloc1dr64
  module procedure dealloc1dr128

  module procedure dealloc2dr32
  module procedure dealloc2dr64
  module procedure dealloc2dr128

  module procedure dealloc3dr32
  module procedure dealloc3dr64
  module procedure dealloc3dr128

  module procedure dealloc1di16
  module procedure dealloc1di32
  module procedure dealloc1di64

  module procedure dealloc1dl

  module procedure dealloc1dc
end interface dealloc

contains

subroutine alloc_error(ename, error, requested_bits)
  implicit none

  !The variable name we are allocating
  character(len=*), intent(in) :: ename

  !The requested allocation size
  integer(kind=int64), intent(in) :: requested_bits

  !The error code
  integer, intent(in) :: error

  write(lout,*) 'Memory allocation error for array ', ename
  write(lout,*) 'Current allocation is: ', allocated_bits, ' bits'
  write(lout,*) 'Requested allocation is: ', requested_bits, ' bits'
  write(lout,*) 'Allocation error code is ', error
  write(lout,*) 'Exiting!'
  stop
end subroutine alloc_error

subroutine print_alloc(ename, type, requested_bits)
  
  use numerical_constants, only : one
  
  implicit none
  
  character(len=*), intent(in) :: ename, type
  integer(kind=int64) :: requested_bits
  
#ifdef SPAMMY
  
  write(lout,"(a)") "ALLOC> Memory allocation for "//type//" array '"//ename//"'"
  
  if((real(allocated_bits,real64)/mbyte < one) .and. (real(allocated_bits,real64)/kbyte > one)) then
    write(lout,"(a,f7.2,a)") "ALLOC>   Current allocation is:   ",real(allocated_bits,real64)/kbyte," kb"
  else if(real(allocated_bits,real64)/mbyte >= one) then
    write(lout,"(a,f7.2,a)") "ALLOC>   Current allocation is:   ",real(allocated_bits,real64)/mbyte," Mb"
  else
    write(lout,"(a,f7.2,a)") "ALLOC>   Current allocation is:   ",real(allocated_bits,real64)/byte," b"
  end if

  if((real(requested_bits,real64)/mbyte < one) .and. (real(requested_bits,real64)/kbyte > one)) then
    write(lout,"(a,f7.2,a)") "ALLOC>   Requested allocation is: ",real(requested_bits,real64)/kbyte," kb"
  else if(real(requested_bits,real64)/mbyte >= one) then
    write(lout,"(a,f7.2,a)") "ALLOC>   Requested allocation is: ",real(requested_bits,real64)/mbyte," Mb"
  else
    write(lout,"(a,f7.2,a)") "ALLOC>   Requested allocation is: ",real(requested_bits,real64)/byte," b"
  end if
#endif
  
  if(allocated_bits > maximum_bits) then
    maximum_bits = allocated_bits
  end if
  
end subroutine print_alloc

subroutine alloc_exit
  implicit none

  logical lopen
  integer, parameter :: memunit = 47894

  !write the maximum allocated memory to a file
  inquire(unit=memunit, opened=lopen)
  if(.not.lopen) then
    open(memunit, file='maximum_memory_allocation_mbytes.txt', form='formatted' )
    write(memunit,*) real(maximum_bits,real64)/real(mbyte,real64)
    close(memunit)
  else
    write(lout,*) 'WARNING: unit already open for memory usage writing!'
  end if
end subroutine alloc_exit

! ================================================================================================ !
!  Allocation and Resizing
! ~~~~~~~~~~~~~~~~~~~~~~~~~
!  Inputs:
!    ename   :: The name of the variable being assigned,for error printing
!    eIdx    :: The new end index (size) that the array should be set to
!    initial :: The inital value for the array to be set to
!    input   :: The input variable to allocated
!    fIdx    :: The index of the first element
! ================================================================================================ !

! ================================================================================================ !
!  INTEGERS
! ================================================================================================ !

! 1D int16 Array
subroutine resize1di16(input, eIdx, initial, ename, fIdxIn)
  
  implicit none
  
  integer(kind=int16), allocatable, intent(inout) :: input(:)
  integer,                          intent(in)    :: eIdx
  integer(kind=int16),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn
  
  integer(kind=int16), allocatable :: buffer(:)   ! Buffer array
  integer                          :: fIdx        ! First index
  integer                          :: oIdx        ! Old end index
  integer(kind=int64)              :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(int16)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx    = size(input)+fIdx-1
    request = (eIdx-oIdx) * storage_size(int16)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D int16",request)
  
end subroutine resize1di16

! 1D int32 Array
subroutine resize1di32(input, eIdx, initial, ename, fIdxIn)
  
  implicit none
  
  integer(kind=int32), allocatable, intent(inout) :: input(:)
  integer,                          intent(in)    :: eIdx
  integer(kind=int32),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn
  
  integer(kind=int32), allocatable :: buffer(:)   ! Buffer array
  integer                          :: fIdx        ! First index
  integer                          :: oIdx        ! Old end index
  integer(kind=int64)              :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(int32)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx    = size(input)+fIdx-1
    request = (eIdx-oIdx) * storage_size(int32)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D int32",request)
  
end subroutine resize1di32

! 2D int32 Array
subroutine resize2di32(input, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  integer(kind=int32), allocatable, intent(inout) :: input(:,:)
  integer,                          intent(in)    :: eIdx1,eIdx2
  integer(kind=int32),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  integer(kind=int32), allocatable :: buffer(:,:) ! Buffer array
  integer                          :: fIdx1,fIdx2 ! First indices
  integer                          :: oIdx1,oIdx2 ! Old end indices
  integer(kind=int64)              :: request     ! Requested size addition
  
  integer i, j, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) * storage_size(int32)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)) * storage_size(int32)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    
    do i=fIdx2,oIdx2
      do j=fIdx1,oIdx1
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"2D int32",request)
  
end subroutine resize2di32

! 3D int32 Array
subroutine resize3di32(input, eIdx1, eIdx2, eIdx3, initial, ename, fIdxIn1, fIdxIn2, fIdxIn3)
  
  implicit none
  
  integer(kind=int32), allocatable, intent(inout) :: input(:,:,:)
  integer,                          intent(in)    :: eIdx1,eIdx2,eIdx3
  integer(kind=int32),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn1,fIdxIn2,fIdxIn3
  
  integer(kind=int32), allocatable :: buffer(:,:,:)     ! Buffer array
  integer                          :: fIdx1,fIdx2,fIdx3 ! First indices
  integer                          :: oIdx1,oIdx2,oIdx3 ! Old end indices
  integer(kind=int64)              :: request           ! Requested size addition
  
  integer i, j, k, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  if(present(fIdxIn3)) then
    fIdx3 = fIdxIn3
  else
    fIdx3 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1) * storage_size(int32)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          input(k,j,i) = initial
        end do
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    oIdx3   = size(input,3)+fIdx3-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)*(oIdx3-fIdx3+1)) * storage_size(int32)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          buffer(k,j,i) = initial
        end do
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    if(oIdx3 > eIdx3) oIdx3 = eIdx3
    
    do i=fIdx3,oIdx3
      do j=fIdx2,oIdx2
        do k=fIdx1,oIdx1
          buffer(k,j,i) = input(k,j,i)
        end do
      end do
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"3D int32",request)
  
end subroutine resize3di32

! 1D int64 Array
subroutine resize1di64(input, eIdx, initial, ename, fIdxIn)
  
  implicit none
  
  integer(kind=int64), allocatable, intent(inout) :: input(:)
  integer,                          intent(in)    :: eIdx
  integer(kind=int64),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn
  
  integer(kind=int64), allocatable :: buffer(:)   ! Buffer array
  integer                          :: fIdx        ! First index
  integer                          :: oIdx        ! Old end index
  integer(kind=int64)              :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(int64)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx, eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx    = size(input)+fIdx-1
    request = (eIdx-oIdx) * storage_size(int64)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D int64",request)
  
end subroutine resize1di64

! 2D int64 Array
subroutine resize2di64(input, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  integer(kind=int64), allocatable, intent(inout) :: input(:,:)
  integer,                          intent(in)    :: eIdx1,eIdx2
  integer(kind=int64),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  integer(kind=int64), allocatable :: buffer(:,:) ! Buffer array
  integer                          :: fIdx1,fIdx2 ! First indices
  integer                          :: oIdx1,oIdx2 ! Old end indices
  integer(kind=int64)              :: request     ! Requested size addition
  
  integer i, j, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) * storage_size(int64)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)) * storage_size(int64)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    
    do i=fIdx2,oIdx2
      do j=fIdx1,oIdx1
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"2D int64",request)
  
end subroutine resize2di64

! 3D int64 Array
subroutine resize3di64(input, eIdx1, eIdx2, eIdx3, initial, ename, fIdxIn1, fIdxIn2, fIdxIn3)
  
  implicit none
  
  integer(kind=int64), allocatable, intent(inout) :: input(:,:,:)
  integer,                          intent(in)    :: eIdx1,eIdx2,eIdx3
  integer(kind=int64),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn1,fIdxIn2,fIdxIn3
  
  integer(kind=int64), allocatable :: buffer(:,:,:)     ! Buffer array
  integer                          :: fIdx1,fIdx2,fIdx3 ! First indices
  integer                          :: oIdx1,oIdx2,oIdx3 ! Old end indices
  integer(kind=int64)              :: request           ! Requested size addition
  
  integer i, j, k, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  if(present(fIdxIn3)) then
    fIdx3 = fIdxIn3
  else
    fIdx3 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1) * storage_size(int64)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          input(k,j,i) = initial
        end do
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    oIdx3   = size(input,3)+fIdx3-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)*(oIdx3-fIdx3+1)) * storage_size(int64)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          buffer(k,j,i) = initial
        end do
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    if(oIdx3 > eIdx3) oIdx3 = eIdx3
    
    do i=fIdx3,oIdx3
      do j=fIdx2,oIdx2
        do k=fIdx1,oIdx1
          buffer(k,j,i) = input(k,j,i)
        end do
      end do
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"3D int64",request)
  
end subroutine resize3di64

! ================================================================================================ !
!  REALS
! ================================================================================================ !

! 1D real32 Array
subroutine resize1dr32(input, eIdx, initial, ename, fIdxIn)
  
  implicit none
  
  real(kind=real32), allocatable, intent(inout) :: input(:)
  integer,                        intent(in)    :: eIdx
  real(kind=real32),              intent(in)    :: initial
  character(len=*),               intent(in)    :: ename
  integer,           optional,    intent(in)    :: fIdxIn
  
  real(kind=real32), allocatable :: buffer(:)   ! Buffer array
  integer                        :: fIdx        ! First index
  integer                        :: oIdx        ! Old end index
  integer(kind=int64)            :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(real32)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx    = size(input)+fIdx-1
    request = (eIdx-oIdx) * storage_size(real32)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D real32",request)
  
end subroutine resize1dr32

! 2D real32 Array
subroutine resize2dr32(input, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  real(kind=real32), allocatable, intent(inout) :: input(:,:)
  integer,                        intent(in)    :: eIdx1,eIdx2
  real(kind=real32),              intent(in)    :: initial
  character(len=*),               intent(in)    :: ename
  integer,           optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  real(kind=real32), allocatable :: buffer(:,:) ! Buffer array
  integer                        :: fIdx1,fIdx2 ! First indices
  integer                        :: oIdx1,oIdx2 ! Old end indices
  integer(kind=int64)            :: request     ! Requested size addition
  
  integer i, j, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) * storage_size(real32)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)) * storage_size(real32)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    
    do i=fIdx2,oIdx2
      do j=fIdx1,oIdx1
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"2D real32",request)
  
end subroutine resize2dr32

! 3D real32 Array
subroutine resize3dr32(input, eIdx1, eIdx2, eIdx3, initial, ename, fIdxIn1, fIdxIn2, fIdxIn3)
  
  implicit none
  
  real(kind=real32), allocatable, intent(inout) :: input(:,:,:)
  integer,                        intent(in)    :: eIdx1,eIdx2,eIdx3
  real(kind=real32),              intent(in)    :: initial
  character(len=*),               intent(in)    :: ename
  integer,           optional,    intent(in)    :: fIdxIn1,fIdxIn2,fIdxIn3
  
  real(kind=real32), allocatable :: buffer(:,:,:)     ! Buffer array
  integer                        :: fIdx1,fIdx2,fIdx3 ! First indices
  integer                        :: oIdx1,oIdx2,oIdx3 ! Old end indices
  integer(kind=int64)            :: request           ! Requested size addition
  
  integer i, j, k, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  if(present(fIdxIn3)) then
    fIdx3 = fIdxIn3
  else
    fIdx3 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1) * storage_size(real32)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          input(k,j,i) = initial
        end do
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    oIdx3   = size(input,3)+fIdx3-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)*(oIdx3-fIdx3+1)) * storage_size(real32)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          buffer(k,j,i) = initial
        end do
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    if(oIdx3 > eIdx3) oIdx3 = eIdx3
    
    do i=fIdx3,oIdx3
      do j=fIdx2,oIdx2
        do k=fIdx1,oIdx1
          buffer(k,j,i) = input(k,j,i)
        end do
      end do
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"3D real32",request)
  
end subroutine resize3dr32

! 4D real32 Array
subroutine resize4dr32(input, eIdx1, eIdx2, eIdx3, eIdx4, initial, ename, fIdxIn1, fIdxIn2, fIdxIn3, fIdxIn4)
  
  implicit none
  
  real(kind=real32), allocatable, intent(inout) :: input(:,:,:,:)
  integer,                        intent(in)    :: eIdx1,eIdx2,eIdx3,eIdx4
  real(kind=real32),              intent(in)    :: initial
  character(len=*),               intent(in)    :: ename
  integer,           optional,    intent(in)    :: fIdxIn1,fIdxIn2,fIdxIn3,fIdxIn4
  
  real(kind=real32), allocatable :: buffer(:,:,:,:)         ! Buffer array
  integer                        :: fIdx1,fIdx2,fIdx3,fIdx4 ! First indices
  integer                        :: oIdx1,oIdx2,oIdx3,oIdx4 ! Old end indices
  integer(kind=int64)            :: request                 ! Requested size addition
  
  integer i, j, k, l, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  if(present(fIdxIn3)) then
    fIdx3 = fIdxIn3
  else
    fIdx3 = 1
  end if
  if(present(fIdxIn4)) then
    fIdx4 = fIdxIn4
  else
    fIdx4 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)*(eIdx4-fIdx4+1) * storage_size(real32)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3,fIdx4:eIdx4), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx4,eIdx4
      do j=fIdx3,eIdx3
        do k=fIdx2,eIdx2
          do l=fIdx1,eIdx1
            input(l,k,j,i) = initial
          end do
        end do
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    oIdx3   = size(input,3)+fIdx3-1
    oIdx4   = size(input,4)+fIdx4-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)*(eIdx4-fIdx4+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)*(oIdx3-fIdx3+1)*(oIdx4-fIdx4+1)) * storage_size(real32)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3,fIdx4:eIdx4), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx4,eIdx4
      do j=fIdx3,eIdx3
        do k=fIdx2,eIdx2
          do l=fIdx1,eIdx1
            buffer(l,k,j,i) = initial
          end do
        end do
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    if(oIdx3 > eIdx3) oIdx3 = eIdx3
    if(oIdx4 > eIdx4) oIdx4 = eIdx4
    
    do i=fIdx4,oIdx4
      do j=fIdx3,oIdx3
        do k=fIdx2,oIdx2
          do l=fIdx1,oIdx1
            buffer(l,k,j,i) = input(l,k,j,i)
          end do
        end do
      end do
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"4D real32",request)
  
end subroutine resize4dr32

! 1D real64 Array
subroutine resize1dr64(input, eIdx, initial, ename, fIdxIn)
  
  implicit none
  
  real(kind=real64), allocatable, intent(inout) :: input(:)
  integer,                        intent(in)    :: eIdx
  real(kind=real64),              intent(in)    :: initial
  character(len=*),               intent(in)    :: ename
  integer,           optional,    intent(in)    :: fIdxIn
  
  real(kind=real64), allocatable :: buffer(:)   ! Buffer array
  integer                        :: fIdx        ! First index
  integer                        :: oIdx        ! Old end index
  integer(kind=int64)            :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(real64)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx   = size(input)+fIdx-1
    request = (eIdx-oIdx) * storage_size(real64)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D real64",request)
  
end subroutine resize1dr64

! 2D real64 Array
subroutine resize2dr64(input, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  real(kind=real64), allocatable, intent(inout) :: input(:,:)
  integer,                        intent(in)    :: eIdx1,eIdx2
  real(kind=real64),              intent(in)    :: initial
  character(len=*),               intent(in)    :: ename
  integer,           optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  real(kind=real64), allocatable :: buffer(:,:) ! Buffer array
  integer                        :: fIdx1,fIdx2 ! First indices
  integer                        :: oIdx1,oIdx2 ! Old end indices
  integer(kind=int64)            :: request     ! Requested size addition
  
  integer i, j, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) * storage_size(real64)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)) * storage_size(real64)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    
    do i=fIdx2,oIdx2
      do j=fIdx1,oIdx1
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"2D real64",request)
  
end subroutine resize2dr64

! 3D real64 Array
subroutine resize3dr64(input, eIdx1, eIdx2, eIdx3, initial, ename, fIdxIn1, fIdxIn2, fIdxIn3)
  
  implicit none
  
  real(kind=real64), allocatable, intent(inout) :: input(:,:,:)
  integer,                        intent(in)    :: eIdx1,eIdx2,eIdx3
  real(kind=real64),              intent(in)    :: initial
  character(len=*),               intent(in)    :: ename
  integer,           optional,    intent(in)    :: fIdxIn1,fIdxIn2,fIdxIn3
  
  real(kind=real64), allocatable :: buffer(:,:,:)     ! Buffer array
  integer                        :: fIdx1,fIdx2,fIdx3 ! First indices
  integer                        :: oIdx1,oIdx2,oIdx3 ! Old end indices
  integer(kind=int64)            :: request           ! Requested size addition
  
  integer i, j, k, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  if(present(fIdxIn3)) then
    fIdx3 = fIdxIn3
  else
    fIdx3 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1) * storage_size(real64)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          input(k,j,i) = initial
        end do
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    oIdx3   = size(input,3)+fIdx3-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)*(oIdx3-fIdx3+1)) * storage_size(real64)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          buffer(k,j,i) = initial
        end do
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    if(oIdx3 > eIdx3) oIdx3 = eIdx3
    
    do i=fIdx3,oIdx3
      do j=fIdx2,oIdx2
        do k=fIdx1,oIdx1
          buffer(k,j,i) = input(k,j,i)
        end do
      end do
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"3D real64",request)
  
end subroutine resize3dr64

! 4d real64 array
subroutine resize4dr64(input, eIdx1, eIdx2, eIdx3, eIdx4, initial, ename, fIdxIn1, fIdxIn2, fIdxIn3, fIdxIn4)
  
  implicit none
  
  real(kind=real64), allocatable, intent(inout) :: input(:,:,:,:)
  integer,                        intent(in)    :: eIdx1,eIdx2,eIdx3,eIdx4
  real(kind=real64),              intent(in)    :: initial
  character(len=*),               intent(in)    :: ename
  integer,           optional,    intent(in)    :: fIdxIn1,fIdxIn2,fIdxIn3,fIdxIn4
  
  real(kind=real64), allocatable :: buffer(:,:,:,:)         ! Buffer array
  integer                        :: fIdx1,fIdx2,fIdx3,fIdx4 ! First indices
  integer                        :: oIdx1,oIdx2,oIdx3,oIdx4 ! Old end indices
  integer(kind=int64)            :: request                 ! Requested size addition
  
  integer i, j, k, l, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  if(present(fIdxIn3)) then
    fIdx3 = fIdxIn3
  else
    fIdx3 = 1
  end if
  if(present(fIdxIn4)) then
    fIdx4 = fIdxIn4
  else
    fIdx4 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)*(eIdx4-fIdx4+1) * storage_size(real64)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3,fIdx4:eIdx4), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx4,eIdx4
      do j=fIdx3,eIdx3
        do k=fIdx2,eIdx2
          do l=fIdx1,eIdx1
            input(l,k,j,i) = initial
          end do
        end do
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    oIdx3   = size(input,3)+fIdx3-1
    oIdx4   = size(input,4)+fIdx4-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)*(eIdx4-fIdx4+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)*(oIdx3-fIdx3+1)*(oIdx4-fIdx4+1)) * storage_size(real64)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3,fIdx4:eIdx4), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx4,eIdx4
      do j=fIdx3,eIdx3
        do k=fIdx2,eIdx2
          do l=fIdx1,eIdx1
            buffer(l,k,j,i) = initial
          end do
        end do
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    if(oIdx3 > eIdx3) oIdx3 = eIdx3
    if(oIdx4 > eIdx4) oIdx4 = eIdx4
    
    do i=fIdx4,oIdx4
      do j=fIdx3,oIdx3
        do k=fIdx2,oIdx2
          do l=fIdx1,oIdx1
            buffer(l,k,j,i) = input(l,k,j,i)
          end do
        end do
      end do
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"4D real64",request)
  
end subroutine resize4dr64

! 1D real128 Array
subroutine resize1dr128(input, eIdx, initial, ename, fIdxIn)
  
  implicit none
  
  real(kind=real128), allocatable, intent(inout) :: input(:)
  integer,                         intent(in)    :: eIdx
  real(kind=real128),              intent(in)    :: initial
  character(len=*),                intent(in)    :: ename
  integer,            optional,    intent(in)    :: fIdxIn
  
  real(kind=real128), allocatable :: buffer(:)   ! Buffer array
  integer                         :: fIdx        ! First index
  integer                         :: oIdx        ! Old end index
  integer(kind=int64)             :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(real128)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx   = size(input)+fIdx-1
    request = (eIdx-oIdx) * storage_size(real128)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D real128",request)
  
end subroutine resize1dr128

! 2D real128 Array
subroutine resize2dr128(input, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  real(kind=real128), allocatable, intent(inout) :: input(:,:)
  integer,                         intent(in)    :: eIdx1,eIdx2
  real(kind=real128),              intent(in)    :: initial
  character(len=*),                intent(in)    :: ename
  integer,            optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  real(kind=real128), allocatable :: buffer(:,:) ! Buffer array
  integer                         :: fIdx1,fIdx2 ! First indices
  integer                         :: oIdx1,oIdx2 ! Old end indices
  integer(kind=int64)             :: request     ! Requested size addition
  
  integer i, j, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) * storage_size(real128)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    request = ((eIdx1-fIdx1+1) * (eIdx2-fIdx2+1) - &
               (eIdx1-oIdx1)  * (eIdx2-oIdx2)) * storage_size(real128)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    
    do i=fIdx2,oIdx2
      do j=fIdx1,oIdx1
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"2D real128",request)
  
end subroutine resize2dr128

! 3D real128 Array
subroutine resize3dr128(input, eIdx1, eIdx2, eIdx3, initial, ename, fIdxIn1, fIdxIn2, fIdxIn3)
  
  implicit none
  
  real(kind=real128), allocatable, intent(inout) :: input(:,:,:)
  integer,                         intent(in)    :: eIdx1,eIdx2,eIdx3
  real(kind=real128),              intent(in)    :: initial
  character(len=*),                intent(in)    :: ename
  integer,            optional,    intent(in)    :: fIdxIn1,fIdxIn2,fIdxIn3
  
  real(kind=real128), allocatable :: buffer(:,:,:)     ! Buffer array
  integer                         :: fIdx1,fIdx2,fIdx3 ! First indices
  integer                         :: oIdx1,oIdx2,oIdx3 ! Old end indices
  integer(kind=int64)             :: request           ! Requested size addition
  
  integer i, j, k, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  if(present(fIdxIn3)) then
    fIdx3 = fIdxIn3
  else
    fIdx3 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1) * storage_size(real128)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          input(k,j,i) = initial
        end do
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    oIdx3   = size(input,3)+fIdx3-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)*(oIdx3-fIdx3+1)) * storage_size(real128)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx3,eIdx3
      do j=fIdx2,eIdx2
        do k=fIdx1,eIdx1
          buffer(k,j,i) = initial
        end do
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    if(oIdx3 > eIdx3) oIdx3 = eIdx3
    
    do i=fIdx3,oIdx3
      do j=fIdx2,oIdx2
        do k=fIdx1,oIdx1
          buffer(k,j,i) = input(k,j,i)
        end do
      end do
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"3D real128",request)
  
end subroutine resize3dr128

! 4D real128 Array
subroutine resize4dr128(input, eIdx1, eIdx2, eIdx3, eIdx4, initial, ename, fIdxIn1, fIdxIn2, fIdxIn3, fIdxIn4)
  
  implicit none
  
  real(kind=real128), allocatable, intent(inout) :: input(:,:,:,:)
  integer,                         intent(in)    :: eIdx1,eIdx2,eIdx3,eIdx4
  real(kind=real128),              intent(in)    :: initial
  character(len=*),                intent(in)    :: ename
  integer,            optional,    intent(in)    :: fIdxIn1,fIdxIn2,fIdxIn3,fIdxIn4
  
  real(kind=real128), allocatable :: buffer(:,:,:,:)         ! Buffer array
  integer                         :: fIdx1,fIdx2,fIdx3,fIdx4 ! First indices
  integer                         :: oIdx1,oIdx2,oIdx3,oIdx4 ! Old end indices
  integer(kind=int64)             :: request                 ! Requested size addition
  
  integer i, j, k, l, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  if(present(fIdxIn3)) then
    fIdx3 = fIdxIn3
  else
    fIdx3 = 1
  end if
  if(present(fIdxIn4)) then
    fIdx4 = fIdxIn4
  else
    fIdx4 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)*(eIdx4-fIdx4+1) * storage_size(real128)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3,fIdx4:eIdx4), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx4,eIdx4
      do j=fIdx3,eIdx3
        do k=fIdx2,eIdx2
          do l=fIdx1,eIdx1
            input(l,k,j,i) = initial
          end do
        end do
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    oIdx3   = size(input,3)+fIdx3-1
    oIdx4   = size(input,4)+fIdx4-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)*(eIdx3-fIdx3+1)*(eIdx4-fIdx4+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)*(oIdx3-fIdx3+1)*(oIdx4-fIdx4+1)) * storage_size(real128)
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2,fIdx3:eIdx3,fIdx4:eIdx4), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx4,eIdx4
      do j=fIdx3,eIdx3
        do k=fIdx2,eIdx2
          do l=fIdx1,eIdx1
            buffer(l,k,j,i) = initial
          end do
        end do
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    if(oIdx3 > eIdx3) oIdx3 = eIdx3
    if(oIdx4 > eIdx4) oIdx4 = eIdx4
    
    do i=fIdx4,oIdx4
      do j=fIdx3,oIdx3
        do k=fIdx2,oIdx2
          do l=fIdx1,oIdx1
            buffer(l,k,j,i) = input(l,k,j,i)
          end do
        end do
      end do
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"4D real128",request)
  
end subroutine resize4dr128

! ================================================================================================ !
!  CHARACTERS
! ================================================================================================ !

! 1D character Array
subroutine resize1dc(input, strlen, eIdx, initial, ename, fIdxIn)
  
  implicit none
  
  character(len=:), allocatable, intent(inout) :: input(:)
  integer,                       intent(in)    :: strlen
  integer,                       intent(in)    :: eIdx
  character(len=*),              intent(in)    :: initial
  character(len=*),              intent(in)    :: ename
  integer,          optional,    intent(in)    :: fIdxIn
  
  character(len=:), allocatable :: buffer(:)   ! Buffer array
  integer                       :: fIdx        ! First index
  integer                       :: oIdx        ! Old end index
  integer(kind=int64)           :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * strlen * CHARACTER_STORAGE_SIZE
    
    allocate(character(strlen) :: input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx    = size(input)+fIdx-1
    request = (eIdx-oIdx) * strlen * CHARACTER_STORAGE_SIZE
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(character(strlen) :: buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D character",request)
  
end subroutine resize1dc

! 2D character Array
subroutine resize2dc(input, strlen, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  character(len=:), allocatable, intent(inout) :: input(:,:)
  integer,                       intent(in)    :: strlen
  integer,                       intent(in)    :: eIdx1,eIdx2
  character(len=*),              intent(in)    :: initial
  character(len=*),              intent(in)    :: ename
  integer,          optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  character(len=:), allocatable :: buffer(:,:) ! Buffer array
  integer                       :: fIdx1,fIdx2 ! First index
  integer                       :: oIdx1,oIdx2 ! Old end index
  integer(kind=int64)           :: request     ! Requested size addition
  
  integer i, j, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) * strlen * CHARACTER_STORAGE_SIZE
    
    allocate(character(strlen) :: input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)  - &
               (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1)) * strlen * CHARACTER_STORAGE_SIZE
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(character(strlen) :: buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    
    do i=fIdx2,oIdx2
      do j=fIdx1,oIdx1
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"2D character",request)
  
end subroutine resize2dc

! ==================================================================== !
!  STRINGS
! ==================================================================== !

subroutine resize1ds(input, eIdx, initial, ename, fIdxIn)
  
  use strings
  
  implicit none
  
  type(string), allocatable, intent(inout) :: input(:)
  integer,                   intent(in)    :: eIdx
  type(string),              intent(in)    :: initial
  character(len=*),          intent(in)    :: ename
  integer,      optional,    intent(in)    :: fIdxIn
  
  type(string), allocatable :: buffer(:)   ! Buffer array
  integer                   :: fIdx        ! First index
  integer                   :: oIdx        ! Old end index
  integer(kind=int64)       :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = eIdx-fIdx+1
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx    = size(input)+fIdx-1
    request = eIdx-oIdx
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D string",request)
  
end subroutine resize1ds

! ================================================================================================ !
!  LOGICALS
! ================================================================================================ !

! 1D logical Array
subroutine resize1dl(input, eIdx, initial, ename, fIdxIn)
  
  implicit none
  
  logical, allocatable, intent(inout) :: input(:)
  integer,              intent(in)    :: eIdx
  logical,              intent(in)    :: initial
  character(len=*),     intent(in)    :: ename
  integer, optional,    intent(in)    :: fIdxIn
  
  logical, allocatable :: buffer(:)   ! Buffer array
  integer              :: fIdx        ! First index
  integer              :: oIdx        ! Old end index
  integer(kind=int64)  :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = eIdx-fIdx+1
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oIdx    = size(input)+fIdx-1
    request = eIdx-oIdx
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oIdx > eIdx) oIdx = eIdx
    do i=fIdx,oIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"1D logical",request)
  
end subroutine resize1dl

! 2D logical Array
subroutine resize2dl(input, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  logical, allocatable, intent(inout) :: input(:,:)
  integer,              intent(in)    :: eIdx1,eIdx2
  logical,              intent(in)    :: initial
  character(len=*),     intent(in)    :: ename
  integer, optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  logical, allocatable :: buffer(:,:) ! Buffer array
  integer              :: fIdx1,fIdx2 ! First indices
  integer              :: oIdx1,oIdx2 ! Old end indices
  integer(kind=int64)  :: request     ! Requested size addition
  
  integer i, j, error
  
  if(present(fIdxIn1)) then
    fIdx1 = fIdxIn1
  else
    fIdx1 = 1
  end if
  if(present(fIdxIn2)) then
    fIdx2 = fIdxIn2
  else
    fIdx2 = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oIdx1   = size(input,1)+fIdx1-1
    oIdx2   = size(input,2)+fIdx2-1
    request = ((eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) - (oIdx1-fIdx1+1)*(oIdx2-fIdx2+1))
    
    if(request == 0.0) then
      write(lout,"(a)") "ALLOC> No additional allocating needed for array '"//ename//"'"
      return
    end if
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oIdx1 > eIdx1) oIdx1 = eIdx1
    if(oIdx2 > eIdx2) oIdx2 = eIdx2
    
    do i=fIdx2,oIdx2
      do j=fIdx1,oIdx1
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,"2D logical",request)
  
end subroutine resize2dl

! ================================================================================================ !
!  DEALLOCTIONS
! ~~~~~~~~~~~~~~
!  These exist to keep track of the total allocated memory
! ================================================================================================ !

subroutine dealloc1dr32(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real32), allocatable, intent(inout) :: input(:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input)*storage_size(real32))
  deallocate(input)
end subroutine dealloc1dr32

subroutine dealloc1dr64(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real64), allocatable, intent(inout) :: input(:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input)*storage_size(real64))
  deallocate(input)
end subroutine dealloc1dr64

subroutine dealloc1dr128(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real128), allocatable, intent(inout) :: input(:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input)*storage_size(real128))
  deallocate(input)
end subroutine dealloc1dr128

subroutine dealloc2dr32(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real32), allocatable, intent(inout) :: input(:,:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input,1)*size(input,2)*storage_size(real32))
  deallocate(input)
end subroutine dealloc2dr32

subroutine dealloc2dr64(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real64), allocatable, intent(inout) :: input(:,:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input,1)*size(input,2)*storage_size(real64))
  deallocate(input)
end subroutine dealloc2dr64

subroutine dealloc2dr128(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real128), allocatable, intent(inout) :: input(:,:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input,1)*size(input,2)*storage_size(real128))
  deallocate(input)
end subroutine dealloc2dr128

subroutine dealloc3dr32(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real32), allocatable, intent(inout) :: input(:,:,:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input,1)*size(input,2)*size(input,3)*storage_size(real32))
  deallocate(input)
end subroutine dealloc3dr32

subroutine dealloc3dr64(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real64), allocatable, intent(inout) :: input(:,:,:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input,1)*size(input,2)*size(input,3)*storage_size(real64))
  deallocate(input)
end subroutine dealloc3dr64

subroutine dealloc3dr128(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  real(kind=real128), allocatable, intent(inout) :: input(:,:,:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input,1)*size(input,2)*size(input,3)*storage_size(real128))
  deallocate(input)
end subroutine dealloc3dr128

subroutine dealloc1di16(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  integer(kind=int16), allocatable, intent(inout) :: input(:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input)*storage_size(int16))
  deallocate(input)
end subroutine dealloc1di16

subroutine dealloc1di32(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  integer(kind=int32), allocatable, intent(inout) :: input(:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input)*storage_size(int32))
  deallocate(input)
end subroutine dealloc1di32

subroutine dealloc1di64(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  integer(kind=int64), allocatable, intent(inout) :: input(:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - (size(input)*storage_size(int64))
  deallocate(input)
end subroutine dealloc1di64

subroutine dealloc1dl(input, ename)
  implicit none
  character(len=*), intent(in) :: ename
  logical, allocatable, intent(inout) :: input(:)
  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  allocated_bits = allocated_bits - size(input)
  deallocate(input)
end subroutine dealloc1dl

subroutine dealloc1dc(input, strlen, ename)
  implicit none
  character(len=*), intent(in) :: ename
  character(len=:), allocatable, intent(inout) :: input(:)
  integer, intent(in) :: strlen

  !Check that we are already allocated
  if(allocated(input) .eqv. .FALSE.) then
    write(lout,*) 'ERROR: Trying to deallocate a NULL pointer: ', ename
    stop
  end if
  !TODO: Check what size(input) actually does
  allocated_bits = allocated_bits - (size(input)*strlen * CHARACTER_STORAGE_SIZE)
  deallocate(input)
end subroutine dealloc1dc

end module mod_alloc
