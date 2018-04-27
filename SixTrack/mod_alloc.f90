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
  
  module procedure resize1di16 ! 1D int16
  module procedure resize1di32 ! 1D int32
  module procedure resize2di32 ! 2D int32
  module procedure alloc3di32
  module procedure resize1di64 ! 1D int64
  module procedure resize2di64 ! 2D int64
  module procedure alloc3di64

!real
!single
  module procedure alloc1dr32
  module procedure alloc2dr32
  module procedure alloc3dr32
  module procedure alloc4dr32

!double
  module procedure alloc1dr64
  module procedure alloc2dr64
  module procedure alloc3dr64
  module procedure alloc4dr64

!quad
  module procedure alloc1dr128
  module procedure alloc2dr128
  module procedure alloc3dr128
  module procedure alloc4dr128

!characters
  module procedure alloc1dc
  module procedure alloc2dc

  ! Strings
  module procedure alloc1ds

!logical
  module procedure alloc1dl
  module procedure alloc2dl
end interface alloc

interface resize

  module procedure resize1di16

  module procedure resize1di32
  module procedure resize2di32
  module procedure resize3di32

  module procedure resize1di64
  module procedure resize2di64
  module procedure resize3di64

  module procedure resize1dr32
  module procedure resize2dr32
  module procedure resize3dr32
  module procedure resize4dr32

  module procedure resize1dr64
  module procedure resize2dr64
  module procedure resize3dr64
  module procedure resize4dr64

  module procedure resize1dr128
  module procedure resize2dr128
  module procedure resize3dr128
  module procedure resize4dr128

  module procedure resize1dc
  module procedure resize2dc
  
  ! Strings
  module procedure resize1ds

  module procedure resize1dl
  module procedure resize2dl
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

subroutine print_alloc(ename, requested_bits)

  use numerical_constants, only : one
  use file_units

  implicit none

  character(len=*), intent(in) :: ename
  integer(kind=int64) :: requested_bits
  
#ifdef spammy
  
  write(lout,*) 'Memory allocation for array ', ename

  if((real(allocated_bits,real64)/mbyte .lt. one) .and. (real(allocated_bits,real64)/kbyte .gt. one)) then
    write(lout,*) 'Current allocation is: ', real(allocated_bits,real64)/kbyte, ' kbytes'
  else if(real(allocated_bits,real64)/mbyte .ge. one) then
    write(lout,*) 'Current allocation is: ', real(allocated_bits,real64)/mbyte, ' Mbytes'
  else
    write(lout,*) 'Current allocation is: ', real(allocated_bits,real64)/byte, ' bytes'
  end if

  if((real(requested_bits,real64)/mbyte .lt. one) .and. (real(requested_bits,real64)/kbyte .gt. one)) then
    write(lout,*) 'Requested allocation is: ', real(requested_bits,real64)/kbyte, ' kbytes'
  else if(real(requested_bits,real64)/mbyte .ge. one) then
    write(lout,*) 'Requested allocation is: ', real(requested_bits,real64)/mbyte, ' Mbytes'
  else
    write(lout,*) 'Requested allocation is: ', real(requested_bits,real64)/byte, ' bytes'
  end if
  write(lout,*) ''

  if(allocated_bits.gt.maximum_bits) then
    maximum_bits = allocated_bits
  end if

#endif

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

!alloc int32
subroutine alloc3di32(input, startsize1, startsize2, startsize3, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  integer(kind=int32), allocatable, intent(inout) :: input(:,:,:)

  !The inital value for the array to be set to
  integer(kind=int32), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2, startsize3

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l,m

  request = startsize1 * startsize2 * startsize3 * storage_size(int32)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2, startsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      do m=1, startsize3
       input(k,l,m) = initial
      end do
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc3di32

!alloc int64
subroutine alloc3di64(input, startsize1, startsize2, startsize3, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  integer(kind=int64), allocatable, intent(inout) :: input(:,:,:)

  !The inital value for the array to be set to
  integer(kind=int64), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2, startsize3

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l,m

  request = startsize1 * startsize2 * startsize3 * storage_size(int64)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2, startsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      do m=1, startsize3
        input(k,l,m) = initial
      end do
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc3di64


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
  integer                          :: oeIdx       ! Old end index
  integer(kind=int64)              :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> INFO Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(int16)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oeIdx = size(input)+fIdx-1
    request = (eIdx-oeIdx) * storage_size(int16)
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oeIdx > eIdx) oeIdx = eIdx
    do i=fIdx,oeIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)
  
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
  integer                          :: oeIdx       ! Old end index
  integer(kind=int64)              :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> INFO Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(int32)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      input(i) = initial
    end do
    
  else
    
    oeIdx = size(input)+fIdx-1
    request = (eIdx-oeIdx) * storage_size(int32)
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx, eIdx
      buffer(i) = initial
    end do
    
    if(oeIdx > eIdx) oeIdx = eIdx
    do i=fIdx,oeIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)
  
end subroutine resize1di32

! 2D int32 Array
subroutine resize2di32(input, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  integer(kind=int32), allocatable, intent(inout) :: input(:,:)
  integer,                          intent(in)    :: eIdx1,eIdx2
  integer(kind=int32),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  integer(kind=int32), allocatable :: buffer(:,:)   ! Buffer array
  integer                          :: fIdx1,fIdx2   ! First indices
  integer                          :: oeIdx1,oeIdx2 ! Old end indices
  integer(kind=int64)              :: request       ! Requested size addition
  
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
    
    write(lout,"(a)") "ALLOC> INFO Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) * storage_size(int32)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oeIdx1 = size(input,1)+fIdx1-1
    oeIdx2 = size(input,2)+fIdx2-1
    request = (eIdx1-oeIdx1)*(eIdx2-oeIdx2) * storage_size(int32)
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oeIdx1 > eIdx1) oeIdx1 = eIdx1
    if(oeIdx2 > eIdx2) oeIdx2 = eIdx2
    
    do i=fIdx2,oeIdx2
      do j=fIdx1,oeIdx2
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)
  
end subroutine resize2di32


!resize 3d int32 array
subroutine resize3di32(input, newsize1, newsize2, newsize3, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  integer(kind=int32), allocatable, intent(inout) :: input(:,:,:)

  !A buffer used to expand the array
  integer(kind=int32), allocatable :: buffer(:,:,:)

  !The inital value for the array to be set to
  integer(kind=int32), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2, newsize3

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2, oldsize3

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j, k

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, newsize3, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)
  oldsize3 = size(input,3)

  !log our request in size change
  request = ((newsize1*newsize2*newsize3) - (oldsize1*oldsize2*oldsize3)) * storage_size(int32)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2, newsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize3
    do j=1,newsize2
      do k=1,newsize1
        buffer(k,j,i) = initial
      end do
    end do
  end do
  !Copy the data over
  do i=1,oldsize3
    do j=1,oldsize2
      do k=1,oldsize1
        buffer(k,j,i)=input(k,j,i)
      end do
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

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
  integer                          :: oeIdx       ! Old end index
  integer(kind=int64)              :: request     ! Requested size addition
  
  integer i, error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .true.) then
    
    write(lout,"(a)") "ALLOC> INFO Allocating array '"//ename//"'"
    
    request = (eIdx-fIdx+1) * storage_size(int64)
    
    allocate(input(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx, eIdx
      input(i) = initial
    end do
    
  else
    
    oeIdx = size(input)+fIdx-1
    request = (eIdx-oeIdx) * storage_size(int64)
    
    allocate(buffer(fIdx:eIdx), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx,eIdx
      buffer(i) = initial
    end do
    
    if(oeIdx > eIdx) oeIdx = eIdx
    do i=fIdx,oeIdx
      buffer(i) = input(i)
    end do
    
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)
  
end subroutine resize1di64

! 2D int64 Array
subroutine resize2di64(input, eIdx1, eIdx2, initial, ename, fIdxIn1, fIdxIn2)
  
  implicit none
  
  integer(kind=int64), allocatable, intent(inout) :: input(:,:)
  integer,                          intent(in)    :: eIdx1,eIdx2
  integer(kind=int64),              intent(in)    :: initial
  character(len=*),                 intent(in)    :: ename
  integer,             optional,    intent(in)    :: fIdxIn1,fIdxIn2
  
  integer(kind=int64), allocatable :: buffer(:,:)   ! Buffer array
  integer                          :: fIdx1,fIdx2   ! First indices
  integer                          :: oeIdx1,oeIdx2 ! Old end indices
  integer(kind=int64)              :: request       ! Requested size addition
  
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
    
    write(lout,"(a)") "ALLOC> INFO Allocating array '"//ename//"'"
    
    request = (eIdx1-fIdx1+1)*(eIdx2-fIdx2+1) * storage_size(int64)
    
    allocate(input(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        input(j,i) = initial
      end do
    end do
    
  else
    
    oeIdx1 = size(input,1)+fIdx1-1
    oeIdx2 = size(input,2)+fIdx2-1
    request = (eIdx1-oeIdx1)*(eIdx2-oeIdx2) * storage_size(int64)
    
    allocate(buffer(fIdx1:eIdx1,fIdx2:eIdx2), stat=error)
    if(error /= 0) call alloc_error(ename, error, request)
    
    do i=fIdx2,eIdx2
      do j=fIdx1,eIdx1
        buffer(j,i) = initial
      end do
    end do
    
    if(oeIdx1 > eIdx1) oeIdx1 = eIdx1
    if(oeIdx2 > eIdx2) oeIdx2 = eIdx2
    
    do i=fIdx2,oeIdx2
      do j=fIdx1,oeIdx2
        buffer(j,i) = input(j,i)
      end do
    end do
        
    call move_alloc(buffer,input)
    
  end if
  
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)
  
end subroutine resize2di64

!resize 3d int64 array
subroutine resize3di64(input, newsize1, newsize2, newsize3, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  integer(kind=int64), allocatable, intent(inout) :: input(:,:,:)

  !A buffer used to expand the array
  integer(kind=int64), allocatable :: buffer(:,:,:)

  !The inital value for the array to be set to
  integer(kind=int64), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2, newsize3

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2, oldsize3

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j, k

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, newsize3, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)
  oldsize3 = size(input,3)

  !log our request in size change
  request = ((newsize1*newsize2*newsize3) - (oldsize1*oldsize2*oldsize3)) * storage_size(int64)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2, newsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize3
    do j=1,newsize2
      do k=1,newsize1
        buffer(k,j,i) = initial
      end do
    end do
  end do

  !Copy the data over
  do i=1,oldsize3
    do j=1,oldsize2
      do k=1,oldsize1
        buffer(k,j,i)=input(k,j,i)
      end do
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize3di64

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            REALS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine alloc1dr32(input, startsize, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real32), allocatable, intent(inout) :: input(:)

  !The inital value for the array to be set to
  real(kind=real32), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k

  request = startsize * storage_size(real32)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize
    input(k) = initial
  end do

  call print_alloc(ename,request)
end subroutine alloc1dr32


subroutine alloc2dr32(input, startsize1, startsize2, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real32), allocatable, intent(inout) :: input(:,:)

  !The inital value for the array to be set to
  real(kind=real32), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l

  request = startsize1 * startsize2 * storage_size(real32)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
     input(k,l) = initial
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc2dr32

subroutine alloc3dr32(input, startsize1, startsize2, startsize3, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real32), allocatable, intent(inout) :: input(:,:,:)

  !The inital value for the array to be set to
  real(kind=real32), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2, startsize3

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l,m

  request = startsize1 * startsize2 * startsize3 * storage_size(real32)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2, startsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      do m=1, startsize3
       input(k,l,m) = initial
      end do
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc3dr32

subroutine alloc4dr32(input, startsize1, startsize2, startsize3, startsize4, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real32), allocatable, intent(inout) :: input(:,:,:,:)

  !The inital value for the array to be set to
  real(kind=real32), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2, startsize3, startsize4

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l,m,n

  request = startsize1 * startsize2 * startsize3 * startsize4 * storage_size(real32)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2, startsize3, startsize4), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      do m=1, startsize3
        do n=1, startsize4
         input(k,l,m,n) = initial
        end do
      end do
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc4dr32

!alloc real64
subroutine alloc1dr64(input, e_index, initial, ename, fIdxIn)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real64), allocatable, intent(inout) :: input(:)

  !The inital value for the array to be set to
  real(kind=real64), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: e_index

  !The initial index
  integer, intent(in), optional :: fIdxIn
  integer :: fIdx

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k

  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if

  request = (e_index-fIdx+1) * storage_size(real64)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(fIdx:e_index), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=fIdx, e_index
    input(k) = initial
  end do

  call print_alloc(ename,request)
end subroutine alloc1dr64


subroutine alloc2dr64(input, startsize1, startsize2, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real64), allocatable, intent(inout) :: input(:,:)

  !The inital value for the array to be set to
  real(kind=real64), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l

  request = startsize1 * startsize2 * storage_size(real64)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
     input(k,l) = initial
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc2dr64

subroutine alloc3dr64(input, startsize1, startsize2, startsize3, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real64), allocatable, intent(inout) :: input(:,:,:)

  !The inital value for the array to be set to
  real(kind=real64), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2, startsize3

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l,m

  request = startsize1 * startsize2 * startsize3 * storage_size(real64)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2, startsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      do m=1, startsize3
       input(k,l,m) = initial
      end do
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc3dr64

subroutine alloc4dr64(input, startsize1, startsize2, startsize3, startsize4, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real64), allocatable, intent(inout) :: input(:,:,:,:)

  !The inital value for the array to be set to
  real(kind=real64), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2, startsize3, startsize4

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l,m,n

  request = startsize1 * startsize2 * startsize3 * startsize4 * storage_size(real64)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2, startsize3, startsize4), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      do m=1, startsize3
        do n=1, startsize4
         input(k,l,m,n) = initial
        end do
      end do
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc4dr64

!resize 1d real32 array
subroutine resize1dr32(input, newsize, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real32), allocatable, intent(inout) :: input(:)

  !A buffer used to expand the array
  real(kind=real32), allocatable :: buffer(:)

  !The inital value for the array to be set to
  real(kind=real32), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize

  !A variable to track the old size of the array
  integer :: oldsize

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize = size(input)

  !log our request in size change
  request = (newsize-oldsize) * storage_size(real32)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Copy the data over
  do i=1,oldsize
    buffer(i) = input(i)
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Set the initial values of the buffer
  if(newsize.gt.oldsize) then
    do i=oldsize+1,newsize
      buffer(i) = initial
    end do
  end if

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize1dr32

!resize 2d real32 array
subroutine resize2dr32(input, newsize1, newsize2, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real32), allocatable, intent(inout) :: input(:,:)

  !A buffer used to expand the array
  real(kind=real32), allocatable :: buffer(:,:)

  !The inital value for the array to be set to
  real(kind=real32), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)

  !log our request in size change
  request = ((newsize1*newsize2) - (oldsize1*oldsize2)) * storage_size(real32)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize2
    do j=1,newsize1
      buffer(j,i) = initial
    end do
  end do

  !Copy the data over
  do i=1,oldsize2
    do j=1,oldsize1
      buffer(j,i)=input(j,i)
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize2dr32


!resize 3d real32 array
subroutine resize3dr32(input, newsize1, newsize2, newsize3, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real32), allocatable, intent(inout) :: input(:,:,:)

  !A buffer used to expand the array
  real(kind=real32), allocatable :: buffer(:,:,:)

  !The inital value for the array to be set to
  real(kind=real32), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2, newsize3

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2, oldsize3

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j, k

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, newsize3, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)
  oldsize3 = size(input,3)

  !log our request in size change
  request = ((newsize1*newsize2*newsize3) - (oldsize1*oldsize2*oldsize3)) * storage_size(real32)

  !Allocate a buffer with the new array size, and fill it with the old input
  allocate(buffer(newsize1, newsize2, newsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize3
    do j=1,newsize2
      do k=1,newsize1
        buffer(k,j,i) = initial
      end do
    end do
  end do

  !Copy the data over
  do i=1,oldsize3
    do j=1,oldsize2
      do k=1,oldsize1
        buffer(k,j,i)=input(k,j,i)
      end do
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize3dr32

!resize 4d real32 array
subroutine resize4dr32(input, newsize1, newsize2, newsize3, newsize4, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real32), allocatable, intent(inout) :: input(:,:,:,:)

  !A buffer used to expand the array
  real(kind=real32), allocatable :: buffer(:,:,:,:)

  !The inital value for the array to be set to
  real(kind=real32), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2, newsize3, newsize4

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2, oldsize3, oldsize4

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j, k, l

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, newsize3, newsize4, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)
  oldsize3 = size(input,3)
  oldsize4 = size(input,4)

  !log our request in size change
  request = ((newsize1*newsize2*newsize3*newsize4) - (oldsize1*oldsize2*oldsize3*oldsize4)) * storage_size(real32)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2, newsize3, newsize4), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize4
    do j=1,newsize3
      do k=1,newsize2
        do l=1,newsize1
          buffer(l,k,j,i) = initial
        end do
      end do
    end do
  end do

  !Copy the data over
  do i=1,oldsize4
    do j=1,oldsize3
      do k=1,oldsize2
        do l=1,oldsize1
          buffer(l,k,j,i)=input(l,k,j,i)
        end do
      end do
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize4dr32


!resize 1d real64 array
subroutine resize1dr64(input, newsize, initial, ename, fIdxIn)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real64), allocatable, intent(inout) :: input(:)

  !A buffer used to expand the array
  real(kind=real64), allocatable :: buffer(:)

  !The inital value for the array to be set to
  real(kind=real64), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize

  !The initial index
  integer, intent(in), optional :: fIdxIn
  integer :: fIdx

  !A variable to track the old size of the array
  integer :: oldsize

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i

  !Set to non-zero in case of an allocation error
  integer :: error

  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize, initial, ename, fIdx)
    return
  end if

  !get the old size of the array
  oldsize = size(input)

  !log our request in size change
  request = (newsize-oldsize) * storage_size(real64)

  !Allocate a buffer with the new array size
  allocate(buffer(fIdx:newsize), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Copy the data over
  do i=fIdx,oldsize-(1-fIdx)
    buffer(i) = input(i)
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Set the initial values of the buffer
  do i=(oldsize-(1-fIdx)+1),newsize
    buffer(i) = initial
  end do

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize1dr64

!resize 2d real64 array
subroutine resize2dr64(input, newsize1, newsize2, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real64), allocatable, intent(inout) :: input(:,:)

  !A buffer used to expand the array
  real(kind=real64), allocatable :: buffer(:,:)

  !The inital value for the array to be set to
  real(kind=real64), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)

  !log our request in size change
  request = ((newsize1*newsize2) - (oldsize1*oldsize2)) * storage_size(real64)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize2
    do j=1,newsize1
      buffer(j,i) = initial
    end do
  end do

  !Copy the data over
  do i=1,oldsize2
    do j=1,oldsize1
      buffer(j,i)=input(j,i)
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize2dr64


!resize 3d real64 array
subroutine resize3dr64(input, newsize1, newsize2, newsize3, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real64), allocatable, intent(inout) :: input(:,:,:)

  !A buffer used to expand the array
  real(kind=real64), allocatable :: buffer(:,:,:)

  !The inital value for the array to be set to
  real(kind=real64), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2, newsize3

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2, oldsize3

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j, k

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, newsize3, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)
  oldsize3 = size(input,3)

  !log our request in size change
  request = ((newsize1*newsize2*newsize3) - (oldsize1*oldsize2*oldsize3)) * storage_size(real64)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2, newsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize3
    do j=1,newsize2
      do k=1,newsize1
        buffer(k,j,i) = initial
      end do
    end do
  end do

  !Copy the data over
  do i=1,oldsize3
    do j=1,oldsize2
      do k=1,oldsize1
        buffer(k,j,i)=input(k,j,i)
      end do
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize3dr64


!resize 4d real64 array
subroutine resize4dr64(input, newsize1, newsize2, newsize3, newsize4, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real64), allocatable, intent(inout) :: input(:,:,:,:)

  !A buffer used to expand the array
  real(kind=real64), allocatable :: buffer(:,:,:,:)

  !The inital value for the array to be set to
  real(kind=real64), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2, newsize3, newsize4

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2, oldsize3, oldsize4

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j, k, l

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, newsize3, newsize4, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)
  oldsize3 = size(input,3)
  oldsize4 = size(input,4)

  !log our request in size change
  request = ((newsize1*newsize2*newsize3*newsize4) - (oldsize1*oldsize2*oldsize3*oldsize4)) * storage_size(real64)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2, newsize3, newsize4), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize4
    do j=1,newsize3
      do k=1,newsize2
        do l=1,newsize1
          buffer(l,k,j,i) = initial
        end do
      end do
    end do
  end do

  !Copy the data over
  do i=1,oldsize4
    do j=1,oldsize3
      do k=1,oldsize2
        do l=1,oldsize1
          buffer(l,k,j,i)=input(l,k,j,i)
        end do
      end do
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize4dr64

!alloc real128
subroutine alloc1dr128(input, startsize, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real128), allocatable, intent(inout) :: input(:)

  !The inital value for the array to be set to
  real(kind=real128), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k

  request = startsize * storage_size(real128)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize
    input(k) = initial
  end do

  call print_alloc(ename,request)
end subroutine alloc1dr128


subroutine alloc2dr128(input, startsize1, startsize2, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real128), allocatable, intent(inout) :: input(:,:)

  !The inital value for the array to be set to
  real(kind=real128), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l

  request = startsize1 * startsize2 * storage_size(real128)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
     input(k,l) = initial
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc2dr128

subroutine alloc3dr128(input, startsize1, startsize2, startsize3, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real128), allocatable, intent(inout) :: input(:,:,:)

  !The inital value for the array to be set to
  real(kind=real128), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2, startsize3

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l,m

  request = startsize1 * startsize2 * startsize3 * storage_size(real128)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2, startsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      do m=1, startsize3
       input(k,l,m) = initial
      end do
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc3dr128

subroutine alloc4dr128(input, startsize1, startsize2, startsize3, startsize4, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real128), allocatable, intent(inout) :: input(:,:,:,:)

  !The inital value for the array to be set to
  real(kind=real128), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2, startsize3, startsize4

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l,m,n

  request = startsize1 * startsize2 * startsize3 * startsize4 * storage_size(real128)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2, startsize3, startsize4), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      do m=1, startsize3
        do n=1, startsize4
         input(k,l,m,n) = initial
        end do
      end do
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc4dr128

!resize 1d real128 array
subroutine resize1dr128(input, newsize, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real128), allocatable, intent(inout) :: input(:)

  !A buffer used to expand the array
  real(kind=real128), allocatable :: buffer(:)

  !The inital value for the array to be set to
  real(kind=real128), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize

  !A variable to track the old size of the array
  integer :: oldsize

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize = size(input)

  !log our request in size change
  request = (newsize-oldsize) * storage_size(real128)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Copy the data over
  do i=1,oldsize
    buffer(i) = input(i)
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Set the initial values of the buffer
  do i=oldsize+1,newsize
    buffer(i) = initial
  end do

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize1dr128

!resize 2d real128 array
subroutine resize2dr128(input, newsize1, newsize2, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real128), allocatable, intent(inout) :: input(:,:)

  !A buffer used to expand the array
  real(kind=real128), allocatable :: buffer(:,:)

  !The inital value for the array to be set to
  real(kind=real128), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)

  !log our request in size change
  request = ((newsize1*newsize2) - (oldsize1*oldsize2)) * storage_size(real128)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize2
    do j=1,newsize1
      buffer(j,i) = initial
    end do
  end do

  !Copy the data over
  do i=1,oldsize2
    do j=1,oldsize1
      buffer(j,i)=input(j,i)
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize2dr128


!resize 3d real128 array
subroutine resize3dr128(input, newsize1, newsize2, newsize3, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real128), allocatable, intent(inout) :: input(:,:,:)

  !A buffer used to expand the array
  real(kind=real128), allocatable :: buffer(:,:,:)

  !The inital value for the array to be set to
  real(kind=real128), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2, newsize3

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2, oldsize3

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j, k

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, newsize3, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)
  oldsize3 = size(input,3)

  !log our request in size change
  request = ((newsize1*newsize2*newsize3) - (oldsize1*oldsize2*oldsize3)) * storage_size(real128)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2, newsize3), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize3
    do j=1,newsize2
      do k=1,newsize1
        buffer(k,j,i) = initial
      end do
    end do
  end do

  !Copy the data over
  do i=1,oldsize3
    do j=1,oldsize2
      do k=1,oldsize1
        buffer(k,j,i)=input(k,j,i)
      end do
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize3dr128


!resize 4d real128 array
subroutine resize4dr128(input, newsize1, newsize2, newsize3, newsize4, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  real(kind=real128), allocatable, intent(inout) :: input(:,:,:,:)

  !A buffer used to expand the array
  real(kind=real128), allocatable :: buffer(:,:,:,:)

  !The inital value for the array to be set to
  real(kind=real128), intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2, newsize3, newsize4

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2, oldsize3, oldsize4

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j, k, l

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, newsize3, newsize4, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)
  oldsize3 = size(input,3)
  oldsize4 = size(input,4)

  !log our request in size change
  request = ((newsize1*newsize2*newsize3*newsize4) - (oldsize1*oldsize2*oldsize3*oldsize4)) * storage_size(real128)

  !Allocate a buffer with the new array size, and fill it with the old input
  allocate(buffer(newsize1, newsize2, newsize3, newsize4), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize4
    do j=1,newsize3
      do k=1,newsize2
        do l=1,newsize1
          buffer(l,k,j,i) = initial
        end do
      end do
    end do
  end do

  !Copy the data over
  do i=1,oldsize4
    do j=1,oldsize3
      do k=1,oldsize2
        do l=1,oldsize1
          buffer(l,k,j,i)=input(l,k,j,i)
        end do
      end do
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize4dr128

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                          CHARACTERS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!alloc char 1d
subroutine alloc1dc(input, strlen, e_index, initial, ename, fIdxIn)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  character(len=:), allocatable, intent(inout) :: input(:)

  !The inital value for the array to be set to
  character(len=*), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: e_index
  integer, intent(in) :: strlen

  !The initial index
  integer, intent(in), optional :: fIdxIn
  integer :: fIdx

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k

  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if

  request = (e_index-fIdx+1) * strlen * CHARACTER_STORAGE_SIZE

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(character(strlen) :: input(fIdx:e_index), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=fIdx, e_index
    input(k) = initial
  end do

  call print_alloc(ename,request)
end subroutine alloc1dc


!alloc char 2d
subroutine alloc2dc(input, strlen, startsize1, startsize2, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  character(len=:), allocatable, intent(inout) :: input(:,:)

  !The inital value for the array to be set to
  character(len=*), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1
  integer, intent(in) :: startsize2
  integer, intent(in) :: strlen

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l

  request = startsize1 * startsize2 * strlen * CHARACTER_STORAGE_SIZE

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(character(strlen) :: input(startsize1,startsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
      input(k,l) = initial
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc2dc

!resize 1d character array
subroutine resize1dc(input, strlen, eIdx, initial, ename, fIdxIn)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  character(len=:), allocatable, intent(inout) :: input(:)

  !The inital value for the array to be set to
  character(len=*), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: eIdx
  integer, intent(in) :: strlen

  !The initial index
  integer, intent(in), optional :: fIdxIn
  integer :: fIdx

  !A buffer used to expand the array
  character(len=:), allocatable :: buffer(:)

  !A variable to track the old size of the array
  integer :: oeIdx

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i

  !Set to non-zero in case of an allocation error
  integer :: error

  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
      call alloc(input, strlen, eIdx, initial, ename, fIdx)
    return
  end if

  !get the old size of the array
  oeIdx = size(input)+fIdx-1
  ! write(lout,"(3(a,i4))") "Resize: old end is ",oeIdx,", new end is ",eIdx,", first index is ",fIdx

  !log our request in size change
  request = (eIdx-oeIdx) * strlen * CHARACTER_STORAGE_SIZE

  !Allocate a buffer with the new array size
  allocate(character(strlen) :: buffer(fIdx:eIdx), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

 !copy old values across
  do i=fIdx,oeIdx
    buffer(i) = input(i)
  end do

  !Set the initial values of the buffer
  if(eIdx > oeIdx) then
    do i=oeIdx+1, eIdx
      buffer(i) = initial
    end do
  end if

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize1dc


!resize 2d character array
subroutine resize2dc(input, strlen, newsize1, newsize2, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  character(len=:), allocatable, intent(inout) :: input(:,:)

  !The inital value for the array to be set to
  character(len=*), intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: newsize1
  integer, intent(in) :: newsize2
  integer, intent(in) :: strlen

  !A buffer used to expand the array
  character(len=:), allocatable :: buffer(:,:)

  !A variable to track the old size of the array
  integer :: oldsize1
  integer :: oldsize2

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i,j

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
      call alloc(input, strlen, newsize1, newsize2, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)

  !log our request in size change
  request = ((newsize1*newsize2)-(oldsize1*oldsize2)) * strlen * CHARACTER_STORAGE_SIZE

  !Allocate a buffer with the new array size
  allocate(character(strlen) :: buffer(newsize1,newsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Set the initial values of the buffer
  do i=1,newsize2
    do j=1,newsize1
      buffer(j,i) = initial
    end do
  end do

 !copy old values across
  do i=1,oldsize2
    do j=1,oldsize1
      buffer(j,i) = input(j,i)
    end do
  end do

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize2dc

! ==================================================================== !
!  STRINGS
! ==================================================================== !

subroutine alloc1ds(input, e_index, initial, ename, fIdxIn)
  
  use strings
  
  implicit none
  
  character(len=*),              intent(in)    :: ename
  type(string),     allocatable, intent(inout) :: input(:)
  type(string),                  intent(in)    :: initial
  integer,                       intent(in)    :: e_index
  integer,          optional,    intent(in)    :: fIdxIn
  
  integer :: fIdx
  integer :: error
  integer(kind=int64) :: request
  integer :: k

  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  request = (e_index-fIdx+1) !*storage_size(type(string))

  ! Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if
  
  ! Do the allocation
  allocate(input(fIdx:e_index), stat=error)
  
  ! Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if
  
  ! Log the number of allocated bits
  allocated_bits = allocated_bits + request
  
  !Initialise the array
  do k=fIdx, e_index
    input(k) = initial
  end do
    
  call print_alloc(ename,request)
  
end subroutine alloc1ds

subroutine resize1ds(input, eIdx, initial, ename, fIdxIn)
  
  use strings
  
  implicit none
  
  character(len=*),              intent(in)    :: ename
  type(string),     allocatable, intent(inout) :: input(:)
  type(string),     allocatable                :: buffer(:)
  type(string),                  intent(in)    :: initial
  integer,                       intent(in)    :: eIdx
  integer,          optional,    intent(in)    :: fIdxIn

  integer :: fIdx
  integer :: oeIdx
  integer(kind=int64) :: request
  integer :: i
  integer :: error
  
  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if
  
  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, eIdx, initial, ename, fIdx)
    return
  end if
  
  ! Get the old end index of the array
  oeIdx = size(input)+fIdx-1
  
  !log our request in size change
  request = (eIdx-oeIdx) !*storage_size(type(string))
  
  !Allocate a buffer with the new array size
  allocate(buffer(fIdx:eIdx), stat=error)
  
  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if
  
  !Copy the data over
  do i=fIdx,oeIdx
    buffer(i) = input(i)
  end do
  
  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)
  
  !Set the initial values of the buffer
  if(eIdx > oeIdx) then
    do i=oeIdx+1,eIdx
      buffer(i) = initial
    end do
  end if
  
  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)
  
end subroutine resize1ds
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                          LOGICAL
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine alloc1dl(input, e_index, initial, ename, fIdxIn)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  logical, allocatable, intent(inout) :: input(:)

  !The inital value for the array to be set to
  logical, intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: e_index

  !The initial index
  integer, intent(in), optional :: fIdxIn
  integer :: fIdx

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k

  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if

  request = e_index-fIdx+1 !* storage_size(logical)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(fIdx:e_index), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=fIdx, e_index
    input(k) = initial
  end do

  call print_alloc(ename,request)
end subroutine alloc1dl

subroutine alloc2dl(input, startsize1, startsize2, initial, ename)
  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  logical, allocatable, intent(inout) :: input(:,:)

  !The inital value for the array to be set to
  logical, intent(in) :: initial

  !The initial size to allocate
  integer, intent(in) :: startsize1, startsize2

  !Set to non-zero in case of an allocation error
  integer :: error

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: k,l

  request = startsize1 * startsize2 !* storage_size(logical)

  !Check that we are not already allocated
  if(allocated(input) .eqv. .TRUE.) then
    write(lout,*) 'ERROR: input array is already allocated for: ', ename
    stop
  end if

  !Do the allocation
  allocate(input(startsize1, startsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Log the number of allocated bits
  allocated_bits = allocated_bits + request

  !Initialise the array
  do k=1, startsize1
    do l=1, startsize2
     input(k,l) = initial
    end do
  end do

  call print_alloc(ename,request)
end subroutine alloc2dl

!resize 1d logical array
subroutine resize1dl(input, eIdx, initial, ename, fIdxIn)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  logical, allocatable, intent(inout) :: input(:)

  !A buffer used to expand the array
  logical, allocatable :: buffer(:)

  !The inital value for the array to be set to
  logical, intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: eIdx

  !The initial index
  integer, intent(in), optional :: fIdxIn
  integer :: fIdx

  !A variable to track the old size of the array
  integer :: oeIdx

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i

  !Set to non-zero in case of an allocation error
  integer :: error

  if(present(fIdxIn)) then
    fIdx = fIdxIn
  else
    fIdx = 1
  end if

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, eIdx, initial, ename, fIdx)
    return
  end if

  !get the old size of the array
  oeIdx = size(input)+fIdx-1

  !log our request in size change
  request = (eIdx-oeIdx) !* storage_size(logical)

  !Allocate a buffer with the new array size
  allocate(buffer(fIdx:eIdx), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Copy the data over
  do i=fIdx,oeIdx
    buffer(i) = input(i)
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Set the initial values of the buffer
  if(eIdx > oeIdx) then
    do i=oeIdx+1,eIdx
      buffer(i) = initial
    end do
  end if

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize1dl

!resize 2d logical array
subroutine resize2dl(input, newsize1, newsize2, initial, ename)

  implicit none

  !The name of the variable being assigned - for error printing
  character(len=*), intent(in) :: ename

  !The input variable to allocated
  logical, allocatable, intent(inout) :: input(:,:)

  !A buffer used to expand the array
  logical, allocatable :: buffer(:,:)

  !The inital value for the array to be set to
  logical, intent(in) :: initial

  !The new size that the array should be set to
  integer, intent(in) :: newsize1, newsize2

  !A variable to track the old size of the array
  integer :: oldsize1, oldsize2

  !To keep track of the requested allocation size
  integer(kind=int64) :: request

  !Loop variable
  integer :: i, j

  !Set to non-zero in case of an allocation error
  integer :: error

  if(allocated(input) .neqv. .TRUE.) then
    write(lout,*) 'INFO: array ', ename, ' is not allocated.'
    call alloc(input, newsize1, newsize2, initial, ename)
    return
  end if

  !get the old size of the array
  oldsize1 = size(input,1)
  oldsize2 = size(input,2)

  !log our request in size change
  request = ((newsize1*newsize2) - (oldsize1*oldsize2)) !* storage_size(logical)

  !Allocate a buffer with the new array size
  allocate(buffer(newsize1, newsize2), stat=error)

  !Print and exit if we have an error
  if(error.ne.0) then
    call alloc_error(ename, error, request)
  end if

  !Set the initial values of the buffer
  do i=1,newsize2
    do j=1,newsize1
      buffer(j,i) = initial
    end do
  end do

  !Copy the data over
  do i=1,oldsize2
    do j=1,oldsize1
      buffer(j,i)=input(j,i)
    end do
  end do

  !update the number of bits allocated (can be negative)
  allocated_bits = allocated_bits + request
  call print_alloc(ename,request)

  !Do a pointer swap and deallocate the buffer
  call move_alloc(buffer,input)

end subroutine resize2dl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                          DEALLOCATIONS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!These exist to keep track of the total allocated memory
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
