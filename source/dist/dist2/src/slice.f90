
subroutine readMatrixFromFile(matrix)
   ! opening the file for reading
         use, intrinsic :: iso_fortran_env
    real(kind=real64), dimension(36) :: reada
    real(kind=real64), dimension(6, 6) :: matrix
   open (2, file = 'data1.dat', status = 'old')

    
      read(2,*) reada
      matrix=RESHAPE(reada, shape(matrix)  )
   
   close(2)
   
end subroutine readMatrixFromFile

program slice
      use, intrinsic :: iso_fortran_env
      implicit none

      integer MAX
      parameter (MAX = 5)
      external cprintarr
      real arr(1:MAX,1:MAX)
      integer i,j,index, maxa


      real(kind=real64), dimension(6) :: coordinates, physcord
      double precision, external :: test
      real(kind=real64) canon(1:6), compare(1:6)
      real(kind=real64) identity2(1:6,1:6)
      real(kind=real64) momentum, mass, one, e1,e2, e3, betx1,zero, angle, energy
      real(kind=real64), dimension(6, 6) :: identity, results, testm, tas, emit
      call readMatrixFromFile(tas)
      e1 = 1.0d0
      e2 = 2.0d0
      e3 = 0.03d0 
      emit(:,:) = 0.0
      physcord(:) = 0.0
      emit(1,1) = e1
      emit(2,2) = e1
      emit(3,3) = e2
      emit(4,4) = e3

      maxa = 6
      one =1.0d0
      zero = 0.0d0
      momentum = 4000.0
      mass = 3000.0
      energy = 5000.0
      !print *, tas
      do i=1,630
        angle = 0.01*i
        call a2c(e1, angle, e2, angle, e3, angle, tas, physcord)
        print *, e1, angle, physcord(1), physcord(2),physcord(3), physcord(4),physcord(5), physcord(6)
        physcord(:)= 0
      end do
      identity(:,:) = 0
      results(:,:) = 0
      testm(:,:) = 0
      testm(2,3) = one
      identity(1,1) = one
      identity(2,2) = one
      identity(3,3) = one
      identity(4,4) = one
      identity(5,5) = one
      identity(6,6) = one
      results(:,:) = 0

      print *, "valllueee", e2
      !call six2canonical(coordinates, momentum,mass, canon)
      !call canonical2six(canon, momentum, mass, compare)

      !print *, canon
      !print *, compare
      !print *, coordinates
      !print*, identity
    
      !call mtrx_mult(maxa,maxa,maxa, identity, testm, results)
      !print *, "ouuttt"
      !print*, results(1,1),results(2,2),results(3,3), results(2,3)
      !call readtas
! Populate my 2D array with some values.
    !  index = 1
    !  do j = 1,MAX
    !      do i = 1,MAX
    !!          arr(i,j) = index
     !         index = index + 1
      !    enddo
      !enddo
 
      
      !call cprintarr(arr, MAX, MAX)
!# Pass an array slice to C
      !call cprintarr(arr(2:MAX-1,2:MAX-1), MAX-2, MAX-2)
 
      end program slice
 
