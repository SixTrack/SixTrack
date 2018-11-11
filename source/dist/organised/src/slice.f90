
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
      real(kind=real64), dimension(2500,2500) :: distribution
      double precision, external :: test, normalcdfinv
      real(kind=real64) canon(1:6), compare(1:6)
      real(kind=real64) identity2(1:6,1:6)
      real(kind=real64) momentum, mass, one, e1,e2, e3, betx1,zero, angle, energy, ppf, pia2, six
      real(kind=real64), dimension(6, 6) :: identity, results, testm, tas, emit
      call readMatrixFromFile(tas)
      e1 = 1.0d0
      e2 = 1.0d0
      e3 = 0.03d0 
      emit(:,:) = 0.0
      physcord(:) = 0.0
      emit(1,1) = e1
      emit(2,2) = e1
      emit(3,3) = e2
      emit(4,4) = e3
      pia2 = 2.00d0*3.1415
      maxa = 6
      one =1.00000d0
      zero = 0.0d0
      momentum = 4000.0
      mass = 1000.0
      energy = 5000.0
      six = 6.000d0
      !do i=1,630
        !angle = 0.01*i
        !call a2c(e1, angle, e2, angle, e3, angle, tas, physcord)
       ! print *, e1, angle, physcord(1), physcord(2),physcord(3), physcord(4),physcord(5), physcord(6)
       !physcord(:)= 0
      !end do
      print *, tas
      maxa = 6;
      index = 1;
      call initializedistribution(index, maxa)
      call setemittance12(one,one)
      call setmassmom(mass, momentum)
      call setparameter(1,zero,six,50,3);
      call setparameter(2,zero,zero,index,0);
      call setparameter(3,zero,six,50,3);
      call setparameter(4,zero,zero,index,0);
      call setparameter(5,zero,zero,index,0);
      call setparameter(6,zero,zero,index,0);
      call settasmatrix(tas)
      !call dist2sixcoord(distribution)

      !call dist2sixcoord(emit)
      call calcualteInverse
      !call testmatrix(tas)
      call setdeltap(0.01d0)
 
      end program slice
 
