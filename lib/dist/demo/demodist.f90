subroutine readMatrixFromFile(matrix)
   ! opening the file for reading
    use, intrinsic :: iso_fortran_env
    real(kind=real64), dimension(36) :: reada
    real(kind=real64), dimension(6, 6) :: matrix
    open (2, file = '../data/tasm.txt', status = 'old')

    
      read(2,*) reada
      matrix=RESHAPE(reada, shape(matrix)  )
   close(2)
   
end subroutine readMatrixFromFile

program demodist
      use, intrinsic :: iso_fortran_env
      implicit none

      integer i ;
      real(kind=real64), dimension(6) :: coordinates 
      real(kind=real64), dimension(2500,6) :: distribution1
      real(kind=real64), dimension(10002) :: x,xp,y, yp, sigma, delta
      real(kind=real64) momentum, mass, one, e1,e2, e3, dp, betx1, zero, pia2, six
      real(kind=real64), dimension(6, 6) :: identity, results, testm, tas
      real(kind=real64) dim
      
      call readMatrixFromFile(tas)
      dim = 6
      e1 = 1.0d0
      e2 = 2.0d0
      e3 = 0.03d0
      dp = 0.0001d0
      pia2 = 2.00d0*3.1415
      zero = 0.0d0
      momentum = 4000.0
      mass = 938.0
      one =1d0
      six = 6.000d0

      


      ! Initialize 3 distributions with dimenstion 6

      call initializedistribution(3, 6)
      ! Set the tas matrix 
      call settasmatrix(tas)  
      ! Set the emittance
      call setemittance12(e1,e2)

      call setemittance3(e3)
      ! Set mass and momentum 
      call setmassmom(mass, momentum)
    
      ! Set the parameters to generate the distribution
      !1. (1=x, 2=x', 3=y 4=y', 5=ds 6=dp)
      !2. Start value of scan
      !3. End value of scan
      !4. Number of points
      !5. type of spacing 0 - constant, 3 - linear spacing  
      call setparameter(1,zero,six,100,3);
      call setparameter(2,zero,zero,1,0);
      call setparameter(3,zero,six,100,3);
      call setparameter(4,zero,zero,1,0);
      call setparameter(5,zero,zero,1,0);
      call setparameter(6,zero,zero,1,0);

      !Print the settings mainly for debugging
      !call printdistsettings()
      !Get the filled vectors 
      !call getcoordvectors(x,xp,y, yp, sigma, delta)
  


    !Distribution 2: a matched distribution

    ! Change the distribution to 1
    call setdistribution(1)
    call settasmatrix(tas)
    call setemittance12(e1,e2)
    call setemittance3(e3)
    call setmassmom(mass, momentum)

    call setparameter(1,zero,one,10000,6);
    call setparameter(2,zero,pia2,10000,4);
    call setparameter(3,zero,one,10000,6);
    call setparameter(4,zero,pia2,10000,4);
    call setparameter(5,zero,zero,1,0);
    call setparameter(6,zero,zero,1,0);


    do i = 0, 9999!00 
        call getcoord(coordinates,i)
        
        print *, i,coordinates
    enddo
    
!e1 from fort.10 1.999999999648927940E+00
!e2 from fort.10 9.999999998426114534E-01
 
      end program demodist