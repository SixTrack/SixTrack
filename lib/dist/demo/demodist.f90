
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
      real(kind=real64), dimension(6) :: acoord, physical 
      real(kind=real64), dimension(2500,6) :: distribution1
      real(kind=real64), dimension(10002) :: x,xp,y, yp, sigma, delta
      real(kind=real64) momentum, mass, one, e1,e2, e3, dp, betx1, zero, pia2, six
      real(kind=real64), dimension(6, 6) :: identity, results, testm, tas
      real(kind=real64) dim
      real(kind=real64) emitt(3)
      
      call readMatrixFromFile(tas)
      dim = 6
      e1 = 3.0d0
      e2 = 1.0d0
      e3 = 1.0d0 
      dp = 0.100d0
      pia2 = 2.00d0*4.D0*DATAN(1.D0)
      zero = 0.0d0
      momentum = 4000.0
      mass = 938.0
      one =1d0
      six = 6.000d0

      
      acoord(1) = 7.0;
      acoord(2) = 0;
      acoord(3) = 5.0;
      acoord(4) = 0;
      acoord(5) = 1.0;
      acoord(6) = 0!1.567718E+00 + 0.005!-pia2/4.0d0+0.0038+100;
!-1.571020E+00

      call initializedistribution(3, 6)
      ! Set the tas matrix 
      call settasmatrix(tas)  
      ! Set the emittance
      call setemittance12(e1,e2)
      call setdeltap(dp)
      !call setemittance3(e3)
      call setmassmom(mass, momentum)
      call action2canonical(acoord,physical)
      print *, "physical", physical
     ! call canonical2emittance(physical, emitt)
     ! print *, "emittance", emitt
      ! Set the parameters to generate the distribution
      !1. (1=x, 2=x', 3=y 4=y', 5=ds 6=dp)
      !2. Start value of scan
      !3. End value of scan
      !4. Number of points
      !5. type of spacing 0 - constant, 3 - linear spacing  
     
 
      end program demodist
 
