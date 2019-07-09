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
      use mod_dist 
      implicit none


      integer i, npart, strlen, a0, z0;
      real(kind=real64), dimension(6) :: coordinates 
      real(kind=real64), dimension(6) :: closed
      real(kind=real64), dimension(2500,6) :: distribution1
      real(kind=real64), dimension(10002) :: x,xp,y, yp, sigma, deltap
      real(kind=real64) momentum, mass, one, e1,e2, e3, dp, betx1, zero, pia2, six
      real(kind=real64), dimension(6, 6) :: identity, results, testm, tas
      real(kind=real64) energy0, mass0

      character(len=256) filename, fileout
      
      call readMatrixFromFile(tas)

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
      closed(1)=10
      

      filename = '../format_example1.txt'
      fileout  = 'format1_out.txt'
      ! Initialize 3 distributions with dimenstion 6

      call dist_initializedistribution(3)
      ! Set the tas matrix 
      strlen = LEN_TRIM(filename) 
      call dist_readfile(filename,strlen)
      call dist_getrefpara(energy0, mass0, a0, z0)
      print *, "energy0, mass0, a0, z0 ", energy0, mass0, a0, z0
      
      call dist_addclosedorbit(closed)
      call dist_get6trackcoord(x,xp,y,yp,sigma,deltap, npart)
      do i=1,npart
        print *, x(i),y(i),sigma(i)
      enddo
      strlen = LEN_TRIM(fileout) 
      call dist_writefile(fileout, strlen)

! Reads the track file (added mass 2 it)
      call dist_initializedistribution(2)
      filename = '../out_test-track_ap_collimator.obs0001.p0001'
      strlen = LEN_TRIM(filename) 
      call dist_readfile(filename,strlen)

      fileout = 'fromMadx_out.txt'
      strlen = LEN_TRIM(fileout) 
      call dist_writefile(fileout, strlen)

    !Distribution 2: a matched distribution

    ! Change the distribution to 1
!   call setdistribution(1)
!    call settasmatrix(tas)
!    call setemittance12(e1,e2)
!    call setemittance3(e3)
!    call setmassmom(mass, momentum)!

    !call setparameter(1,zero,one,100,6);
    !call setparameter(2,zero,pia2,100,4);
    !call setparameter(3,zero,one,100,6);
    !call setparameter(4,zero,pia2,100,4);
    !call setparameter(5,zero,zero,1,0);
    !call setparameter(6,zero,zero,1,0);



!e1 from fort.10 1.999999999648927940E+00
!e2 from fort.10 9.999999998426114534E-01
 
      end program demodist