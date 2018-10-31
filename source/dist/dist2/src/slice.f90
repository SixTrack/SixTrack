subroutine readtas
      use, intrinsic :: iso_fortran_env
     ! integer, parameter :: sp = REAL32
      !integer, parameter :: real64a = REAL64
      real(kind=real64) ta_tmp(6,6)
      real(kind=real64) dmmac_tmp,dnms_tmp,dizu0_tmp,dnumlr_tmp,sigcor_tmp,dpscor_tmp

      !For the actual tracking data
      real(kind=real64) b_tmp,c_tmp,d_tmp,e_tmp,f_tmp,g_tmp,h_tmp,p_tmp
      real c1_tmp,d1_tmp,e1_tmp,f1_tmp,g1_tmp,h1_tmp,p1_tmp

      character(len=80) title(20),chxtit(20),chytit(20)
      character(len=8) cdate,ctime,progrm ! Note: Keep in sync with maincr
      character(len=80) sixtit,commen     ! Note: Keep in sync with mod_common
                                          ! DANGER: If the len changes, CRCHECK will break.
      real qwc_tmp(3), clo_tmp(3), clop_tmp(3)
      character(len=11) hvs
      character(len=8192) ch
      character(len=25) ch1
      integer errno,l1,l2
      logical rErr
      real(kind=real64) dummy64, dam_tmp
      real(kind=real64) di0_tmp(2), dip0_tmp(2)
      real(kind=real64) sigcor64, dpscor64, zero64
      integer i,j,ia,idummy,ierro,ifipa,ihalf,ilapa,ipa,ipa1,itopa,numl, nfile, icode
      nfile =90
      ierro=0
      ia=0

      !OPEN(nfile, action='read',iostat=ierro, file="fort.90", form='unknown', access='unknown',position='asis')
       
      read(nfile,end=510, iostat=ierro) &
     &     sixtit,commen,cdate,ctime,progrm, &
     &     ifipa,ilapa,itopa,icode,numl, &
     &     qwc_tmp(1),qwc_tmp(2),qwc_tmp(3), &
     &     clo_tmp(1),clop_tmp(1),clo_tmp(2),clop_tmp(2), &
     &     clo_tmp(3),clop_tmp(3), &
     &     di0_tmp(1),dip0_tmp(1),di0_tmp(2),dip0_tmp(2), &
     &     dummy64,dummy64, &
     &     ta_tmp(1,1),ta_tmp(1,2),ta_tmp(1,3), &
     &     ta_tmp(1,4),ta_tmp(1,5),ta_tmp(1,6), &
     &     ta_tmp(2,1),ta_tmp(2,2),ta_tmp(2,3), &
     &     ta_tmp(2,4),ta_tmp(2,5),ta_tmp(2,6), &
     &     ta_tmp(3,1),ta_tmp(3,2),ta_tmp(3,3), &
     &     ta_tmp(3,4),ta_tmp(3,5),ta_tmp(3,6), &
     &     ta_tmp(4,1),ta_tmp(4,2),ta_tmp(4,3), &
     &     ta_tmp(4,4),ta_tmp(4,5),ta_tmp(4,6), &
     &     ta_tmp(5,1),ta_tmp(5,2),ta_tmp(5,3), &
     &     ta_tmp(5,4),ta_tmp(5,5),ta_tmp(5,6), &
     &     ta_tmp(6,1),ta_tmp(6,2),ta_tmp(6,3), &
     &     ta_tmp(6,4),ta_tmp(6,5),ta_tmp(6,6), &
     &     dmmac_tmp,dnms_tmp,dizu0_tmp,dnumlr_tmp, &
     &     sigcor_tmp,dpscor_tmp
   
     print *, ierro, icode, numl
     print *, ta_tmp

510 continue 
  print *, "errrorrr"
     end subroutine readtas 
program slice
      use, intrinsic :: iso_fortran_env
      implicit none

      integer MAX
      parameter (MAX = 5)
      external cprintarr
      real arr(1:MAX,1:MAX)
      integer i,j,index, maxa


      real(kind=real64), dimension(6) :: coordinates 
      real(kind=real64) canon(1:6), compare(1:6)
      real(kind=real64) identity2(1:6,1:6)
      real(kind=real64) momentum, mass, one
      real(kind=real64), dimension(6, 6) :: identity, results
      maxa = 6
      one =1.0d0
      momentum = 5000.0
      mass = 3000.0
      coordinates(1) = 1.0 
      coordinates(2) = 2.0  
      coordinates(3) = 3.0  
      coordinates(4) = 4.0 
      coordinates(5) = 10.0
      coordinates(6) = 0.1
      identity(:,:) = 0
      results(:,:) = 0
      identity(1,1) = one
      identity(2,2) = one
      identity(3,3) = one
      identity(4,4) = one
      identity(5,5) = one
      identity(6,6) = one
      results(:,:)=0
      !call six2canonical(coordinates, momentum,mass, canon)
      !call canonical2six(canon, momentum, mass, compare)

      !print *, canon
      !print *, compare
      !print *, coordinates
      print*, identity
    
      call mtrx_mult(maxa,maxa,maxa, identity, identity, results)
      print *, "ouuttt"
      print*, results(1,1),results(2,2),results(3,3)
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
 
