      double precision function acos_rn(x)
      implicit none
      double precision atan_rn,x,pi,pi2
      logical myisnan
      data pi  /3.1415926535897932d0/
      data pi2 /1.5707963267948966d0/
      if (myisnan(x,x)) then
        acos_rn=x
      elseif (abs(x).eq.0.0d0) then
        acos_rn=pi2
      else
!       acos_rn=atan_rn(sqrt(1.0d0-x*x)/x)
! Try using (1-x)*(1+x) in case x is very small.........
! or close to 1.....write a test program!!!
         acos_rn=atan_rn(sqrt((1.0d0-x)*(1.0d0+x))/x)
        if (x.lt.0.0d0) then
          acos_rn=pi+acos_rn
        endif
      endif
      end
