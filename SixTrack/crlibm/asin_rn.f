      double precision function asin_rn(x)
      implicit none
      double precision atan_rn,x,pi2
      logical myisnan
      data pi2 /1.5707963267948966d0/
      if (myisnan(x,x)) then
        asin_rn=x
        return
      endif
      if (abs(x).eq.1.0d0) then
        asin_rn=sign(pi2,x)
      else 
!       asin_rn=atan_rn(x/sqrt(1.0d0-x*x))
! Try using (1-x)*(1+x) in case x is very small.........
! or close to 1.....write a test program!!!
        asin_rn=atan_rn(x/sqrt((1.0d0-x)*(1.0d0+x)))
      endif
      end
