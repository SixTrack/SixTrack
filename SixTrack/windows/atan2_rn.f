      double precision function atan2_rn(y,x)
      implicit none
      double precision atan_rn,x,y,pi,pi2
      logical myisnan
      data pi  /3.1415926535897932d0/
      data pi2 /1.5707963267948966d0/
      if (x.eq.0d0) then
         if (y.eq.0d0) then
C Should get me a NaN 
           atan2_rn=atan_rn(y/x)
         else
           atan2_rn=sign(pi2,y)
         endif
      else
        if (y.eq.0d0) then
          if (x.gt.0d0) then
            atan2_rn=0d0
          else
            atan2_rn=pi
          endif
        else          
          atan2_rn=atan_rn(y/x)
          if (x.lt.0d0) then
            atan2_rn=sign(pi,y)+atan2_rn
          endif
        endif
      endif
      end
