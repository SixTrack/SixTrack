      integer function dtostr(x,results)
! Uses the dtoa_c.c version of dtoa via the dtoaf.c interface in
! crlibm
      implicit none
      double precision x
      character*(*) results
      integer dtoaf 
      integer ilen,mode,ndigits,decpoint,mysign
      integer i,l,d,e
      character*1 str(999)
      character*24 lstr
      character*3 e3
      mode=2
      ndigits=17
      ilen=dtoaf(x,mode,ndigits,decpoint,mysign,str(1))
      if (ilen.le.0.or.ilen.gt.17) then
! Always returns 17 or less characters as requested
      write (*,10000)
      write (*,*) 'Routine dtoa[f] returned string length ',ilen
      call abend('Error writing fort.10                             ')
10000 format(5x///t10,'++++++++++++++++++++++++'/ t10,                  &
     &'+++++ERROR DETECTED+++++'/ t10,'++++++++++++++++++++++++'/ t10)
! Never returns
      endif
      lstr=' '
      do i=1,ilen
        lstr(i:i)=str(i)
      enddo
! Now try my formatting
      d=decpoint
      e=0
      l=1
      lstr=' '
      if (mysign.ne.0) then
        lstr(l:l)='-'
      endif
      if (decpoint.eq.9999) then
! Infinity or Nan
        do i=1,ilen
          lstr(l+i:l+i)=str(i)
        enddo
      else
! Pad with zeros
        do i=ilen+1,17
          str(i)='0'
        enddo
        if (decpoint.le.0) then
          e=decpoint-1
          d=1
        else
! I am using 17 as decision point to avoid dddd.e+eee
! but rather d.ddde+eee
          if (decpoint.ge.17) then
            e=decpoint-1
            d=1
          else
            d=decpoint
          endif
        endif
! and copy with the decimal point
        do i=1,17
          lstr(l+i:l+i)=str(i)
          if (i.eq.d) then
            l=l+1
            lstr(l+i:l+i)='.'
          endif
        enddo
! and add exponent e+/-nnn
        l=20
        lstr(l:l)='e'
        l=21
        lstr(l:l)='+'
        if (e.lt.0) then
          lstr(l:l)='-'
          e=-e
        endif
        l=22
        write (e3,'(I3.3)') e
        lstr(l:l+2)=e3(1:3)
      endif  
      results=lstr
      dtostr=24
      return
      end
