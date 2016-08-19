      program old_converter
!--   July 2016
!--   Author: Vikas Gupta,email-id: vikasgupta.rmps@gmail.com
!--   fortran program to convert  ouput data of non-stf version of sixtrack
!--   from binary to ASCII format into one single older.dat file.
!--   this requires all the fort.# files and not just fort.90
      character*8 cdate,ctime,progrm
      character*80 sixtit,commen
      integer i,rph,ifipa,ilapa,itopa,imax,icode,numl
      dimension qwc(3),clo(3),clop(3),di0(2),dip0(2)
      dimension ta(6,6)
      double precision dummy,dmmac,dnms,dizu0,dnumlr,sigcor,dpscor
      read(90,iostat=ierro) sixtit,commen,cdate,ctime,                   &
     &progrm,ifipa,ilapa,itopa,icode,numl
       if(ifipa.ne.ilapa) then
       imax=itopa/2
       rph=2
       else
       imax=itopa
       rph=1
       endif
      open(unit=13,file='older.dat',form='formatted',status='unknown')
       rewind 90
       do 20 i=1,imax
       read(91-i,iostat=ierro) sixtit,commen,cdate,ctime,               &
     &progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
     &clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0 &
     &(2),dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
     &ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
     &ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
     &ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
     &ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
     &dnumlr,sigcor,dpscor
      write(13,*) sixtit,                                               &
     &progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
     &clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0 &
     &(2),dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
     &ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
     &ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
     &ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
     &ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
     &dnumlr,sigcor,dpscor
      if(rph.eq.2) then 
 35     read(91-i,iostat=ierro) ia,ifipa,b,c,d,e,f,g,h,p,ilapa,b,       &
     &c1,d1,e1,f1,g1,h1,p1
        if(ierro.eq.0) then
          write(13,*) ia,ifipa,b,c,d,e,f,g,h,p,ilapa,                   &
     &b,c1,d1,e1,f1,g1,h1,p1
          goto 35
        endif
      endif
      if(rph.eq.1) then
 40     read(91-i,iostat=ierro) ia,ifipa,b,c,d,e,f,g,h,p
        if(ierro.eq.0) then
          write(13,*) ia,ifipa,b,c,d,e,f,g,h,p
          goto 40
         endif
       endif
 20   continue
      end program
