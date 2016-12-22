      program new_converter
!--   July 2016
!--   Author: Vikas Gupta, email-id:vikasgupta.rmps@gmail.com
!--   fortran program to convert ouput data of stf version of sixtrack
!--   from binary to ASCII format
!--   input=fort.90 and output=output.dat 
      character*8 cdate,ctime,progrm
      character*80 sixtit,commen
      integer i,k,l,rph,ifipa,ilapa,itopa,imax,icode,numl,turns
      dimension qwc(3),clo(3),clop(3),di0(2),dip0(2)
      dimension ta(6,6)
      double precision dummy,dmmac,dnms,dizu0,dnumlr,sigcor,dpscor
       open(unit=12,file='output.dat',form='formatted',status='unknown')
      read(90,iostat=ierro) sixtit,commen,cdate,ctime,                  &
     &progrm,ifipa,ilapa,itopa,icode,numl
       imax=itopa
       turns=numl+1
       if(ifipa.eq.ilapa) then
          rph=1              !--rph is record per header 
       else 
          rph=2
       endif
       rewind 90
       do 20 i=1,imax,rph
          rewind 90
 10       read(90,iostat=ierro) sixtit,commen,cdate,ctime,              &
     &progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
     &clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0 &
     &(2),dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
     &ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
     &ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
     &ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
     &ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
     &dnumlr,sigcor,dpscor
       if(ifipa.ne.i) then
         goto 10
       endif
       if(rph.eq.1) then
         write(12,*) sixtit,commen,cdate,ctime,                         &
     &progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
     &clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0 &
     &(2),dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
     &ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
     &ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
     &ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
     &ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
     &dnumlr,sigcor,dpscor
       else
           write(12,*) sixtit,commen,cdate,ctime,                       &
     &progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
     &clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0 &
     &(2),dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
     &ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
     &ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
     &ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
     &ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
     &dnumlr,sigcor,dpscor
       endif
      rewind 90
      do k=1,imax,rph
         read(90)
      enddo 
      if(rph.eq.1) then
        do l=1,turns
 30     read(90,iostat=ierro) ia,ifipa,b,c,d,e,f,g,h,p
        if(ifipa.ne.i) then
          goto 30
        endif
        write(12,*) ia,ifipa,b,c,d,e,f,g,h,p
        enddo
      endif
 
      if(rph.eq.2) then
        do l=1,turns
 40     read(90,iostat=ierro) ia,ifipa,b,c,d,e,f,g,h,p,ilapa,b,          &
     &c1,d1,e1,f1,g1,h1,p1
        if(ifipa.ne.i) then
          goto 40
        endif
          write(12,*) ia,ifipa,b,c,d,e,f,g,h,p,ilapa,                    &
     &b,c1,d1,e1,f1,g1,h1,p1
        enddo
      endif 
 20   continue
      end program
