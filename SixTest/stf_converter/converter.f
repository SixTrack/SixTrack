      program converter
      character*8 cdate,ctime,progrm
      character*80 sixtit,commen
      integer i,j,k,l,rph,ifipa,ilapa,itopa,imax,icode,numl,turns
      dimension qwc(3),clo(3),clop(3),di0(2),dip0(2)
      dimension ta(6,6)
      double precision dummy,dmmac,dnms,dizu0,dnumlr,sigcor,dpscor

C      open(90,form='unformatted',status='unknown')
      read(1,iostat=ierro) sixtit,commen,cdate,ctime,                   &
     &progrm,ifipa,ilapa,itopa,icode,numl
       imax=itopa
       turns=numl
        write(*,*) 'total particles and number of turns ',imax,' ',turns
       if(ifipa.eq.ilapa) then
          rph=1              !--rph is record per header 
       else 
          rph=2
       endif
       write(*,*) 'particle per record',rph
       rewind 1
       do 20 i=1,imax,rph
          rewind 1
 10       read(1,iostat=ierro) sixtit,commen,cdate,ctime,               &
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
         open(91-i,form='unformatted',status='unknown')
         write(91-i,iostat=ierro) sixtit,commen,cdate,ctime,            &
     &progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
     &clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0 &
     &(2),dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
     &ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
     &ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
     &ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
     &ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
     &dnumlr,sigcor,dpscor,zero,zero,zero,zero,                         &
     &zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,                &
     &zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,                &
     &zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,                &
     &zero,zero,zero,zero,zero,zero,zero,zero,zero,zero
       else
         j=(i+1)/2
         open(91-j,form='unformatted',status='unknown')
         write(91-j,iostat=ierro) sixtit,commen,cdate,ctime,            &
     &progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
     &clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0 &
     &(2),dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
     &ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
     &ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
     &ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
     &ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
     &dnumlr,sigcor,dpscor,zero,zero,zero,zero,                         &
     &zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,                &
     &zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,                &
     &zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,                &
     &zero,zero,zero,zero,zero,zero,zero,zero,zero,zero
       endif
      rewind 1
      do k=1,imax,rph
         read(1)
      enddo 

      if(rph.eq.1) then
        do l=1,turns
        write(*,*) 'turn number',l
 30     read(1,iostat=ierro) ia,ifipa,b,c,d,e,f,g,h,p
        if(ifipa.ne.i) then
          goto 30
        endif
        write(91-i,iostat=ierro) ia,ifipa,b,c,d,e,f,g,h,p
        enddo
      endif
 
      if(rph.eq.2) then
        do l=1,turns
 40     read(1,iostat=ierro) ia,ifipa,b,c,d,e,f,g,h,p,ilapa,b,           &
     &c1,d1,e1,f1,g1,h1,p1
        if(ifipa.ne.i) then
          goto 40
        endif
        write(91-j,iostat=ierro) ia,ifipa,b,c,d,e,f,g,h,p,ilapa,         &
     &b,c1,d1,e1,f1,g1,h1,p1
        enddo
      endif 
 20   continue
      end program
