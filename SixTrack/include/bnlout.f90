!GRD-2007 ADDITIONAL OUTPUT FILE FOR MULTI-PARTICLES BEAM-BEAM STUDIES
do j=1,napx
  if(pstop(nlostp(j))) goto 599
#ifdef DEBUG
!           if (n.ge.990) then
!             write(99,*) 'before j ',j,xv(1,j),xv(2,j),yv(1,j),yv(2,j)
!             endfile (99,iostat=ierro)
!             backspace (99,iostat=ierro)
!           endif
#endif
    x_temp=(xv(1,j)-torbx(i-1))*c1m3                             !hr03
    y_temp=(xv(2,j)-torby(i-1))*c1m3                             !hr03
    xp_temp=(yv(1,j)-torbxp(i-1))*c1m3                           !hr03
    yp_temp=(yv(2,j)-torbyp(i-1))*c1m3                           !hr03
    twojx = (tbetax(i-1)*(xp_temp**2)+((two*talphax(i-1))*x_temp)*xp_temp)+((one+talphax(i-1)**2)/tbetax(i-1))*(x_temp**2)
    twojy = (tbetay(i-1)*(yp_temp**2)+((two*talphay(i-1))*y_temp)*yp_temp)+((one+talphay(i-1)**2)/tbetay(i-1))*(y_temp**2)
    twojr = sqrt(twojx**2+twojy**2)
    if(n.eq.1) then
      if(j.eq.1) then
        limit_twojx = 25.0_fPrec*(2.5e-6_fPrec/(e0/pma))                   !hr03
        limit_twojy = 25.0_fPrec*(2.5e-6_fPrec/(e0/pma))                   !hr03
        limit_twojr = 25.0_fPrec*(2.5e-6_fPrec/(e0/pma))                   !hr03
      endif
    endif
    if(twojr.le.limit_twojr) then
      sumtwojx=sumtwojx+twojx
      sumtwojy=sumtwojy+twojy
      sumsquarex=sumsquarex+x_temp**2
      sumsquarey=sumsquarey+y_temp**2
      n_nocut=n_nocut +1
    else
      n_cut=n_cut+1
#ifndef BOINC
      write(53,'(i8,1x,i8)') n,namepart(j)
#else
      write(10,'(a10,i8,1x,i8)') 'lostID    ',n,namepart(j)
#endif
#ifdef CR
#ifndef BOINC
      endfile (53,iostat=ierro)
      backspace (53,iostat=ierro)
      bllrec=bllrec+1
#else
      endfile (10,iostat=ierro)
      backspace (10,iostat=ierro)
      bnlrec=bnlrec+1
#endif
#endif
    endif
599 continue
  enddo
#ifdef DEBUG
!     write(99,*) 'after  update bnl n ',n
!     write(99,*)                                                       &
!    &n_cut,                                                            &
!    &n_nocut,                                                          &
!    &sumsquarex,                                                       &
!    &sumsquarey,                                                       &
!    &sumtwojx,                                                         &
!    &sumtwojy,                                                         &
!    &limit_twojx,limit_twojy,limit_twojr,                              &
!    &totals,                                                           &
!    &(namepart(j),j=1,napx)
!     endfile (99,iostat=ierro)
!     backspace (99,iostat=ierro)
#endif
  if(mod(n,nwr(3)).eq.0) then
    write(lout,*) 'dumping stats at turn number ',n
#ifdef CRLIBM
    ! Use dtostr for correct binary decimal conversion
    l1=1
#ifndef BOINC
    ! use Unit 52 
#else
    ! use Unit 10 and initialise string with header
    ch(l1:l1+10)='output    '           
    l1=l1+11
#endif
    ! Now do conversions
    ! First the 3 integers using internal read
    write(ch(l1:l1+8),'(i8)') n
    l1=l1+9
    write(ch(l1:l1+9),'(1x,i8)') n_cut
    l1=l1+10
    write(ch(l1:l1+9),'(1x,i8)') n_nocut
    l1=l1+10
    ! and now the four real(kind=fPrec)
    ! We return the length of the string (always 24)
    errno=dtostr(sumsquarex,ch1)
    ch(l1:l1+errno)=' '//ch1(1:errno)
    l1=l1+errno+1
    errno=dtostr(sumsquarey,ch1)
    ch(l1:l1+errno)=' '//ch1(1:errno)
    l1=l1+errno+1
    errno=dtostr(sumtwojx,ch1)
    ch(l1:l1+errno)=' '//ch1(1:errno)
    l1=l1+errno+1
    errno=dtostr(sumtwojy,ch1)
    ch(l1:l1+errno)=' '//ch1(1:errno)
    l1=l1+errno+1
#ifndef BOINC
    ! write string to 52
    write(52,'(a)') ch(1:l1-1)
#else
    ! write string to 10
    write(10,'(a)') ch(1:l1-1)
#endif
#else
#ifndef BOINC
    write(52,'(i8,2(1x,i8),4(1x,e15.8))')                       &
#else
    write(10,'(a10,i8,2(1x,i8),4(1x,e15.8))')                   &
      'output    ',&
#endif
      n,n_cut,n_nocut, &
      sumsquarex,      &
      sumsquarey,      &
      sumtwojx,        &
      sumtwojy
#endif

#ifdef CR
#ifndef BOINC
    endfile (52,iostat=ierro)
    backspace (52,iostat=ierro)
#else
    endfile (10,iostat=ierro)
    backspace (10,iostat=ierro)
#endif
    bnlrec=bnlrec+1
#endif
    sumsquarex=zero
    sumsquarey=zero
    sumtwojx=zero
    sumtwojy=zero
    n_cut=0
    n_nocut=0
  endif
