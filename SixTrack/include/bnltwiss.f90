if(n.eq.1.and.lhc.eq.9) then
  totals=totals+strack(i)
  sampl(i)=totals
#ifdef CRLIBM
  ! Use dtostr for correct binary decimal conversion
  l1=1
#ifndef BOINC
  ! Use Unit 51
#else
  ! Use Unit 10 and initialise string with header
  ch(l1:l1+10)='SixTwiss  '
  l1=l1+11
#endif
  ! Now do the conversions
  ! A 5-digit integer, followed by 7 real(kind=fPrec) numbers
  write(ch(l1:l1+5),'(i5)') i
  l1=l1+6
  ! and now the 7 real(kind=fPrec)
  ! We return the length of the string (always 24)
  errno=dtostr(sampl(i),ch1)
  ch(l1:l1+errno)=' '//ch1(1:errno)
  l1=l1+errno+1
  errno=dtostr(tbetax(i),ch1)
  ch(l1:l1+errno)=' '//ch1(1:errno)
  l1=l1+errno+1
  errno=dtostr(tbetay(i),ch1)
  ch(l1:l1+errno)=' '//ch1(1:errno)
  l1=l1+errno+1
  errno=dtostr(talphax(i),ch1)
  ch(l1:l1+errno)=' '//ch1(1:errno)
  l1=l1+errno+1
  errno=dtostr(talphay(i),ch1)
  ch(l1:l1+errno)=' '//ch1(1:errno)
  l1=l1+errno+1
  errno=dtostr(torbx(i),ch1)
  ch(l1:l1+errno)=' '//ch1(1:errno)
  l1=l1+errno+1
  errno=dtostr(torby(i),ch1)
  ch(l1:l1+errno)=' '//ch1(1:errno)
  l1=l1+errno+1
#ifndef BOINC
  ! write string to 51
  write(51,'(a)') ch(1:l1-1)
#else
  ! write string to 10
  write(10,'(a)') ch(1:l1-1)
#endif
#else
#ifndef BOINC
  write(51,'(i5,(1x,f15.10),6(1x,f20.13))')     &
#else
  write(10,'(a10,i5,(1x,f15.10),6(1x,f20.13))') &
    'SixTwiss  ',                               &
#endif
    i,sampl(i),tbetax(i),tbetay(i),talphax(i),talphay(i),torbx(i),torby(i)
#endif

#ifdef BOINC
  bnlrec=bnlrec+1
  endfile (10,iostat=ierro)
  backspace (10,iostat=ierro)
#else
  endfile (51,iostat=ierro)
  backspace (51,iostat=ierro)
#endif
endif
