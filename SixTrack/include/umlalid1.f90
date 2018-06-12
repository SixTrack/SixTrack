ibb=ibb+1
if(ibb.gt.nbb) call prror(102)
imbb(i)=ibb
!FOX  YP(1)=Y(1)*(ONE+DPDA)/MTCDA ;
!FOX  YP(2)=Y(2)*(ONE+DPDA)/MTCDA ;
!FOX  DPDA1=DPDA*C1E3 ;
!FOX  DELTAS=SIGMDA/RV ;
call dacop(x(1),damap(1))
call dacop(yp(1),damap(2))
call dacop(x(2),damap(3))
call dacop(yp(2),damap(4))
do j=1,2
  ii=2*j
  call dapek(damap(ii-1),jj,c(j))
  call dapek(damap(ii),jj,cp(j))
enddo
call dacsu(damap(1),c(1),damap(1))
call dacsu(damap(2),cp(1),damap(2))
call dacsu(damap(3),c(2),damap(3))
call dacsu(damap(4),cp(2),damap(4))
if(ndimf.eq.3) then
  call dacop(deltas,damap(5))
  call dacop(dpda1,damap(6))
  call dapek(damap(5),jj,c(3))
  call dapek(damap(6),jj,cp(3))
  call dacsu(damap(5),c(3),damap(5))
  call dacsu(damap(6),cp(3),damap(6))
  if(iflag2.eq.1.and.ithick.eq.1) then
    call dacct(damap,nvar,corrnew,nvar,damap,nvar)
  endif
endif
call dainv(damap,nvar,damapi,nvar)
call dacct(damap,nvar,aa2,nvar,aa2r,nvar)
call dacct(damap,nvar,damap1,nvar,damap,nvar)
call dacct(damap,nvar,damapi,nvar,damap,nvar)
#ifdef DEBUG
!         write(*,*) 'bbcu set in umlalid1'
!     call warr('umlid1bbcu',bbcu(ibb,1),ibb,1,0,0)
!     call warr('umlid1bbcu',bbcu(ibb,2),ibb,2,0,0)
!     call warr('umlid1bbcu',bbcu(ibb,3),ibb,3,0,0)
#endif
  do ii=1,3
    call damul(damap(i4(ii,1)),damap(i4(ii,2)),angno)
    call averaged(angno,aa2r,.false.,angno,rv)
    do j=1,ndimf
      j1=2*j
      jj(j1-1)=1
      jj(j1)=1
      call dapek(angno,jj,angnoe(j))
      jj(j1-1)=0
      jj(j1)=0
    enddo
    if(ndimf.eq.3) then
      bbcu(ibb,ii)=two*((emitx*angnoe(1)+emity*angnoe(2))+emitz*angnoe(3))
#ifdef DEBUG
!     call warr('umlid1bbcii',bbcu(ibb,ii),ibb,ii,1,0)
#endif
    else
      bbcu(ibb,ii)=two*(emitx*angnoe(1)+emity*angnoe(2))
#ifdef DEBUG
!     call warr('umlid1bbcii',bbcu(ibb,ii),ibb,ii,2,0)
#endif
    endif
  enddo
  if (beam_expflag .eq. 0) then !Old-style input
    if(parbe(ix,2).gt.zero) then
      do ii=4,10
        call damul(damap(i4(ii,1)),damap(i4(ii,2)),angno)
        call averaged(angno,aa2r,.false.,angno, rv)
        do j=1,ndimf
          j1=2*j
          jj(j1-1)=1
          jj(j1)=1
          call dapek(angno,jj,angnoe(j))
          jj(j1-1)=0
          jj(j1)=0
        enddo
        if(ndimf.eq.3) then
          bbcu(ibb,ii) = two * ((emitx*angnoe(1)+emity*angnoe(2))+emitz*angnoe(3))
        else
          bbcu(ibb,ii)=two*(emitx*angnoe(1)+emity*angnoe(2))
        endif
      enddo
    endif
    if(lhc.eq.1) then
      dummy=bbcu(ibb,1)
      bbcu(ibb,1)=bbcu(ibb,2)
      bbcu(ibb,2)=dummy
      dummy=bbcu(ibb,4)
      bbcu(ibb,4)=bbcu(ibb,9)
      bbcu(ibb,9)=dummy
      dummy=bbcu(ibb,5)
      bbcu(ibb,5)=bbcu(ibb,7)
      bbcu(ibb,7)=dummy
      dummy=bbcu(ibb,6)
      bbcu(ibb,6)=bbcu(ibb,10)
      bbcu(ibb,10)=dummy
    endif
    if(lhc.eq.2) then
      bbcu(ibb,1)=bbbx(ix)
      bbcu(ibb,2)=bbby(ix)
      bbcu(ibb,3)=bbbs(ix)
    endif

  !Indentation break, sorry :(

else if (beam_expflag .eq. 1) then !New style input
  if(parbe(ix,2).gt.zero) then
     bbcu(ibb,1)=parbe(ix,7)
     bbcu(ibb,4)=parbe(ix,8)
     bbcu(ibb,6)=parbe(ix,9)
     bbcu(ibb,2)=parbe(ix,10)
     bbcu(ibb,9)=parbe(ix,11)
     bbcu(ibb,10)=parbe(ix,12)
     bbcu(ibb,3)=parbe(ix,13)
     bbcu(ibb,5)=parbe(ix,14)
     bbcu(ibb,7)=parbe(ix,15)
     bbcu(ibb,8)=parbe(ix,16)
  endif
  if(parbe(ix,2).eq.zero) then
     bbcu(ibb,1)=parbe(ix,1)
     bbcu(ibb,2)=parbe(ix,3)
  endif
else
   write(lout,'(a)') "ERROR in +cd umlalid1"
   write(lout,'(a)') "beam_expflag was", beam_expflag
   write(lout,'(a)') " expected 0 or 1. This is a BUG!"
   call prror(-1)
end if

if (.not.beam_expfile_open) then
  inquire(unit=600,opened=lopen)
  if (lopen) then
    write(lout,'(a)') "Error when opening beam_expert.txt"
    write(lout,'(a)') "Unit 600 already taken."
    call prror(-1)
  endif
#ifdef BOINC
  call boincrf("beam_expert.txt",filename)
  open(600,file=filename,status='replace',action="write")
#else
  open(600,file="beam_expert.txt",status='replace',action="write")
#endif
  beam_expfile_open = .true.
  !This line will be a comment if copy-pasted into fort.3
  write(600,'(a,g13.6,a,g13.6,a,g13.6,a)')"/ ******* USING emitx=",emitx,", emity=",emity,", emitz=",emitz," ******"
endif

if(parbe(ix,2).eq.0.0) then !4D
  !Note: One should always use the CRLIBM version when converting,
  ! in order to guarantee the exact same results from the converted input file.
#ifndef CRLIBM
  write(600,"(a16,1x,a1,1x,5g30.20)") bez(ix), "0", bbcu(ibb,1),bbcu(ibb,2),parbe(ix,5), parbe(ix,6), ptnfac(ix)
#else
  l1 = 1
  write(ch,'(a16,1x,a1)') bez(ix), "0"
  l1 = len(trim(ch))+1

  errno=dtostr(bbcu(ibb,1),ch1) ! Return value is the string length (always 24)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,2),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(parbe(ix,5),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(parbe(ix,6),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(ptnfac(ix),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  write(600,*) trim(ch)
#endif
else                      ! 6D
#ifndef CRLIBM
  write(600,"(a16,1x,i4,1x,4g30.20)") bez(ix), int(parbe(ix,2)),parbe(ix,1), parbe(ix,3),parbe(ix,5), parbe(ix,6)
  write(600,"(5g30.20)") bbcu(ibb,1), bbcu(ibb,4), bbcu(ibb,6), bbcu(ibb,2), bbcu(ibb,9)
  write(600,"(6g30.20)") bbcu(ibb,10), bbcu(ibb,3), bbcu(ibb,5), bbcu(ibb,7), bbcu(ibb,8), ptnfac(ix)
#else
  l1 = 1
  write(ch,'(a16,1x,i4)') bez(ix), int(parbe(ix,2))
  l1 = len(trim(ch))+1

  errno=dtostr(parbe(ix,1),ch1) ! Return value is the string length (always 24)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(parbe(ix,3),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(parbe(ix,5),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(parbe(ix,6),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  write(600,*) trim(ch)

  l1 = 1
  ch = ' '

  errno=dtostr(bbcu(ibb,1),ch1) ! Return value is the string length (always 24)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,4),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,6),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,2),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,9),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  write(600,*) trim(ch)

  l1 = 1
  ch = ' '

  errno=dtostr(bbcu(ibb,10),ch1) ! Return value is the string length (always 24)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,3),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,5),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,7),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(bbcu(ibb,8),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  errno=dtostr(ptnfac(ix),ch1)
  ch(l1:l1+errno) = ch1(1:errno)
  l1 = l1+errno+1

  write(600,*) trim(ch)
#endif
endif !END if(parbe(ix,2).eq.0.0)

if((bbcu(ibb,1).le.pieni).or.(bbcu(ibb,2).le.pieni)) then
  call prror(88)
endif
if(ibbc.eq.1) then
  sfac1=bbcu(ibb,1)+bbcu(ibb,2)
  sfac2=bbcu(ibb,1)-bbcu(ibb,2)
  sfac2s=one
  if(sfac2.lt.zero) sfac2s=-one                            !hr08
  sfac3=sqrt(sfac2**2+(four*bbcu(ibb,3))*bbcu(ibb,3))          !hr03
  if(sfac3.gt.sfac1) call prror(103)
  sfac4=(sfac2s*sfac2)/sfac3                                   !hr03
  sfac5=(((-one*sfac2s)*two)*bbcu(ibb,3))/sfac3                !hr03
  sigman(1,ibb)=sqrt(((sfac1+sfac2*sfac4)+(two*bbcu(ibb,3))*sfac5)*half)    !hr03
  sigman(2,ibb)=sqrt(((sfac1-sfac2*sfac4)-(two*bbcu(ibb,3))*sfac5)*half)    !hr03
  bbcu(ibb,11)=sqrt(half*(one+sfac4))
  bbcu(ibb,12)=(-one*sfac2s)*sqrt(half*(one-sfac4))            !hr03
  if(bbcu(ibb,3).lt.zero) bbcu(ibb,12)=-one*bbcu(ibb,12)       !hr03
else
  bbcu(ibb,11)=one
  sigman(1,ibb)=sqrt(bbcu(ibb,1))
  sigman(2,ibb)=sqrt(bbcu(ibb,2))
endif
if(parbe(ix,2).gt.zero) then !6D -> convert units
  do ii=1,10
    bbcu(ibb,ii)=bbcu(ibb,ii)*c1m6
  enddo
endif
