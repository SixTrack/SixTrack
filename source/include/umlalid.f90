iwrite=0
if(nlin.eq.0) then
  iwrite=1
else
  do ii=1,nlin
    if(typ.eq.bezl(ii)) iwrite=1
  enddo
endif
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
  if(iflag1.eq.1.and.ithick.eq.1) then
    call dacct(damap,nvar,corrnew,nvar,damap,nvar)
  endif
else
  call dacop(dpda1,damap(5))
  do j1=1,4
    do ii=1,4
      jj(ii)=1
      call dapek(damap(j1),jj,rdd(j1,ii))
      jj(ii)=0
    enddo
  enddo
  jj(5)=1
  do j1=1,4
    call dapek(damap(j1),jj,rdd(j1,5))
  enddo
  jj(5)=0
  do j1=1,2
    ii=2*j1
    d(j1)=(((rdd(ii-1,1)*dicu(1)+rdd(ii-1,2)*dicu(2))+rdd(ii-1,3)*dicu(3))+rdd(ii-1,4)*dicu(4))+rdd(ii-1,5)              !hr03
    dp(j1)=(((rdd(ii,1)*dicu(1)+rdd(ii,2)*dicu(2))+rdd(ii,3)*dicu(3))+rdd(ii,4)*dicu(4))+rdd(ii,5)                    !hr03
  enddo
endif
call dacct(damap,nvar,aa2,nvar,damap,nvar)
!         calculate linear 6D optics parameter for each element
!         by calculating the matrix of eigenvectors (tas)
do j=1,ndimf
  ii=2*j
  if(j.eq.1) then
    i2=4
    i3=6
  elseif(j.eq.2) then
    i2=2
    i3=6
  elseif(j.eq.3) then
    i2=2
    i3=4
  endif
  jj(ii-1)=1
  call dapek(damap(ii-1),jj,angp(1,ii-1))
  call dapek(damap(ii),jj,au(ii,ii-1))
  jj(ii-1)=0
  jj(ii)=1
  call dapek(damap(ii-1),jj,angp(1,ii))
  call dapek(damap(ii),jj,au(ii,ii))
  jj(ii)=0
  jj(i2-1)=1
  call dapek(damap(ii-1),jj,au(i2-1,i2-1))
  call dapek(damap(ii),jj,au(i2,i2-1))
  jj(i2-1)=0
  jj(i2)=1
  call dapek(damap(ii-1),jj,au(i2-1,i2))
  call dapek(damap(ii),jj,au(i2,i2))
  jj(i2)=0
  jj(i3-1)=1
  call dapek(damap(ii-1),jj,au(i3-1,i3-1))
  call dapek(damap(ii),jj,au(i3,i3-1))
  jj(i3-1)=0
  jj(i3)=1
  call dapek(damap(ii-1),jj,au(i3-1,i3))
  call dapek(damap(ii),jj,au(i3,i3))
  jj(i3)=0

!     Store tas matrix (normalisation of phase space) and closed orbit for FMA and DUMP normalization.
!     Variable added to DUMP block module variables;
!     units dumptas: mm,mrad,mm,mrad,mm,1.e-3 -> convert later to 1.e3
  if(ic(i)-nblo.gt.0) then !check if structure element is a block
    if(ldump(ic(i)-nblo)) then !check if particles are dumped at this element
      dumptas(ic(i)-nblo,ii-1,ii-1)=angp(1,ii-1)
      dumptas(ic(i)-nblo,ii-1,ii  )=angp(1,ii)
      dumptas(ic(i)-nblo,ii  ,ii-1)=au(ii,ii-1)
      dumptas(ic(i)-nblo,ii  ,ii  )=au(ii,ii  )
      dumptas(ic(i)-nblo,ii-1,i2-1)=au(i2-1,i2-1)
      dumptas(ic(i)-nblo,ii  ,i2-1)=au(i2  ,i2-1)
      dumptas(ic(i)-nblo,ii-1,i2  )=au(i2-1,i2  )
      dumptas(ic(i)-nblo,ii  ,i2  )=au(i2  ,i2  )
      dumptas(ic(i)-nblo,ii-1,i3-1)=au(i3-1,i3-1)
      dumptas(ic(i)-nblo,ii  ,i3-1)=au(i3  ,i3-1)
      dumptas(ic(i)-nblo,ii-1,i3  )=au(i3-1,i3  )
      dumptas(ic(i)-nblo,ii  ,i3  )=au(i3  ,i3  )
!    closed orbit in canonical variables x,px,y,py,sig,delta [mm,mrad,mm,mrad,mm,1.e-3]
!    convert to x,xp,y,yp,sig,delta [mm,mrad,mm,mrad,mm,1]
!     -> check units used in dumpclo (is x' or px used?)
      dumpclo(ic(i)-nblo,2*j-1)=c(j)
      if (j.eq.3) then !dp/p
        dumpclo(ic(i)-nblo,2*j)  =cp(j)*c1m3
      else ! xp,yp
        dumpclo(ic(i)-nblo,2*j)  =cp(j)/(one+cp(3)*c1m3)
      endif
    endif
  endif

  b1(j)=angp(1,ii-1)**2+angp(1,ii)**2                          !hr08
  b2(j)=au(i2-1,i2-1)**2+au(i2-1,i2)**2                        !hr08
  b3(j)=au(i3-1,i3-1)**2+au(i3-1,i3)**2                        !hr08
  al1(j)=-one*(angp(1,ii-1)*au(ii,ii-1)+angp(1,ii)*au(ii,ii))  !hr03
  al2(j)=-one*(au(i2-1,i2-1)*au(i2,i2-1)+au(i2-1,i2)*au(i2,i2)) !hr03
  al3(j)=-one*(au(i3-1,i3-1)*au(i3,i3-1)+au(i3-1,i3)*au(i3,i3)) !hr03
  g1(j)=au(ii,ii-1)**2+au(ii,ii)**2                            !hr04
  g2(j)=au(i2,i2-1)**2+au(i2,i2)**2                            !hr04
  g3(j)=au(i3,i3-1)**2+au(i3,i3)**2                            !hr04
  if(ndimf.eq.3) then
    call dainv(damap,nvar,damapi,nvar)
    jj(6)=1
    call dapek(damapi(5),jj,aui(1))
    call dapek(damapi(6),jj,aui(2))
    jj(6)=0
    if(j.lt.3) then
      d(j)=au(i3-1,i3-1)*aui(1)+au(i3-1,i3)*aui(2)
      dp(j)=au(i3,i3-1)*aui(1)+au(i3,i3)*aui(2)
    else
      d(j)=angp(1,ii-1)*aui(1)+angp(1,ii)*aui(2)
      dp(j)=au(ii,ii-1)*aui(1)+au(ii,ii)*aui(2)
    endif
  endif
  sx=angp(2,ii-1)*angp(1,ii)-angp(1,ii-1)*angp(2,ii)
  cx=angp(1,ii-1)*angp(2,ii-1)+angp(1,ii)*angp(2,ii)
  if(abs(sx).gt.c1m15.or.abs(cx).gt.c1m15) then
    dphi(j)=atan2_mb(sx,cx)/x2pi
  else
    dphi(j)=zero
  endif
  phi(j)=phi(j)+dphi(j)
enddo ! end include/of optics calculation

if(ic(i)-nblo.gt.0) then !check if structure element is a block
   if(ldump(ic(i)-nblo)) then !check if particles are dumped at this element

!     do the unit conversion + inversion of dumptas
!     convert from units [mm,mrad,mm,mrad,1.e-3] to [mm,mrad,mm,mrad,1] as needed for normalization

     dumptas(ic(i)-nblo,1:5,6)=dumptas(ic(i)-nblo,1:5,6)*c1e3
     dumptas(ic(i)-nblo,6,1:5)=dumptas(ic(i)-nblo,6,1:5)*c1m3

!     invert the tas matrix
      call invert_tas(dumptasinv(ic(i)-nblo,:,:),dumptas(ic(i)-nblo,:,:))
!     dumptas and dumptasinv are now in units [mm,mrad,mm,mrad,1]

   endif
endif

do j=1,ndimf
  ii=2*j
  angp(2,ii-1)=angp(1,ii-1)
  angp(2,ii)=angp(1,ii)
enddo
!         write optics parameter for each element (LINE block)
if(iwrite.eq.1) then
  iii=i
  if(typ(:8).eq.'START   ') iii=0
  write(lout,10030) iii,typ(:8),tl,phi(1),b1(1),al1(1),g1(1),d(1),dp(1),c(1),cp(1)
  if(ndimf.eq.3) then
    write(lout,10040) b2(1),al2(1),g2(1)
    write(lout,10050) typ(9:16),b3(1),al3(1),g3(1)
  else
    write(lout,10055) typ(9:16),b2(1),al2(1),g2(1)
  endif
  write(lout,10060)
  write(lout,10070) phi(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),cp(2)
  write(lout,10080) b2(2),al2(2),g2(2)
  if(ndimf.eq.3) then
    write(lout,10090) b3(2),al3(2),g3(2)
    write(lout,10060)
    write(lout,10100) -phi(3),b1(3),al1(3),g1(3),d(3),dp(3),c(3),cp(3)
    write(lout,10080) b2(3),al2(3),g2(3)
    write(lout,10040) b3(3),al3(3),g3(3)
  endif
  write(lout,10010)
endif
