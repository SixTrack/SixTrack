!GRD-2007 REQUIRED FOR RHIC BEAM-BEAM STUDIES
!GRD-2007 OPTION IS ACTIVATED IF "LHC" PARAMETER IN BEAM-BEAM
!GRD-2007 BLOCK IN FORT.3 IS SET TO 9
if(lhc.eq.9) then
  ! Write a line to fort.10 for all platforms (including LSF)
#ifdef BOINC
  if (.not.restart) then
    write(10,'(a10,a60)') 'title     ',sixtit(1:60)
    endfile (10,iostat=ierro)
    backspace (10,iostat=ierro)
    bnlrec=bnlrec+1
  endif
#else
  write(10,'(a60)') sixtit(1:60)
  endfile (10,iostat=ierro)
  backspace (10,iostat=ierro)
#endif
  call bnlrdis(20000)
  write(lout,*) 'Sample number 1'
  do j = 1, napx
    pstop(nlostp(j))=.false.
  enddo
  do j = 1, napx
    xv(1,j)  = c1e3*myx(j) +torbx(1)
    yv(1,j)  = c1e3*myxp(j)+torbxp(1)
    xv(2,j)  = c1e3*myy(j) +torby(1)
    yv(2,j)  = c1e3*myyp(j)+torbyp(1)
    sigmv(j) = mys(j)
    ejv(j)   = myp(j)
    ejfv(j)=sqrt(ejv(j)**2-nucm0**2)                         !hr03
    rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
    dpsv(j)=(ejfv(j)-e0f)/e0f
    oidpsv(j)=one/(one+dpsv(j))
    moidpsv(j)=mtc(j)/(one+dpsv(j))
    dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)*mtc(j)                      !hr03
    namepart(j) = samplenumber*100 + j
  enddo
endif
