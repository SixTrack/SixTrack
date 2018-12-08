
  !---- Zero the arrays
  NORMAL = zero
  SKEW = zero
  PNL = zero
  PSL = zero


fact=1   
do j=1,nordm-1   
  fact=fact*j   
end do


  do j=1,napx 
  
  crabamp=ed(ix)*nzz(j)
  NORMAL = zero
  SKEW = zero
  PNL = zero
  PSL = zero
  
  !nord = max(nn, ns, n_ferr/2-1)


  if(isSkew .eq. 0) then 
    if(nordm .eq. 1) then
      print *, "are we hereeee..?"
      pnl(0) = pi/2 - crabpase_t
      normal(0) = crabamp/e0f
    else
      pnl(nordm-1) =  -crabpase_t
      normal(nordm-1) = -crabamp*fact
    endif
  else
    if(nordm .eq. 1) then
      psl(0) = -pi/2 - crabpase_t
      skew(0) = crabamp/e0f
    else
      psl(nordm-1) = - crabpase_t
      skew(nordm-1) = crabamp*fact
    endif
  endif
  krf = (((one/(clight*(e0f/e0)))*crabfreq)*two)*pi
  


  print * , "krf",  "crabfreq", "crabamp",  "phase", pnl(0), normal(0)
  print * , krf, crabfreq, crabamp, crabph(ix), pnl(1), normal(1), "ddddd"
  print *, psl(1), skew(1), "skeeewww11"

    ! if zero:th order prepare
    ! if not zero change phase with a minus sign and no multiply
  !---- Prepare to calculate the kick and the matrix elements
    x_t = xv1(j)*c1m3 
    y_t = xv2(j)*c1m3
    !---- Vector with strengths 
    do iord = 0, nordm
      print *,"skewww",   krf, skew(0), psl(0)
      field_cos(1,iord) = (normal(iord) * cos(pnl(iord)  - krf * sigmv(j)))
      field_sin(1,iord) = (normal(iord) * sin(pnl(iord)  - krf * sigmv(j)))
      field_cos(2,iord) = (skew(iord)   * cos(psl(iord)  - krf * sigmv(j)))
      field_sin(2,iord) = (skew(iord)   * sin(psl(iord)  - krf * sigmv(j)))
    enddo
    print *, "is fieeeld", field_cos(2,0), field_sin(2,0), nord
    Cm2 = zero; Sm2a = zero; Cm1 = zero; Sm1a = zero;
    Cp0 = zero; Sp0 = zero; Cp1 = zero; Sp1 = zero;

    do iord = nordm, 0, -1
    print * , "iiooord", iord, nordm

      Cp0 = Cp0 * (x_t+imag*y_t) / (iord+1)   + field_cos(1,iord)+imag*field_cos(2,iord);
      Sp1 = Sp1 * (x_t+imag*y_t) / (iord+2)   + field_sin(1,iord)+imag*field_sin(2,iord);
      !Cp0 = Cp0 * (x_t+imag*y_t) / (iord+1)  ! + field_cos(1,iord)+imag*field_cos(2,iord);
      !Sp1 = Sp1 * (x_t+imag*y_t) / (iord+2)  ! + field_sin(1,iord)+imag*field_sin(2,iord);

      
    enddo
    Sp1 = Sp1 * (x_t+imag*y_t);
    

    dpx = -REAL(Cp0)*c1e3*moidpsv(j);
    dpy = AIMAG(Cp0)*c1e3*moidpsv(j);
    dpt = - krf * REAL(Sp1)*c1e3*e0f;
    yv1(j) = yv1(j) + dpx
    yv2(j) = yv2(j) + dpy
    ejv(j)=ejv(j)+dpt
    

    ejf0v(j)=ejfv(j)
    ejfv(j)=sqrt(ejv(j)**2-nucm(j)**2)
    rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
    dpsv(j)=(ejfv(j)*(nucm0/nucm(j))-e0f)/e0f
    oidpsv(j)=one/(one+dpsv(j))
    moidpsv(j)=mtc(j)/(one+dpsv(j))
    omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
    dpsv1(j)=(dpsv(j)*c1e3)*oidpsv(j)
    yv1(j)=(ejf0v(j)/ejfv(j))*yv1(j)
    yv2(j)=(ejf0v(j)/ejfv(j))*yv2(j)
    if(ithick.eq.1) call envarsv(dpsv,moidpsv,rvv,ekv)

    print *, "aaaaa", dpx, dpy, dpt, Sp1, imag
  end do



