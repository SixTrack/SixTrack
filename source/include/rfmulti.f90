
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
  if(newstyle_rf .eq. 1) then
    irrtr = irm_rf(ix)
    NORMAL = norrfamp(irrtr,:)
    SKEW = skrfamp(irrtr,:)
    PNL = norrfph(irrtr,:)
    PSL = skrfph(irrtr,:)
    nordm=nmu_rf(irrtr)
    crabfreq = freq_rfm(irrtr)

  else
    crabamp=ed(ix)*nzz(j)
    NORMAL = zero
    SKEW = zero
    PNL = zero
    PSL = zero

    !nord = max(nn, ns, n_ferr/2-1)
    if(isSkew .eq. 0) then 
      if(nordm .eq. 1) then
        pnl(1) = pi/2 - crabpase_t
        normal(1) = crabamp/e0f
      else
        pnl(nordm) =  -crabpase_t
        normal(nordm) = -crabamp*fact
      endif
    else
      if(nordm .eq. 1) then
        psl(1) = -pi/2 - crabpase_t
        skew(1) = crabamp/e0f
      else
        psl(nordm) = - crabpase_t
        skew(nordm) = crabamp*fact
      endif
    endif
  endif
  
  krf = (((one/(clight*(e0f/e0)))*crabfreq)*two)*pi
  



  print * , krf, crabfreq, crabamp, crabph(ix), pnl(1), normal(1), "ddddd"
  print *, psl(1), skew(1), "skeeewww11"

    ! if zero:th order prepare
    ! if not zero change phase with a minus sign and no multiply
  !---- Prepare to calculate the kick and the matrix elements
#include "alignva.f90"
    x_t = crkve*c1m3 
    y_t = cikve*c1m3
    !---- Vector with strengths

    do iord = 1, nordm
      print *, "sigmmmma", sigmv(j)
      print *, "pnnnllll", pnl(iord)
      print *, "streeeeength", normal(iord)
      field_cos(1,iord) = (normal(iord) * cos(pnl(iord)*twopi  - krf * sigmv(j)))
      field_sin(1,iord) = (normal(iord) * sin(pnl(iord)*twopi  - krf * sigmv(j)))
      field_cos(2,iord) = (skew(iord)   * cos(psl(iord)*twopi  - krf * sigmv(j)))
      field_sin(2,iord) = (skew(iord)   * sin(psl(iord)*twopi  - krf * sigmv(j)))
      print *, "streeeeength", field_cos(1,iord)
    enddo
    Cm2 = zero; Sm2a = zero; Cm1 = zero; Sm1a = zero;
    Cp0 = zero; Sp0 = zero; Cp1 = zero; Sp1 = zero;

    do iord = nordm, 1, -1

      Cp0 = Cp0 * (x_t+imag*y_t) / (iord)   + field_cos(1,iord)+imag*field_cos(2,iord);
      !Sp0 = Sp0 * (x_t+imag*y_t) / (iord+1)   + field_sin(1,iord)+imag*field_sin(2,iord);
      !Cp1 = Cp1 * (x_t+imag*y_t) / (iord+2)   + field_cos(1,iord)+imag*field_cos(2,iord);
      Sp1 = Sp1 * (x_t+imag*y_t) / (iord+1)   + field_sin(1,iord)+imag*field_sin(2,iord);
    enddo
    Sp1 = Sp1 * (x_t+imag*y_t);
    

    dpx = -REAL(Cp0)*c1e3*moidpsv(j);
    dpy = AIMAG(Cp0)*c1e3*moidpsv(j);
    dpt = - krf * REAL(Sp1)*c1e3*e0f;
    yv1(j) = yv1(j) + dpx
    yv2(j) = yv2(j) + dpy
    print *, "ejjjbefore", ejv(j)
    ejv(j) = ejv(j) + dpt
    print *, "ejjjafter", ejv(j)

    

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



