
  !---- Zero the arrays

do j=1,napx 

    irrtr = irm_rf(ix)
    nordm=nmu_rf(irrtr)
    crabfreq = freq_rfm(irrtr)*c1e3

  
  krf = (((one/(clight*(e0f/e0)))*crabfreq)*two)*pi
  

    ! if zero:th order prepare
    ! if not zero change phase with a minus sign and no multiply
  !---- Prepare to calculate the kick and the matrix elements
#include "alignva.f90"
    x_t = crkve*c1m3 
    y_t = cikve*c1m3
    !---- Vector with strengths
    do iord = 1, nordm
      field_cos(1,iord) = (nzz(j)*(norrfamp(irrtr,iord)) * cos((norrfph(irrtr,iord)*twopi)  - krf * sigmv(j)))
      field_sin(1,iord) = (nzz(j)*(norrfamp(irrtr,iord)) * sin((norrfph(irrtr,iord)*twopi)  - krf * sigmv(j)))
      field_cos(2,iord) = (nzz(j)*(skrfamp(irrtr,iord))  * cos((skrfph(irrtr,iord)*twopi)   - krf * sigmv(j)))
      field_sin(2,iord) = (nzz(j)*(skrfamp(irrtr,iord))  * sin((skrfph(irrtr,iord)*twopi)   - krf * sigmv(j)))
    
    enddo
    
    Cp0 = zero 
    Sp1 = zero

    do iord = nordm, 1, -1
      Cp0 = (Cp0 * (x_t+(imag*y_t))) / (iord)     + field_cos(1,iord) + (imag*field_cos(2,iord));
      Sp1 = (Sp1 * (x_t+(imag*y_t))) / (iord+1)   + field_sin(1,iord) + (imag*field_sin(2,iord));
    enddo
    
    Sp1 = Sp1 * (x_t+imag*y_t);
    

    yv1(j) = yv1(j) -((REAL(Cp0)*c1e3)*moidpsv(j))
    yv2(j) = yv2(j) + ((AIMAG(Cp0)*c1e3)*moidpsv(j)) 
    ejv(j) = ejv(j)  - ((REAL(Sp1)*(c1e3*(e0f*(crabfreq*(two*pi))))))/clight


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

  end do
