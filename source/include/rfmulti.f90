
do j=1,napx

  irrtr = irm_rf(ix)
  nordm = nmu_rf(irrtr)
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
    field_cos(1,iord) = (nzz(j)*(norrfamp(irrtr,iord)) * cos_mb(((norrfph(irrtr,iord)*twopi))  - (krf * sigmv(j))))
    field_sin(1,iord) = (nzz(j)*(norrfamp(irrtr,iord)) * sin_mb(((norrfph(irrtr,iord)*twopi))  - (krf * sigmv(j))))
    field_cos(2,iord) = (nzz(j)*(skrfamp(irrtr,iord))  * cos_mb(((skrfph(irrtr,iord)*twopi))   - (krf * sigmv(j))))
    field_sin(2,iord) = (nzz(j)*(skrfamp(irrtr,iord))  * sin_mb(((skrfph(irrtr,iord)*twopi))   - (krf * sigmv(j))))
  end do

  Cp0 = zero
  Sp1 = zero

  do iord = nordm, 1, -1
    Cp0 = (((Cp0 * (x_t+(imag*y_t))) / (iord))     + (field_cos(1,iord)) + (imag*field_cos(2,iord)));
    Sp1 = (((Sp1 * (x_t+(imag*y_t))) / (iord+1))   + (field_sin(1,iord)) + (imag*field_sin(2,iord)));
  end do

  Sp1 = Sp1 * (x_t+imag*y_t);

  yv1(j) = yv1(j) - (((REAL(Cp0))*tiltc(i))  + ((AIMAG(Cp0)))*tilts(i))*(moidpsv(j)*c1e3)
  yv2(j) = yv2(j) + (((-REAL(Cp0))*tilts(i)) + ((AIMAG(Cp0)))*tiltc(i))*(moidpsv(j)*c1e3)
  ejv(j)   = ejv(j) - ((real(Sp1)*(c1e3*(e0f*(crabfreq*(two*pi))))))/clight

end do

call part_updatePartEnergy(1,.true.)
if(ithick == 1) call envarsv(dpsv,moidpsv,rvv,ekv)
