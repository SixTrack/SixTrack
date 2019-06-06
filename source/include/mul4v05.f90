! start include/mul4v05.f90

yv1j = (bbiv(1,i) + bbiv(2,i)*xlvj) + aaiv(2,i)*zlvj
yv2j = (aaiv(1,i) - bbiv(2,i)*zlvj) + aaiv(2,i)*xlvj

if(curveff .and. abs(dki(ix,3)) > pieni) then !check for length and that flag is set
  if (abs(dki(ix,1)) > pieni) then ! Effect from the horizontal bend
    yv1j = yv1j + bbiv(2,i) * (strack(i) * ((xlvj**2 - 0.5*zlvj**2)*c1m3))
    yv2j = yv2j - bbiv(2,i) * (strack(i) * (((xlvj*zlvj)*c1m3) + (((dki(ix,1)*c1m6)*zlvj**3) / 6.0)))
  else ! Effect from vertical bend
    yv1j = yv1j + bbiv(2,i) * (strack(i) * (((xlvj*zlvj)*c1m3) + (((dki(ix,1)*c1m6)*xlvj**3) / 6.0)))
    yv2j = yv2j - bbiv(2,i) * (strack(i) * ((zlvj**2 - 0.5*xlvj**2)*c1m3))
  endif
endif

crkve=xlvj
cikve=zlvj
! end include/mul4v05.f90
