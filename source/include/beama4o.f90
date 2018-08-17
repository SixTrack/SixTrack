! start include/beama4o.f90
  beamoff(4,imbb(i))=(rkb(j)*(crzb(j)-exp_mb(-one*tkb(j))*cbzb(j)))*sign(one,crkveb(j))
  beamoff(5,imbb(i))=(rkb(j)*(crxb(j)-exp_mb(-one*tkb(j))*cbxb(j)))*sign(one,cikveb(j))
end do
! end include/beama4o.f90
