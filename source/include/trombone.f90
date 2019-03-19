! The values are stored in the temp vector which are used for the multiplication.
temptr(1)=xv1(j)
temptr(2)=yv1(j)/moidpsv(j)
temptr(3)=xv2(j)
temptr(4)=yv2(j)/moidpsv(j)
temptr(5)=sigmv(j)
temptr(6)=((mtc(j)*ejv(j)-e0)/e0f)*c1e3*(e0/e0f)
! Adding the closed orbit. The previous values are stored in the temptr vector.
xv1(j)  = cotr(irrtr,1)
yv1(j)  = cotr(irrtr,2)
xv2(j)  = cotr(irrtr,3)
yv2(j)  = cotr(irrtr,4)
sigmv(j) = cotr(irrtr,5)
pttemp   = cotr(irrtr,6)

! Multiplying the arbitrary matrix to the coordinates.
do kxxa=1,6
  xv1(j)   =  xv1(j)+temptr(kxxa)*rrtr(irrtr,1,kxxa)
  yv1(j)   =  yv1(j)+temptr(kxxa)*rrtr(irrtr,2,kxxa)
  xv2(j)   =  xv2(j)+temptr(kxxa)*rrtr(irrtr,3,kxxa)
  yv2(j)   =  yv2(j)+temptr(kxxa)*rrtr(irrtr,4,kxxa)
  sigmv(j)  =  sigmv(j)+temptr(kxxa)*rrtr(irrtr,5,kxxa)
  pttemp    =  pttemp+temptr(kxxa)*rrtr(irrtr,6,kxxa)
enddo
! Transforming back to the tracked coordinates of Sixtrack...
ejv(j)  = (e0f*pttemp/(c1e3*(e0/e0f))+e0)/mtc(j)
call part_updatePartEnergy(1,.false.)

! We have to go back to angles after we updated the energy.
yv1(j) = yv1(j)*moidpsv(j)
yv2(j) = yv2(j)*moidpsv(j)
