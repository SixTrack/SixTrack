! start include/wirekick.f90
!----------------------------------------------------------------------
! Wire element.
! MODEL OF STRAIGHT CURRENT WIRE
!
!     The model provides a transfer map of a straight current wire.
!     Description:
!     1. Infinitly thin wire with arbitrary orientation.
!     2. Thin element in SixTrack (L)=0
!     3. Parameters:
!     dx, dy: horizontal and vertical distances between wire midpoint
!     and closed orbit [mm]
!     (parameters are given by dx and dy in WIRE block)
!     tx, ty: tilt of the wire w.r.t the closed orbit in the
!     horizontal and vertical planes (in degrees)
!     (parameters are given by tiltx and tilty in WIRE block)
!     L - physical length of the wire element [m]
!     cur - current of the wire [Amperes]
!     embl - embedding drift (integrated length or integration interval) [m]
!     4. The transport map is given for canonical variables (x,px...)
!
! The MAP is constructed out of the following steps:
!     1. Declaration of shifted canonical variables:
!          rx = x+dx; ry = y+dy  in the same way as for the BEAM-BEAM element
!     2. Symplectic Rotation by the tilt angles tx, ty (in 4D space: px, rx, py, ry)
!     3. Wire kick for a longitudinally aligned wire (= kick for tx=ty=0)
!     4. Symplectic Rotation back by the tilt angles -ty, -yx (in 4D space: ...taking only PX, PY)
!--------------------------------------------------------------
!     Normalization factor (in SI) NNORM = (mu0*I*e)/(4*Pi*P0)
!     e -> 1; m0/4Pi -> 1.0e-7; N -> 1.0e-7*I

tx   = wire_tiltx(ix) ! tilt x [degrees]
ty   = wire_tilty(ix) ! tilt y [degrees]
tx   = tx*rad         ! [rad]
ty   = ty*rad         ! [rad]
dx   = wire_dispx(ix) ! displacement x [mm]
dy   = wire_dispy(ix) ! displacement y [mm]
embl = wire_lint(ix)  ! integrated length [m]
l    = wire_lphys(ix) ! physical length [m]
cur  = wire_current(ix)

if(abs(wire_flagco(ix)) /= 1) then
  write(lerr,"(a)") "WIRE> ERROR Wirekick in WIRE block must be either 1 or -1 for '"//trim(bez(ix))//"'"
  call prror
end if

if(wire_flagco(ix) == 1) then
  dxi = (dx+wire_clo(1,wire_num(i)))*c1m3
  dyi = (dy+wire_clo(2,wire_num(i)))*c1m3
else if(wire_flagco(ix) == -1) then
  dxi = dx*c1m3
  dyi = dy*c1m3
end if

do j=1,napx
  chi = (sqrt(e0**2-nucm(j)**2)*c1e6)/clight
  NNORM=c1m7/chi
  yv1(j) = yv1(j)*c1m3 ! [m]
  yv2(j) = yv2(j)*c1m3 ! [m]

  ! 1 shift
  if(wire_flagco(ix) == 1) then
    xi = (xv1(j)+dx)*c1m3 ! [m]
    yi = (xv2(j)+dy)*c1m3 ! [m]
  else if(wire_flagco(ix) == -1) then
    xi = (xv1(j)+( dx-wire_clo(1,wire_num(i)) ))*c1m3 ! [m]
    yi = (xv2(j)+( dy-wire_clo(2,wire_num(i)) ))*c1m3 ! [m]
  end if

  ! x'-> px; y'->py
  yv1(j) = yv1(j)*(one + dpsv(j))/mtc(j)
  yv2(j) = yv2(j)*(one + dpsv(j))/mtc(j)

  if(ibeco == 0) then
    ! 2 symplectic rotation of coordinate system (tx, ty)
    yi = yi-(((xi*sin_mb(tx))*yv2(j))/sqrt((one+dpsv(j))**2-yv2(j)**2))/&
          cos_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx)
    xi = xi*(cos_mb(tx)-sin_mb(tx)*tan_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx))
    yv1(j) = sqrt((one+dpsv(j))**2-yv2(j)**2)*sin_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx)

    xi = xi-(((yi*sin_mb(ty))*yv1(j))/sqrt((one+dpsv(j))**2-yv1(j)**2))/&
          cos_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty)
    yi = yi*(cos_mb(ty)-sin_mb(ty)*tan_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty))
    yv2(j) = sqrt((one+dpsv(j))**2-yv1(j)**2)*sin_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty)

    ! 3 apply wire kick
    RTWO = xi**2+yi**2
    yv1(j) = yv1(j)-(((CUR*NNORM)*xi)*(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO)))/RTWO
    yv2(j) = yv2(j)-(((CUR*NNORM)*yi)*(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO)))/RTWO

  elseif(ibeco == 1) then
    ! 2 symplectic rotation of coordinate system (tx, ty)

    dyi = dyi-(((dxi*sin_mb(tx))*yv2(j))/sqrt((one+dpsv(j))**2-yv2(j)**2))/&
          cos_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx)
    dxi = dxi*(cos_mb(tx)-sin_mb(tx)*tan_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx))

    yi = yi-(((xi*sin_mb(tx))*yv2(j))/sqrt((one+dpsv(j))**2-yv2(j)**2))/&
          cos_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx)
    xi = xi*(cos_mb(tx)-sin_mb(tx)*tan_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx))

    yv1(j) = sqrt((one+dpsv(j))**2-yv2(j)**2)*sin_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx)

    dxi = dxi-(((dyi*sin_mb(ty))*yv1(j))/sqrt((one+dpsv(j))**2-yv1(j)**2))/&
          cos_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty)
    dyi = dyi*(cos_mb(ty)-sin_mb(ty)*tan_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty))

    xi = xi-(((yi*sin_mb(ty))*yv1(j))/sqrt((one+dpsv(j))**2-yv1(j)**2))/&
          cos_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty)
    yi = yi*(cos_mb(ty)-sin_mb(ty)*tan_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty))

    yv2(j) = sqrt((one+dpsv(j))**2-yv1(j)**2)*sin_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty)

    ! 3 apply wire kick
    RTWO = xi**2+yi**2
    yv1(j) = yv1(j)-(((CUR*NNORM)*xi)*(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO)))/RTWO
    yv2(j) = yv2(j)-(((CUR*NNORM)*yi)*(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO)))/RTWO

    ! subtract closed orbit kick
    ! wire kick is negative px -> px - wirekick - (-closed orbit kick)
    RTWO = dxi**2+dyi**2
    yv1(j) = yv1(j)+(((CUR*NNORM)*dxi)*(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO)))/RTWO
    yv2(j) = yv2(j)+(((CUR*NNORM)*dyi)*(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO)))/RTWO

  end if

  ! 4 SYMPLECTIC ROTATION OF COORDINATE SYSTEM (-ty, -tx)
  yv2(j) = sqrt((one+dpsv(j))**2-yv1(j)**2)*sin_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))+ty)
  yv1(j) = sqrt((one+dpsv(j))**2-yv2(j)**2)*sin_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))+tx)

  ! px -> x'; py -> y'
  yv1(j) = yv1(j)*mtc(j)/(one + dpsv(j))
  yv2(j) = yv2(j)*mtc(j)/(one + dpsv(j))

  ! END OF WIRE MAP
  yv1(j) = yv1(j)*c1e3
  yv2(j) = yv2(j)*c1e3
end do

! end include/wirekick.f90
