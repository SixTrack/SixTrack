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




      tx = wire_tiltx(ix) !tilt x [degrees]
      ty = wire_tilty(ix) !tilt y [degrees]
      tx = tx*(pi/c180e0) ![rad]
      ty = ty*(pi/c180e0) ![rad]
      dx = wire_dispx(ix) !displacement x [mm]
      dy = wire_dispy(ix) !displacement y [mm]
      embl = wire_lint(ix) !integrated length [m]
      l = wire_lphys(ix) !physical length [m]
      cur = wire_current(ix)

      if (abs(wire_flagco(ix)).ne.1) then
        write(lout,                                                     &
     &fmt='((A,A,/),(A,I0,A,A,/),(A,I0,A,I0,/))')                       &
     &'ERROR: in wirekick -  wire_flagco defined in WIRE block must ',  &
     &'be either 1 or -1! Did you define all wires in the WIRE block?', &
     &'bez(',ix,') = ',bez(ix),                                         &
     &'wire_flagco(',ix,') = ',wire_flagco(ix)
        call prror
      endif


      IF (wire_flagco(ix).eq.1) THEN
         dxi = (dx+wire_clo(1,wire_num(i)) )*c1m3
         dyi = (dy+wire_clo(2,wire_num(i)) )*c1m3
      ELSE IF (wire_flagco(ix).eq.-1) THEN
         dxi = (dx)*c1m3
         dyi = (dy)*c1m3
      END IF

      do j=1, napx
      chi = (sqrt(e0**2-nucm(j)**2)*c1e6)/clight
      NNORM=c1m7/chi
      yv1(j) = yv1(j) * c1m3 !SI
      yv2(j) = yv2(j) * c1m3 !SI

! 1 shift
      IF (wire_flagco(ix).eq.1) THEN
         xi = (xv1(j)+dx)*c1m3 !SI
         yi = (xv2(j)+dy)*c1m3 !SI
      ELSE IF (wire_flagco(ix).eq.-1) THEN
         xi = (xv1(j)+( dx-wire_clo(1,wire_num(i)) ))*c1m3 !SI
         yi = (xv2(j)+( dy-wire_clo(2,wire_num(i)) ))*c1m3 !SI
      END IF

! x'-> px; y'->py
      yv1(j) = yv1(j)*(one + dpsv(j))/mtc(j)
      yv2(j) = yv2(j)*(one + dpsv(j))/mtc(j)

!ibeco = 0
      if(ibeco.eq.0) then
! 2 symplectic rotation of coordinate system (tx, ty)
          yi = yi-(((xi*sin_mb(tx))*yv2(j))/                           &
     &sqrt((one+dpsv(j))**2-yv2(j)**2))/                               &
     &cos_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-        &
     &yv2(j)**2))-tx)
          xi = xi*(cos_mb(tx)-sin_mb(tx)*tan_mb(atan_mb(yv1(j)/        &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx))
          yv1(j) = sqrt((one+dpsv(j))**2-yv2(j)**2)*                  &
     &sin_mb(atan_mb(yv1(j)/                                           &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx)

          xi = xi-(((yi*sin_mb(ty))*yv1(j))/                           &
     &sqrt((one+dpsv(j))**2-yv1(j)**2))/                               &
     &cos_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-        &
     &yv2(j)**2))-ty)
          yi = yi*(cos_mb(ty)-sin_mb(ty)*tan_mb(atan_mb(yv2(j)/        &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty))
          yv2(j) = sqrt((one+dpsv(j))**2-yv1(j)**2)*                  &
     &sin_mb(atan_mb(yv2(j)/                                           &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty)

! 3 apply wire kick
          RTWO = xi**2+yi**2
          yv1(j) = yv1(j)-(((CUR*NNORM)*xi)*                          &
     &(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO) ))/RTWO
          yv2(j) = yv2(j)-(((CUR*NNORM)*yi)*                          &
     &(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO) ))/RTWO

! ibeco = 1
      elseif(ibeco.eq.1) then
! 2 symplectic rotation of coordinate system (tx, ty)

          dyi = dyi-(((dxi*sin_mb(tx))*yv2(j))/                        &
     &sqrt((one+dpsv(j))**2-yv2(j)**2))/                               &
     &cos_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-        &
     &yv2(j)**2))-tx)
          dxi = dxi*(cos_mb(tx)-sin_mb(tx)*tan_mb(atan_mb(yv1(j)/      &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx))

          yi = yi-(((xi*sin_mb(tx))*yv2(j))/                           &
     &sqrt((one+dpsv(j))**2-yv2(j)**2))/                               &
     &cos_mb(atan_mb(yv1(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-        &
     &yv2(j)**2))-tx)
          xi = xi*(cos_mb(tx)-sin_mb(tx)*tan_mb(atan_mb(yv1(j)/        &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx))

          yv1(j) = sqrt((one+dpsv(j))**2-yv2(j)**2)*                  &
     &sin_mb(atan_mb(yv1(j)/                                           &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-tx)

          dxi = dxi-(((dyi*sin_mb(ty))*yv1(j))/                        &
     &sqrt((one+dpsv(j))**2-yv1(j)**2))/                               &
     &cos_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-        &
     &yv2(j)**2))-ty)
          dyi = dyi*(cos_mb(ty)-sin_mb(ty)*tan_mb(atan_mb(yv2(j)/      &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty))

          xi = xi-(((yi*sin_mb(ty))*yv1(j))/                           &
     &sqrt((one+dpsv(j))**2-yv1(j)**2))/                               &
     &cos_mb(atan_mb(yv2(j)/sqrt(((one+dpsv(j))**2-yv1(j)**2)-        &
     &yv2(j)**2))-ty)
          yi = yi*(cos_mb(ty)-sin_mb(ty)*tan_mb(atan_mb(yv2(j)/        &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty))

         yv2(j) = sqrt((one+dpsv(j))**2-yv1(j)**2)*                   &
     &sin_mb(atan_mb(yv2(j)/                                           &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))-ty)

! 3 apply wire kick
          RTWO = xi**2+yi**2
          yv1(j) = yv1(j)-(((CUR*NNORM)*xi)*                          &
     &(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO) ))/RTWO
          yv2(j) = yv2(j)-(((CUR*NNORM)*yi)*                          &
     &(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO) ))/RTWO
! subtract closed orbit kick
! wire kick is negative px -> px - wirekick - (-closed orbit kick)
          RTWO = dxi**2+dyi**2
          yv1(j) = yv1(j)+(((CUR*NNORM)*dxi)*                         &
     &(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO) ))/RTWO
          yv2(j) = yv2(j)+(((CUR*NNORM)*dyi)*                         &
     &(sqrt((embl+L)**2+four*RTWO)-sqrt((embl-L)**2+four*RTWO) ))/RTWO

      endif

! 4 SYMPLECTIC ROTATION OF COORDINATE SYSTEM (-ty, -tx)
      yv2(j) = sqrt((one+dpsv(j))**2-yv1(j)**2)*                      &
     &sin_mb(atan_mb(yv2(j)/                                           &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))+ty)
      yv1(j) = sqrt((one+dpsv(j))**2-yv2(j)**2)*                      &
     &sin_mb(atan_mb(yv1(j)/                                           &
     &sqrt(((one+dpsv(j))**2-yv1(j)**2)-yv2(j)**2))+tx)

! px -> x'; py -> y'
      yv1(j) = yv1(j)*mtc(j)/(one + dpsv(j))
      yv2(j) = yv2(j)*mtc(j)/(one + dpsv(j))
!-----------------------------------------------------------------------
! END OF WIRE MAP
!-----------------------------------------------------------------------
      yv1(j) = yv1(j) * c1e3
      yv2(j) = yv2(j) * c1e3
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
! end include/wirekick.f90
