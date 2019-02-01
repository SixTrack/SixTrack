! start include/multl07e.f90
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr03
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr03
qu1=tiltck*qu-tiltsk*qv
qv=tiltck*qv+tiltsk*qu
qu=qu1
dyy11=tiltc(k)*dyy1-tilts(k)*dyy2
dyy2=tiltc(k)*dyy2+tilts(k)*dyy1
dyy1=dyy11
! end include/multl07e.f90
