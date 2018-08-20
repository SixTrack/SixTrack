! start include/kickq04v.f90
#ifndef TILT
qu=(three*ekk)*cikve                                             !hr02
qv=((-one*three)*ekk)*crkve                                      !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(three*ekk)*(tiltck*cikve-tiltsk*crkve)                       !hr02
qv=((-one*three)*ekk)*(tiltck*crkve+tiltsk*cikve)                !hr02
#endif
! end include/kickq04v.f90
